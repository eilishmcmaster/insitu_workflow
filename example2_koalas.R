library(ozmaps)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(sf)
library(ggpubr)
library(units)

# Load data ####
koalas <- read.csv2('data/ALA_koala_records.tsv', sep='\t')
koalas <- koalas[koalas$state=='New South Wales',]
koalas_sf <- st_as_sf(koalas, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>% st_transform(crs=3577)

nsw_map <- ozmaps::ozmap_states %>%
  filter(NAME == "New South Wales") %>%
  st_transform(crs = 3577)

## Define input variables #### 
## Set parameters for cluster creation ####
min_n <- 100  # Minimum number of observations in a cluster
max_area <-  2e+08  # Maximum area (m2) for clusters = 200 km2
eps_value <- 5000  # DBSCAN epsilon (maximum distance for a point to be from other points to be included in a cluster)
expand_m <- 100  # Buffer distance for convex hulls
# because this data is in GDA94, one unit is equivalent to 1 meter

# Set the default theme to theme_bw and remove x and y axis text and titles
theme_set(theme_bw())  # Set default theme to theme_bw
theme_update(
  axis.title = element_blank(),       # Remove x and y axis titles
  theme(legend.key.size = unit(0.5, 'lines'), legend.key.height = unit(0.5, 'lines'))
)

# Initial plot of individual distrubutions ####
# Plot vegetation map with points
initial_map_plot <- ggplot() +
  geom_sf(data = ozmaps::ozmap_states, fill = "gray", color = "gray", size = 0.2, alpha=0.5) +
  geom_point(
    color='steelblue3',
    data = koalas, 
    aes(x = as.numeric(decimalLongitude), y = as.numeric(decimalLatitude)),
    alpha = 0.05, 
    size = 0.2
  ) +
  xlim(113, 153)+coord_sf()+
  geom_sf(data = nsw_map, fill = "transparent", color = "black", size = 0.2) 

# Spatial clustering method: defining in situ management sites  #################################################################
# Step 1: DBSCAN clusters points based on density.
# Step 2: Large DBSCAN clusters are refined using K-means clustering if their area exceeds max_area
## 1. Extract coordinates and calculate distance matrix
coords_df <- st_coordinates(koalas_sf) # Extract coordinates from spatial data
pairwise_geo_dist <- dist(coords_df, method = 'maximum') # Calculate pairwise geographic distances

# ## 2. Plot k-NN distances for clustering
# kNNdistplot(pairwise_geo_dist, k = 6)
# abline(h = 125, lty = 2)

## 3. Perform DBSCAN clustering
# DBSCAN identifies density-based clusters using the distance matrix
db <- dbscan::dbscan(pairwise_geo_dist, eps = eps_value, minPts = 100)
cluster_df <- as.data.frame(coords_df)
cluster_df$cluster <- as.factor(db$cluster)
cluster_df$cluster[cluster_df$cluster == "0"] <- NA # Assign NA to noise points

## 4. Filter small clusters and renumber them
# Remove clusters smaller than min_n and renumber remaining clusters
cluster_df <- cluster_df %>%
  mutate(cluster = ifelse(cluster %in% names(which(table(cluster) < min_n)), NA, cluster)) %>%
  mutate(cluster = as.numeric(as.factor(cluster)))

# Add species and sample metadata
cluster_df$species <- koalas_sf$species
cluster_df$sample <- koalas_sf$sample

## 5. Compute convex hulls for visualization
# Generate convex hulls for each cluster for visualization purposes
hull_sf <- cluster_df[!is.na(cluster_df$cluster), ] %>%
  group_by(cluster) %>%
  slice(chull(X, Y)) %>%
  st_as_sf(coords = c("X", "Y"),crs=3577) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_convex_hull() %>%
  filter(!st_geometry_type(geometry) %in% c("POINT")) %>%
  st_buffer(dist = expand_m) %>%
  mutate(perimeter = st_perimeter(geometry), area = st_area(geometry))

## 6. Save initial plot with convex hulls
bbox <- st_bbox(hull_sf)
plot_list <- list()
plot_list[[length(plot_list) + 1]] <- ggplot() +
  geom_sf(data = nsw_map, color = "black", fill='grey', size = 0.1, alpha=0.5) +
  geom_point(data = cluster_df, aes(X, Y), color = 'steelblue3', alpha = 0.05, size = 0.2) +
  geom_sf(data = st_transform(hull_sf, crs = 4326),
          fill = "transparent", color = 'red', linewidth = 0.3)+
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), expand = 20000)

## 7. Identify large clusters and refine with K-means
# Identify clusters with a area exceeding max_area
large_clusters <- hull_sf$cluster[which(hull_sf$area >= set_units(max_area, m^2))] %>% as.vector()
new_cluster_hull_sf <- hull_sf

while (length(large_clusters) > 0) {
  large_cluster_df <- cluster_df[cluster_df$cluster %in% large_clusters, ]
  new_cluster_df <- data.frame()
  
  # Apply K-means to split each large cluster
  for (original_big_cluster in large_clusters) {
    subset_df <- large_cluster_df[large_cluster_df$cluster == original_big_cluster, ]
    p <- new_cluster_hull_sf$area[new_cluster_hull_sf$cluster == original_big_cluster]
    start_i <- ceiling(as.numeric(p / max_area))  # Ensure numeric value for kmeans centers
    
    # Ensure start_i does not exceed the number of points
    start_i <- min(start_i, nrow(subset_df))
    
    # Proceed only if there are enough points for clustering
    if (start_i > 1) {
      kmeans_result <- kmeans(subset_df[, 1:2], centers = start_i)
      subset_df$cluster <- paste(subset_df$cluster, kmeans_result$cluster, sep = '.')
    }
    new_cluster_df <- rbind(new_cluster_df, subset_df)
  }
  
  # Update cluster data and recompute convex hulls
  cluster_df <- rbind(cluster_df[!cluster_df$cluster %in% large_clusters, ], new_cluster_df)
  new_cluster_hull_sf <- cluster_df[!is.na(cluster_df$cluster), ] %>%
    group_by(cluster) %>%
    slice(chull(X, Y)) %>%
    st_as_sf(coords = c("X", "Y"), crs = 3577) %>%
    summarise(geometry = st_combine(geometry)) %>%
    st_convex_hull() %>%
    filter(!st_geometry_type(geometry) %in% c("POINT")) %>%
    st_buffer(dist = expand_m) %>%
    mutate(perimeter = st_perimeter(geometry), area = st_area(geometry))
  
  plot_list[[length(plot_list) + 1]] <- ggplot() +
    geom_sf(data = nsw_map, color = "black", fill='grey', size = 0.1, alpha=0.5) +
    geom_point(data = cluster_df, aes(X, Y), color = 'steelblue3', alpha = 0.05, size = 0.2) +
    geom_sf(data = st_transform(new_cluster_hull_sf, crs = 4326), fill = "transparent", color = 'red', linewidth = 0.3)+
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), expand = 20000)
  
  large_clusters <- new_cluster_hull_sf$cluster[which(new_cluster_hull_sf$area >= set_units(max_area, m^2))] %>% as.vector()
}

## 9. Finalize clusters and handle remaining NA values
# Assign new cluster IDs to remaining NA values
final_cluster_df <- cluster_df
max_cluster <- final_cluster_df %>%
  mutate(first_part = as.numeric(sub("\\..*", "", cluster))) %>%
  summarize(max_cluster = max(first_part, na.rm = TRUE)) %>%
  pull(max_cluster)

num_na <- sum(is.na(final_cluster_df$cluster))
new_cluster_values <- seq(from = max_cluster + 1, length.out = num_na)
final_cluster_df$cluster[is.na(final_cluster_df$cluster)] <- as.character(new_cluster_values)

# Recalculate convex hulls for final clusters
final_hull_sf <- final_cluster_df[!is.na(final_cluster_df$cluster), ] %>%
  group_by(cluster) %>%
  slice(chull(X, Y)) %>%
  st_as_sf(coords = c("X", "Y")) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_convex_hull() %>%
  filter(!st_geometry_type(geometry) %in% c("POINT")) %>%
  st_buffer(dist = expand_m) %>%
  mutate(perimeter = st_perimeter(geometry), area = st_area(geometry))

# Set CRS for final convex hull data
final_hull_sf <- st_set_crs(final_hull_sf, value = 3577)

# Arrange and display clustering plots
clustering_process_plots <- ggarrange(
  initial_map_plot + ggtitle('Koala observations in NSW 2022-25'),
  plot_list[[1]] + ggtitle('HDBSCAN clusters'),
  plot_list[[2]] + ggtitle('K-means clusters round 1'),
  # plot_list[[3]] + ggtitle('K-means clusters round 2'),
  common.legend = TRUE, legend = 'bottom',
  widths=c(1,0.8,0.8),
  ncol = 3, labels = 'auto'
)

clustering_process_plots

ggsave('outputs/example2_map_process.png', clustering_process_plots, width=25, height=10, units = 'cm', bg='white')
