# Load necessary libraries ####
library(fastDiversity) #https://github.com/eilishmcmaster/fastDiversity
library(dplyr)
library(sf)
library(openxlsx)
library(ggplot2)
library(raster)
library(plotly)
library(tidyverse)
library(dbscan)
library(ggpubr)
library(ggthemes)


# Load data ####
## Genotype SNP data ####
# this has been pre-filtered for quality and subset to 1000 loci for this example
# data is biallelic (0=aa, 1=Aa, 2=AA)
eidothea_gt <- read.csv('data/eidothea_snp_data.tsv', row.names = 1)
elaeocarpus_gt <- read.csv('data/elaeocarpus_snp_data.tsv', row.names = 1)
uromyrtus_gt <- read.csv('data/uromyrtus_snp_data.tsv', row.names = 1)

## Other data ####
veg_sf <- read_sf('data/vegetation_map.gpkg') # vegetation map polygons for risk assessment
fire_risk_df <- read.xlsx('data/risk_data.xlsx') # risk values for vegetation types (burn area ratios)
meta_sf <- read_sf('data/spatial_data.gpkg') # spatial metadata of target individuals in genotype data
# spatial data is in GDA94, and locations have been modified to mask the location of the endangered species

## Define input variables #### 
## Define weights for G-value components ####
# here all factors have equal weight, but relative importance could be changed on a case by case basis
weights <- c(
  "dist_within_clusters" = 1,
  "dist_between_clusters" = 1,
  "ta" = 1, # proportion of total alleles
  "n_normalized" = 1, # proportion of maximum n per cluster
  "ho_normalized" = 1 # proportion of maximum Ho per cluster
)

## Set parameters for cluster creation ####
min_n <- 3  # Minimum cluster size
max_perimeter <- 700  # Maximum perimeter for clusters
eps_value <- 100  # DBSCAN epsilon (maximum distance for a point to be from other points to be included in a cluster)
expand_m <- 10  # Buffer distance for convex hulls
# because this data is in GDA94, one unit is equivalent to 1 meter

## Define color schemes for plots ####
vegetation_colors <- c(
  `Rainforests (Combined)` = "#3288BD", 
  `Rainforests (Pyrophyte 0-<1%)` = "#3288BD", 
  `Rainforests (Pyrophyte 1%-10%)` = "#66C2A5", 
  `Rainforests (Pyrophyte 11%-30%)` = "#ABDDA4", 
  `Wet Sclerophyll Forests (Shrubby sub-formation)` = "#FFFFBF", 
  `Wet Sclerophyll Forests (Grassy sub-formation)` = "#FEE08B", 
  `Dry Sclerophyll Forests (Shrubby sub-formation)` = "#F46D43", 
  Heathlands = "#D53E4F",
  `Grassy Woodlands` = 'grey90'
)

gvalue_bin_cols <- c(">1" = "steelblue2", "0.9—1" = "green", "0.8—0.9" = "orange", '<0.8'='red')

# Set the default theme to theme_bw and remove x and y axis text and titles
theme_set(theme_bw())  # Set default theme to theme_bw
theme_update(
  axis.title = element_blank(),       # Remove x and y axis titles
  axis.text = element_blank(),         # Remove x and y axis text
  theme(legend.key.size = unit(0.5, 'lines'), legend.key.height = unit(0.5, 'lines'))
)

# Initial plot of individual distrubutions ####
# Plot vegetation map with points
initial_map_plot <- ggplot() +
  geom_sf(data = veg_sf, aes(fill = Vegetation), color = 'transparent', alpha = 0.5) +
  scale_fill_manual(values = vegetation_colors) +
  geom_sf(data = meta_sf, aes(color = species, shape = species))

initial_map_plot

# Spatial clustering method: defining in situ management sites  #################################################################
# Step 1: DBSCAN clusters points based on density.
# Step 2: Large DBSCAN clusters are refined using K-means clustering if their perimeter exceeds max_perimeter.

## 1. Extract coordinates and calculate distance matrix
coords_df <- st_coordinates(meta_sf) # Extract coordinates from spatial data
pairwise_geo_dist <- dist(coords_df, method = 'maximum') # Calculate pairwise geographic distances

## 2. Plot k-NN distances for clustering
kNNdistplot(pairwise_geo_dist, k = 6)
abline(h = 125, lty = 2)

## 3. Perform DBSCAN clustering
# DBSCAN identifies density-based clusters using the distance matrix
db <- dbscan::dbscan(pairwise_geo_dist, eps = eps_value, minPts = 3)
cluster_df <- as.data.frame(coords_df)
cluster_df$cluster <- as.factor(db$cluster)
cluster_df$cluster[cluster_df$cluster == "0"] <- NA # Assign NA to noise points

## 4. Filter small clusters and renumber them
# Remove clusters smaller than min_n and renumber remaining clusters
cluster_df <- cluster_df %>%
  mutate(cluster = ifelse(cluster %in% names(which(table(cluster) < min_n)), NA, cluster)) %>%
  mutate(cluster = as.numeric(as.factor(cluster)))

# Add species and sample metadata
cluster_df$species <- meta_sf$species
cluster_df$sample <- meta_sf$sample

## 5. Compute convex hulls for visualization
# Generate convex hulls for each cluster for visualization purposes
hull_sf <- cluster_df[!is.na(cluster_df$cluster), ] %>%
  group_by(cluster) %>%
  slice(chull(X, Y)) %>%
  st_as_sf(coords = c("X", "Y")) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_convex_hull() %>%
  filter(!st_geometry_type(geometry) %in% c("POINT")) %>%
  st_buffer(dist = expand_m) %>%
  mutate(perimeter = st_perimeter(geometry), area = st_area(geometry))

## 6. Save initial plot with convex hulls
plot_list <- list()
plot_list[[length(plot_list) + 1]] <- ggplot() +
  geom_sf(data = hull_sf, fill = 'grey', alpha = 0.3, color = 'black') +
  geom_point(data = cluster_df, aes(X, Y, shape = species, color = species), size = 0.5) +
  scale_color_brewer(palette = "Set1")

## 7. Identify large clusters and refine with K-means
# Identify clusters with a perimeter exceeding max_perimeter
large_clusters <- hull_sf$cluster[hull_sf$perimeter >= max_perimeter] %>% as.vector()
new_cluster_hull_sf <- hull_sf

## 8. Loop through large clusters and apply K-means clustering
while(length(large_clusters) > 0) {
  large_cluster_df <- cluster_df[cluster_df$cluster %in% large_clusters, ]
  new_cluster_df <- data.frame()
  
  # Apply K-means to split each large cluster
  for(original_big_cluster in large_clusters) {
    subset_df <- large_cluster_df[large_cluster_df$cluster == original_big_cluster, ]
    p <- new_cluster_hull_sf$perimeter[new_cluster_hull_sf$cluster == original_big_cluster]
    start_i <- ceiling(p / max_perimeter)
    kmeans_result <- kmeans(subset_df[, 1:2], centers = start_i)
    subset_df$cluster <- paste(subset_df$cluster, kmeans_result$cluster, sep = '.')
    new_cluster_df <- rbind(new_cluster_df, subset_df)
  }
  
  # Update cluster data and recompute convex hulls
  cluster_df <- rbind(cluster_df[!cluster_df$cluster %in% large_clusters, ], new_cluster_df)
  new_cluster_hull_sf <- cluster_df[!is.na(cluster_df$cluster), ] %>%
    group_by(cluster) %>%
    slice(chull(X, Y)) %>%
    st_as_sf(coords = c("X", "Y")) %>%
    summarise(geometry = st_combine(geometry)) %>%
    st_convex_hull() %>%
    filter(!st_geometry_type(geometry) %in% c("POINT")) %>%
    st_buffer(dist = expand_m) %>%
    mutate(perimeter = st_perimeter(geometry), area = st_area(geometry))
  
  plot_list[[length(plot_list) + 1]] <- ggplot() +
    geom_sf(data = new_cluster_hull_sf, fill = 'grey', alpha = 0.3, color = 'black') +
    scale_color_brewer(palette = "Set1") +
    geom_point(data = cluster_df, aes(X, Y, shape = species, color = species), size = 0.5)
  
  large_clusters <- new_cluster_hull_sf$cluster[new_cluster_hull_sf$perimeter >= max_perimeter] %>% as.vector()
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
  plot_list[[1]] + ggtitle('DBSCAN clusters'),
  plot_list[[2]] + ggtitle('K-means clusters round 1'),
  plot_list[[3]] + ggtitle('K-means clusters round 2'),
  common.legend = TRUE, legend = 'bottom',
  ncol = 3, labels = 'auto'
)

clustering_process_plots

# G-value calculation: genetic value for in situ clusters #################################################################
# Calculates G-value variables for each cluster of the target species using distances, heterozygosity, sample size, and allele counts
## 1. Prepare input and output lists for loop
species_names <- c("Eidothea hardeniana", "Elaeocarpus sedentarius", "Uromyrtus australis")
gt_list <- list(eidothea_gt, elaeocarpus_gt, uromyrtus_gt)
names(gt_list) <- species_names
gvalue_list <- list()

## 2. Loop through calculations for each species
for(species in species_names) {
  
  ### 2.1 Subset data to target species
  subset_final_cluster_df <- final_cluster_df[final_cluster_df$species == species,]
  subset_clusters <- unique(subset_final_cluster_df$cluster)
  
  ### 2.2 Get genetic data for the species
  gt0 <- gt_list[[species]]
  gt <- gt0[subset_final_cluster_df$sample,]
  
  ### 2.3 Calculate heterozygosity for each individual
  subset_final_cluster_df$ho <- as.numeric(rowSums((gt == 1), na.rm = TRUE) / rowSums(!is.na(gt)))
  
  ### 2.4 Get mean Ho and number of individuals (n) per cluster
  cluster_summary_df <- subset_final_cluster_df %>%
    group_by(cluster) %>%
    summarise(mean_ho = mean(ho), n = n())
  
  ### 2.5 Calculate total alleles
  # Compute total common alleles (maf > 0.05) per cluster
  total_alleles_df <- calculate_total_alleles(gt, subset_final_cluster_df$cluster, 0.05)
  total_alleles_df$ta <- as.numeric(total_alleles_df$total_allele_count) / 
    max(as.numeric(total_alleles_df$total_allele_count)) # Normalize total alleles
  
  # Merge allele data with cluster summary
  cluster_summary_df <- merge(cluster_summary_df,
                              total_alleles_df[1:(nrow(total_alleles_df) - 1),
                                               c('population', 'total_allele_count', 'ta')], by.x = 'cluster', by.y = 'population')
  
  ### 2.6 Compute Euclidean distances
  pairwise_euclidean_dist <- as.matrix(dist(gt, diag = TRUE)) 
  pairwise_euclidean_dist_normalized <- pairwise_euclidean_dist / max(pairwise_euclidean_dist, na.rm = TRUE) # Normalize distances
  
  #### 2.6.1 Within clusters
  dist_within_clusters <- sapply(subset_clusters, function(cluster) {
    samples <- subset_final_cluster_df$sample[subset_final_cluster_df$cluster == cluster]
    sub_matrix <- pairwise_euclidean_dist_normalized[samples, samples, drop = FALSE]
    mean(sub_matrix[upper.tri(sub_matrix)])
  })
  dist_within_clusters[is.na(dist_within_clusters)] <- 0
  
  #### 2.6.2 Between clusters
  dist_between_clusters <- numeric(length(subset_clusters))
  names(dist_between_clusters) <- subset_clusters
  
  for (i in seq_along(subset_clusters)) {
    samples_within <- subset_final_cluster_df$sample[subset_final_cluster_df$cluster == subset_clusters[i]]
    samples_other <- subset_final_cluster_df$sample[subset_final_cluster_df$cluster != subset_clusters[i]]
    dist_between_clusters[subset_clusters[i]] <- mean(pairwise_euclidean_dist_normalized[samples_within, samples_other, drop = FALSE])
  }
  
  # Add diversity metrics to cluster summary
  cluster_summary_df$dist_between_clusters <- dist_between_clusters[cluster_summary_df$cluster]
  cluster_summary_df$dist_within_clusters <- dist_within_clusters[cluster_summary_df$cluster]
  cluster_summary_df$n_normalized <- cluster_summary_df$n / max(cluster_summary_df$n)
  cluster_summary_df$ho_normalized <- cluster_summary_df$mean_ho / max(cluster_summary_df$mean_ho)
  
  final_cluster_summary_df <- cluster_summary_df[, c('cluster', 'ta',
                                                     'dist_between_clusters', 'dist_within_clusters', 
                                                     'n_normalized', 'ho_normalized')]
  
  ### 2.7 Calculate G-values
  # Calculate raw and normalized G-values
  final_cluster_summary_df$gvalue_raw <- Reduce(`+`, lapply(names(weights), 
                                                            function(col) final_cluster_summary_df[[col]] * weights[col]))
  final_cluster_summary_df$gvalue_final <- final_cluster_summary_df$gvalue_raw / max(final_cluster_summary_df$gvalue_raw) 
  
  final_cluster_summary_df$species <- species
  
  # Save output for the species
  gvalue_list[[species]] <- final_cluster_summary_df
}

## 3. Finalize diversity dataframe
# Combine G-value data for all species
gvalues_df <- bind_rows(gvalue_list)

# Risk calculation ############################################################################################################
# here we calculate the relative risk of each cluster based on what proportion is made of different vegetation types
# the burn area ratio (calculation not shown here) is used to calculate the final risk score


## 1. Compute spatial overlap between clusters and vegetation types ####
# overlap final hulls with vegetation map
intersect_sf <- st_intersection(final_hull_sf, veg_sf)

# Add area of overlap for intersected features
intersect_sf$area_overlap <- st_area(intersect_sf)

# Retain relevant columns and drop geometry
intersect_sf <- intersect_sf %>%
  dplyr::select(cluster, perimeter, area, Vegetation, area_overlap)

intersect_sf$geometry <- NULL

## 2. Merge intersected data with final_hull_sf to retain attributes ####
intersect_sf <- merge(final_hull_sf, intersect_sf, 
                         by = c('cluster', 'area'), all = TRUE) %>%
  units::drop_units()  # Drop units

## 3. Calculate proportion of overlap for each feature ####
intersect_sf$proportion_overlap <- intersect_sf$area_overlap / intersect_sf$area
intersect_sf$geometry <- NULL  # Drop geometry column

## 4. Summarize total proportion overlap by cluster and Vegetation ####
cluster_veg_type_df <- intersect_sf %>%
  as.data.frame() %>%
  dplyr::select(cluster, Vegetation, proportion_overlap) %>%
  group_by(cluster, Vegetation) %>%
  summarize(proportion_overlap = sum(proportion_overlap, na.rm = TRUE)) %>%
  pivot_wider(names_from = Vegetation, values_from = proportion_overlap)

# Remove unnecessary 'NA' column
cluster_veg_type_df$`NA` <- NULL
numeric_columns <- colnames(cluster_veg_type_df)[2:ncol(cluster_veg_type_df)]
# Initialize a risk column
cluster_veg_type_df$risk <- NA

## 5. Calculate risk for each cluster ####
for (i in 1:nrow(cluster_veg_type_df)) {
  row_values <- cluster_veg_type_df[i, numeric_columns]
  row_values[is.na(row_values)] <- 0  # Replace NA with 0
  total_sum <- 0
  
  # Compute risk based on overlap and burn area ratio
  for (col in numeric_columns) {
    if (col %in% fire_risk_df$Vegetation) {
      burn_area_ratio_total <- fire_risk_df$Total[fire_risk_df$Vegetation == col]
      total_sum <- total_sum + (row_values[[col]] * burn_area_ratio_total)
    }
  }
  
  # Assign risk to the current cluster
  cluster_veg_type_df$risk[i] <- total_sum
}

## 6. Merge risk values with gvalues_df ####
gvalues_df <- merge(gvalues_df, 
                          cluster_veg_type_df[, c('cluster', 'risk')], 
                          by = 'cluster')

gvalues_df <- gvalues_df %>%
  group_by(cluster) %>%
  mutate(multispecies = ifelse(n() >= 2, 'y', 'n')) %>%
  ungroup()

gvalues_df <- gvalues_df %>%
  mutate(gvalue_bin = factor(case_when(
    gvalue_final > 1 ~ ">1",
    gvalue_final >= 0.9 & gvalue_final <= 1 ~ "0.9—1",
    gvalue_final >= 0.8 & gvalue_final < 0.9 ~ "0.8—0.9",  
    gvalue_final < 0.8 ~ "<0.8"),
    levels = c("<0.8", "0.8—0.9", "0.9—1", ">1")))

gvalues_df$gvalue_bin <- factor(gvalues_df$gvalue_bin, levels = c(">1", "0.9—1", "0.8—0.9", "<0.8"))


## 7. calculate additive G-values across species ####
gvalues_combined_df <- gvalues_df %>%
  group_by(cluster) %>%
  reframe(
    multispecies = multispecies,
    species = species,
    gvalue_final = sum(gvalue_final),
    risk = risk
  )

gvalues_combined_df <- gvalues_combined_df%>%
  mutate(gvalue_bin = factor(case_when(
    gvalue_final > 1 ~ ">1",
    gvalue_final >= 0.9 & gvalue_final <= 1 ~ "0.9—1",
    gvalue_final >= 0.8 & gvalue_final < 0.9 ~ "0.8—0.9",  
    gvalue_final < 0.8 ~ "<0.8"),
    levels = c("<0.8", "0.8—0.9", "0.9—1", ">1")))

gvalues_combined_df$gvalue_bin <- factor(gvalues_combined_df$gvalue_bin, levels = c(">1", "0.9—1", "0.8—0.9", "<0.8"))


# Final plots  ############################################################################################################
# Faceted plot for risk vs G-value
risk_gvalue_faceted_plot <- ggplot(gvalues_df, aes(y = risk, x = gvalue_final, color = gvalue_bin, shape = multispecies, size = multispecies)) +
  geom_hline(yintercept = 1, colour = 'red', linetype = 'dotted') +
  geom_point() +
  facet_grid(~species) +
  scale_size_manual(values = c(1, 2)) +
  scale_color_manual(values = gvalue_bin_cols) +
  labs(x = "G-value", y = "Fire risk", color = "G-value bin", shape = "Multi-species site", size = "Multi-species site") +
  theme_bw() + theme(strip.text = element_text(face = "italic"))

# Combined plot for risk vs G-value
risk_gvalue_combined_plot <- ggplot(gvalues_combined_df, aes(y = risk, x = gvalue_final, color = gvalue_bin, shape = multispecies, size = multispecies)) +
  geom_hline(yintercept = 1, colour = 'red', linetype = 'dotted') +
  geom_point() +
  scale_size_manual(values = c(1, 2)) +
  scale_color_manual(values = gvalue_bin_cols) +
  labs(x = "G-value", y = "Fire risk", color = "G-value bin", shape = "Multi-species site", size = "Multi-species site") +
  theme_bw() + theme(legend.position = 'bottom')



combined_gvalue_risk_plots <- ggarrange(risk_gvalue_faceted_plot, risk_gvalue_combined_plot,
                                               common.legend=TRUE, legend.grob=get_legend(risk_gvalue_combined_plot),
                                               legend='bottom', 
                                               labels='auto', widths=c(1,0.6))


## Export plots ####

initial_map_plot
ggsave('outputs/example1_map.png', initial_map_plot, width=20, height=10, units = 'cm', bg='white')

clustering_process_plots
ggsave('outputs/example1_clustering_process_plots.png', clustering_process_plots, width=20, height=10, units = 'cm', bg='white')

combined_gvalue_risk_plots
ggsave('outputs/example1_combined_gvalue_risk_plots.png', combined_gvalue_risk_plots, width=20, height=10, units = 'cm', bg='white')
