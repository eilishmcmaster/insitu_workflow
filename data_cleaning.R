library(dplyr)
library(sf)
library(openxlsx)
library(ggplot2)
library(raster)
library(plotly)




veg <- st_read('/Users/eilishmcmaster/Documents/insitu_workflow/data/Nightcap_reserves_vegetation.gpkg') %>%
  st_transform(crs='EPSG:3577')

meta_all <- read.xlsx('/Users/eilishmcmaster/Documents/insitu_workflow/data/spatial_data.xlsx')
meta_sf <- st_as_sf(meta_all, coords = c("long", "lat"), crs = 'EPSG:4326') %>%
  st_transform(crs='EPSG:3577')
# meta_sf <- st_as_sf(meta_all, coords = c("long", "lat"), crs = 'EPSG:3577')

# set subset axes
subset_xlim <- c(2053000, 2055500)
subset_ylim <- c(-3276700, -3273500)

meta_cropped <- meta_sf %>%
  filter(st_coordinates(.)[, 1] >= subset_xlim[1] & 
           st_coordinates(.)[, 1] <= subset_xlim[2] &
           st_coordinates(.)[, 2] >= subset_ylim[1] & 
           st_coordinates(.)[, 2] <= subset_ylim[2])

bbox <- st_bbox(c(xmin = subset_xlim[1], xmax = subset_xlim[2],
                  ymin = subset_ylim[1], ymax = subset_ylim[2]),
                crs = st_crs(veg))

# Crop the vegetation data to the bounding box
veg_cropped <- st_crop(veg, bbox)


vegetation_colors <- c(`Rainforests (Combined)` = "#3288BD", 
                       `Rainforests (Pyrophyte 0-<1%)` = "#3288BD", 
                       `Rainforests (Pyrophyte 1%-10%)` = "#66C2A5", 
                       `Rainforests (Pyrophyte 11%-30%)` = "#ABDDA4", 
                       `Wet Sclerophyll Forests (Shrubby sub-formation)` = "#FFFFBF", 
                       `Wet Sclerophyll Forests (Grassy sub-formation)` = "#FEE08B", 
                       `Dry Sclerophyll Forests (Shrubby sub-formation)` = "#F46D43", 
                       Heathlands = "#D53E4F",
                       `Grassy Woodlands`='grey90')


map_plot <- ggplot() +
  geom_sf(data= veg_cropped, 
          aes(fill=Vegetation), 
          color='transparent', alpha=0.5)+
  scale_fill_manual(values=vegetation_colors)+
  geom_sf(data=meta_cropped,
          aes(color=sp, shape=sp),
          inherit.aes = FALSE, size=1)+
  theme_bw()+ 
  theme(legend.key.size = unit(0.5, 'lines'),legend.key.height = unit(0.5, 'lines'), axis.title = element_blank())

map_plot


##
library(sf)
library(dplyr)

# Function to translate geometries
translate_geometry <- function(sf_obj, dx, dy) {
  # Translate geometries by adding the specified offsets
  sf_obj <- sf_obj %>%
    mutate(geom = st_geometry(.) + c(dx, dy))
  return(sf_obj)
}

# Apply to veg_cropped and meta_cropped
veg_translated <- veg_cropped %>% mutate(geom = st_geometry(.) + c(-5845, -23891)) %>% st_sf(crs = 3577)
meta_translated <- meta_cropped %>% mutate(geometry = st_geometry(.) + c(-5845, -23891)) %>% st_sf(crs = 3577)


ggplot() +
  geom_sf(data= veg_translated, 
          aes(fill=Vegetation), 
          color='transparent', alpha=0.5)+
  scale_fill_manual(values=vegetation_colors)+
  geom_sf(data=meta_translated,
          aes(color=sp, shape=sp),
          inherit.aes = FALSE, size=1)+
  theme_bw()+ 
  theme(legend.key.size = unit(0.5, 'lines'),legend.key.height = unit(0.5, 'lines'), axis.title = element_blank())


write_sf(veg_translated, '/Users/eilishmcmaster/Documents/insitu_workflow/data/Vegetation_transformed.gpkg')
write_sf(meta_translated, '/Users/eilishmcmaster/Documents/insitu_workflow/data/spatial_data.gpkg')
