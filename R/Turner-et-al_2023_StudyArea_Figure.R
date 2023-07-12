### Load Packages ----
libs <- c('ggplot2', "dplyr", "sf", "data.table",
          "rnaturalearth", "rnaturalearthdata", "rnaturalearthhires",
          "ggspatial", "raster", "cowplot")
lapply(libs, require, character.only = TRUE)

### Load map shapefiles ----
# Read in Canada base map
canada_map <- ne_states(country = 'canada', returnclass = "sf")

# Study areas
GHA26 <- st_read("Input/GHA26/Study_area_Eastern_MB.shp")
RMNP <- st_read("Input/RMNP/rmnp2.shp")

# Reproject to UTM14N coordinates
utm14N <- "+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

GHA26 <- st_transform(GHA26, utm14N)
RMNP <- st_transform(RMNP, utm14N)
canada_map <- st_transform(canada_map, utm14N)

### Generate subset of wolf steps ----
steps <- fread("Input/RMNP/map_pts.csv")

GHA.steps <- steps[pop == "GHA26" & id == "W01" & 
                     step_id_ %in% 2239:2255,]

RMNP.steps <- steps[pop == "RMNP" & id == "W03" & 
                      step_id_ %in% 2460:2479,]

### Mapping Rasters ----
## This part will take a while to run - after you have done it once, you can skip to "START FROM HERE"
mb.lc <- raster("GIS_Layers/Landcover/Landcover_2015_30.tif")
mb.lc <- projectRaster(mb.lc, crs = utm14N)

rmnp.lc <- raster("Input/RMNP/RMNPlandcover_wgs84.tif")
rmnp.lc <- projectRaster(rmnp.lc, crs = utm14N)

# Crop raster as small as you can 
lc.GHA26 <- crop(mb.lc, GHA26)
lc.RMNP <- crop(rmnp.lc, RMNP)

# Turn raster into as spatial pixels data frame and save as a data.frame
GHA26_spdf <- as(lc.GHA26, "SpatialPixelsDataFrame")
GHA26_df <- as.data.frame(GHA26_spdf)
colnames(GHA26_df) <- c("Value", "x", "y")

RMNP_spdf <- as(lc.RMNP, "SpatialPixelsDataFrame")
RMNP_df <- as.data.frame(RMNP_spdf)
colnames(RMNP_df) <- c("Value", "x", "y")

# Merge legend with raster data.frame 
legend <- read.csv("GIS_Layers/Landcover/CoverType_2015_Legend.csv")
RMNP.landcover <- merge(RMNP_df, legend, by = "Value")
GHA.landcover <- merge(GHA26_df, legend, by = "Value")

### Build figure
## Part (a)
# Plot the raster data.frame using geom_tile
RMNP.plot <- ggplot() +  
  geom_tile(data=RMNP.landcover, aes(x=x, y=y, fill= Cover)) + 
  scale_fill_grey(start = 0.5, end = 1, 
                  labels=c('Coniferous', 'Deciduous', 'Mixed', 
                           'Open', 'Urban', 'Wet')) +
  geom_sf(data=RMNP, fill= "grey", color="black", 
          size=0.25, alpha = 0.2) + 
  geom_rect(aes(xmin = 382500, xmax = 391200, ymin = 5632000, ymax = 5638000),
            colour = 'black', fill = NA, size = 2) +
  annotation_scale(location = "br", width_hint = 0.25) +
  theme_bw() +
  theme(legend.title = element_text(size=14),
        legend.text = element_text(size=14))

# Create inset showing subset of steps across habitats
RMNP.inset <- ggplot() +  
  geom_tile(data=RMNP.landcover, aes(x=x, y=y, fill= Cover)) + 
  scale_fill_grey(start = 0.5, end = 1) + 
  geom_point(data = RMNP.steps, aes(x = x1_, y = y1_)) +
  geom_path(data = RMNP.steps, aes(x = x1_, y = y1_)) +
  coord_sf(xlim = c(382500, 391200), ylim = c(5632000, 5638000), expand = TRUE) +
  annotation_scale(location = "br", width_hint = 0.25) +
  theme_void() + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill= NA, size = 1),
        legend.position = "none") 

# Add inset to study map
RMNP.full <- ggdraw(RMNP.plot) + 
  draw_plot(RMNP.inset, width = 0.3, height = 0.3, x = 0.1, y = 0.15)

## Part (b)
# Create map showing study areas within MB
mb.study.areas <- ggplot() +  
  geom_sf(data = canada_map, colour = "black", fill = "white") + 
  geom_sf(data=GHA26, fill= "grey", color="black", size=0.25, 
          alpha = 0.9) +
  geom_sf(data=RMNP, fill= "grey", color="black", size=0.25, 
          alpha = 0.9) +
  annotation_scale(location = "br", width_hint = 0.25) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         pad_x = unit(0.25, "cm"), pad_y = unit(0.25, "cm"),
                         style = north_arrow_fancy_orienteering) +
  xlim(c(150000, 1104200)) + 
  ylim(c(5300449, 6743390)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", size = 1),
        panel.border = element_rect(colour = "black", fill= NA, size = 1)) 

## Part (c)
GHA.plot <- ggplot() +  
  geom_tile(data=GHA.landcover, aes(x=x, y=y, fill= Cover)) + 
  scale_fill_grey(start = 0.5, end = 1, labels=c('Coniferous', 'Deciduous', 'Mixed', 
                                                 'Open', 'Urban', 'Wet')) +
  geom_sf(data=GHA26, fill= "grey", color="black", 
          size=0.25, alpha = 0.2) +
  geom_rect(aes(xmin = 750000, xmax = 768000, ymin = 5600000, ymax = 5619000),
            fill = NA, colour = "black", size = 2) +
  annotation_scale(location = "br", width_hint = 0.25) +
  theme_bw() +
  theme(legend.title = element_text(size=14),
        legend.text = element_text(size=14))

# Create inset showing subset of steps across habitats
GHA.inset <- ggplot() +  
  geom_tile(data=GHA.landcover, aes(x=x, y=y, fill= Cover)) + 
  scale_fill_grey(start = 0.5, end = 1) + 
  geom_point(data = GHA.steps, aes(x = x1_, y = y1_)) +
  geom_path(data = GHA.steps, aes(x = x1_, y = y1_)) +
  coord_sf(xlim = c(750000, 768000), ylim = c(5600000, 5619000), expand = TRUE) + 
  annotation_scale(location = "br", width_hint = 0.25) +
  theme_void() + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill= NA, size = 1),
        legend.position = "none") 

# Add inset to study maps
GHA.full <- ggdraw(GHA.plot) + 
  draw_plot(GHA.inset, width = 0.3, height = 0.3, x = 0.1, y = 0.15)

## Create full map
MB.full <- ggdraw() +
  draw_plot(RMNP.full, x = 0, y = 0.125, width = .4, height = .75) +
  draw_plot(inset, x = 0.4, y = 0.1, width = .2, height = 0.8) +
  draw_label("Manitoba", x = 0.49, y = 0.6, size = 12) +
  draw_label("Ontario", x = 0.55, y = 0.39, size = 12) + 
  draw_label("United States of America", x = 0.48, y = 0.15, size = 12) + 
  draw_plot(GHA.full, x = 0.65, y = 0.05, width = .35, height = 0.9) +
  draw_plot_label(label = c("a)", "b)", "c)"), size = 15,
                  x = c(0, 0.4, 0.6), y = c(1, 1, 1))

ggsave("JWT_2023_manuscript_FigS1_Update.jpeg", plot = MB.full, bg = "white", 
       width = 21, height = 7, units = "in", dpi = "print")


