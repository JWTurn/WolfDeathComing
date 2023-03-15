# Load Packages
libs <- c('ggplot2', "dplyr", "sf", "rgdal", "data.table",
          "rnaturalearth", "rnaturalearthdata", "rnaturalearthhires",
          "ggspatial", "raster", "dplyr", "cowplot")
lapply(libs, require, character.only = TRUE)

# Read in MB base map
canada_map <- ne_states(country = 'canada', returnclass = "sf")
map_mb <- canada_map[5,] 

### Mapping Shapefiles without aesthetic ----
GHA26 <- readOGR("GIS_Layers/Simple_Polygon/Study_area_Eastern_MB.shp")
RMNP <- readOGR("Input/RMNP/RMNP_StudyArea/rmnp2.shp")

## Reproject to lat/long coordinates
utm14N <- "+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
GHA26 <- spTransform(GHA26, utm14N)
RMNP <- spTransform(RMNP, utm14N)
map_mb <- st_transform(map_mb, utm14N)

# Collect subset of wolf steps
steps <- fread("Input/RMNP/map_pts.csv")
GHA.steps <- steps[id == "W01" & step_id_ %in% 2242:2246,]
#GHA.steps <- steps[id == "W01" & step_id_ %in% 2144:2248,]

#RMNP.steps <- steps[pop == "RMNP" & id == "W03" & 
#                      step_id_ %in% 2490:2494,]

RMNP.steps <- steps[pop == "RMNP" & id == "W03" & 
                      step_id_ %in% 2492:2496,]

## Because we have not transformed our shapefile into a `sf`, we can use geom_polygon
# Can be a bit slower this way if using a more complicated or larger polygon
ggplot() +  
  geom_polygon(data=GHA26, mapping = aes(x=long, y=lat, group=group), 
               fill= "forestgreen", color="black", size=0.25, alpha = 0.2) +  
  geom_polygon(data=RMNP, mapping = aes(x=long, y=lat, group=group), 
               fill= "forestgreen", color="black", size=0.25, alpha = 0.2) 

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

# Merge legend with raster data.frame and save for the future
#reclass <- read.csv("Input/RMNP/rcl_fine.csv")
legend <- read.csv("GIS_Layers/Landcover/CoverType_2015_Legend.csv")
RMNP.dat <- merge(RMNP_df, legend, by = "Value")
fwrite(RMNP.dat, "RMNP_Landcover_SPDF_easting-northing.Rds")

# START FROM HERE
GHA.landcover <- fread("Input/MB_Landcover_SPDF_easting-northing.Rds")
RMNP.landcover <- fread("Input/RMNP/RMNP_Landcover_SPDF_easting-northing.Rds")

# Plot the raster data.frame using geom_tile
GHA.plot <- ggplot() +  
  geom_tile(data=GHA.landcover, aes(x=x, y=y, fill= Cover)) + 
  scale_fill_grey(start = 0.5, end = 1, labels=c('Coniferous', 'Deciduous', 'Mixed', 
                                                 'Open', 'Urban', 'Wet')) +
  geom_polygon(data=GHA26, mapping = aes(x=long, y=lat, group=group), 
               fill= "grey", color="black", size=0.25, alpha = 0.2) +
  geom_rect(aes(xmin = 755500, xmax = 752100, ymin = 5615000, ymax = 5604100),
            fill = NA, colour = "black", size = 2) +
  annotation_scale(location = "br", width_hint = 0.25) + 
  theme_bw() +
  theme(legend.title = element_text(size=14),
        legend.text = element_text(size=14))

RMNP.plot <- ggplot() +  
  geom_tile(data=RMNP.landcover, aes(x=x, y=y, fill= Cover)) + 
  scale_fill_grey(start = 0.5, end = 1, 
                  labels=c('Coniferous', 'Deciduous', 'Mixed', 
                           'Open', 'Urban', 'Wet')) +
  geom_polygon(data=RMNP, mapping = aes(x=long, y=lat, group=group), 
               fill= "grey", color="black", size=0.25, alpha = 0.2) + 
  geom_rect(aes(xmin = 389900, xmax = 392000, ymin = 5633200, ymax = 5636000),
            colour = 'black', fill = NA, size = 2) +
  annotation_scale(location = "br", width_hint = 0.25) +
  theme_bw() +
  theme(legend.title = element_text(size=14),
        legend.text = element_text(size=14))

# Create first inset showing GHA26 within MB
inset <- ggplot() +  
  geom_sf(data = map_mb, colour = "black", fill = "white") + 
  geom_polygon(data=GHA26, mapping = aes(x=long, y=lat, group=group), 
               fill= "grey", color="black", size=0.25, alpha = 0.9) +
  geom_polygon(data=RMNP, mapping = aes(x=long, y=lat, group=group), 
               fill= "grey", color="black", size=0.25, alpha = 0.9) +
  annotation_scale(location = "br", width_hint = 0.25) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         pad_x = unit(0.25, "cm"), pad_y = unit(0.25, "cm"),
                         style = north_arrow_fancy_orienteering) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", size = 1),
        panel.border = element_rect(colour = "black", fill= NA, size = 1)) 

# Create second inset showing subset of steps across habitats
GHA.inset <- ggplot() +  
  geom_tile(data=GHA.landcover, aes(x=x, y=y, fill= Cover)) + 
  scale_fill_grey(start = 0.5, end = 1) + 
  geom_point(data = GHA.steps, aes(x = x1_, y = y1_)) +
  geom_path(data = GHA.steps, aes(x = x1_, y = y1_)) +
  coord_sf(xlim = c(752000, 755500), ylim = c(5604000, 5615000), expand = TRUE) + 
  annotation_scale(location = "bl", width_hint = 0.25) +
  theme_void() + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill= NA, size = 1),
        legend.position = "none") 

RMNP.inset <- ggplot() +  
  geom_tile(data=RMNP.landcover, aes(x=x, y=y, fill= Cover)) + 
  scale_fill_grey(start = 0.5, end = 1) + 
  geom_point(data = RMNP.steps, aes(x = x1_, y = y1_)) +
  geom_path(data = RMNP.steps, aes(x = x1_, y = y1_)) +
  coord_sf(xlim = c(391900, 389900), ylim = c(5635900, 5633200), expand = TRUE) +
  annotation_scale(location = "br", width_hint = 0.25) +
  theme_void() + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill= NA, size = 1),
        legend.position = "none") 


## Add insets to study maps
GHA.full <- ggdraw(GHA.plot) + 
  draw_plot(GHA.inset, width = 0.3, height = 0.3, x = 0.1, y = 0.15)


RMNP.full <- ggdraw(RMNP.plot) + 
  draw_plot(RMNP.inset, width = 0.3, height = 0.3, x = 0.1, y = 0.15)

## Create full map
MB.full <- ggdraw() +
  draw_plot(RMNP.full, x = 0, y = 0.125, width = .4, height = .75) +
  draw_plot(inset, x = 0.4, y = 0.1, width = .2, height = 0.8) +
  draw_plot(GHA.full, x = 0.65, y = 0.05, width = .35, height = 0.9) +
  draw_plot_label(label = c("a)", "b)", "c)"), size = 15,
                  x = c(0, 0.4, 0.6), y = c(1, 1, 1))

ggsave("JWT_2023_manuscript_FigS1.jpeg", plot = MB.full, bg = "white", 
       width = 21, height = 7, units = "in", dpi = "print")

ggsave("JWT_2023_Manuscript_FigS1.png", plot = MB.full, bg = "white", 
       width = 8, height = 3, units = "in", dpi = "screen")
