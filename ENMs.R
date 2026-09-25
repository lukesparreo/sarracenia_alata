#Code for examining potential niche differences occupied by eastern and western Sarracenia alata populations

#Load necessary packages
library(maxnet)
library(terra)
library(dismo)
library(tidyverse)
library(pROC)
library(geodata)
library(rnaturalearth)
library(rnaturalearthdata)
library(ENMeval)
library(usdm)
library(readr)
library(sf)
library(tidyterra)
library(ggplot2)
library(patchwork)
library(ecospat)
library(ade4)

#1 Load occurrences from GBIF

occ_west <- read_tsv("s_alata_west.csv", show_col_types = FALSE,
                     col_types = cols(eventDate = col_character()))
occ_east <- read_tsv("s_alata_east.csv", show_col_types = FALSE,
                     col_types = cols(eventDate = col_character()))

pres_west <- occ_west %>%
  mutate(scientific_name = "Sarracenia alata") %>%
  rename(longitude = decimalLongitude, latitude = decimalLatitude)

pres_east <- occ_east %>%
  mutate(scientific_name = "Sarracenia alata") %>%
  rename(longitude = decimalLongitude, latitude = decimalLatitude)

occ_vect_west <- vect(pres_west, geom = c("longitude", "latitude"), crs = "EPSG:4326")
occ_vect_east <- vect(pres_east, geom = c("longitude", "latitude"), crs = "EPSG:4326")

#2 Preparing variable rasters, cropping study region, generating background points on cropped regions

#Load WorldClim rasters — tiles covering each study region only
hull_vect_west   <- convHull(occ_vect_west)
hull_ea_west     <- project(hull_vect_west, "EPSG:6933")
hull_buf_ea_west <- buffer(hull_ea_west, width = 111320)
hull_buf_west    <- project(hull_buf_ea_west, "EPSG:4326")

hull_vect_east   <- convHull(occ_vect_east)
hull_ea_east     <- project(hull_vect_east, "EPSG:6933")
hull_buf_ea_east <- buffer(hull_ea_east, width = 111320)
hull_buf_east    <- project(hull_buf_ea_east, "EPSG:4326")

#Get a representative point inside each buffered study region to pick the right tile
centroid_west <- geom(centroids(hull_buf_west))
centroid_east <- geom(centroids(hull_buf_east))

bio_tile_west1 <- worldclim_tile(
  var  = "bio",
  lon  = centroid_west[, "x"],
  lat  = 35,    
  path = "climate"
)

bio_tile_west2 <- worldclim_tile(
  var  = "bio",
  lon  = centroid_west[, "x"],
  lat  = 25,     
  path = "climate"
)

bio_tile_east1 <- worldclim_tile(
  var  = "bio",
  lon  = -91,
  lat  = centroid_east[, "y"],
  path = "climate"
)

bio_tile_east2 <- worldclim_tile(
  var  = "bio",
  lon  = -89,
  lat  = centroid_east[, "y"],
  path = "climate"
)

bio_tile_east <- mosaic(bio_tile_east1, bio_tile_east2)
bio_tile_west <- mosaic(bio_tile_west1, bio_tile_west2)

elev_tile_west1 <- worldclim_tile(
  var  = "elev",
  lon  = centroid_west[, "x"],
  lat  = 35,
  path = "climate"
)

elev_tile_west2 <- worldclim_tile(
  var  = "elev",
  lon  = centroid_west[, "x"],
  lat  = 25,
  path = "climate"
)

elev_tile_east1 <- worldclim_tile(
  var  = "elev",
  lon  = -91,
  lat  = centroid_east[, "y"],
  path = "climate"
)

elev_tile_east2 <- worldclim_tile(
  var  = "elev",
  lon  = -89,
  lat  = centroid_east[, "y"],
  path = "climate"
)

elev_tile_west <- mosaic(elev_tile_west1, elev_tile_west2)
elev_tile_east <- mosaic(elev_tile_east1, elev_tile_east2)

#Country outlines for plotting
countries  <- ne_countries(scale = "medium", returnclass = "sf")
study_area <- countries %>% filter(admin %in% c("United States of America", "Canada"))
study_vect <- vect(study_area)

#Create convex hull buffer around points, crop study region

#West
bio_crop_west <- crop(bio_tile_west, hull_buf_west)
bio_crop_west <- mask(bio_crop_west, hull_buf_west)

elev_crop_west <- crop(elev_tile_west, hull_buf_west)
elev_crop_west <- mask(elev_crop_west, hull_buf_west)
elev_crop_west <- resample(elev_crop_west, bio_crop_west, method = "bilinear")

elev_focal_west <- focal(elev_crop_west, w = 5, fun = "mean", na.rm = TRUE)
tpi_broad_west  <- elev_crop_west - elev_focal_west
names(tpi_broad_west) <- "TPI_broad"

bio_crop_west <- c(bio_crop_west, tpi_broad_west)

#East
bio_crop_east <- crop(bio_tile_east, hull_buf_east)
bio_crop_east <- mask(bio_crop_east, hull_buf_east)

elev_crop_east <- crop(elev_tile_east, hull_buf_east)
elev_crop_east <- mask(elev_crop_east, hull_buf_east)
elev_crop_east <- resample(elev_crop_east, bio_crop_east, method = "bilinear")

elev_focal_east <- focal(elev_crop_east, w = 5, fun = "mean", na.rm = TRUE)
tpi_broad_east  <- elev_crop_east - elev_focal_east
names(tpi_broad_east) <- "TPI_broad"

bio_crop_east <- c(bio_crop_east, tpi_broad_east)

#Generate background points — west & east
bg_xy_west <- spatSample(bio_crop_west, size = 10000, method = "random",
                         na.rm = TRUE, xy = TRUE)
bg_west <- data.frame(
  scientific_name = "background",
  longitude       = bg_xy_west$x,
  latitude        = bg_xy_west$y
)

bg_xy_east <- spatSample(bio_crop_east, size = 10000, method = "random",
                         na.rm = TRUE, xy = TRUE)
bg_east <- data.frame(
  scientific_name = "background",
  longitude       = bg_xy_east$x,
  latitude        = bg_xy_east$y
)

#3 Extract environmental values

#West grid
presvals_west <- terra::extract(bio_crop_west, pres_west[, c("longitude", "latitude")])
presvals_west <- presvals_west[, -1]
presvals_west <- presvals_west[complete.cases(presvals_west), ]

bgvals_west <- terra::extract(bio_crop_west, bg_west[, c("longitude", "latitude")])
bgvals_west <- bgvals_west[, -1]
bgvals_west <- bgvals_west[complete.cases(bgvals_west), ]

env_all_west <- rbind(presvals_west, bgvals_west)
p_all_west   <- c(rep(1, nrow(presvals_west)), rep(0, nrow(bgvals_west)))

#East grid
presvals_east <- terra::extract(bio_crop_east, pres_east[, c("longitude", "latitude")])
presvals_east <- presvals_east[, -1]
presvals_east <- presvals_east[complete.cases(presvals_east), ]

bgvals_east <- terra::extract(bio_crop_east, bg_east[, c("longitude", "latitude")])
bgvals_east <- bgvals_east[, -1]
bgvals_east <- bgvals_east[complete.cases(bgvals_east), ]

env_all_east <- rbind(presvals_east, bgvals_east)
p_all_east   <- c(rep(1, nrow(presvals_east)), rep(0, nrow(bgvals_east)))

#4a ENMEval west
enmeval_west <- ENMevaluate(
  occs = pres_west[, c("longitude", "latitude")],
  envs       = bio_crop_west,
  bg         = bg_west[, c("longitude", "latitude")],
  algorithm  = "maxnet",
  partitions = "block",
  tune.args  = list(
    fc = c("L", "LQ", "LQH", "LQHP"),
    rm = c(1, 2, 4)
  ),
  parallel = FALSE
)

eval_tbl_west <- enmeval_west@results

print(eval_tbl_west[, c(
  "fc",
  "rm",
  "AICc",
  "auc.val.avg",
  "or.10p.avg"
)])

#Select best east model based on AICc
idx_west <- which.min(eval_tbl_west$AICc)
best_west <- eval_tbl_west[idx_west, ]
mod_west <- enmeval_west@models[[idx_west]]
cat(
  "West best model — FC:",
  best_west$fc,
  "| rm:",
  best_west$rm,
  "\n"
)

#4b ENMEval east
enmeval_east <- ENMevaluate(
  occs = pres_east[, c("longitude", "latitude")],
  envs       = bio_crop_east,
  bg         = bg_east[, c("longitude", "latitude")],
  algorithm  = "maxnet",
  partitions = "block",
  tune.args  = list(
    fc = c("L", "LQ", "LQH", "LQHP"),
    rm = c(1, 2, 4)
  ),
  parallel = FALSE
)

eval_tbl_east <- enmeval_east@results

print(eval_tbl_east[, c(
  "fc",
  "rm",
  "AICc",
  "auc.val.avg",
  "or.10p.avg"
)])

#Select best east model based on AICc
idx_east <- which.min(eval_tbl_east$AICc)
best_east <- eval_tbl_east[idx_east, ]
mod_east <- enmeval_east@models[[idx_east]]
cat(
  "East best model — FC:",
  best_east$fc,
  "| rm:",
  best_east$rm,
  "\n"
)

#5 view response curves and predictions for the two best models

#West prediction
par(mfrow = c(1,1))
pred_west <- terra::predict(
  bio_crop_west,
  mod_west,
  type = "cloglog",
  na.rm = TRUE
)

plot(
  pred_west,
  main = paste0(
    "West suitability (FC:",
    best_west$fc,
    ", rm:",
    best_west$rm,
    ")"
  )
)

points(
  occ_vect_west,
  pch = 16,
  cex = 0.4,
  col = "white"
)

#West curves
bio_vars_west <- names(bio_crop_west)

means_all_west <- colMeans(
  env_all_west,
  na.rm = TRUE
)

par(mfrow = c(4,5))

for (v in bio_vars_west) {
  var_seq <- seq(
    min(env_all_west[[v]], na.rm = TRUE),
    max(env_all_west[[v]], na.rm = TRUE),
    length.out = 100
  )
  
  newdata <- as.data.frame(
    matrix(
      rep(means_all_west, 100),
      nrow = 100,
      byrow = TRUE
    )
  )
  
  colnames(newdata) <- names(means_all_west)
  
  newdata[[v]] <- var_seq
  
  pred <- predict(
    mod_west,
    newdata,
    type = "cloglog"
  )
  
  plot(
    var_seq,
    pred,
    type = "l",
    lwd = 2,
    col = "steelblue",
    ylim = c(0,1),
    main = v,
    xlab = "",
    ylab = "Suitability"
  )
}

#East prediction
par(mfrow = c(1,1))
pred_east <- terra::predict(
  bio_crop_east,
  mod_east,
  type = "cloglog",
  na.rm = TRUE
)

plot(
  pred_east,
  main = paste0(
    "East suitability (FC:",
    best_east$fc,
    ", rm:",
    best_east$rm,
    ")"
  )
)

points(
  occ_vect_east,
  pch = 16,
  cex = 0.4,
  col = "white"
)

#East curves
bio_vars_east <- names(bio_crop_east)

means_all_east <- colMeans(
  env_all_east,
  na.rm = TRUE
)

par(mfrow = c(4,5))

for (v in bio_vars_east) {
  
  var_seq <- seq(
    min(env_all_east[[v]], na.rm = TRUE),
    max(env_all_east[[v]], na.rm = TRUE),
    length.out = 100
  )
  
  newdata <- as.data.frame(
    matrix(
      rep(means_all_east, 100),
      nrow = 100,
      byrow = TRUE
    )
  )
  
  colnames(newdata) <- names(means_all_east)
  
  newdata[[v]] <- var_seq
  
  pred <- predict(
    mod_east,
    newdata,
    type = "cloglog"
  )
  
  plot(
    var_seq,
    pred,
    type = "l",
    lwd = 2,
    col = "steelblue",
    ylim = c(0,1),
    main = v,
    xlab = "",
    ylab = "Suitability"
  )
}

#Compare predictions in one plot
par(mfrow = c(2,1))
plot(pred_west,
     main = "West predicted suitability")
plot(pred_east,
     main = "East predicted suitability")

#6 Make PCA and compare niche differences

#Combine environmental rasters
bio_combined <- mosaic(
  bio_crop_west,
  bio_crop_east
)

#Extract available environmental space
env_background <- terra::values(
  bio_combined,
  na.rm = TRUE
)

env_background <- as.data.frame(env_background)
env_background <- env_background[complete.cases(env_background), ]

#Run PCA
pca_env <- dudi.pca(
  env_background,
  scannf = FALSE, #will not ask for # of axes
  nf = 2 #number of axes
)

#Variance explained
pca_env$eig / sum(pca_env$eig)

#Extract occurrence environments
west_env <- terra::extract(
  bio_combined,
  occ_vect_west
)

west_env <- west_env[, -1]

east_env <- terra::extract(
  bio_combined,
  occ_vect_east
)

east_env <- east_env[, -1]

west_env <- west_env[complete.cases(west_env), ]
east_env <- east_env[complete.cases(east_env), ]

#Project occurrences into PCA space
west_scores <- suprow(
  pca_env,
  west_env
)$li

east_scores <- suprow(
  pca_env,
  east_env
)$li

#Plot PCA
plot(
  west_scores[,1],
  west_scores[,2],
  pch = 16,
  col = "blue",
  xlab = paste0(
    "PC1 (",
    round(100*pca_env$eig[1]/sum(pca_env$eig),1),
    "%)"
  ),
  ylab = paste0(
    "PC2 (",
    round(100*pca_env$eig[2]/sum(pca_env$eig),1),
    "%)"
  ),
  xlim = range(c(west_scores[,1], east_scores[,1])),
  ylim = range(c(west_scores[,2], east_scores[,2]))
)

points(
  east_scores[,1],
  east_scores[,2],
  pch = 16,
  col = "red"
)

legend(
  "topright",
  legend = c("West", "East"),
  col = c("blue", "red"),
  pch = 16
)

#Calculating niche overlap
#Extract scores from PCA
env_scores <- pca_env$li #where each occurance falls in PCA space

#Create niche grids
z_west <- ecospat.grid.clim.dyn(
  glob = env_scores,
  glob1 = west_scores,
  sp = west_scores,
  R = 100 #resolution of PCA space
)

z_east <- ecospat.grid.clim.dyn(
  glob = env_scores,
  glob1 = east_scores,
  sp = east_scores,
  R = 100 #resolution of PCA space
)

#Niche equivalency test for divergence
equivalency <- ecospat.niche.equivalency.test(
  z_west, z_east,
  rep = 1000,
  overlap.alternative    = "lower", #testing divergence, not conservatism
  expansion.alternative  = "higher",
  stability.alternative  = "lower",
  unfilling.alternative  = "higher"
)

equivalency

#Niche similarity tests (bidirectional)
#West -> East
similarity_west_east <- ecospat.niche.similarity.test(
  z_west,
  z_east,
  rep = 1000
)

similarity_west_east

#East -> West
similarity_east_west <- ecospat.niche.similarity.test(
  z_east,
  z_west,
  rep = 1000
)

similarity_east_west

#Plot West -> East null distribution
hist(similarity_west_east$sim$D,
     breaks = 30,
     main = "Niche similarity test: West to East",
     xlab = "Schoener's D")

abline(v = similarity_west_east$obs$D,
       col = "red",
       lwd = 3)

# Plot East -> West null distribution
hist(similarity_east_west$sim$D,
     breaks = 30,
     main = "Niche similarity test: East to West",
     xlab = "Schoener's D")

abline(v = similarity_east_west$obs$D,
       col = "red",
       lwd = 3)
