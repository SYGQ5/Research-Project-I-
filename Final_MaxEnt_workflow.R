# ------------------------------------------
# Simple MaxEnt workflow in R using maxnet
# with bias-weighted background sampling
# No deduplication of presences
# ------------------------------------------

library(terra)
library(readr)
library(dplyr)
library(maxnet)

# ------------------------------------------
# 1. File paths
# ------------------------------------------
occ_file   <- file.choose() # Species presence csv file e.g."Gactornis_enarratus_presence.csv"
env_file   <- file.choose() #"madagascar_sdm_predictors_1929_1931_50km.tif"
bias_file  <- file.choose() #"AVE_sampling_bias.tif"

# ------------------------------------------
# 2. Read presence data
# ------------------------------------------
occ <- read_csv(occ_file, show_col_types = FALSE)

# Keep duplicates
occ <- occ %>%
  filter(!is.na(longitude), !is.na(latitude))

# ------------------------------------------
# 3. Read environmental predictors
# ------------------------------------------
env <- rast(env_file)

# ------------------------------------------
# 4. Read and align bias raster
# ------------------------------------------
bias <- rast(bias_file)

# Skip this if bias and env rasters are already matched to same grid size
# Project / resample bias raster to match predictors
if (!same.crs(env, bias)) {
  bias <- project(bias, env[[1]])
}

bias <- resample(bias, env[[1]], method = "bilinear")

# Mask bias raster to predictor extent
bias <- mask(bias, env[[1]])

# Replace NA with 0 for weighted sampling
bias[is.na(bias)] <- 0

# --------------------------------------------
# 5. Extract environmental values at presences
# --------------------------------------------
occ_vect <- vect(occ, geom = c("longitude", "latitude"), crs = "EPSG:4326")

# Reproject presences if predictors are not lon/lat
if (!same.crs(occ_vect, env)) {
  occ_vect <- project(occ_vect, crs(env))
}

pres_env <- extract(env, occ_vect, ID = FALSE)

# Drop rows with missing predictor values
keep_pres <- complete.cases(pres_env)
occ_vect  <- occ_vect[keep_pres, ]
pres_env  <- pres_env[keep_pres, ]

# -----------------------------------------------------------------
# 6. Sample background points weighted by bias raster or randomly
# -----------------------------------------------------------------
set.seed(42)

# Number of available non-NA cells
n_available <- sum(!is.na(values(env[[1]])))

# Use all available cells, or cap at a smaller number if preferred
n_bg <- n_available

# Do you want to run bias-weighted or random background selection?
# Use for selection of background points weighted by bias raster
bg_xy <- spatSample(
  bias,
  size = n_bg,
  method = "weights",
  xy = TRUE,
  values = FALSE,
  na.rm = TRUE
)
bg_xy <- as.data.frame(bg_xy)
bg_vect <- vect(bg_xy, geom = c("x", "y"), crs = crs(bias))
bg_env  <- extract(env, bg_vect, ID = FALSE)
bg_env <- bg_env[complete.cases(bg_env), , drop = FALSE]

# Use for random selection of background points instead of points weighted by bias raster
bg_xy <- spatSample(
  env[[1]],
  size = n_bg,
  method = "random",
  xy = TRUE,
  values = FALSE,
  na.rm = TRUE
)
bg_xy <- as.data.frame(bg_xy)
bg_vect <- vect(bg_xy, geom = c("x", "y"), crs = crs(env))
bg_env  <- extract(env, bg_vect, ID = FALSE)
bg_env <- bg_env[complete.cases(bg_env), , drop = FALSE]

# Drop background rows with missing predictor values
keep_bg <- complete.cases(bg_env)
bg_xy   <- bg_xy[keep_bg, , drop = FALSE]
bg_env  <- bg_env[keep_bg, ]

# ------------------------------------------
# 7. Build modelling table
# ------------------------------------------
pres_df <- as.data.frame(pres_env)
bg_df   <- as.data.frame(bg_env)

pb <- c(rep(1, nrow(pres_df)), rep(0, nrow(bg_df)))
model_df <- bind_rows(pres_df, bg_df)

# ------------------------------------------
# 8. Fit MaxEnt model with maxnet
# ------------------------------------------
# Start simple for this MRes project:
# linear + quadratic + hinge features
form <- maxnet.formula(
  p = pb,
  data = model_df,
  classes = "lqh"
)

mx <- maxnet(
  p = pb,
  data = model_df,
  f = form,
  regmult = 1
)

# ------------------------------------------
# 9. Predict across study area
# ------------------------------------------
names(env)
names(model_df)

# Do you want predictions for the bias-weighted or random models?
# Use this routine if using bias raster
pred_bias <- terra::predict(
  env,
  mx,
  fun = function(model, data) {
    predict(model, data, type = "cloglog")
  },
  na.rm = TRUE
)
# Change file name according to species 
writeRaster(
  pred_bias,
  "Gactornis_enarratus_prediction_bias.tif",
  overwrite = TRUE
)

# Use this routine if using random background selection
pred_random <- terra::predict(
  env,
  mx,
  fun = function(model, data) {
    predict(model, data, type = "cloglog")
  },
  na.rm = TRUE
)
# Change file name according to species 
writeRaster(
  pred_random,
  "Gactornis_enarratus_prediction_random.tif",
  overwrite = TRUE
)

# Comparison of models using Boyce Index
library(modEvA)

# Presence coordinates for Boyce Index
# Use the projected / aligned presence vector (occ_vect) already used in modelling
pres_xy <- crds(occ_vect)

# Bias-weighted model
boyce_bias <- modEvA::Boyce(
  obs = crds(occ_vect),
  pred = pred_bias,
  n.bins = NA,
  rm.dup.points = TRUE, #change to FALSE if you want to keep duplicate observations
  rm.dup.classes = FALSE,
  plot = FALSE,
  na.rm = TRUE
)

print(boyce_bias)

# Random model
boyce_random <- modEvA::Boyce(
  obs = crds(occ_vect),
  pred = pred_random,
  n.bins = NA,
  rm.dup.points = TRUE, #change to FALSE if you want to keep duplicate observations
  rm.dup.classes = FALSE,
  plot = FALSE,
  na.rm = TRUE
)

print(boyce_random)

# ------------------------------------------
# 10. Quick plots
# ------------------------------------------

# Plot of suitability (bias) - change "main" title dependent on species
plot(pred_bias, main = "Suitability (bias): Gactornis enarratus")
points(crds(occ_vect), pch = 16, cex = 0.5)

# Plot of suitability (random) - change "main" title dependent on species
plot(pred_random, main = "Suitability (random): Gactornis enarratus")
points(crds(occ_vect), pch = 16, cex = 0.5)

# ------------------------------------------
# 11. Comparisons bias vs random
# ------------------------------------------
diff_map <- pred_bias - pred_random
plot(diff_map, main = "Difference in suitability (bias - random)")

cor(
  values(pred_bias, na.rm = TRUE),
  values(pred_random, na.rm = TRUE)
)
# ------------------------------------------
# 12. Final plots
# ------------------------------------------

library(rnaturalearth)

library(sf)

# Creating vector for Madagascar outline

mad_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

mad_sf <- mad_sf[mad_sf$admin == "Madagascar", ]

mad_vect <- vect(mad_sf)

# To create one figure with three plots horizontally, so one row with three columns
# Also to create line spacing for overall shared title and x axis label 
par(mfrow = c(1, 3), oma = c(0.5, 0, 2, 0))

# For plotting suitability (bias)

plot(pred_bias,
     
     ylab = expression("Latitude ("*degree*")"),
     
     main = "(a) Bias", 
     
     cex.main = 1,
     
     col = hcl.colors(100, "viridis"))

plot(mad_vect, add = TRUE, border = "white", lwd = 2)

points(occ_vect,
       
       pch = 21,
       
       bg = "white",
       
       col = "black",
       
       cex = 1.2,
       
       lwd = 0.5)

# For plotting suitability (random)

plot(pred_random,
     
     main = "(b) Random",
     
     cex.main = 1,
     
     col = hcl.colors(100, "viridis"))

plot(mad_vect, add = TRUE, border = "white", lwd = 2)

points(occ_vect,
       
       pch = 21,
       
       bg = "white",
       
       col = "black",
       
       cex = 1.2,
       
       lwd = 0.5)

# For plotting the Difference plot
diff_map <- pred_bias - pred_random
plot(diff_map,
     main = "(c) Difference \n(bias - random)",
     cex.main = 1
)

cor(
  values(pred_bias, na.rm = TRUE),
  values(pred_random, na.rm = TRUE)
)

plot(mad_vect, add = TRUE, border = "white", lwd = 2)

# Shared Titles
# Adding one shared main title for the three tile plots

# Change text in "italic()" code for the species being plotted
mtext(
  expression("Suitability: " * italic("Gactornis enarratus")),
  side = 3, outer = TRUE, adj = 0.5, cex = 1, line = -3
)
# Adding shared x axis label - longitude and degrees unit
mtext(
  expression("Longitude ("*degree*")"),
  side = 1, outer = TRUE, adj = 0.5, cex = 0.5, line = -5
)

