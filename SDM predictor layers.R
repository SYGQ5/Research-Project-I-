# =========================================================
# Build SDM predictor layers for Madagascar
# CRU TS 4.09 (1929-1931) + elevation
# =========================================================

library(terra)
library(sf)
library(rnaturalearth)
library(geodata)

# ---------------------------------------------------------
# 1. File paths
# ---------------------------------------------------------
pre_file <- file.choose()
tmp_file <- file.choose()
out_dir  <- "sdm_predictors"

dir.create(out_dir, showWarnings = FALSE)
file.exists(pre_file)
file.exists(tmp_file)

# ---------------------------------------------------------
# 2. Madagascar polygon
# ---------------------------------------------------------
madagascar <- ne_countries(
  scale = "medium",
  country = "Madagascar",
  returnclass = "sf"
)

mad_vect <- vect(madagascar)

# ---------------------------------------------------------
# 3. Read CRU monthly data
# ---------------------------------------------------------
tmp <- rast(tmp_file)
pre <- rast(pre_file)

# Crop and mask first for speed
tmp <- mask(crop(tmp, mad_vect), mad_vect)
pre <- mask(crop(pre, mad_vect), mad_vect)

time(tmp)[1:3]
time(tmp)[(nlyr(tmp)-2):nlyr(tmp)]

# ---------------------------------------------------------
# 4. Subset to 1929-01 through 1931-12 by layer index
# ---------------------------------------------------------

# CRU starts in January 1901
start_year  <- 1901
target_year <- 1929

# First layer for Jan 1929
first_layer <- (target_year - start_year) * 12 + 1

# Last layer for Dec 1931
last_layer  <- first_layer + 36 - 1

tmp_2931 <- tmp[[first_layer:last_layer]]
pre_2931 <- pre[[first_layer:last_layer]]

# Sanity checks
stopifnot(nlyr(tmp_2931) == 36)
stopifnot(nlyr(pre_2931) == 36)

time(tmp_2931)[1:3]
time(tmp_2931)[34:36]

# ---------------------------------------------------------
# 5. Helper indices for years and months
# ---------------------------------------------------------
dates <- time(tmp_2931)
yrs   <- as.integer(format(dates, "%Y"))
mons  <- as.integer(format(dates, "%m"))

# ---------------------------------------------------------
# 6. Predictor 1: annual mean temperature
# mean of all 36 monthly temperature rasters
# ---------------------------------------------------------
bio_temp_mean <- app(tmp_2931, mean, na.rm = TRUE)
names(bio_temp_mean) <- "temp_mean_1929_1931"

plot(bio_temp_mean, main = "Mean temperature 1929-1931")
crs(bio_temp_mean)
ext(bio_temp_mean)
res(bio_temp_mean)

# ---------------------------------------------------------
# 7. Predictor 2: temperature seasonality
# SD of monthly temp within each year, then mean across years
# ---------------------------------------------------------
tmp_sd_years <- lapply(1929:1931, function(y) {
  app(tmp_2931[[yrs == y]], sd, na.rm = TRUE)
})

tmp_sd_stack <- rast(tmp_sd_years)
bio_temp_seasonality <- app(tmp_sd_stack, mean, na.rm = TRUE)
names(bio_temp_seasonality) <- "temp_seasonality_1929_1931"

# ---------------------------------------------------------
# 8. Predictor 3: annual precipitation
# Sum within each year, then mean across years
# ---------------------------------------------------------
pre_sum_years <- lapply(1929:1931, function(y) {
  app(pre_2931[[yrs == y]], sum, na.rm = TRUE)
})

pre_sum_stack <- rast(pre_sum_years)
bio_pre_annual <- app(pre_sum_stack, mean, na.rm = TRUE)
names(bio_pre_annual) <- "pre_annual_1929_1931"

plot(bio_pre_annual, main = "Annual precipitation 1929-1931")

# ---------------------------------------------------------
# 9. Predictor 4: precipitation seasonality
# Coefficient of variation within each year, then mean across years
# CV = sd / mean
# ---------------------------------------------------------
pre_cv_years <- lapply(1929:1931, function(y) {
  pre_y <- pre_2931[[yrs == y]]
  pre_sd <- app(pre_y, sd, na.rm = TRUE)
  pre_mean <- app(pre_y, mean, na.rm = TRUE)
  pre_sd / pre_mean
})

pre_cv_stack <- rast(pre_cv_years)
bio_pre_seasonality <- app(pre_cv_stack, mean, na.rm = TRUE)
names(bio_pre_seasonality) <- "pre_seasonality_1929_1931"

# ---------------------------------------------------------
# 10. Predictor 5: dry-season precipitation
# Use May-October as dry season
# Sum dry-season months within each year, then mean across years
# ---------------------------------------------------------
dry_idx <- mons %in% 5:10

pre_dry_years <- lapply(1929:1931, function(y) {
  idx <- yrs == y & mons %in% 5:10
  app(pre_2931[[idx]], sum, na.rm = TRUE)
})

pre_dry_stack <- rast(pre_dry_years)
bio_pre_dry <- app(pre_dry_stack, mean, na.rm = TRUE)
names(bio_pre_dry) <- "pre_dry_MayOct_1929_1931"

plot(bio_pre_dry, main = "Dry-season precipitation 1929-1931")

# ---------------------------------------------------------
# 11. Elevation
# Download at 30 arc-sec, then aggregate/resample to CRU grid
# ---------------------------------------------------------
elev <- geodata::elevation_30s(
  country = "MDG",
  path = out_dir
)

# Use one CRU-derived layer as template
template <- bio_temp_mean

# Aggregate elevation to roughly match CRU resolution
# CRU is ~0.5 degrees; elevation_30s is ~30 arc-sec = 0.008333 deg
fact_x <- round(res(template)[1] / res(elev)[1])
fact_y <- round(res(template)[2] / res(elev)[2])

elev_agg <- aggregate(
  elev,
  fact = c(fact_x, fact_y),
  fun = mean,
  na.rm = TRUE
)

# Align exactly to template
elev_resampled <- resample(elev_agg, template, method = "bilinear")
elev_resampled <- mask(crop(elev_resampled, mad_vect), mad_vect)
names(elev_resampled) <- "elevation_mean"

# ---------------------------------------------------------
# 12. Combine predictors into one stack
# ---------------------------------------------------------
predictors <- c(
  bio_temp_mean,
  bio_temp_seasonality,
  bio_pre_annual,
  bio_pre_seasonality,
  bio_pre_dry,
  elev_resampled
)

values_df <- as.data.frame(values(predictors, na.rm = TRUE))
cor(values_df, use = "complete.obs")

# Just including predictors with low correlations
predictors <- c(
  bio_temp_mean,
  bio_temp_seasonality,
  bio_pre_dry
)

values_df <- as.data.frame(values(predictors, na.rm = TRUE))
cor(values_df, use = "complete.obs")

# ---------------------------------------------------------------------
# Optional: Resample predictors to match an existing bias raster
# ---------------------------------------------------------------------

# Choose the bias raster
bias_file <- file.choose()
bias_rast <- rast(bias_file)

# Make sure Madagascar polygon is in the same CRS as the bias raster
mad_vect_bias <- project(mad_vect, crs(bias_rast))

# Project/resample predictors directly onto the bias raster grid
predictors_biasgrid <- project(
  predictors,
  bias_rast,
  method = "bilinear"
)

# Mask to Madagascar
predictors_biasgrid <- mask(crop(predictors_biasgrid, mad_vect_bias), mad_vect_bias)

# Keep names
names(predictors_biasgrid) <- names(predictors)

# Replace predictors object
predictors <- predictors_biasgrid

# Quick checks
print(predictors)
res(predictors)
crs(predictors)
ext(predictors)

# ---------------------------------------------------------
# 13. Write outputs
# ---------------------------------------------------------
writeRaster(
  predictors,
  file.path(out_dir, "madagascar_sdm_predictors_1929_1931_50km.tif"),
  overwrite = TRUE
)

# Also write each layer separately if useful
for (i in 1:nlyr(predictors)) {
  writeRaster(
    predictors[[i]],
    file.path(out_dir, paste0(names(predictors)[i], ".tif")),
    overwrite = TRUE
  )
}

# ---------------------------------------------------------
# 14. Quick checks
# ---------------------------------------------------------
print(predictors)
plot(predictors)