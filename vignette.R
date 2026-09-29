
# https://ncss-tech.github.io/rosettaPTF/

remotes::install_github("ncss-tech/rosettaPTF")
library(rosettaPTF)
library(terra)
library(soilDB)

# find python binaries
rosettaPTF::find_python()

# download latest python 3.10.x
# reticulate::install_python(version = "3.10:latest")
# reticulate::virtualenv_create("r-reticulate")


# Install rosetta-soil Python Module
rosettaPTF::install_rosetta()


# =======================================================================================
# vignette examples --------------------------------------------------------------------
# =======================================================================================

# obtain mukey map from SoilWeb Web Coverage Service (800m resolution SSURGO derived)
res <- mukey.wcs(aoi = list(aoi = c(-114.16, 47.65,-114.08, 47.68), crs = 'EPSG:4326'))
# request input data from SDA
varnames <- c("sandtotal_r", "silttotal_r", "claytotal_r", "dbthirdbar_r")
resprop <- get_SDA_property(property = varnames,
                            method = "Dominant Component (numeric)",
                            mukeys = unique(values(res$mukey)))
# keep only those where we have a complete set of 4 parameters (sand, silt, clay, bulk density; model code #3)
soildata <- resprop[complete.cases(resprop), c("mukey", varnames)]

# run Rosetta on the mapunit-level aggregate data
system.time(resrose <- run_rosetta(soildata[,varnames]))
# transfer mukey to result
resprop$mukey <- as.numeric(resprop$mukey)
resrose$mukey <- as.numeric(soildata$mukey)

# merge property (input) and rosetta parameters (output) into RAT
levels(res) <- merge(cats(res)[[1]], resprop, by.x = "ID", by.y = "mukey", all.x = TRUE, sort = FALSE)
levels(res) <- merge(cats(res)[[1]], resrose, by.x = "ID", by.y = "mukey", all.x = TRUE, sort = FALSE)

# convert categories based on mukey to numeric values
res2 <- catalyze(res)

# make a plot of the predicted Ksat
plot(res2, "log10_Ksat_mean")


# test with spat rast ------------------------------------------------------------------
res3 <- rast(list(
  res2[["sandtotal_r"]],
  res2[["silttotal_r"]],
  res2[["claytotal_r"]],
  res2[["dbthirdbar_r"]]
))




# SpatRaster to data.frame interface (one call on all cells)
system.time(test2 <- run_rosetta(res3))
# make a plot of the predicted Ksat (identical to mukey-based results)
plot(test2, "log10_Ksat_mean")




# ========================================================================================
### read in solus rasters ---------------------------------------------------------------

library(tidyverse)
library(sf)
library(terra)
library(tigris)

# get Levy County boundary
levy <- counties(state = "FL", cb = TRUE, resolution = "20m", year = 2025) %>%
  filter(NAME == "Levy") %>%
  st_transform(crs(rast(l[1])))

# select SOLUS layers needed for initial Rosetta testing
(l = list.files("D:\\Projects\\fbrc\\data\\solus", full.names=TRUE, pattern = "\\.tif$"))
(solus_files <- l[str_detect(basename(l), "^(sandtotal|claytotal|dbovendry)_")] )

# read SOLUS layers
solus <- rast(solus_files)

# crop and mask to Levy County
solus_levy <- crop(solus, vect(levy)) %>%
  mask(vect(levy))

plot(solus_levy["sandtotal_15_cm_p"])


# summarize distributions and missing data
global(solus_levy, c("min", "mean", "max", "notNA"), na.rm = TRUE)

# read SOLUS soil-depth rasters
solus_depth <-
  rast(c("D:/Projects/fbrc/data/solus/resdept_all_cm_p.tif",
         "D:/Projects/fbrc/data/solus/anylithicdpt_cm_p.tif")) %>%
  crop(levy) %>%
  mask(levy)

names(solus_depth) <- c("resdept", "lithic_depth")

# align soil-depth rasters with the existing Levy SOLUS grid
solus_depth <- resample(solus_depth, solus_levy[[1]], method = "bilinear")

# extract SOLUS restriction depths by grid cell
depth_dat <-
  as.data.frame(solus_depth, cells = TRUE, na.rm = FALSE) %>%
  select(cell, resdept, lithic_depth)

# compare prediction availability with SOLUS restriction depth
solus_rosetta %>%
  left_join(depth_dat, by = "cell") %>%
  mutate(available = !is.na(sandtotal) & !is.na(claytotal) & !is.na(dbovendry)) %>%
  group_by(depth) %>%
  summarise(
    n = n(),
    pct_available = mean(available) * 100,
    resdept_available = median(resdept[available], na.rm = TRUE),
    resdept_missing = median(resdept[!available], na.rm = TRUE),
    lithic_available = median(lithic_depth[available], na.rm = TRUE),
    lithic_missing = median(lithic_depth[!available], na.rm = TRUE)
  )

# test whether deep SOLUS NAs occur below the predicted restrictive depth
solus_rosetta %>%
  left_join(depth_dat, by = "cell") %>%
  filter(depth %in% c(60, 100, 150), !is.na(resdept)) %>%
  mutate(
    available = !is.na(sandtotal) & !is.na(claytotal) & !is.na(dbovendry),
    below_restriction = depth >= resdept
  ) %>%
  count(depth, available, below_restriction) %>%
  group_by(depth, available) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()


# convert SOLUS rasters to cell-by-depth Rosetta inputs
solus_rosetta <-
  as.data.frame(solus_levy, cells = TRUE, xy = TRUE, na.rm = FALSE) %>%
  pivot_longer(
    -c(cell, x, y),
    names_to = c("variable", "depth"),
    names_pattern = "(sandtotal|claytotal|dbovendry)_(\\d+)_cm_p",
    values_to = "value"
  ) %>%
  mutate(depth = as.integer(depth)) %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate(
    dbovendry = dbovendry / 100,
    silttotal = 100 - sandtotal - claytotal
  ) %>%
  arrange(cell, depth)

# inspect Rosetta input distributions by depth
solus_rosetta %>%
  group_by(depth) %>%
  summarise(
    n = n(),
    n_complete = sum(complete.cases(sandtotal, silttotal, claytotal, dbovendry)),
    across(c(sandtotal, silttotal, claytotal, dbovendry),
           list(min = ~min(.x, na.rm = TRUE), mean = ~mean(.x, na.rm = TRUE),
                max = ~max(.x, na.rm = TRUE)))
  )



# check derived silt and texture sums
solus_rosetta %>%
  summarise(
    silt_min = min(silttotal, na.rm = TRUE),
    silt_max = max(silttotal, na.rm = TRUE),
    sum_min = min(sandtotal + silttotal + claytotal, na.rm = TRUE),
    sum_max = max(sandtotal + silttotal + claytotal, na.rm = TRUE),
    n_silt_negative = sum(silttotal < 0, na.rm = TRUE)
  )

# count complete and unique Rosetta input combinations
solus_rosetta %>%
  drop_na(sandtotal, silttotal, claytotal, dbovendry) %>%
  summarise(
    n_predictions = n(),
    n_unique = n_distinct(pick(sandtotal, silttotal, claytotal, dbovendry))
  )


# ======================================================================================
# build unique Rosetta input combinations
rosetta_inputs <-
  solus_rosetta %>%
  drop_na(sandtotal, silttotal, claytotal, dbovendry) %>%
  distinct(sandtotal, silttotal, claytotal, dbovendry) %>%
  mutate(rosetta_id = row_number())

# run Rosetta on unique SOLUS soil combinations
system.time(
  rosetta_results <- run_rosetta(
    rosetta_inputs %>% select(sandtotal, silttotal, claytotal, dbovendry)
  )
)

# combine Rosetta inputs and outputs
rosetta_results <-
  bind_cols(rosetta_inputs, rosetta_results)

glimpse(rosetta_results)

# calculate Rosetta AWC for each unique soil combination
rosetta_results <-
  rosetta_results %>%
  mutate(
    alpha = 10^log10_alpha_mean,
    npar = 10^log10_npar_mean,
    m = 1 - 1 / npar,
    theta_33 = theta_r_mean + (theta_s_mean - theta_r_mean) /
      (1 + (alpha * (33 * 10.197))^npar)^m,
    theta_1500 = theta_r_mean + (theta_s_mean - theta_r_mean) /
      (1 + (alpha * (1500 * 10.197))^npar)^m,
    awc = theta_33 - theta_1500
  )

# plot awc ----
# attach Rosetta AWC to 15-cm SOLUS cells
awc_15 <-
  solus_rosetta %>%
  filter(depth == 15) %>%
  left_join(
    rosetta_results %>%
      select(sandtotal, silttotal, claytotal, dbovendry, awc),
    by = c("sandtotal", "silttotal", "claytotal", "dbovendry")
  )

# create 15-cm Rosetta AWC raster
awc_15_rast <- solus_levy[[1]]
values(awc_15_rast) <- NA_real_
awc_15_rast[awc_15$cell] <- awc_15$awc
names(awc_15_rast) <- "awc_15cm"

# plot Rosetta AWC at 15-cm depth
plot(awc_15_rast, main = "Volumetric water capcity at 15 cm")


# select contrasting 15-cm SOLUS sites based on sand content
test_sites <-
  solus_rosetta %>%
  filter(depth == 15, complete.cases(sandtotal, silttotal, claytotal, dbovendry)) %>%
  crossing(target_sand = c(40, 70, 95)) %>%
  group_by(target_sand) %>%
  slice_min(abs(sandtotal - target_sand), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(site = factor(paste0(target_sand, "% sand"),
                       levels = paste0(c(40, 70, 95), "% sand")))

test_sites %>%
  select(site, x, y, sandtotal, silttotal, claytotal, dbovendry)



# attach existing Rosetta predictions to selected sites
test_sites <-
  test_sites %>%
  left_join(rosetta_results,
            by = c("sandtotal", "silttotal", "claytotal", "dbovendry"))


# compare soil properties and Rosetta hydraulic parameters
test_sites %>%
  select(site, sandtotal, silttotal, claytotal, dbovendry,
         theta_r_mean, theta_s_mean, log10_alpha_mean,
         log10_npar_mean, log10_Ksat_mean)


# define pressure heads spanning near-saturation through very dry soil
pressure_head <- 10^seq(0, 5, length.out = 300)

# calculate van Genuchten water-retention curves
retention_curves <-
  test_sites %>%
  select(site, theta_r_mean, theta_s_mean, log10_alpha_mean, log10_npar_mean) %>%
  crossing(h_cm = pressure_head) %>%
  mutate(
    alpha = 10^log10_alpha_mean,
    n = 10^log10_npar_mean,
    m = 1 - 1 / n,
    theta = theta_r_mean + (theta_s_mean - theta_r_mean) /
      (1 + (alpha * h_cm)^n)^m,
    kpa = h_cm / 10.197
  )

# plot soil-water retention curves
ggplot(retention_curves, aes(kpa, theta, color = site)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = c(33, 1500), linetype = "dashed") +
  scale_x_log10() +
  labs(x = "Soil water suction (-kPa)", y = expression(theta~(cm^3~cm^-3)),
       color = "15-cm soil") +
  theme_bw()



# calculate water content at conventional field-capacity and wilting potentials
retention_points <-
  test_sites %>%
  select(site, theta_r_mean, theta_s_mean, log10_alpha_mean, log10_npar_mean) %>%
  crossing(kpa = c(33, 1500)) %>%
  mutate(
    h_cm = kpa * 10.197,
    alpha = 10^log10_alpha_mean,
    n = 10^log10_npar_mean,
    m = 1 - 1 / n,
    theta = theta_r_mean + (theta_s_mean - theta_r_mean) /
      (1 + (alpha * h_cm)^n)^m
  ) %>%
  select(site, kpa, theta)

retention_points


#
# define soil water potentials for comparison
water_potentials <- c(1, 5, 10, 33, 100, 300, 500, 1000, 1500)

# calculate water content across selected potentials
retention_table <-
  test_sites %>%
  select(site, theta_r_mean, theta_s_mean, log10_alpha_mean, log10_npar_mean) %>%
  crossing(kpa = water_potentials) %>%
  mutate(
    h_cm = kpa * 10.197,
    alpha = 10^log10_alpha_mean,
    n = 10^log10_npar_mean,
    m = 1 - 1 / n,
    theta = theta_r_mean + (theta_s_mean - theta_r_mean) /
      (1 + (alpha * h_cm)^n)^m
  )

# display water content by site and matric potential
retention_table %>%
  select(site, kpa, theta) %>%
  pivot_wider(names_from = kpa, values_from = theta, names_prefix = "kPa_")


# calculate water lost as soil dries from -10 kPa
retention_table %>%
  group_by(site) %>%
  mutate(theta_10 = theta[kpa == 10], water_lost_from_10 = theta_10 - theta) %>%
  ungroup() %>%
  filter(kpa >= 10) %>%
  select(site, kpa, theta, water_lost_from_10)



# ======================================================================================
# estimate whole-profile water storage  
# ======================================================================================
# attach Rosetta hydraulic parameters to every SOLUS cell and depth
solus_hydraulic <-
  solus_rosetta %>%
  left_join(
    rosetta_results %>%
      select(sandtotal, silttotal, claytotal, dbovendry,
             theta_r_mean, theta_s_mean, log10_alpha_mean, log10_npar_mean),
    by = c("sandtotal", "silttotal", "claytotal", "dbovendry")
  )

# attach restrictive depth and retain accessible soil depths
solus_theta_depth <-
  solus_theta %>%
  left_join(depth_dat %>% select(cell, resdept), by = "cell") %>%
  filter(!is.na(resdept), depth <= pmin(resdept, 150))

# integrate water storage to the deepest valid accessible depth
profile_storage <-
  solus_theta_depth %>%
  group_by(cell, x, y, kpa) %>%
  arrange(depth, .by_group = TRUE) %>%
  filter(n() >= 2) %>%
  summarise(
    profile_depth = max(depth),
    storage_mm = sum(diff(depth) * (head(theta, -1) + tail(theta, -1)) / 2) * 10,
    .groups = "drop"
  )

# summarize effective profile depths across Levy County
profile_storage %>%
  distinct(cell, profile_depth) %>%
  count(profile_depth) %>%
  mutate(pct = n / sum(n) * 100)

# calculate profile AWC between -33 and -1500 kPa
profile_awc <-
  profile_storage %>%
  filter(kpa %in% c(33, 1500)) %>%
  select(cell, x, y, profile_depth, kpa, storage_mm) %>%
  pivot_wider(names_from = kpa, values_from = storage_mm, names_prefix = "storage_") %>%
  mutate(awc_mm = storage_33 - storage_1500)

# map Rosetta profile AWC
awc_profile_rast <- solus_levy[[1]]
values(awc_profile_rast) <- NA_real_
awc_profile_rast[profile_awc$cell] <- profile_awc$awc_mm
names(awc_profile_rast) <- "awc_mm"

# plot
plot(awc_profile_rast, main = "AWC (mm), 0-150 cm")
