
pacman::p_load(here, tidyverse, terra, sf, tigris, rosettaPTF, soiltexture)


# select SOLUS layers needed for initial Rosetta testing ----
(l = list.files("D:\\Projects\\fbrc\\data\\solus", full.names=TRUE, pattern = "\\.tif$"))
(solus_files <- l[str_detect(basename(l), "^(sandtotal|claytotal|dbovendry)_(0|5|15|30|60)_cm_p\\.tif$")] )
# read SOLUS layers as stack
solus <- rast(solus_files)

# get Levy County boundary
levy <- counties(state = "FL", cb = TRUE, resolution = "20m", year = 2025) %>%
  filter(NAME == "Levy") %>%
  st_transform(crs(rast(l[1])))


# crop and mask to Levy County
solus_levy <- crop(solus, vect(levy)) %>%
  mask(vect(levy))


# Read Noah topsoil texture raster ======================================================
soil_top <- terra::rast("D:/Projects/RAPID/data/statsgo/STATSGO/topsoil30snew/topsoil30snew")
soil_top


# Crop Noah topsoil classes to southeastern US
soil_top_seus <- crop(soil_top, ext(-107, -74, 24, 40))

# Inspect soil classes and spatial pattern
freq(soil_top_seus)
plot(soil_top_seus)


# Remove missing/water cells and plot soil texture classes
soil_top_seus <- classify(soil_top_seus, rbind(c(0, NA), c(14, NA)))
plot(soil_top_seus)


# Crop and mask Noah topsoil classes to Levy County
soil_top_levy <- soil_top_seus %>%
  crop(vect(st_transform(levy, 4326))) %>%
  mask(vect(st_transform(levy, 4326)))

# Inspect soil classes represented in Levy County
freq(soil_top_levy)
plot(soil_top_levy)
# Add Noah soil texture labels
soil_top_levy_texture <- as.factor(soil_top_levy)
levels(soil_top_levy_texture) <- data.frame(
  ID = c(1, 3, 12, 13),
  texture = c("Sand", "Sandy loam", "Clay", "Organic material")
)

# Plot soil texture classes
plot(soil_top_levy_texture, main = "Noah soil texture")

# ====================================================================================
# Prepare 0-cm Rosetta inputs
solus_levy_0 <- c(
  solus_levy[["sandtotal_0_cm_p"]],
  100 - solus_levy[["sandtotal_0_cm_p"]] - solus_levy[["claytotal_0_cm_p"]],
  solus_levy[["claytotal_0_cm_p"]],
  solus_levy[["dbovendry_0_cm_p"]] / 100
)
names(solus_levy_0) <- c("sandtotal", "silttotal", "claytotal", "dbovendry")
# plot 
plot(solus_levy_0)

# Retain unique SOLUS soil combinations
solus_levy_0_unique <- as.data.frame(solus_levy_0, na.rm = TRUE) %>%
  distinct() %>%
  mutate(soil_id = row_number())

# Run Rosetta once per unique soil combination
rosetta_levy_0 <- run_rosetta(
  solus_levy_0_unique %>% select(-soil_id),
  vars = c("sandtotal", "silttotal", "claytotal", "dbovendry")
) %>%
  bind_cols(solus_levy_0_unique, .)



# Calculate Rosetta PWP at -1500 kPa
rosetta_levy_0 <- rosetta_levy_0 %>%
  mutate(
    alpha = 10^log10_alpha_mean,
    npar = 10^log10_npar_mean,
    m = 1 - 1 / npar,
    theta_1500 = theta_r_mean + (theta_s_mean - theta_r_mean) /
      (1 + (alpha * (1500 * 10.197))^npar)^m
  )


# =======================================================================================
# Match each SOLUS cell to its Rosetta prediction
cells <- as.data.frame(solus_levy_0, cells = TRUE, na.rm = TRUE) %>%
  left_join(rosetta_levy_0 %>%
              select(sandtotal, silttotal, claytotal, dbovendry, soil_id),
            by = c("sandtotal", "silttotal", "claytotal", "dbovendry"))

# Create soil-ID raster
soil_id <- solus_levy_0[[1]]
values(soil_id) <- NA
soil_id[cells$cell] <- cells$soil_id

# Map Rosetta PWP and saturation to SOLUS grid
rosetta_levy_rast <- c(
  subst(soil_id, rosetta_levy_0$soil_id, rosetta_levy_0$theta_1500),
  subst(soil_id, rosetta_levy_0$soil_id, rosetta_levy_0$theta_s_mean)
)
names(rosetta_levy_rast) <- c("pwp_rosetta", "sat_rosetta")

# Inspect local hydraulic variation
global(rosetta_levy_rast, c("min", "mean", "max"), na.rm = TRUE)
plot(rosetta_levy_rast)


# =======================================================================================
# https://ral.ucar.edu/sites/default/files/public/product-tool/noah-multiparameterization-land-surface-model-noah-mp-lsm/SOILPARM.TBL_.txt

# Noah STAS hydraulic parameters by soil class
noah_soil <- tibble(
  class = 1:19,
  wltsmc = c(.010,.028,.047,.084,.084,.066,.067,.120,.103,.100,.126,.138,.066,
             0,.006,.028,.030,.006,.010),
  refsmc = c(.236,.383,.383,.360,.383,.329,.314,.387,.382,.338,.404,.412,.329,
             0,.170,.283,.454,.170,.236),
  maxsmc = c(.339,.421,.434,.476,.476,.439,.404,.464,.465,.406,.468,.468,.439,
             1,.200,.421,.468,.200,.339)
)

# Map Noah hydraulic parameters to Levy soil classes
noah_levy <- c(
  subst(soil_top_levy, noah_soil$class, noah_soil$wltsmc),
  #subst(soil_top_levy, noah_soil$class, noah_soil$refsmc),
  subst(soil_top_levy, noah_soil$class, noah_soil$maxsmc)
)
names(noah_levy) <- c("pwp_noah", 
                      #"ref_noah", 
                      "sat_noah")

# Inspect Noah hydraulic properties
global(noah_levy, c("min", "mean", "max"), na.rm = TRUE)

# plot
par(mfrow = c(1,2))
plot(noah_levy[[1]], main="pwp_noah")
plot(noah_levy[[2]], main="sat_noah")
par(mfrow = c(1, 1))


# ======================================================================================
# Project Noah classes to the SOLUS/Rosetta grid
soil_top_levy_solus <- project(soil_top_levy, rosetta_levy_rast, method = "near")

# Retain Rosetta cells Noah classifies as sand
rosetta_noah_sand <- mask(rosetta_levy_rast, soil_top_levy_solus == 1)

# Summarize Rosetta properties within Noah "sand"
global(rosetta_noah_sand, c("min", "mean", "max"), na.rm = TRUE)
global(rosetta_noah_sand, median, na.rm = TRUE)

# Examine SOLUS variability within Noah "sand"
solus_noah_sand <- mask(solus_levy_0, soil_top_levy_solus == 1)

# Summarize SOLUS properties within Noah "sand"
global(solus_noah_sand, c("min", "mean", "max"), na.rm = TRUE)
global(solus_noah_sand, median, na.rm = TRUE)



# Summarize Noah hydraulic properties
global(noah_levy, c("mean", "sd", "min", "max"), na.rm = TRUE)
# Summarize Rosetta hydraulic properties
global(rosetta_levy_rast, c("mean", "sd", "min", "max"), na.rm = TRUE)

# Map Noah hydraulic properties onto the SOLUS grid
noah_levy_solus <- c(
  subst(soil_top_levy_solus, noah_soil$class, noah_soil$wltsmc),
  subst(soil_top_levy_solus, noah_soil$class, noah_soil$maxsmc)
)
names(noah_levy_solus) <- c("pwp_noah", "sat_noah")

# Combine Noah texture class and Rosetta hydraulic properties
hydraulic_compare <- as.data.frame(c(soil_top_levy_solus, rosetta_levy_rast), na.rm = TRUE) %>%
  rename(class = 1) %>%
  mutate(texture = factor(class, levels = c(1, 3, 12, 13),
                          labels = c("Sand", "Sandy loam", "Clay", "Organic material")))

# Plot Rosetta PWP within each Noah soil class
ggplot(hydraulic_compare, aes(texture, pwp_rosetta)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "Noah soil texture", y = "Rosetta PWP (m³/m³)") +
  theme_bw()

# Plot Rosetta saturation within each Noah soil class
ggplot(hydraulic_compare, aes(texture, sat_rosetta)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "Noah soil texture", y = "Rosetta saturation (m³/m³)") +
  theme_bw()



# Add Noah/STAS hydraulic values to each texture class
noah_plot <- noah_soil %>%
  filter(class %in% c(1, 3, 12, 13)) %>%
  mutate(texture = factor(class, levels = c(1, 3, 12, 13),
                          labels = c("Sand", "Sandy loam", "Clay", "Organic material")))


# Compare Rosetta PWP distributions with Noah wilting thresholds
ggplot(hydraulic_compare, aes(texture, pwp_rosetta)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(data = noah_plot, aes(texture, wltsmc, shape = "Noah/STAS"), color='red', size = 3) +
  scale_fill_discrete(name = NULL) +
  scale_shape_manual(name = NULL, values = 18) +
  labs(x = "Noah soil texture", y = "Wilting point/threshold (m³/m³)") +
  theme_bw() +
  theme(legend.position = 'top')


# Compare Rosetta and Noah saturated water content
ggplot(hydraulic_compare, aes(texture, sat_rosetta)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(data = noah_plot, aes(texture, maxsmc, shape = "Noah/STAS"), color='red', size = 3) +
  scale_fill_discrete(name = NULL) +
  scale_shape_manual(name = NULL, values = 18) +
  labs(x = "Noah soil texture", y = "Saturated water content (m³/m³)") +
  theme_bw() + 
  theme(legend.position = 'top')

# At this point, I think we have a very testable approach: 
# take an actual SPoRT VSM/RSM day over Levy, reconstruct RSM ourselves from WLTSMC/MAXSMC, 
# transfer that relative state to Rosetta PWP/θsat, and see what the resulting ~100-m moisture 
# surface looks like. That would tell us very quickly whether this approach produces physically 
# reasonable spatial variation before we invest further in it.


