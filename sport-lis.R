
# install packages and libs w/ pacman ------------------------------------------
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(here, tidyverse, terra, sf, tigris, rosettaPTF, soiltexture)
# install rosetta from github; https://ncss-tech.github.io/rosettaPTF/
# remotes::install_github("ncss-tech/rosettaPTF")

# read SOLUS layers needed for Rosetta ----
(l = list.files("D:\\Projects\\fbrc\\data\\solus", full.names=TRUE, pattern = "\\.tif$"))
(solus_files <- l[str_detect(basename(l), "^(sandtotal|claytotal|dbovendry)_(0|5|15|30|60)_cm_p\\.tif$")] )
# read SOLUS layers as stack
solus <- rast(solus_files)


# get Levy County boundary for testing ---------------------------------------------------
levy <- 
  counties(state = "FL", cb = TRUE, resolution = "20m", year = 2025) %>%
  filter(NAME == "Levy") %>%
  st_transform(crs(rast(l[1])))


# crop and mask solus to Levy County
solus_levy <- 
  crop(solus, vect(levy)) %>%
  mask(vect(levy))


# read topsoil raster used by noah-lsm for sport-lis -------------------------------------
# https://ral.ucar.edu/model/noah-multiparameterization-land-surface-model-noah-mp-lsm
# soil look up table: https://ral.ucar.edu/sites/default/files/public/product-tool/noah-multiparameterization-land-surface-model-noah-mp-lsm/SOILPARM.TBL_.txt
(soil_top <- terra::rast("D:/Projects/RAPID/data/statsgo/STATSGO/topsoil30snew/topsoil30snew"))

# crop noah topsoil classes to southeastern US
soil_top_seus <- crop(soil_top, ext(-107, -74, 24, 40))

# remove missing/water cells and plot soil texture classes
soil_top_seus <- classify(soil_top_seus, rbind(c(0, NA), c(14, NA)))
plot(soil_top_seus)


# Crop and mask Noah topsoil classes to Levy County
soil_top_levy <- soil_top_seus %>%
  crop(vect(st_transform(levy, 4326))) %>%
  mask(vect(st_transform(levy, 4326)))

# add noah soil texture labels
soil_top_levy_texture <- as.factor(soil_top_levy)
levels(soil_top_levy_texture) <- data.frame(
  ID = c(1, 3, 12, 13),
  texture = c("Sand", "Sandy loam", "Clay", "Organic material"))

# plot soil texture classes
plot(soil_top_levy_texture, main = "Noah soil texture")

# ====================================================================================
# Prepare 0-cm Rosetta inputs
solus_levy_0 <- c(
  solus_levy[["sandtotal_0_cm_p"]],
  100 - solus_levy[["sandtotal_0_cm_p"]] - solus_levy[["claytotal_0_cm_p"]],
  solus_levy[["claytotal_0_cm_p"]],
  solus_levy[["dbovendry_0_cm_p"]] / 100)
# add layer names
names(solus_levy_0) <- c("sandtotal", "silttotal", "claytotal", "dbovendry")
# plot 
plot(solus_levy_0)

# Retain unique SOLUS soil combinations
solus_levy_0_unique <- as.data.frame(solus_levy_0, na.rm = TRUE) %>%
  distinct() %>%
  mutate(soil_id = row_number())

# run rosetta once per unique soil combination
rosetta_levy_0 <- run_rosetta(
  solus_levy_0_unique %>% select(-soil_id),
  vars = c("sandtotal", "silttotal", "claytotal", "dbovendry")) %>%
  bind_cols(solus_levy_0_unique, .) %>%
  mutate(alpha = 10^log10_alpha_mean,
         npar = 10^log10_npar_mean,
         m = 1 - 1 / npar,
         ksat = 10^log10_Ksat_mean * 10,
         theta_1500 = theta_r_mean + (theta_s_mean - theta_r_mean) /
           (1 + (alpha * (1500 * 10.197))^npar)^m)

# calculate unsaturated hydraulic conductivity
k_vgm <- function(h, alpha, n, m, ksat) {
  se <- (1 + (alpha * abs(h))^n)^(-m)
  ksat * sqrt(se) * (1 - (1 - se^(1 / m))^m)^2
}

# find pressure head where hydraulic conductivity reaches 1 mm/day
find_h_fc <- function(alpha, n, m, ksat, target = 1) {
  if (ksat <= target) return(0)
  uniroot(\(h) k_vgm(h, alpha, n, m, ksat) - target, c(-1e6, 0))$root
}

# calculate field capacity and available water capacity
rosetta_levy_0 <- rosetta_levy_0 %>%
  mutate(h_fc = pmap_dbl(list(alpha, npar, m, ksat), find_h_fc),
         theta_fc = theta_r_mean + (theta_s_mean - theta_r_mean) /
           (1 + (alpha * abs(h_fc))^npar)^m,
         awc = theta_fc - theta_1500)


# ========================================================================================
# match each solus cell to its rosetta prediction
cells <- as.data.frame(solus_levy_0, cells = TRUE, na.rm = TRUE) %>%
  left_join(rosetta_levy_0 %>%
              select(sandtotal, silttotal, claytotal, dbovendry, soil_id),
            by = c("sandtotal", "silttotal", "claytotal", "dbovendry"))

# create soil-id raster
soil_id <- solus_levy_0[[1]]
values(soil_id) <- NA
soil_id[cells$cell] <- cells$soil_id

# map rosetta hydraulic properties to solus grid
rosetta_levy_rast <- c(
  subst(soil_id, rosetta_levy_0$soil_id, rosetta_levy_0$theta_1500),
  subst(soil_id, rosetta_levy_0$soil_id, rosetta_levy_0$theta_fc),
  subst(soil_id, rosetta_levy_0$soil_id, rosetta_levy_0$theta_s_mean),
  subst(soil_id, rosetta_levy_0$soil_id, rosetta_levy_0$awc))
names(rosetta_levy_rast) <- c("pwp_rosetta", "fc_rosetta", "sat_rosetta", "awc_rosetta")

# inspect local hydraulic variation
global(rosetta_levy_rast, c("min", "mean", "max"), na.rm = TRUE)
plot(rosetta_levy_rast)




# check hydraulic ordering and awc distribution
global(rosetta_levy_rast[["awc_rosetta"]],
       c("min", "mean", "max"), na.rm = TRUE)

global(rosetta_levy_rast[["fc_rosetta"]] <= rosetta_levy_rast[["pwp_rosetta"]],
       "sum", na.rm = TRUE)

global(rosetta_levy_rast[["fc_rosetta"]] >= rosetta_levy_rast[["sat_rosetta"]],
       "sum", na.rm = TRUE)




# read in a sport-lis RSM tile ---------------------------------------------------------
# read sport-lis 0-10 cm relative soil moisture
(sport_rsm <- rast("D:/Projects/RAPID/data/sport/20260923_0000_sport_lis_rsm0-10cm_conus3km_float_wgs84.tif"))

# inspect raster and values
sport_rsm
global(sport_rsm, c("min", "mean", "max"), na.rm = TRUE)
plot(sport_rsm)

# remove sport-lis fill values
sport_rsm[sport_rsm == 9999] <- NA

# inspect valid rsm values
global(sport_rsm, c("min", "mean", "max"), na.rm = TRUE)
# inspect rsm distribution
quantile(values(sport_rsm, na.rm = TRUE),
         probs = c(0, .01, .25, .5, .75, .99, 1))


# crop and mask sport-lis rsm to levy county
sport_rsm_levy <- sport_rsm %>%
  crop(vect(st_transform(levy, crs(sport_rsm)))) %>%
  mask(vect(st_transform(levy, crs(sport_rsm))))

# inspect levy county rsm
global(sport_rsm_levy, c("min", "mean", "max"), na.rm = TRUE)
values(sport_rsm_levy, na.rm = TRUE)

# project sport rsm to the solus/rosetta grid --------------------------------------------
sport_rsm_solus <- project(sport_rsm_levy, rosetta_levy_rast, method = "near") / 100
names(sport_rsm_solus) <- "rsm_sport"

# calculate localized volumetric soil moisture
theta_localized <- rosetta_levy_rast[["pwp_rosetta"]] +
  sport_rsm_solus * (rosetta_levy_rast[["sat_rosetta"]] -
                       rosetta_levy_rast[["pwp_rosetta"]])
names(theta_localized) <- "vsm_localized"

# calculate fraction of plant-available water
paw_fraction <- (theta_localized - rosetta_levy_rast[["pwp_rosetta"]]) /
  rosetta_levy_rast[["awc_rosetta"]]
paw_fraction <- clamp(paw_fraction, 0, 1)
names(paw_fraction) <- "paw_fraction"

# combine and inspect moisture outputs
moisture_levy <- c(sport_rsm_solus, theta_localized, paw_fraction)
global(moisture_levy, c("min", "mean", "max"), na.rm = TRUE)
plot(moisture_levy)


# calculate relative wetness corresponding to local field capacity ----------------------
rsm_at_fc <- (rosetta_levy_rast[["fc_rosetta"]] -
                rosetta_levy_rast[["pwp_rosetta"]]) /
  (rosetta_levy_rast[["sat_rosetta"]] -
     rosetta_levy_rast[["pwp_rosetta"]])
names(rsm_at_fc) <- "rsm_at_fc"

# inspect local field-capacity position within the pwp-saturation range
global(rsm_at_fc, c("min", "mean", "max"), na.rm = TRUE)
plot(rsm_at_fc)



# =======================================================================================
# use the noah lookup table 
# https://ral.ucar.edu/sites/default/files/public/product-tool/noah-multiparameterization-land-surface-model-noah-mp-lsm/SOILPARM.TBL_.txt

# noah stas hydraulic parameters by soil class
noah_soil <- tibble(
  class = 1:19,
  bb = c(2.79,4.26,4.74,5.33,5.33,5.25,6.77,8.72,8.17,10.73,10.39,11.55,5.25,
         0,2.79,4.26,11.55,2.79,2.79),
  wltsmc = c(.010,.028,.047,.084,.084,.066,.067,.120,.103,.100,.126,.138,.066,
             0,.006,.028,.030,.006,.010),
  refsmc = c(.236,.383,.383,.360,.383,.329,.314,.387,.382,.338,.404,.412,.329,
             0,.170,.283,.454,.170,.236),
  maxsmc = c(.339,.421,.434,.476,.476,.439,.404,.464,.465,.406,.468,.468,.439,
             1,.200,.421,.468,.200,.339),
  satpsi = c(.069,.036,.141,.759,.759,.355,.135,.617,.263,.098,.324,.468,.355,
             0,.069,.036,.468,.069,.069))



# map noah parameters to the solus grid
noah_params <- c(
  subst(soil_top_levy_solus, noah_soil$class, noah_soil$wltsmc),
  subst(soil_top_levy_solus, noah_soil$class, noah_soil$maxsmc),
  subst(soil_top_levy_solus, noah_soil$class, noah_soil$satpsi),
  subst(soil_top_levy_solus, noah_soil$class, noah_soil$bb))
names(noah_params) <- c("wltsmc", "maxsmc", "satpsi", "bb")

# reconstruct noah volumetric soil moisture from sport rsm
vsm_noah <- noah_params[["wltsmc"]] + sport_rsm_solus *
  (noah_params[["maxsmc"]] - noah_params[["wltsmc"]])
names(vsm_noah) <- "vsm_noah"

# convert noah moisture state to pressure head in meters
h_noah <- noah_params[["satpsi"]] *
  (vsm_noah / noah_params[["maxsmc"]])^(-noah_params[["bb"]])
names(h_noah) <- "h_noah_m"

# convert pressure head to cm and evaluate the rosetta retention curve
h_noah_cm <- h_noah * 100

theta_localized_h <- rosetta_levy_rast[["pwp_rosetta"]]
theta_localized_h <- rosetta_levy_0$theta_r_mean[1] # don't run this yet

# inspect inferred noah pressure head
global(h_noah, c("min", "mean", "max"), na.rm = TRUE)
quantile(values(h_noah, na.rm = TRUE),
         probs = c(0, .05, .25, .5, .75, .95, 1))

# map rosetta retention parameters to the solus grid -------------------------------------
rosetta_vg <- c(
  subst(soil_id, rosetta_levy_0$soil_id, rosetta_levy_0$theta_r_mean),
  subst(soil_id, rosetta_levy_0$soil_id, rosetta_levy_0$alpha),
  subst(soil_id, rosetta_levy_0$soil_id, rosetta_levy_0$npar),
  subst(soil_id, rosetta_levy_0$soil_id, rosetta_levy_0$m))
names(rosetta_vg) <- c("theta_r", "alpha", "n", "m")

# convert noah pressure head to cm for the rosetta curve
h_noah_cm <- h_noah * 100

# evaluate local rosetta retention curve at the noah pressure head
theta_localized_h <- rosetta_vg[["theta_r"]] +
  (rosetta_levy_rast[["sat_rosetta"]] - rosetta_vg[["theta_r"]]) /
  (1 + (rosetta_vg[["alpha"]] * h_noah_cm)^rosetta_vg[["n"]])^rosetta_vg[["m"]]
names(theta_localized_h) <- "vsm_localized_h"

# calculate raw fraction of plant-available water
paw_fraction_h <- (theta_localized_h - rosetta_levy_rast[["pwp_rosetta"]]) /
  rosetta_levy_rast[["awc_rosetta"]]
names(paw_fraction_h) <- "paw_fraction_h"

# inspect pressure-head localization before bounding paw
moisture_h <- c(sport_rsm_solus, h_noah, theta_localized_h, paw_fraction_h)
global(moisture_h, c("min", "mean", "max"), na.rm = TRUE)
plot(moisture_h)


# classify localized moisture relative to pwp and field capacity ----------------------
moisture_class <- ifel(theta_localized_h < rosetta_levy_rast[["pwp_rosetta"]], 1,
                       ifel(theta_localized_h < rosetta_levy_rast[["fc_rosetta"]], 2, 3))
names(moisture_class) <- "moisture_class"

# summarize landscape moisture status
freq(moisture_class) %>%
  as_tibble() %>%
  mutate(status = recode(as.character(value),
                         `1` = "below pwp",
                         `2` = "pwp to fc",
                         `3` = "above fc"),
         percent = 100 * count / sum(count)) %>%
  select(status, count, percent)

plot(moisture_class)
# inspect raw fraction of plant-available water
quantile(values(paw_fraction_h, na.rm = TRUE),
         probs = c(0, .05, .25, .5, .75, .95, .99, 1))


# summarize rsm and pressure head by noah texture class
as.data.frame(c(soil_top_levy_solus, sport_rsm_solus, h_noah), na.rm = TRUE) %>%
  rename(class = 1) %>%
  group_by(class) %>%
  summarise(n = n(),
            rsm_mean = mean(rsm_sport),
            h_median = median(h_noah_m),
            h_p95 = quantile(h_noah_m, .95),
            h_max = max(h_noah_m),
            .groups = "drop") %>%
  left_join(noah_soil %>% select(class, bb, wltsmc, maxsmc, satpsi), by = "class")

# ====================================================================================
# Map Noah hydraulic parameters to Levy soil classes
noah_levy <- 
  c(subst(soil_top_levy, noah_soil$class, noah_soil$wltsmc),
    #subst(soil_top_levy, noah_soil$class, noah_soil$refsmc),
    subst(soil_top_levy, noah_soil$class, noah_soil$maxsmc))
# name layers
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
# project Noah classes to the SOLUS/Rosetta grid
soil_top_levy_solus <- project(soil_top_levy, rosetta_levy_rast, method = "near")

# retain Rosetta cells Noah classifies as sand
rosetta_noah_sand <- mask(rosetta_levy_rast, soil_top_levy_solus == 1)

# summarize Rosetta properties within Noah "sand"
global(rosetta_noah_sand, c("min", "mean", "max"), na.rm = TRUE)
global(rosetta_noah_sand, median, na.rm = TRUE)

# examine SOLUS variability within Noah "sand"
solus_noah_sand <- mask(solus_levy_0, soil_top_levy_solus == 1)

# summarize SOLUS properties within Noah "sand"
global(solus_noah_sand, c("min", "mean", "max"), na.rm = TRUE)
global(solus_noah_sand, median, na.rm = TRUE)

# summarize Noah hydraulic properties
global(noah_levy, c("mean", "sd", "min", "max"), na.rm = TRUE)
# Summarize Rosetta hydraulic properties
global(rosetta_levy_rast, c("mean", "sd", "min", "max"), na.rm = TRUE)

# map Noah hydraulic properties onto the SOLUS grid
noah_levy_solus <- c(
  subst(soil_top_levy_solus, noah_soil$class, noah_soil$wltsmc),
  subst(soil_top_levy_solus, noah_soil$class, noah_soil$maxsmc))
names(noah_levy_solus) <- c("pwp_noah", "sat_noah")

# combine Noah texture class and Rosetta hydraulic properties
hydraulic_compare <- as.data.frame(c(soil_top_levy_solus, rosetta_levy_rast), na.rm = TRUE) %>%
  rename(class = 1) %>%
  mutate(texture = factor(class, levels = c(1, 3, 12, 13),
                          labels = c("Sand", "Sandy loam", "Clay", "Organic material")))


# add Noah/STAS hydraulic values to each texture class
noah_plot <- noah_soil %>%
  filter(class %in% c(1, 3, 12, 13)) %>%
  mutate(texture = factor(class, levels = c(1, 3, 12, 13),
                          labels = c("Sand", "Sandy loam", "Clay", "Organic material")))


# compare Rosetta PWP distributions with Noah wilting thresholds
ggplot(hydraulic_compare, aes(texture, pwp_rosetta)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(data = noah_plot, aes(texture, wltsmc, shape = "Noah/STAS"), color='red', size = 3) +
  scale_fill_discrete(name = NULL) +
  scale_shape_manual(name = NULL, values = 18) +
  labs(x = "Noah soil texture", y = "Wilting point/threshold (m³/m³)") +
  theme_bw() +
  theme(legend.position = 'top')


# compare Rosetta and Noah saturated water content
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


