# install packages and libraries -------------------------------------------------------


if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(here, tidyverse, terra, sf, tigris, rosettaPTF)


# read solus layers -------------------------------------------------------------------
l <- list.files("D:/Projects/fbrc/data/solus", full.names = TRUE, pattern = "\\.tif$")
solus_files <- l[str_detect(basename(l),
                            "^(sandtotal|claytotal|dbovendry)_(0|5|15|30|60)_cm_p\\.tif$")]
solus <- rast(solus_files)

# get levy county boundary
levy <- counties(state = "FL", cb = TRUE, resolution = "20m", year = 2025) %>%
  filter(NAME == "Levy") %>%
  st_transform(crs(solus))

# crop and mask solus to levy county
solus_levy <- solus %>%
  crop(vect(levy)) %>%
  mask(vect(levy))


# prepare 0-cm rosetta inputs ---------------------------------------------------------
solus_levy_0 <- c(
  solus_levy[["sandtotal_0_cm_p"]],
  100 - solus_levy[["sandtotal_0_cm_p"]] - solus_levy[["claytotal_0_cm_p"]],
  solus_levy[["claytotal_0_cm_p"]],
  solus_levy[["dbovendry_0_cm_p"]] / 100
)
names(solus_levy_0) <- c("sandtotal", "silttotal", "claytotal", "dbovendry")

# retain unique solus soil combinations
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

# find pressure head where conductivity reaches 1 mm/day
find_h_fc <- function(alpha, n, m, ksat, target = 1) {
  if (ksat <= target) return(0)
  uniroot(\(h) k_vgm(h, alpha, n, m, ksat) - target, c(-1e6, 0))$root
}

# derive field capacity and available water capacity
rosetta_levy_0 <- rosetta_levy_0 %>%
  mutate(h_fc = pmap_dbl(list(alpha, npar, m, ksat), find_h_fc),
         theta_fc = theta_r_mean + (theta_s_mean - theta_r_mean) /
           (1 + (alpha * abs(h_fc))^npar)^m,
         awc = theta_fc - theta_1500)


# map rosetta properties back to solus grid ------------------------------------------
cells <- as.data.frame(solus_levy_0, cells = TRUE, na.rm = TRUE) %>%
  left_join(rosetta_levy_0 %>%
              select(sandtotal, silttotal, claytotal, dbovendry, soil_id),
            by = c("sandtotal", "silttotal", "claytotal", "dbovendry"))

# create soil-id raster
soil_id <- solus_levy_0[[1]]
values(soil_id) <- NA
soil_id[cells$cell] <- cells$soil_id

# map local hydraulic properties
rosetta_levy_rast <- c(
  subst(soil_id, rosetta_levy_0$soil_id, rosetta_levy_0$theta_1500),
  subst(soil_id, rosetta_levy_0$soil_id, rosetta_levy_0$theta_fc),
  subst(soil_id, rosetta_levy_0$soil_id, rosetta_levy_0$theta_s_mean),
  subst(soil_id, rosetta_levy_0$soil_id, rosetta_levy_0$awc)
)
names(rosetta_levy_rast) <- c("pwp_rosetta", "fc_rosetta",
                              "sat_rosetta", "awc_rosetta")

# inspect local hydraulic properties
global(rosetta_levy_rast, c("min", "mean", "max"), na.rm = TRUE)


# read noah soil classes --------------------------------------------------------------
soil_top <- rast("D:/Projects/RAPID/data/statsgo/STATSGO/topsoil30snew/topsoil30snew")

# crop and mask noah soil classes to levy county
soil_top_levy <- soil_top %>%
  crop(vect(st_transform(levy, crs(soil_top)))) %>%
  mask(vect(st_transform(levy, crs(soil_top)))) %>%
  classify(rbind(c(0, NA), c(14, NA)))

# add noah soil texture labels
soil_top_levy_texture <- as.factor(soil_top_levy)
levels(soil_top_levy_texture) <- data.frame(
  ID = c(1, 3, 12, 13),
  texture = c("Sand", "Sandy loam", "Clay", "Organic material")
)
# plot texture
plot(soil_top_levy_texture, main='Noah texture classes')

# define noah stas moisture parameters ------------------------------------------------
noah_soil <- tibble(
  class = 1:19,
  wltsmc = c(.010,.028,.047,.084,.084,.066,.067,.120,.103,.100,.126,.138,.066,
             0,.006,.028,.030,.006,.010),
  refsmc = c(.236,.383,.383,.360,.383,.329,.314,.387,.382,.338,.404,.412,.329,
             0,.170,.283,.454,.170,.236),
  maxsmc = c(.339,.421,.434,.476,.476,.439,.404,.464,.465,.406,.468,.468,.439,
             1,.200,.421,.468,.200,.339)
) %>%
  mutate(rsm_ref = (refsmc - wltsmc) / (maxsmc - wltsmc))

# map reference rsm on the native noah soil grid
rsm_ref_noah <- subst(soil_top_levy, noah_soil$class, noah_soil$rsm_ref)
names(rsm_ref_noah) <- "rsm_ref"

# project noah reference rsm to the solus grid
rsm_ref <- project(rsm_ref_noah, rosetta_levy_rast, method = "near")

# inspect noah reference-moisture position
global(rsm_ref, c("min", "mean", "max"), na.rm = TRUE)
plot(rsm_ref, main="Reference soil moisture")

# plot noah texture classes and reference soil moisture ---
par(mfrow = c(1, 2))
plot(soil_top_levy_texture, main = "Noah texture classes")
plot(rsm_ref, main = "Reference soil moisture")
par(mfrow = c(1, 1))

# read wet and dry sport-lis periods --------------------------------------------------
sport_files <- list.files(here("data", "sport"), pattern = "\\.tif$",
                          full.names = TRUE)

wet_files <- sport_files[str_extract(basename(sport_files), "^\\d{8}") %in%
                           c("20260921", "20260922", "20260923")]
dry_files <- sport_files[str_extract(basename(sport_files), "^\\d{8}") %in%
                           c("20260927", "20260928", "20260929")]

# read and name wet-period rsm
wet_rsm <- rast(wet_files)
names(wet_rsm) <- str_extract(basename(wet_files), "^\\d{8}")
wet_rsm[wet_rsm == 9999] <- NA

# read and name dry-period rsm
dry_rsm <- rast(dry_files)
names(dry_rsm) <- str_extract(basename(dry_files), "^\\d{8}")
dry_rsm[dry_rsm == 9999] <- NA

# crop, mask, and project wet rsm to the solus grid
wet_rsm_solus <- wet_rsm %>%
  crop(vect(st_transform(levy, crs(wet_rsm)))) %>%
  mask(vect(st_transform(levy, crs(wet_rsm)))) %>%
  project(rosetta_levy_rast, method = "near") / 100

# crop, mask, and project dry rsm to the solus grid
dry_rsm_solus <- dry_rsm %>%
  crop(vect(st_transform(levy, crs(dry_rsm)))) %>%
  mask(vect(st_transform(levy, crs(dry_rsm)))) %>%
  project(rosetta_levy_rast, method = "near") / 100

# inspect wet and dry rsm
global(wet_rsm_solus, c("min", "mean", "max"), na.rm = TRUE)
global(dry_rsm_solus, c("min", "mean", "max"), na.rm = TRUE)

plot(wet_rsm_solus, range = c(0, 1))
plot(dry_rsm_solus, range = c(0, 1))


# localize wet-period rsm -------------------------------------------------------------
theta_wet_low <- rosetta_levy_rast[["pwp_rosetta"]] +
  (wet_rsm_solus / rsm_ref) *
  (rosetta_levy_rast[["fc_rosetta"]] - rosetta_levy_rast[["pwp_rosetta"]])

theta_wet_high <- rosetta_levy_rast[["fc_rosetta"]] +
  ((wet_rsm_solus - rsm_ref) / (1 - rsm_ref)) *
  (rosetta_levy_rast[["sat_rosetta"]] - rosetta_levy_rast[["fc_rosetta"]])

# combine wet-period piecewise localized vsm
vsm_wet <- ifel(wet_rsm_solus <= rsm_ref, theta_wet_low, theta_wet_high)
names(vsm_wet) <- names(wet_rsm_solus)

# calculate and bound wet-period fraction paw
paw_wet <- (vsm_wet - rosetta_levy_rast[["pwp_rosetta"]]) /
  rosetta_levy_rast[["awc_rosetta"]]
paw_wet <- clamp(paw_wet, 0, 1)
names(paw_wet) <- names(wet_rsm_solus)


# localize dry-period rsm -------------------------------------------------------------
theta_dry_low <- rosetta_levy_rast[["pwp_rosetta"]] +
  (dry_rsm_solus / rsm_ref) *
  (rosetta_levy_rast[["fc_rosetta"]] - rosetta_levy_rast[["pwp_rosetta"]])

theta_dry_high <- rosetta_levy_rast[["fc_rosetta"]] +
  ((dry_rsm_solus - rsm_ref) / (1 - rsm_ref)) *
  (rosetta_levy_rast[["sat_rosetta"]] - rosetta_levy_rast[["fc_rosetta"]])

# combine dry-period piecewise localized vsm
vsm_dry <- ifel(dry_rsm_solus <= rsm_ref, theta_dry_low, theta_dry_high)
names(vsm_dry) <- names(dry_rsm_solus)

# calculate and bound dry-period fraction paw
paw_dry <- (vsm_dry - rosetta_levy_rast[["pwp_rosetta"]]) /
  rosetta_levy_rast[["awc_rosetta"]]
paw_dry <- clamp(paw_dry, 0, 1)
names(paw_dry) <- names(dry_rsm_solus)


# compare wet and dry periods ---------------------------------------------------------
global(vsm_wet, c("min", "mean", "max"), na.rm = TRUE)
global(vsm_dry, c("min", "mean", "max"), na.rm = TRUE)

global(paw_wet, c("min", "mean", "max"), na.rm = TRUE)
global(paw_dry, c("min", "mean", "max"), na.rm = TRUE)

plot(paw_wet, range = c(0, 1))
plot(paw_dry, range = c(0, 1))


# summarize fraction of cells above the noah reference point --------------------------
wet_segment <- ifel(wet_rsm_solus <= rsm_ref, 0, 1)
dry_segment <- ifel(dry_rsm_solus <= rsm_ref, 0, 1)

global(wet_segment, "mean", na.rm = TRUE) * 100
global(dry_segment, "mean", na.rm = TRUE) * 100


# inspect wet and dry sport-lis rsm
global(wet_rsm_solus, c("min", "mean", "max"), na.rm = TRUE)
global(dry_rsm_solus, c("min", "mean", "max"), na.rm = TRUE)





# compare representative wet and dry days
compare_days <- c(
  wet_rsm_solus[["20260922"]],
  vsm_wet[["20260922"]],
  paw_wet[["20260922"]],
  dry_rsm_solus[["20260929"]],
  vsm_dry[["20260929"]],
  paw_dry[["20260929"]]
)

names(compare_days) <- c(
  "wet_rsm", "wet_vsm", "wet_paw",
  "dry_rsm", "dry_vsm", "dry_paw"
)

# inspect representative wet and dry conditions
global(compare_days, c("min", "mean", "max"), na.rm = TRUE)
plot(compare_days, nc = 3)

# compare wet and dry rsm with a common scale
plot(compare_days[[c("wet_rsm", "dry_rsm")]], range = c(0, 1))

# compare wet and dry localized vsm with a common scale
plot(compare_days[[c("wet_vsm", "dry_vsm")]],
     range = range(values(compare_days[[c("wet_vsm", "dry_vsm")]]), na.rm = TRUE))

# compare wet and dry paw with a common scale
plot(compare_days[[c("wet_paw", "dry_paw")]], range = c(0, 1))



# create unique ids for native sport cells ----------------------------------------------
sport_id <- wet_rsm[["20260922"]]
values(sport_id) <- seq_len(ncell(sport_id))

# crop and project sport cell ids to the solus grid
sport_id_solus <- sport_id %>%
  crop(vect(st_transform(levy, crs(sport_id)))) %>%
  mask(vect(st_transform(levy, crs(sport_id)))) %>%
  project(rosetta_levy_rast, method = "near")

names(sport_id_solus) <- "sport_id"

# summarize localized vsm within each sport cell
vsm_by_sport <- as.data.frame(
  c(sport_id_solus, wet_rsm_solus[["20260922"]], vsm_wet[["20260922"]]),
  na.rm = TRUE
) %>%
  rename(sport_id = 1, rsm = 2, vsm = 3) %>%
  group_by(sport_id) %>%
  summarise(
    n = n(),
    rsm = first(rsm),
    vsm_min = min(vsm),
    vsm_mean = mean(vsm),
    vsm_max = max(vsm),
    vsm_range = vsm_max - vsm_min,
    .groups = "drop"
  ) %>%
  arrange(desc(vsm_range))

print(vsm_by_sport, n=Inf)


#---------------------------------------------------------------------------------------
# isolate representative sport cell
sport_cell_id <- 1523378
sport_cell <- ifel(sport_id_solus == sport_cell_id, 1, NA)

# get extent of selected sport cell
cell_ext <- ext(as.polygons(sport_cell, na.rm = TRUE))

# extract and crop rsm and localized vsm
cell_compare <- c(
  mask(wet_rsm_solus[["20260922"]], sport_cell),
  mask(vsm_wet[["20260922"]], sport_cell)
) %>%
  crop(cell_ext)

names(cell_compare) <- c("rsm", "vsm_localized")

# inspect within-cell variation
global(cell_compare, c("min", "mean", "max"), na.rm = TRUE)
plot(cell_compare, nc = 2)

# extract and crop local hydraulic properties
cell_hydraulics <- rosetta_levy_rast %>%
  mask(sport_cell) %>%
  crop(cell_ext)

# inspect local hydraulic variation
global(cell_hydraulics, c("min", "mean", "max"), na.rm = TRUE)
plot(cell_hydraulics, nc = 2)


#-----------------------------------------------------------------------------------------
# summarize robust within-sport-cell vsm variation
vsm_cell_summary <- as.data.frame(
  c(sport_id_solus, wet_rsm_solus[["20260922"]], vsm_wet[["20260922"]]),
  na.rm = TRUE
) %>%
  rename(sport_id = 1, rsm = 2, vsm = 3) %>%
  group_by(sport_id) %>%
  summarise(
    n = n(),
    rsm = first(rsm),
    vsm_mean = mean(vsm),
    vsm_sd = sd(vsm),
    vsm_iqr = IQR(vsm),
    vsm_p10 = quantile(vsm, .10),
    vsm_p90 = quantile(vsm, .90),
    vsm_p90_p10 = vsm_p90 - vsm_p10,
    vsm_range = max(vsm) - min(vsm),
    .groups = "drop"
  ) %>%
  filter(n >= 800)

# summarize robust within-cell variation
vsm_cell_summary %>%
  summarise(
    n_cells = n(),
    median_iqr = median(vsm_iqr),
    median_p90_p10 = median(vsm_p90_p10),
    median_range = median(vsm_range),
    median_sd = median(vsm_sd)
  )

# inspect distribution of robust within-cell spread
quantile(vsm_cell_summary$vsm_p90_p10,
         probs = c(0, .1, .25, .5, .75, .9, .95, 1))

# plot within-cell vsm range distribution
ggplot(vsm_cell_summary, aes(vsm_range)) +
  geom_histogram(bins = 30) +
  labs(x = "Within-SPoRT-cell localized VSM range (m³/m³)",
       y = "Number of SPoRT cells") +
  theme_bw()

# examine whether within-cell variation depends on rsm
ggplot(vsm_cell_summary, aes(rsm, vsm_range)) +
  geom_point(alpha = 0.6) +
  labs(x = "SPoRT RSM",
       y = "Within-SPoRT-cell localized VSM range (m³/m³)") +
  theme_bw()



# where does localization add meaningful spatial detail ---------------------------------------
# map robust within-cell vsm spread back to sport cells
vsm_spread_rast <- subst(
  sport_id_solus,
  vsm_cell_summary$sport_id,
  vsm_cell_summary$vsm_p90_p10
)
names(vsm_spread_rast) <- "vsm_p90_p10"

# inspect and plot robust sub-grid variation
global(vsm_spread_rast, c("min", "mean", "max"), na.rm = TRUE)
plot(vsm_spread_rast)



# map robust within-cell vsm spread back to sport cells
vsm_spread_rast <- subst(
  sport_id_solus,
  vsm_cell_summary$sport_id,
  vsm_cell_summary$vsm_p90_p10,
  others = NA
)
names(vsm_spread_rast) <- "vsm_p90_p10"

# inspect robust sub-grid variation
global(vsm_spread_rast, c("min", "mean", "max"), na.rm = TRUE)
plot(vsm_spread_rast)


# --------------------------------------------------------------------------------------
# how much of the landscape gains min, moderate, sub-stantial sub-grid differentiation
# quantify proportion of cells for spread classes

# classify robust within-cell vsm spread
vsm_cell_summary <- vsm_cell_summary %>%
  mutate(
    spread_class = case_when(
      vsm_p90_p10 < 0.01 ~ "<0.01",
      vsm_p90_p10 < 0.03 ~ "0.01-0.03",
      vsm_p90_p10 < 0.06 ~ "0.03-0.06",
      TRUE ~ ">0.06"
    )
  )

# summarize proportion of sport cells by spread class
vsm_cell_summary %>%
  count(spread_class) %>%
  mutate(percent = 100 * n / sum(n))


# order spread classes for summaries and plots
vsm_cell_summary <- vsm_cell_summary %>%
  mutate(
    spread_class = factor(
      spread_class,
      levels = c("<0.01", "0.01-0.03", "0.03-0.06", ">0.06")
    )
  )

