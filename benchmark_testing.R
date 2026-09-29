# https://ncss-tech.github.io/rosettaPTF/

remotes::install_github("ncss-tech/rosettaPTF")

# load packages
pacman::p_load(here, tidyverse, terra, sf, tigris, rosettaPTF, soiltexture)

# find python binaries
rosettaPTF::find_python()

# install these once ------------------------------
# download latest python 3.10.x
# reticulate::install_python(version = "3.10:latest")
# reticulate::virtualenv_create("r-reticulate")
# install rosetta-soil Python Module
# rosettaPTF::install_rosetta()

# read data ============================================================================


# select SOLUS layers needed for initial Rosetta testing ----
(l = list.files("D:\\Projects\\fbrc\\data\\solus", full.names=TRUE, pattern = "\\.tif$"))
(solus_files <- l[str_detect(basename(l), "^(sandtotal|claytotal|dbovendry)_(0|5|15|30|60)_cm_p\\.tif$")] )
# read SOLUS layers as stack
solus <- rast(solus_files)

# mask rasters to FL---------------------------------------------------------
# get Florida boundary
fl <- states(cb = TRUE, resolution = "20m", year = 2025) %>%
  filter(STUSPS == "FL") %>%
  st_transform(crs(solus))

# crop and mask SOLUS to Florida
solus_fl <- crop(solus, vect(fl)) %>%
  mask(vect(fl))

# smaller test for levy county
# get Levy County boundary
levy <- counties(state = "FL", cb = TRUE, resolution = "20m", year = 2025) %>%
  filter(NAME == "Levy") %>%
  st_transform(crs(rast(l[1])))

# mask to levy
# crop and mask to Levy County
solus_levy <- crop(solus, vect(levy)) %>%
  mask(vect(levy))

# ========================================================================================
# testing on 0cm layer
### prepare Rosetta inputs ---------------------------------------------------------------

# prepare depth-specific Rosetta stacks
solus_fl_0 <- c(solus_fl[["sandtotal_0_cm_p"]],
                100 - solus_fl[["sandtotal_0_cm_p"]] - solus_fl[["claytotal_0_cm_p"]],
                solus_fl[["claytotal_0_cm_p"]], solus_fl[["dbovendry_0_cm_p"]] / 100)

names(solus_fl_0) <- c("sandtotal", "silttotal", "claytotal", "dbovendry")


# prepare 0-cm Levy County Rosetta inputs
solus_levy_0 <- c(
  solus_levy[["sandtotal_0_cm_p"]],
  100 - solus_levy[["sandtotal_0_cm_p"]] - solus_levy[["claytotal_0_cm_p"]],
  solus_levy[["claytotal_0_cm_p"]],
  solus_levy[["dbovendry_0_cm_p"]] / 100
)

names(solus_levy_0) <- c("sandtotal", "silttotal", "claytotal", "dbovendry")




# convert rast to data frame -------------------------------------------------------------
# extract Florida Rosetta inputs and retain unique combinations
# solus_fl_0_unique <- as.data.frame(solus_fl_0, na.rm = TRUE) %>%
#   distinct()





# extract unique Levy County Rosetta input combinations
solus_levy_0_unique <- as.data.frame(solus_levy_0, na.rm = TRUE) %>%
  distinct()


# run Rosetta for each unique soil combination
rosetta_levy_0_unique <- run_rosetta(
  solus_levy_0_unique,
  vars = c("sandtotal", "silttotal", "claytotal", "dbovendry")
)

# classify unique combinations using USDA soil texture classes
usda_texture <- soiltexture::TT.points.in.classes(
  tri.data = solus_levy_0_unique %>%
    transmute(CLAY = claytotal, SILT = silttotal, SAND = sandtotal),
  class.sys = "USDA.TT",
  PiC.type = "t"
)

# combine soil inputs, texture class, and Rosetta predictions
rosetta_levy_0_lookup <- bind_cols(
  solus_levy_0_unique %>% mutate(texture = usda_texture),
  rosetta_levy_0_unique
)



# run Rosetta once for each unique Florida soil combination
system.time(
  rosetta_levy_0_unique <- run_rosetta(
    solus_levy_0_unique,
    vars = c("sandtotal", "silttotal", "claytotal", "dbovendry")
  )
)





# classify unique combinations using USDA soil texture classes
usda_texture <- soiltexture::TT.points.in.classes(
  tri.data = solus_levy_0_unique %>%
    transmute(CLAY = claytotal, SILT = silttotal, SAND = sandtotal),
  class.sys = "USDA.TT", PiC.type = "t"
)

# combine inputs and predictions and calculate AWC at both field-capacity definitions
rosetta_levy_0_lookup <- bind_cols(
  solus_levy_0_unique %>% mutate(texture = usda_texture),
  rosetta_levy_0_unique
) %>%
  mutate(
    alpha = 10^log10_alpha_mean,
    npar = 10^log10_npar_mean,
    m = 1 - 1 / npar,
    theta_10 = theta_r_mean + (theta_s_mean - theta_r_mean) /
      (1 + (alpha * (10 * 10.197))^npar)^m,
    theta_33 = theta_r_mean + (theta_s_mean - theta_r_mean) /
      (1 + (alpha * (33 * 10.197))^npar)^m,
    theta_1500 = theta_r_mean + (theta_s_mean - theta_r_mean) /
      (1 + (alpha * (1500 * 10.197))^npar)^m,
    awc_10 = (theta_10 - theta_1500) * 1000,
    awc_33 = (theta_33 - theta_1500) * 1000,
    awc_diff = awc_10 - awc_33
  )

# compare AWC definitions among USDA texture classes
rosetta_levy_0_lookup %>%
  group_by(texture) %>%
  summarise(
    n = n(),
    sand = mean(sandtotal),
    awc_10 = mean(awc_10),
    awc_33 = mean(awc_33),
    difference = mean(awc_diff),
    pct_difference = 100 * mean(awc_diff / awc_33),
    .groups = "drop"
  ) %>%
  arrange(desc(sand))




# =======================================================================================
# calculate VG-Mualem conductivity at pressure head h (cm)
vg_k <- function(h, alpha, npar, K0, L) {
  m <- 1 - 1 / npar
  Se <- (1 + (alpha * h)^npar)^(-m)
  K0 * Se^L * (1 - (1 - Se^(1 / m))^m)^2
}

# solve for pressure head where conductivity equals drainage threshold
find_h_fc <- function(alpha, npar, K0, L, q_fc) {
  f <- \(log_h) vg_k(10^log_h, alpha, npar, K0, L) - q_fc
  10^uniroot(f, interval = c(-3, 7))$root
}


# calculate flux-defined FC for each unique Levy soil
rosetta_levy_0_flux <- rosetta_levy_0_lookup %>%
  mutate(
    K0 = 10^log10_K0_mean,
    h_fc_05 = pmap_dbl(list(alpha, npar, K0, lpar_mean),
                       \(a, n, k, l) find_h_fc(a, n, k, l, 0.05)),
    h_fc_1 = pmap_dbl(list(alpha, npar, K0, lpar_mean),
                      \(a, n, k, l) find_h_fc(a, n, k, l, 0.10)),
    theta_fc_05 = theta_r_mean + (theta_s_mean - theta_r_mean) /
      (1 + (alpha * h_fc_05)^npar)^m,
    theta_fc_1 = theta_r_mean + (theta_s_mean - theta_r_mean) /
      (1 + (alpha * h_fc_1)^npar)^m,
    awc_flux_05 = (theta_fc_05 - theta_1500) * 1000,
    awc_flux_1 = (theta_fc_1 - theta_1500) * 1000
  )




# compare flux-defined FC and AWC among USDA texture classes
rosetta_levy_0_flux %>%
  mutate(psi_fc_05 = -h_fc_05 / 10.197,
         psi_fc_1 = -h_fc_1 / 10.197) %>%
  group_by(texture) %>%
  summarise(
    n = n(),
    sand = mean(sandtotal),
    psi_fc_05 = mean(psi_fc_05),
    psi_fc_1 = mean(psi_fc_1),
    awc_flux_05 = mean(awc_flux_05),
    awc_flux_1 = mean(awc_flux_1),
    awc_10 = mean(awc_10),
    awc_33 = mean(awc_33),
    .groups = "drop"
  ) %>%
  arrange(desc(sand))



# calculate Assouline-Or field capacity and AWC
rosetta_levy_0_flux <- rosetta_levy_0_flux %>%
  mutate(
    theta_fc_ao = theta_r_mean + (theta_s_mean - theta_r_mean) *
      (1 + (1 - 1 / npar) * (1 - 2 * npar))^(-1 + 1 / npar),
    awc_ao = (theta_fc_ao - theta_1500) * 1000
  )


# calculate equivalent pressure head for Assouline-Or field capacity
rosetta_levy_0_flux <- rosetta_levy_0_flux %>%
  mutate(
    se_fc_ao = (theta_fc_ao - theta_r_mean) / (theta_s_mean - theta_r_mean),
    h_fc_ao = ((se_fc_ao^(-1 / m) - 1)^(1 / npar)) / alpha,
    psi_fc_ao = -h_fc_ao / 10.197
  )

# compare field-capacity definitions by USDA texture
rosetta_levy_0_flux %>%
  mutate(psi_fc_05 = -h_fc_05 / 10.197,
         psi_fc_1 = -h_fc_1 / 10.197) %>%
  group_by(texture) %>%
  summarise(
    n = n(),
    sand = mean(sandtotal),
    psi_ao = mean(psi_fc_ao),
    psi_flux_05 = mean(psi_fc_05),
    psi_flux_1 = mean(psi_fc_1),
    awc_ao = mean(awc_ao),
    awc_flux_05 = mean(awc_flux_05),
    awc_flux_1 = mean(awc_flux_1),
    awc_10 = mean(awc_10),
    awc_33 = mean(awc_33),
    .groups = "drop"
  ) %>%
  arrange(desc(sand))



# check flux-derived field capacity against physical water-content bounds
rosetta_levy_0_flux %>%
  summarise(
    n = n(),
    fc05_below_pwp = sum(theta_fc_05 <= theta_1500),
    fc05_above_sat = sum(theta_fc_05 >= theta_s_mean),
    fc1_below_pwp = sum(theta_fc_1 <= theta_1500),
    fc1_above_sat = sum(theta_fc_1 >= theta_s_mean)
  )

# inspect ranges of flux-derived FC and AWC
rosetta_levy_0_flux %>%
  summarise(across(
    c(theta_fc_05, theta_fc_1, awc_flux_05, awc_flux_1),
    list(min = min, median = median, max = max)
  ))



# verify solved FC heads reproduce the specified drainage fluxes
rosetta_levy_0_flux %>%
  summarise(
    max_error_05 = max(abs(
      pmap_dbl(list(h_fc_05, alpha, npar, K0, lpar_mean),
               \(h, a, n, k, l) vg_k(h, a, n, k, l)) - 0.05
    )),
    max_error_1 = max(abs(
      pmap_dbl(list(h_fc_1, alpha, npar, K0, lpar_mean),
               \(h, a, n, k, l) vg_k(h, a, n, k, l)) - 0.10
    ))
  )



# ========================================================================================
# now scale to all depth layers 

# calculate Rosetta hydraulic properties and FC for one SOLUS depth
run_depth <- function(depth) {
  x <- c(
    solus_levy[[paste0("sandtotal_", depth, "_cm_p")]],
    100 - solus_levy[[paste0("sandtotal_", depth, "_cm_p")]] -
      solus_levy[[paste0("claytotal_", depth, "_cm_p")]],
    solus_levy[[paste0("claytotal_", depth, "_cm_p")]],
    solus_levy[[paste0("dbovendry_", depth, "_cm_p")]] / 100
  )
  names(x) <- c("sandtotal", "silttotal", "claytotal", "dbovendry")
  
  # create unique soil combinations and integer lookup ID
  unique_x <- as.data.frame(x, na.rm = TRUE) %>%
    distinct() %>%
    mutate(soil_id = row_number())
  
  # estimate hydraulic properties and field capacity
  bind_cols(
    unique_x,
    run_rosetta(unique_x %>% select(-soil_id),
                vars = c("sandtotal", "silttotal", "claytotal", "dbovendry"))
  ) %>%
    mutate(
      depth = depth,
      alpha = 10^log10_alpha_mean,
      npar = 10^log10_npar_mean,
      m = 1 - 1 / npar,
      K0 = 10^log10_K0_mean,
      h_fc_05 = pmap_dbl(list(alpha, npar, K0, lpar_mean),
                         \(a, n, k, l) find_h_fc(a, n, k, l, 0.05)),
      h_fc_1 = pmap_dbl(list(alpha, npar, K0, lpar_mean),
                        \(a, n, k, l) find_h_fc(a, n, k, l, 0.10)),
      theta_fc_05 = theta_r_mean + (theta_s_mean - theta_r_mean) /
        (1 + (alpha * h_fc_05)^npar)^m,
      theta_fc_1 = theta_r_mean + (theta_s_mean - theta_r_mean) /
        (1 + (alpha * h_fc_1)^npar)^m,
      theta_1500 = theta_r_mean + (theta_s_mean - theta_r_mean) /
        (1 + (alpha * (1500 * 10.197))^npar)^m,
      awc_flux_05 = (theta_fc_05 - theta_1500) * 1000,
      awc_flux_1 = (theta_fc_1 - theta_1500) * 1000
    )
}

# run Rosetta for each SOLUS depth
rosetta_levy <- c(0, 5, 15, 30, 60) %>%
  set_names() %>%
  map(run_depth)


# summarize hydraulic properties across depth
rosetta_levy %>%
  list_rbind() %>%
  group_by(depth) %>%
  summarise(
    n = n(),
    psi_fc_05 = median(-h_fc_05 / 10.197),
    psi_fc_1 = median(-h_fc_1 / 10.197),
    awc_05 = median(awc_flux_05),
    awc_1 = median(awc_flux_1),
    .groups = "drop"
  )


# map Rosetta lookup values back to the SOLUS raster grid
rasterize_depth <- function(depth) {
  x <- c(
    solus_levy[[paste0("sandtotal_", depth, "_cm_p")]],
    100 - solus_levy[[paste0("sandtotal_", depth, "_cm_p")]] -
      solus_levy[[paste0("claytotal_", depth, "_cm_p")]],
    solus_levy[[paste0("claytotal_", depth, "_cm_p")]],
    solus_levy[[paste0("dbovendry_", depth, "_cm_p")]] / 100
  )
  names(x) <- c("sandtotal", "silttotal", "claytotal", "dbovendry")
  
  # match raster-cell combinations to the Rosetta lookup
  lookup <- rosetta_levy[[as.character(depth)]]
  cells <- as.data.frame(x, cells = TRUE, na.rm = TRUE) %>%
    left_join(lookup %>%
                select(sandtotal, silttotal, claytotal, dbovendry, soil_id),
              by = c("sandtotal", "silttotal", "claytotal", "dbovendry"))
  
  # create soil-ID raster
  soil_id <- x[[1]]
  values(soil_id) <- NA
  soil_id[cells$cell] <- cells$soil_id
  
  # substitute hydraulic properties onto the grid
  c(
    subst(soil_id, lookup$soil_id, lookup$theta_fc_05),
    subst(soil_id, lookup$soil_id, lookup$theta_fc_1),
    subst(soil_id, lookup$soil_id, lookup$theta_1500)
  ) %>%
    setNames(c(
      paste0("theta_fc_05_", depth),
      paste0("theta_fc_1_", depth),
      paste0("theta_1500_", depth)
    ))
}


# map hydraulic properties at all SOLUS depths
hydraulic_levy <- c(0, 5, 15, 30, 60) %>%
  map(rasterize_depth) %>%
  rast()



# inspect hydraulic raster
hydraulic_levy
plot(hydraulic_levy)


# calculate volumetric AWC at each SOLUS depth
awc_05 <- hydraulic_levy[[paste0("theta_fc_05_", c(0, 5, 15, 30, 60))]] -
  hydraulic_levy[[paste0("theta_1500_", c(0, 5, 15, 30, 60))]]

awc_1 <- hydraulic_levy[[paste0("theta_fc_1_", c(0, 5, 15, 30, 60))]] -
  hydraulic_levy[[paste0("theta_1500_", c(0, 5, 15, 30, 60))]]



# calculate TAW to any rooting depth from 0-60 cm
calc_taw <- function(awc, drz) {
  depths <- c(0, 5, 15, 30, 60)
  stopifnot(drz > 0, drz <= max(depths))
  
  full <- which(depths[-1] <= drz)
  
  if (length(full)) {
    taw_layers <- (awc[[full]] + awc[[full + 1]]) / 2 *
      (depths[full + 1] - depths[full]) * 10
    taw <- sum(taw_layers)
  } else {
    taw <- awc[[1]] * 0
  }
  
  if (drz < max(depths)) {
    i <- findInterval(drz, depths)
    awc_tip <- awc[[i]] + (awc[[i + 1]] - awc[[i]]) *
      (drz - depths[i]) / (depths[i + 1] - depths[i])
    taw <- taw + (awc[[i]] + awc_tip) / 2 * (drz - depths[i]) * 10
  }
  
  taw
}


# calculate TAW at selected rooting depths
taw_levy <- c(
  calc_taw(awc_05, 5),
  calc_taw(awc_05, 15),
  calc_taw(awc_05, 30),
  calc_taw(awc_05, 60),
  calc_taw(awc_1, 5),
  calc_taw(awc_1, 15),
  calc_taw(awc_1, 30),
  calc_taw(awc_1, 60)
)

names(taw_levy) <- c(
  "taw_05_5cm", "taw_05_15cm", "taw_05_30cm", "taw_05_60cm",
  "taw_1_5cm", "taw_1_15cm", "taw_1_30cm", "taw_1_60cm"
)


# inspect root-zone available water
taw_levy
global(taw_levy, c("min", "mean", "max"), na.rm = TRUE)
plot(taw_levy)


# now test TAW at arbitrary rooting depths
taw_levy_test <- c(
  calc_taw(awc_05, 10),
  calc_taw(awc_05, 20),
  calc_taw(awc_05, 40),
  calc_taw(awc_1, 10),
  calc_taw(awc_1, 20),
  calc_taw(awc_1, 40)
)

names(taw_levy_test) <- c(
  "taw_05_10cm", "taw_05_20cm", "taw_05_40cm",
  "taw_1_10cm", "taw_1_20cm", "taw_1_40cm"
)

# check arbitrary-depth TAW
global(taw_levy_test, c("min", "mean", "max"), na.rm = TRUE)

# summarize mean TAW by rooting depth and FC threshold
global(taw_levy, "mean", na.rm = TRUE) %>%
  rownames_to_column("variable") %>%
  mutate(
    fc = if_else(str_detect(variable, "^taw_05"), "0.5 mm/day", "1.0 mm/day"),
    root_depth = parse_number(str_extract(variable, "\\d+cm"))
  ) %>%
  select(root_depth, fc, mean) %>%
  pivot_wider(names_from = fc, values_from = mean) %>%
  arrange(root_depth)


# The next step is not more raster processing. 
# Before deciding exactly what soil quantity to feed CWB, 
# we should inspect CropWaterBalance itself—specifically how it uses AWC and Drz, 
# its units, whether AWC can vary through time, 
# and whether CWB internally calculates something equivalent to:
#   TAW=AWC×D











# raster based approach ==================================================================
# # corrected SpatRaster wrapper for rosettaPTF -------------------------------------------
# run_rosetta_rast <- function(soildata, vars = NULL, rosetta_version = 3,
#                              estimate_type = "log", cores = 1,
#                              core_thresh = 20000L, file = paste0(tempfile(), ".tif"),
#                              nrows = 100L, overwrite = TRUE) {
#   
#   if (any(!terra::inMemory(soildata))) {
#     terra::readStart(soildata)
#     on.exit(try(terra::readStop(soildata), silent = TRUE), add = TRUE)
#   }
#   
#   # initialize Rosetta output raster
#   out <- terra::rast(soildata)
#   sample_res <- rosettaPTF:::run_rosetta.default(
#     list(c(33, 33, 34)), rosetta_version = rosetta_version,
#     estimate_type = estimate_type
#   )
#   cnm <- colnames(sample_res)
#   terra::nlyr(out) <- length(cnm)
#   names(out) <- cnm
#   
#   # initialize output file
#   terra::writeStart(out, filename = file, overwrite = overwrite)
#   on.exit(try(out <- terra::writeStop(out), silent = TRUE), add = TRUE)
#   
#   # define processing blocks from input raster dimensions
#   start_row <- seq(1L, nrow(soildata), by = nrows)
#   n_row <- pmin(nrows, nrow(soildata) - start_row + 1L)
#   
#   # process raster blocks in parallel
#   if (cores > 1 && ncell(soildata) > core_thresh) {
#     cls <- parallel::makeCluster(cores)
#     on.exit(parallel::stopCluster(cls), add = TRUE)
#     
#     for (i in seq_along(start_row)) {
#       blockdata <- terra::readValues(
#         soildata, row = start_row[i], nrows = n_row[i], dataframe = TRUE
#       )
#       
#       # n <- max(cores, ceiling(nrow(blockdata) / core_thresh))
#       # X <- split(blockdata, rep(seq_len(n), length.out = nrow(blockdata)))
#       # split block into contiguous chunks while preserving cell order
#       n <- min(nrow(blockdata), max(cores, ceiling(nrow(blockdata) / core_thresh)))
#       X <- split(blockdata, cut(seq_len(nrow(blockdata)), breaks = n, labels = FALSE))
#       
#       r <- do.call("rbind", parallel::clusterApply(
#         cls, X, \(x) rosettaPTF::run_rosetta(
#           x, vars = vars, rosetta_version = rosetta_version,
#           estimate_type = estimate_type
#         )
#       ))
#       
#       terra::writeValues(out, as.matrix(r), start_row[i], nrows = n_row[i])
#     }
#     
#     # process raster blocks sequentially
#   } else {
#     for (i in seq_along(start_row)) {
#       r <- rosettaPTF::run_rosetta(
#         terra::readValues(
#           soildata, row = start_row[i], nrows = n_row[i], dataframe = TRUE
#         ),
#         vars = vars, rosetta_version = rosetta_version,
#         estimate_type = estimate_type
#       )
#       
#       terra::writeValues(out, as.matrix(r), start_row[i], nrows = n_row[i])
#     }
#   }
#   
#   out <- terra::writeStop(out)
#   out
# }
# 
# 
# # test corrected raster method on levy county
# system.time(
#   rosetta_levy_0 <- run_rosetta_rast(solus_levy_0, cores = 1, nrows = 100)
# )
# 
# # ========================================================================================
# ### benchmark Rosetta cores --------------------------------------------------------------
# 
# bench_cores <- c(4, 8, 16, 22)
# 
# bench_time <- map_dfr(bench_cores, \(n_cores) {
#   tm <- system.time(
#     run_rosetta_rast(solus_levy_0, cores = n_cores, nrows = 100)
#   )
#   
#   tibble(cores = n_cores, elapsed_sec = unname(tm["elapsed"]))
# })
# 
# # add existing single-core result
# bench_time <- bind_rows(
#   tibble(cores = 1, elapsed_sec = 116.71),
#   bench_time
# )
# 
# bench_time
# 
# # benchmark raster block size with four cores
# bench_nrows <- c(50, 100, 250, 500)
# 
# bench_blocks <- map_dfr(bench_nrows, \(n_rows) {
#   tm <- system.time(
#     run_rosetta_rast(solus_levy_0, cores = 4, nrows = n_rows)
#   )
#   
#   tibble(nrows = n_rows, elapsed_sec = unname(tm["elapsed"]))
# })
# 
# bench_blocks
# 
# 
# # run optimized Levy County Rosetta test
# rosetta_levy_0_4core <- run_rosetta_rast(
#   solus_levy_0, cores = 4, nrows = 500
# )
# 
# # compare single- and four-core predictions
# compareGeom(rosetta_levy_0, rosetta_levy_0_4core)
# global(abs(rosetta_levy_0 - rosetta_levy_0_4core), "max", na.rm = TRUE)
# 
# 
# 
# # run 0-cm Rosetta model for Florida rasters 
# system.time(
#   rosetta_fl_0 <- run_rosetta_rast(
#     solus_fl_0, cores = 4, nrows = 500
#   )
# )











