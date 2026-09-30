# ==============================================================================
# TITLE
# Association between mosquito-net use and childhood anaemia in India:
# District-level Bayesian spatial analysis using NFHS-5 (2019-21)
#
# SCRIPT:
# 01_bayesian_spatial_anemia_analysis.R
#
# PURPOSE:
# Reproducible district-level spatial analysis including:
#   1. Data validation and preparation
#   2. Global Moran's I
#   3. Local Indicators of Spatial Association (LISA)
#   4. Non-spatial beta regression
#   5. Bayesian spatial beta regression using BYM2-INLA
#   6. Model comparison using DIC and WAIC
#   7. Posterior fixed effects
#   8. BYM2 hyperparameters
#   9. Posterior spatial effects
#  10. Posterior exceedance probabilities
#  11. Reproducibility/session information
#
# IMPORTANT:
# This script assumes that:
#   - `df` contains the district-level analytical dataset.
#   - `districts` is an sf object containing district polygons.
#   - `net_use` has already been constructed from the NFHS-5 raw variables.
#
# Required net_use categories:
#   no_net
#   untreated_net
#   itn
#
# No absolute/local file paths are used in this script.
# ==============================================================================


# ==============================================================================
# 1. PACKAGES
# ==============================================================================

required_packages <- c(
  "sf",
  "spdep",
  "INLA",
  "dplyr",
  "ggplot2",
  "tmap"
)

missing_packages <- setdiff(
  required_packages,
  rownames(installed.packages())
)

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "The following packages are missing: ",
      paste(missing_packages, collapse = ", "),
      "\nPlease install them before running this script."
    )
  )
}

library(sf)
library(spdep)
library(INLA)
library(dplyr)
library(ggplot2)
library(tmap)


# ==============================================================================
# 2. REPRODUCIBILITY SETTINGS
# ==============================================================================

set.seed(12345)

options(
  scipen = 999,
  stringsAsFactors = FALSE
)


# ==============================================================================
# 3. OUTPUT DIRECTORIES
# ==============================================================================

dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/tables", showWarnings = FALSE)
dir.create("outputs/figures", showWarnings = FALSE)
dir.create("outputs/spatial", showWarnings = FALSE)


# ==============================================================================
# 4. REQUIRED INPUT OBJECTS
# ==============================================================================

# ----------------------------------------------------------------------
# IMPORTANT:
#
# The following two objects must already exist:
#
#   df         = district-level analytical dataset
#   districts  = district polygon sf object
#
# The data preparation/import stage should create these objects.
# No machine-specific paths are included in this GitHub script.
# ----------------------------------------------------------------------

if (!exists("df")) {
  stop(
    "Object `df` was not found. Please load the district-level analytical dataset first."
  )
}

if (!exists("districts")) {
  stop(
    "Object `districts` was not found. Please load the district polygon sf object first."
  )
}

if (!inherits(districts, "sf")) {
  stop("`districts` must be an sf spatial object.")
}


# ==============================================================================
# 5. REQUIRED VARIABLES
# ==============================================================================

required_variables <- c(
  "dist_code",
  "Anemic",
  "net_use",
  "male",
  "LBW",
  "child_12_17_months",
  "single_birth",
  "Stunted",
  "bord_1",
  "size_large",
  "no_edu",
  "age_30plus",
  "poorest",
  "SC",
  "rural",
  "Hindu"
)

missing_variables <- setdiff(
  required_variables,
  names(df)
)

if (length(missing_variables) > 0) {
  stop(
    paste0(
      "The following required variables are missing from `df`: ",
      paste(missing_variables, collapse = ", ")
    )
  )
}

if (!"REGCODE" %in% names(districts)) {
  stop("`districts` must contain the district identifier `REGCODE`.")
}


# ==============================================================================
# 6. CLEAN DISTRICT IDENTIFIERS
# ==============================================================================

df <- df %>%
  mutate(
    dist_code = as.character(dist_code)
  )

districts <- districts %>%
  mutate(
    REGCODE = as.character(REGCODE)
  )


# ==============================================================================
# 7. PREPARE NET-USE EXPOSURE
# ==============================================================================

# ----------------------------------------------------------------------
# The reviewer specifically requested three categories:
#
#   1. no_net
#   2. untreated_net
#   3. itn
#
# The script does not guess the original NFHS numeric coding.
# The variable must therefore already be constructed from the
# underlying NFHS variables.
# ----------------------------------------------------------------------

df <- df %>%
  mutate(
    net_use = factor(
      net_use,
      levels = c(
        "no_net",
        "untreated_net",
        "itn"
      )
    )
  )

if (any(is.na(df$net_use))) {
  stop(
    paste0(
      "`net_use` contains missing or unrecognised categories.\n",
      "Required categories are: no_net, untreated_net, itn."
    )
  )
}

# Explicit reference category
df$net_use <- relevel(
  df$net_use,
  ref = "no_net"
)

# Verify exposure
net_use_check <- as.data.frame(
  table(
    df$net_use,
    useNA = "ifany"
  )
)

names(net_use_check) <- c(
  "Net_use_category",
  "N"
)

write.csv(
  net_use_check,
  "outputs/tables/Net_use_distribution.csv",
  row.names = FALSE
)


# ==============================================================================
# 8. PREPARE NUMERIC VARIABLES
# ==============================================================================

numeric_variables <- c(
  "Anemic",
  "male",
  "LBW",
  "child_12_17_months",
  "single_birth",
  "Stunted",
  "bord_1",
  "size_large",
  "no_edu",
  "age_30plus",
  "poorest",
  "SC",
  "rural",
  "Hindu"
)

df[numeric_variables] <- lapply(
  df[numeric_variables],
  function(x) as.numeric(as.character(x))
)


# ==============================================================================
# 9. CHECK OUTCOME SCALE
# ==============================================================================

if (any(is.na(df$Anemic))) {
  stop("The outcome variable `Anemic` contains missing values.")
}

# Convert percentage to proportion if necessary
if (max(df$Anemic, na.rm = TRUE) > 1) {
  
  df$Anemic <- df$Anemic / 100
  
}

if (min(df$Anemic, na.rm = TRUE) < 0 ||
    max(df$Anemic, na.rm = TRUE) > 1) {
  
  stop(
    "The outcome `Anemic` must be between 0 and 1."
  )
  
}


# ==============================================================================
# 10. BETA-REGRESSION BOUNDARY ADJUSTMENT
# ==============================================================================

# Beta likelihood requires observations strictly between 0 and 1.
# The transformation is applied ONLY to the outcome.

if (any(
  df$Anemic <= 0 |
  df$Anemic >= 1,
  na.rm = TRUE
)) {
  
  n <- nrow(df)
  
  df$Anemic <- (
    df$Anemic * (n - 1) + 0.5
  ) / n
  
}


# ==============================================================================
# 11. MERGE DISTRICT DATA WITH POLYGONS
# ==============================================================================

districts <- districts %>%
  left_join(
    df,
    by = c(
      "REGCODE" = "dist_code"
    )
  ) %>%
  filter(
    !is.na(Anemic)
  )

if (nrow(districts) == 0) {
  stop(
    "No observations remained after joining the analytical data to the district polygons."
  )
}


# ==============================================================================
# 12. CREATE DISTRICT INDEX
# ==============================================================================

districts <- districts %>%
  mutate(
    district_index = seq_len(n())
  )


# ==============================================================================
# 13. CHECK ANALYTICAL DATA
# ==============================================================================

model_variables <- c(
  "Anemic",
  "net_use",
  "male",
  "LBW",
  "child_12_17_months",
  "single_birth",
  "Stunted",
  "bord_1",
  "size_large",
  "no_edu",
  "age_30plus",
  "poorest",
  "SC",
  "rural",
  "Hindu"
)

missing_model_data <- colSums(
  is.na(
    st_drop_geometry(districts)[model_variables]
  )
)

missing_model_data <- data.frame(
  Variable = names(missing_model_data),
  Missing = as.numeric(missing_model_data)
)

write.csv(
  missing_model_data,
  "outputs/tables/Missing_data_check.csv",
  row.names = FALSE
)


# ==============================================================================
# 14. DESCRIPTIVE DISTRICT-LEVEL SUMMARY
# ==============================================================================

descriptive_summary <- data.frame(
  Statistic = c(
    "Number of districts",
    "Mean anaemia prevalence",
    "SD anaemia prevalence",
    "Minimum anaemia prevalence",
    "Maximum anaemia prevalence"
  ),
  Value = c(
    nrow(districts),
    mean(districts$Anemic, na.rm = TRUE),
    sd(districts$Anemic, na.rm = TRUE),
    min(districts$Anemic, na.rm = TRUE),
    max(districts$Anemic, na.rm = TRUE)
  )
)

write.csv(
  descriptive_summary,
  "outputs/tables/Descriptive_summary.csv",
  row.names = FALSE
)


# ==============================================================================
# 15. SPATIAL NEIGHBOURHOOD STRUCTURE
# ==============================================================================

# First-order queen contiguity
nb <- poly2nb(
  districts,
  queen = TRUE,
  snap = 0.001
)

# Spatial weights
lw <- nb2listw(
  nb,
  style = "W",
  zero.policy = TRUE
)

# Neighbour summary
neighbour_summary <- data.frame(
  District = seq_along(nb),
  Number_of_neighbours = card(nb)
)

write.csv(
  neighbour_summary,
  "outputs/tables/District_neighbour_summary.csv",
  row.names = FALSE
)

# Identify disconnected districts
isolated_districts <- which(
  card(nb) == 0
)

if (length(isolated_districts) > 0) {
  
  warning(
    paste(
      length(isolated_districts),
      "district(s) have no spatial neighbours."
    )
  )
  
}


# ==============================================================================
# 16. CREATE INLA GRAPH
# ==============================================================================

graph_file <- file.path(
  tempdir(),
  "india_district.graph"
)

nb2INLA(
  graph_file,
  nb
)

g <- inla.read.graph(
  graph_file
)


# ==============================================================================
# 17. GLOBAL MORAN'S I
# ==============================================================================

moran_result <- moran.test(
  districts$Anemic,
  lw,
  zero.policy = TRUE
)

moran_table <- data.frame(
  Moran_I = unname(moran_result$estimate[1]),
  Expected_I = unname(moran_result$estimate[2]),
  Variance = unname(moran_result$estimate[3]),
  P_value = moran_result$p.value
)

write.csv(
  moran_table,
  "outputs/tables/Global_Morans_I.csv",
  row.names = FALSE
)


# ==============================================================================
# 18. LOCAL INDICATORS OF SPATIAL ASSOCIATION (LISA)
# ==============================================================================

lisa_result <- localmoran(
  districts$Anemic,
  lw,
  zero.policy = TRUE
)

lisa_table <- data.frame(
  Local_I = lisa_result[, "Ii"],
  Expected_I = lisa_result[, "E.Ii"],
  Variance = lisa_result[, "Var.Ii"],
  Z_score = lisa_result[, "Z.Ii"],
  P_value = lisa_result[, "Pr(z != E(Ii))"]
)

districts <- bind_cols(
  districts,
  lisa_table
)

# LISA cluster classification
anaemia_mean <- mean(
  districts$Anemic,
  na.rm = TRUE
)

lag_anaemia <- lag.listw(
  lw,
  districts$Anemic,
  zero.policy = TRUE
)

districts$LISA_cluster <- "Not significant"

districts$LISA_cluster[
  districts$Anemic > anaemia_mean &
    lag_anaemia > anaemia_mean &
    districts$P_value < 0.05
] <- "High-High"

districts$LISA_cluster[
  districts$Anemic < anaemia_mean &
    lag_anaemia < anaemia_mean &
    districts$P_value < 0.05
] <- "Low-Low"

districts$LISA_cluster[
  districts$Anemic > anaemia_mean &
    lag_anaemia < anaemia_mean &
    districts$P_value < 0.05
] <- "High-Low"

districts$LISA_cluster[
  districts$Anemic < anaemia_mean &
    lag_anaemia > anaemia_mean &
    districts$P_value < 0.05
] <- "Low-High"

write.csv(
  st_drop_geometry(districts) %>%
    select(
      REGCODE,
      Anemic,
      Local_I,
      Z_score,
      P_value,
      LISA_cluster
    ),
  "outputs/tables/LISA_results.csv",
  row.names = FALSE
)


# ==============================================================================
# 19. NON-SPATIAL BETA REGRESSION
# ==============================================================================

formula_baseline <- Anemic ~
  net_use +
  male +
  LBW +
  child_12_17_months +
  single_birth +
  Stunted +
  bord_1 +
  size_large +
  no_edu +
  age_30plus +
  poorest +
  SC +
  rural +
  Hindu


model_baseline <- inla(
  formula_baseline,
  family = "beta",
  data = districts,
  
  control.predictor = list(
    compute = TRUE
  ),
  
  control.compute = list(
    dic = TRUE,
    waic = TRUE,
    cpo = TRUE
  )
)


# ==============================================================================
# 20. BYM2 PRIOR SPECIFICATION
# ==============================================================================

# IMPORTANT:
# These are the documented R-INLA BYM2 defaults.
#
# Precision:
#   P(sigma > 1) = 0.01
#
# Spatial mixing:
#   P(phi <= 0.5) = 0.5

bym2_hyper <- list(
  
  prec = list(
    prior = "pc.prec",
    param = c(1, 0.01)
  ),
  
  phi = list(
    prior = "pc",
    param = c(0.5, 0.5)
  )
  
)


# ==============================================================================
# 21. BAYESIAN SPATIAL BYM2 MODEL
# ==============================================================================

formula_bym2 <- Anemic ~
  net_use +
  male +
  LBW +
  child_12_17_months +
  single_birth +
  Stunted +
  bord_1 +
  size_large +
  no_edu +
  age_30plus +
  poorest +
  SC +
  rural +
  Hindu +
  f(
    district_index,
    model = "bym2",
    graph = g,
    hyper = bym2_hyper
  )


# ------------------------------------------------------------------------------
# INLA computational settings
#
# simplified.laplace and auto are used explicitly here so that the
# computational specification is transparent in the GitHub version.
#
# NOTE:
# The originally reported model did not manually specify control.inla.
# Therefore, if reproducing the ORIGINAL published estimates exactly,
# remove this control.inla block and rely on the INLA defaults.
#
# For the revised reviewer model, the settings are made explicit.
# ------------------------------------------------------------------------------

model_bym2 <- inla(
  formula_bym2,
  family = "beta",
  data = districts,
  
  control.predictor = list(
    compute = TRUE
  ),
  
  control.compute = list(
    dic = TRUE,
    waic = TRUE,
    cpo = TRUE
  ),
  
  control.inla = list(
    strategy = "simplified.laplace",
    int.strategy = "auto"
  )
)


# ==============================================================================
# 22. MODEL COMPARISON
# ==============================================================================

model_comparison <- data.frame(
  
  Model = c(
    "Non-spatial Beta Regression",
    "Bayesian Spatial BYM2 Beta Regression"
  ),
  
  DIC = c(
    model_baseline$dic$dic,
    model_bym2$dic$dic
  ),
  
  Effective_parameters_DIC = c(
    model_baseline$dic$p.eff,
    model_bym2$dic$p.eff
  ),
  
  WAIC = c(
    model_baseline$waic$waic,
    model_bym2$waic$waic
  ),
  
  Effective_parameters_WAIC = c(
    model_baseline$waic$p.eff,
    model_bym2$waic$p.eff
  )
)

write.csv(
  model_comparison,
  "outputs/tables/Model_Comparison.csv",
  row.names = FALSE
)


# ==============================================================================
# 23. BAYESIAN FIXED-EFFECT RESULTS
# ==============================================================================

fixed_effects <- model_bym2$summary.fixed

fixed_effect_table <- data.frame(
  
  Covariate = rownames(fixed_effects),
  
  Posterior_Mean = fixed_effects$mean,
  
  SD = fixed_effects$sd,
  
  Lower_95_CrI = fixed_effects$`0.025quant`,
  
  Upper_95_CrI = fixed_effects$`0.975quant`
  
)

write.csv(
  fixed_effect_table,
  "outputs/tables/BYM2_Fixed_Effects.csv",
  row.names = FALSE
)


# ==============================================================================
# 24. BYM2 HYPERPARAMETERS
# ==============================================================================

hyperparameters <- model_bym2$summary.hyperpar

hyperparameter_table <- data.frame(
  
  Parameter = rownames(hyperparameters),
  
  Posterior_Mean = hyperparameters$mean,
  
  SD = hyperparameters$sd,
  
  Lower_95_CrI = hyperparameters$`0.025quant`,
  
  Upper_95_CrI = hyperparameters$`0.975quant`
  
)

write.csv(
  hyperparameter_table,
  "outputs/tables/BYM2_Hyperparameters.csv",
  row.names = FALSE
)


# ==============================================================================
# 25. POSTERIOR SPATIAL RANDOM EFFECT
# ==============================================================================

spatial_effect <- model_bym2$summary.random$district_index

districts$posterior_spatial_mean <-
  spatial_effect$mean[
    seq_len(nrow(districts))
  ]

districts$posterior_spatial_lower <-
  spatial_effect$`0.025quant`[
    seq_len(nrow(districts))
  ]

districts$posterior_spatial_upper <-
  spatial_effect$`0.975quant`[
    seq_len(nrow(districts))
  ]


# ==============================================================================
# 26. POSTERIOR PROBABILITY OF POSITIVE SPATIAL EFFECT
# ==============================================================================

# Calculate P(u_i > 0 | data) using the posterior mean and SD.
#
# This is preferable to simply counting positive posterior means.
# A district is considered strongly supported as positive when:
#
#   P(u_i > 0 | data) > 0.80
#
# This threshold is reported explicitly and should not be confused
# with simply having a positive posterior mean.

districts$posterior_probability_positive <-
  1 - pnorm(
    0,
    mean = spatial_effect$mean[
      seq_len(nrow(districts))
    ],
    sd = spatial_effect$sd[
      seq_len(nrow(districts))
    ]
  )

districts$posterior_probability_negative <-
  pnorm(
    0,
    mean = spatial_effect$mean[
      seq_len(nrow(districts))
    ],
    sd = spatial_effect$sd[
      seq_len(nrow(districts))
    ]
  )

districts$spatial_exceedance_class <- case_when(
  
  districts$posterior_probability_positive > 0.80 ~
    "Positive spatial effect: P(u > 0) > 0.80",
  
  districts$posterior_probability_negative > 0.80 ~
    "Negative spatial effect: P(u < 0) > 0.80",
  
  TRUE ~
    "Uncertain spatial effect"
  
)


# ==============================================================================
# 27. EXCEEDANCE SUMMARY
# ==============================================================================

exceedance_summary <- data.frame(
  
  Category = c(
    "Positive: P(u > 0) > 0.80",
    "Negative: P(u < 0) > 0.80",
    "Uncertain"
  ),
  
  Number_of_districts = c(
    
    sum(
      districts$posterior_probability_positive > 0.80,
      na.rm = TRUE
    ),
    
    sum(
      districts$posterior_probability_negative > 0.80,
      na.rm = TRUE
    ),
    
    sum(
      districts$posterior_probability_positive <= 0.80 &
        districts$posterior_probability_negative <= 0.80,
      na.rm = TRUE
    )
    
  )
)

write.csv(
  exceedance_summary,
  "outputs/tables/Posterior_Exceedance_Summary.csv",
  row.names = FALSE
)


# ==============================================================================
# 28. EXPORT DISTRICT-LEVEL SPATIAL RESULTS
# ==============================================================================

spatial_results <- districts %>%
  st_drop_geometry() %>%
  select(
    REGCODE,
    Anemic,
    net_use,
    posterior_spatial_mean,
    posterior_spatial_lower,
    posterior_spatial_upper,
    posterior_probability_positive,
    posterior_probability_negative,
    spatial_exceedance_class,
    LISA_cluster,
    Local_I,
    Z_score,
    P_value
  )

write.csv(
  spatial_results,
  "outputs/spatial/District_Spatial_Posterior_Results.csv",
  row.names = FALSE
)


# ==============================================================================
# 29. POSTERIOR SPATIAL EFFECT MAP
# ==============================================================================

tmap_mode("plot")

spatial_effect_map <-
  
  tm_shape(districts) +
  
  tm_polygons(
    col = "posterior_spatial_mean",
    palette = "-RdYlBu",
    style = "quantile",
    n = 5,
    border.col = "grey60",
    lwd = 0.1,
    title = "Posterior spatial effect"
  ) +
  
  tm_layout(
    frame = FALSE,
    legend.position = c(
      "right",
      "bottom"
    )
  )


tmap_save(
  spatial_effect_map,
  "outputs/figures/Posterior_Spatial_Effect_Map.png",
  width = 8.5,
  height = 9,
  units = "in",
  dpi = 300
)


# ==============================================================================
# 30. POSTERIOR EXCEEDANCE MAP
# ==============================================================================

exceedance_map <-
  
  tm_shape(districts) +
  
  tm_polygons(
    col = "spatial_exceedance_class",
    palette = "Set2",
    border.col = "grey60",
    lwd = 0.1,
    title = "Posterior spatial classification"
  ) +
  
  tm_layout(
    frame = FALSE,
    legend.position = c(
      "right",
      "bottom"
    )
  )


tmap_save(
  exceedance_map,
  "outputs/figures/Posterior_Exceedance_Map.png",
  width = 8.5,
  height = 9,
  units = "in",
  dpi = 300
)


# ==============================================================================
# 31. LISA MAP
# ==============================================================================

lisa_map <-
  
  tm_shape(districts) +
  
  tm_polygons(
    col = "LISA_cluster",
    palette = "Set3",
    border.col = "grey60",
    lwd = 0.1,
    title = "LISA cluster"
  ) +
  
  tm_layout(
    frame = FALSE,
    legend.position = c(
      "right",
      "bottom"
    )
  )


tmap_save(
  lisa_map,
  "outputs/figures/LISA_Clusters.png",
  width = 8.5,
  height = 9,
  units = "in",
  dpi = 300
)


# ==============================================================================
# 32. SAVE ANALYTICAL DATA WITH SPATIAL RESULTS
# ==============================================================================

st_write(
  districts,
  "outputs/spatial/Bayesian_Spatial_Anemia_Results.gpkg",
  delete_dsn = TRUE,
  quiet = TRUE
)


# ==============================================================================
# 33. MODEL OBJECTS
# ==============================================================================

saveRDS(
  model_baseline,
  "outputs/model_baseline.rds"
)

saveRDS(
  model_bym2,
  "outputs/model_bym2.rds"
)


# ==============================================================================
# 34. SESSION INFORMATION
# ==============================================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  "outputs/session_info.txt"
)


# ==============================================================================
# 35. ANALYSIS SUMMARY
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("BAYESIAN SPATIAL ANAEMIA ANALYSIS COMPLETED\n")
cat("============================================================\n")

cat(
  "Number of districts:",
  nrow(districts),
  "\n"
)

cat(
  "Moran's I:",
  round(
    moran_result$estimate[1],
    4
  ),
  "\n"
)

cat(
  "Moran's I p-value:",
  format.pval(
    moran_result$p.value,
    digits = 4
  ),
  "\n"
)

cat(
  "BYM2 DIC:",
  round(
    model_bym2$dic$dic,
    3
  ),
  "\n"
)

cat(
  "BYM2 WAIC:",
  round(
    model_bym2$waic$waic,
    3
  ),
  "\n"
)

cat(
  "Districts with P(u > 0) > 0.80:",
  sum(
    districts$posterior_probability_positive > 0.80,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Districts with P(u < 0) > 0.80:",
  sum(
    districts$posterior_probability_negative > 0.80,
    na.rm = TRUE
  ),
  "\n"
)

cat("============================================================\n")
cat("All outputs saved under the 'outputs/' directory.\n")
cat("============================================================\n")