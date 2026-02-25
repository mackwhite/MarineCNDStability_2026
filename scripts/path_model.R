librarian::shelf(tidyverse, vegan, readxl, dplyr, splitstackshape, codyn, lavaan,
                 MuMIn, corrplot, performance, ggeffects, ggpubr, parameters, ggstats,
                 brms, mixedup, rstatix, sf, ggspatial, semPlot, tidySEM)

vars_to_plot <- dat_ready[, c("s_rich_mean", "s_div_mean", "t_rich_mean", "t_div_mean", 
                              "spp_turnover", "troph_turnover", "spp_synchrony", "comm_n_stability")]
m <- cor(vars_to_plot, use = "complete.obs")
corrplot(m, 
         method = "circle",      
         type = "upper",         
         addCoef.col = "black",  
         tl.col = "black",       
         diag = FALSE)           

### step one
path_model <- '
  # 1. Stability Direct Effects
  comm_n_stability ~ s_rich_mean + s_div_mean + t_rich_mean + t_div_mean + 
                     troph_turnover + spp_turnover + spp_synchrony

  # 2. Mediators
  spp_synchrony  ~ s_rich_mean + s_div_mean + t_rich_mean + t_div_mean
  spp_turnover   ~ s_rich_mean + s_div_mean + t_rich_mean + t_div_mean
  troph_turnover ~ s_rich_mean + s_div_mean + t_rich_mean + t_div_mean
  
  # 3. Predictor Covariances (FREE estimation - no multipliers)
  # We just use ~~ to say "these are related"
  s_rich_mean ~~ s_div_mean
  t_rich_mean ~~ t_div_mean
  s_rich_mean ~~ t_rich_mean
  s_div_mean  ~~ t_div_mean
'

fit <- sem(path_model, data = dat_ready)
summary(
      fit,
      standardized = TRUE,
      fit.measures = TRUE,
      rsquare = TRUE
)

### step two
path_model <- '
  # 1. Stability Direct Effects
  comm_n_stability ~ s_rich_mean + troph_turnover + spp_turnover + spp_synchrony

  # 2. Mediators
  spp_synchrony  ~ s_rich_mean + s_div_mean + t_rich_mean + t_div_mean
  spp_turnover   ~ s_rich_mean + s_div_mean + t_rich_mean + t_div_mean
  troph_turnover ~ s_rich_mean + s_div_mean + t_rich_mean + t_div_mean
  
  # 3. Predictor Covariances (FREE estimation - no multipliers)
  # We just use ~~ to say "these are related"
  s_rich_mean ~~ s_div_mean
  t_rich_mean ~~ t_div_mean
  s_rich_mean ~~ t_rich_mean
  s_div_mean  ~~ t_div_mean
'

fit <- sem(path_model, data = dat_ready)
summary(
      fit,
      standardized = TRUE,
      fit.measures = TRUE,
      rsquare = TRUE
)

### step three
path_model <- '
  # 1. Stability Direct Effects
  comm_n_stability ~ s_rich_mean + troph_turnover + spp_turnover + spp_synchrony

  # 2. Mediators
  spp_synchrony  ~ s_rich_mean + s_div_mean + t_div_mean
  spp_turnover   ~ t_rich_mean
  troph_turnover ~ s_rich_mean + s_div_mean + t_div_mean
  
  # 3. Predictor Covariances (FREE estimation - no multipliers)
  # We just use ~~ to say "these are related"
  s_rich_mean ~~ s_div_mean
  t_rich_mean ~~ t_div_mean
  s_rich_mean ~~ t_rich_mean
  s_div_mean  ~~ t_div_mean
'

fit <- sem(path_model, data = dat_ready)
summary(
      fit,
      standardized = TRUE,
      fit.measures = TRUE,
      rsquare = TRUE
)

### step four
path_model <- '
  # 1. Stability Direct Effects
  comm_n_stability ~ s_rich_mean + spp_synchrony + spp_turnover + troph_turnover

  # 2. Mediators
  spp_synchrony  ~ s_div_mean 
  spp_turnover   ~ t_rich_mean
  troph_turnover ~ t_div_mean
  
  # 3. Predictor Covariances (FREE estimation - no multipliers)
  # We just use ~~ to say "these are related"
  s_rich_mean ~~ s_div_mean
  t_rich_mean ~~ t_div_mean
  s_rich_mean ~~ t_rich_mean
  s_div_mean  ~~ t_div_mean
'

fit <- sem(path_model, data = dat_ready)
summary(
      fit,
      standardized = TRUE,
      fit.measures = TRUE,
      rsquare = TRUE
)
modificationIndices(fit, sort = TRUE, maximum.number = 10)

path_model_v3 <- '
  # 1. Regressions (with labels for mediation)
  comm_n_stability ~ cp*s_rich_mean + b1*spp_synchrony + b2*spp_turnover
  
  spp_synchrony  ~ a1*s_rich_mean
  spp_turnover   ~ a2*t_rich_mean + d1*troph_turnover
  troph_turnover ~ a3*t_rich_mean

  # 2. Covariances
  s_rich_mean ~~ t_rich_mean

  # 3. Indirect Effects (Mediation)
  # Effect of richness via synchrony
  ind_rich_sync := a1 * b1
  
  # Effect of trophic richness via turnover
  ind_troph_turn := a3 * d1 * b2
  
  # Total Effect of s_rich_mean
  total_rich_effect := cp + (a1 * b1)
'

fit3 <- sem(path_model_v3, data = dat_ready)
summary(fit3, standardized = TRUE, fit.measures = TRUE)
modificationIndices(fit3, sort = TRUE, maximum.number = 5)

path_model_final <- '
   # 1. Regressions
   comm_n_stability ~ cp*s_rich_mean + b1*spp_synchrony + b2*spp_turnover
   
   spp_synchrony  ~ a1*s_rich_mean
   
   # ADDED: s_rich_mean now also predicts turnover (from MI 40)
   spp_turnover   ~ a2*t_rich_mean + d1*troph_turnover + a4*s_rich_mean
   
   troph_turnover ~ a3*t_rich_mean

   # 2. Covariances
   s_rich_mean ~~ t_rich_mean

   # 3. Indirect Effects
   ind_rich_sync := a1 * b1
   ind_rich_turn := a4 * b2        # New indirect path
   ind_troph_turn := a3 * d1 * b2
   
   total_rich_effect := cp + (a1 * b1) + (a4 * b2)
'
fit_final <- sem(path_model_final, data = dat_ready)
summary(fit_final, standardized = TRUE, fit.measures = TRUE)
modificationIndices(fit_final, sort = TRUE, maximum.number = 3)

cols <- c("royalblue", "firebrick")
semPaths(fit_final, 
         whatLabels = "std", 
         layout = "tree", 
         edge.label.cex = 1.2, 
         edge.color = cols,      # Applies colors based on sign
         posCol = "royalblue",   # Positive paths
         negCol = "firebrick",   # Negative paths
         fade = FALSE, 
         residuals = FALSE,
         sizeMan = 12, 
         edge.width = 1.5,
         style = "lisrel")
saveRDS(fit_final, "output/fit_final_path.rds")

path_model_ultimate <- '
   # 1. Regressions
   comm_n_stability ~ cp*s_rich_mean + b1*spp_synchrony + b2*spp_turnover
   
   spp_synchrony  ~ a1*s_rich_mean
   
   spp_turnover   ~ a2*t_rich_mean + d1*troph_turnover + a4*s_rich_mean
   
   # ADDED: s_rich_mean now also predicts trophic turnover (from MI 45)
   troph_turnover ~ a3*t_rich_mean + a5*s_rich_mean

   # 2. Covariances
   s_rich_mean ~~ t_rich_mean

   # 3. Indirect Effects (Updated to include the new chain)
   ind_rich_sync := a1 * b1
   ind_rich_turn := a4 * b2
   ind_troph_turn := a3 * d1 * b2
   
   # New indirect path: s_rich -> troph_turn -> spp_turn -> stability
   ind_rich_troph_chain := a5 * d1 * b2 
   
   total_rich_effect := cp + (a1 * b1) + (a4 * b2) + (a5 * d1 * b2)
'

# Use the Robust Estimator (MLR) to give the Chi-Square its best chance
fit_ultimate <- sem(path_model_ultimate, data = dat_ready, estimator = "MLR")
summary(fit_ultimate, standardized = TRUE, fit.measures = TRUE)

# Customizing node labels for clarity
node_labels <- c("S. Richness", "T. Richness", "Synchrony", 
                 "T. Turnover", "S. Turnover", "Stability")

# Define colors: Blue for positive, Red for negative
# Path coefficients: cp(+), a1(-), a2(-), d1(+), a4(+), a3(-), a5(+), b1(-), b2(-)
semPaths(fit_ultimate, 
         whatLabels = "std", 
         layout = "tree",       # Organized hierarchy from top to bottom
         rotation = 2,          # Rotates to flow left-to-right
         edge.label.cex = 1.1, 
         sizeMan = 12, 
         sizeLat = 10,
         color = "lightgrey", 
         edge.color = "black", 
         edge.width = 1.5,
         posCol = "blue",       # Custom positive color
         negCol = "red",        # Custom negative color
         fade = FALSE, 
         residuals = FALSE, 
         intercepts = FALSE)
inspect(fit_ultimate, "rsquare")
