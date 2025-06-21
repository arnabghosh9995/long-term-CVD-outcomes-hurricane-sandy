
# Load data file
all_data_quarters_wshape2 <- readRDS("cvd_hurricane.rds")
regions_data <- all_data_quarters_wshape2 %>%
  filter(state_flood %in% c("CT", "NYC", "NJ"))
newdata<- regions_data %>% # Create the adjacency matrix
  filter(QUARTER == 1) %>%
  select_('indx_zip2', 'geometry', 'state_flood') %>%
  mutate(zip_code = as.factor(indx_zip2)) %>%
  arrange(zip_code)
alldata.adj <- poly2nb(newdata$geometry) 

nb2INLA("map.adj", alldata.adj)
g <- inla.read.graph(filename = "map.adj")


regions_data$idarea <- as.numeric(as.factor(regions_data$indx_zip2)) # Create an alternative variable for area and time
regions_data$idarea1 <- regions_data$idarea
regions_data$idtime <- 1 + regions_data$year - min(regions_data$year)

# Moran's I calculations
library(spdep)
library(dplyr)

newdata<- regions_data %>%
  filter(QUARTER == 1) %>%
  select_('indx_zip2', 'geometry', 'state_flood') %>%
  mutate(zip_code = as.factor(indx_zip2)) %>%
  arrange(zip_code)

nb <- poly2nb(newdata$geometry, queen = TRUE) # Islands found
no_neighbors <- which(card(nb) == 0) #Check for islands
newdata <- newdata[-no_neighbors, ] # Removed islands


nb <- poly2nb(newdata$geometry, row.names = newdata$indx_zip2, queen = TRUE) # Rerun with names of list added from ZIP codes
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

regions_data$indx_zip2 <- as.character(regions_data$indx_zip2)
lw_ids <- attr(lw, "region.id")

morans_results <- regions_data %>%
  group_by(QUARTER) %>%
  summarise(
    moran_test = list(
      tryCatch(
        {
          values_named <- setNames(CVD, indx_zip2)[lw_ids] # change the outcome
          moran.test(values_named, lw)
        },
        error = function(e) NA
      )
    ),
    .groups = "drop"
  ) %>%
  mutate(
    moran_I = sapply(moran_test, function(x) if (is.list(x)) x$estimate["Moran I statistic"] else NA),
    p_value = sapply(moran_test, function(x) if (is.list(x)) x$p.value else NA)
  ) %>%
  dplyr::select(QUARTER, moran_I, p_value)

print(morans_results)

# Matching with propensity scores
regions_data_imputed <- regions_data %>% # Imputing with median value based on distribution
  mutate(fraction_rent_occupied_unt = ifelse(is.na(fraction_rent_occupied_unt), median(fraction_rent_occupied_unt, na.rm = TRUE), fraction_rent_occupied_unt),
         median_hshld_incm_norm = ifelse(is.na(median_hshld_incm_norm), median(median_hshld_incm_norm, na.rm = TRUE), median_hshld_incm_norm),
         fraction_racial_minor = ifelse(is.na(fraction_racial_minor), median(fraction_racial_minor, na.rm = TRUE), fraction_racial_minor))


saveRDS(regions_data_imputed, file = "regions_data_imputed.rds") # Saved version for sensitivity of matching strategy

pre_regional_analysis <- regions_data_imputed %>%
  filter(hurricane == 0)
match_model <- matchit(stormsurgedummy ~ mean_age + fraction_age65over + fraction_rent_occupied_unt + Charlson_adj_avrg + ADI_NATRANK + fraction_female + fraction_overcrowd_housng + median_hshld_incm_norm + fraction_racial_minor + state_flood, data = pre_regional_analysis, method = "nearest",ratio = 3, distance = "logit", replace = TRUE, distance.options = list())

summary(match_model, un = FALSE) # Assess the average absolute within-pair difference between each covariate
plot(match_model, type = "jitter", interactive = FALSE) # Checking the matching using SMD
matched_df <- match.data(match_model) # Applying the modeling to the entire dataset # Extract model
matched_ids <- matched_df %>% pull(indx_zip2)
matched_dataset <- regions_data_imputed %>%
  filter(indx_zip2 %in% matched_ids)

saveRDS(matched_dataset, "matched_dataset")

# Checking the parallel trends assumption
pre_treatment_data <- subset(matched_dataset, hurricane == 0)  # assuming hurricane == 0 is pre-treatment
mod_pretrend <- CVD ~ as.factor(year)*stormsurgedummy + f(idarea, model = 'bym', graph = g) + mean_age + fraction_age65over + fraction_rent_occupied_unt + Charlson_adj_avrg + ADI_NATRANK + fraction_female + as.factor(year) + fraction_overcrowd_housng + median_hshld_incm_norm + fraction_racial_minor + state_flood + offset(log(n_bene_in_zip))
model_pretrend <- inla(mod_pretrend, 
                       data = pre_treatment_data, 
                       family = "nbinomial", 
                       E = n_bene_in_zip) 
summary(model_pretrend)

###############################
# Setting up output functions #
###############################
# Extraction for the event study
extract_exp_last_es_variables <- function(model) {
  # Extract the summary of fixed effects
  coefficients <- model$summary.fixed
  
  # Get the last 5 rows (last 5 variables) and exponentiate the mean and CI
  last_5_variables <- tail(coefficients, 31)
  exp_last_variables <- data.frame(
    Variable = rownames(last_5_variables),
    Mean = exp(last_5_variables$mean),
    `2.5%` = exp(last_5_variables$`0.025quant`),
    `97.5%` = exp(last_5_variables$`0.975quant`)
  )
  
  # Return the exponentiated results for the last 5 variables
  return(exp_last_variables)
}

# For extracting the DiD estimator
extract_exp_last_variable_did <- function(model) {
  # Extract the mean and confidence interval (0.025 and 0.975 quantiles) for fixed effects
  coefficients <- model$summary.fixed
  
  # Get the last row (last variable) and exponentiate the mean and CI
  last_variable <- tail(coefficients, 1)
  exp_last_variable <- data.frame(
    Variable = rownames(last_variable),
    Mean = exp(last_variable$mean),
    `2.5%` = exp(last_variable$`0.025quant`),
    `97.5%` = exp(last_variable$`0.975quant`)
  )
  
  # Return the exponentiated results for the last variable
  return(exp_last_variable)
}

# Extraction for the event study (all relevant coefficients)
extract_exp_last_es_variables <- function(model) {
  # Extract the summary of fixed effects
  coefficients <- model$summary.fixed
  
  # Get the last 31 rows and exponentiate the mean and CIs
  last_31_variables <- tail(coefficients, 31)
  exp_last_variables <- data.frame(
    Variable = rownames(last_31_variables),
    Mean = exp(last_31_variables$mean),
    Lower = exp(last_31_variables$`0.025quant`),
    Upper = exp(last_31_variables$`0.975quant`)
  )
  
  # Create a new column 'Time' by extracting numbers from 'Variable'
  exp_last_variables$Time <- gsub("[^0-9]", "", exp_last_variables$Variable)
  exp_last_variables$Time <- as.character(exp_last_variables$Time)
  
  # Add a reference row (neutral effect) after the 9th row
  library(tibble)
  exp_last_variables <- add_row(exp_last_variables,
                                Variable = "0",
                                Time = "10",
                                Mean = 1,
                                Lower = 1,
                                Upper = 1,
                                .after = 9)
  
  return(exp_last_variables)

  ######################
  # INLA outputs - DID #
  ######################
  # CVD
  modelformula.allregions.CVD <- CVD ~ stormsurgedummy*hurricane + f(idarea, model = 'bym', graph = g) + f(idtime, model = 'rw1') + mean_age + fraction_age65over + fraction_rent_occupied_unt + Charlson_adj_avrg + ADI_NATRANK + fraction_female + as.factor(year) + fraction_overcrowd_housng + median_hshld_incm_norm + fraction_racial_minor + state_flood + offset(log(n_bene_in_zip))
  
  modelINLA.allregions.CVD.matched <-  inla(modelformula.allregions.CVD, family = "nbinomial", 
                                            data = matched_dataset, 
                                            control.family = list(link = 'log'), 
                                            inla.mode = "experiment", 
                                            control.compute=list(openmp.strategy = "huge", dic=TRUE, cpo=TRUE), 
                                            control.predictor = list(link=1),
                                            control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                            verbose = TRUE,
                                            debug = TRUE, 
                                            safe = TRUE)
  
  modelINLA.allregions.CVD.matched.50 <-  inla(modelformula.allregions.CVD, family = "nbinomial", 
                                               data = matched_dataset %>% filter(n_bene_in_zip >= 50), 
                                               control.family = list(link = 'log'), 
                                               inla.mode = "experiment", 
                                               control.compute=list(openmp.strategy = "huge", dic=TRUE, cpo=TRUE), 
                                               control.predictor = list(link=1),
                                               control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                               verbose = FALSE,
                                               debug = TRUE, 
                                               safe = TRUE)
  
  modelINLA.allregions.CVD.matched.100 <-  inla(modelformula.allregions.CVD, family = "nbinomial", 
                                                data = matched_dataset %>% filter(n_bene_in_zip >= 100), 
                                                control.family = list(link = 'log'), 
                                                inla.mode = "experiment", 
                                                control.compute=list(openmp.strategy = "huge", dic=TRUE, cpo=TRUE), 
                                                control.predictor = list(link=1),
                                                control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                                verbose = FALSE,
                                                debug = TRUE, 
                                                safe = TRUE)
  
  extract_exp_last_variable_did(modelINLA.allregions.CVD.matched)
  extract_exp_last_variable_did(modelINLA.allregions.CVD.matched.50)
  extract_exp_last_variable_did(modelINLA.allregions.CVD.matched.100)
  
  # AMI
  modelformula.allregions.AMI <- AMI ~ stormsurgedummy*hurricane + f(idarea, model = 'bym', graph = g) + f(idtime, model = 'rw1') + mean_age + fraction_age65over + fraction_rent_occupied_unt + Charlson_adj_avrg + ADI_NATRANK + fraction_female + as.factor(year) + fraction_overcrowd_housng + median_hshld_incm_norm + fraction_racial_minor + state_flood + offset(log(n_bene_in_zip))
  
  modelINLA.allregions.ami.matched <-  inla(modelformula.allregions.AMI, family = "nbinomial", 
                                            data = matched_dataset, 
                                            control.family = list(link = 'log'), 
                                            inla.mode = "experiment", 
                                            control.compute=list(openmp.strategy = "huge", dic=TRUE, cpo=TRUE), 
                                            control.predictor = list(link=1),
                                            control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                            verbose = FALSE,
                                            debug = TRUE, 
                                            safe = TRUE)

  modelINLA.allregions.ami.matched.50 <-  inla(modelformula.allregions.AMI, family = "nbinomial", 
                                               data = matched_dataset %>% filter(n_bene_in_zip >= 50), 
                                               control.family = list(link = 'log'), 
                                               inla.mode = "experiment", 
                                               control.compute=list(openmp.strategy = "huge", dic=TRUE, cpo=TRUE), 
                                               control.predictor = list(link=1),
                                               control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                               verbose = FALSE,
                                               debug = TRUE, 
                                               safe = TRUE)
  
  modelINLA.allregions.ami.matched.100 <-  inla(modelformula.allregions.AMI, family = "nbinomial", 
                                                data = matched_dataset %>% filter(n_bene_in_zip >= 100), 
                                                control.family = list(link = 'log'), 
                                                inla.mode = "experiment", 
                                                control.compute=list(openmp.strategy = "huge", dic=TRUE, cpo=TRUE), 
                                                control.predictor = list(link=1),
                                                control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                                verbose = FALSE,
                                                debug = TRUE, 
                                                safe = TRUE)

  extract_exp_last_variable_did(modelINLA.allregions.ami.matched)
  extract_exp_last_variable_did(modelINLA.allregions.ami.matched.50)
  extract_exp_last_variable_did(modelINLA.allregions.ami.matched.100)
  
  # HF
  modelformula.allregions.chf <- CHF ~ stormsurgedummy*hurricane + f(idarea, model = 'bym', graph = g) + f(idtime, model = 'rw1') + mean_age + fraction_age65over + fraction_rent_occupied_unt + Charlson_adj_avrg + ADI_NATRANK + fraction_female + as.factor(year) + fraction_overcrowd_housng + median_hshld_incm_norm + fraction_racial_minor + state_flood + offset(log(n_bene_in_zip))
  
  modelINLA.allregions.chf.matched <-  inla(modelformula.allregions.chf, family = "nbinomial", 
                                            data = matched_dataset, 
                                            control.family = list(link = 'log'), 
                                            inla.mode = "experiment", 
                                            control.compute=list(openmp.strategy = "huge", dic=TRUE, cpo=TRUE), 
                                            control.predictor = list(link=1),
                                            control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                            verbose = FALSE,
                                            debug = TRUE, 
                                            safe = TRUE)
  
  modelINLA.allregions.chf.matched.50 <-  inla(modelformula.allregions.chf, family = "nbinomial", 
                                               data = matched_dataset %>% filter(n_bene_in_zip >= 50), 
                                               control.family = list(link = 'log'), 
                                               inla.mode = "experiment", 
                                               control.compute=list(openmp.strategy = "huge", dic=TRUE, cpo=TRUE), 
                                               control.predictor = list(link=1),
                                               control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                               verbose = FALSE,
                                               debug = TRUE, 
                                               safe = TRUE)
  
  modelINLA.allregions.chf.matched.100 <-  inla(modelformula.allregions.chf, family = "nbinomial", 
                                                data = matched_dataset %>% filter(n_bene_in_zip >= 100), 
                                                control.family = list(link = 'log'), 
                                                inla.mode = "experiment", 
                                                control.compute=list(openmp.strategy = "huge", dic=TRUE, cpo=TRUE), 
                                                control.predictor = list(link=1),
                                                control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                                verbose = FALSE,
                                                debug = TRUE, 
                                                safe = TRUE)
  
  extract_exp_last_variable_did(modelINLA.allregions.chf.matched)
  extract_exp_last_variable_did(modelINLA.allregions.chf.matched.50)
  extract_exp_last_variable_did(modelINLA.allregions.chf.matched.100)
  
  # Stroke
  modelformula.allregions.stroke <- stroke ~ stormsurgedummy*hurricane + f(idarea, model = 'bym', graph = g) + f(idtime, model = 'rw1') + mean_age + fraction_age65over + fraction_rent_occupied_unt + Charlson_adj_avrg + ADI_NATRANK + fraction_female + as.factor(year) + fraction_overcrowd_housng + median_hshld_incm_norm + fraction_racial_minor + state_flood + offset(log(n_bene_in_zip))
  
  modelINLA.allregions.stroke.matched <-  inla(modelformula.allregions.stroke, family = "nbinomial", 
                                               data = matched_dataset, 
                                               control.family = list(link = 'log'), 
                                               inla.mode = "experiment", 
                                               control.compute=list(openmp.strategy = "huge", dic=TRUE, cpo=TRUE), 
                                               control.predictor = list(link=1),
                                               control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                               verbose = FALSE,
                                               debug = TRUE, 
                                               safe = TRUE)

  modelINLA.allregions.stroke.matched.50 <-  inla(modelformula.allregions.stroke, family = "nbinomial", 
                                                  data = matched_dataset %>% filter(n_bene_in_zip >= 50), 
                                                  control.family = list(link = 'log'), 
                                                  inla.mode = "experiment", 
                                                  control.compute=list(openmp.strategy = "huge", dic=TRUE, cpo=TRUE), 
                                                  control.predictor = list(link=1),
                                                  control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                                  verbose = FALSE,
                                                  debug = TRUE, 
                                                  safe = TRUE)
  
  modelINLA.allregions.stroke.matched.100 <-  inla(modelformula.allregions.stroke, family = "nbinomial", 
                                                   data = matched_dataset %>% filter(n_bene_in_zip >= 100), 
                                                   control.family = list(link = 'log'), 
                                                   inla.mode = "experiment", 
                                                   control.compute=list(openmp.strategy = "huge", dic=TRUE, cpo=TRUE), 
                                                   control.predictor = list(link=1),
                                                   control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                                   verbose = FALSE,
                                                   debug = TRUE, 
                                                   safe = TRUE)
  
  
  extract_exp_last_variable_did(modelINLA.allregions.stroke.matched)
  extract_exp_last_variable_did(modelINLA.allregions.stroke.matched.50)
  extract_exp_last_variable_did(modelINLA.allregions.stroke.matched.100)
  
  ##############################
  # INLA outputs - Event Study #
  ##############################
  matched_dataset %>%
    group_by(QUARTER) %>%
    summarise(count = n())
  
  matched_dataset$QUARTER <- as.factor(matched_dataset$QUARTER) 
  matched_dataset$QUARTER <- factor(matched_dataset$QUARTER, levels = sort(as.numeric(levels(matched_dataset$QUARTER))))
  matched_dataset$QUARTER <- relevel(matched_dataset$QUARTER, ref = "10") # one quarter before the hurricane strikes
  
  # CVD
  modelformula.allregions.CVDya <- CVD ~ QUARTER*stormsurgedummy + f(idarea1, model = 'bym', graph = g) + f(idtime, model = 'rw1') + mean_age + fraction_age65over + fraction_rent_occupied_unt + Charlson_adj_avrg + ADI_NATRANK + fraction_female + fraction_overcrowd_housng +  median_hshld_incm_norm + fraction_racial_minor + state_flood + offset(log(n_bene_in_zip))
  
  es.INLA.allregions.CVD.match <-  inla(modelformula.allregions.CVDya, family = "nbinomial", 
                                        data = matched_dataset, 
                                        control.family = list(link = 'log'), 
                                        inla.mode = "experiment", 
                                        control.compute=list(openmp.strategy = "huge", waic = TRUE, dic=TRUE, cpo=TRUE), 
                                        control.predictor = list(link=1),
                                        control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                        verbose = FALSE,
                                        debug = TRUE, 
                                        safe = TRUE)
  
  es_cvd <- extract_exp_last_es_variables(es.INLA.allregions.CVD.match)
  
  create_coef_plot_with_band(es_cvd, sheet = 1, x_var = 'Time', coef_var = 'Mean', lower_var = 'Lower', upper_var = 'Upper')
  
  modelformula.allregions.AMIya <- AMI ~ QUARTER*stormsurgedummy + f(idarea1, model = 'bym', graph = g) + f(idtime, model = 'rw1') + mean_age + fraction_age65over + fraction_rent_occupied_unt + Charlson_adj_avrg + ADI_NATRANK + fraction_female + fraction_overcrowd_housng +  median_hshld_incm_norm + fraction_racial_minor + state_flood + offset(log(n_bene_in_zip))
  
  # AMI
  es.INLA.allregions.AMI.match <-  inla(modelformula.allregions.AMIya, family = "nbinomial", 
                                        data = matched_dataset, 
                                        control.family = list(link = 'log'), 
                                        inla.mode = "experiment", 
                                        control.compute=list(openmp.strategy = "huge", waic = TRUE, dic=TRUE, cpo=TRUE), 
                                        control.predictor = list(link=1),
                                        control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                        verbose = FALSE,
                                        debug = TRUE, 
                                        safe = TRUE)
  
  es_ami <- extract_exp_last_es_variables(es.INLA.allregions.AMI.match)
  
  create_coef_plot_with_band(es_ami, sheet = 1, x_var = 'Time', coef_var = 'Mean', lower_var = 'Lower', upper_var = 'Upper')
  
  # HF
  modelformula.allregions.CHFya <- CHF ~ QUARTER*stormsurgedummy + f(idarea1, model = 'bym', graph = g) + f(idtime, model = 'rw1') + mean_age + fraction_age65over + fraction_rent_occupied_unt + Charlson_adj_avrg + ADI_NATRANK + fraction_female + fraction_overcrowd_housng +  median_hshld_incm_norm + fraction_racial_minor + state_flood + offset(log(n_bene_in_zip))

  es.INLA.allregions.CHF.match <-  inla(modelformula.allregions.CHFya, family = "nbinomial", 
                                        data = matched_dataset, 
                                        control.family = list(link = 'log'), 
                                        inla.mode = "experiment", 
                                        control.compute=list(openmp.strategy = "huge", waic = TRUE, dic=TRUE, cpo=TRUE), 
                                        control.predictor = list(link=1),
                                        control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                        verbose = FALSE,
                                        debug = TRUE, 
                                        safe = TRUE)
  
  es_chf <- extract_exp_last_es_variables(es.INLA.allregions.CHF.match)
  
  create_coef_plot_with_band(es_chf, sheet = 1, x_var = 'Time', coef_var = 'Mean', lower_var = 'Lower', upper_var = 'Upper')
  
# Stroke
  modelformula.allregions.strokeya <- stroke ~ QUARTER*stormsurgedummy + f(idarea1, model = 'bym', graph = g) + f(idtime, model = 'rw1') + mean_age + fraction_age65over + fraction_rent_occupied_unt + Charlson_adj_avrg + ADI_NATRANK + fraction_female + fraction_overcrowd_housng +  median_hshld_incm_norm + fraction_racial_minor + state_flood + offset(log(n_bene_in_zip))

  es.INLA.allregions.stroke.match <-  inla(modelformula.allregions.strokeya, family = "nbinomial", 
                                           data = matched_dataset, 
                                           control.family = list(link = 'log'), 
                                           inla.mode = "experiment", 
                                           control.compute=list(openmp.strategy = "huge", waic = TRUE, dic=TRUE, cpo=TRUE), 
                                           control.predictor = list(link=1),
                                           control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                           verbose = FALSE,
                                           debug = TRUE, 
                                           safe = TRUE)
  
  es_stroke <- extract_exp_last_es_variables(es.INLA.allregions.stroke.match)

  create_coef_plot_with_band(es_stroke, sheet = 1, x_var = 'Time', coef_var = 'Mean', lower_var = 'Lower', upper_var = 'Upper')
  
  
  #####################################
  ## Sensitivity analysis - Matching ##
  #####################################
  
  # Nearest neighbor 2:1 with replacement
  match_model_2to1 <- matchit(stormsurgedummy ~ mean_age + fraction_age65over + fraction_rent_occupied_unt + Charlson_adj_avrg + ADI_NATRANK + fraction_female + fraction_overcrowd_housng + median_hshld_incm_norm + fraction_racial_minor + state_flood, data = pre_regional_analysis, method = "nearest",ratio = 1, distance = "logit", replace = TRUE, distance.options = list())
  
  
  summary(match_model_2to1, un = FALSE) # Assess the average absolute within-pair difference between each covariate
  plot(match_model_2to1, type = "jitter", interactive = FALSE) # Checking the matching using SMD
  
  matched_df <- match.data(match_model_2to1) # Extract model
  matched_ids <- matched_df %>% pull(indx_zip2)
  matched_dataset2to1 <- regions_data_imputed %>%
    filter(indx_zip2 %in% matched_ids)
  
  modelINLA.allregions.CVD.matched2to1 <-  inla(modelformula.allregions.CVD, family = "nbinomial", 
                                            data = matched_dataset2to1, 
                                            control.family = list(link = 'log'), 
                                            inla.mode = "experiment", 
                                            control.compute=list(openmp.strategy = "huge", dic=TRUE, cpo=TRUE), 
                                            control.predictor = list(link=1),
                                            control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                            verbose = TRUE,
                                            debug = TRUE, 
                                            safe = TRUE)
  
  modelINLA.allregions.CVD.matched2to1.50 <-  inla(modelformula.allregions.CVD, family = "nbinomial", 
                                               data = matched_dataset2to1 %>% filter(n_bene_in_zip >= 50), 
                                               control.family = list(link = 'log'), 
                                               inla.mode = "experiment", 
                                               control.compute=list(openmp.strategy = "huge", dic=TRUE, cpo=TRUE), 
                                               control.predictor = list(link=1),
                                               control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                               verbose = FALSE,
                                               debug = TRUE, 
                                               safe = TRUE)
  
  modelINLA.allregions.CVD.matched2to1.100 <-  inla(modelformula.allregions.CVD, family = "nbinomial", 
                                                data = matched_dataset2to1 %>% filter(n_bene_in_zip >= 100), 
                                                control.family = list(link = 'log'), 
                                                inla.mode = "experiment", 
                                                control.compute=list(openmp.strategy = "huge", dic=TRUE, cpo=TRUE), 
                                                control.predictor = list(link=1),
                                                control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                                verbose = FALSE,
                                                debug = TRUE, 
                                                safe = TRUE)
  
  extract_exp_last_variable_did(modelINLA.allregions.CVD.matched2to1)
  extract_exp_last_variable_did(modelINLA.allregions.CVD.matched2to1.50)
  extract_exp_last_variable_did(modelINLA.allregions.CVD.matched2to1.100)
  
  # Nearest neighbor 1:1 with replacement
  match_model_1to1 <- matchit(stormsurgedummy ~ mean_age + fraction_age65over + fraction_rent_occupied_unt + Charlson_adj_avrg + ADI_NATRANK + fraction_female + fraction_overcrowd_housng + median_hshld_incm_norm + fraction_racial_minor + state_flood, data = pre_regional_analysis, method = "nearest", distance = "logit", replace = TRUE, ratio = 3)

  summary(match_model_1to1, un = FALSE) # Assess the average absolute within-pair difference between each covariate
  plot(match_model_1to1, type = "jitter", interactive = FALSE)
  
  matched_df <- match.data(match_model_1to1) # Extract model
  matched_ids <- matched_df %>% pull(indx_zip2)
  matched_dataset3to1 <- regions_data_imputed %>%
    filter(indx_zip2 %in% matched_ids)
  
  modelformula.allregions.CVD <- CVD ~ stormsurgedummy*hurricane + f(idarea, model = 'bym', graph = g) + f(idtime, model = 'rw1') + mean_age + fraction_age65over + fraction_rent_occupied_unt + Charlson_adj_avrg + ADI_NATRANK + fraction_female + as.factor(year) + fraction_overcrowd_housng + median_hshld_incm_norm + fraction_racial_minor + state_flood + offset(log(n_bene_in_zip))
  
  modelINLA.allregions.CVD.matched1to1 <-  inla(modelformula.allregions.CVD, family = "nbinomial", 
                                            data = match_model_1to1, 
                                            control.family = list(link = 'log'), 
                                            inla.mode = "experiment", 
                                            control.compute=list(openmp.strategy = "huge", dic=TRUE, cpo=TRUE), 
                                            control.predictor = list(link=1),
                                            control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                            verbose = TRUE,
                                            debug = TRUE, 
                                            safe = TRUE)
  

  modelINLA.allregions.CVD.matched1to1.50 <-  inla(modelformula.allregions.CVD, family = "nbinomial", 
                                               data = match_model_1to1 %>% filter(n_bene_in_zip >= 50), 
                                               control.family = list(link = 'log'), 
                                               inla.mode = "experiment", 
                                               control.compute=list(openmp.strategy = "huge", dic=TRUE, cpo=TRUE), 
                                               control.predictor = list(link=1),
                                               control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                               verbose = FALSE,
                                               debug = TRUE, 
                                               safe = TRUE)
  
  modelINLA.allregions.CVD.matched1to1.100 <-  inla(modelformula.allregions.CVD, family = "nbinomial", 
                                                data = match_model_1to1 %>% filter(n_bene_in_zip >= 100), 
                                                control.family = list(link = 'log'), 
                                                inla.mode = "experiment", 
                                                control.compute=list(openmp.strategy = "huge", dic=TRUE, cpo=TRUE), 
                                                control.predictor = list(link=1),
                                                control.inla(strategy = "adaptive", int.strategy = "ccd"),
                                                verbose = FALSE,
                                                debug = TRUE, 
                                                safe = TRUE)
  
  extract_exp_last_variable_did(modelINLA.allregions.CVD.matched1to1)
  extract_exp_last_variable_did(modelINLA.allregions.CVD.matched1to1.50)
  extract_exp_last_variable_did(modelINLA.allregions.CVD.matched1to1.100)