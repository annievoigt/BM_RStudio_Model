run_model <- function(data, output_path) {
  # set up brms model formula
  model_formula <- brmsformula(hole_number_rad ~ x2*y2 + s(target_distance) + target_distance:treatment + treatment + s(previous_target_distance) + target_hole + Trial + (Trial|Day) + (x2*y2+target_distance:treatment|MouseID),
                               kappa ~ x2*y2 + s(target_distance) + target_distance:treatment + treatment + s(previous_target_distance)+ target_hole + Trial + running_trial + (Trial|Day) + (running_trial+x2*y2+target_distance:treatment|MouseID),
                               family = von_mises(), center = TRUE)
  
  # get prior for model formula --> check that you set all required cooefficient priors
  #get_prior(model_formula, data = exp_data, family = von_mises(link_kappa = "log", link = ), center = TRUE)
  
  # set priors
  new_priors <- c(set_prior("normal(0, 1)", class = "b", coef = c(    "running_trial", "sprevious_target_distance_1", "starget_distance_1", "target_hole12", "target_hole16", "target_hole8", "treatmentCNO",
                                                                      "treatmentCNO:target_distance", "treatmentM:target_distance", "treatmentSaline", "treatmentSaline:target_distance",
                                                                      "Trial", "x2", "x2:y2","y2"), dpar = "kappa"),
                  set_prior("normal(0, 1)", class = "b", coef = c("sprevious_target_distance_1", "starget_distance_1", "target_hole12", "target_hole16", "target_hole8", "treatmentCNO",
                                                                  "treatmentCNO:target_distance", "treatmentM:target_distance", "treatmentSaline", "treatmentSaline:target_distance",
                                                                  "Trial", "x2", "x2:y2", "y2")), 
                  set_prior("normal(10, 5)", class = "Intercept"),
                  set_prior("cauchy(0, 0.1)", class = "Intercept", dpar = "kappa"))
  
  fit <- brm(formula = model_formula, data = exp_data, family = von_mises(), 
             prior = new_priors, chains = 4, iter = 20, warmup = 1,
             cores = 4, init = "0", control = list(adapt_delta = 0.95), save_all_pars = T, file = "R/results/model_fit.rds")
  summary(fit)

  # plot model summary
  plot(fit)
}
