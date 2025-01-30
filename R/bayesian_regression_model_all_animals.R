#!/usr/bin/env Rscript

# Check if running from Snakemake or command line
if (exists("snakemake")) {
  # Running from Snakemake
  #output_plot_path <- snakemake@output[["plot"]]
  csv_files <- snakemake@input[["csv_files"]]
  model_fit_path <- snakemake@output[["model_fit"]]
  trial_info_path <- snakemake@config[["trial_info_xlsx"]]
} else {
  # Running from command line or RStudio
  #output_plot_path <- "R/results/results.pdf"
  csv_files <- list.files("R/all_csv/", pattern = "*.csv", full.names = TRUE)
  model_fit_path <- "R/results/model_fit.rds"
  trial_info_path <- "R/trial_info.xlsx"
}

# load libraries
library(ggplot2)
library(data.table)
library(brms)
library(bayesplot)
library(readxl)  # for reading Excel files

# Read trial info from Excel file
trial_info <- data.table(read_excel(trial_info_path))[!is.na(`Animal No.`),]
setkey(trial_info, Filename)

# Load Barnes Maze files
exp_data <- rbindlist(lapply(X = csv_files, FUN = function(file_name) {
  tmp_dt <- fread(file_name)
  # Get the basename of the file for lookup in trial_info
  base_filename <- gsub(pattern = ".csv", x = basename(file_name), replacement = "")
  tmp_dt[, Filename := base_filename]
  # Get trial info for this file
  return(tmp_dt)
}))
exp_data <- merge.data.table(exp_data, trial_info, by = "Filename")
setnames(exp_data, old = c("Animal No.", "Correct hole"), new = c("MouseID", "target_hole"))
exp_data <- exp_data[, TrialUnique := paste0(MouseID,"_",Round, "_", Day, "_", Trial)]

# add target hole according to round
exp_data[Round==1,target_hole := 2,]
exp_data[Round==2,target_hole := 16,]
exp_data[Round==3,target_hole := 12,]
exp_data[Round==4,target_hole := 8,]

# sort key by hierarchy
setkeyv(x = exp_data, c("Round", "Day", "Trial", "start_frame"))

# Transform data to radians and calculate x, y coordinates on a unit circle
exp_data[, target_rad := as.numeric(target_hole) * (2 * pi / 20) - (pi / 20) - pi,][
  ,`:=`(  x_target = sin(target_rad),
          y_target = -cos(target_rad))]

exp_data[, `:=`(
  previous_hole_rad = as.numeric(previous_hole) * (2 * pi / 20) - (pi / 20) - pi - pi,
  hole_number_rad = as.numeric(hole_number) * (2 * pi / 20) - (pi / 20) - pi
)]

exp_data[, `:=`(
  x1 = sin(hole_number_rad),
  y1 = -cos(hole_number_rad),
  x2 = sin(previous_hole_rad),
  y2 = -cos(previous_hole_rad)
)]

# define starting point as 0,0 (center of unit circle)
exp_data[previous_hole==0,`:=`(x2=0,y2=0),]

# calculate distance to target hole and distance to target hole at time t-1
exp_data[,`:=`(distance=sqrt((x1-x2)^2+(y1-y2)^2), target_distance=sqrt((x2-x_target)^2+(y2-y_target)^2)),]
exp_data[,previous_target_distance:=shift(target_distance, fill = 1),by=TrialUnique]

# set target hole and day as factor (make it discrete) and add levels
exp_data[,`:=`(target_hole=factor(target_hole), Day=factor(Day, levels=c("1", "2", "3", "4", "5", "12"))),]
exp_data[,index := .I]

exp_data[,poke_count:=1:.N, by = .(MouseID, Round, Day, Trial)]
# get the running trial number for each animal
exp_data[,running_trial:=0,][poke_count==1,running_trial:=1, by = .(MouseID, Round, Day)][,running_trial:=cumsum(running_trial),by=MouseID]




# Save all plots to the specified output path
#pdf(output_plot_path)

# Plot Barnes Maze data

# plot overview of Barnes Maze data as graphs

unique_animal_rounds <- exp_data[,.N,by=.(MouseID, Round)]

# plotting animals
# lapply(seq_along(unique_animal_rounds$N), FUN = function(index) {
#   ggplot(data = exp_data[MouseID==unique_animal_rounds$MouseID[index]&Round==unique_animal_rounds$Round[index]], aes(x = x1, y = y1, color = start_frame, group = interaction(MouseID, TrialUnique))) +
#     ggtitle(paste(unique_animal_rounds$MouseID[index], unique_animal_rounds$Round[index])) +
#     geom_point() +
#     geom_point(data = exp_data[MouseID==unique_animal_rounds$MouseID[index]&Round==unique_animal_rounds$Round[index]], aes(x = x_target, y = y_target, group=Trial), color = "red", shape = 1, inherit.aes = F)+
#     geom_path() +
#     facet_wrap(paste("Day", Day) ~ paste("Trial", Trial), ncol = 4, nrow = 6) +
#     guides(colour = guide_colorbar(title = "Frame")) +
#     coord_fixed() +
#     theme_void()
#   
#   ggplot(data = exp_data, aes(x = x1, y = y1, group=Round)) +
#     geom_density_2d_filled(contour_var = "ndensity")+ 
#     ggtitle(paste(unique_animal_rounds$MouseID[index], unique_animal_rounds$Round[index])) +
#     geom_point(data = exp_data[MouseID==unique_animal_rounds$MouseID[index]&Round==unique_animal_rounds$Round[index]], aes(x = x_target, y = y_target, group=Trial), color = "red", shape = 1, inherit.aes = F)+
#     facet_wrap(paste("Day", Day) ~ paste("Trial", Trial), ncol = 4, nrow = 6) +
#     coord_fixed() +
#     theme_void()
# })

# y ~ 1 + X + (1|Day) example of formula for hierarchical model

setnames(exp_data, old = c("CNO or Saline"), new = c("treatment"))

model_formula <- brmsformula(hole_number_rad ~ x2*y2 + s(target_distance) + target_distance:treatment + s(previous_target_distance) + target_hole*treatment + Trial + log(num_frames)*treatment + (Trial|Day) + (treatment+x2*y2+target_distance:treatment|MouseID),
                             kappa ~ x2*y2 + s(target_distance) + target_distance:treatment + s(previous_target_distance)+ target_hole*treatment + Trial + running_trial + log(num_frames)*treatment + (Trial|Day) + (treatment+running_trial+x2*y2+target_distance:treatment|MouseID),
                             family = von_mises(), center = TRUE)

# get prior for model formula --> check that you set all required cooefficient priors
get_prior(model_formula, data = exp_data, family = von_mises(link_kappa = "log", link = ), center = TRUE)

# set priors
new_priors <- c(set_prior("normal(0, 1)", class = "b", coef = c("lognum_frames", "running_trial", "sprevious_target_distance_1", "starget_distance_1", "target_hole12", "target_hole12:treatmentCNO",
                                                                "target_hole12:treatmentSaline", "target_hole16", "target_hole16:treatmentCNO", "target_hole16:treatmentSaline", "target_hole8",
                                                                "target_hole8:treatmentCNO", "target_hole8:treatmentSaline", "treatmentCNO", "treatmentCNO:lognum_frames", "treatmentCNO:target_distance",
                                                                "treatmentM:target_distance", "treatmentSaline", "treatmentSaline:lognum_frames", "treatmentSaline:target_distance", "Trial", "x2", "x2:y2", "y2"), dpar = "kappa"),
                set_prior("normal(0, 1)", class = "b", coef = coef_vector <- c("lognum_frames", "sprevious_target_distance_1", "starget_distance_1", "target_hole12", "target_hole12:treatmentCNO",
                                                                               "target_hole12:treatmentSaline", "target_hole16", "target_hole16:treatmentCNO", "target_hole16:treatmentSaline",
                                                                               "target_hole8", "target_hole8:treatmentCNO", "target_hole8:treatmentSaline", "treatmentCNO", "treatmentCNO:lognum_frames",
                                                                               "treatmentCNO:target_distance", "treatmentM:target_distance", "treatmentSaline", "treatmentSaline:lognum_frames",
                                                                               "treatmentSaline:target_distance", "Trial", "x2", "x2:y2", "y2")), 
                set_prior("normal(10, 5)", class = "Intercept"),
                set_prior("cauchy(0, 0.1)", class = "Intercept", dpar = "kappa"))

available_threads <- snakemake@threads
chains <- 4
threads <- floor(available_threads/chains)
fit <- brm(formula = model_formula, data = exp_data, family = von_mises(), 
           prior = new_priors, chains = chains, iter = 2000, warmup = 1000, threads = threading(threads),
           cores = 4, init = "0", control = list(adapt_delta = 0.9), save_all_pars = T, file = model_fit_path)
summary(fit)

# plot model summary
plot(fit)