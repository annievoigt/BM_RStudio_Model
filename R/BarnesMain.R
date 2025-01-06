# load libraries
library(ggplot2)
library(data.table)
library(brms)
library(bayesplot)

# load data files
all_files <- list.files("testing/distances_between_holes_poked_csv", pattern = ".csv")
test_me <- unlist(strsplit(all_files[1], split = "_"))[5]
grep(x = test_me, pattern = "R[[1-9]]*")
sub("*D[[1-9]]*", "", test_me)
substr(test_me, 1, 2)

unlist(strsplit(test_me, split = "[A-Z]"))

# Load Barnes Maze files
exp_data <- rbindlist(lapply(X = all_files, FUN = function(file_name) {
  tmp_dt <- fread(paste0("testing/distances_between_holes_poked_csv/", file_name))
  split_string <- unlist(strsplit(file_name, split = "_"))
  conditions_string <- split_string[6]
  condition_vec <- unlist(strsplit(conditions_string, split = "[A-Z]"))
  mouse_id <- gsub(x = tail(split_string, 1), pattern = ".csv", replacement = "")
  tmp_dt[,`:=`("MouseID"= mouse_id, "Round" = as.integer(condition_vec[2]), "Day" = as.integer(condition_vec[3]), "Trial" = as.integer(condition_vec[4])),]
}))
exp_data <- exp_data[,TrialUnique := paste0(Round, "_", Day, "_", Trial),]

# add target hole according to round
exp_data[Round==1,target_hole := 2,]
exp_data[Round==2,target_hole := 16,]

# sort key by hierarchy
setkeyv(x = exp_data, c("Round", "Day", "Trial","start_frame"))

exp_data[,`:=`(hole_number=hole_number+1,previous_hole=previous_hole+1),] # remove when data offset is corrected in files

# transform data to radians and x y coordinates on a unit circle
exp_data[,`:=`(target_rad = as.numeric(target_hole) * pi / 20 * 2 - pi/20 -pi),][,`:=`(x_target = sin(target_rad), y_target = -cos(target_rad)),]
exp_data[,`:=`(previous_hole_rad = as.numeric(previous_hole) * pi / 20 * 2 - pi/20 -pi - pi, hole_number_rad = as.numeric(hole_number) * pi / 20 * 2 - pi/20 -pi),]
exp_data[,`:=`(x1 = sin(hole_number_rad), y1 = -cos(hole_number_rad), x2 = sin(previous_hole_rad), y2 = -cos(previous_hole_rad)),]

# define starting point as 0,0 (center of unit circle)
exp_data[previous_hole==0,`:=`(x2=0,y2=0),]

# calculate distance to target hole and distance to target hole at time t-1
exp_data[,`:=`(distance=sqrt((x1-x2)^2+(y1-x2)^2), target_distance=sqrt((x2-x_target)^2+(y2-y_target)^2)),]
exp_data[,previous_target_distance:=shift(target_distance, fill = 1),by=TrialUnique]

# set target hole and day as factor (make it discrete) and add levels
exp_data[,`:=`(target_hole=factor(target_hole), Day=factor(Day, levels=c("1", "2", "3", "4", "5", "12"))),]

# Plot Barnes Maze data

# plot overview of Barnes Maze data as graphs
ggplot(data = exp_data, aes(x = x1, y = y1, color = start_frame, group = interaction(MouseID, TrialUnique))) +
  geom_point() +
  geom_point(data = exp_data, aes(x = x_target, y = y_target), color = "red", shape = 1, inherit.aes = F)+
  geom_path() +
  facet_wrap(paste("Round", Round) + Day ~ paste("Trial", Trial)) +
  coord_fixed() +
  theme_void()

# Plot density of Barnes Maze data
ggplot(data = exp_data, aes(x = x1, y = y1, group=Round)) +
  geom_density_2d_filled(contour_var = "ndensity")+ 
  facet_wrap(paste("Round", Round) + Day ~ paste("Trial", Trial)) +
  coord_fixed() +
  theme_void()

# y ~ 1 + X + (1|Day) example of formula for hierarchical model

# set up brms model formula
model_formula <- brmsformula(hole_number_rad ~ x2*y2 + s(target_distance) + s(previous_target_distance) + Trial  + (Trial|Day) + target_hole,
                             kappa ~ x2*y2 + s(target_distance) + s(previous_target_distance)+ (Trial|Day) + Trial + target_hole,
                             family = von_mises(), center = TRUE)

# get prior for model formula --> check that you set all required cooefficient priors
get_prior(model_formula, data = exp_data, family = von_mises(link_kappa = "log", link = ), center = TRUE)

# set priors
new_priors <- c(set_prior("normal(0, 1)", class = "b", coef = c("starget_distance_1", "sprevious_target_distance_1", "x2", "y2", "x2:y2","Trial", "target_hole16"), dpar = "kappa"),
                set_prior("normal(0, 1)", class = "b", coef = c("starget_distance_1", "sprevious_target_distance_1","x2", "y2", "x2:y2", "Trial", "target_hole16")), 
                set_prior("normal(10, 5)", class = "Intercept"),
                set_prior("cauchy(0, 0.1)", class = "Intercept", dpar = "kappa"))

# fit brms model with 4 chains, 2000 iterations, 1000 warmup, 4 cores and init = 0 (starting point of iteration) and adapt_delta = 0.95
fit_prior <- brm(formula = model_formula, data = exp_data, family = von_mises(),
                 prior = new_priors, chains = 4, iter = 2000, warmup = 1000,
                 cores = 4, init = "0", control = list(adapt_delta = 0.95),
                 sample_prior = "only")
summary(fit_prior)

# plot model summary
plot(fit_prior)

# get y_rep from brms fit
y_rep_prior <- posterior_predict(fit_prior, newdata = exp_data)

# get hole number by back transforming radians to hole number
y_rep_prior_hole <- round((y_rep_prior+pi + pi/20) / 2 / pi * 20)

# plot posterior predictive checks of hole number
pp_check(fit_prior, ndraws = 1, type = "dens_overlay")

fit <- brm(formula = model_formula, data = exp_data, family = von_mises(), 
           prior = new_priors, chains = 4, iter = 2000, warmup = 1000,
           cores = 4, init = "0", control = list(adapt_delta = 0.95))
summary(fit)

# plot model summary
plot(fit)

# get y_rep from brms fit
y_rep <- posterior_predict(fit, newdata = exp_data)

# get hole number by back transforming radians to hole number
y_rep_hole <- round((y_rep+pi + pi/20) / 2 / pi * 20)

# plot posterior predictive checks of hole number
bayesplot::ppc_dens_overlay_grouped(y = exp_data$hole_number, yrep = y_rep_hole[1:50,], group = exp_data$Round, bw = .1)+
  geom_vline(data = exp_data[,.(target_hole=as.integer(as.character(unique(target_hole)))),by=Round], aes(xintercept = target_hole, group = Round), colour="red")

# Calculate percentage of holes
round_i <- exp_data[Round==1,index,]
round_i2 <- exp_data[Round==2,index,]
sum(y_rep_hole[,round_i]==exp_data[Round==1,unique(target_hole),])/length(y_rep_hole[,round_i])
sum(y_rep_hole[,round_i]==exp_data[Round==2,unique(target_hole),])/length(y_rep_hole[,round_i])
sum(y_rep_hole[,round_i2]==exp_data[Round==2,unique(target_hole),])/length(y_rep_hole[,round_i2])
sum(y_rep_hole[,round_i2]==exp_data[Round==1,unique(target_hole),])/length(y_rep_hole[,round_i2])

# plot posterior graphs from model for specific days (1) and compare rounds
index_day <- exp_data[Trial==1&Day==1&Round==1,index,]
index_day2 <- exp_data[Trial==1&Day==1&Round==2,index,]
plot(c(0, sin(y_rep[1,index_day])), c(0,-cos(y_rep[1,index_day])), type = "l", col = scales::alpha("blue", 0.1), xlim = c(-1, 1), ylim = c(-1, 1))
for (i in 2:500) {
  lines(c(0,sin(y_rep[i,index_day])), c(0,-cos(y_rep[i,index_day])), col = scales::alpha("blue", 0.1))
}

for (i in 1:500) {
  lines(c(0,sin(y_rep[i,index_day2])), c(0,-cos(y_rep[i,index_day2])), col = scales::alpha("red", 0.1))
}

points(sin(exp_data[Trial==1&Day==1&Round==1,target_rad][1]), -cos(exp_data[Trial==1&Day==1&Round==1,target_rad][1]), col = "black")
points(sin(exp_data[Trial==1&Day==1&Round==2,target_rad][1]), -cos(exp_data[Trial==1&Day==1&Round==2,target_rad][1]), col = "green")

# plot posterior predictive checks of radians
bayesplot::ppc_intervals_grouped(y = exp_data$hole_number_rad, yrep = y_rep, group = exp_data$Round)

# plot kappa for different conditions
brms::conditional_effects(fit, dpar = "kappa", robust = T)

# print stan code of model (can be used to compile from scratch)
brms::make_stancode(model_formula, prior=new_priors, family=von_mises(), data=exp_data)