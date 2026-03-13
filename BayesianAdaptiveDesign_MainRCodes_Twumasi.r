options(repr.plot.res=300,repr.plot.weight=8,repr.plot.height=8)#setting plot size

# setting working directory
setwd("/Users/clementtwumasi/Desktop/ICR_CTSU_T_BayesianDesign_task")


# Loading some relevant R codes
library(ggplot2)
library(tidyverse)

# Part 1: Phase II Bayesian Adaptive Design – Simulation Task

###############################################################
# Bayesian Adaptive Phase II Trial Simulation
# Endpoint: Objective Response Rate (ORR)
# Model: Beta-Binomial
###############################################################

# set.seed(123)

###############################################################
# Trial design parameters
###############################################################

Nmax <- 30            # maximum sample size
Ninterim <- 12        # interim analysis point

prior_a <- 1          # Beta(a=1,b=1) prior parameter or assuming prior ~ Uniform(0,1)
prior_b <- 1

futility_prob <- 0.10 # stop if P(ORR>0.3) < 0.10
success_prob <- 0.80  # declare success if P(ORR>0.3) > 0.80

target_ORR <- 0.30

nsim <- 10000          # number of simulated trials


###############################################################
# Function to simulate ONE trial
###############################################################

simulate_trial <- function(true_ORR){

 # Assuming objective responses follow a binomial distribution
  responses <- rbinom(n=Nmax,size=1,prob=true_ORR)

  r_interim <- sum(responses[1:Ninterim])

  # posterior after interim
  post_a <- prior_a + r_interim
  post_b <- prior_b + Ninterim - r_interim
  # Posterior probability
  prob_promising <- 1 - pbeta(q=target_ORR,shape1=post_a,shape2=post_b)

  # Futility stopping
  if(prob_promising < futility_prob){

    return(list(
      stop_early=TRUE,
      success=FALSE,
      sample_size=Ninterim
    ))
  }

  # Continue to final analysis
  r_final <- sum(responses)

  post_a_final <- prior_a + r_final
  post_b_final <- prior_b + Nmax - r_final
  # Posterior probability
  prob_promising_final <- 1 - pbeta(q=target_ORR,shape1=post_a_final,shape2=post_b_final)

  success <- prob_promising_final > success_prob

  return(list(
    stop_early=FALSE,
    success=success,
    sample_size=Nmax
  ))
}


###############################################################
# Function to run simulations
###############################################################

run_sim <- function(true_ORR){

  results <- replicate(nsim,simulate_trial(true_ORR),simplify=FALSE)

  # proportion of times prob_promising < futility_prob
  stop_rate <- mean(sapply(results,function(x) x$stop_early))

 # proportion of times prob_promising_final > success_prob
  success_rate <- mean(sapply(results,function(x) x$success))

                              
  mean_sample <- mean(sapply(results,function(x) x$sample_size))

  return(data.frame(
    True_ORR=true_ORR,
    Early_Stop=stop_rate,
    Success=success_rate,
    Avg_Sample_Size=mean_sample
  ))
}


###############################################################
# Run simulation scenarios
###############################################################
# set.seed(1)
scenarios <- seq(0.1, 0.70, by=0.05)

sim_results <- map_df(scenarios,run_sim)

print(sim_results)

names(sim_results)

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

###############################################################
# Prepare data for Panel A (Success + Early Stop)
###############################################################

plot_data1 <- sim_results %>%
  select(True_ORR, Success, Early_Stop) %>%
  pivot_longer(cols = c(Success, Early_Stop),
               names_to = "Metric",
               values_to = "Probability")

plot_data1$Metric <- factor(plot_data1$Metric,
                            levels = c("Success", "Early_Stop"),
                            labels = c("Declare Success",
                                       "Early Stop (Futility)"))

###############################################################
# Panel A: Decision Performance
###############################################################

p1 <- ggplot(plot_data1,
             aes(x = True_ORR, y = Probability, color = Metric)) +

  geom_line(size = 1.3) +
  geom_point(size = 2.5) +

  scale_color_manual(values = c(
    "Declare Success" = "#1f78b4",
    "Early Stop (Futility)" = "#e31a1c"
  )) +

  labs(
    x = "",
    y = "Probability",
    color = "Decision Metric:",
    title="A)"
      # , title = "Operating Characteristics: Decision Performance"
  ) +

  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
     panel.border = element_rect(colour = "black", fill = NA)
      )+  geom_vline(xintercept = .3, linetype="dashed")+
  annotate("text", x = .22, y = .8, label = "Target ORR=0.30")


###############################################################
# Panel B: Expected Sample Size
###############################################################

p2 <- ggplot(sim_results,
             aes(x = True_ORR, y = Avg_Sample_Size)) +

  geom_line(size = 1.3, color = "#33a02c") +
  geom_point(size = 2.5, color = "#33a02c") +

  labs(
    x = "True Objective Response Rate (ORR)",
    y = "Expected Sample Size",
      title="B)"
      # , title = "Operating Characteristics: Trial Efficiency"
  ) +

  theme_minimal(base_size = 12)  +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
     panel.border = element_rect(colour = "black", fill = NA)
      )+  geom_vline(xintercept = .3, linetype="dashed")+
  annotate("text", x = .38, y = 25, label = "Target ORR=0.30")

###############################################################
# Combine Panels
###############################################################

final_OC_plot <- p1 / p2

final_OC_plot

###############################################################
# Predictive Probability of Success
###############################################################

predictive_success <- function(r, n, Nmax = 30,
                               target = 0.30,
                               success_prob = 0.80,
                               a0 = 1, b0 = 1,
                               nsim = 10000){

  # posterior parameters
  post_a <- a0 + r
  post_b <- b0 + n - r

  successes <- 0

  for(i in 1:nsim){

    # draw true ORR from posterior
    theta <- rbeta(1, post_a, post_b)

    # simulate future responses
    future <- rbinom(1, Nmax - n, theta)

    r_final <- r + future

    # posterior at final analysis
    post_a_final <- a0 + r_final
    post_b_final <- b0 + Nmax - r_final

    # posterior probability of success
    PP_final <- 1 - pbeta(target, post_a_final, post_b_final)

    if(PP_final > success_prob){
      successes <- successes + 1
    }
  }

  PPoS <- successes/nsim

  return(PPoS)
}

# Example
predictive_success(r = 3, n = 12)

# Compute Posterior Probability for All Possible Outcomes

###############################################################
# Decision boundary calculations
###############################################################

# Trial parameters
n_interim <- 12
target <- 0.30
futility_threshold <- 0.10

prior_a <- 1
prior_b <- 1

# possible response counts
responses <- 0:n_interim

results <- data.frame()

for(r in responses){

  post_a <- prior_a + r
  post_b <- prior_b + n_interim - r

  PP <- 1 - pbeta(target, post_a, post_b)

  decision <- ifelse(PP < futility_threshold,
                     "Stop for Futility",
                     "Continue")

  results <- rbind(results,
                   data.frame(
                     Responses = r,
                     PosteriorProb = PP,
                     Decision = decision
                   ))
}

results

# library(ggplot2)

ggplot(results,
       aes(x = Responses,
           y = PosteriorProb,
           colour = Decision)) +

  geom_point(size = 4) +

  geom_line(aes(group = 1), linewidth = 1) +

  geom_hline(yintercept = futility_threshold,
             linetype = "dashed",
             colour = "black") +

  scale_colour_manual(values = c(
    "Stop for Futility" = "red",
    "Continue" = "darkgreen"
  )) +

  labs(
    x = "Number of Responses at Interim (n = 12)",
    y = "Posterior Probability ORR > 30%"
      # ,title = "Bayesian Interim Decision Boundary"
      ,colour = "Decision:"
  ) +

  theme_minimal(base_size = 13)+
theme(
  legend.position = "top",
  axis.line = element_line(),
  axis.ticks.x = element_line(),   # restore x ticks
  axis.text.x = element_text(size = 8),  # ensure labels show
  panel.border = element_rect(colour = "black", fill = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

# R Code to Generate the Operating Characteristic (OC) Surface
# NB: The OC surface shows how trial decisions behave under different true treatment effects.

###############################################################
# Operating Characteristic Surface
###############################################################

true_ORR <- seq(0.05,0.70,by=0.05)

n_interim <- 12
Nmax <- 30

grid_results <- data.frame()

for(theta in true_ORR){

  for(r in 0:n_interim){

    post_a <- 1 + r
    post_b <- 1 + n_interim - r

    nsim <- 3000
    success <- 0

    for(i in 1:nsim){

      theta_draw <- rbeta(1,post_a,post_b)

      future <- rbinom(1,Nmax-n_interim,theta_draw)

      r_final <- r + future

      post_a_final <- 1 + r_final
      post_b_final <- 1 + Nmax - r_final

      PP_final <- 1 - pbeta(0.30,post_a_final,post_b_final)

      if(PP_final > 0.80){
        success <- success + 1
      }

    }

    PPoS <- success/nsim

    grid_results <- rbind(grid_results,
                          data.frame(
                            True_ORR=theta,
                            Interim_Response=r,
                            PPoS=PPoS
                          ))
  }
}

# library(ggplot2)

ggplot(grid_results,
       aes(x=True_ORR,
           y=Interim_Response,
           fill=PPoS))+

geom_tile()+

scale_fill_gradient(
  low="white",
  high="darkblue"
)+

labs(
x="True Response Rate",
y="Responses at Interim (n=12)",
fill="PPoS:"
# ,title="Operating Characteristic Surface for Bayesian Adaptive Trial"
)+

theme_minimal(base_size=13)+
theme(
  # legend.position = "top",
  axis.line = element_line(),
  axis.ticks.x = element_line(),   # restore x ticks
  axis.text.x = element_text(size = 8),  # ensure labels show
  panel.border = element_rect(colour = "black", fill = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)


#Part 2: Bayesian Adaptive Phase II Design – Interim Analysis Reporting

##Main automated summary stats function to estimate summary statistics of the data (irrespective of variable type)
#For both categorical and numeric data adaptively
setwd("/Users/clementtwumasi/Desktop/ICR_CTSU_T_BayesianDesign_task")
source("Novel_summary_computation_script.R")

# importing the empirical dataset 

data_obs<- read.csv("ICR-CTSU Senior Statistician - Bayesian Adaptive Design Task - Simulation and Reporting Exercise - Data.csv")

# view the observed empirical data
head(data_obs, n=10)

# Number of patients (sample size)

cat("The number of patients in this cohort is:", "", length(unique(data_obs$Patient_ID)), "\n")

# Check for missing values in each variable
missing_value_count<- function(x) sum(is.na(x))

paste("Number of missingness out of N=",dim(data_obs)[1],"","participants:")

# Print the missing counts
print(apply(data_obs,2,missing_value_count))

# visualising missingness
naniar::gg_miss_var(data_obs[,-1])

naniar::mcar_test(data_obs)

# setting some character variables as factor variables (to be able to summarise the data using my novel summary function

data_obs$Visit<- as.factor(data_obs$Visit)
data_obs$RECIST_Status<- as.factor(data_obs$RECIST_Status)
data_obs$Comments<- as.factor(data_obs$Comments)

head(data_obs)

#View categorical variables
cat("Categorical variables of the data:","\n")
Categorical_covariates<- Categorical_variables(data=data_obs)
print(Categorical_covariates[-c(1,3)])  #ignoring patient ID, Dates

# Summarising the qualitative variables of the data
#Extract all qualitative data
# AT BASELINE
Data_categorical<- subset(data_obs, Visit=="Baseline")[ ,Categorical_covariates]

categorical_summary<- list()

for(index in 1:dim(Data_categorical)[2]){
  categorical_summary[[index]]<- Summary_stats(data=Data_categorical, variable_index=index,plot.graph=F)
}

Categorical_summaries<- do.call("rbind",categorical_summary)
subset(Categorical_summaries, 	n_freq!=0)

# Summarising the qualitative variables of the data
#Extract all qualitative data
# AT 3-Month
Data_categorical<- subset(data_obs, Visit=="3-Month")[ ,Categorical_covariates]

categorical_summary<- list()

for(index in 1:dim(Data_categorical)[2]){
  categorical_summary[[index]]<- Summary_stats(data=Data_categorical, variable_index=index,plot.graph=F)
}

Categorical_summaries<- do.call("rbind",categorical_summary)
subset(Categorical_summaries, 	n_freq!=0)

# Summarising the qualitative variables of the data
#Extract all qualitative data
# AT 3-Month
Data_categorical<- subset(data_obs, Visit=="6-Month")[ ,Categorical_covariates]

categorical_summary<- list()

for(index in 1:dim(Data_categorical)[2]){
  categorical_summary[[index]]<- Summary_stats(data=Data_categorical, variable_index=index,plot.graph=F)
}

Categorical_summaries<- do.call("rbind",categorical_summary)
subset(Categorical_summaries, 	n_freq!=0)

subset(Categorical_summaries, Variable_name=="RECIST_Status")[, 1:2]

#View quantitative variables
cat("Quantitative/numerical variables of the data:","\n")
Quantitative_covariates<- Quantitative_variables(data=data_obs)

print(Quantitative_covariates)

#Extract all quantitative data (AT 3-MONTH)
Data_quantitative<- data_obs[ ,Quantitative_covariates]
#Data_quantitative

numeric_summary<- list()

numeric_summary[[index]]<-  Summary_stats(data=Data_quantitative, 
                                    variable_index=index,plot.graph=FALSE) 



numeric_summaries<- do.call("rbind",numeric_summary)
numeric_summaries

#Extract all quantitative data (AT 3-MONTH)
Data_quantitative<- subset(data_obs, Visit=="3-Month")[ ,Quantitative_covariates]
#Data_quantitative

numeric_summary<- list()

numeric_summary[[index]]<-  Summary_stats(data=Data_quantitative, 
                                    variable_index=index,plot.graph=FALSE) 



numeric_summaries<- do.call("rbind",numeric_summary)
numeric_summaries

#Extract all quantitative data (AT 6-MONTH)
Data_quantitative<- subset(data_obs, Visit=="6-Month")[ ,Quantitative_covariates]
#Data_quantitative

numeric_summary<- list()

numeric_summary[[index]]<-  Summary_stats(data=Data_quantitative, 
                                    variable_index=index,plot.graph=FALSE) 



numeric_summaries<- do.call("rbind",numeric_summary)
numeric_summaries

subset(data_obs, Patient_ID=="P-01")

subset(data_obs, Patient_ID=="P-05")

subset(data_obs, Patient_ID=="P-15")

subset(data_obs, Patient_ID=="P-16")

levels(data_obs$RECIST_Status)


table(data_obs$RECIST_Status)


# reordering the factor levels of RECIST Status
data_obs$RECIST_Status<- factor(data_obs$RECIST_Status, levels=c('CR', 'PR', 'SD', 'PD', 'NE', 'Pending', 'Target Lesions Identified'))

levels(data_obs$RECIST_Status)

# Selecting categories of interest
data_obs2 <- data_obs %>%
filter(RECIST_Status %in% c("CR","PR","SD","PD"))

head(data_obs2)

unique(data_obs2$Visit)


data_obs2$Visit<- as.character(data_obs2$Visit)

visit_counts2 <- data_obs %>%
  group_by(Visit) %>%
  summarise(n = n())

visit_counts2

visit_counts <- data_obs2 %>%
  group_by(Visit) %>%
  summarise(n = n())

visit_counts 

visit_labels <- setNames(
  paste0(visit_counts$Visit, " (n = ", visit_counts$n, ")"),
  visit_counts$Visit
)


visit_labels

###############################################################
# Waterfall plot
###############################################################
library(ggplot2)
library(RColorBrewer)
library(tidytext)   # for reorder_within()

color_selected <- brewer.pal(7, "Set1")

recist_palette <- c(
  "CR" = color_selected[1],
  "PR" = color_selected[2],
  "SD" = "chartreuse",
  "PD" = color_selected[4]
)

ggplot(data_obs2,
       aes(x = reorder_within(gsub("P-", "", Patient_ID), Percent_Change, Visit),
           y = Percent_Change,
           fill = RECIST_Status)) +

  geom_bar(stat = "identity") +

  # Percentage labels (slanted)
  geom_text(
    aes(label = paste0(round(Percent_Change,1), "%"),
        y = ifelse(Percent_Change > 0,
                   Percent_Change + 3,
                   Percent_Change - 3)),
    size = 2.1,
    angle = 45,        # slanted labels
    hjust = .5
  ) +

  scale_fill_manual(values = recist_palette) +

  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ Visit, labeller = labeller(Visit = visit_labels), scales = "free_x")+
  # facet_wrap(~ Visit, scales = "free_x") +

  scale_x_reordered() +

  labs(
    x = "Patient ID",
    y = "Change in tumour size from Baseline (%)",
    fill = "RECIST Status:"
  ) +

  theme_minimal(base_size = 12) +
scale_x_reordered() +
theme(
  legend.position = "top",
  axis.line = element_line(),
  axis.ticks.x = element_line(),   # restore x ticks
  axis.text.x = element_text(size = 8),  # ensure labels show
  panel.border = element_rect(colour = "black", fill = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

###############################################################
# Keep only RECIST visits relevant for response
###############################################################

data_clean <- data_obs %>%
filter(RECIST_Status %in% c("CR","PR","SD","PD"))

###############################################################
# Convert to binary response
###############################################################

data_clean <- data_clean %>%
mutate(
Response = ifelse(RECIST_Status %in% c("CR","PR"),1,0)
)

dim(data_clean)
dim(data_obs)

data_clean

###############################################################
# Posterior summary function (Posterior Probability and 95% Credible Interval)
###############################################################

posterior_summary <- function(responses){

n <- length(responses)
r <- sum(responses)

post_a <- prior_a + r
post_b <- prior_b + n - r

median_est <- qbeta(0.5,post_a,post_b)
lower <- qbeta(0.025,post_a,post_b)
upper <- qbeta(0.975,post_a,post_b)

# i.e., prob_above_target=Prob_ORR_above_30
prob_above_target <- 1 - pbeta(q=target_ORR,shape1=post_a,shape2=post_b)

data.frame(
Responses=r,
Sample=n,
Median_ORR=median_est,
CI_Lower=lower,
CI_Upper=upper,
Prob_ORR_above_30=prob_above_target
)

}

# Interim Posterior Calculation

###############################################################
# Compute interim ORR
###############################################################

responses <- data_clean$Response

posterior_summary(responses)

###############################################################
# Interim decision rule
###############################################################

posterior <- posterior_summary(responses)

if(posterior$Prob_ORR_above_30 < 0.10){

decision <- "Stop for futility"

}else{

decision <- "Continue recruitment"

}

print(decision)

set.seed(123)
theta <- seq(0,1,length=1000)

prior <- dbeta(theta,prior_a ,prior_b)

post_a <- prior_a + sum(responses)
post_b <- prior_b + length(responses) - sum(responses)

posterior <- dbeta(theta,post_a,post_b)

# density_posterior<- density(posterior, from=0, to=5, n=256)$y
# density_prior<- density(prior, from=0, to=1, n=256)$y


plot(theta,prior,type="l",lwd=2,ylim=c(0,max(prior,posterior)),
     ylab="Density",xlab="True Response Rate (θ)")

lines(theta,posterior,col="red",lwd=3)

legend("topright",
       legend=c("Prior","Posterior"),
       col=c("black","red"),
       lwd=2)


# Additional Posterior Predictive Check 

theta_post <- rbeta(10000,prior_a+sum(responses),prior_b+length(responses)-sum(responses))


predicted_responses <- rbinom(10000,30,theta_post)

hist(predicted_responses, col="blue", main="Posterior distribution of responses", xlab="Predicted responses")












