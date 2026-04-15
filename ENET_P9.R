
# Load libraries. 
# These contain functions such as 'goric' that will be used in this R code. 
# Each time you re-open this R file you have to execute this step.
#
# Install the packages (once) 
# If you want to use restriktor from github:
#if (!require("devtools")) install.packages("devtools")
#library(devtools) 
#install_github("LeonardV/restriktor")
#library(restriktor) # for goric function
# If from CRAN:
if (!require("restriktor")) install.packages("restriktor")
#
# Load 
library(restriktor) # for goric function

################################################################################


# Notably, it is only possible to load the data if you are using the correct 
# working directory (with both your R script and data file). 
# The command `getwd()` shows you your current working directory.
#
# You can change the working directory to the one you prefer using the function 
# `setwd()` by specifying the correct location between parentheses. 
# Alternatively, in Rstudio, you can use the "Session" tab (on top) or you can 
# use the "Files"-pane (on top of probably the right lower box of Rstudio, this 
# pane is located next to the panes for "Plots", "Packages", "Help", and "Viewer").


### Analysis E-NET per person ### 

## Data preparation ##

# Read/Load the Data # 
# If you open the data file in a text editor, 
# you can see that the variable labels have been inserted in the first line of the file, which is called a header. 
# Therefore, you have to specify 'header = TRUE' when loading the data:
data_ENET_P9 <- read.table("data/ENET09.dat", header=TRUE)
data_ENET_P9

# Missings #
# Declare -999 as missing
data_ENET_P9$event[data_ENET_P9$event == -999] <- NA 
data_ENET_P9

# Make the variable 'event' a factor #
#
# Since a .txt file was loaded, R does not know the measurement levels of the 
# variables and assumes all of them to be continuous, meaning that there
# measurement level is interval or ratio. 
# Hence, especially when there are more than two groups, 
# one has to tell R that the variable `event` is a factor by using the 
# `factor()` function on the `event` variable (i.e., a grouping / categorical / 
# nominal variable):
#
data_ENET_P9$event <- factor(data_ENET_P9$event) 
data_ENET_P9
# this command tells R that it is a factor 
# and not a continuous variable like PCL5 and CLL

# Make the time variables #
#
# Center time such that at 1/3th the value is 0, 
# then the maximum/minimum will take place there:
percentile <- 1/3
data_ENET_P9$Time <- data_ENET_P9$Time - quantile(data_ENET_P9$Time, percentile) 
# If you want the variable time to have a 0 value, then you could use:
#data_ENET_P9$Time <- data_ENET_P9$Time - data_ENET_P9$Time[round(length(data_ENET_P9$Time)*percentile)]
#
# Make the quadratic term of time:
data_ENET_P9$sqTime <- data_ENET_P9$Time^2 # quadratic term


# Fit the (unconstrained) MANOVA model #
lm_fit_P9 <- lm(cbind(PCL5, CLL) ~ Time + sqTime + event, data = data_ENET_P9)
#
# You should not look at the results before specifying your hypotheses.
# If you, at some point, do want to inspect the results, then use:
#lm_fit_P9
#summary(lm_fit_P9)


# Extract parameter estimates of interest and their covariance matrix #
#
# The goric function can take some fit objects and it can use the unconstrained
# parameter estimates of interest and their covariance matrix.
# We will use the latter.
#
# parameter estimates with self-assigned labels:
est <- c(coef(lm_fit_P9))
names(est) <- c(paste0("beta_PCL5_", 0:3), paste0("beta_CLL_", 0:3)) 
# Specify restrictions in hypotheses using those names/labels.
# 0 is then intercept, 
# 1 linear relationship with Time, 
# 2 linear relationship with sqTime, thus the quadratic relationship with Time, 
# 3 the relationship with event.
#
# the covariance matrix of the parameter estimates from above
VCOV <- vcov(lm_fit_P9)
#
# Again, you should not look at the results before specifying your hypotheses.


## Hypotheses specification ##
#
# These are known before inspecting the results.
# We have to use the labels/names of the estimates.
# Restrictions in a hypothesis can be separated by a ',', ';' or '&',
# here, we will use the ',' to separate restrictions for one outcome
# and '&' to separate the restrictions for the two outcomes.
#
# No trend
H_NoTrend <- 'beta_PCL5_1 = 0, beta_PCL5_2 = 0 & beta_CLL_1 = 0, beta_CLL_2 = 0'
# beta_PCL5_1 = 0 & beta_CLL_1 = 0 -> no linear effect of Time
# beta_PCL5_2 = 0 & beta_CLL_2 = 0 -> no quadratic effect of Time
#
# For PCL5 linear downward trend and for CLL a linear upward trend
H_Linear <- 'beta_PCL5_1 < 0, beta_PCL5_2 = 0 & beta_CLL_1 > 0; beta_CLL_2 = 0'
# beta_PCL5_1 < 0 & beta_CLL_1 > 0 -> negative/positive linear effect of Time
# beta_PCL5_2 = 0 & beta_CLL_2 = 0 -> no quadratic effect of Time
#
# For PCL5 a hyperbolic trend and for CLL a linear upward trend
H_Hyperbolic <- 'beta_PCL5_2 < 0  & beta_CLL_1 > 0; beta_CLL_2 = 0'
# beta_PCL5_2 < 0 (and beta_PCL5_1 free) -> 
#                       negative quadratic (hyperbolic) trend for PCL5: 
#                       first an increase in symptoms and afterwards a decrease
# beta_CLL_2 = 0 -> no quadratic effect of Time
# beta_CLL_1 > 0 -> linear effect of Time and a positive one
#
# Unconstrained as failsafe, this will be included by default when there is
# more than one informative hypothesis (like here).


## GORICA evaluation ##
#
# Calculate GORIC values and weights
#
#In the calculation of the GORIC, an iterative process is needed to calculate 
#the penalty / complexity part. Therefore, one needs to set a seed value:
#1. Then, you will obtain the same penalty value every time you run this code.
#2. Then, you can change the seed value to check the sensitivity of the penalty 
#   value.
#   If it is sensitive, then increase number of iterations used in calculation 
#   of the penalty.
#
set.seed(123) # Set seed value
output_P9 <- goric(est, VCOV = VCOV, 
                   hypotheses = list(H_NoTrend = H_NoTrend, 
                                     H_Linear = H_Linear, 
                                     H_Hyperbolic = H_Hyperbolic))
output_P9
#summary(output_P9)
#output_P9$ratio.gw # matrix with ratios of GORICA weights


################################################################################


# Extra

# 'Null' populations
est_no <- est
est_no[c(2,3,6,7)] <- 0
est_lin <- est
est_lin[c(3,7)] <- 0
my_pop_est <- rbind("Null" = rep(0, length(est)), 
                    "NoTrend" = est_no, 
                    "Linear" = est_lin, 
                    "Observed" = est)
#
# Benchmarks
set.seed(123) # Set seed value
iter <- 100 # to decrease computation time
benchmarks_P9 <- benchmark(output_P9, pop_est = my_pop_est, iter = iter)
benchmarks_P9
#
plot(benchmarks_P9, log_scale = T, x_lim = c(0.01,13))$plots[1]
plot(benchmarks_P9, log_scale = T, x_lim = c(0.01,13))$plots[2]
#
# Our findings are likely under all the populations specified,
# thus including the null and no-trend population.
