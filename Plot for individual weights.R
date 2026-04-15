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


# evidence synthesis ------------------------------------------------------

study_names <- c('P1', 'P3', 'P4', 'P5', 'P6', 'P7', 'P9', 'P10', 'P11')


GORICA <- list(
  output_P1$result[,4],
  output_P3$result[,4],
  output_P4$result[,4],
  output_P5$result[,4],
  output_P6$result[,4],
  output_P7$result[,4],
  output_P9$result[,4],
  output_P10$result[,4],
  output_P11$result[,4]
  )

evSyn_Added <- evSyn(GORICA, hypo_names = output_P1$result[,1],
                     study_names = study_names,
                     order_studies = c("ascending"))
evSyn_Added

plot(evSyn_Added)


#evSyn_Added_desc <- evSyn(GORICA, hypo_names = output_P1$result[,1],
#                     study_names = study_names,
#                     order_studies = c("descending"))
#evSyn_Added_desc
#plot(evSyn_Added_desc)


###


LL <- list(
  output_P1$result[,2],
  output_P3$result[,2],
  output_P4$result[,2],
  output_P5$result[,2],
  output_P6$result[,2],
  output_P7$result[,2],
  output_P9$result[,2],
  output_P10$result[,2],
  output_P11$result[,2]
)

PT <- list(
  output_P1$result[,3],
  output_P3$result[,3],
  output_P4$result[,3],
  output_P5$result[,3],
  output_P6$result[,3],
  output_P7$result[,3],
  output_P9$result[,3],
  output_P10$result[,3],
  output_P11$result[,3]
)


## Added-evidence approach, but with LL and PT
#evSyn_Added <- evSyn(object = LL, PT = PT, hypo_names = output_P1$result[,1],
#                     study_names = study_names,
#                     order_studies = c("ascending"))
#evSyn_Added
#plot(evSyn_Added)


# Equal-evidence approach
evSyn_Equal <- evSyn(object = LL, PT = PT, type = "equal", 
                     hypo_names = output_P1$result[,1],
                     study_names = study_names,
                     order_studies = c("ascending"))
evSyn_Equal

plot(evSyn_Equal)





################################################################################


# Exploratory, post-hoc analysis for future research

GORICA_post <- list(
  output_P1_post$result[,4],
  output_P3_post$result[,4],
  output_P4_post$result[,4],
  output_P5_post$result[,4],
  output_P6_post$result[,4],
  output_P7_post$result[,4],
  output_P9_post$result[,4],
  output_P10_post$result[,4],
  output_P11_post$result[,4]
)

evSyn_Added_post <- evSyn(GORICA_post, hypo_names = output_P1_post$result[,1],
                          study_names = study_names,
                          order_studies = c("ascending"))
evSyn_Added_post

plot(evSyn_Added_post)


###

LL_post <- list(
  output_P1_post$result[,2],
  output_P3_post$result[,2],
  output_P4_post$result[,2],
  output_P5_post$result[,2],
  output_P6_post$result[,2],
  output_P7_post$result[,2],
  output_P9_post$result[,2],
  output_P10_post$result[,2],
  output_P11_post$result[,2]
)

PT_post <- list(
  output_P1_post$result[,3],
  output_P3_post$result[,3],
  output_P4_post$result[,3],
  output_P5_post$result[,3],
  output_P6_post$result[,3],
  output_P7_post$result[,3],
  output_P9_post$result[,3],
  output_P10_post$result[,3],
  output_P11_post$result[,3]
)


## Added-evidence approach, but with LL and PT
#evSyn_Added_post <- evSyn(object = LL_post, PT = PT_post, 
#                          hypo_names = output_P1_post$result[,1],
#                          study_names = study_names,
#                          order_studies = c("ascending"))
#evSyn_Added_post
#plot(evSyn_Added_post)



# Equal-evidence approach
evSyn_Equal_post <- evSyn(object = LL_post, PT = PT_post, 
                          hypo_names = output_P1_post$result[,1],
                          study_names = study_names,
                          order_studies = c("ascending"))
evSyn_Equal_post

plot(evSyn_Equal_post)


