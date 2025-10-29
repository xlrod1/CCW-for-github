library(dplyr)
library(purrr)
library(data.table)
library(cmprsk)
library(tidyr)
library(ggplot2)
library(survival)
library(foreach)
library(doParallel)
library(progressr)
library(slider)
library(zoo)


#num_cores <- parallel::detectCores() - 5


source("boot_SW_CIF_TV.R")
source("f_compute_survival_for_clone_SW_TV.R")
df.only.base<- read.csv("data_base_line_new.csv") 
measur_df<-fread("data_opening.csv") 

#select only the women who has time-varying confounders
df.only.base<-df.only.base %>% right_join(.,measur_df %>% select(id) %>% distinct()) 
 
#Complete cse analysis: select only the women who has complete information
list.of.confounder<-age ~nulliparity+preg.week+smoking+GBS
vars <- all.vars(list.of.confounder)
df.only.base<-df.only.base %>%
  # select only the vars in the formula, then drop any rows with NA
  filter(
    if_all(all_of(vars), ~ !is.na(.))
  )


#--------------------------------------------------------------------------------------------------------------
#Step 1:cloning+censoring. You can specify any treatment windows that you like.
#--------------------------------------------------------------------------------------------------------------
brks <- c(0,8,16,24,32,Inf)
brks.e <- c(0,8,16,24,32)

make_clone <- function(data, LB, UB = Inf) {
  # build the clone name, e.g. "window_0_12" or "window_12_Inf"
  clone_name <- paste0("window_", LB, "_", ifelse(is.infinite(UB), "Inf", UB))
  
  data %>%
    mutate(
      clone       = clone_name,
      # 1) start everyone at their true delivery time
      fup_outcome = time.prom.to.delivery,
      
      # 2) too‐early induction: censor at oxytocin time if < LB
      fup_outcome = if_else(
        !is.na(time.prom.oxitocin.hour) & time.prom.oxitocin.hour < LB,
        time.prom.oxitocin.hour,
        fup_outcome
      ),
      
      # 3) too‐late induction: censor at oxytocin time if > UB
      fup_outcome = if_else(
        !is.na(time.prom.oxitocin.hour) & time.prom.oxitocin.hour > UB,
        UB,
        fup_outcome
      ),
      
      # 4) no‐oxytocin & delivery after UB: censor at UB
      fup_outcome = if_else(
        is.na(time.prom.oxitocin.hour) & time.prom.to.delivery > UB,
        UB,
        fup_outcome
      ),
      
      # 5) flag event vs. censor
      outcome = as.integer(fup_outcome == time.prom.to.delivery),
      censor  = as.integer(fup_outcome <  time.prom.to.delivery)
    )
}

# combine all the clones in one data‐frame
df_clones <- map2_dfr(
  brks[-length(brks)],    # LB
  brks[-1],               # UB
  ~ make_clone(df.only.base, .x, .y)
)
#adding the variable for cometing risks:
df_clones <- df_clones%>%
  mutate(
    status = factor(
      case_when(
        censor == 1                   ~ 0L,
        outcome == 1 & NVD == 1       ~ 1L,
        outcome == 1 & NVD == 0       ~ 2L
      ),
      levels = 0:2,
      labels = c("censored","vaginal","other")
    )
  )


#------------------------------------------------------------------------------------------------------------
#Step 3+4:Wighting+estimating CIF and survival curves
#-----------------------------------------------------------------------------------------------------------

#Weighting models
glm_formula_first <- (censor == 0) ~opening_clean_cf+age+nulliparity+preg.week+smoking+GBS
glm_formula_rest <- (censor == 0) ~ ns(t, df =3) +opening_clean_cf +age+nulliparity+preg.week+smoking+GBS
cox_formula <- Surv(t_start, t_stop, censor)~opening_clean_cf+age+nulliparity+preg.week+smoking+GBS



#Weighting loop: In each iteration, we will calculate the weighted and estimated the CIF
#The lopp uses the function " bootstrap_cif_TV_with_CI" which calculates the weights, the estimators and the bootstrap 95% CI
clone_info <- tibble(
  LB    = brks[-length(brks)],
  UB    = brks[-1],
  clone.n = ifelse(is.finite(UB),
                 paste0("window_", LB, "_", UB),
                 paste0("window_", LB, "_Inf"))
)

handlers(global = TRUE)  # Use default progress bar
clone_info_list <- split(clone_info, seq_len(nrow(clone_info)))



## For a parallel calculation 
registerDoParallel(cl)
cl <- makeCluster(10)


results_list <- list()
for (i in seq_along(clone_info_list)) {
  ci <- clone_info_list[[i]]
  clone.n <- ci$clone.n[[1]]
  df_sub  <- df_clones[df_clones$clone == clone.n, , drop = FALSE]
  
  message("Running bootstrap for clone: ", clone.n)
  
  df_clone      =      df_sub
  only.measur.full.time=measur_df
  
  results_list[[i]] <- bootstrap_cif_TV_with_CI(
    df_clone      =      df_sub,
    ci                 = ci,
    only.measur.full.time=measur_df,
    glm_formula_first  = glm_formula_first,
    glm_formula_rest   = glm_formula_rest,
    cox_formula        = cox_formula,
    B                  = 250,
    seed               = 2024 + i  # different seed per clone
  )
}
stopCluster(cl)


#This is the output of the estimated CIFs and survival curves
cif_by_clone <- bind_rows(results_list)

write.csv()
