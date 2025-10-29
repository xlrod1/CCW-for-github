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



source("boot_SW_CIF_TV.R")
source("cloning_functions.R")
source("f_compute_survival_for_clone_SW_TV_no_glm.R")
df.only.base<- read.csv("data_base_line.csv") 
measur_df<-fread("opening_measurment.csv") 

df.only.base<-df.only.base %>% right_join(.,measur_df %>% select(id) %>% distinct()) 



#Crating the first time each women reached opening if 4. this can be changed
c_o<-4
first.time.open<-measur_df%>%
  filter(opening_clean_cf >=c_o) %>%
  group_by(id) %>%
  summarize(first_time_reach_c_o = min(time.prom.measurement)) %>%
  ungroup()
df.only.base<- df.only.base %>% left_join(first.time.open)


#for now we will do complete case analysis
list.of.confounder<-age ~nulliparity+preg.week+smoking+GBS
vars <- all.vars(list.of.confounder)
df.only.base<-df.only.base %>%
  # select only the vars in the formula, then drop any rows with NA
  filter(
    if_all(all_of(vars), ~ !is.na(.))
  )




#------------------------------------------------------------------------------------------------------------------
#Step 1: cloning+censoring. We use the function make_clone_strategy
#------------------------------------------------------------------------------------------------------------------

#These are the treatemnt windows which can be changed
brks<- c(0,12,Inf)
brks.e <- c(0,12)



# build all the clones in one data‐frame
df_clones <- purrr::map_dfr(1:4, ~ make_clone_strategy(df.only.base, strategy = .x, c_0 = 4, T = 12))
df_clones <-clone_censor (df.only.base, treatment_threshold = 12) 


epsilon <- 1e-4
df_clones<- df_clones %>%
  mutate(fup_outcome = if_else(fup_outcome == 0, epsilon, fup_outcome)) %>% 
  rename(clone=clone_id)


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


#---------------------------------------------------------------------------------------------------------------
#Step 3+4: weighting+estimation of CID and survival curves


glm_formula <- (censor == 0) ~ ns(t, df =3) +opening_clean_cf+age+nulliparity+preg.week+smoking+GBS
cox_formula <- Surv(t_start, t_stop, censor)~opening_clean_cf+age+nulliparity+preg.week+smoking+GBS


#For parallel computation
handlers(global = TRUE)  # Use default progress bar
cl <- makeCluster(10)
registerDoParallel(cl)

clone_info_list <- unique(df_clones$clone)
results_list <- list()
#i=1
for (i in seq_along(clone_info_list)) {
  ci <- clone_info_list[[i]]
  #clone.n <- ci$clone.n[[1]]
  df_sub  <- df_clones[df_clones$clone == ci, , drop = FALSE]
  
  message("Running bootstrap for clone: ", ci)
  
  df_clone      =      df_sub
  only.measur.full.time=measur_df
  
  results_list[[i]] <- bootstrap_cif_TV_with_CI(
    df_clone      =      df_sub,
    ci                 = ci,
    only.measur.full.time=measur_df,
   # glm_formula_first  = glm_formula_first,
    glm_formula= glm_formula,
    cox_formula        = cox_formula,
    B                  = 250,
    seed               = 2024 + i  # different seed per clone
  )
}

stopCluster(cl)


# The output: the estimated CIF under each treatment
cif_by_clone <- bind_rows(results_list)











