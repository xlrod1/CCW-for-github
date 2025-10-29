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
60*250


source("f_compute_survival_for_clone_SW_TV_b.R")


df.only.base<- read.csv("data_base_line_new.csv") 
measur_df<-fread("opening_height2.csv") 
measur_df<-measur_df %>% select(-height_clean_cf) %>% ###take only the opening data
  mutate(opening_clean_cf=ifelse(time.prom.measurement==0&is.na(opening_clean_cf),0,opening_clean_cf)
  ) %>% 
  filter(!is.na(opening_clean_cf))

df.only.base<-df.only.base %>% right_join(.,measur_df %>% select(id) %>% distinct()) %>% 
  mutate(GBS.pos=ifelse(is.na(GBS)|GBS==0,0,1),
         nulliparity=ifelse(parity==0,1,0))

#for now we will do complete cases analysis
list.of.confounder<-age ~nulliparity+preg.week+smoking+GBS


vars <- all.vars(list.of.confounder)
df.only.base<-df.only.base %>%
  filter(
    if_all(all_of(vars), ~ !is.na(.))
  )

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

# build all the clones in one data‐frame
df_clones <- map2_dfr(
  brks[-length(brks)],    # LB
  brks[-1],               # UB
  ~ make_clone(df.only.base, .x, .y)
)

# jj<-df_clones %>% filter(id==2)


#add the variabel foe competenf risks:
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




glm_formula_first <- (censor == 0) ~opening_clean_cf+age+nulliparity+preg.week+smoking+GBS
glm_formula_rest <- (censor == 0) ~ ns(t, df =3) +opening_clean_cf +age+nulliparity+preg.week+smoking+GBS
cox_formula <- Surv(t_start, t_stop, censor)~opening_clean_cf+age+nulliparity+preg.week+smoking+GBS


clone_info <- tibble(
  LB    = brks[-length(brks)],
  
  
  UB    = brks[-1],
  clone.n = ifelse(is.finite(UB),
                   paste0("window_", LB, "_", UB),
                   paste0("window_", LB, "_Inf"))
)

handlers(global = TRUE)  # Use default progress bar


clone_info_list <- split(clone_info, seq_len(nrow(clone_info)))


start=Sys.time()
cl <- makeCluster(10)
registerDoParallel(cl)
# Not in parallel — just loop over windows
results_list <-data_list<- list()
#i=1
for (i in seq_along(clone_info_list)) {
  ci <- clone_info_list[[i]]
  clone.n <- ci$clone.n[[1]]
  df_sub  <- df_clones[df_clones$clone == clone.n, , drop = FALSE]
  
  message("Running clone: ", clone.n)
  
  df_clone      =      df_sub
  only.measur.full.time=measur_df
  
  
  LB <- ci$LB[[1]]
  UB <- ci$UB[[1]]
  is_first_window  <- (LB == 0)
  glm_formula_use  <- if (is_first_window) glm_formula_first else glm_formula_rest


  # Point estimate (original)
  result <- compute_survival_for_clone_TV_balance(
    data                  = df_clone,
    LB=LB,
    UB=UB,
    first.w               = is_first_window,
    first.cut.off         = ifelse(is_first_window, UB, NA_real_),
    only.measur.full.time=measur_df,
    glm_censor_formula    = glm_formula_use,
    cox_censor_formula    = cox_formula,
    time_grid             = seq(0, 72, by = 0.5),
    competing.risks       = TRUE
  )
  
  results_list[[i]]<-result$surv_dat_long
  data_list[[i]]<-result$data
}



stop=Sys.time()
stopCluster(cl)

print(min(stop-start))
# row-bind them all into one data frame
cif_by_clone <- bind_rows(results_list)
data_by_clone <- bind_rows(data_list)

