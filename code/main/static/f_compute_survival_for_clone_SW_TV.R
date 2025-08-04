library(survival)
library(ggplot2)
library(dplyr)
library(splines)
library(data.table)
library(slider)
library(zoo)
# Function: compute_survival_for_clone
# Description: Computes survival probabilities for a given clone dataset, including 
# Kaplan-Meier estimates with different weighting methods.
#
# Parameters:
# - data: A dataframe containing 'fup_outcome', 'outcome', 'censor', and relevant covariates.
#
# Returns:
# - A dataframe with estimated survival probabilities for different methods.
#data<-boot_data

compute_survival_for_clone_TV <- function(data,
                                          LB,
                                          UB,
                                           first.w=FALSE,
                                           first.cut.off=NA,
                                           #event_times,
                                           only.measur.full.time=NULL,
                                           glm_censor_formula,
                                           cox_censor_formula,
                                           time_grid =seq(0, 72, by = 0.5),
                                           competing.risks=FALSE
                                           ) {
 
  data$id.new<-seq(1:nrow(data))
   #work only with complete cases:when i'll have the mesurment it will be differents
  # if(is.null(only.measur.full.time)){
  #   
  #   
  #   vars <- all.vars(glm_censor_formula)[!all.vars(glm_censor_formula)=="t"]
  #   data<-data %>%
  #     # select only the vars in the formula, then drop any rows with NA
  #     filter(
  #       if_all(all_of(vars), ~ !is.na(.))
  #     )
  #   
  # }
  #data=jj
  max.fup<-max(data$fup_outcome)
  event_times<-seq(0,max.fup,by=0.5)
  
  # 1) Expand to long via survSplit
  #data=jj %>% filter(clone=="window_0_12")
  sd0 <- data %>% mutate(t_start = 0)
  
  long_evt   <- survSplit(sd0, cut = event_times, end = "fup_outcome",
                          start = "t_start", event = "outcome", id = "ID")
  long_cens  <- survSplit(sd0, cut = event_times, end = "fup_outcome",
                          start = "t_start", event = "censor",  id = "ID")
  long_evt$censor <- long_cens$censor
  long_evt$t_stop <- long_evt$fup_outcome
  df_long <- long_evt %>% arrange(ID) %>% 
    mutate(t=t_start)
  
  #f<-df_long %>% filter(id==1)
  #jj<-data_long %>% filter(id==2)
  # 2) join your measurement‐times and fill any missing weights- by id not ID
  if(!is.null(only.measur.full.time)){
    #only.measur.full.time
    # Step 1: Create time intervals
    # Ensure valid internal self-reference
    only.measur.full.time <- as.data.table(only.measur.full.time)
    #ff<-  only.measur.full.time %>% filter(id==1)
    # Step 1: Create time intervals
    only.measur.full.time[, time_bin_start := floor(time.prom.measurement * 2) / 2]
    only.measur.full.time[, time_bin_end := time_bin_start + 0.5]
    
    #Step 2: Keep one measurement per bin (e.g., latest per bin)
    only.measur.filtered <- only.measur.full.time[
      order(id,time.prom.measurement),  # sort before slicing
      .SD[.N],  # keep last (latest) per bin
      by = .(id, time_bin_start)
    ]
    
    only.measur.filtered <- only.measur.filtered %>%
      rename(t_start = time_bin_start) 
    
    df_long_with_tv <- df_long %>%
      left_join(only.measur.filtered, by = c("id", "t_start"))
    
    df_long_with_tv <- df_long_with_tv %>%
      arrange(ID, t_start) %>%
      group_by(ID) %>%
      mutate(
        opening_clean_cf = na.locf(opening_clean_cf, na.rm = FALSE)
        #,
        #height_clean_cf = na.locf(height_clean_cf, na.rm = FALSE)
      )
   # if I wnat to bulid a rolling weights of the itme depdendent variables 
    # df_long_with_tv <- df_long_with_tv %>%
    #   arrange(id, t_start) %>%
    #   group_by(id) %>%
    #   mutate(
    #     opening_max_3h = slide_dbl(opening_clean_cf, ~max(.x, na.rm = TRUE),
    #                                .before = 6 - 1, .complete = FALSE)  # 6*0.5 = 3 hours
    #   )
    
    
    df_long<-df_long_with_tv %>% ungroup()
    
   # f<-df_long %>% filter(id==1)
  }
  
  #assign df_long in the Global envirment- this is for the Cox predict
 # assign("df_long", df_long, envir = .GlobalEnv)
  
  # 3) fit the two censoring models. 
  # 1. Create a private environment that ONLY contains df_long
  
  local_env <- new.env(parent = baseenv())  # <-- not emptyenv(), allows base functions
  
  local_env$Surv <- survival::Surv  # 🔴 This is the key fix
  local_env$df_long <- df_long
  
  
  
  fit_cox  <- coxph(cox_censor_formula,
                    data  =df_long,
                    ties  = "efron")
  
   
  if(!is.null(only.measur.full.time)){
    cox_censor_num<- update(cox_censor_formula, . ~ .-opening_clean_cf)
    
    fit_cox_n<-coxph(cox_censor_num,
                     data  =df_long,
                     ties  = "efron")
  }
  
  if(is.null(only.measur.full.time)){
    cox_censor_num<- update(cox_censor_formula, . ~ 1)
    
    fit_cox_n<-coxph(cox_censor_num,
                     data  =df_long,
                     ties  = "efron")
  }
  
  
  
  # summary(fit_cox)
  # summary(summary(fit_cox))
  # Inject df_long into the environment the model will look in
  # 3. Inject the custom environment into BOTH the formula and terms objects
  environment(fit_cox$formula) <- local_env
  environment(fit_cox$terms) <- local_env
  environment(fit_cox_n$formula) <- local_env
  environment(fit_cox_n$terms) <- local_env
  
  expected_vec = predict( fit_cox, type = "expected")
  expected0 = predict( fit_cox_n, type = "expected")
  
  # 3. Now build your expected.data without any scoping magic- this are the predcitiond from teh cox model
  expected.data <- df_long %>%
    select(ID,t_start) %>%
    mutate(expected = expected_vec,
           expected0=expected0) %>% 
            distinct()
  
 # jj<-expected.data %>% filter(id==1)
  # 4) linear predictors & predictions
  if (first.w==TRUE){
    #if its the first window than the logistic regresion shoukd be estimated only on the the cebsoring point,
    #all other pints shoud get a weiught of 1
    #Second, we will create weights using logistoc regression
    #first.cut.off=4
    data.at.f<-df_long %>%filter(t_stop==first.cut.off)
    # estimate teh censoring weights only at time=12.its a simple logistic regression because the censoring can accure only at one time point
    
    
    if(!is.null(only.measur.full.time)){
      glm_censor_formula_num<- update(glm_censor_formula, . ~ . -opening_clean_cf)
      fit_glm0  <- glm(glm_censor_formula_num,
                       family = binomial(link = "logit"),
                       data   = data.at.f)
      
    }
    
    if(is.null(only.measur.full.time)){
       
      glm_censor_formula_num<- update(glm_censor_formula, . ~  1)
      fit_glm0  <- glm(glm_censor_formula_num,
                       family = binomial(link = "logit"),
                       data   = data.at.f)
      
    }
    
    fit_glm  <- glm(glm_censor_formula,
                    family = binomial(link = "logit"),
                    data   = data.at.f)
    #glm_censor_forlmua_num<- update(glm_censor_formula, . ~ .-opening_clean_cf -height_clean_cf)
    # fit_glm0  <- glm(glm_censor_forlmua_num,
    #                 family = binomial(link = "logit"),
    #                 data   = data.at.f)
    #print(summary(fit_glm))
    #predication data
    
    pred.data<-data.at.f %>% select(ID,t_start) %>% 
      mutate(pred=predict(fit_glm , type = "response",newdata=data.at.f),
             pred0=predict(fit_glm0 , type = "response",newdata=data.at.f)) %>% 
      distinct()
    
    
    #join the prediction from the Cox and the log model- want also the SW and not SW  
    df_long <- df_long  %>%
      left_join(.,pred.data) %>% 
      left_join(.,expected.data) %>% 
      group_by(ID) %>% 
      mutate(cumhaz = cumsum(expected),
             S_hat_cox = exp(-cumhaz),
             W_cox = 1 / S_hat_cox,
            prob_uncensor_log= ifelse(t_stop!=first.cut.off,1,pred),
             W_t_log=1/prob_uncensor_log,
            cumhaz_SW = cumsum(expected-expected0),
            S_hat_cox_SW = exp(-cumhaz_SW),
            W_cox_SW = 1 / S_hat_cox_SW,
            prob_uncensor_log_SW= ifelse(t_stop!=first.cut.off,1,pred/pred0),
            W_t_log_SW=1/prob_uncensor_log_SW) %>% 
      mutate(W_log= cumprod( W_t_log),
             W_log_SW= cumprod( W_t_log_SW))  %>% 
      mutate(
        censoring_weight_log     = W_log,
        censoring_weight_cox = W_cox ,
        censoring_weight_log_SW     = W_log_SW,
        censoring_weight_cox_SW = W_cox_SW ,
      )
    
    # Truncate censoring weights at 1st and 99th percentiles
    truncate_weights <- function(w, lower = 0.01, upper = 0.99) {
      q <- quantile(w, probs = c(lower, upper), na.rm = TRUE)
      pmax(pmin(w, q[2]), q[1])
    }
    
    df_long <- df_long %>%
      ungroup() %>%  # remove ID grouping to avoid per-ID truncation
      mutate(
        censoring_weight_log      = truncate_weights(censoring_weight_log),
        censoring_weight_cox      = truncate_weights(censoring_weight_cox),
        censoring_weight_log_SW   = truncate_weights(censoring_weight_log_SW),
        censoring_weight_cox_SW   = truncate_weights(censoring_weight_cox_SW)
      )
    
    # jj<-df_long %>% filter(!is.na(pred))
    # jj2<-df_long %>% filter(X==4467)
      #4467
  }else{
    
    # Define censoring times: early (t_start < LB) or late (t_start == UB), unless UB is Inf
    if (is.finite(UB)) {
      data_censor_model <- df_long %>%
        filter(t_start < LB | t_start == UB)
    } else {
      data_censor_model <- df_long %>%
        filter(t_start < LB)  # Only early treatment matters
    }
    
    
    
    if(!is.null(only.measur.full.time)){
      glm_censor_formula_num<- update(glm_censor_formula, . ~ . -opening_clean_cf -height_clean_cf)
      fit_glm0  <- glm(glm_censor_formula_num,
                       family = binomial(link = "logit"),
                       data   = data_censor_model)
      
    }
    
    if(is.null(only.measur.full.time)){
      
      glm_censor_formula_num<- update(glm_censor_formula, . ~  1)
      fit_glm0  <- glm(glm_censor_formula_num,
                       family = binomial(link = "logit"),
                       data   = data_censor_model)
      
    }
    
    
    
    
    fit_glm  <- glm(glm_censor_formula,
                    family = binomial(link = "logit"),
                    data   = data_censor_model)
    
    
    #summary(fit_glm)
    df_long <-df_long %>% 
      left_join(.,expected.data) %>% 
      mutate(pred = if_else(t_start < LB | (is.finite(UB) & t_start == UB),
                            predict(fit_glm,  type = "response", newdata = df_long),
                            1),
             pred0 = if_else(t_start < LB | (is.finite(UB) & t_start == UB),
                             predict(fit_glm0, type = "response", newdata = df_long),
                             1)
    ) %>% #this the logistic prediction
      group_by(ID) %>%
      mutate(cumhaz = cumsum(expected),
             S_hat_cox = exp(-cumhaz),
             W_cox = 1 / S_hat_cox,
             prob_uncensor_log= pred,
             W_t_log=1/prob_uncensor_log,
             cumhaz_SW = cumsum(expected-expected0),
             S_hat_cox_SW = exp(-cumhaz_SW),
             W_cox_SW = 1 / S_hat_cox_SW,
             prob_uncensor_log_SW= pred/pred0,
             W_t_log_SW=1/prob_uncensor_log_SW) %>% 
      mutate(W_log= cumprod( W_t_log),
             W_log_SW= cumprod( W_t_log_SW))  %>% 
      mutate(
        censoring_weight_log     = W_log,
        censoring_weight_cox = W_cox ,
        censoring_weight_log_SW     = W_log_SW,
        censoring_weight_cox_SW = W_cox_SW ,
      )
    
    
    # Truncate censoring weights at 1st and 99th percentiles
    truncate_weights <- function(w, lower = 0.01, upper = 0.99) {
      q <- quantile(w, probs = c(lower, upper), na.rm = TRUE)
      pmax(pmin(w, q[2]), q[1])
    }
    
    df_long <- df_long %>%
      ungroup() %>%  # remove ID grouping to avoid per-ID truncation
      mutate(
        censoring_weight_log      = truncate_weights(censoring_weight_log),
        censoring_weight_cox      = truncate_weights(censoring_weight_cox),
        censoring_weight_log_SW   = truncate_weights(censoring_weight_log_SW),
        censoring_weight_cox_SW   = truncate_weights(censoring_weight_cox_SW)
      )
    
  }
 
  #df.complete<- df_long %>% filter(!is.na(censoring_weight_log),!is.na(censoring_weight_cox))
  if (competing.risks==FALSE){
    df.complete<-df_long
  # 5) build weighted survival curves
  km_w_logit <- survfit(Surv(t_start, t_stop, outcome) ~ 1,
                        data    = df.complete,
                        weights = df.complete$censoring_weight_log)
  km_w_cox   <- survfit(Surv(t_start, t_stop, outcome) ~ 1,
                        data    = df.complete,
                        weights = df.complete$censoring_weight_cox)
  km_unw     <- survfit(Surv(t_start, t_stop, outcome) ~ 1,
                        data    = df.complete)
  # 6) extract on a fine grid
  time_grid <-time_grid
  
  # first get the summaries
  sum_logit <- summary(km_w_logit, times = time_grid, extend = TRUE)
  sum_cox   <- summary(km_w_cox,   times = time_grid, extend = TRUE)
  sum_unw   <- summary(km_unw,     times = time_grid, extend = TRUE)
  
  # then build your table, including n.risk
  surv_dat <- tibble::tibble(
    time        = time_grid,
    logit_w     = sum_logit$surv,
    logit_risk  = sum_logit$n.risk,
    cox_w       = sum_cox$surv,
    cox_risk    = sum_cox$n.risk,
    unweighted  = sum_unw$surv,
    unw_risk    = sum_unw$n.risk
  )
  }
  
  if (competing.risks) {
    # 1) build a 3‐level “status”:
    #    0 = censored
    #    1 = vaginal delivery  (NVD==1 & outcome==1)
    #    2 = other modes       (NVD==0 & outcome==1)
   
    df2 <- df_long %>%
      mutate(
        status = case_when(
          # 1) vaginal delivery event on that interval
          outcome == 1 & NVD == 1  ~ 1L,
          # 2) competing (CS/instrumental) on that interval
          outcome == 1 & NVD == 0  ~ 2L,
          # 3) otherwise (including final‐interval censoring!), stay at risk → 0
          TRUE                      ~ 0L
        )
      )
    #jj<-df2 %>% filter(is.na(status))
    # convert to factor so survfit() auto-detects mstate
    df2 <- df2 %>%
      mutate(
        status_fac = factor(
          status,
          levels = 0:2,
          labels = c("censored", "vaginal", "other")
        )
      )
    # 2) fit three A–J curves via survfit(); “id” clusters by subject
    # 2) fit three A–J curves via survfit(); “id” clusters by subject
    fit_unw <- survfit(
      Surv(t_start, t_stop, status_fac) ~ 1,
      data    = df2,
      id      = ID
    )
    
    fit_cox <- survfit(
      Surv(t_start, t_stop, status_fac) ~ 1,
      data    = df2,
      id      = ID,
      weights = df2$censoring_weight_cox
    )
    fit_log <- survfit(
      Surv(t_start, t_stop,status_fac) ~ 1,
      data    = df2,
      id      = ID,
      weights = df2$censoring_weight_log
    )
    fit_cox_SW <- survfit(
      Surv(t_start, t_stop, status_fac) ~ 1,
      data    = df2,
      id      = ID,
      weights = df2$censoring_weight_cox_SW
    )
    fit_log_SW <- survfit(
      Surv(t_start, t_stop,status_fac) ~ 1,
      data    = df2,
      id      = ID,
      weights = df2$censoring_weight_log_SW
    )
    
    
    # 1. Extract summaries at *native time points*
    times<-time_grid
    # sum_u <- summary(fit_unw, times = fit_unw$time, extend = TRUE)
    # sum_c <- summary(fit_cox, times = fit_cox$time, extend = TRUE)
    # sum_l <- summary(fit_log, times = fit_log$time, extend = TRUE)
    
    sum_u <- summary(fit_unw, times = times, extend = TRUE)
    sum_c <- summary(fit_cox, times = times, extend = TRUE)
    sum_l <- summary(fit_log, times = times, extend = TRUE)
    sum_c_SW <- summary(fit_cox_SW, times = times, extend = TRUE)
    sum_l_SW <- summary(fit_log_SW, times = times, extend = TRUE)
    
    
    # 2. Turn each into a tidy data frame
    extract_cif_df <- function(sum_obj, method_label) {
      tibble(
        time = sum_obj$time,
        survival = sum_obj$pstate[, 1],      # 🟢 overall survival
        
        vaginal = sum_obj$pstate[, 2],
        competing = sum_obj$pstate[, 3]
      ) %>%
        pivot_longer(cols = c(survival, vaginal, competing),
                     names_to = "event", values_to = "prob") %>%
        mutate(method = method_label)
    }
    
    df_unw <- extract_cif_df(sum_u, "No weights")
    df_cox <- extract_cif_df(sum_c, "Cox‐model weights")
    df_log <- extract_cif_df(sum_l, "Logistic‐model weights")
    df_cox_SW <- extract_cif_df(sum_c_SW, "Cox.SW‐model weights")
    df_log_SW <- extract_cif_df(sum_l_SW, "Logistic.SW‐model weights")
    
    # 3. Combine everything
    surv_dat_long <- bind_rows(df_unw, df_cox, df_log,df_cox_SW,df_log_SW)
   
    
    return(surv_dat_long)
  }
  
  
  
  
  
  return(surv_dat)
}

# Function: plot_survival_curves
# Description: Plots Kaplan-Meier survival curves for all clones and estimation methods.
#
# Parameters:
# - survival_list: A list of survival dataframes (one per clone).
#
# Returns:
# - A ggplot2 object with survival curves.

plot_survival_curves <- function(survival_list) {
  survival_combined <- bind_rows(survival_list, .id = "clone_id")
  
  ggplot(survival_combined, aes(x = time)) +
   # geom_line(aes(y = manual, color = "Manual Weights")) +
    geom_line(aes(y = log, color = "Log Model Weights")) +
    geom_line(aes(y = cox, color = "Cox Model Weights")) +
    geom_line(aes(y = not_w, color = "No Weights")) +
    facet_wrap(~clone_id) +
    labs(title = "Survival Curves for Different Clones", x = "Time", y = "Survival Probability") +
    theme_minimal() +
    scale_color_manual(values = c("Manual Weights" = "blue", "Log Model Weights" = "red", "Cox Model Weights" = "green", "No Weights" = "black"))
}

# Example usage:
# survival_clone1 <- compute_survival_for_clone(clone_1_data)
# survival_clone2 <- compute_survival_for_clone(clone_2_data)
# survival_clone3 <- compute_survival_for_clone(clone_3_data)
# survival_clone4 <- compute_survival_for_clone(clone_4_data)

# survival_list <- list(clone_1 = survival_clone1, clone_2 = survival_clone2, clone_3 = survival_clone3, clone_4 = survival_clone4)
# plot_survival_curves(survival_list)


# jj.c<-df2 %>% select(id,NVD) %>% distinct() 
# jj.c<-df2n.c %>% select(id,NVD,censor) 
# jj<-df2%>% filter(id==20906)
# 
# table(jj.c$NVD)
