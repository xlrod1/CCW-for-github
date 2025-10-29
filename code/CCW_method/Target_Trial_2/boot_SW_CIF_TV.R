bootstrap_cif_TV_with_CI <- function(df_clone, ci,only.measur.full.time=NULL,
glm_formula, cox_formula, B = 200, seed = 123) {
  set.seed(seed)

  # 
  # data                  = df_clone
  # LB=LB
  # UB=LB
  # first.w               = is_first_window
  # first.cut.off         = ifelse(is_first_window, UB, NA_real_)
  # only.measur.full.time=only.measur.full.time
  # glm_censor_formula    = glm_formula_use
  # cox_censor_formula    = cox_formula
  # time_grid             = seq(0, 72, by = 0.5)
  # competing.risks       = TRUE
  # 
  # 
  # Point estimate (original)
  original <- compute_survival_for_clone_TV(
    data                  = df_clone,
    only.measur.full.time=only.measur.full.time,
    glm_censor_formula    = glm_formula,
    cox_censor_formula    = cox_formula,
    time_grid             = seq(0, 72, by = 0.5),
    competing.risks       = TRUE
  )
  #B=20
  # Bootstrap in parallel
  boot_list <- foreach(b = 1:B, .packages = c("dplyr", "survival", "tidyr","splines","zoo","slider","data.table"),
                       .export = c("compute_survival_for_clone_TV")) %dopar% {
    boot_data <- df_clone[sample(nrow(df_clone), replace = TRUE), ]
    
    #data=boot_data
    
    compute_survival_for_clone_TV(
      data                  = boot_data,
      only.measur.full.time=only.measur.full.time,
      glm_censor_formula    = glm_formula,
      cox_censor_formula    = cox_formula,
      time_grid             = seq(0, 72, by = 0.5),
      competing.risks       = TRUE
    )
  }
  
  
  boot_results <- bind_rows(boot_list[!sapply(boot_list, is.null)])
  if (!all(c("time", "event", "prob") %in% names(boot_results))) {
    stop("Missing expected columns in boot_results. Possibly all bootstraps failed.")
  }
  # Summarize bootstrap distribution
  ci_summary <- boot_results %>%
    group_by(time, event,method) %>%
    summarise(
      cif_lower = quantile(prob, 0.025, na.rm = TRUE),
      cif_upper = quantile(prob, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Join with point estimate
  output <- original %>%
    left_join(ci_summary, by = c("time", "event","method")) %>%
    mutate(clone = ci)
  
  return(output)
}