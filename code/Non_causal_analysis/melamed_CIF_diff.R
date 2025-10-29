library(riskRegression)
library(dplyr)
library(ggplot2)
library(tidyr)

# Time grid for prediction
prediction_times <- seq(0, 120, by = 1)

# Empty list to store all predicted CIFs
cif_curves <- list()
#interval_label<-"0–8"
# Loop over intervals
for (interval_label in unique(pseudodata$interval)) {
  message("Processing interval: ", interval_label)
  
  # Prepare data
  data_interval <- pseudodata %>%
    filter(interval == interval_label) %>%
    mutate(
      treatment = factor(treatment, levels = c(0, 1)),
      event_type = case_when(
        cs == 1 ~ 2,
        cs == 0 ~ 1,
        TRUE ~ 0
      )
    )
  
  if (length(unique(data_interval$treatment)) < 2) next
  
  # Drop constant covariates
  vars <- c("treatment", "age", "nulliparity","preg.week", "smoking")
  var_levels <- sapply(data_interval[, vars], function(x) length(unique(x)))
  covariates <- names(var_levels[var_levels > 1 & names(var_levels) != "treatment"])
  model_formula <- as.formula(paste(
    "Hist(delivery_time, event_type) ~ treatment",
    if (length(covariates) > 0) paste("+", paste(covariates, collapse = " +")) else ""
  ))
  
  # Fit model
  tryCatch({
    csc_model <- CSC(model_formula, data = data_interval)
    
    # Create newdata with standard covariates
    newdata <- data.frame(treatment = factor(c(0, 1), levels = c(0, 1)))
    for (v in covariates) {
      if (is.numeric(data_interval[[v]])) {
        newdata[[v]] <- mean(data_interval[[v]], na.rm = TRUE)
      } else {
        newdata[[v]] <- as.character(sort(unique(data_interval[[v]]))[1])
      }
    }
    
    # Predict CIFs for vaginal delivery (cause = 1)
    pred <- predict(csc_model, newdata = newdata, times = prediction_times,
                     cause = 1, se = TRUE, band = TRUE, type =  "absRisk")
    
    
    # Extract predictions for each group
    group_labels <- c("expectant", "induced")
    # Extract and reshape predictions for both groups
    pred_list <- lapply(1:2, function(i) {
      data.frame(
        time = pred$times,
        CIF = pred$absRisk[i, ],
        lower = pred$absRisk.lower[i, ],
        upper = pred$absRisk.upper[i, ],
        se   = pred$absRisk.se[i, ],  # ← ADD THIS
        group = group_labels[i],
        interval = interval_label
      )
    })
    
    df_long <- bind_rows(pred_list)
    cif_curves[[interval_label]] <- df_long
    
    
  }, error = function(e) {
    message("Skipping interval ", interval_label, " due to error: ", e$message)
  })
}

# Combine all intervals
cif_all_1<- bind_rows(cif_curves)

# Define desired order
interval_order <- c("0–8", "8–16", "16–24", "24–32", "32–36")

# Force interval to be an ordered factor
cif_all_1 $interval <- factor(cif_all_1$interval, levels = interval_order)


#same but with cause=2
# Empty list to store all predicted CIFs
cif_curves <- list()
#interval_label<-"0–8"
# Loop over intervals
for (interval_label in unique(pseudodata$interval)) {
  message("Processing interval: ", interval_label)
  
  # Prepare data
  data_interval <- pseudodata %>%
    filter(interval == interval_label) %>%
    mutate(
      treatment = factor(treatment, levels = c(0, 1)),
      event_type = case_when(
        cs == 1 ~ 2,
        cs == 0 ~ 1,
        TRUE ~ 0
      )
    )
  
  if (length(unique(data_interval$treatment)) < 2) next
  
  # Drop constant covariates
  vars <- c("treatment", "age", "nulliparity","preg.week", "smoking")
  var_levels <- sapply(data_interval[, vars], function(x) length(unique(x)))
  covariates <- names(var_levels[var_levels > 1 & names(var_levels) != "treatment"])
  model_formula <- as.formula(paste(
    "Hist(delivery_time, event_type) ~ treatment",
    if (length(covariates) > 0) paste("+", paste(covariates, collapse = " +")) else ""
  ))
  
  # Fit model
  tryCatch({
    csc_model <- CSC(model_formula, data = data_interval)
    
    # Create newdata with standard covariates
    newdata <- data.frame(treatment = factor(c(0, 1), levels = c(0, 1)))
    for (v in covariates) {
      if (is.numeric(data_interval[[v]])) {
        newdata[[v]] <- mean(data_interval[[v]], na.rm = TRUE)
      } else {
        newdata[[v]] <- as.character(sort(unique(data_interval[[v]]))[1])
      }
    }
    
    # Predict CIFs for vaginal delivery (cause = 1)
    pred <- predict(csc_model, newdata = newdata, times = prediction_times,
                    cause = 2, se = TRUE, band = TRUE, type =  "absRisk")
    
    
    # Extract predictions for each group
    group_labels <- c("expectant", "induced")
    # Extract and reshape predictions for both groups
    pred_list <- lapply(1:2, function(i) {
      data.frame(
        time = pred$times,
        CIF = pred$absRisk[i, ],
        lower = pred$absRisk.lower[i, ],
        upper = pred$absRisk.upper[i, ],
        se   = pred$absRisk.se[i, ],  # ← ADD THIS
        group = group_labels[i],
        interval = interval_label
      )
    })
    
    df_long <- bind_rows(pred_list)
    cif_curves[[interval_label]] <- df_long
    
    
  }, error = function(e) {
    message("Skipping interval ", interval_label, " due to error: ", e$message)
  })
}

# Combine all intervals
cif_all_2 <- bind_rows(cif_curves)

# Force interval to be an ordered factor
cif_all_2 $interval <- factor(cif_all_2$interval, levels = interval_order)



cif_all_1$cause=1
cif_all_2$cause=2
cif_all<-rbind(cif_all_1,cif_all_2)


# Compute composite (cause = 3)
cif_composite <- cif_all %>%
  filter(cause %in% c(1, 2)) %>%
  group_by(time, group, interval) %>%
  summarise(
    CIF = sum(CIF),
    se = sqrt(sum(se^2)),  # Assuming independence
    .groups = "drop"
  ) %>%
  mutate(
    cause = 3,
    lower = CIF - qnorm(0.975) * se,
    upper = CIF + qnorm(0.975) * se
  )
cif_all_full <- bind_rows(cif_all, cif_composite)
cif_all_full <- cif_all_full %>%
  mutate(cause_label = factor(cause,
                              levels = c(1, 2, 3),
                              labels = c("Non-Assisted", "Assisted", "Composite")
  ))


# Compute CIF difference + CI for each cause + interval
cif_diff <- cif_all %>%
  filter(group %in% c("expectant", "induced")) %>%
  select(time, group, CIF, se, interval, cause) %>%
  group_by(cause, interval, time) %>%
  pivot_wider(names_from = group, values_from = c(CIF, se)) %>%
  mutate(
    diff    = CIF_induced - CIF_expectant,
    se_diff = sqrt(se_induced^2 + se_expectant^2),
    z       = qnorm(0.975),
    lower   = diff - z * se_diff,
    upper   = diff + z * se_diff
  ) %>%
  ungroup()

# Plot
# Define colors
interval_colors <- c(
  "0–8"   = "red",
  "8–16"  = "gray70",
  "16–24" = "gray70",
  "24–32" = "gray70",
  "32–36" = "blue"
)

ggplot(cif_diff, aes(x = time, y = diff, color = interval, fill = interval)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ cause, labeller = as_labeller(c(`1` = "Non assisted vaginal delivery", `2` = "Assisted delivery"))) +
  scale_color_manual(values = interval_colors, drop = FALSE) +
  scale_fill_manual(values = interval_colors, drop = FALSE) +
  labs(
   # title = "Difference in CIF (Induced − Expectant) by Cause and Interval",
    x = "Time from PROM (hours)",
    y = "Cumulative Incidence Difference",
    color = "Induction Interval",
    fill = "Induction Interval"
  ) +
  theme_minimal(base_size = 14)+
  theme(legend.position = "bottom")  # ⬅️ Legend below the plot



# ggplot(cif_all, aes(x = time, y = CIF, color = group, linetype = group, fill = group)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
#   facet_wrap(~ interval) +
#   labs(
#     #  title = "Adjusted CIF for Vaginal Delivery by Induction Interval",
#     x = "Time interval from PROM (hours)",
#     y = "Cumulative Incidence function",
#     color = "Group", linetype = "Group", fill = "Group"
#   ) +
#   theme_minimal()


write.csv(cif_diff ,"cif_diff_all_cause.csv",row.names=FALSE)
write.csv(cif_all ,"cif_all_cause.csv",row.names=FALSE)

ggsave("cif_all_cause.jpeg", 
       width = 10, height = 4.5,
       plot = last_plot(),
       units = "in")

#Make good table 
# 1. Define target time points and desired row order
summary_times <- c(12, 24, 48, 72)
treatment_order <- c("0–8", "8–16", "16–24", "24–32", "32–36")
outcome_order <- c("Non-Assisted", "Assisted")

# 2. Format CIF estimates with 2 decimal digits
cif_table <- cif_diff %>%
  filter(time %in% summary_times, cause %in% c(1, 2)) %>%
  mutate(
    time = paste0(time, "h"),
    Outcome = recode(cause, `1` = "Non-Assisted", `2` = "Assisted"),
    interval = factor(interval, levels = treatment_order),
    Outcome = factor(Outcome, levels = outcome_order),
    estimate_str = sprintf("%.2f%% (%.2f%%, %.2f%%)",
                           100 * diff,
                           100 * lower,
                           100 * upper)
    
  ) %>%
  select(Outcome, interval, time, estimate_str) %>%
  pivot_wider(
    names_from = time,
    values_from = estimate_str
  ) %>%
  arrange(Outcome, interval)

# 3. View the result
View(cif_table)  # or print(cif_table)
write.csv(cif_table, "cif_summary_table.csv", row.names = FALSE)
