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
  vars <- c("treatment", "age", "parity", "GBS.pos", "smoking")
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
cif_all <- bind_rows(cif_curves)
# Define desired order
interval_order <- c("0–8", "8–16", "16–24", "24–32", "32-36")

# Force interval to be an ordered factor
cif_all $interval <- factor(cif_all $interval, levels = interval_order)
# Plot
ggplot(cif_all, aes(x = time, y = CIF, color = group, linetype = group, fill = group)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  facet_wrap(~ interval) +
  labs(
  #  title = "Adjusted CIF for Vaginal Delivery by Induction Interval",
    x = "Time interval from PROM (hours)",
    y = "Cumulative Incidence function",
    color = "Group", linetype = "Group", fill = "Group"
  ) +
  theme_minimal()


write.csv(cif_all ,"G:\\My Drive\\Phd4\\R code\\final_code\\SW\\final\\code\\alternative_methods\\cif_all_not_vaginal.csv",row.names=FALSE)
ggsave("G:\\My Drive\\Phd4\\R code\\final_code\\SW\\final\\code\\alternative_methods\\cif_not_vaginal.jpeg", 
       width = 10, height = 4.5,
       plot = last_plot(),
       units = "in")

