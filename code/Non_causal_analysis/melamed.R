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
library(broom)
library(tibble)
df.only.base<- read.csv("data_base_line.csv") 

df<-df.only.base %>% 
  select(id,
         time.prom.to.delivery,
         time.prom.oxitocin.hour,
         NVD,
         age ,parity,GWG,smoking,GBS.pos,preg.week,GBS) %>% 
  mutate(spontaneous_labor_time=ifelse(is.na(time.prom.oxitocin.hour),time.prom.to.delivery,NA),
         cs=1-NVD,
         nulliparity=ifelse(parity==0,1,0))

brks <- c(0, 8, 16, 24, 32, 36)
labels <- c("0–8", "8–16", "16–24", "24–32", "32–36")


#Create Risk Sets by Interval
# Clean variable names
df_clean <- df %>%
  rename(
    delivery_time = time.prom.to.delivery,
    induction_time = time.prom.oxitocin.hour
  )
pseudodata_list <- list()

for (i in 1:(length(brks) - 1)) {
  t_start <- brks[i]
  t_end <- brks[i + 1]
  label <- labels[i]
  
  subset_df <- df_clean %>%
    filter(delivery_time > t_end,
           is.na(induction_time) | induction_time > t_start,
           is.na(spontaneous_labor_time) | spontaneous_labor_time > t_start) %>%
    mutate(treatment = ifelse(!is.na(induction_time) & induction_time > t_start & induction_time <= t_end, 1, 0),
           interval = label)
  
  pseudodata_list[[i]] <- subset_df
}

# Combine all intervals
pseudodata <- bind_rows(pseudodata_list)



#Step 1: Estimate Risk Difference / Relative Risk per Interval
results_rr <- pseudodata %>%
  group_by(interval) %>%
  group_split() %>%
  lapply(function(subset_df) {
    if (length(unique(subset_df$treatment)) < 2) return(NULL)  # skip if no variation
    print(table(subset_df$treatment, subset_df$nulliparity))
    model <- glm(cs ~ treatment + age + nulliparity+preg.week+smoking,
                 family = binomial(link = "log"), data = subset_df)
    
    broom::tidy(model, conf.int = TRUE) %>%
      filter(term == "treatment") %>%
      mutate(interval = unique(subset_df$interval))
  }) %>%
  bind_rows()


results_rr <- results_rr %>%
  mutate(RR = exp(estimate),
         RR_low = exp(conf.low),
         RR_high = exp(conf.high))


interval_order <- c("0–8", "8–16", "16–24", "24–32", "32–36")
results_rr$interval <- factor(results_rr$interval, levels = interval_order)
ggplot(results_rr, aes(x = interval, y = RR)) +
  geom_point() +
  geom_errorbar(aes(ymin = RR_low, ymax = RR_high), width = 0.1) +
  labs(x = "Time Interval from PROM (hours)",
       y = "Risk Ratio (Induction vs Expectant) for assisted delivery"
    #   title = "Effect of Induction by Time Since PROM") +
  )+
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  )

write.csv(results_rr,"RR_table.csv",row.names=FALSE)

library(quantreg)

results_time <- pseudodata %>%
  group_by(interval) %>%
  group_split() %>%
  lapply(function(subset_df) {
    if (length(unique(subset_df$treatment)) < 2) return(NULL)
    
    model <- rq(delivery_time ~ treatment + age + nulliparity+preg.week+smoking,
                tau = 0.5, data = subset_df)
    
    broom::tidy(model, conf.int = TRUE) %>%
      filter(term == "treatment") %>%
      mutate(interval = unique(subset_df$interval))
  }) %>%
  bind_rows()


results_time <- results_time %>%
  filter(term == "treatment") %>%
  mutate(
    interval = factor(interval, levels = unique(interval)),
    median_diff = estimate,
    lower_ci = conf.low,
    upper_ci = conf.high
  )

interval_order <- c("0–8", "8–16", "16–24", "24–32", "32–36")
results_time $interval <- factor(results_time $interval, levels = interval_order)


ggplot(results_time, aes(x = interval, y = median_diff)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    x = "Time Interval from PROM (hours)",
    y = "Difference (Induction vs Expectant) in Median Time to Delivery (hours)"
  #  title = "Effect of Induction on Time to Delivery by Interval"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  )



# #####################################Cometing risks####################

#First plot if hazard ration of delivery in each induction group
library(riskRegression)
library(dplyr)
library(broom)
library(ggplot2)

# Storage for results
hr_results <- list()

# Loop over each interval
for (interval_label in unique(pseudodata$interval)) {
  message("Processing interval: ", interval_label)

  # Prepare data
  data_interval <- pseudodata %>%
    filter(interval == interval_label) %>%
    mutate(
      treatment = factor(treatment, levels = c(0, 1)),
      event_type = case_when(
        cs == 1 ~ 2,  # cesarean
        cs == 0 ~ 1,  # vaginal
        TRUE ~ 0      # censored
      )
    )

  # Skip if not enough treatment variation
  if (length(unique(data_interval$treatment)) < 2) {
    message("  Skipped: only one treatment group.")
    next
  }

  # Skip if no vaginal deliveries
  if (sum(data_interval$event_type == 1) == 0) {
    message("  Skipped: no vaginal deliveries.")
    next
  }

  # Drop constant covariates
  vars <- c("treatment", "age", "nulliparity", "preg.week", "smoking")
  var_levels <- sapply(data_interval[, vars], function(x) length(unique(x)))
  covariates <- names(var_levels[var_levels > 1 & names(var_levels) != "treatment"])

  # Build formula
  model_formula <- as.formula(paste(
    "Hist(delivery_time, event_type) ~ treatment",
    if (length(covariates) > 0) paste("+", paste(covariates, collapse = " +")) else ""
  ))

  # Define model data
  model_vars <- c("delivery_time", "event_type", "treatment", covariates)
  data_model <- data_interval[, model_vars, drop = FALSE] %>% na.omit()

  if (nrow(data_model) == 0) {
    message("  Skipped: no complete data after NA removal.")
    next
  }

  # Fit model
  tryCatch({
    model <- CSC(model_formula, data = data_model)

    if (!"models" %in% names(model) || !"Cause 2" %in% names(model$models)) {
      message("  Skipped: cause-specific model for vaginal delivery not found.")
      next
    }

    model1 <- model$models$`Cause 2`

    if (!inherits(model1, "coxph")) {
      message("  Skipped: model1 is not a Cox model.")
      next
    }

    coef_table <- summary(model1)$coefficients

    if (!"treatment1" %in% rownames(coef_table)) {
      message("  Skipped: treatment1 not in coefficient table.")
      next
    }

    hr_est <- coef_table["treatment1", "exp(coef)"]
    hr_lci <- exp(coef_table["treatment1", "coef"] - 1.96 * coef_table["treatment1", "se(coef)"])
    hr_uci <- exp(coef_table["treatment1", "coef"] + 1.96 * coef_table["treatment1", "se(coef)"])

    hr_results[[interval_label]] <- data.frame(
      interval = interval_label,
      HR = hr_est,
      HR_lower = hr_lci,
      HR_upper = hr_uci
    )
  }, error = function(e) {
    message("  Skipped: CSC model failed — ", e$message)
  })
}

# Combine results
hr_df <- bind_rows(hr_results)
# Define desired order
interval_order <- c("0–8", "8–16", "16–24", "24–32", "32–36")

# Force interval to be an ordered factor
hr_df$interval <- factor(hr_df$interval, levels = interval_order)
# Plot HRs with confidence intervals
ggplot(hr_df, aes(x = interval, y = HR)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = HR_lower, ymax = HR_upper), width = 0.1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
  labs(
  #  title = "Adjusted Hazard Ratio for Vaginal Delivery by Induction Interval",
    x = "Time interval from PROM (hours)",
    y = "Hazard Ratio (Induced vs Expectant) for assited delivery"
  )  +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  )



write.csv(hr_df,"hr_df_table_non_vaginal.csv",row.names=FALSE)
ggsave("hr_plot_non_vaginal.jpeg",
       width = 10, height = 7,
       plot = last_plot(),
       units = "in")









