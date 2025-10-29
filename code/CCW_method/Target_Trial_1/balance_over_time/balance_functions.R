# Function to compute SMD for binary
smd_binary <- function(x, treat, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, length(x))
  p1 <- weighted.mean(x[treat == 1], w = weights[treat == 1])
  p0 <- weighted.mean(x[treat == 0], w = weights[treat == 0])
  smd <- (p1 - p0) / sqrt((p1*(1 - p1) + p0*(1 - p0)) / 2)
  return(smd)
}

# Function to compute SMD for continuous
smd_continuous <- function(x, treat, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, length(x))
  m1 <- weighted.mean(x[treat == 1], w = weights[treat == 1])
  m0 <- weighted.mean(x[treat == 0], w = weights[treat == 0])
  s1 <- sqrt(weighted.var(x[treat == 1], w = weights[treat == 1]))
  s0 <- sqrt(weighted.var(x[treat == 0], w = weights[treat == 0]))
  sd_pool <- sqrt((s1^2 + s0^2) / 2)
  return((m1 - m0) / sd_pool)
}

# Define weighted.var manually if not available
weighted.var <- function(x, w) {
  m <- weighted.mean(x, w)
  sum(w * (x - m)^2) / sum(w)
}

# Loop over time and covariates
calc_smd_over_time <- function(data, covariates, time_var, treat_var, weight_var) {
  time_points <- sort(unique(data[[time_var]]))
  
  smd_results <- expand.grid(
    time = time_points,
    covariate = covariates,
    stringsAsFactors = FALSE
  )
  
  smd_results <- smd_results %>%
    rowwise() %>%
    mutate(
      unadj = {
        sub <- data[data[[time_var]] == time, ]
        x <- sub[[covariate]]
        treat <- sub[[treat_var]]
        if (is.numeric(x)) smd_continuous(x, treat) else smd_binary(x, treat)
      },
      adj = {
        sub <- data[data[[time_var]] == time, ]
        x <- sub[[covariate]]
        treat <- sub[[treat_var]]
        w <- sub[[weight_var]]
        if (is.numeric(x)) smd_continuous(x, treat, w) else smd_binary(x, treat, w)
      }
    ) %>% ungroup()
  
  return(smd_results)
}