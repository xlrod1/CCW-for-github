# Function to compute over time for a pair
compute_smd_overtime <- function(dat,
                                 pair,
                                 time_var = "t_start",
                                 weight_var = "censoring_weight_cox_SW",
                                 covars) {
  d <- dat %>%
    dplyr::filter(.data[["clone"]] %in% pair) %>%
    mutate(clone = factor(.data[["clone"]], levels = pair))  # preserve order
  
  d %>%
    group_by(.data[[time_var]]) %>%
    summarise(
      across(all_of(covars),
             ~ smd_num(.x, g = .data[["clone"]]),
             .names = "SMD_unw_{.col}"),
      across(all_of(covars),
             ~ smd_num(.x, g = .data[["clone"]], w = .data[[weight_var]]),
             .names = "SMD_w_{.col}"),
      .groups = "drop"
    ) %>%
    mutate(comparison = paste(pair, collapse = " vs "))
}


# Numeric/binary SMD function
smd_num <- function(x, g, w = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  g <- droplevels(factor(g))
  i1 <- g == levels(g)[1]; i0 <- g == levels(g)[2]
  m1 <- weighted.mean(x[i1], w[i1], na.rm=TRUE)
  m0 <- weighted.mean(x[i0], w[i0], na.rm=TRUE)
  v1 <- sum(w[i1] * (x[i1] - m1)^2, na.rm=TRUE) / sum(w[i1], na.rm=TRUE)
  v0 <- sum(w[i0] * (x[i0] - m0)^2, na.rm=TRUE) / sum(w[i0], na.rm=TRUE)
  (m1 - m0) / sqrt((v1 + v0) / 2)
}


plot_smd_comparison <- function(df_long, comparison_value,
                                abs_smd = TRUE,
                                free_y = FALSE) {
  
  dfp <- df_long %>% filter(comparison == comparison_value)
  
  p <- ggplot(
    dfp,
    aes(x = t_start, y = if (abs_smd) abs(SMD) else SMD, color = type)
  ) +
    geom_hline(yintercept = 0.1, linetype = "dashed") +
    geom_hline(yintercept = 0.2, linetype = "dashed") +
    geom_line(linewidth = 0.9, alpha = 0.9) +
    facet_wrap(~ covariate, scales = if (free_y) "free_y" else "fixed") +
    labs(
      x = "Time (t_start, hours)",
      y = if (abs_smd) "|SMD|" else "SMD",
      color = "Estimation",
      title = paste0("Balance over time: ", comparison_value)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
  
  p
}

# --- Helpers (only if not already defined) ---
upper_bound_from_clone <- function(cl) {
  if (grepl("Inf$", cl)) return(Inf)
  nums <- as.numeric(unlist(regmatches(cl, gregexpr("\\d+", cl))))
  nums[2]
}
cap_time_for_pair <- function(pair) min(sapply(pair, upper_bound_from_clone))

pair_cap <- function(dat, pair, time_var = "t_start") {
  cap <- cap_time_for_pair(pair)        # your helper from earlier
  dat %>% dplyr::filter(.data[[time_var]] < cap)
}
# small runner
run_pairs <- function(pairs, dat, covars) {
  lapply(pairs, function(p)
    compute_smd_overtime(
      dat        = dat,
      pair       = p,
      covars     = covars,
      time_var   = "t_start",
      weight_var = "censoring_weight_cox_SW",
      cap_time   = cap_time_for_pair(p)
    )
  ) %>% bind_rows()
}
