
plot_standardized_covariates_trimmed <- function(df, covariates, time_var = "t_start", clone_var = "clone",
                                                 weight_var = NULL, clone_bounds) {
  df_long <- df %>%
    pivot_longer(cols = starts_with("std_"), names_to = "covariate", values_to = "value") %>%
    mutate(covariate = gsub("^std_", "", covariate))
  
  # Filter by clone-specific max time
  df_long <- df_long %>%
    filter(.data[[time_var]] <= clone_bounds[.data[[clone_var]]])
  
  # Unweighted
  unweighted <- df_long %>%
    group_by(.data[[clone_var]], .data[[time_var]], covariate) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), type = "Unweighted", .groups = "drop")
  
  # Weighted
  if (!is.null(weight_var)) {
    weighted <- df_long %>%
      filter(!is.na(.data[[weight_var]])) %>%
      group_by(.data[[clone_var]], .data[[time_var]], covariate) %>%
      summarise(mean_value = weighted.mean(value, .data[[weight_var]], na.rm = TRUE), type = "Weighted", .groups = "drop")
    
    plot_df <- bind_rows(unweighted, weighted)
  } else {
    plot_df <- unweighted
  }
  
  ggplot(plot_df, aes_string(x = time_var, y = "mean_value", color = clone_var)) +
    geom_line() +
    facet_grid(type ~ covariate, scales = "fixed") +
    labs(
      title = "Covariate balance over time by clone (standardized + truncated)",
      x = "Time (hours)",
      y = "Standardized Mean",
      color = "Clone"
    ) +
    theme_minimal()
}
standardize_covariates <- function(df, covariate_names) {
  df_std <- df
  
  for (cov in covariate_names) {
    mean_val <- mean(df[[cov]], na.rm = TRUE)
    sd_val <- sd(df[[cov]], na.rm = TRUE)
    df_std[[paste0("std_", cov)]] <- (df[[cov]] - mean_val) / sd_val
  }
  
  return(df_std)
}


plot_var_mean_over_time_trimmed <- function(df, covariate,
                                            time_var = "t_start",
                                            weight_var = "censoring_weight_cox",
                                            clone_var = "clone") {
  cov_sym <- sym(covariate)
  
  # Grace period end times for each clone
  grace_bounds <- tibble::tibble(
    clone = c("window_0_8", "window_8_16", "window_16_24", "window_24_32", "window_32_Inf"),
    UB = c(8-time.diff, 16-time.diff, 24-time.diff, 32-time.diff, 32)
  )
  
  # Merge and truncate df
  df_trimmed <- df %>%
    left_join(grace_bounds, by = setNames("clone", clone_var)) %>%
    filter(.data[[time_var]] <= UB)
  
  # 1. Unweighted means
  unweighted <- df_trimmed %>%
    group_by(!!sym(clone_var), !!sym(time_var)) %>%
    summarise(
      mean_value = mean(!!cov_sym, na.rm = TRUE),
      type = "Unweighted",
      .groups = "drop"
    )
  
  # 2. Weighted means
  weighted <- df_trimmed %>%
    filter(!is.na(.data[[weight_var]])) %>%
    group_by(!!sym(clone_var), !!sym(time_var)) %>%
    summarise(
      mean_value = weighted.mean(!!cov_sym, .data[[weight_var]], na.rm = TRUE),
      type = "Weighted",
      .groups = "drop"
    )
  
  # 3. Combine
  plot_df <- bind_rows(unweighted, weighted)
  
  #4. Plot
  ggplot(plot_df, aes_string(x = time_var, y = "mean_value", color = clone_var)) +
    geom_line() +
    facet_wrap(~ type, ncol = 1) +
    labs(
      title = paste("Mean of", covariate, "over time by clone"),
      x = "Time (hours)",
      y = paste("Mean", covariate),
      color = "Clone"
    ) +
    theme_minimal()
  
  #return(plot_df)
  
}


plot_multiple_covariates_combined <- function(df, covariates,
                                              time_var = "t_start",
                                              weight_var = "censoring_weight_cox",
                                              clone_var = "clone",
                                              time.diff = 0) {
  grace_bounds <- tibble::tibble(
    clone = c("window_0_8", "window_8_16", "window_16_24", "window_24_32", "window_32_Inf"),
    UB = c(8 - time.diff, 16 - time.diff, 24 - time.diff, 32 - time.diff, 32)
  )
  
  all_plot_df <- purrr::map_dfr(covariates, function(covariate) {
    cov_sym <- rlang::sym(covariate)
    
    df_trimmed <- df %>%
      left_join(grace_bounds, by = setNames("clone", clone_var)) %>%
      filter(.data[[time_var]] <= UB)
    
    # Unweighted
    unweighted <- df_trimmed %>%
      group_by(.data[[clone_var]], .data[[time_var]]) %>%
      summarise(
        mean_value = mean(!!cov_sym, na.rm = TRUE),
        type = "Unweighted",
        .groups = "drop"
      )
    
    # Weighted
    weighted <- df_trimmed %>%
      filter(!is.na(.data[[weight_var]])) %>%
      group_by(.data[[clone_var]], .data[[time_var]]) %>%
      summarise(
        mean_value = weighted.mean(!!cov_sym, .data[[weight_var]], na.rm = TRUE),
        type = "Weighted",
        .groups = "drop"
      )
    
    bind_rows(unweighted, weighted) %>%
      mutate(covariate = factor(covariate, levels = covariates))  # << FIXED
  })
  
  # Create a combined facet label
  all_plot_df <- all_plot_df %>%
    mutate(facet_label = paste(type, covariate, sep = ": "))
  covariate_labels <- c(
    "opening_clean_cf" = "Cervical Dilation (cm)",
    "height_clean_cf" = "Fetal Station (cm)"
  )
  
  all_plot_df <- all_plot_df %>%
    mutate(
      cov_label = recode(covariate, !!!covariate_labels),
      facet_label = paste(type, cov_label, sep = ": ")
    )
  clone_levels <- c("window_0_8", "window_8_16", "window_16_24", "window_24_32", "window_32_Inf")
  
  all_plot_df <- all_plot_df %>%
    mutate(!!clone_var := factor(.data[[clone_var]], levels = clone_levels))
  
  clone_levels <- c("window_0_8", "window_8_16", "window_16_24", "window_24_32", "window_32_Inf")
  clone_labels <- c("0–8", "8–16", "16–24", "24–32", "32+")
  
  
  ggplot(all_plot_df, aes_string(x = time_var, y = "mean_value", color = clone_var)) +
    geom_line() +
    facet_wrap(~ facet_label, scales = "free_y", ncol = length(unique(all_plot_df$covariate))) +
    labs(
  #    title = "Covariate trends over time by treatment window",
      x = "Time (hours)",
      y = "Mean value",
      color = "Treatment Strategy"
    ) +
    scale_color_manual(
      values = RColorBrewer::brewer.pal(5, "Set1"),
      breaks = clone_levels,
      labels = clone_labels
    ) +
    theme_minimal(base_size = 14)+
    theme(legend.position = "bottom")
  
  
}