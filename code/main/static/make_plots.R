library(tidyverse)

#lead the data
cif_plot_df_TV<-read.csv( "G:\\My Drive\\Phd4\\R code\\final_code\\SW\\final\\results\\main\\cif_TV.csv")
cif_plot_df_not_TV<-read.csv( "G:\\My Drive\\Phd4\\R code\\final_code\\SW\\final\\results\\main\\cif_no_TV.csv")

cif_plot_df_TV<-cif_plot_df_TV %>% 
  mutate(TV="yes")

cif_plot_df_not_TV<-cif_plot_df_not_TV %>% 
  mutate(TV="no") %>% 
  filter(method!="No weights")


cif_plot_df<-rbind(cif_plot_df_TV,cif_plot_df_not_TV)
highlight_pts <- cif_plot_df %>% filter(time %in% c(12, 24, 48,72))

h.cox.SW<-highlight_pts %>% filter(TV=="yes",
  method %in% c( "Cox.SW‐model weights")) %>% 
  arrange(event,clone,time) 
  # mutate(
  #   prob = sprintf("%.2f", prob),
  #   cif_lower= sprintf("%.2f",  cif_lower),
  #   cif_upper = sprintf("%.2f", cif_upper)
  # )

#plot unweighted, not_Tv and TV side by side
plot_data <- cif_plot_df %>%
  filter(#clone %in% c("window_0_8", "window_32_Inf"),
         method %in% c("No weights", "Cox.SW‐model weights")) %>%
  mutate(
    # Label treatment strategy
    clone_label = recode(clone,
                         "window_0_8" = "Induction: 0–8h",
                         "window_32_Inf" = "Induction: 32+h"),
    
    # Define method_label using method + TV
    method_label = case_when(
      method == "No weights" ~ "Unweighted",
      method == "Cox.SW‐model weights" & TV == "no" ~ "Weighted\n(Baseline confounders only)",
      method == "Cox.SW‐model weights" & TV == "yes" ~ "Weighted\n(Including time-varying confounders)",
      TRUE ~ NA_character_
    ),
    
    # Recode event label
    event_label = recode(event,
                         "Vaginal delivery" = "Non-assisted vaginal delivery",
                         "Other modes" = "Assisted delivery",
                         "Overall survival" = "Overall survival"),
    
    # Create unique curve ID for grouping
    curve_id = paste(clone_label, method_label, event_label, sep = "_"),
    
    # Optional: cleaner legend
    legend_label = paste(clone_label, "|", method_label, "|", event_label)
  )

plot_data <- plot_data %>%
  mutate(
    event_label = factor(event_label,
                         levels = c("Non-assisted vaginal delivery", "Assisted delivery", "Overall survival"))
  )
ggplot(plot_data, aes(x = time, y = prob,
                      color = clone_label,
                      linetype = method_label,
                      group = curve_id)) +
  geom_step(size = 1) +
  facet_wrap(~ event_label, nrow = 1) +
  labs(
    x = "Time from PROM (hours)",
    y = "Cumulative Incidence",
    linetype = "Estimation Method",
    color = "Treatment Strategy"
  ) +
  scale_linetype_manual(
    name = "Estimation Method",
    values = c(
      "Unweighted" = "dashed",
      "Weighted\n(Baseline confounders only)" = "dotdash",
      "Weighted\n(Including time-varying confounders)" = "solid"
    )
  ) +
  scale_color_manual(
    name = "Treatment Strategy",  # ✅ rename legend title
    values = c("Induction: 0–8h" = "red", "Induction: 32+h" = "blue")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    #legend.justification = "left",      # ⬅️ anchor the box to the left
    #legend.box.just = "left",           # ⬅️ left-align legend items
    #plot.margin = margin(10, -100, 10, 10)  # ⬅️ extra space at bottom
  )+
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper, fill = clone_label),
              alpha = 0.2, color = NA)+
  scale_fill_manual(
    name = "Treatment Strategy",  # ✅ rename legend title
    values = c("Induction: 0–8h" = "red", "Induction: 32+h" = "blue")
  )

ggsave("G:\\My Drive\\Phd4\\R code\\final_code\\SW\\final\\results\\main\\plots\\jpg\\plot_cif_main.jpeg", 
       width = 10, height = 4.5,
       plot = last_plot(),
       units = "in")


#For appendix- include all groups

# Load and prepare full dataset
plot_data_all <- df%>%
  filter(method %in% c("Cox.SW‐model weights", "No weights")) %>%
  mutate(
    clone_label = recode(clone,
                         "window_0_8" = "Pitocin: 0–8h",
                         "window_8_16" = "Pitocin: 8–16h",
                         "window_16_24" = "Pitocin: 16–24h",
                         "window_24_32" = "Pitocin: 24–32h",
                         "window_32_Inf" = "Pitocin: 32+h"),
    method_label = ifelse(method == "Cox.SW‐model weights", "Weighted", "Unweighted"),
    event_label = recode(event,
                         "Vaginal delivery" = "Vaginal",
                         "Other modes" = "Non-vaginal"),
    curve_id = paste(clone_label, method_label, event_label, sep = "_")
  )
ggplot(plot_data_all , aes(x = time, y = cif,
                           color = clone_label,
                           linetype = method_label,
                           group = curve_id)) +
  geom_step(size = 1) +
  facet_wrap(~ event_label, nrow = 1) +
  labs(
    #    title = "Cumulative Incidence by Mode of Delivery, Timing, and Estimation Method",
    x = "Time from PROM (hours)",
    y = "Cumulative Incidence",
    #   color = "Delivery Mode",
    linetype = "Estimation Method"
  ) +
  #scale_color_manual(values = c("Vaginal" = "blue", "Non-vaginal" = "red")) +
  scale_linetype_manual(values = c("Weighted" = "solid", "Unweighted" = "dashed")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = 5, b = 5),
    plot.margin = margin(10, 10, 20, 10)  # top, right, bottom, left
  )+
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper, fill = clone_label),
              alpha = 0.2, color = NA)

ggsave("G:\\My Drive\\Phd4\\overleaf\\plots\\CIF\\SW_CI\\TV\\cif_plot_appendix_tv.pdf", 
       width = 15, height = 8,
       plot = last_plot(),
       units = "in", device = cairo_pdf)

# scale_color_manual(name = "Treatment Strategy", values = c("Pitocin: 0–8h" = "red", "Pitocin: 32+h" = "blue")) +
#scale_fill_manual(name = "Treatment Strategy", values = c("Pitocin: 0–8h" = "red", "Pitocin: 32+h" = "blue"))


# 

# Plot with all treatment groups
# 
# ggplot(plot_data, aes(x = time, y = cif,
#                       color = clone_label,
#                       linetype = method_label,
#                       group = curve_id)) +
#   geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper, fill = clone_label),
#               alpha = 0.2, color = NA) +
#   geom_step(size = 1) +
#   labs(
#     title = "Cumulative Incidence by Treatment, Estimation, and Delivery Mode",
#     x = "Hours since PROM",
#     y = "Cumulative Incidence",
#     color = "Pitocin Timing",
#     fill = "Pitocin Timing",
#     linetype = "Estimation Method"
#   ) +
#   scale_color_manual(values = c(
#     "Pitocin: 0–8h" = "red",
#     "Pitocin: 32+h" = "blue"
#   )) +
#   scale_fill_manual(values = c(
#     "Pitocin: 0–8h" = "red",
#     "Pitocin: 32+h" = "blue"
#   )) +
#   scale_linetype_manual(values = c("Weighted" = "solid", "Unweighted" = "dashed")) +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "bottom") +
#   guides(
#     fill = guide_legend(order = 1),
#     color = guide_legend(order = 1),
#     linetype = guide_legend(order = 2)
#   )
# 
# 
# ggplot(plot_data, aes(
#   x = time, y = cif,
#   color = legend_label,
#   fill = legend_label,
#   group = curve_id
# )) +
#   geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper), alpha = 0.2, color = NA) +
#   geom_step(aes(linetype = legend_label), size = 1) +
#   labs(
#     title = "Cumulative Incidence by Treatment, Estimation, and Delivery Mode",
#     x = "Hours since PROM",
#     y = "Cumulative Incidence",
#     color = "Group",
#     fill = "Group",
#     linetype = "Group"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "bottom")
# 
# 
# 
# #
# ggplot(plot_data, aes(x = time, y = cif,
#                       color = clone_label,
#                       linetype = linetype,
#                       group = curve_id)) +
#   geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper, fill = clone_label),
#               alpha = 0.2, color = NA) +
#   geom_step(size = 1) +
#   labs(
#     title = "Cumulative Incidence by Treatment Timing, Weighting, and Delivery Mode",
#     x = "Hours since PROM",
#     y = "Cumulative Incidence",
#     color = "Pitocin Timing",
#     fill = "Pitocin Timing",
#     linetype = "Estimation Method"
#   ) +
#   scale_color_manual(values = c("Pitocin: 0–8h" = "red", "Pitocin: 32+h" = "blue")) +
#   scale_fill_manual(values = c("Pitocin: 0–8h" = "red", "Pitocin: 32+h" = "blue")) +
#   scale_linetype_manual(values = c("Weighted" = "solid", "Unweighted" = "dashed")) +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "bottom")
