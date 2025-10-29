base_path <- "balance\\"

comapre_0_8<-read.csv( file.path(base_path, "compare_0_8_compare.csv"))
comapre_8_16<-read.csv( file.path(base_path, "compare_8_16_set.csv"))
comapre_16_24<-read.csv( file.path(base_path, "compare_16_24_set.csv"))
comapre_24_32<-read.csv( file.path(base_path, "compare_24_32_set.csv"))



#PLOT
#truncate until 8 
comapre_0_8 <-comapre_0_8%>% filter(t_start<7.5)

smd_long <- comapre_0_8 %>%
  pivot_longer(
    cols = starts_with("SMD_"),
    names_to   = c("type", "covariate"),
    names_pattern = "SMD_(unw|w)_(.*)",
    values_to  = "SMD"
  ) %>%
  mutate(
    type = ifelse(type == "unw", "Unweighted", "Weighted"),
    # optional: nicer covariate labels
    covariate = recode(covariate,
                       `preg.week` = "Pregnancy week",
                       `opening_clean_cf` = "Cervical dilation",
                       `nulliparity` = "Nulliparity",
                       `GBS` = "GBS status",
                       .default = covariate
    )
  )
# choose the exact order of columns (left → right)
comp_order <- c(
  "window_8_16 vs window_0_8",
  "window_16_24 vs window_0_8",
  "window_24_32 vs window_0_8",
  "window_32_Inf vs window_0_8"
)

# pretty labels for the strip headers (optional)
comp_labels <- c(
  "window_8_16 vs window_0_8"   = "0–8 vs 8–16",
  "window_16_24 vs window_0_8"  = "0–8 vs 16–24",
  "window_24_32 vs window_0_8"  = "0–8 vs 24–32",
  "window_32_Inf vs window_0_8" = "0–8 vs 32+"
)

smd_long <- smd_long %>%
  mutate(comparison = factor(comparison, levels = comp_order))


p<-ggplot(smd_long, aes(x = t_start, y = abs(SMD),
                     linetype = type)) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50") +
#  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.9, alpha = 0.9, color = "black") +
  facet_grid(covariate ~ comparison,
             scales = "free_y",
             labeller = labeller(comparison = comp_labels),
             ) +
  labs(x = "Time from PROM (hour)", y = "|SMD|",
       linetype = "Estimation") +
  scale_linetype_manual(values = c("Unweighted" = "dashed",
                                   "Weighted"   = "solid")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        strip.text.y = element_text(angle = 0, face = "bold"),
        strip.text.x = element_text(face = "bold"),
        panel.spacing.y = unit(1.5, "lines")
        ) # increase vertical space between rows)
p
# Save as PNG
ggsave(
  filename = "SMD_plot.jpeg",
  plot = p,
  path = "balance",
  width = 12,
  height = 8,
  dpi = 300)



#PLOT
#truncate until 8 
comapre_8_16 <-comapre_8_16 %>% filter(t_start<15.5)

smd_long <- comapre_8_16 %>%
  pivot_longer(
    cols = starts_with("SMD_"),
    names_to   = c("type", "covariate"),
    names_pattern = "SMD_(unw|w)_(.*)",
    values_to  = "SMD"
  ) %>%
  mutate(
    type = ifelse(type == "unw", "Unweighted", "Weighted"),
    # optional: nicer covariate labels
    covariate = recode(covariate,
                       `preg.week` = "Pregnancy week",
                       `opening_clean_cf` = "Cervical dilation",
                       `nulliparity` = "Nulliparity",
                       `GBS` = "GBS status",
                       .default = covariate
    )
  )
# choose the exact order of columns (left → right)
comp_order <- c(
  "window_8_16 vs window_16_24",
  "window_8_16 vs window_24_32",
  "window_8_16 vs window_32_Inf"
)

# pretty labels for the strip headers (optional)
comp_labels <- c(
  
  "window_8_16 vs window_16_24" = "8-16 vs 16–24",
  "window_8_16 vs window_24_32" = "8-16 vs 24–32",
  "window_8_16 vs window_32_Inf"= "8-16 vs 32+"
)

smd_long <- smd_long %>%
  mutate(comparison = factor(comparison, levels = comp_order))

unique( smd_long $comparison)


p<-ggplot(smd_long, aes(x = t_start, y = abs(SMD),
                        linetype = type)) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50") +
  #  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.9, alpha = 0.9, color = "black") +
  facet_grid(covariate ~ comparison,
   #          scales = "free_y",
             labeller = labeller(comparison = comp_labels)) +
  labs(x = "Time from PROM (hour)", y = "|SMD|",
       linetype = "Estimation") +
  scale_linetype_manual(values = c("Unweighted" = "dashed",
                                   "Weighted"   = "solid")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        strip.text.y = element_text(angle = 0, face = "bold"),
        strip.text.x = element_text(face = "bold"),
        panel.spacing.y = unit(1.5, "lines"))
# Save as PNG
ggsave(
  filename = "SMD_plot_comapre_8_16.jpeg",
  plot = p,
  path = "balance",
  width = 12,
  height = 8,
  dpi = 300)



#PLOT
#truncate until 8 
comapre_16_24 <-comapre_16_24 %>% filter(t_start<23.5)

smd_long <-comapre_16_24%>%
  pivot_longer(
    cols = starts_with("SMD_"),
    names_to   = c("type", "covariate"),
    names_pattern = "SMD_(unw|w)_(.*)",
    values_to  = "SMD"
  ) %>%
  mutate(
    type = ifelse(type == "unw", "Unweighted", "Weighted"),
    # optional: nicer covariate labels
    covariate = recode(covariate,
                       `preg.week` = "Pregnancy week",
                       `opening_clean_cf` = "Cervical dilation",
                       `nulliparity` = "Nulliparity",
                       `GBS` = "GBS status",
                       .default = covariate
    )
  )
# choose the exact order of columns (left → right)
comp_order <- c(
  "window_16_24 vs window_24_32",
  "window_16_24 vs window_32_Inf"
)

# pretty labels for the strip headers (optional)
comp_labels <- c(
  
  "window_16_24 vs window_24_32" = "16-24 vs 24–32",
  "window_16_24 vs window_32_Inf"= "16-24 vs 32+"
)

smd_long <- smd_long %>%
  mutate(comparison = factor(comparison, levels = comp_order))

unique( smd_long $comparison)


p<-ggplot(smd_long, aes(x = t_start, y = abs(SMD),
                        linetype = type)) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50") +
  #  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.9, alpha = 0.9, color = "black") +
  facet_grid(covariate ~ comparison,
#             scales = "free_y",
             labeller = labeller(comparison = comp_labels)) +
  labs(x = "Time from PROM (hour)", y = "|SMD|",
       linetype = "Estimation") +
  scale_linetype_manual(values = c("Unweighted" = "dashed",
                                   "Weighted"   = "solid")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        strip.text.y = element_text(angle = 0, face = "bold"),
        strip.text.x = element_text(face = "bold"),
        panel.spacing.y = unit(1.5, "lines"))
# Save as PNG
ggsave(
  filename = "SMD_plot_comapre_16_24.jpeg",
  plot = p,
  path = "balance",
  width = 12,
  height = 8,
  dpi = 300)




comapre_24_32 <-comapre_24_32%>% filter(t_start<31.5)

smd_long <-comapre_24_32%>%
  pivot_longer(
    cols = starts_with("SMD_"),
    names_to   = c("type", "covariate"),
    names_pattern = "SMD_(unw|w)_(.*)",
    values_to  = "SMD"
  ) %>%
  mutate(
    type = ifelse(type == "unw", "Unweighted", "Weighted"),
    # optional: nicer covariate labels
    covariate = recode(covariate,
                       `preg.week` = "Pregnancy week",
                       `opening_clean_cf` = "Cervical dilation",
                       `nulliparity` = "Nulliparity",
                       `GBS` = "GBS status",
                       .default = covariate
    )
  )
# choose the exact order of columns (left → right)
comp_order <- c(
 
  "window_24_32 vs window_32_Inf"
)

# pretty labels for the strip headers (optional)
comp_labels <- c(
  
 
  "window_24_32 vs window_32_Inf"= "24-32 vs 32+"
)

smd_long <- smd_long %>%
  mutate(comparison = factor(comparison, levels = comp_order))

unique( smd_long $comparison)


p<-ggplot(smd_long, aes(x = t_start, y = abs(SMD),
                        linetype = type)) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50") +
  #  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.9, alpha = 0.9, color = "black") +
  facet_grid(covariate ~ comparison,
 #            scales = "free_y",
             labeller = labeller(comparison = comp_labels)) +
  labs(x = "Time from PROM (hour)", y = "|SMD|",
       linetype = "Estimation") +
  scale_linetype_manual(values = c("Unweighted" = "dashed",
                                   "Weighted"   = "solid")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        strip.text.y = element_text(angle = 0, face = "bold"),
        strip.text.x = element_text(face = "bold"),
        panel.spacing.y = unit(1.5, "lines"))
p
# Save as PNG
ggsave(
  filename = "SMD_plot_comapre_24_32.jpeg",
  plot = p,
  path = "balance",
  width = 12,
  height = 8,
  dpi = 300)
