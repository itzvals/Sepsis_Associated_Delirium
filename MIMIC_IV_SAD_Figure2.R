#########################################################################################################
# ZHANG PAPER REPRODUCTION
# FIGURE 2A/B ICU Mortality Kaplan-Meier Curves and Boxplot
#########################################################################################################

# Clean work space
rm(list = ls())

# Load packages
library(survival)
library(survminer)
library(ggplot2)
library(data.table)
library(cowplot)

# Load finalized SAD cohort (post MICE) 
load("sepsis_cohort_complete.RData")  

# Rename DF
sepsis_cohort <- complete_cohort 

# Check number patients/variables
cat(sprintf("✓ Loaded: %d patients, %d variables\n", 
            nrow(sepsis_cohort), ncol(sepsis_cohort)))

# Convert to data.table if not already
setDT(sepsis_cohort)

# Ensure datetime columns are in correct format
sepsis_cohort[, icu_admit_time := as.POSIXct(icu_admit_time)]
sepsis_cohort[, icu_discharge_time := as.POSIXct(icu_discharge_time)]
sepsis_cohort[, death_time := as.POSIXct(death_time)]

#########################################################################################################
# CALCULATE ICU DEATH RELATED VARIABLES
#########################################################################################################

# Calculate days to death from ICU admission
sepsis_cohort[, days_to_death := as.numeric(
  difftime(death_time, icu_admit_time, units = "days")
)]

# Calculate ICU length of stay
if(!"icu_los" %in% names(sepsis_cohort)) {
  sepsis_cohort[, icu_los := as.numeric(
    difftime(icu_discharge_time, icu_admit_time, units = "days")
  )]
}

#########################################################################################################
# CREATE VARIABLES FOR 28-DAY ICU MORTALITY
#########################################################################################################

# Event indicator: Death in ICU within 28 days
sepsis_cohort[, icu_mortality_28d := as.integer(
  !is.na(death_time) &                          
    death_time <= icu_discharge_time &            # Death occurred in ICU
    days_to_death <= 28                    
)]

# Time variable that follows all patients for full 28 days
sepsis_cohort[, time_28d_icu := pmin(
  ifelse(icu_mortality_28d == 1, days_to_death, 28),  # If patient died in the ICU: use death time, else 28
  28,                                                   
  na.rm = TRUE
)]

# Safety checks for time variable
sepsis_cohort[is.na(time_28d_icu) | time_28d_icu < 0, time_28d_icu := 0.01]
sepsis_cohort[time_28d_icu > 28, time_28d_icu := 28]

# Summary statistics
cat(sprintf("Total patients: %d\n", nrow(sepsis_cohort)))
cat(sprintf("Deaths in ICU within 28 days: %d (%.1f%%)\n", 
            sum(sepsis_cohort$icu_mortality_28d), 
            100 * mean(sepsis_cohort$icu_mortality_28d)))
cat(sprintf("\nBy group:\n"))
sepsis_cohort[, .(
  N = .N,
  Deaths = sum(icu_mortality_28d),
  Mortality_Rate = sprintf("%.1f%%", 100 * mean(icu_mortality_28d))
), by = sad_group]

#########################################################################################################
# SURVIVAL ANALYSIS - 28-DAY ICU MORTALITY
#########################################################################################################

# Create survival object
surv_obj_icu <- Surv(time  = sepsis_cohort$time_28d_icu,
                     event = sepsis_cohort$icu_mortality_28d)

# Fit survival model
fit_icu <- survfit(surv_obj_icu ~ sad_group, data = sepsis_cohort)

# Statistical test (Log-rank test)
surv_diff_icu <- survdiff(surv_obj_icu ~ sad_group, data = sepsis_cohort)
# View
surv_diff_icu

# Calculate p-value
pval_icu <- 1 - pchisq(surv_diff_icu$chisq, length(surv_diff_icu$n) - 1)

# View Log rank test output
cat(sprintf("Chi-square = %.2f\n", surv_diff_icu$chisq))
cat(sprintf("P-value = %.4g\n", pval_icu))

#########################################################################################################
# FIGURE 2A - KAPLAN-MEIER SURVIVAL CURVES
#########################################################################################################

# Plot
fig_icu <- ggsurvplot(
  fit_icu,
  data = sepsis_cohort,
  conf.int = TRUE,                    
  pval = FALSE,                      
  risk.table = TRUE,                  
  risk.table.height = 0.25,          
  risk.table.title = "Number at risk",
  xlim = c(0, 28),                    
  break.time.by = 7,                  
  censor = FALSE,
  legend.title = "",                  
  legend.labs = c("Non-SAD", "SAD"),  
  palette = c("#0d9cf2", "#f0b70f"),
  xlab = "Time (days)",
  ylab = "Survival Probability",
  ggtheme = theme_bw(base_size = 16)
)

# Adjust p value
fig_icu$plot <- fig_icu$plot +
  scale_y_continuous(limits = c(0.5, 1.0), 
                     breaks = seq(0.5, 1.0, 0.1)) +
  # Custom log-rank test annotation
  annotate("text", 
           x = 1,              
           y = 0.62,           
           label = paste0("Log-rank test\nP < 0.01"),
           hjust = 0.2,         
           vjust = 1,         
           size = 5)        

# View the plot
print(fig_icu)

# Save 
png("Figure2A_28d_ICU_Mortality.png", width = 2400, height = 2100, res = 300)
print(fig_icu)
dev.off()


################################################################################
# Figure 2B: ICU Length of Stay Boxplots
################################################################################

# Create 50-day capped version for visualization 
sepsis_cohort[, icu_los_50d := pmin(icu_los, 50)]

# Mann-Whitney U test - use full ICU LOS data
wilcox_result <- wilcox.test(
  icu_los ~ sad_group,
  data = sepsis_cohort,
  exact = FALSE
)
# View results
print(wilcox_result)

# Summary statistics
los_stats <- sepsis_cohort[, .(
  N = .N,
  median = median(icu_los, na.rm = TRUE),
  q25 = quantile(icu_los, 0.25, na.rm = TRUE),
  q75 = quantile(icu_los, 0.75, na.rm = TRUE),
  mean = mean(icu_los, na.rm = TRUE),
  sd = sd(icu_los, na.rm = TRUE)
), by = sad_group]

# View length of stay by group
print(los_stats)

# Create p-value label
if(wilcox_result$p.value < 0.001) {
  pval_label <- "P < 0.001"
} else if(wilcox_result$p.value < 0.01) {
  pval_label <- "P < 0.01"
} else {
  pval_label <- sprintf("P = %.3f", wilcox_result$p.value)
}

# Create boxplot with 50-day cap
fig_boxplot <- ggplot(sepsis_cohort, 
                aes(x = sad_group, y = icu_los_50d, fill = sad_group)) +
  geom_boxplot(
    outlier.shape = 16,
    outlier.size = 1,
    outlier.alpha = 0.3,
    width = 0.6
  ) +
  scale_fill_manual(values = c("Non-SAD" = "#0d9cf2", "SAD" = "#f0b70f")) +
  scale_y_continuous(
    limits = c(0, 55),
    breaks = seq(0, 50, 10)
  ) +
  scale_x_discrete(labels = c("Non-SAD" = "Non-SAD", "SAD" = "SAD")) + 
  labs(
    x = "",
    y = "ICU Length of Stay (days)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  annotate(
    "text",
    x = 2.35,
    y = 52,
    label = pval_label,
    size = 5
  )

print(fig_boxplot)

# Save
ggsave(
  "Figure2B_LOS_boxplot.png", fig_boxplot, width = 8, height = 6, dpi = 300)

################################################################################
# Combine figures 2A/B into one figure (Report)
################################################################################

# Extract components from Figure 2A
survival_plot <- fig_icu$plot
risk_table_plot <- fig_icu$table

# Combine survival plot with risk table vertically
fig2a_complete <- plot_grid(
  survival_plot,
  risk_table_plot,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(3, 1)
)

# Combine Figure 2A (with risk table) and Figure 2B horizontally
fig2_combined <- plot_grid(
  fig2a_complete,
  fig_boxplot,
  labels = c("A", "B"),
  label_size = 18,
  label_fontface = "bold",
  ncol = 2,
  rel_widths = c(1.2, 1)
)

# Display
print(fig2_combined)

# Save combined figure
ggsave(
  "Figure2_Combined.png", fig2_combined, width = 14, height = 7,dpi = 300)
