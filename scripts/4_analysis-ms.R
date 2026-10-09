# Libraries

library(glmmTMB)
library(ggplot2)
library(dplyr)
library(cowplot)
library(here)
library(data.table)
library(patchwork)
library(stringr)
library(colorspace)
library(scales)
library(dplyr)
library(ggplot2)
library(scales)
library(reshape2)
library(tidyverse)
library(tidyr)
library(flextable)
library(corrplot)

# Colours

cols <- c(
  Insecta   = "#F2C400",
  Mammalia  = "#8B4513",
  Aves      = "#ACA0DC",
  Squamata  = "#98FB98",
  Arachnida = "#2B0078",
  Amphibia  = "#41B6C4"
)


# Working directory
setwd("/workdir/")

# setwd(here("model-simulation_march_2026"))
# resf=NULL
# for(i in 1:100){
#   load(paste0(i,"_model.RData"))
#   resf=rbind(resf,as.data.table(res_simdf))
# }
# 
# save(resf, file="model-simulation-join.RData")

load("model-simulation-join.RData")

# Filter non-convergent models

resf$conv = 1
resf$conv[which(
  resf$convergence == FALSE |
    is.na(resf$aic) |
    is.na(resf$mu_zero_estimated_se) | resf$hessian_pos_def == FALSE)] = 0
resf_filtered = subset(resf, conv == 1)
resf_filtered$r = round(resf_filtered$r, 2)
resf_filtered$phi = round(resf_filtered$phi, 2)

# Create db for 5 reliability metrics

stat_df = resf_filtered %>%
  filter(between(r_estimated, -2, 2), duration < 40, sdev < 2) %>%
  # between(phi_estimated, -.9, .9),
  # between(sdev_estimated, 0, 1.5)) %>%
  group_by(duration, sdev, phi) %>%
  summarise(
    bias = mean((r_estimated - r)),
    error = mean(sqrt(((
      r_estimated - r
    )) ^ 2), na.rm = TRUE),
    bias_phi = mean((phi_estimated - phi)),
    bias_sdev_sim = mean((sdev_sim - sdev)),
    bias_phi_sim = mean((phi_sim - phi)),
    bias_sdev = mean((sdev_estimated - sdev)),
    bias_phi2 = mean(abs(phi_estimated - phi)),
    bias_sdev2 = mean(abs(sdev_estimated - sdev)),
    n_pow = sum(r == (-0.02), na.rm = TRUE),
    power = sum(slope_pval[r == (-0.02) & r_estimated < 0] < 0.05, na.rm = TRUE) / sum(r == (-0.02)),
    power2 = sum(slope_pval[r == (-0.2) & r_estimated < 0] < 0.05, na.rm = TRUE) / sum(r == (-0.2)),
    false_pos = sum(slope_pval[r == 0] < 0.05, na.rm = TRUE) / sum(r == 0, na.rm = TRUE),
    n_fp = sum(r == 0.00, na.rm = TRUE),
    type_M=mean(abs(r_estimated[slope_pval<0.05 & r!=0])/abs(r[slope_pval<0.05 & r!=0])),
    #type_M = mean(abs(r_estimated[slope_pval < 0.05 & abs(r) > 1e-6] / r[slope_pval < 0.05 & abs(r) > 1e-6]), na.rm = TRUE),
    type_S = sum(slope_pval < 0.05 & r != 0 & sign(r) != sign(r_estimated), na.rm = TRUE) / sum(slope_pval < 0.05 & r != 0, na.rm = TRUE),
    n_models = length(model_id)
  )


fwrite(stat_df, here("errors_simulations.csv"))
stat_df_or= stat_df
####

# Pipeline to work with the db directly

stat_df_or =fread("CODE/data/errors_simulations.csv")
# Calculate min n of years for 80% power (to detect -2%/year)

min_duration_80 <- stat_df_or %>%
  filter(power >= 0.80) %>%
  group_by(sdev, phi) %>%
  summarise(min_duration_80 = min(duration, na.rm = TRUE), .groups = "drop")

mean(min_duration_80$min_duration_80)

# Formula to convert standard deviation in Coefficient of Variation (CV)
stat_df_or$sdev_original = stat_df_or$sdev
stat_df_or$cv = round(sqrt(exp(stat_df_or$sdev ^ 2) - 1), 1)
stat_df_or$cv_factor = factor(stat_df_or$cv, levels = sort(unique(stat_df_or$cv)))


# Figure 2. RELIABILITY METRICS

n_unconvergent_models=max(stat_df_or$n_models)*0.2

stat_df=stat_df_or %>% filter(n_models>n_unconvergent_models) %>% ungroup() 


# 72 combinations
colors = scales::rescale(c(log(min(stat_df$error)), log(0.02), log(max(stat_df$error))))
colnames(stat_df)
p1 = stat_df  %>%
  filter(phi %in% c(-0.75, 0.00, 0.75)) %>%
  tidyr::complete(duration, sdev_original, phi) %>%  
  ggplot(aes(x = duration, y = sdev_original, fill = error * 100)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = c("#6e9fc1", "white", "#e87061"),
    values = colors,
    trans = "log",
    labels = scales::label_number() ,
    na.value = "grey80") +
  scale_x_continuous(breaks = c(3, 12, 21, 30, 39)) +
  labs(fill = "Estimation error (%/year)", title = "A)") +
  xlab("Duration (years)") +
  ylab("Inter-annual variability")  + theme_minimal()  +
  facet_wrap( ~ phi, labeller = label_bquote(phi == .(phi))) +
  coord_cartesian(expand = FALSE) +  theme_minimal(base_size = 16)

p2 = stat_df  %>%
  filter(phi %in% c(-0.75, 0.00, 0.75), n_pow != 0) %>%
  tidyr::complete(duration, sdev_original, phi)%>%  
  ggplot(aes(x = duration, y = sdev_original, fill = power * 100)) +
  geom_raster() +
  scale_x_continuous(breaks = c(3, 21, 39)) +
  scale_fill_gradientn(colours = c("#e87061", "white", "#6e9fc1"),
                       values = c(0, .8, 1),
                       na.value = "grey80") +
  labs(fill = "Power (%)\nto detect\na 2%/year\ndecline", title = "D)") +
  xlab("Duration (years)") +
  ylab("Inter-annual variability")  + theme_minimal()  +
  facet_wrap( ~ phi, drop = TRUE, labeller = label_bquote(phi == .(phi))) +
  coord_cartesian(expand = FALSE)


p3 = stat_df  %>%
  filter(phi %in% c(-0.75, 0.00, 0.75), !n_fp %in% c(0, 1)) %>%
  tidyr::complete(duration, sdev_original, phi)%>%  
  ggplot(aes(x = duration, y = sdev_original, fill = false_pos * 100)) +
  scale_x_continuous(breaks = c(3, 21, 39)) +
  geom_raster() +
  scale_fill_gradientn(
    colours = c("#6e9fc1", "white", "#e87061"),
    values = rescale(c(0, 5, 100)),
    limits = c(0, 100) ,
    na.value = "grey80" ) +
  labs(fill = "False positives (%)", title = "E)") +
  xlab("Duration (years)") +
  ylab("Inter-annual variability")  + theme_minimal()  +
  facet_wrap( ~ phi, drop = TRUE, labeller = label_bquote(phi == .(phi))) +
  coord_cartesian(expand = FALSE)


p4 = stat_df %>%
  filter(phi %in% c(-0.75, 0.00, 0.75)) %>%
  tidyr::complete(duration, sdev_original, phi)%>%  
  ggplot(aes(x = duration, y = sdev_original, fill = bias * 100)) +
  geom_raster() +
  scale_fill_gradientn(colours = c("#8B4513", "white", "#ACA0DC"),
                       values = c(0, .4, 1),
                       na.value = "grey80") +
  scale_x_continuous(breaks = c(3, 21, 39)) +
  labs(fill = "Bias (%/year)", title = "D)") +
  xlab("Duration (years)") +
  ylab("Inter-annual variability")  + theme_minimal()  +
  facet_wrap( ~ phi, drop = TRUE, labeller = label_bquote(phi == .(phi))) +
  coord_cartesian(expand = FALSE)


p5 = stat_df   %>%
  filter(phi %in% c(-0.75, 0.00, 0.75)) %>%
  tidyr::complete(duration, sdev_original, phi)%>%  
  ggplot(aes(x = duration, y = sdev_original, fill = type_S * 100)) +
  geom_raster() +
  scale_fill_gradientn(
    colours = c("#6e9fc1", "white", "#e87061"),
    values = rescale(c(0, 1, 100)),
    limits = c(0, 100) ,
    na.value = "grey80" ) +
  scale_x_continuous(breaks = c(3, 21, 39)) +
  labs(fill = "Type S error (%)", title = "C)") +
  xlab("Duration (years)") +
  ylab("Inter-annual variability")  + theme_minimal()  +
  facet_wrap( ~ phi, drop = TRUE, labeller = label_bquote(phi == .(phi))) +
  coord_cartesian(expand = FALSE)


p6 = stat_df  %>%
  filter(phi %in% c(-0.75, 0.00, 0.75)) %>%
  tidyr::complete(duration, sdev_original, phi)%>%  
  ggplot(aes(x = duration, y = sdev_original, fill = type_M)) +
  geom_raster() +
  scale_x_continuous(breaks = c(3, 21, 39)) +
  scale_fill_gradientn(
    colours = c("#6e9fc1", "white", "#e87061"),
    values = rescale(c(1, 1.2, 2, 20)),
    limits = c(1, 20) ,
    na.value = "grey80" ) +
  labs(fill = "Type M error", title = "B)") +
  xlab("Duration (years)") +
  ylab("Inter-annual variability")  + theme_minimal()  +
  facet_wrap( ~ phi, drop = TRUE, labeller = label_bquote(phi == .(phi))) +
  coord_cartesian(expand = FALSE)

bottom <- (p6 + p5) / (p2 + p3)

Fig2 = p1 / bottom + plot_layout(heights = c(1.7, 2))

ggsave("RESULTS/Fig_2.png", Fig2, scale = .8, width=15, height = 10)

# Read BIOTIME RESULTS

# setwd(here("RESULTS/results-biotime_march_2026"))
# list_models= list.files()
# 
# table=NULL
# 
# for (i in list_models) {
#   load(i)
#   table = rbind(table, res)
# }
# 
# fwrite(table, "biotime_table.csv")

table = fread("CODE/data/biotime_table.csv")

# filter
table2 = table %>% filter(convergence == T,
                          hessian_pos_def == T,
                          !is.na(log_lik),
                          !is.na(phi_estimated_lwr))

table2 = table2 %>% mutate(phi_estimated = if_else(sdev_estimated < 1e-5, 0, phi_estimated))
table2 = table2 %>% mutate(phi_estimated = if_else(
  abs(phi_estimated) > 0.89 &
    sdev_estimated < 1e-3, 0, phi_estimated))
table2 = table2  %>% mutate(duration = n_years)

# Explore empirical data from biotime

table2 <- table2 %>%
  mutate(study_id = word(model_id, 1, sep = "-"))

table2 %>% group_by (class) %>%
  summarize(
    duration_mean = mean(duration, na.rm = T),
    sdev_mean = mean(sdev_estimated, na.rm = T),
    phi_mean = mean(phi_estimated, na.rm = T),
    models_sum = n_distinct(model_id),
    spp_sum = n_distinct(species),
    study_id = n_distinct(study_id)
  )

sum_study = table2 %>% group_by (study_id) %>%
  reframe(
    duration_mean = mean(duration, na.rm = T),
    sdev_mean = mean(sdev_estimated, na.rm = T),
    phi_mean = mean(phi_estimated, na.rm = T),
    models_sum = n_distinct(model_id),
    spp_sum = n_distinct(species),
    class = class
  ) %>% distinct()

table2 <- table2 %>%
  group_by(class) %>% 
  mutate(study_index = as.numeric(as.factor(study_id))) %>%
  ungroup()

# Explore metadata
# meta_data=fread("/Users/juliag.dealedo/ONE/Postdoc/ANTENNA/WP3/DATA/BIOTIME/biotime_v2_full_2025/data/biotime_v2_metadata_2025.csv")

# Figure 3. DATA BIOTIME

table2$sdev_estimated_original = table2$sdev_estimated
table2$cv_estimated = round(sqrt(exp(table2$sdev_estimated ^ 2) - 1), 1)
table2$cv_estimated_factor = factor(table2$cv_estimated, levels = sort(unique(table2$cv_estimated)))

table2  %>% group_by(class) %>% summarise(mean = mean(sdev_estimated))

p2.1 = table2 %>% filter(!is.na(duration)) %>% ggplot(aes(y = n_years, x = class, fill = class)) +
  geom_jitter(
    alpha = 0.1,
    color = "grey70",
    height = 0,
    size = 2
  ) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = cols) + ylab("Duration (years)") + xlab("") +
  guides(fill = "none") +
  theme_minimal(base_size = 22) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))


p2.2 = table2  %>% ggplot(aes(y = sdev_estimated, x = class, fill = class)) +
  geom_jitter(
    alpha = 0.1,
    color = "grey70",
    height = 0,
    size = 2
  ) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = cols) + ylab("Inter-annual variability") + xlab("") +
  guides(fill = "none") +
  theme_minimal(base_size = 22) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p2.3 = table2 %>% ggplot(aes(y = phi_estimated, x = class, fill = class)) +
  geom_jitter(
    alpha = 0.1,
    color = "grey70",
    height = 0,
    size = 2
  ) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = cols) + ylab("Temporal autocorrelation") + xlab("") +
  guides(fill = "none") +
  theme_minimal(base_size = 22) +  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p2.1 + p2.2 + p2.3 + plot_annotation(tag_levels = c("A", ")"))

ggsave("RESULTS/Fig_3.png", scale = 1, width=15, height = 10)


## Figure 4. Reliability + biotime
# Extract values by rounding

categories_sdev <- unique(stat_df$sdev)
categories_dur <- unique(stat_df$duration)
categories_phi <- unique(stat_df$phi)

table3 <- table2 %>%
  mutate(sdev_rounded = categories_sdev [apply(outer(sdev_estimated, categories_sdev, function(x, y)
    abs(x - y)), 1, which.min)]) %>%
  mutate(duration_rounded = categories_dur [apply(outer(duration, categories_dur, function(x, y)
    abs(x - y)), 1, which.min)]) %>%
  mutate(phi_rounded = categories_phi [apply(outer(phi_estimated, categories_phi, function(x, y)
    abs(x - y)), 1, which.min)])

table3$combi = paste(table3$duration_rounded,
                     table3$phi_rounded,
                     table3$sdev_rounded)

# BIAS BOXPLOT

# Load raw data

stat_df_biotime = resf_filtered %>%
  group_by(duration, sdev, phi) %>%
  filter(between(r_estimated, -2, 2)) %>%
  summarise(
    bias = mean((r_estimated - r), na.rm = TRUE),
    error = mean(sqrt(((r_estimated - r)) ^ 2), na.rm = TRUE),
    bias_phi = mean((phi_estimated - phi), na.rm = TRUE),
    bias_sdev = mean((sdev_estimated - sdev), na.rm = TRUE),
    bias_phi2 = mean(abs(phi_estimated - phi)),
    bias_sdev2 = mean(abs(sdev_estimated - sdev)),
    n_pow = sum(r == 0, na.rm = TRUE),
    power = sum(slope_pval[r == (-0.02) & r_estimated < 0] < 0.05, na.rm = TRUE) / sum(r == (-0.02), na.rm = TRUE),
    false_pos = sum(slope_pval[r == 0] < 0.05, na.rm = TRUE) / sum(r == 0, na.rm = TRUE),
    type_M = mean(abs(r_estimated[slope_pval < 0.05 & abs(r) > 1e-6] / r[slope_pval < 0.05 & abs(r) > 1e-6]), na.rm = TRUE),
    type_S = sum(slope_pval < 0.05 & r != 0 & sign(r) != sign(r_estimated), na.rm = TRUE) / sum(slope_pval < 0.05 & r != 0, na.rm = TRUE),
    n_models = length(model_id)
  )


stat_df_biotime$combi = paste(stat_df_biotime$duration,
                              stat_df_biotime$phi,
                              stat_df_biotime$sdev)
error_biotime = merge(stat_df_biotime, table3, by = "combi")

error_biotime_filtered = error_biotime %>% filter(!is.na(error))
#fwrite(error_biotime_filtered, "error_biotime_filtered.csv")


error_biotime_filtered = fread("CODE/data/error_biotime_filtered.csv")

nrow(error_biotime_filtered %>% filter(n_models<=1020))/nrow(error_biotime_filtered)*100
16000*0.03


error_biotime_filtered=error_biotime_filtered %>% filter(n_models>1020)


# GRAPHS
perc_error <- error_biotime_filtered %>%
  group_by(class) %>%
  summarise(perc_reliable_error = 100 * mean(error <= 0.02, na.rm = TRUE),
            .groups = "drop")

p_error <- error_biotime_filtered %>%
  ggplot(aes(x = class, y = error * 100, fill = class)) +
  geom_jitter(
    alpha = 0.2,
    size = 2,
    color = "grey70",
    height = 0.00
  ) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_text(
    data = perc_error,
    aes(x = class,
        y = .7,
        label = paste0(round(perc_reliable_error, 1), "%")),
    vjust = 1.5,
    inherit.aes = FALSE
  ) +
  ylim(-10, 100) +
  scale_fill_manual(values = cols) +
  scale_y_log10() +
  theme_minimal(base_size = 22) +
  labs(title = "Expected error in trend estimation", #subtitle = "Text shows % of time series with expected error ≤ 2%/year",
       x = "Class", y = "Expected error in growth rate (%/year)") +  guides(fill = "none")


perc_power <- error_biotime_filtered %>%
  group_by(class) %>%
  summarise(perc_reliable_power = 100 * mean(power >= 0.80, na.rm = TRUE),
            .groups = "drop")

p_power <- error_biotime_filtered %>%
  ggplot(aes(x = class, y = power * 100, fill = class)) +
  geom_jitter(
    alpha = 0.2,
    size = 2,
    color = "grey70",
    height = 0.00
  ) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_hline(yintercept = 80, linetype = "dashed") +
  geom_text(
    data = perc_power,
    aes(
      x = class,
      y = 100,
      label = paste0(round(perc_reliable_power, 1), "%")
    ),
    vjust = 1.5,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = cols) +
  theme_minimal(base_size = 22) +
  labs(title = "Expected power to detect -2%/year", #subtitle = "Text shows % of time series with power ≥ 80% to detect a −2%/year trend",
       x = "Class", y = "Power (%)") +  guides(fill = "none") + theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_power
perc_bias <- error_biotime_filtered %>%
  group_by(class) %>%
  summarise(perc_reliable_bias = 100 * mean(abs(bias) <= 0.02, na.rm = TRUE),
            .groups = "drop")

p_bias <- error_biotime_filtered %>%
  filter(bias > -0.15) %>%
  ggplot(aes(x = class, y = bias * 100, fill = class)) +
  geom_jitter(
    alpha = 0.2,
    size = 2,
    color = "grey70",
    height = 0.00
  ) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -2, linetype = "dashed") +
  geom_text(
    data = perc_bias,
    aes(
      x = class,
      y = 4,
      label = paste0(round(perc_reliable_bias, 1), "%")
    ),
    vjust = 1.5,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = cols) +
  theme_minimal(base_size = 22) +
  labs(
    title = "Expected bias in trend estimation",
    subtitle = "Text shows % of time series with expected bias between −2% and +2%/year",
    x = "Class",
    y = "Expected bias in growth rate (%/year)"
  ) +  guides(fill = "none") + theme(axis.text.x = element_text(angle = 20, hjust = 1))


perc_falsepositives <- error_biotime_filtered %>%
  group_by(class) %>%
  summarise(
    perc_reliable_fp = 100 * mean(false_pos <= 0.05, na.rm = TRUE),
    .groups = "drop"
  )

p_false <- error_biotime_filtered %>%
  ggplot(aes(x = class, y = false_pos * 100, fill = class)) +
  geom_jitter(
    alpha = 0.2,
    size = 2,
    color = "grey70",
    height = 0.00
  ) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_hline(yintercept = 5, linetype = "dashed") +
  geom_text(
    data = perc_falsepositives,
    aes(
      x = class,
      y = 1.6,
      label = paste0(round(perc_reliable_fp, 1), "%")
    ),
    vjust = 1.5,
    inherit.aes = FALSE
  ) +
  ylim(-10, 100) +
  scale_fill_manual(values = cols) +
  theme_minimal(base_size = 22) +
  labs(title = "Expected false positive rate", #subtitle = "Text shows % of time series with expected false positive rate ≤ 5%",
       x = "Class", y = "Expected false positive rate (%)") +  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))


error_significant <- error_biotime_filtered #%>% filter(slope_pval < 0.05)

perc_typeM <- error_significant %>%
  group_by(class) %>%
  summarise(perc_reliable_typeM = 100 * mean(type_M <= 1.2, na.rm = TRUE),
            .groups = "drop")

p_type_m <- error_significant %>%
  ggplot(aes(x = class, y = type_M, fill = class)) +
  geom_jitter(
    alpha = 0.2,
    size = 2,
    color = "grey70",
    height = 0.00
  ) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_hline(yintercept = 1.2, linetype = "dashed") +
  geom_text(
    data = perc_typeM,
    aes(
      x = class,
      y = 1.1,
      label = paste0(round(perc_reliable_typeM, 1), "%")
    ),
    vjust = 1.5,
    inherit.aes = FALSE
  ) +
  scale_y_log10() +
  ylim(-.4, 14) +
  scale_fill_manual(values = cols) +
  theme_minimal(base_size = 22) +
  labs(title = "Expected Type M error", #subtitle = "Text shows % of significant time series with expected Type M ≤ 1.2",
       x = "Class", y = "Expected Type M error") + guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))


perc_typeS <- error_significant %>%
  group_by(class) %>%
  summarise(
    perc_reliable_typeS = 100 * mean(type_S <= 0.01, na.rm = TRUE),
    .groups = "drop"  )

p_type_s <- error_significant %>%
  ggplot(aes(x = class, y = type_S, fill = class)) +
  geom_jitter(
    alpha = 0.2,
    size = 2,
    color = "grey70",
    height = 0.00
  ) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_hline(yintercept = 0.01, linetype = "dashed") +
  geom_text(
    data = perc_typeS,aes(x = class,y = 0.0015,
                          label = paste0(round(perc_reliable_typeS, 1), "%")    ),
    vjust = 1.5,
    inherit.aes = FALSE ) +
  ylim(-.0003, .15) +
  scale_y_log10() +
  scale_fill_manual(values = cols) +
  theme_minimal(base_size = 22) +
  labs(title = "Expected Type S error", # subtitle = "Text shows % of significant time series with expected Type S ≤ 1%",
       x = "Class", y = "Expected Type S error") +
  guides(fill = "none")  + theme(axis.text.x = element_text(angle = 20, hjust = 1))


bottom <- (p_type_m + p_type_s) / (p_power + p_false)

p_error / bottom +
  plot_layout(heights = c(1.4, 2.2))  + plot_annotation(tag_levels = c("A", ")"))

 ggsave("RESULTS/Fig_4.png", scale = 1.1, width=15, height=15)



############################
## Supplementary material
############################

stat_df =fread("CODE/data/errors_simulations.csv")



p_models = stat_df %>%
  filter(phi %in% c(-0.75, 0.00, 0.75)) %>%
  ggplot(aes(x = duration, y = sdev, fill = n_models)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(c( "#127475", "#0E9594", "#F5DFBB", "orange", "#F2542D")))+
  #   limits = c(300,5100),
  #   breaks = c(300, 1000, 3000, 5000)) +
  scale_x_continuous(breaks = c(3, 9, 15, 21, 27, 33, 39)) +
  labs(fill = "Number of models per cell") +
  xlab("Duration (years)") +
  ylab("Inter-annual variability") +
  theme_minimal() +
  facet_wrap( ~ phi, drop = TRUE, labeller = label_bquote(phi == .(phi))) + coord_cartesian(expand =
                                                                                              FALSE)

p_models

ggsave("RESULTS/Fig_S1.png", p_models, width=7, height=2)

#

## Bias

p_bias_phi = ggplot(data = stat_df, aes(x = duration, y = sdev, fill = bias_phi_sim)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "#e87061",
    mid = "grey95",
    high = "#6e9fc1",
    midpoint = 0.0
  ) +
  labs(fill = "Bias in autocorrelation estimation") +
  xlab("Duration (years)") +
  ylab("Inter-annual variability")  +
  facet_wrap( ~ phi, drop = TRUE, labeller = label_bquote(phi == .(phi))) +
  coord_cartesian(expand = FALSE)+ theme_bw() + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  ) 
p_bias_phi

p_bias_sdev = ggplot(data = stat_df, aes(x = duration, y = sdev, fill = bias_sdev_sim)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "#e87061",
    mid = "grey95",
    high = "#6e9fc1",
    midpoint = 0.0
  )  +
  labs(fill = "Bias in inter-annual variability estimation") +
  xlab("Duration (years)") +
  ylab("Inter-annual variability")  +
  facet_wrap( ~ phi, drop = TRUE, labeller = label_bquote(phi == .(phi))) +
  coord_cartesian(expand = FALSE)+ theme_bw() + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  ) 

p_bias_slope = ggplot(data = stat_df, aes(x = duration, y = sdev, fill = bias)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "#e87061",
    mid = "grey95",
    high = "#6e9fc1",
    midpoint = 0.0
  ) +
  labs(fill = "Bias in slope estimation") +
  xlab("Duration (years)") +
  ylab("Inter-annual variability")  +
  facet_wrap( ~ phi, drop = TRUE, labeller = label_bquote(phi == .(phi))) +
  coord_cartesian(expand = FALSE)+ theme_bw() + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  ) 

FigS1 = p_bias_phi / p_bias_sdev / p_bias_slope + plot_annotation(tag_levels = c("A", ")"))

ggsave("RESULTS/Fig_S2.png", FigS1, width = 8, height = 13)


## Correlation among metrics and characteristics


M <- stat_df_or %>%
  ungroup() %>%  select(power, type_M,  type_S,  error,  false_pos,  duration,  sdev,  phi)

# Better names for the figure
names(M) <- c(
  "Power",
  "Type M",
  "Type S",
  "Estimation error",
  "False positives",
  "Duration",
  "Inter-annual variability",
  "Autocorrelation"
)

# Full Spearman correlation matrix
M_cor <- cor( M,
              use = "pairwise.complete.obs",
              method = "pearson")

png("RESULTS/Fig_S3.3.png", width = 2400, height = 2100, res = 300)
corrplot(
  M_cor,
  method = "color",
  type = "upper",
  order = "original",
  addCoef.col = "black",
  col = colorRampPalette(
    c("#e87061", "white", "lightblue")
  )(200),
  number.cex = 0.75,
  number.digits = 2,
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 0.9,  diag = FALSE)
dev.off()

# Figure MAP


table = fread("CODE/data/biotime_table.csv")

# filter
table2 = table %>% filter(convergence == T,
                          hessian_pos_def == T,!is.na(log_lik),!is.na(phi_estimated_lwr))

table2 = table2 %>% mutate(phi_estimated = if_else(sdev_estimated < 1e-5, 0, phi_estimated))
table2 = table2 %>% mutate(phi_estimated = if_else(
  abs(phi_estimated) > 0.89 &
    sdev_estimated < 1e-3,
  0,
  phi_estimated
))
table2 = table2  %>% mutate(duration = n_years)

# Map: Figure S3
conversor = fread("CODE/data/conversor.csv")

colnames(conversor) = c("LATITUDE", "LONGITUDE", "model_id")
table_map = table2 %>% left_join(conversor, by = "model_id")
table_map = table_map %>% group_by(model_id, LATITUDE, LONGITUDE, class) %>%
  summarise(n_models = n_distinct(model_id))

world <- map_data("world")

map=ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "white",
    color = "gray80"
  ) +
  geom_point(
    data = table_map,
    aes(x = LONGITUDE, y = LATITUDE, color = class),
    size = 1,
    alpha = .2
  ) +
  scale_color_manual(values = cols) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  coord_fixed(1.3) +
  theme_void()

ggsave("RESULTS/Fig_S3.png", map, width=15, height=10)


# Summary of the models
table2 <- table2 %>%
  mutate(study_id = word(model_id, 1, sep = "-"))

table2 %>% group_by (class) %>% 
  summarize(duration_mean = mean(duration, na.rm = T),
            sdev_mean = mean(sdev_estimated, na.rm = T),
            phi_mean = mean(phi_estimated, na.rm=T),
            models_sum = n_distinct(model_id),
            spp_sum = n_distinct(species),
            study_id = n_distinct(study_id))

table2 <- table2 %>%
  group_by(class) %>% #filter(class=="Amphibia") %>%
  mutate(study_index = as.numeric(as.factor(study_id))) %>%
  ungroup()
n_studies = table2 %>% group_by(class) %>% summarise (n_studies = n_distinct(study_id))
max_studies <- max(table2$study_index)
cols_local <- qualitative_hcl(max_studies, palette = "Dynamic")

p_sup_mat <- table2 %>%
  ggplot(aes(x = class, y = n_years)) +
  geom_jitter(
    aes(color = factor(study_index)),
    alpha = .7,
    size = 1.2,
    width = .3,
    height = 0.00
  ) +
  geom_text(
    data = n_studies,
    aes(x = class, y = 0, label = n_studies),
    # <-- important
    inherit.aes = FALSE,
    vjust = 1.5
  ) +
  geom_boxplot(aes(fill = class), alpha = 0.2, outlier.shape = NA) +
  scale_color_manual(values = cols_local, guide = "none") +
  scale_fill_manual(values = cols) +
  theme_minimal() +
  theme(legend.position = "none") +
  ylab("Duration (years)") +
  xlab("Class (labels and colours indicate number of studies)")


ggsave("RESULTS/Fig_S4.png", p_sup_mat, scale = 1, width=7, height=4.5)


# Power and time series characteristics

load("CODE/data/model-simulation-join.RData")

resf$conv=1
resf$conv[which(resf$convergence==FALSE | is.na(resf$aic) | is.na(resf$mu_zero_estimated_se) | resf$hessian_pos_def==FALSE)]=0
resf_filtered=subset(resf, conv==1)
rm(resf)

length(unique(resf_filtered$phi))*length(unique(resf_filtered$sdev))*length(unique(resf_filtered$r))*length(unique(resf_filtered$duration))*100

resf_filtered$r=round(resf_filtered$r, 2)
resf_filtered$phi=round(resf_filtered$phi, 2)

power_df <- resf_filtered %>%
  filter(between(r_estimated, -2, 2), duration < 40, sdev < 2, r != 0) %>%
  group_by(duration, sdev, phi, r) %>%
  summarise(
    power = mean(slope_pval < 0.05, na.rm = TRUE),
    n_models = n(),
    .groups = "drop"
  )



# Plot power


duration_df <- power_df %>%
  filter(n_models > 20) %>%
  group_by(duration) %>%
  summarise(
    mean_power = mean(power, na.rm = TRUE),
    n = sum(!is.na(power)),
    sd_power = sd(power, na.rm = TRUE),
    se = sd_power / sqrt(n),
    lower_95 = mean_power -  1.96 * se,
    upper_95 = mean_power + 1.96 * se,
    .groups = "drop"
  )

duration <- duration_df %>%
  ggplot(aes(x = duration, y = mean_power)) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 2) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  theme_minimal() +
  labs(x = "Duration (years)", y = "Mean power")

auto = power_df %>%
  filter(n_models > 20) %>%
  group_by(phi) %>%
  summarise(
    mean_power = mean(power, na.rm = TRUE),
    n = sum(!is.na(power)),
    sd_power = sd(power, na.rm = TRUE),
    se = sd_power / sqrt(n),
    lower_95 = mean_power - 1.96 * se,
    upper_95 = mean_power + 1.96 * se,
    .groups = "drop"
  ) %>%
  ggplot(aes(x = phi, y = mean_power)) + geom_line() + ylim(0, 1) +
  geom_point() +
  ylim(0, 1) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = .11) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  theme_minimal() + labs(x = "Autocorrelation", y = "Mean power")

variab = power_df %>%
  filter(n_models > 20) %>%
  group_by(sdev) %>%
  summarise(
    mean_power = mean(power, na.rm = TRUE),
    n = sum(!is.na(power)),
    sd_power = sd(power, na.rm = TRUE),
    se = sd_power / sqrt(n),
    lower_95 = mean_power - 1.96 * se,
    upper_95 = mean_power + 1.96 * se,
    .groups = "drop"
  ) %>%
  ggplot(aes(x = sdev, y = mean_power)) +
  geom_line() +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = .11) +
  ylim(0, 1) +
  geom_point() +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  theme_minimal() +
  labs(x = "Inter-annual variability ", y = "Mean power")

slope = power_df %>%
  filter(n_models > 20) %>%
  group_by(r) %>%
  summarise(
    mean_power = mean(power, na.rm = TRUE),
    n = sum(!is.na(power)),
    sd_power = sd(power, na.rm = TRUE),
    se = sd_power / sqrt(n),
    lower_95 = mean_power - 1.96 * se,
    upper_95 = mean_power + 1.96 * se,
    .groups = "drop"
  ) %>%
  ggplot(aes(x = r, y = mean_power)) +
  geom_line() +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = .011) +
  geom_point() +
  ylim(0, 1) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  theme_minimal() +
  labs(x = "Temporal trend", y = "Mean power")


library(patchwork)
duration + auto + variab + slope + plot_annotation(tag_levels = "A")
ggsave("RESULTS/Fig_S5.png",  scale = 1, width=7, height=4.5)



# Fig S6

biotime_error_duration <- error_biotime_filtered %>%
  group_by(duration.x, class) %>%
  summarise(
    prop_reliable = mean(error <= 0.020, na.rm = TRUE),
    n = n(),
    .groups = "drop")

figS9 = ggplot(biotime_error_duration,
               aes(x = duration.x, y = prop_reliable * 100, color = class)) +
  geom_line() +
  geom_point(aes(size = (n)), alpha = 0.8) +
  geom_hline(yintercept = 50, linetype = "dashed") +
  theme_minimal() +
  scale_color_manual(values = cols) +
  labs(x = "Duration (years)",
       y = "% time series with expected error ≤ 2%/year",
       color = "Class",
       size = "N time series")

ggsave("RESULTS/Fig_S6.png", figS9,  height = 4.5, width=7)


# Partial pooling

library(glmmTMB)
library(dplyr)
library(data.table)
library(tidyverse)

gm_mean <- function(x){exp(mean(log(x)))}

setwd(here("CODE/data/data_simulation_march_2026/"))
list_models= list.files()

table_sim=NULL

for (i in list_models) {
  load(i)
  table_sim = rbind(table_sim, resf)
}


d=table_sim
rm(resf, table_sim)

#we will draw three lines, fixing the rest of the parameters. 
pooling <- c(1,5,10,15,20,25,30,35,40,45,50) 
d2 = d %>% filter (phi==0, year<11, sdev %in% c(0.15, 1.5))
d2$sim_id <- paste(d2$repi, d2$r, d2$sdev, d2$phi, sep = "_")

#loop
once <- function(x, pooling, sdev){
  x <- subset(x, sdev == sdev)
  out <- data.frame(pool_n = NA, error = NA, variance = NA)
  for(i in 1:length(pooling)){
    #select the right simulations
    subsampled_sim <- sample(x$sim_id, pooling[i], replace = F)
    temp <- x[which(x$sim_id %in% subsampled_sim),]
    #calculate real trend
    real_trend <- gm_mean(temp$r + 1) - 1 #WARNING: Need to use better metric? Geometric mean might be ok.
    
    #calculate observed trend
    temp$year2 <- temp$year - min(temp$year) #to start at year 0
    temp$YEAR3 <- as.factor(temp$year2) #autocor strusture
    if(is.numeric(temp$abund)){
      m <- glmmTMB(abund ~ year2 + ar1(YEAR3 + 0 | sim_id),
                   data = temp,
                   family = poisson)
      obs_estimate <- coef(m)$cond$sim_id$year2
    }else{
      print("ABUNDANCE NOT NUMERIC, STOP THE MACHINES")
    }
    #calculate error
    obs_trend = exp(obs_estimate) - 1
    error <- mean(sqrt((real_trend-obs_trend)^2))
    #store values
    out[i,"pool_n"] <- pooling[i]
    out[i,"error"] <- error
    out[i,"variance"] <- sdev
  }
  out
}

out <- once(x = d2, pooling = pooling, sdev = 0.15)
scatter.smooth(out$error ~ out$pool_n, las = 1)

#Now lets do 100 iterations (a bit slow)

for(j in 1:100){
  temp <- once(x = d2, pooling = pooling, sdev = 0.15)
  out <- rbind(out, temp)
}
#fix warning message when model don't converge? Ignoring it for now.

scatter.smooth(out$error ~ out$pool_n, las = 1)
boxplot(out$error ~ out$pool_n, las = 1)
#cleaner plotting only means:
mean_out07 <- tapply(out$error, out$pool_n, mean)
plot(mean_out07, las = 1, type = "b", col = "red")
unique(d2$sdev)
#Now we can do for var = 0.1 and 1.5
out_01 <- once(x = d2, pooling = pooling, sdev = 0.15)
for(j in 1:100){
  temp <- once(x = d2, pooling = pooling, sdev = 0.15)
  out_01 <- rbind(out_01, temp)
}
mean_out01 <- tapply(out_01$error, out_01$pool_n, mean)
sd_out01 <- tapply(out_01$error, out_01$pool_n, sd)

out_15 <- once(x = d2, pooling = pooling, sdev = 1.5)
for(j in 1:100){
  temp <- once(x = d2, pooling = pooling, sdev = 1.5)
  out_15 <- rbind(out_15, temp)
}

out_total = rbind(out_15, out_01)

fwrite(out_total, "RESULTS/pp.csv")

mean_out15 <- tapply(out_15$error, out_15$pool_n, mean)
sd_out15 <- tapply(out_15$error, out_15$pool_n, sd)
n_out01 <- tapply(
  out_01$error,
  out_01$pool_n,
  function(x) sum(!is.na(x))
)

n_out15 <- tapply(
  out_15$error,
  out_15$pool_n,
  function(x) sum(!is.na(x))
)

round(sqrt(exp(0.15^2)-1),1)


df_plot <- data.frame(
  pooling = as.numeric(names(mean_out01)),
  var_015 = mean_out01,
  var_15 = mean_out15,
  n_15 = n_out15,
  n_01 = n_out01,
  sd_01 = sd_out01,
  sd_15 = sd_out15
) %>%
  pivot_longer(cols = starts_with("var_"),
               names_to = "variance",
               values_to = "mean_error")  %>%
  mutate(
    sd = ifelse(variance == "var_015", sd_01, sd_15),
    n  = ifelse(variance == "var_015", n_01, n_15),
    se = sd / sqrt(n),
    lower_95 = mean_error - 1.96 * se,
    upper_95 = mean_error + 1.96 * se
  )




ggplot(df_plot,
       aes(x = pooling,
           y = mean_error * 100,
           color = variance,
           group = variance)) +
  geom_errorbar(aes(ymin = lower_95 * 100, ymax = upper_95 * 100), width = 0.4) +
  geom_point() +
  geom_line() +
  scale_color_manual(
    values = c(
      "var_015" = "darkred",
      "var_15"  = "goldenrod2"
    ),
    labels = c("var_015" = "0.15", "var_15"  = "1.5")
  ) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  labs(x = "Partial pooling", y = "Mean error (%/year)", color = "Inter-annual variability") +
  theme_minimal()



ggsave("RESULTS/Fig_S7.png", scale=0.5)



# Read filtered table. Code to TableS1

error_biotime_filtered = fread("CODE/data/error_biotime_filtered.csv")


error_biotime_filtered <- error_biotime_filtered %>%
  mutate(
    safe_error = error <= 0.02,
    safe_power = power >= 0.80,
    safe_fp = false_pos <= 0.05,
    safe_typeS = type_S <= 0.01,
    safe_typeM = type_M < 1.1
  )

safe_summary_wide <- error_biotime_filtered %>%
  select(class,
         slope_estimated,
         safe_typeM,
         safe_typeS,
         safe_fp,
         safe_power,
         safe_error) %>%
  pivot_longer(cols = starts_with("safe_"),
               names_to = "metric",
               values_to = "safe") %>%
  mutate(
    metric = gsub("safe_", "", metric),
    safe = ifelse(safe, "reliable", "unreliable")
  ) %>%
  group_by(class, metric, safe) %>%
  summarise(
    n = n(),
    mean = mean(slope_estimated, na.rm = TRUE),
    median = median(slope_estimated, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = safe,
    values_from = c(n, mean, median),
    names_glue = "{.value}_{safe}")

all_summary <- error_biotime_filtered %>%
  select(class,
         slope_estimated,
         safe_typeM,
         safe_typeS,
         safe_fp,
         safe_power,
         safe_error) %>%
  pivot_longer(cols = starts_with("safe_"),
               names_to = "metric",
               values_to = "safe") %>%
  mutate(metric = gsub("safe_", "", metric)) %>%
  group_by(class, metric) %>%
  summarise(
    n_all = n(),
    mean_all = mean(slope_estimated, na.rm = TRUE),
    median_all = median(slope_estimated, na.rm = TRUE),
    .groups = "drop"
  )

safe_summary_wide <- safe_summary_wide %>%
  left_join (all_summary, by = c("class", "metric")) %>%
  select (class, metric, n_reliable, n_unreliable, n_all, mean_reliable, mean_unreliable, 
          mean_all, median_reliable, median_unreliable, median_all) %>%
  mutate (across(c (mean_reliable, mean_unreliable, mean_all, median_reliable, 
                    median_unreliable, median_all)))

# To print in a table
fwrite(safe_summary_wide, "RESULTS/tableS1.csv")
# print(qflextable(safe_summary_wide), preview = "docx")




