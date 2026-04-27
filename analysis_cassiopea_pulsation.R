########################################################### #
# Dawn peaks and individual variation shape functional pulsation rhythms in the 
#photosymbiotic jellyfish Cassiopea 
# Authors: Noémie Buratto and Olivia Bleeckx

# This script reproduces the analyses and figures reported in the manuscript.
# Additional / exploratory analyses not reported in the main text are provided
# at the end of this file under "ADDITIONAL CODE".
########################################################### #

########################### #
# 0) Packages + settings ####
########################### #
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)       
library(rptR)
library(janitor)
library(lubridate)
library(lomb)      
library(mgcv) 
library(MASS)

theme_set(theme_minimal(base_size = 13))
options(stringsAsFactors = FALSE)

geommean <- function(x) exp(mean(log(x)))

pal_tp  <- c("Dawn"="#8B5E3C","Noon"="#556B2F","Dusk"="#3B6FB6","Midnight"="#7E5AA7")
pal_day <- c("D1"="#6B8E23","D2"="#5F9EA0","D3"="#8B4513")

########################### #
# 1) Import + cleaning ####
########################### #

path <- "cassiopea_pulsation_dataset.xlsx"

df <- read_excel(path) |>
  janitor::clean_names()

required <- c(
  "individual_id","day","time_point","time_med",
  "pulse_count_rep1","pulse_count_rep2","pulsation_frequency",
  "lux_mean","ph","temperature_c","oxygen_mg_l",
  "conductivity","salinity","umbrella_diameter_cm"
)
missing <- setdiff(required, names(df))
if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

df <- df |>
  mutate(
    individual_id = factor(individual_id),
    
    # standardize day as D1/D2/D3
    day = as.character(day),
    day = if_else(grepl("^D", day), day, paste0("D", day)),
    day = factor(day),
    
    time_point = factor(time_point, levels = c("Dawn","Noon","Dusk","Midnight")),
    
    pulsation_frequency  = as.numeric(pulsation_frequency),
    umbrella_diameter_cm = as.numeric(umbrella_diameter_cm),
    temperature_c        = as.numeric(temperature_c),
    lux_mean             = as.numeric(lux_mean),
    conductivity         = as.numeric(conductivity),
    oxygen_mg_l           = as.numeric(oxygen_mg_l),
    ph                   = as.numeric(ph),
    
    z_umbrella = as.numeric(scale(umbrella_diameter_cm)),
    z_temp     = as.numeric(scale(temperature_c)),
    z_lux      = as.numeric(scale(lux_mean)),
    z_cond     = as.numeric(scale(conductivity)),
    z_ox     = as.numeric(scale(oxygen_mg_l)),
    z_ph     = as.numeric(scale(ph))
  )

# log-model requires > 0
if (any(df$pulsation_frequency <= 0, na.rm = TRUE)) {
  stop("pulsation_frequency contains non-positive values; log-transform is not defined.")
}

########################### #
# 2) Time variables for chronobiology sections ####
########################### #

# Parse time_med robustly:
parse_time_to_hours <- function(x) {
  if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) {
    return(hour(x) + minute(x)/60 + second(x)/3600)
  }
  if (inherits(x, "difftime")) {
    # unlikely; treat as hours
    return(as.numeric(x, units = "hours"))
  }
  # try character parsing
  xx <- as.character(x)
  suppressWarnings({
    tt <- hms::as_hms(xx)
    if (all(is.na(tt))) tt <- lubridate::hms(xx)
    if (all(is.na(tt))) return(rep(NA_real_, length(xx)))
    return(lubridate::hour(tt) + lubridate::minute(tt)/60 + lubridate::second(tt)/3600)
  })
}

df <- df |>
  mutate(
    time_med_h = parse_time_to_hours(time_med)
  )

# fallback if time_med missing: typical hours by time_point
tp_codes <- c(Dawn = 6, Noon = 12, Dusk = 18, Midnight = 0)
df <- df |>
  mutate(
    phase24 = if_else(is.na(time_med_h),
                      unname(tp_codes[as.character(time_point)]),
                      time_med_h),
    # keep in [0,24)
    phase24 = phase24 %% 24
  )

# Continuous time across 72 h (needed for 72 h figure + Lomb detrending)
df <- df |>
  mutate(
    day_num = as.integer(gsub("^D","", as.character(day))) - 1L,
    t72 = day_num*24 + phase24
  )


########################### #
# 3) Figure 2: Pulsation frequency over 72 h (trajectories + mean ± SE) ####
########################### #

df_plot <- df |>
  arrange(day, time_point) |>
  mutate(tp72 = factor(paste(day, time_point, sep = "-"),
                       levels = unique(paste(day, time_point, sep = "-"))))

df_sum72 <- df_plot |>
  group_by(tp72) |>
  summarise(
    mean_freq = mean(pulsation_frequency, na.rm = TRUE),
    se_freq   = sd(pulsation_frequency, na.rm = TRUE)/sqrt(sum(!is.na(pulsation_frequency))),
    .groups = "drop"
  )

p_fig72 <- ggplot() +
  geom_line(data = df_plot,
            aes(x = tp72, y = pulsation_frequency, group = individual_id),
            color = "grey75", alpha = 0.5, linewidth = 0.5) +
  geom_point(data = df_plot,
             aes(x = tp72, y = pulsation_frequency),
             color = "grey60", size = 0.9, alpha = 0.35) +
  geom_line(data = df_sum72,
            aes(x = tp72, y = mean_freq, group = 1),
            color = "#8F9779", linewidth = 1.6) +
  geom_point(data = df_sum72,
             aes(x = tp72, y = mean_freq),
             color = "#556B2F", size = 2.4) +
  geom_errorbar(data = df_sum72,
                aes(x = tp72, ymin = mean_freq - se_freq, ymax = mean_freq + se_freq),
                width = 0.12, color = "#8F9779", linewidth = 0.6) +
  labs(
    title    = expression("Pulsation frequency of "*italic("Cassiopea")*" spp. over 72 h"),
    subtitle = "Green line = mean ± SE · Grey lines = individual trajectories",
    x = "Time point (Day - Period)",
    y = "Pulsation frequency (pulse/min)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_fig72)

########################### #
# 4) Collinearity checks (temp, lux, conductivity) ####
########################### #
cor_matrix <- cor(dplyr::select(df, z_temp, z_lux, z_cond, z_ox, z_ph, z_umbrella),
                  use = "pairwise.complete.obs")
print(round(cor_matrix, 3))


anova(lm(z_temp~time_point, df))
summary(lm(z_temp ~ time_point, df))$r.squared
plot(df$time_point, df$z_cond)

anova(lm(z_lux~time_point, df))
summary(lm(z_lux ~ time_point, df))$r.squared
plot(df$time_point, df$z_lux)

anova(lm(z_cond~time_point, df))
plot(df$time_point, df$z_cond)

anova(lm(oxygen_mg_l~time_point, df))
plot(df$time_point, df$oxygen_mg_l)

anova(lm(z_ph~time_point, df))
plot(df$time_point, df$z_ph)

########################### #
# 5) Main model: LMM time_point + umbrella + random effects (individual, day) ####
########################### #

m_cycle_ML <- lmer(log(pulsation_frequency) ~ time_point +
                     (1|individual_id) + (1|day),
                   data = df, REML = FALSE)

m_cycle_umb_ML <- lmer(log(pulsation_frequency) ~ time_point + z_umbrella +
                     (1|individual_id) + (1|day),
                   data = df, REML = FALSE)

m_cycle_cond_ML <- lmer(log(pulsation_frequency) ~ time_point + z_umbrella + z_cond +
                          (1|individual_id) + (1|day),
                        data = df, REML = FALSE)

m_cycle_lux_ML <- lmer(log(pulsation_frequency) ~ time_point + z_umbrella + z_lux +
                          (1|individual_id) + (1|day),
                        data = df, REML = FALSE)

m_cycle_temp_ML <- lmer(log(pulsation_frequency) ~ time_point + z_umbrella + z_temp +
                         (1|individual_id) + (1|day),
                       data = df, REML = FALSE)

m_cycle_umb_interaction_ML <- lmer(log(pulsation_frequency) ~ time_point * z_umbrella +
                         (1|individual_id) + (1|day),
                       data = df, REML = FALSE)

anova(m_cycle_ML, m_cycle_umb_ML)
anova(m_cycle_umb_ML, m_cycle_cond_ML)
anova(m_cycle_umb_ML, m_cycle_lux_ML)
anova(m_cycle_umb_ML, m_cycle_temp_ML)
anova(m_cycle_umb_ML, m_cycle_umb_interaction_ML)


# final model (as in paper: conductivity not retained)
m_final <- update(m_cycle_umb_ML, REML = TRUE)
print(summary(m_final))
anova(m_final)

emm_tp <- emmeans(m_final, ~ time_point)
print(pairs(emm_tp, adjust = "tukey", type = "response"))

########################### #
# 6) Diagnostics (LMM) ####
########################### #
car::qqPlot(resid(m_final))
plot(fitted(m_final), resid(m_final),
     xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs fitted")
abline(h = 0, col = "red")

########################### #
# 7) Figure 3: violins + Tukey letters ####
########################### #
emm_df <- as.data.frame(emm_tp) |>
  mutate(
    emmean   = exp(emmean),
    lower.CL = exp(lower.CL),
    upper.CL = exp(upper.CL)
  )

# Letters used in the manuscript (edit if you change grouping)
letters_manual <- c("Dawn"="a","Noon"="b","Dusk"="b","Midnight"="c")
emm_df$group <- letters_manual[as.character(emm_df$time_point)]

span <- diff(range(df$pulsation_frequency, na.rm = TRUE))
emm_df$y_letters <- emm_df$upper.CL + 0.10 * span

emm_df$x_pos <- as.numeric(emm_df$time_point) + 0.08 

p_fig_lmm <- ggplot(df, aes(x = time_point, y = pulsation_frequency, fill = time_point)) +
  geom_violin(trim = FALSE, width = 0.8, alpha = 0.25, colour = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.90, colour = "grey25") +
  geom_point(aes(colour = time_point),
             position = position_jitter(width = 0.08, height = 0),
             size = 0.8, alpha = 0.25, show.legend = FALSE) +
  geom_text(data = emm_df,
            aes(x = x_pos, y = y_letters, label = group),
            colour = "red", size = 6, fontface = "bold",
            inherit.aes = FALSE) +
  scale_fill_manual(values = pal_tp, guide = "none") +
  scale_colour_manual(values = pal_tp, guide = "none") +
  labs(title = "Effect of the nychthemeral cycle on pulsation frequency",
       x = "", y = "Pulsation frequency (pulse/min)")
print(p_fig_lmm)

########################### #
# 8) Cyclic GAMM (24 h) + peak time CI + MESOR + amplitude (Figure 4) ####
########################### #

m_gam24 <- mgcv::gam(
  log(pulsation_frequency) ~ s(phase24, bs = "cc", # Cyclic cubic regression spline
                               k = 5) + # Number of basis functions
    z_umbrella +
    s(individual_id, 
      bs = "re") + # Random effect
    s(day, 
      bs = "re"), # Random effect
  data  = df,
  method = "REML",
  knots = list(phase24 = c(0, 24)) # Ensuring continuity at the 0–24 h boundary
)

# Diagnostics
summary(m_gam24)
mgcv::gam.check(m_gam24)

## ---- Population mean curve + 95% CI (exclude random effects) ----
pred24_gam <- data.frame(
  phase24    = seq(0, 24, length.out = 721),
  z_umbrella = mean(df$z_umbrella, na.rm = TRUE),
  individual_id = levels(df$individual_id)[1],
  day           = levels(df$day)[1]
)

pr24_gam <- predict(
  m_gam24, newdata = pred24_gam, se.fit = TRUE,
  exclude = c("s(individual_id)", "s(day)"),
  type = "link"
)

pred24_gam <- pred24_gam %>%
  mutate(
    fit = exp(pr24_gam$fit),
    lo  = exp(pr24_gam$fit - 1.96 * pr24_gam$se.fit),
    hi  = exp(pr24_gam$fit + 1.96 * pr24_gam$se.fit)
  )

## ---- Derived rhythm descriptors (from the fitted curve) ----
idx_max_gam <- which.max(pred24_gam$fit)
idx_min_gam <- which.min(pred24_gam$fit)

phase_peak_gam   <- pred24_gam$phase24[idx_max_gam]
phase_trough_gam <- pred24_gam$phase24[idx_min_gam]
phase_amplitude_gam    <- (pred24_gam$fit[idx_max_gam] - pred24_gam$fit[idx_min_gam]) / 2
phase_mesor_gam        <- mean(pred24_gam$fit)

cat(sprintf(
  "\nCircadian GAM (24 h, k=5): Peak=%.1f h, Trough=%.1f h, Amplitude=%.2f, MESOR=%.2f\n",
  phase_peak_gam, phase_trough_gam, phase_amplitude_gam, phase_mesor_gam
))

## ---- Coefficient simulation  ----
newdata_gam <- data.frame(
  phase24 = seq(0, 24, length.out = 900)
) %>%
  mutate(
    individual_id = df$individual_id[1],
    day           = df$day[1],
    z_umbrella    = mean(df$z_umbrella, na.rm = TRUE)
  )

# Population prediction (exclude random effects) on response scale
newdata_gam$fit <- predict(
  m_gam24, newdata = newdata_gam,
  exclude = c("s(individual_id)", "s(day)")
) |> exp()

peak_time_gam <- newdata_gam$phase24[which.max(newdata_gam$fit)]

# Coefficients + covariance matrix
beta_gam <- coef(m_gam24)
Vb_gam   <- vcov(m_gam24)

# Simulation of parameters
set.seed(123)     # ensures the same CI across runs
n_sim <- 1000
sim_beta_gam <- MASS::mvrnorm(n = n_sim, mu = beta_gam, Sigma = Vb_gam)

# Design matrix for the spline basis
Xp_gam <- predict(m_gam24, newdata_gam, type = "lpmatrix")


sim_fits_gam <- Xp_gam %*% t(sim_beta_gam) |> exp()

## ---- Morning peak 95% CI by coefficient simulation  ----
get_peak_gam <- function(sim) newdata_gam$phase24[which.max(sim)]
peak_times_gam <- apply(sim_fits_gam, 2, get_peak_gam)

peak_CI_gam <- quantile(peak_times_gam, c(0.025, 0.975))

c(
  peak  = mean(peak_times_gam),
  lower = peak_CI_gam[1],
  upper = peak_CI_gam[2]
)


df <- df %>% mutate(phase24_plot = if_else(phase24 == 0, 24, phase24))
pred24_gam <- pred24_gam %>% mutate(phase24_plot = if_else(phase24 == 0, 24, phase24))

peak_time_plot <- if_else(peak_time_gam == 0, 24, peak_time_gam)
peak_CI_plot   <- c(if_else(peak_CI_gam[1] == 0, 24, peak_CI_gam[1]),
                    if_else(peak_CI_gam[2] == 0, 24, peak_CI_gam[2]))

## ---- Nocturnal minimum 95% CI by coefficient simulation  ----
low_time_gam <- newdata_gam$phase24[which.min(newdata_gam$fit)]

get_low_gam <- function(sim) newdata_gam$phase24[which.min(sim)]
low_times_gam <- apply(sim_fits_gam, 2, get_low_gam)

low_CI_gam <- quantile(low_times_gam, c(0.025, 0.975))

c(
  low  = mean(low_times_gam),
  lower = low_CI_gam[1],
  upper = low_CI_gam[2]
)

low_time_plot <- if_else(low_time_gam == 0, 24, low_time_gam)
low_CI_plot   <- c(if_else(low_CI_gam[1] == 0, 24, low_CI_gam[1]),
                    if_else(low_CI_gam[2] == 0, 24, low_CI_gam[2]))

## ---- MESOR 95% CI by coefficient simulation  ----
mesor_gam <- mean(newdata_gam$fit)

mesors_gam <- apply(sim_fits_gam, 2, mean)

mesor_CI_gam <- quantile(mesors_gam, c(0.025, 0.975))

c(
  MESOR  = mean(mesors_gam),
  lower = mesor_CI_gam[1],
  upper = mesor_CI_gam[2]
)

## ---- Amplitude 95% CI by coefficient simulation  ----
idx_max_gam <- which.max(newdata_gam$fit)
idx_min_gam <- which.min(newdata_gam$fit)
amplitude_gam    <- (newdata_gam$fit[idx_max_gam] - newdata_gam$fit[idx_min_gam]) / 2


get_amplitude_gam <- function(sim) (max(sim) - min(sim)) / 2
amplitudes_gam <- apply(sim_fits_gam, 2, get_amplitude_gam)

amplitude_CI_gam <- quantile(amplitudes_gam, c(0.025, 0.975))

c(
  amplitude  = mean(amplitudes_gam),
  lower = amplitude_CI_gam[1],
  upper = amplitude_CI_gam[2]
)


## ---- Final figure  ----
ggplot(pred24_gam, aes(phase24_plot, fit)) +
  
  # 95% CI for peak time
  annotate(
    "rect",
    xmin = peak_CI_plot[1], xmax = peak_CI_plot[2],
    ymin = 0, ymax = max(df$pulsation_frequency, na.rm = TRUE),
    alpha = 0.13, fill = "pink"
  ) +
  
  # Peak (vertical line + label)
  geom_vline(xintercept = peak_time_plot, color = "red", linetype = "dashed") +
  annotate(
    "text",
    x = peak_time_plot, y = max(pred24_gam$fit) * 0.25,
    label = sprintf("Peak = %.1f h", peak_time_gam),
    vjust = 1, color = "red", size = 4
  ) +
  
  # 95% CI for low time
  annotate(
    "rect",
    xmin = low_CI_plot[1], xmax = low_CI_plot[2],
    ymin = 0, ymax = max(df$pulsation_frequency, na.rm = TRUE),
    alpha = 0.13, fill = "pink"
  ) +
  
  # Low (vertical line + label)
  geom_vline(xintercept = low_time_plot, color = "red", linetype = "dashed") +
  annotate(
    "text",
    x = low_time_plot-1, # Minus 1 to be in the figure frame
    y = max(pred24_gam$fit) * 0.25,
    label = sprintf("Low = %.1f h", low_time_gam),
    vjust = 1, color = "red", size = 4
  ) +
  
  # GAM curve + 95% CI
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#8F9779", alpha = 0.25) +
  geom_line(colour = "#556B2F", linewidth = 1.6) +
  
  # Raw data + geometric means (red)
  geom_point(data = df, aes(y = pulsation_frequency, x = phase24_plot), alpha = 0.45) +
  stat_summary(data = df, aes(y = pulsation_frequency, x = phase24_plot),
               fun = geommean, color = "red") +
  
  # MESOR line + label
  geom_hline(yintercept = phase_mesor_gam, colour = "#8B5E3C", linewidth = 1.5) +
  annotate(
    "text",
    x = 6.2,
    y = phase_mesor_gam - 4,
    label = sprintf("MESOR = %.2f", phase_mesor_gam),
    colour = "#8B5E3C",
    hjust = 0,
    size = 4
  ) +
  
  # Amplitude arrow + label
  annotate(
    "segment",
    x = peak_time_plot + 10, xend = peak_time_plot + 10,
    y = phase_mesor_gam - phase_amplitude_gam,
    yend = phase_mesor_gam + phase_amplitude_gam,
    colour = "#7E5AA7", linewidth = 1.2,
    arrow = arrow(ends = "both", length = unit(0.18, "cm"))
  ) +
  annotate(
    "text",
    x = peak_time_plot + 10,
    y = phase_mesor_gam + phase_amplitude_gam,
    label = sprintf("Amplitude = %.2f", phase_amplitude_gam),
    vjust = -1, colour = "#7E5AA7", size = 4
  ) +
  
  # Axis: display time points in the paper order
  scale_x_continuous(
    breaks = c(6, 12, 18, 24),
    labels = c("Dawn", "Noon", "Dusk", "Midnight"),
    limits = c(0, 24)
  ) +
  
  labs(
    title = "Circadian smooth (24 h, cyclic GAM k=5)",
    subtitle = sprintf("Peak = %.1f h (95%% CI %.1f–%.1f) | MESOR = %.2f | Amplitude = %.2f",
                       peak_time_gam, peak_CI_gam[1], peak_CI_gam[2], phase_mesor_gam, phase_amplitude_gam),
    x = "Circadian phase (h)",
    y = "Pulsation frequency (pulse/min)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

########################### #
# 9) Frequency analysis: Lomb–Scargle on detrended mean series (Figure 5) ####
########################### #

# rebuild t72 in the same way as above (phase24 via time_point)
df <- df |>
  mutate(
    phase24 = unname(tp_codes[as.character(time_point)]),
    day_num = as.integer(gsub("^D", "", as.character(day))) - 1L,
    t72 = day_num*24 + phase24,
    t72 = if_else(phase24 == 0, t72 + 24, t72)
  )

df_mean <- df |>
  group_by(t72) |>
  summarise(mean_freq = geommean(pulsation_frequency), .groups = "drop") |>
  arrange(t72)

# detrend first
y_dt <- resid(lm(mean_freq ~ t72, data = df_mean))

# Lomb–Scargle periodogram (period in hours)
ls <- lomb::lsp(
  x = y_dt,
  times = df_mean$t72,
  type = "period",
  from = 12, to = 72,
  ofac = 12,
  plot = FALSE
)

spec_df <- data.frame(period = ls$scanned, power = ls$power) |>
  mutate(power_rel = power / max(power, na.rm = TRUE))

peak <- spec_df |> slice_max(power, n = 1)
print(peak)

ggplot(spec_df, aes(x = period, y = power_rel)) +
  geom_line(linewidth = 1.2, colour = "#556B2F") +
  geom_point(size = 1.2, colour = "#556B2F") +
  geom_vline(xintercept = c(12, 24, 48), linetype = 3) +
  scale_x_continuous(limits = c(8, 72), breaks = seq(12, 72, 12)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0, 1, 0.2)) +
  labs(title = "", x = "Period (hours)", y = "Relative power")


########################### #
# 10) Interindividual variability: variance decomposition + rptR ####
########################### #

vc <- as.data.frame(VarCorr(m_final))
var_ind <- vc$vcov[vc$grp == "individual_id"]
var_day <- vc$vcov[vc$grp == "day"]
var_res <- attr(VarCorr(m_final), "sc")^2
var_tot <- var_ind + var_day + var_res

cat("\n--- Variance proportions (final LMM) ---\n")
cat(sprintf("Between individuals: %.3f\n", var_ind/var_tot))
cat(sprintf("Between days:        %.3f\n", var_day/var_tot))
cat(sprintf("Within individuals:  %.3f\n", var_res/var_tot))

df_rep <- df |>
  dplyr::select(pulsation_frequency, individual_id, day, time_point, z_umbrella) |>
  tidyr::drop_na()

rpt_res <- rptR::rptGaussian(
  formula = pulsation_frequency ~ time_point + z_umbrella + (1|individual_id) + (1|day),
  grname  = c("individual_id","day"),
  data    = df_rep,
  nboot   = 500,
  npermut = 0
)
print(summary(rpt_res))





########################################################### #
# 02_additional_code.R ####
# Additional analyses not shown in main text (redundant/robustness)
########################################################### #

library(readxl)
library(dplyr)
library(tidyr)
library(janitor)
library(ggplot2)
library(mgcv)
library(MASS)
library(lomb)

theme_set(theme_minimal(base_size = 13))
geommean <- function(x) exp(mean(log(x)))

path <- "cassiopea_pulsation_dataset.xlsx"
df <- read_excel(path) |> clean_names()

# (Re)build needed variables exactly as in main script
tp_codes <- c(Dawn = 6, Noon = 12, Dusk = 18, Midnight = 0)

df <- df |>
  mutate(
    individual_id = factor(individual_id),
    day = as.character(day),
    day = if_else(grepl("^D", day), day, paste0("D", day)),
    day = factor(day),
    time_point = factor(time_point, levels = c("Dawn","Noon","Dusk","Midnight")),
    pulsation_frequency  = as.numeric(pulsation_frequency),
    umbrella_diameter_cm = as.numeric(umbrella_diameter_cm),
    z_umbrella = as.numeric(scale(umbrella_diameter_cm)),
    phase24 = unname(tp_codes[as.character(time_point)]),
    day_num = as.integer(gsub("^D","", as.character(day))) - 1L,
    t72 = day_num*24 + phase24,
    t72 = if_else(phase24 == 0, t72 + 24, t72)
  )

########################### #
## A) GAM 72 h exploratory (not essential to manuscript figures) ####
########################### #
m_gam72 <- mgcv::gam(
  log(pulsation_frequency) ~ s(t72, 
                               k = 12, 
                               bs = "cr") +  # Cubic Regression Spline
    z_umbrella + s(individual_id, 
                   bs = "re") + # Random effect
    s(day, 
      bs = "re"),# Random effect
  data = df, 
  method = "REML" 
)
summary(m_gam72)
mgcv::gam.check(m_gam72)

########################### #
## B) Cosinor models (24h ± 12h harmonic) — redundant with GAMM ####
########################### #
df <- df |>
  mutate(
    w = 2*pi*(phase24 %% 24)/24,
    sin_w = sin(w), cos_w = cos(w)
  )

m_cosinor_24 <- mgcv::gam(
  log(pulsation_frequency) ~ sin_w + cos_w + z_umbrella +
    s(individual_id, bs = "re") + s(day, bs = "re"),
  data = df, method = "REML"
)
m_cosinor_24_12 <- mgcv::gam(
  log(pulsation_frequency) ~ sin_w + cos_w + sin(2*w) + cos(2*w) + z_umbrella +
    s(individual_id, bs = "re") + s(day, bs = "re"),
  data = df, method = "REML"
)

AIC(m_cosinor_24, m_cosinor_24_12)
anova(m_cosinor_24, m_cosinor_24_12, test = "Chisq")

########################### #
## C) Alternative spectral approach: Fourier on regularized mean series ####
########################### #
df_mean <- df |>
  group_by(t72) |>
  summarise(mean_freq = geommean(pulsation_frequency), .groups = "drop") |>
  arrange(t72)

y_dt <- resid(lm(mean_freq ~ t72, data = df_mean))
y_ts <- ts(y_dt) # (irregular sampling -> Lomb is preferred; this is just a check)
sp <- spec.pgram(y_ts, plot = FALSE)
sp

########################### #
## D) Exploratory scatterplots (environmental data) ####
########################### #

scatter_env <- function(df, xvar, xlab) {
  ggplot(df, aes(x = .data[[xvar]], y = pulsation_frequency, color = time_point)) +
    geom_point(alpha = 0.6, size = 1.6, na.rm = TRUE) +
    scale_color_manual(values = pal_tp) +
    labs(
      x = xlab,
      y = "Pulsation frequency (pulse/min)",
      color = "Time point"
    )
}

print(scatter_env(df, "oxygen_mg_l", "Oxygen (mg/L)"))
print(scatter_env(df, "ph", "pH"))
print(scatter_env(df, "temperature_c", "Temperature (°C)"))

