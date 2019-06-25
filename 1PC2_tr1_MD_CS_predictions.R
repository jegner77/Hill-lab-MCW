# 
# Chemical shift predictions (Sparta+) of MD simulations
# experiment = 100-ns 1PC2_121SKY, tr1
# force field = AMBER99SB
# MD simulation time = 100 ns
# GROMACS
# 

# Load libraries, Import data, and tidy data  ----
library(tidyverse)
library(broom)
library(readxl)
library(sheetr)

data_path <- "./pred_shifts/" # Update to path of tab files

files <- dir(data_path, pattern = "*.tab")

pred <- as_tibble(data_frame(snapshot = files)) %>%
  mutate(data = map(snapshot, 
                             function(x) read_table2(file.path(data_path, x),
                                                     col_names = F, skip = 26)),
         data = map(data, ~slice(., 14:nrow(.))),
         data = map(data, ~rename(., residue = X1,
                                  aa = X2,
                                  nuclei = X3,
                                  ss_shift = X4,
                                  predicted_shift = X5,
                                  rc_shift = X6,
                                  hm_shift = X7,
                                  ef_shift = X8,
                                  sigma = X9)),
         data = map(data, ~mutate(., residue = as.numeric(residue),
                                  ss_shift = as.numeric(ss_shift),
                                  predicted_shift = as.numeric(predicted_shift),
                                  rc_shift = as.numeric(rc_shift),
                                  hm_shift = as.numeric(hm_shift),
                                  ef_shift = as.numeric(ef_shift),
                                  sigma = as.numeric(sigma))),
         data = map(data, ~select(., -X10, -X11)),
         snapshot = sapply(strsplit(snapshot, split = "_"),
                           function(x) (x[2])),
         snapshot = parse_number(snapshot)) %>%
  arrange(snapshot)

# Analyze HN predicted chemical shifts to experimental values ----

exp_ref <- read_excel("apo_hFis1_aa1-125_chemical_shifts.xlsx") %>%
  gather(., key = "nuclei", value = "exp_shift", 2:3) %>%
  mutate(., label = parse_number(label),
         nuclei = if_else(exp_shift < 20, "HN", "N")) %>%
  rename(., residue = label)

HN_shifts <- pred %>%
  unnest() %>%
  left_join(., exp_ref, by = c("residue", "nuclei")) %>%
  filter(., nuclei == "HN") %>%
  nest(-residue) %>%
  mutate(CS_by_time = map2(.x = data, .y = residue,
                     ~ggplot(data = .x, aes(x = snapshot/100, y = predicted_shift)) +
                       geom_point(size = 0.5) +
                       geom_line(stat = "smooth", method = "loess", span = 0.2,
                                 color = "purple", size = 1.5) +
                       geom_hline(aes(yintercept = exp_shift), color = "cyan",
                                  size = 1) +
                       labs(title = .y,
                            x = "Time (ns)",
                            y = "Predicted HN CS (ppm)",
                            caption = "purple = local regression of predicted CS,
                            cyan = experimental CS in HEPES") +
                       theme_bw()),
         population = map2(.x = data, .y = residue,
                           ~ggplot(data = .x, aes(x = predicted_shift)) +
                             geom_density(size = 1) +
                             geom_point(aes(x = mean(predicted_shift), y = 0),
                                        size = 5, color = "purple") +
                             geom_point(aes(x = exp_shift, y = 0),
                                        size = 5, color = "cyan") +
                             labs(title = .y,
                                  x = "CS (ppm)",
                                  y = "Population",
                                  caption = "purple = mean value of predicted CS,
                                            cyan = experimental CS in HEPES") +
                             theme_bw())
  )

map2(.x = paste0(HN_shifts$residue, "_HN_time.pdf"), .y = HN_shifts$CS_by_time,
     ggsave)
map2(.x = paste0(HN_shifts$residue, "_HN_population.pdf"), .y = HN_shifts$population,
     ggsave)

# Analyze N predicted chemical shifts to experimental values ----

N_shifts <- pred %>%
  unnest() %>%
  left_join(., exp_ref, by = c("residue", "nuclei")) %>%
  filter(., nuclei == "N") %>%
  nest(-residue) %>%
  mutate(CS_by_time = map2(.x = data, .y = residue,
                           ~ggplot(data = .x, aes(x = snapshot/100, 
                                                  y = predicted_shift)) +
                             geom_point(size = 0.5) +
                             geom_line(stat = "smooth", method = "loess", span = 0.2,
                                       color = "purple", size = 1.5) +
                             geom_hline(aes(yintercept = exp_shift), color = "cyan",
                                        size = 1) +
                             scale_x_continuous(limits = c(0, 100),
                                        breaks = c(0, 20, 40, 60, 80, 100)) +
                             scale_y_continuous(limits = c(105, 130),
                                        breaks = c(105, 110, 115, 120, 125, 130)) +
                             labs(title = .y,
                               x = "Time (ns)",
                               y = "Predicted N CS (ppm)",
                               caption = "purple = local regression of predicted CS,
                                          cyan = experimental CS in HEPES") +
                             theme_bw()),
         population = map2(.x = data, .y = residue,
                           ~ggplot(data = .x, aes(x = predicted_shift)) +
                             geom_density(size = 1) +
                             geom_point(aes(x = mean(predicted_shift), y = 0),
                                        size = 5, color = "purple") +
                             geom_point(aes(x = exp_shift, y = 0),
                                        size = 5, color = "cyan") +
                             scale_x_continuous(limits = c(min(.$predicted_shift)-4, 
                                                           max(.$predicted_shift)+4),
                                                breaks = c(106, 108, 110, 112,
                                                           114, 116, 118, 120,
                                                           122, 124, 126, 128)) +
                             labs(title = .y,
                                  x = "N CS (ppm)",
                                  y = "Population",
                                  caption = "purple = mean value of predicted CS,
                                  cyan = experimental CS in HEPES") +
                             theme_bw())
         )

map2(.x = paste0(N_shifts$residue, "_N_time.pdf"), .y = N_shifts$CS_by_time,
     ggsave)
map2(.x = paste0(N_shifts$residue, "_N_population.pdf"), .y = N_shifts$population,
     ggsave)

# Calculate mean and sd for predicted chemical shifts ----

comp_1pc2 <- pred %>%
  unnest() %>%
  left_join(., exp_ref, by = c("residue", "nuclei")) %>%
  group_by(., nuclei, residue) %>%
  summarise(., avg_1pc2 = mean(predicted_shift),
            sd_1pc2 = sd(predicted_shift),
            avg_exp = mean(exp_shift),
            diff_1pc2 = abs(avg_1pc2 - avg_exp)) %>%
  left_join(., exp_ref, by = c("residue", "nuclei")) %>%
  select(-exp_shift) %>%
  filter(., nuclei == "HN" | nuclei == "N")

write_csv(comp_1pc2, "comparison_1pc2_to_exp.csv")

# Plotting chemical shift differences for HN nuclei ----
comp_1pc2 %>%
  filter(., nuclei == "HN") %>%
  ggplot(., aes(x = residue, y = diff_1pc2)) +
  geom_bar(aes(fill = secondary), stat = "identity") +
  geom_hline(aes(yintercept = 0.49), color = "red", linetype = 2) +
  labs(title = "HN CS Differences of 100 ns 1PC2 MD & Experimental",
       fill = "Secondary Structure",
       x = "Residue",
       y = expression(paste("|CS"["1PC2"], " - CS"["Exp"],"|"))) + 
  scale_fill_manual(breaks = c("arm", "helix", "loop"),
                    labels = c("N-arm", "Helix", "Loop"),
                    values = c("#1b9e77", "#7570b3", "#d95f02")) +
  theme_bw() +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "top")

ggsave("comparison_1pc2_exp_HN.pdf",
       width = 14, height = 10, dpi = 300, units = "cm")

# Plotting chemical shift differences for N nuclei ----
comp_1pc2 %>%
  filter(., nuclei == "N") %>%
  ggplot(., aes(x = residue, y = diff_1pc2)) +
  geom_bar(aes(fill = secondary), stat = "identity") +
  geom_hline(aes(yintercept = 2.45), color = "red", linetype = 2) +
  labs(title = "N CS Differences of 100 ns 1PC2 MD & Experimental",
       fill = "Secondary Structure",
       x = "Residue",
       y = expression(paste("|CS"["1PC2"], " - CS"["Exp"],"|"))) + 
  scale_fill_manual(breaks = c("arm", "helix", "loop"),
                    labels = c("N-arm", "Helix", "Loop"),
                    values = c("#1b9e77", "#7570b3", "#d95f02")) +
  theme_bw() +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "top")

ggsave("comparison_1pc2_exp_N.pdf",
       width = 14, height = 10, dpi = 300, units = "cm")

# Calculate difference of differences between 1PC2-EXP and 1IYG-EXP ----
  ## Is experimental chemical shifts closer to arm-out (1PC2) or arm-in (1IYG)?

comparison <- read_csv("comparison_1iyg_to_exp.csv") %>%
  left_join(., comp_1pc2, by = c("residue", "nuclei", "secondary", "avg_exp")) %>%
  group_by(., nuclei, residue) %>%
  mutate(., arm_in_out = diff_1pc2 - diff_1iyg)
  # arm_in_out = absolute difference(1PC2-EXP shifts) minus 
  # absolute difference(1IYG-EXP shifts)
    # thus, positive values indicate experimental shift closer to Arm-in (or 1IYG)
    # and negative values indicate experimental shift closer to Arm-out (or 1PC2)

# Plotting chemical shift differences between 1PC2/1IYG and Exp for HN nuclei ----
comparison %>%
  filter(., nuclei == "HN") %>%
  ggplot(., aes(x = residue, y = arm_in_out)) +
  geom_bar(aes(fill = secondary), stat = "identity") +
  geom_hline(aes(yintercept = -0.49), color = "red", linetype = 2) +
  geom_hline(aes(yintercept = 0.49), color = "red", linetype = 2) +
  labs(title = "HN CS Differences",
       fill = "Secondary Structure",
       x = "Residue",
       y = expression(paste("|CS"["1PC2"], " - CS"["Exp"],"|", " - ",
                            "|CS"["1IYG"], " - CS"["Exp"],"|"))) + 
  scale_fill_manual(breaks = c("arm", "helix", "loop"),
                    labels = c("N-arm", "Helix", "Loop"),
                    values = c("#1b9e77", "#7570b3", "#d95f02")) +
  theme_bw() +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "top")
  # red dotted line marks Sparta+ prediction accuracy for HN nuclei

ggsave("comparison_1pc2_1iyg_exp_HN.pdf",
       width = 14, height = 10, dpi = 300, units = "cm")

# Plotting chemical shift differences between 1PC2/1IYG and Exp for N nuclei ----
comparison %>%
  filter(., nuclei == "N") %>%
  ggplot(., aes(x = residue, y = arm_in_out)) +
  geom_bar(aes(fill = secondary), stat = "identity") +
  geom_hline(aes(yintercept = -2.45), color = "red", linetype = 2) +
  geom_hline(aes(yintercept = 2.45), color = "red", linetype = 2) +
  labs(title = "N CS Differences",
       fill = "Secondary Structure",
       x = "Residue",
       y = expression(paste("|CS"["1PC2"], " - CS"["Exp"],"|", " - ",
                            "|CS"["1IYG"], " - CS"["Exp"],"|"))) +
  scale_fill_manual(breaks = c("arm", "helix", "loop"),
                    labels = c("N-arm", "Helix", "Loop"),
                    values = c("#1b9e77", "#7570b3", "#d95f02")) +
  theme_bw() +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "top")
  # red dotted line marks Sparta+ prediction accuracy for N nuclei

ggsave("comparison_1pc2_1iyg_exp_N.pdf",
       width = 14, height = 10, dpi = 300, units = "cm")
