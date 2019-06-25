#
# Data analysis of hetNOE of hFis1 aa1-125 in 1IYG sample conditions
#


# Load libraries ---------------------------------------------------------------

library(tidyverse)
library(broom)
library(readxl)
theme_set(theme_bw() +
            theme(axis.text = element_text(size = 12, color = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())
          )

# Import, tidy, and calculate hetNOE value -------------------------------------

hetNOE <- read_tsv("hetnoe_1iyg_condition", col_names = F) %>%
  select(-X1) %>%
  t() %>%
  as_tibble(.) %>%
  slice(1:118) %>%
  rename(residue = V1,
         noe = V2,
         reference = V3) %>%
  mutate(het = noe/reference)

hetNOE

apo_hepes <- read_excel("apo_hFis1_aa1-125_chemical_shifts.xlsx") %>%
  mutate(residue = parse_number(label)) %>%
  select(-label, -H_ppm, -N_ppm)

apo_hepes

hetNOE <- hetNOE %>%
  left_join(apo_hepes, by = "residue")
hetNOE

# Plot hetNOE values for hFis1 in 1IYG condition -------------------------------

hetNOE %>%
  ggplot(aes(x = residue, y = het, fill = secondary)) +
  geom_col() +
  labs(title = "hetNOE for hFis1 in 1IYG conditions",
       x = "Residue",
       y = "1H, 15N hetNOE") +
  scale_fill_manual(breaks = c("arm", "helix", "loop"),
                    labels = c("N-arm", "Helix", "Loop"),
                    values = c("#1b9e77", "#7570b3", "#d95f02")) +
  scale_y_continuous(limits = c(-1.1, 1.1),
                     breaks = c(-1, -0.5, 0, 0.5, 1)) +
  theme(legend.position = "top")

ggsave("hetNOE_hFis1_1iyg_condition.pdf",
       width = 12, height = 8, units = "cm", dpi = 300)

# Data transformations to set-up for chemcial shift ∆∆∆ plot -------------------

# import chemical shifts from HEPES and 1PC2 conditions
shifts_hepes <- read_excel("apo_hFis1_aa1-125_chemical_shifts.xlsx") %>%
  mutate(residue = parse_number(label)) %>%
  rename(H_hepes = H_ppm,
         N_hepes = N_ppm)
shifts_hepes

shifts_1pc2 <- read_csv("hFis1_1pc2_condition_chemical_shifts.csv") %>%
  rename(H_1pc2 = H,
         N_1pc2 = N)
shifts_1pc2

# Generate H shift table
H_shift <- read_delim("H_1iyg.txt", delim = ";") %>%
  rename(label = `Spin label`,
         H_1iyg = hetnoe_noe.3D) %>%
  select(-hetnoe_ref.3D) %>%
  mutate(residue = parse_number(label),
         H_1iyg = parse_number(H_1iyg)) %>%
  left_join(shifts_hepes, by = "residue") %>%
  left_join(shifts_1pc2, by = "residue") %>%
  select(residue, secondary, H_hepes, H_1iyg, H_1pc2) %>%
  gather(3:5, key = "tmp1", value = "H_shift") %>%
  separate(tmp1, into = c("tmp2", "type"), sep = "_") %>%
  select(-tmp2)
  
H_shift

# Generate N shift table and combine with H shift table
shifts <- read_delim("N_1iyg.txt", delim = ";") %>%
  rename(label = `Spin label`,
         N_1iyg = hetnoe_noe.3D) %>%
  select(-hetnoe_ref.3D) %>%
  mutate(residue = parse_number(label)) %>%
  left_join(shifts_hepes, by = "residue") %>%
  left_join(shifts_1pc2, by = "residue") %>%
  select(residue, secondary, N_hepes, N_1iyg, N_1pc2) %>%
  gather(3:5, key = "tmp1", value = "N_shift") %>%
  separate(tmp1, into = c("tmp2", "type"), sep = "_") %>%
  left_join(H_shift, by = c("residue", "type", "secondary")) %>%
  select(residue, secondary, H_shift, N_shift, type)

shifts

write_csv(shifts, "hFis1_all_3_conditions_chemical_shifts.csv")

# Calculate and plot chemical shift ∆∆∆ ----------------------------------------

delta3 <- shifts %>%
  group_by(residue) %>%
  mutate(human_delta = sqrt(((5*(H_shift[type == "1pc2"] - H_shift[type == "hepes"]))^2) + 
                              (N_shift[type == "1pc2"] - N_shift[type == "hepes"])^2),
         mouse_delta = sqrt(((5*(H_shift[type == "1iyg"] - H_shift[type == "hepes"]))^2) + 
                              (N_shift[type == "1iyg"] - N_shift[type == "hepes"])^2),
         delta_delta = human_delta - mouse_delta) 

delta3

delta3 %>%
  select(residue, secondary, delta_delta) %>%
  ungroup() %>%
  distinct(delta_delta, .keep_all = T) %>%
  ggplot(aes(x = residue, y = delta_delta, fill = secondary)) +
  geom_col() +
  scale_y_continuous(limits = c(-0.65, 0.4),
                     breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4)) +
  labs(fill = "Secondary Structure:",
       x = "Residue",
       y = expression(paste("|CS"["OUT"], " - CS"["Exp"],"| -","|CS"["IN"], " - CS"["Exp"],"|"))) + 
  scale_fill_manual(breaks = c("arm", "helix", "loop"),
                    labels = c("N-arm", "Helix", "Loop"),
                    values = c("#1b9e77", "#7570b3", "#d95f02")) +
  theme(legend.position = "top")

ggsave("delta_delta_delta_CS_conditions.pdf",
       width = 12, height = 8, units = "cm", dpi = 300)

delta3 %>%
  select(residue, secondary, human_delta, mouse_delta) %>%
  ungroup() %>%
  distinct(human_delta, .keep_all = T) %>%
  ggplot(aes(x = residue)) +
  geom_line(aes(y = human_delta), color = "purple") +
  geom_line(aes(y = mouse_delta), color = "green") +
  labs(title = "CS Difference between Human/Mouse & HEPES",
       subtitle = "∆Human/Hepes = purple \n∆Mouse/Hepes = green",
       fill = "Secondary Structure:",
       x = "Residue",
       y = expression(paste("|CS"["Human or Mouse"], " - CS"["HEPES"],"|"))) + 
  scale_fill_manual(breaks = c("arm", "helix", "loop"),
                    labels = c("N-arm", "Helix", "Loop"),
                    values = c("#1b9e77", "#7570b3", "#d95f02")) +
  theme(legend.position = "top")

ggsave("delta_delta_delta_CS_conditions_line_plot.pdf",
       width = 12, height = 8, units = "cm", dpi = 300)

# Plot human and mouse deltas against EXP, individually ------------------------

delta3 %>%
  select(residue, secondary, human_delta, mouse_delta) %>%
  ungroup() %>%
  distinct(human_delta, .keep_all = T) %>%
  ggplot(aes(x = residue, y = human_delta, fill = secondary)) +
  geom_col() +
  labs(title = "CS Difference between Human & HEPES",
       fill = "Secondary Structure:",
       x = "Residue",
       y = expression(paste("|CS"["Human"], " - CS"["HEPES"],"|"))) +
  scale_fill_manual(breaks = c("arm", "helix", "loop"),
                    labels = c("N-arm", "Helix", "Loop"),
                    values = c("#1b9e77", "#7570b3", "#d95f02")) +
  theme(legend.position = "top")

ggsave("delta_1pc2_hepes_CS_conditions.pdf",
       width = 12, height = 8, units = "cm", dpi = 300)

delta3 %>%
  select(residue, secondary, human_delta, mouse_delta) %>%
  ungroup() %>%
  distinct(human_delta, .keep_all = T) %>%
  ggplot(aes(x = residue, y = mouse_delta, fill = secondary)) +
  geom_col() +
  labs(title = "CS Difference between Mouse & HEPES",
       fill = "Secondary Structure:",
       x = "Residue",
       y = expression(paste("|CS"["Mouse"], " - CS"["HEPES"],"|"))) +
  scale_fill_manual(breaks = c("arm", "helix", "loop"),
                    labels = c("N-arm", "Helix", "Loop"),
                    values = c("#1b9e77", "#7570b3", "#d95f02")) +
  theme(legend.position = "top")

ggsave("delta_1iyg_hepes_CS_conditions.pdf",
       width = 12, height = 8, units = "cm", dpi = 300)


# Import hetNOE from HEPES condition and plot ----------------------------------

hetNOE_hepes <- read_tsv("hetNOE_hFis1_hepes", col_names = F) %>%
  select(-X1) %>%
  t() %>%
  as_tibble(.) %>%
  slice(1:118) %>%
  rename(residue = V1,
         noe = V3,
         reference = V2) %>%
  mutate(het_hepes = noe/reference) %>%
  left_join(apo_hepes, by = "residue") %>%
  filter(residue != 19) # residue 19 excluded due to xpk overlap

hetNOE_hepes

hetNOE_hepes %>%
  ggplot(aes(x = residue, y = het_hepes, fill = secondary)) +
  geom_col() +
  labs(title = "hetNOE for hFis1 in HEPES",
       x = "Residue",
       y = "1H, 15N hetNOE") +
  scale_fill_manual(breaks = c("arm", "helix", "loop"),
                    labels = c("N-arm", "Helix", "Loop"),
                    values = c("#1b9e77", "#7570b3", "#d95f02")) +
  scale_y_continuous(limits = c(-1.1, 1.1),
                     breaks = c(-1, -0.5, 0, 0.5, 1)) +
  theme(legend.position = "top")

ggsave("hetNOE_hFis1_HEPES_condition.pdf",
       width = 12, height = 8, units = "cm", dpi = 300)

# Calculate and Plot delta-delta-delta for hetNOE values -----------------------

hetNOE_1pc2 <- read_tsv("hFis1_1pc2_hetNOE_noe2.txt", col_names = F) %>%
  select(-X1) %>%
  t() %>%
  as_tibble(.) %>%
  slice(1:118) %>%
  rename(residue = V1,
         noe = V2,
         reference = V3) %>%
  mutate(het = noe/reference) %>%
  filter(residue != 47) %>%
  rename(het_1pc2 = het) %>%
  select(residue, het_1pc2)
hetNOE_1pc2

hetNOE_1pc2 %>%
  ggplot(aes(x = residue, y = het_1pc2)) +
  geom_col() +
  labs(title = "hetNOE for hFis1 in 1PC2",
       x = "Residue",
       y = "1H, 15N hetNOE") +
  scale_fill_manual(breaks = c("arm", "helix", "loop"),
                    labels = c("N-arm", "Helix", "Loop"),
                    values = c("#1b9e77", "#7570b3", "#d95f02")) +
  scale_y_continuous(limits = c(-1.1, 1.1),
                     breaks = c(-1, -0.5, 0, 0.5, 1)) +
  theme(legend.position = "top")

ggsave("hetNOE_hFis1_1pc2_condition.pdf",
       width = 12, height = 8, units = "cm", dpi = 300)

hetNOE_delta <- hetNOE_hepes %>%
  select(residue, secondary, het_hepes) %>%
  left_join(hetNOE, by = c("residue", "secondary")) %>%
  select(-noe, -reference) %>%
  rename(het_1iyg = het) %>%
  left_join(hetNOE_1pc2, by = "residue") %>%
  mutate(human_het_delta = abs(het_1pc2 - het_hepes),
         mouse_het_delta = abs(het_1iyg - het_hepes),
         het_delta = human_het_delta - mouse_het_delta)

hetNOE_delta

hetNOE_delta %>%
  ggplot(aes(x = residue, y = het_delta, fill = secondary)) +
  geom_col() +
  scale_y_continuous(limits = c(-0.25, 0.4),
                     breaks = c(-0.2, 0, 0.2, 0.4)) +
  labs(fill = "Secondary Structure:",
       x = "Residue",
       y = expression(paste("|NOE"["OUT"], " - NOE"["Exp"],"| -","|NOE"["IN"], " - NOE"["Exp"],"|"))) + 
  scale_fill_manual(breaks = c("arm", "helix", "loop"),
                    labels = c("N-arm", "Helix", "Loop"),
                    values = c("#1b9e77", "#7570b3", "#d95f02")) +
  theme(legend.position = "top")

ggsave("delta_delta_delta_CS_conditions.pdf",
       width = 12, height = 8, units = "cm", dpi = 300)

# Calculate average hetNOE values by secondary structure ------------------

avg_hetNOE <- hetNOE_delta %>%
  gather(het_hepes, het_1iyg, het_1pc2, key = "condition", value = "hetNOE") %>%
  group_by(secondary, condition) %>%
  summarise(avg = mean(hetNOE, na.rm = T),
            stdev = sd(hetNOE, na.rm = T)) %>%
  ungroup()

avg_hetNOE

write_csv(avg_hetNOE, "average_hetNOE_by_secondary_structure.csv")

avg_hetNOE %>%
  ggplot(aes(x = secondary, y = avg, fill = condition)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = avg - stdev, ymax = avg + stdev), 
                position = position_dodge(width = 0.9), width = 0.2) +
  scale_fill_brewer(palette = "Dark2",
                    labels = c("1IYG", "1PC2", "PHYS")) +
  labs(title = "hetNOE by secondary structure",
       x = "",
       y = "hetNOE value")

ggsave("average_hetNOE_by_secondary.pdf",
       width = 12, height = 8, units = "cm", dpi = 300)


stats_hetNOE <- hetNOE_delta %>%
  gather(het_hepes, het_1iyg, het_1pc2, key = "condition", value = "hetNOE") %>%
  do(tidy(aov(hetNOE ~ condition*secondary, data = .)))
  
stats_hetNOE  

tukey_hetNOE <- hetNOE_delta %>%
  gather(het_hepes, het_1iyg, het_1pc2, key = "condition", value = "hetNOE") %>%
  do(tidy(TukeyHSD(aov(hetNOE ~ condition*secondary, data = .))))

tukey_hetNOE  

write_csv(stats_hetNOE, "anova_hetNOE_values.csv")

write_csv(tukey_hetNOE, "tukey_posthoc_hetNOE_values.csv")
