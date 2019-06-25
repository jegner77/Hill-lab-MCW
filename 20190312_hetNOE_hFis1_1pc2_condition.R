#
# Data analysis of hetNOE of hFis1 aa1-125 in 1PC2 sample conditions
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

hetNOE <- read_tsv("hFis1_1pc2_hetNOE_noe2.txt", col_names = F) %>%
  select(-X1) %>%
  t() %>%
  as_tibble(.) %>%
  slice(1:118) %>%
  rename(residue = V1,
         noe = V2,
         reference = V3) %>%
  mutate(het = noe/reference) %>%
  filter(residue != 47)

hetNOE

apo_hepes <- read_excel("apo_hFis1_aa1-125_chemical_shifts.xlsx") %>%
  mutate(residue = parse_number(label)) %>%
  select(-label, -H_ppm, -N_ppm)

apo_hepes

hetNOE <- hetNOE %>%
  left_join(apo_hepes, by = "residue")
hetNOE

# Plot hetNOE values for hFis1 in 1PC2 condition -------------------------------

hetNOE %>%
  ggplot(aes(x = residue, y = het, fill = secondary)) +
  geom_col() +
  labs(title = "hetNOE for hFis1 in 1PC2 conditions",
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

# Chemical Shift Analysis of hFis1 in 1PC2 condition ---------------------------

H_shift <- read_delim("H_1pc2_conditions_chemical_shifts.txt", delim = ";") %>%
  rename(label = `Spin label`,
         H = hetnoe_noe.3D) %>%
  select(-hetnoe_ref.3D) %>%
  mutate(residue = parse_number(label),
         H = parse_number(H))
H_shift

shifts <- read_delim("N_1pc2_conditions_chemical_shifts.txt", delim = ";") %>%
  rename(label = `Spin label`,
         N = hetnoe_noe.3D) %>%
  select(-hetnoe_ref.3D) %>%
  mutate(residue = parse_number(label)) %>%
  left_join(H_shift, by = "residue") %>%
  select(residue, label.x, H, N) %>%
  rename(label = label.x)

shifts

write_csv(shifts, "hFis1_1pc2_condition_chemical_shifts.csv")

