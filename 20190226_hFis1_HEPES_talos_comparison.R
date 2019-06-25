# 
# Comparing ramachandran angles of hFis1 from experimental HEPES condition
# against NMR solution structures of N-arm in (1IYG) and out (1PC2)
# mouse_ref = 1IYG = N-arm in
# human_ref = 1PC2 = N-arm out
#


# Load libraries and set theme --------------------------------------------

library(tidyverse)
library(readxl)
library(bio3d)
theme_set(theme_bw() +
            theme(axis.text = element_text(size = 12, color = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())
          )

# Import prot file and aa sequence ----------------------------------------
hepes_rama <- read_table2("hFis1_HEPES_rama.txt") %>%
  mutate(type = "hepes") %>%
  select(RESID, RESNAME, PHI, PSI, type) %>%
  rename(residue = RESID,
         aa = RESNAME,
         phi = PHI,
         psi = PSI)
hepes_rama

mouse_ref_tmp <- torsion.pdb(read.pdb("1IYG_121SKY.pdb"))
mouse_ref <- as_tibble(mouse_ref_tmp$tbl) %>%
  mutate(residue = 1:121,
         type = "mouse",
         aa = hepes_rama$aa)
mouse_ref

human_ref_tmp <- torsion.pdb(read.pdb("1PC2_121SKY.pdb"))
human_ref <- as_tibble(human_ref_tmp$tbl) %>%
  mutate(residue = 1:121,
         type = "human",
         aa = hepes_rama$aa)
human_ref

human_exp <- read_table2("youle_pred.tab") %>% # talos of 1PC2 CS
  mutate(type = "human_exp") %>%
  select(RESID, RESNAME, PHI, PSI, type) %>%
  rename(residue = RESID,
         aa = RESNAME,
         phi = PHI,
         psi = PSI) 
human_exp

reference <- mouse_ref %>%
  union(human_ref) %>%
  select(residue, aa, phi, psi, type) %>%
  union(hepes_rama) %>%
  union(human_exp) %>%
  filter(psi != 9999 & phi != 9999)
reference

# Ramachandran Plot -------------------------------------------------------

### Ramachandran plot for all residues
reference %>%
  ggplot(aes(x = phi, y = psi, color = type)) +
  geom_jitter(alpha = 0.4, size = 1.5, shape = 15) +
  scale_color_brewer(palette = "Dark2",
                     labels = c("EXP", "HUMAN ('OUT')", 
                                "HUMAN (∂1PC2)", "MOUSE ('IN')")) +
  scale_x_continuous(limits = c(-180, 180),
                     breaks = c(-180, -90, 0, 90, 180)) +
  scale_y_continuous(limits = c(-180, 180),
                     breaks = c(-180, -90, 0, 90, 180)) +
  labs(color = "",
       x = "PHI",
       y = "PSI") +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(legend.position = "top",
        axis.text = element_text(size = 10, color = "black"))

ggsave("all_rama_comparison.pdf",
       width = 10, height = 10, units = "cm", dpi = 300)

### Ramachandran plot for N-arm residues only
reference %>%
  filter(residue <= 13) %>% # residue 11 is where helix 1 starts in both
  ggplot(aes(x = phi, y = psi, color = type, shape = type)) +
  geom_jitter(size = 2, alpha = 1) +
  facet_wrap(~ residue) + 
  scale_color_brewer(palette = "Dark2",
                     labels = c("EXP", "'OUT' PDB Ref.", 
                                "'OUT' CS", "'IN' PDB Ref."),
                     guide = guide_legend(override.aes = list(size = 2))) +
  scale_shape_manual(labels = c("EXP", "'OUT' PDB Ref.",
                                "'OUT' CS", "'IN' PDB Ref."),
                     values = c(4, 15, 16, 24)) +
  scale_x_continuous(limits = c(-180, 180),
                     breaks = c(-120, 0, 120)) +
  scale_y_continuous(limits = c(-180, 180),
                     breaks = c(-120, 0, 120)) +
  labs(title = "Facet by Residue #",
       color = "",
       shape = "",
       x = "PHI",
       y = "PSI") +
  theme(legend.position = "top",
        axis.text = element_text(size = 10, color = "black"),
        strip.background = element_rect(fill = "white"))

ggsave("residues2-13_rama_comparison_by_colorSHAPE.svg",
       width = 12, height = 12, units = "cm", dpi = 300)

# Ramachandran plot for each residue --------------------------------------

reference %>%
  ggplot(aes(x = phi, y = psi, color = type)) +
  geom_jitter(shape = 15) +
  facet_wrap(~ residue) + 
  scale_color_brewer(palette = "Dark2",
                     labels = c("EXP", "'OUT' PDB Ref.", 
                                "'OUT' CS", "'IN' PDB Ref.")) +
  scale_x_continuous(limits = c(-180, 180),
                     breaks = c(-120, 0, 120)) +
  scale_y_continuous(limits = c(-180, 180),
                     breaks = c(-120, 0, 120)) +
  labs(title = "Facet by Residue #",
       color = "",
       x = "PHI",
       y = "PSI") +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(legend.position = "top",
        axis.text = element_text(size = 10, color = "black"),
        strip.background = element_rect(fill = "white"))

ggsave("all_rama_comparison_facet_by_residue.pdf",
       width = 30, height = 30, units = "cm", dpi = 300)


# PHI/PSI difference plot between EXP and reference structures ------------
apo <- read_excel("apo_hFis1_aa1-125_chemical_shifts.xlsx", sheet = 1) %>%
  mutate(residue = parse_number(label))
apo

phi_psi <- reference %>%
  group_by(residue) %>%
  filter(type != "human_exp") %>%
  mutate(human_delta = sqrt(((phi[type == "human"] - phi[type == "hepes"])^2) + 
                              (psi[type == "human"] - psi[type == "hepes"])^2),
         mouse_delta = sqrt(((phi[type == "mouse"] - phi[type == "hepes"])^2) + 
                              (psi[type == "mouse"] - psi[type == "hepes"])^2),
         delta_delta = human_delta - mouse_delta) %>%
  left_join(apo, by = "residue") %>%
  select(-H_ppm, -N_ppm, -label)
phi_psi

# delta-delta-delta plot of phi-psi angles between human and mouse
phi_psi %>%
  select(residue, aa, delta_delta, secondary) %>%
  distinct() %>%
  drop_na() %>%
  ggplot(aes(x = residue, y = delta_delta, fill = secondary)) +
  geom_col() +
  labs(fill = "Secondary Structure:",
       x = "Residue",
       y = expression(paste("|PHI/PSI"["OUT"], " - PHI/PSI"["Exp"],"| -","|PHI/PSI"["IN"], " - PHI/PSI"["Exp"],"|"))) + 
  scale_fill_manual(breaks = c("arm", "helix", "loop"),
                    labels = c("N-arm", "Helix", "Loop"),
                    values = c("#1b9e77", "#7570b3", "#d95f02")) +
  theme(legend.position = "top")

ggsave("delta_delta_phi_psi_human_mouse_exp.pdf",
       width = 12, height = 8, units = "cm")

phi_psi %>%
  select(residue, aa, human_delta, mouse_delta) %>%
  distinct() %>%
  ggplot(aes(x = residue)) +
  geom_point(aes(y = human_delta), color = "green") +
  geom_point(aes(y = mouse_delta), color = "purple") +
  geom_line(aes(y = human_delta), color = "green") +
  geom_line(aes(y = mouse_delta), color = "purple") +
  labs(x = "Residue",
       y = "Delta PHI/PSI")

ggsave("delta_phi_psi_human_mouse_exp.pdf",
       width = 12, height = 8, units = "cm")

# PHI/PSI difference plot between EXP & 1PC2 chemical shift derived angles -----
phi_psi_youle <- reference %>%
  group_by(residue) %>%
  filter(type != "mouse" & residue > 3) %>%
  mutate(human_delta = sqrt(((phi[type == "human"] - phi[type == "hepes"])^2) + 
                              (psi[type == "human"] - psi[type == "hepes"])^2),
         human_exp_delta = sqrt(((phi[type == "human_exp"] - phi[type == "hepes"])^2) + 
                              (psi[type == "human_exp"] - psi[type == "hepes"])^2),
         delta_delta = human_delta - human_exp_delta) %>%
  left_join(apo, by = "residue") %>%
  select(-H_ppm, -N_ppm, -label)
phi_psi_youle

# delta-delta-delta plot of phi-psi angles between human pdb and human CS
phi_psi_youle %>%
  select(residue, aa, delta_delta, secondary) %>%
  distinct() %>%
  drop_na() %>%
  ggplot(aes(x = residue, y = delta_delta, fill = secondary)) +
  geom_col() +
  labs(fill = "Secondary Structure:",
       x = "Residue",
       y = expression(paste("|PHI/PSI"["OUT"], " - PHI/PSI"["Exp"],"| -","|PHI/PSI"["OUT-CS"], " - PHI/PSI"["Exp"],"|"))) + 
  scale_fill_manual(breaks = c("arm", "helix", "loop"),
                    labels = c("N-arm", "Helix", "Loop"),
                    values = c("#1b9e77", "#7570b3", "#d95f02")) +
  theme(legend.position = "top")

ggsave("delta_delta_phi_psi_human_humanSHIFTS_exp.pdf",
       width = 12, height = 8, units = "cm")

# PHI/PSI difference plot for 1PC2 PDB and experimental CS from 1PC2  -----
phi_psi_1pc2_only <- reference %>%
  group_by(residue) %>%
  filter(type != "mouse" & residue > 3) %>%
  mutate(delta = sqrt(((phi[type == "human"] - phi[type == "human_exp"])^2) + 
                              (psi[type == "human"] - psi[type == "human_exp"])^2)) %>%
  left_join(apo, by = "residue") %>%
  select(-H_ppm, -N_ppm, -label)
phi_psi_1pc2_only

# delta-delta plot of phi-psi angles between human pdb and human CS
phi_psi_1pc2_only %>%
  select(residue, aa, delta, secondary) %>%
  distinct() %>%
  drop_na() %>%
  ggplot(aes(x = residue, y = delta, fill = secondary)) +
  geom_col() +
  labs(fill = "Secondary Structure:",
       x = "Residue",
       y = expression(paste("|PHI/PSI"["Human"], " - PHI/PSI"["Human_Exp"],"|"))) + 
  scale_fill_manual(breaks = c("arm", "helix", "loop"),
                    labels = c("Arm", "Helix", "Loop"),
                    values = c("#1b9e77", "#7570b3", "#d95f02")) +
  theme(legend.position = "top")

ggsave("delta_phi_psi_1PC2_torsion_differences.pdf",
       width = 12, height = 8, units = "cm")

# PHI difference plot for 1PC2 PDB and experimental CS from 1PC2  -----
phi_1pc2_only <- reference %>%
  group_by(residue) %>%
  filter(type != "mouse" & residue > 3) %>%
  mutate(delta = (phi[type == "human"] - phi[type == "human_exp"])) %>%
  left_join(apo, by = "residue") %>%
  select(-H_ppm, -N_ppm, -label)
phi_1pc2_only

# delta-delta plot of phi angles between human pdb and human CS
phi_1pc2_only %>%
  select(residue, aa, delta, secondary) %>%
  distinct() %>%
  drop_na() %>%
  ggplot(aes(x = residue, y = delta, fill = secondary)) +
  geom_col() +
  labs(fill = "Secondary Structure:",
       x = "Residue",
       y = expression(paste("|PHI"["Human"], " - PHI"["Human_Exp"],"|"))) + 
  scale_fill_manual(breaks = c("arm", "helix", "loop"),
                    labels = c("Arm", "Helix", "Loop"),
                    values = c("#1b9e77", "#7570b3", "#d95f02")) +
  theme(legend.position = "top")

ggsave("delta_PHI_1PC2_torsion_differences.pdf",
       width = 12, height = 8, units = "cm")

# PSI difference plot for 1PC2 PDB and experimental CS from 1PC2  -----
psi_1pc2_only <- reference %>%
  group_by(residue) %>%
  filter(type != "mouse" & residue > 3) %>%
  mutate(delta = (psi[type == "human"] - psi[type == "human_exp"])) %>%
  left_join(apo, by = "residue") %>%
  select(-H_ppm, -N_ppm, -label)
psi_1pc2_only

# delta-delta plot of psi angles between human pdb and human CS
psi_1pc2_only %>%
  select(residue, aa, delta, secondary) %>%
  distinct() %>%
  drop_na() %>%
  ggplot(aes(x = residue, y = delta, fill = secondary)) +
  geom_col() +
  labs(fill = "Secondary Structure:",
       x = "Residue",
       y = expression(paste("|PSI"["Human"], " - PSI"["Human_Exp"],"|"))) + 
  scale_fill_manual(breaks = c("arm", "helix", "loop"),
                    labels = c("Arm", "Helix", "Loop"),
                    values = c("#1b9e77", "#7570b3", "#d95f02")) +
  theme(legend.position = "top")

ggsave("delta_PHI_1PC2_torsion_differences.pdf",
       width = 12, height = 8, units = "cm")



# Density plots -----------------------------------------------------------

phi_1pc2_only %>%
  select(residue, aa, delta, secondary) %>%
  distinct() %>%
  drop_na() %>%
  ggplot(aes(x = delta, color = secondary)) +
  geom_density()


psi_1pc2_only %>%
  select(residue, aa, delta, secondary) %>%
  distinct() %>%
  drop_na() %>%
  ggplot(aes(x = delta, color = secondary)) +
  geom_density()

phi_psi_1pc2_only %>%
  select(residue, aa, delta, secondary) %>%
  distinct() %>%
  drop_na() %>%
  ggplot(aes(x = delta, color = secondary)) +
  geom_density()

