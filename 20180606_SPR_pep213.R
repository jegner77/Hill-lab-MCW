#
# Data analysis of SPR data collected on 20180606 and 20180607
# CM5 chip with ligands attached by primary amine coupling
# flow cell 1 = reference cell
# flow cell 2 = Fis1
# flow cell 3 = ∆N-Fis1
# flow cell 4 = hMiD51
# titrated pep213 and negative control peptide GPTIEEVD over flow cells
#

# Load libraries, Import data, and tidy data  ----

library(tidyverse)
library(broom)
library(minpack.lm)
library(readxl)
library(gridExtra)

fis1_tr1 <- read_tsv("20180606_pep213tr1_negative_peptide_tr2_tr3_hFis1_results_RU_curves.txt") %>%
  select(., 1, ends_with("Y")) %>%
  gather(., "tmp1", "RU", 2:19) %>%
  select(., 1, 13, 14) %>%
  separate(., col = "tmp1", into = c("conc", "tr"), sep = " ") %>%
  rename(., time_sec = `0 tr1_X`) %>%
  mutate(., tr = parse_number(tr),
         conc = parse_number(conc),
         protein = "Fis1")
  
fis1_data <- read_tsv("20180607_pep213_titration_tr2_tr3_results_hFis1_RU_curves.txt") %>%
  select(., 1, ends_with("Y")) %>%
  gather(., "tmp1", "RU", 2:31) %>%
  separate(., col = "tmp1", into = c("conc", "tr"), sep = " ") %>%
  rename(., time_sec = `0 tr1_X`) %>%
  mutate(., tr = parse_number(tr),
         conc = parse_number(conc),
         protein = "Fis1") %>%
  union(., fis1_tr1) %>%
  filter(., time_sec <= 450)

dNfis1_tr1 <- read_tsv("20180606_pep213tr1_negative_peptide_tr2_tr3_dNhFis1_results_RU_curves.txt") %>%
  select(., 1, ends_with("Y")) %>%
  gather(., "tmp1", "RU", 2:19) %>%
  select(., 1, 13, 14) %>%
  separate(., col = "tmp1", into = c("conc", "tr"), sep = " ") %>%
  rename(., time_sec = `0 tr1_X`) %>%
  mutate(., tr = parse_number(tr),
         conc = parse_number(conc),
         protein = "dN-Fis1")

dNfis1_data <- read_tsv("20180607_pep213_titration_tr2_tr3_results_dN_hFis1_RU_curves.txt") %>%
  select(., 1, ends_with("Y")) %>%
  gather(., "tmp1", "RU", 2:31) %>%
  separate(., col = "tmp1", into = c("conc", "tr"), sep = " ") %>%
  rename(., time_sec = `0 tr1_X`) %>%
  mutate(., tr = parse_number(tr),
         conc = parse_number(conc),
         protein = "dN-Fis1") %>%
  union(., dNfis1_tr1) %>%
  filter(., time_sec <= 450)

mid51_tr1 <- read_tsv("20180606_pep213tr1_negative_peptide_tr2_tr3_hMiD51_results_RU_curves.txt") %>%
  select(., 1, ends_with("Y")) %>%
  gather(., "tmp1", "RU", 2:19) %>%
  select(., 1, 13, 14) %>%
  filter(., RU != "NA") %>%
  separate(., col = "tmp1", into = c("conc", "tr"), sep = " ") %>%
  rename(., time_sec = `0 tr1_X`) %>%
  mutate(., tr = parse_number(tr),
         conc = parse_number(conc),
         protein = "MiD51")

mid51_data <- read_tsv("20180607_pep213_titration_tr2_tr3_results_hMiD51_RU_curves.txt") %>%
  select(., 1, ends_with("Y")) %>%
  gather(., "tmp1", "RU", 2:31) %>%
  separate(., col = "tmp1", into = c("conc", "tr"), sep = " ") %>%
  rename(., time_sec = `0 tr1_X`) %>%
  mutate(., tr = parse_number(tr),
         conc = parse_number(conc),
         protein = "MiD51") %>%
  union(., mid51_tr1) %>%
  filter(., time_sec <= 450)

data <- fis1_data %>%
  union(., dNfis1_data) %>%
  union(., mid51_data) %>%
  filter(., tr == 1 | tr == 2 | tr == 3)

TRENDdata_Fis1 <- data %>%
  filter(., protein == "Fis1" & tr == 1) %>%
  split(.$conc)

TRENDdata_Fis1 %>%
  map(~distinct(., time_sec, tr, RU, .keep_all=TRUE)) %>% 
  walk(~.x %>%
         write.csv(file = paste0(unique(.x$conc),"_Fis1.csv"),
                   row.names = FALSE))

# Visualize Response curves and fit curves for each protein and tr  ----

data %>%
  filter(., time_sec != 184.5  & time_sec != 185.5  & time_sec != 277.5 & 
           time_sec != 278.5 & time_sec != 279.5 & time_sec != 280.5) %>%
  ggplot(., aes(x = time_sec, y = RU, color = as.factor(conc))) +
  geom_line() +
  facet_grid(tr ~ protein) +
  scale_x_continuous(limits = c(0, 450),
                     breaks = c(0, 100, 200, 300, 400)) +
  scale_y_continuous(limits = c(-5, 300),
                     breaks = c(0, 100, 200, 300)) +
  labs(title = "SPR Response Curves: Facet by Protein and TR",
       color = "[pep213] (µM)",
       x = "Time (sec)",
       y = "Response (RU)") +
  scale_color_hue(h = c(180, 225)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid')
  )
ggsave("Response_curves_facet_by_protein_tr.pdf", 
       width = 20, height = 16, units = "cm")

# TR1, dN-Fis1 only
data %>%
  filter(., time_sec != 184.5  & time_sec != 185.5  & time_sec != 277.5 & 
           time_sec != 278.5 & time_sec != 279.5 & time_sec != 280.5 &
           protein == "dN-Fis1" & tr == 1) %>%
  mutate(., time_sec = time_sec - 160) %>%
  ggplot(., aes(x = time_sec, y = RU, color = as.factor(conc))) +
  geom_line() +
  scale_x_continuous(limits = c(0, 200),
                     breaks = c(0, 50, 100, 150, 200)) +
  scale_y_continuous(limits = c(-5, 300),
                     breaks = c(0, 100, 200, 300)) +
  labs(title = "SPR Response Curve: dN-Fis1",
       color = "[pep213] (µM)",
       x = "Time (sec)",
       y = "Response (RU)") +
  scale_color_hue(h = c(0, 360)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid')
  )
ggsave("Response_curves_TR1_dNFis1.pdf", 
       width = 12, height = 10, units = "cm")

# TR1, Fis1 only
data %>%
  filter(., time_sec != 184.5  & time_sec != 185.5  & time_sec != 277.5 & 
           time_sec != 278.5 & time_sec != 279.5 & time_sec != 280.5 &
           protein == "Fis1" & tr == 1) %>%
  mutate(., time_sec = time_sec - 160) %>%
  ggplot(., aes(x = time_sec, y = RU, color = as.factor(conc))) +
  geom_line() +
  scale_x_continuous(limits = c(0, 200),
                     breaks = c(0, 50, 100, 150, 200)) +
  scale_y_continuous(limits = c(-5, 300),
                     breaks = c(0, 100, 200, 300)) +
  labs(title = "SPR Response Curve: Fis1",
       color = "[pep213] (µM)",
       x = "Time (sec)",
       y = "Response (RU)") +
  scale_color_hue(h = c(0, 360)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid')
  )
ggsave("Response_curves_TR1_Fis1.pdf", 
       width = 12, height = 10, units = "cm")

# TR1, MiD51 only
data %>%
  filter(., time_sec != 184.5  & time_sec != 185.5  & time_sec != 277.5 & 
           time_sec != 278.5 & time_sec != 279.5 & time_sec != 280.5 &
           protein == "MiD51" & tr == 1 & time_sec > 180) %>%
  mutate(., time_sec = time_sec - 160) %>%
  ggplot(., aes(x = time_sec, y = RU, color = as.factor(conc))) +
  geom_line() +
  scale_x_continuous(limits = c(0, 200),
                     breaks = c(0, 50, 100, 150, 200)) +
  scale_y_continuous(limits = c(-5, 300),
                     breaks = c(0, 100, 200, 300)) +
  labs(title = "SPR Response Curve: MiD51",
       color = "[pep213] (µM)",
       x = "Time (sec)",
       y = "Response (RU)") +
  scale_color_hue(h = c(0, 360)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid')
  )
ggsave("Response_curves_TR1_MiD51.pdf", 
       width = 12, height = 10, units = "cm")

# Kd determination of pep213 by fitting at equilibrium ----

# Calculate average value at equilibrium of association phase
equil_tr <- data %>%
  filter(., time_sec >= 240 & time_sec <= 260) %>%
  group_by(., protein, conc, tr) %>%
  summarise(., Req = mean(RU),
            stdev_Req = sd(RU))

equil <- data %>%
  filter(., time_sec >= 240 & time_sec <= 260) %>%
  group_by(., protein, conc) %>%
  summarise(., Req = mean(RU),
            stdev_Req = sd(RU))

# Plot Req vs concentration at equilibrium with Kd fitting by tr and protein
equil_tr %>%
  ggplot(., aes(x = conc, y = Req, color = as.factor(protein))) +
  geom_point(shape = 15) +
  geom_errorbar(aes(ymin = Req - stdev_Req, ymax = Req + stdev_Req), width = 50) +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ (rmax*x) / (Kd + x),
            method.args = list(start = c(rmax = 200,
                                         Kd = 10),
                               control = nls.control(maxiter = 100, tol = 1e-6)),
            se = FALSE,
            fullrange = TRUE,
            size = 1,
            alpha = 0.6) +
  facet_grid(tr ~ protein) +
  scale_x_continuous(limits = c(0, 1200),
                     breaks = c(0, 200, 400, 600, 800, 1000)) +
  scale_y_continuous(limits = c(-5, 280),
                     breaks = c(0, 50, 100, 150, 200, 250)) +
  labs(title = "Kd fitting: Req vs. [pep213]",
       color = "Protein",
       x = "[pep213] (µM)",
       y = "Req (RU)") +
  scale_color_brewer(palette = "Set2") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid')
  )

ggsave("Req_vs_pep213_Kd_fits_facet_by_protein_tr.pdf", 
       width = 20, height = 16, units = "cm")
  
# Kd fitting from equilibria at association phase
equil_Kd <- equil_tr %>%
  group_by(., protein) %>%
  do(tidy(nlsLM(
      formula = Req ~ (rmax*conc) / (Kd + conc),
      start = list(rmax = 200, Kd = 10),
      trace = TRUE,
      control = nls.control(maxiter = 100, tol = 1e-6),
      data = .))) %>%
      spread(., key = term, value = estimate)

write_csv(equil_Kd, "pep213_Kd_fits_by_protein_tr.csv")

equil_Kd

# Residuals of Kd fitting from equilibria at association phase
equil_Kd_residuals <- equil_tr %>%
  group_by(., protein) %>%
  do(augment(nlsLM(
    formula = Req ~ (rmax*conc) / (Kd + conc),
    start = list(rmax = 200, Kd = 10),
    trace = TRUE,
    control = nls.control(maxiter = 100, tol = 1e-6),
    data = .)))

write_csv(equil_Kd_residuals, "pep213_Kd_fits_by_protein_tr_residuals.csv")

# Plot Req vs concentration at equilibrium with Kd fitting
a <- equil %>%
  ggplot(., aes(x = conc, y = Req, color = as.factor(protein))) +
  geom_point(size = 1, shape = 15) +
  geom_errorbar(aes(ymin = Req - stdev_Req, ymax = Req + stdev_Req), width = 50) +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ (rmax*x) / (kd + x),
            method.args = list(start = c(rmax = 200,
                                         kd = 10),
                               control = nls.control(maxiter = 100, tol = 1e-6)),
            se = FALSE,
            fullrange = F,
            size = 1,
            alpha = 0.6) +
  scale_x_continuous(limits = c(-50, 1050),
                     breaks = c(0, 200, 400, 600, 800, 1000)) +
  scale_y_continuous(limits = c(-5, 280),
                     breaks = c(0, 50, 100, 150, 200, 250)) +
  labs(title = "Kd fitting: Req vs. [pep213]",
       color = "Protein",
       x = "[pep213] (µM)",
       y = "Req (RU)") +
  scale_color_brewer(palette = "Set2") +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.length = unit(0.2, "cm")
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid')
  )
a

# Plot Req vs concentration at equilibrium with Kd fitting
b <- equil_Kd_residuals %>%
  ggplot(., aes(x = conc, y = .resid, color = as.factor(protein))) +
  geom_point(size = 1, shape = 15) +
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  scale_x_continuous(limits = c(-50, 1050),
                     breaks = c(0, 200, 400, 600, 800, 1000)) +
  scale_y_continuous(limits = c(-45, 45),
                     breaks = c(-30, 0, 30)) +
  labs(title = "Kd fitting residuals",
       color = "Protein",
       x = "",
       y = "Residuals") +
  scale_color_brewer(palette = "Set2") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text.x = element_blank()
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid')
  )
b

# Plotting regression with residuals, model 2 ----
lay <- cbind(c(1, 2, 2))

req_plot <- grid.arrange(grobs = list(b, a), layout_matrix = lay,
                      top = "Model 2: Req KD fitting")

ggsave("Req_SPR_KD_fitting_pep213.pdf", req_plot,
       width = 10, height = 10, units = "cm")
ggsave("Req_SPR_KD_fitting_pep213.png", req_plot,
       width = 10, height = 10, units = "cm")
 
# Log (log10) transformation of concentration and plotting/fitting sigmoid ----
 
equil %>%
  filter(., conc != 0) %>%
  mutate(., log_conc = log10(conc)) %>%
  ggplot(., aes(x = log_conc, y = Req, color = protein)) +
  geom_point() +
  geom_errorbar(aes(ymin = Req - stdev_Req, ymax = Req + stdev_Req), width = 0.2) +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ a + ((b - a) / (1 + exp((Kd - x)/c))),
            method.args = list(start = c(a = 0, b = 200, c = 10, Kd = 30),
                              control = nls.control(maxiter = 200, tol = 1e-6)),
            se = F,
            fullrange = F,
            size = 0.5,
            alpha = 0.6) +
  scale_x_continuous(limits = c(-0.8, 3.2),
                     breaks = c(0, 1, 2, 3)) +
  scale_y_continuous(limits = c(-10, 300),
                     breaks = c(0, 100, 200, 300)) +
  labs(title = "Kd fitting: Log10-Scale vs Response, n = 3, Mean +/- SD",
      color = "Protein",
      x = "log10[pep213] (µM)",
      y = "Req (RU)") +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = "black"),
       axis.title.x = element_text(size = 12),
       axis.title.y = element_text(size = 12),
       panel.grid.major = element_blank(), 
       panel.grid.minor = element_blank()
       # panel.border = element_blank()
       # panel.background = element_blank(),
       # axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
       # axis.line.y = element_line(colour = 'black', size=1, linetype='solid')
  )
 
 ggsave("Req_vs_pep213_sigmoid_fit.pdf", 
        width = 16, height = 12, units = "cm")

# Kd determination by sigmoidal fit
sigmoid_fits_n3 <- equil %>%
  filter(., conc != 0) %>%
  mutate(., log_conc = log10(conc)) %>%
  group_by(., protein) %>%
  do(tidy(nlsLM(Req ~ a + ((b - a) / (1 + exp((Kd - log_conc)/c))),
                start = list(a = 0, b = 200, c = 10, Kd = 30), trace = TRUE, 
                control = nls.control(maxiter = 200, tol = 1e-6), .))) %>%
  spread(., key = term, value = estimate) %>%
  mutate(., Kd_correct = 10^(Kd))
 
 write_csv(sigmoid_fits_n3, "pep213_Kd_sigmoid_fits_by_protein.csv")
 
# Residual calculation and plotting of Kd determination by sigmoidal fit
 sigmoid_fits_n3_residuals <- equil %>%
   filter(., conc != 0) %>%
   mutate(., log_conc = log10(conc)) %>%
   group_by(., protein) %>%
   do(augment(nlsLM(Req ~ a + ((b - a) / (1 + exp((Kd - log_conc)/c))),
                 start = list(a = 0, b = 200, c = 10, Kd = 30), trace = TRUE, 
                 control = nls.control(maxiter = 200, tol = 1e-6), .)))
 
 write_csv(sigmoid_fits_n3_residuals,
           "pep213_Kd_sigmoid_fits_by_protein_residuals.csv")
 
sigmoid_fits_n3_residuals %>%
  ggplot(., aes(x = log_conc, y = .resid, color = protein)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linetype = 2) +
  scale_x_continuous(limits = c(-0.8, 3.2),
                    breaks = c(0, 1, 2, 3)) +
  scale_y_continuous(limits = c(-20, 20),
                    breaks = c(-15, 0, 15)) +
  labs(title = "Residuals of Kd fitting: Log10-Scale vs Response, n = 3, Mean +/- SD",
      color = "Protein",
      x = "log10[pep213] (µM)",
      y = "Residuals") +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = "black"),
       axis.title.x = element_text(size = 12),
       axis.title.y = element_text(size = 12),
       panel.grid.major = element_blank(), 
       panel.grid.minor = element_blank()
       # panel.border = element_blank()
       # panel.background = element_blank(),
       # axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
       # axis.line.y = element_line(colour = 'black', size=1, linetype='solid')
  )
 
 ggsave("Req_vs_pep213_sigmoid_fit_residuals.pdf", 
        width = 16, height = 6, units = "cm")
 
# Kd determination of pep213 by association and dissociation rate fitting ----

# Kd fitting by using association equation
association_fits <- data %>%
  filter(., time_sec != 184.5  & time_sec != 185.5  & time_sec != 277.5 & 
           time_sec != 278.5 & time_sec != 279.5 & time_sec != 280.5 &
           time_sec >= 180 & time_sec <= 260 & conc != 0 & time_sec > 180) %>%
  mutate(., time_sec = time_sec - 180) %>%
  group_by(., protein, tr) %>%
  do(tidy(nlsLM(
    formula = RU ~ (((Rmax * conc) / ((kd / ka) + conc)) * (1 - exp(-(ka * conc + kd) * time_sec)) + Ro),
    start = list(Rmax = 200,
                 kd = 1,
                 ka = 1e-4,
                 Ro = -5),
    trace = TRUE,
    control = nls.control(maxiter = 100, tol = 1e-6),
    data = .))) %>%
  spread(., key = term, value = estimate)

write_csv(association_fits, "pep213_Kd_fits_by_association_rate_equation.csv")

association_augment <- data %>%
  filter(., time_sec != 184.5  & time_sec != 185.5  & time_sec != 277.5 & 
           time_sec != 278.5 & time_sec != 279.5 & time_sec != 280.5 &
           time_sec >= 180 & time_sec <= 260 & conc != 0 & tr == 1 & protein == "Fis1" &
           time_sec > 182) %>%
  mutate(., time_sec = time_sec - 180) %>%
  group_by(., protein, tr) %>%
  do(augment(nlsLM(
    formula = RU ~ (((Rmax * conc) / ((kd / ka) + conc)) * (1 - exp(-(ka * conc + kd) * time_sec)) + Ro),
    start = list(Rmax = 10,
                 kd = 0.02,
                 ka = 5e-5,
                 Ro = 1),
    trace = TRUE,
    control = nls.control(maxiter = 100, tol = 1e-6),
    data = .)))

data %>%
  filter(., time_sec != 184.5  & time_sec != 185.5  & time_sec != 277.5 & 
           time_sec != 278.5 & time_sec != 279.5 & time_sec != 280.5 &
           time_sec >= 180 & time_sec <= 260 & conc != 0 & tr == 1 & protein == "Fis1" &
           time_sec > 182) %>%
  mutate(., time_sec = time_sec - 180) %>%
  ggplot(., aes(x = time_sec, y = RU, color = as.factor(conc))) +
  geom_jitter(size = 1) +
  geom_line(aes(x = time_sec, y = association_augment$.fitted)) +
  facet_grid(tr ~ protein) +
  theme_bw()

# dissociation fitting
data %>%
  filter(., time_sec != 184.5  & time_sec != 185.5  & time_sec != 277.5 & 
           time_sec != 278.5 & time_sec != 279.5 & time_sec != 280.5 &
           time_sec >= 275 & time_sec <= 400 & conc > 0) %>%
  ggplot(., aes(x = time_sec, y = RU, color = as.factor(conc))) +
  geom_point(size = 0.2) +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ Ro*exp(-kd * x),
            method.args = list(start = c(Ro = 200,
                                         kd = 0.01),
                               control = nls.control(maxiter = 300, tol = 1e-6)),
            se = FALSE,
            fullrange = TRUE,
            size = 1,
            alpha = 0.6) +
  facet_grid(tr ~ protein) +
  scale_x_continuous(limits = c(275, 400),
                     breaks = c(280, 300, 320, 340, 360, 380, 400)) +
  scale_y_continuous(limits = c(-5, 300),
                     breaks = c(0, 100, 200, 300)) +
  labs(title = "SPR Response Curves: Facet by Protein and TR",
       color = "[pep213] (µM)",
       x = "Time (sec)",
       y = "Response (RU)") +
  scale_color_hue(h = c(180, 225)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid')
  )
ggsave("Dis_rate_fitting_facet_by_protein_tr.pdf",
       width = 20, height = 16, units = "cm")

# # Kd fitting by using dissociation equation
# diss_Kd <- data %>%
#   filter(., time_sec != 184.5  & time_sec != 185.5  & time_sec != 277.5 &
#            time_sec != 278.5 & time_sec != 279.5 & time_sec != 280.5 &
#            time_sec >= 180 & time_sec <= 260 & conc > 0) %>%
#   group_by(., protein, conc) %>%
#   do(augment(nlsLM(
#     formula = RU ~ Ro * exp(-kd * time_sec),
#     start = list(Ro = 200, kd = 0.01),
#     trace = TRUE,
#     control = nls.control(maxiter = 100, tol = 1e-6),
#     data = .)))
# 
# write_csv(diss_Kd, "pep213_Kd_fits_by_dissociation_rate_equation.csv")
# 
# diss_Kd %>%
#   ggplot(., aes(x = time_sec, y = .fitted, color = as.factor(conc))) +
#   geom_point() +
#   facet_grid(~ protein)
