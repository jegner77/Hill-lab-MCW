# 
# Data analysis of pep213 HSQC titration against 15N-hFis1
# For Kd determination, subset of residues analyzed by CSPs and
# due to intermediate exchange,
# subset of residues also analyzed by change in Intensity
#

# Load libraries, import data, and tidy data    --------------------------------

library(tidyverse)
library(broom)
library(minpack.lm)
library(readxl)

# Intensity data with true apo peaks selected
int_files <- list.files(path = "./intensities/")

setwd(dir = "./intensities/") #set path to csv files with shifts only

int <- int_files %>%
  map(., function(x) read_delim(x, delim = ";")) %>%
  reduce(., left_join) %>%
  select(., -X5)

setwd(dir = "../.") # set path back to working directory 
getwd()

tidy_int <- int %>%
  gather(., "tmp", "tmp2", 3:18) %>%
  separate(., col = tmp, into = c("tmp3", "conc", "tmp5", "tmp6")) %>%
  select(., -tmp3, -tmp5, -`peak no`) %>%
  separate(., col = label, into = c("nuclei", "label"), sep = " ") %>%
  spread(., key = tmp6, value = tmp2) %>%
  mutate(., residue = parse_number(label)) %>%
  arrange(., residue) %>%
  mutate(., conc = parse_number(conc),
         Amp = parse_number(Amp),
         Vol = parse_number(Vol)) %>%
  select(., nuclei, residue, label, conc, Amp, Vol)

# # Intensity data with bound position selected for subset of apo peaks
# slowint_files <- list.files(path = "./pep213_intensity_fitting/")
# 
# setwd(dir = "./pep213_intensity_fitting/") #set path to csv files with shifts only
# 
# slowint <- slowint_files %>%
#   map(., function(x) read_delim(x, delim = ";")) %>%
#   reduce(., left_join) %>%
#   select(., -X5)
# 
# setwd(dir = "../.") # set path back to working directory 
# getwd()
# 
# slowtidy_int <- slowint %>%
#   gather(., "tmp", "tmp2", 3:18) %>%
#   separate(., col = tmp, into = c("tmp3", "conc", "tmp5", "tmp6")) %>%
#   select(., -tmp3, -tmp5, -`peak no`) %>%
#   separate(., col = label, into = c("nuclei", "label"), sep = " ") %>%
#   spread(., key = tmp6, value = tmp2) %>%
#   mutate(., residue = parse_number(label)) %>%
#   arrange(., residue) %>%
#   mutate(., conc = parse_number(conc),
#          Amp = parse_number(Amp),
#          Vol = parse_number(Vol)) %>%
#   select(., nuclei, residue, label, conc, Amp, Vol) %>%
#   filter(., residue == 13 | residue == 36 | residue == 40 | residue == 49 | 
#            residue == 57 | residue == 75 | residue == 76 | residue == 82 | 
#            residue == 84 | residue == 117 | residue == 119 | residue == 121)

# load-in manually edited excel file for Kd fitting of fast exchange residues

raw_cs <- read_excel("20180525_hFis1_pep213_shifts.xlsx",
                     sheet = 1) %>%
  mutate(., residue = parse_number(label))

# All Residues: Visualize intensity (Vol) changes   ----------------------

tidy_int %>%
  filter(., conc > 0) %>%
  ggplot(., aes(x = conc, y = Vol, color = factor(residue), shape = factor(nuclei))) +
  geom_point() +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ dmax*((kd+50+x) - ((kd+50+x)^2 - 4*50*x)^0.5)/(2*50),
            method.args = list(start = c(dmax = 1000,
                                         kd = 10),
                               control = nls.control(maxiter = 100, tol = 1e-6)),
            se = FALSE,        
            fullrange = TRUE, 
            size = 1,
            alpha = 0.6) +        # p held constant at 50 (µM) in model
  facet_wrap(~ residue) +
  labs(x = "[pep213] (µM)",
       y = "Volume (Crosspeak Intensity)") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("hFis1_pep213_Vol_by_res_all.pdf",
       width = 40, height = 40, units = "cm")

# Select testing w/ Model 1: Kd Model fitting from intensity (Amp) changes  ----

tidy_int %>%
  filter(., nuclei == "H/N" | nuclei == "HE1/NE1") %>%
  ggplot(., aes(x = conc, y = Amp, color = factor(residue),
                shape = factor(nuclei))) +
  geom_point() +
  facet_wrap(~ residue) +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ dmax*((kd+50+x) - ((kd+50+x)^2 - 4*50*x)^0.5)/(2*50),
            method.args = list(start = c(dmax = 2000,
                                         kd = 10),
                               control = nls.control(maxiter = 100, tol = 1e-6)),
            se = FALSE,        
            fullrange = TRUE,  
            size = 1,
            alpha = 0.6) +        # p held constant at 50 (µM) in model
  scale_x_continuous(limits = c(-10, 2000),
                     breaks = c(0, 500, 1000, 1500, 2000)) +
  labs(color = "Residue",
       shape = "Nuclei",
       x = "[pep213] (µM)",
       y = "Crosspeak Intensity") +
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
  
ggsave("hFis1_pep213_Amp_by_res.pdf",
       width = 40, height = 40, units = "cm")

# Kd determination after excluding residues 13, 75, 119, and 121 (poor fits)
tidy_int %>%
  filter(., residue == 3 & conc > 0 | residue == 40 & nuclei == "HE1/NE1"
         & conc > 0 | residue == 51 & conc > 0 | residue == 57 & conc > 0) %>%
  ggplot(., aes(x = conc, y = Amp, color = factor(residue),
                shape = factor(nuclei))) +
  geom_point() +
  facet_wrap(~ residue) +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ dmax*((kd+50+x) - ((kd+50+x)^2 - 4*50*x)^0.5)/(2*50),
            method.args = list(start = c(dmax = 2000,
                                         kd = 10),
                              control = nls.control(maxiter = 100, tol = 1e-6)),
           se = FALSE,
           fullrange = FALSE,
           size = 1,
           alpha = 0.6) +        # p held constant at 50 (µM) in model
 scale_x_continuous(limits = c(-10, 2000),
                    breaks = c(0, 500, 1000, 1500, 2000)) +
 scale_y_continuous(limits = c(-100, 5000),
                      breaks = c(0, 1000, 2000, 3000, 4000, 5000)) +
 labs(color = "Residue",
      shape = "Nuclei",
      x = "[pep213] (µM)",
      y = "Crosspeak Intensity") +
 theme_bw() +
 theme(axis.text = element_text(size = 10, color = "black"),
       axis.title.x = element_text(size = 12),
       axis.title.y = element_text(size = 12),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       legend.position = "none"
       #panel.border = element_blank(),
       #panel.background = element_blank(),
       #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
       #axis.line.y = element_line(colour = 'black', size=1, linetype='solid')
       )

ggsave("hFis1_pep213_Amp_by_res_select_good_fits_facet.pdf",
        width = 24, height = 14, units = "cm")
 
# Kd determination by gloabl fitting select residue xpk intensities
tidy_int %>%
  filter(., residue == 3 & conc > 0 | residue == 40 & nuclei == "HE1/NE1"
         & conc > 0 | residue == 51 & conc > 0 | residue == 57 & conc > 0) %>%
 ggplot(., aes(x = conc, y = Amp, color = factor(residue),
               shape = factor(nuclei))) +
 geom_point() +
 geom_line(stat = "smooth", method = "nlsLM",
           formula = y ~ dmax*((kd+50+x) - ((kd+50+x)^2 - 4*50*x)^0.5)/(2*50),
           method.args = list(start = c(dmax = 2000,
                                        kd = 10),
                              control = nls.control(maxiter = 100, tol = 1e-6)),
           se = FALSE,        
           fullrange = FALSE,  
           size = 1,
           alpha = 0.6) +        # p held constant at 50 (µM) in model
 #scale_x_continuous(limits = c(-10, 2000),
  #                  breaks = c(0, 500, 1000, 1500, 2000)) +
 #scale_y_continuous(limits = c(-100, 5000),
  #                  breaks = c(0, 1000, 2000, 3000, 4000, 5000)) +
 scale_color_brewer(name = "Residue", palette = "Dark2") +
 labs(color = "Residue",
      shape = "Nuclei",
      x = "[pep213] (µM)",
      y = "Crosspeak Intensity") +
 theme_bw() +
 theme(axis.text = element_text(size = 10, color = "black"),
       axis.title.x = element_text(size = 12),
       axis.title.y = element_text(size = 12),
       panel.grid.major = element_blank(), 
       panel.grid.minor = element_blank(),
       panel.border = element_blank(),
       panel.background = element_blank(),
       axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
       axis.line.y = element_line(colour = 'black', size=1, linetype='solid')
 )

ggsave("hFis1_pep213_Amp_by_res_select_good_fits.pdf",
      width = 24, height = 14, units = "cm")

# Calculate Kd value from subset of data
select_kd_int <- tidy_int %>% 
  filter(., residue == 3 & conc > 0 | residue == 40 & nuclei == "HE1/NE1"
         & conc > 0 | residue == 51 & conc > 0 | residue == 57 & conc > 0) %>%
 group_by(., residue, nuclei) %>%
 do(tidy(nlsLM(
   formula = Vol ~ dmax*((kd+50+conc) - ((kd+50+conc)^2 - 4*50*conc)^0.5)/(2*50),
   start = list(dmax = 2000, kd = 10), 
   trace = TRUE,
   control = nls.control(maxiter = 100, tol = 1e-6),
   data = .))) %>%
 spread(., key = term, value = estimate)

select_kd_int

write_csv(select_kd_int, "hFis1_pep213_Kd_select_kd_by_intensities.csv")

# calculate residuals with augment() and plot, for indicated Residues
model_int_select_residuals <- tidy_int %>%
  filter(., residue == 3 & conc > 0 | residue == 40 & nuclei == "HE1/NE1"
         & conc > 0 | residue == 51 & conc > 0 | residue == 57 & conc > 0) %>%
 group_by(., residue, nuclei) %>%
 do(augment(nlsLM(
   formula = Vol ~ dmax*((kd+50+conc) - ((kd+50+conc)^2 - 4*50*conc)^0.5)/(2*50),
   start = list(dmax = 2000, kd = 10), 
   trace = TRUE,
   control = nls.control(maxiter = 100, tol = 1e-6),
   data = .)))

write_csv(model_int_select_residuals,
         "hFis1_pep213_Kd_select_kd_by_intensities_residuals.csv")

model_int_select_residuals %>%
  filter(., residue == 3 & conc > 0 | residue == 40 & nuclei == "HE1/NE1"
         & conc > 0 | residue == 51 & conc > 0 | residue == 57 & conc > 0) %>%
  ggplot(., aes(x = conc, y = .resid, color = factor(residue),
                shape = factor(nuclei))) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linetype = 2) +
  facet_wrap(~ residue) +
  scale_x_continuous(limits = c(-10, 2000),
                     breaks = c(0, 400, 800, 1200, 1600, 2000)) +
  scale_y_continuous(limits = c(-800, 800),
                     breaks = c(-800, -400, 0, 400, 800)) +
  labs(shape = "Nuclei",
       x = "[pep213] (µM)", y = "Residual") +
 theme_bw()
ggsave("hFis1_pep213_Vol_by_res_select_facet_residuals.pdf",
       width = 24, height = 14, units = "cm")

# Select testing w/ Model 2: Kd Model fitting from intensity (Amp) changes  ----

tidy_int %>%
  filter(., residue == 3 & conc > 0 | residue == 40 & nuclei == "HE1/NE1"
         & conc > 0 | residue == 51 & conc > 0 | residue == 57 & conc > 0) %>%
  ggplot(., aes(x = conc, y = Amp, color = factor(residue),
                shape = factor(nuclei))) +
  geom_point() +
  facet_wrap(~ residue) +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ (dmax*x) / (kd + x),
            method.args = list(start = c(dmax = 2000,
                                         kd = 10),
                               control = nls.control(maxiter = 100, tol = 1e-6)),
            se = FALSE,
            fullrange = F,
            size = 1,
            alpha = 0.6) +
  scale_x_continuous(limits = c(-10, 2000),
                     breaks = c(0, 500, 1000, 1500, 2000)) +
  labs(title = "Model 2",
       color = "Residue",
       shape = "Nuclei",
       x = "[pep213] (µM)",
       y = "Crosspeak Intensity") +
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

ggsave("hFis1_pep213_Amp_by_res_select_model2.pdf",
       width = 24, height = 16, units = "cm")

# Kd determination by global fit of select residue xpk intensities
tidy_int %>%
  filter(., residue == 3 & conc > 0 | residue == 40 & nuclei == "HE1/NE1"
         & conc > 0 | residue == 51 & conc > 0 | residue == 57 & conc > 0) %>%
  ggplot(., aes(x = conc, y = Amp, color = factor(residue),
                shape = factor(nuclei))) +
  geom_point() +
  facet_wrap(~ residue) +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ (dmax*x) / (kd + x),
            method.args = list(start = c(dmax = 2000,
                                         kd = 10),
                               control = nls.control(maxiter = 100, tol = 1e-6)),
            se = FALSE,
            fullrange = F,
            size = 1,
            alpha = 0.6) +
  scale_x_continuous(limits = c(-10, 2000),
                     breaks = c(0, 500, 1000, 1500, 2000)) +
  scale_y_continuous(limits = c(-100, 5000),
                     breaks = c(0, 1000, 2000, 3000, 4000, 5000)) +
  scale_color_brewer(name = "Residue", palette = "Dark2") +
  labs(title = "Model 2",
       color = "Residue",
       shape = "Nuclei",
       x = "[pep213] (µM)",
       y = "Crosspeak Intensity") +
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

ggsave("hFis1_pep213_Amp_by_res_select_good_fits_facet_model2.pdf",
      width = 24, height = 14, units = "cm")

# Calculate Kd value from subset of data
select_kd_int_model2 <- tidy_int %>%
  filter(., residue == 3 & conc > 0 | residue == 40 & nuclei == "HE1/NE1"
         & conc > 0 | residue == 51 & conc > 0 | residue == 57 & conc > 0) %>%
  group_by(., residue, nuclei) %>%
  do(tidy(nlsLM(
    formula = Amp ~ (dmax*conc) / (kd + conc),
    start = list(dmax = 2000, kd = 10),
    trace = TRUE,
    control = nls.control(maxiter = 100, tol = 1e-6),
    data = .))) %>%
  spread(., key = term, value = estimate)

select_kd_int_model2

write_csv(select_kd_int_model2,
          "hFis1_pep213_Kd_select_kd_by_intensities_model2.csv")

# calculate residuals with augment() and plot, for indicated Residues
model_int_select_residuals_model2 <- tidy_int %>%
  filter(., residue == 3 & conc > 0 | residue == 40 & nuclei == "HE1/NE1"
         & conc > 0 | residue == 51 & conc > 0 | residue == 57 & conc > 0) %>%
  group_by(., residue, nuclei) %>%
  do(augment(nlsLM(
    formula = Amp ~ (dmax*conc) / (kd + conc),
    start = list(dmax = 2000, kd = 10),
    trace = TRUE,
    control = nls.control(maxiter = 100, tol = 1e-6),
    data = .)))

write_csv(model_int_select_residuals_model2,
          "hFis1_pep213_Kd_select_kd_by_intensities_residuals_model2.csv")

model_int_select_residuals_model2 %>%
  filter(., residue == 3 & conc > 0 | residue == 40 & nuclei == "HE1/NE1"
         & conc > 0 | residue == 51 & conc > 0 | residue == 57 & conc > 0) %>%
  ggplot(., aes(x = conc, y = .resid, color = factor(residue),
                shape = factor(nuclei))) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linetype = 2) +
  facet_wrap(~ residue) +
  scale_color_brewer(name = "Residue", palette = "Dark2") +
  scale_x_continuous(limits = c(-10, 2000),
                     breaks = c(0, 400, 800, 1200, 1600, 2000)) +
  scale_y_continuous(limits = c(-800, 1200),
                     breaks = c(-800, -400, 0, 400, 800, 1200)) +
  labs(title = "Model 2",
       shape = "Nuclei",
       x = "[pep213] (µM)", y = "Residual") +
  theme_bw()

ggsave("hFis1_pep213_Vol_by_res_select_facet_residuals_model2.pdf",
       width = 24, height = 14, units = "cm")

# Calculate total chemical shift perturbation and plot against residue ---------

total_shift <- raw_cs %>%
  mutate(., CSP = sqrt(((5*(H_0 - H_2000))^2) + (N_0 - N_2000)^2))

mean_sd <- total_shift %>%
  summarise(., mean = mean(CSP),
            sd = sd(CSP)) %>%
  mutate(twoSD = 2*sd, sigma = mean + sd, twosigma = sigma*2)

mean_sd

total_shift <- total_shift %>%
  mutate(., STDEV = 0.506, TwoSTDEV = 1.01)

write_csv(total_shift, "hFis1_pep213_hsqc_CSPs.csv")

# CSP by Residue # for apo and 2 mM pep213
total_shift %>% 
  filter(., residue != "W40") %>%
  ggplot(., aes(x = residue, y = CSP)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = c(0.506, 1.01), color = c("green", "purple"), 
             alpha = 0.6) + 
  scale_x_continuous(limits = c(0, 126),
                     breaks = c(0, 10, 20, 30, 40, 50, 60, 
                                70, 80, 90, 100, 110, 120)) +
  scale_y_continuous(limits = c(0, 3.4),
                     breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  labs(x = "Residue Number", y = "∆∂ ppm") +
  theme_bw() + # green line = 1 SD and purple line = 2 SD from zero 
  theme(axis.text = element_text(color = "black", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color = "black", size = 1, linetype = 1),
        axis.line.y = element_line(color = "black", size = 1, linetype = 1)
  )

suppressWarnings(ggsave("hFis1_pep213_CSP_plot.pdf", 
       width = 14, height = 10, units = "cm")) 
# delta sign does not print onto pdf and instead generates error message.
# suppressWarnings() hides this error and cleans up compiled R script
ggsave("hFis1_pep213_CSP_plot.png", 
       width = 14, height = 10, units = "cm")

# Write out CSV file for CSP mapping in PyMol   --------------------------------

CSP_color1 <- read_csv("hFis1_pycolor_template.csv")

CSP_color2 <- total_shift %>%
  rename(., "number" = residue) %>%
  select(., number, CSP)

CSP_color <- CSP_color1 %>%
  left_join(., CSP_color2, by = "number") %>%
  replace_na(., list(CSP = 0))

write_csv(CSP_color, "hFis1_pep213_pycolor.csv")

CSP_threshold_sigma <- CSP_color %>%
  filter(., CSP > 0.506)
write_csv(CSP_threshold_sigma, "CSP_both_1_2_sigma.csv")
CSP_threshold_1sigma <- CSP_color %>%
  filter(., CSP > 0.506 & CSP < 1.01)
write_csv(CSP_threshold_1sigma, "CSP_1sigma.csv")
CSP_threshold_2sigma <- CSP_color %>%
  filter(., CSP > 1.01)
write_csv(CSP_threshold_2sigma, "CSP_2sigma.csv")

# Calculate Kd values from non-linear fitting of CSP - all residues   ----------

# calculate CSPs at each [ligand]
kd_fitting <- total_shift %>%
  mutate(., CSP0 = sqrt(((5*(H_0 - H_0))^2) + (N_0 - N_0)^2),
         CSP25 = sqrt(((5*(H_0 - H_25))^2) + (N_0 - N_25)^2),
         CSP50 = sqrt(((5*(H_0 - H_50))^2) + (N_0 - N_50)^2),
         CSP150 = sqrt(((5*(H_0 - H_150))^2) + (N_0 - N_150)^2),
         CSP400 = sqrt(((5*(H_0 - H_400))^2) + (N_0 - N_400)^2),
         CSP800 = sqrt(((5*(H_0 - H_800))^2) + (N_0 - N_800)^2),
         CSP1600 = sqrt(((5*(H_0 - H_1600))^2) + (N_0 - N_1600)^2),
         CSP2000 = sqrt(((5*(H_0 - H_2000))^2) + (N_0 - N_2000)^2)) %>%
  select(., 1, 18, 22:29) %>%
  gather(., "conc", "csp", 3:10) %>%
  mutate(., conc = parse_number(conc))

# visualize CSPs by [ligand] plot; W40 indole used as reporter
kd_fitting %>% 
  filter(., label != "W40") %>%
  ggplot(., aes(x = conc, y = csp)) +
  geom_point() +
  facet_wrap(~ residue) +
  scale_x_continuous(limits = c(0, 2000),
                     breaks = c(0, 500, 1000, 1500, 2000)) +
  labs(x = "[pep213] (µM)", y = "CSP (ppm)") +
  theme_bw()

ggsave("hFis1_pep213_csp_by_res.pdf",
       width = 40, height = 40, units = "cm")

# plot Kd fits onto data; W40 indole used as reporter
kd_fitting %>% 
  filter(., label != "W40") %>%
  ggplot(., aes(x = conc, y = csp, color = as.factor(residue))) +
  geom_point() +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ dmax*((kd+50+x) - ((kd+50+x)^2 - 4*50*x)^0.5)/(2*50),
            method.args = list(start = c(dmax = 0.5,
                                         kd = 5),
                               control = nls.control(maxiter = 100, tol = 1e-6)),
            se = FALSE,       
            fullrange = TRUE,  
            size = 1,
            alpha = 0.6) +        # p held constant at 50 (µM) in model
  facet_wrap(~ residue) +
  labs(x = "[pep213] (µM)", y = "CSP (ppm)") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("hFis1_pep213_csp_by_res_fits_facet.pdf",
       width = 40, height = 40, units = "cm")

# Calculate Kd values from non-linear fitting of CSP - select residues   -------

# Visualize Kd fits onto (select residues) data
kd_fitting %>% 
  filter(., residue == 6 | residue == 8 | residue == 9 |
           residue == 36 | residue == 44 | residue == 45 |
           residue == 47 | residue == 51 | residue == 74 | residue == 76 |
           residue == 75 | residue == 79 | residue == 80 | residue == 83 |
           residue == 111 | residue == 119) %>%
  ggplot(., aes(x = conc, y = csp, color = as.factor(residue))) +
  geom_point() +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ dmax*((kd+50+x) - ((kd+50+x)^2 - 4*50*x)^0.5)/(2*50),
            method.args = list(start = c(dmax = 0.5,
                                         kd = 5),
                               control = nls.control(maxiter = 100, tol = 1e-6)),
            se = FALSE,        
            fullrange = TRUE,  
            size = 1,
            alpha = 0.6) +        # p held constant at 50 (µM) in model
  scale_x_continuous(limits = c(0, 2200),
                     breaks = c(0, 800, 1600)) +
  scale_y_continuous(limits = c(0, 2),
                     breaks = c(0, 0.5, 1.0, 1.5, 2.0)) +
  facet_wrap(~ residue) +
  labs(color = "Residue",
       x = "[pep213] (µM)",
       y = "CSP (ppm)") +
  #scale_color_brewer("Residue", palette = "Dark2") +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none"
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid')
  )

ggsave("hFis1_pep213_Kd_CSP_select_fitting_plot_facet.pdf",
       width = 20, height = 18, units = "cm")

# calculate Kd value from subset of data
kd_select_csp <- kd_fitting %>% 
  filter(., residue == 6 | residue == 8 | residue == 9 |
           residue == 36 | residue == 44 | residue == 45 |
           residue == 47 | residue == 51 | residue == 74 | residue == 76 |
           residue == 75 | residue == 79 | residue == 80 | residue == 83 |
           residue == 111 | residue == 119) %>%
  group_by(., residue) %>%
  do(tidy(nlsLM(
    formula = csp ~ dmax*((kd+50+conc) - ((kd+50+conc)^2 - 4*50*conc)^0.5)/(2*50),
    start = list(dmax = 0.5, kd = 5), 
    trace = TRUE,
    control = nls.control(maxiter = 100, tol = 1e-6),
    data = .))) %>%
  spread(., key = term, value = estimate)

write_csv(kd_select_csp, "hFis1_pep213_Kd_CSP_select_values.csv")

# Evaluate non-linear fitting by inspection of residuals - select residues  ----

# calculate residuals using augment from broom package of same residues for Kd's

# Calculate residuals for select residues
kd_csp_select_residuals <- kd_fitting %>%
  filter(., residue == 6 | residue == 8 | residue == 9 |
           residue == 36 | residue == 44 | residue == 45 |
           residue == 47 | residue == 51 | residue == 74 | residue == 76 |
           residue == 75 | residue == 79 | residue == 80 | residue == 83 |
           residue == 111 | residue == 119) %>%
  group_by(., residue) %>%
  do(augment(nlsLM(
    formula = csp ~ dmax*((kd+50+conc) - ((kd+50+conc)^2 - 4*50*conc)^0.5)/(2*50),
    start = list(dmax = 0.5, kd = 5), 
    trace = TRUE,
    control = nls.control(maxiter = 100, tol = 1e-6),
    data = .))) 

write_csv(kd_csp_select_residuals, "hFis1_pep213_Kd_csp_select_residuals.csv")

kd_csp_select_residuals %>%
  ggplot(., aes(x = conc, y = .resid, color = as.factor(residue))) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linetype = 2) +
  facet_wrap(~ residue) +
  scale_x_continuous(limits = c(0, 2200),
                     breaks = c(0, 800, 1600)) +
  scale_y_continuous(limits = c(-1, 1),
                     breaks = c(-0.6, 0, 0.6)) +
  labs(color = "Residue",
       x = "[pep213] (µM)",
       y = "Residuals") +
  scale_color_hue(h = c(180, 270)) +
  theme_bw() +
  theme(legend.position = "none")
ggsave("hFis1_pep213_Kd_csp_residuals_facet.pdf", 
       width = 20, height = 18, units = "cm")
# 
# kd_csp_select_residuals %>%
#   ggplot(., aes(x = conc, y = .resid, color = as.factor(residue))) +
#   geom_point() +
#   geom_hline(yintercept = 0, color = "red", linetype = 2) +
#   scale_x_continuous(breaks = c(0, 500, 1000, 1500, 2000)) +
#   scale_y_continuous(limits = c(-0.16, 0.16),
#                      breaks = c(-0.15, 0, 0.15)) +
#   labs(x = "[pep213] (µM)", y = "Residuals") +
#   theme_bw() +
#   theme(axis.text = element_text(size = 14, color = "black"),
#         axis.title.x = element_text(size = 12),
#         axis.title.y = element_text(size = 12),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
#         axis.line.y = element_line(colour = 'black', size=1, linetype='solid'),
#         legend.position = "none")
# ggsave("hFis1_pep213_Kd_CSP_select_residuals.pdf", 
#        width = 18, height = 6, units = "cm")

# Calculate Kd from CSPs, using simplified equation (model 2) - select residues ----

# Using model 2: calculate Kd value from select CSP data (only used Kd values > 0)
# Kd values less than 0 are not possible
kd_select_csp_model2 <- kd_fitting %>% 
  filter(., residue == 6 | residue == 8 | residue == 9 |
           residue == 36 | residue == 44 | residue == 45 |
           residue == 47 | residue == 51 | residue == 74 | residue == 76 |
           residue == 75 | residue == 79 | residue == 80 | residue == 83 |
           residue == 111 | residue == 119) %>%
  group_by(., residue) %>%
  do(tidy(nlsLM(
    formula = csp ~ (dmax*conc) / (kd + conc),
    start = list(dmax = 0.5, kd = 5), 
    trace = TRUE,
    control = nls.control(maxiter = 100, tol = 1e-6),
    data = .))) %>%
  spread(., key = term, value = estimate)

write_csv(kd_select_csp_model2, "hFis1_pep213_Kd_CSP_select_values_model2.csv")

# plot resulting Kd fits onto select CSP data
kd_fitting %>% 
  filter(., residue == 6 | residue == 8 | residue == 9 |
           residue == 36 | residue == 44 | residue == 45 |
           residue == 47 | residue == 51 | residue == 74 | residue == 76 |
           residue == 75 | residue == 79 | residue == 80 | residue == 83 |
           residue == 111 | residue == 119) %>%
  ggplot(., aes(x = conc, y = csp, color = as.factor(residue))) +
  geom_point() +
  facet_wrap(~ residue) +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ (dmax*x) / (kd + x),
            method.args = list(start = c(dmax = 0.5,
                                         kd = 5),
                               control = nls.control(maxiter = 100, tol = 1e-6)),
            se = FALSE,        
            fullrange = F,  
            size = 1,
            alpha = 0.6) +        
  scale_x_continuous(breaks = c(0, 500, 1000, 1500, 2000)) +
  scale_y_continuous(limits = c(0, 3),
                     breaks = c(0, 1, 2, 3)) +
  labs(color = "Residue", x = "[pep213] (µM)", y = "CSP (ppm)") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("hFis1_pep213_Kd_CSP_select_fitting_plot_model2_facet.pdf",
       width = 20, height = 18, units = "cm")

# # plot resulting Kd fits onto select CSP data - single plot (no facet)
 # kd_fitting %>%
 #   filter(., residue == 6 | residue == 8 | residue == 9 |
 #         residue == 36 | residue == 44 | residue == 45 |
 #         residue == 47 | residue == 51 | residue == 74 | residue == 76 |
 #         residue == 75 | residue == 79 | residue == 80 | residue == 83 |
 #         residue == 111 | residue == 119) %>%
#   ggplot(., aes(x = conc, y = csp, color = as.factor(residue))) +
#   geom_point() +
#   geom_line(stat = "smooth", method = "nlsLM",
#             formula = y ~ (dmax*x) / (kd + x),
#             method.args = list(start = c(dmax = 0.5,
#                                          kd = 5),
#                                control = nls.control(maxiter = 100, tol = 1e-6)),
#             se = FALSE,        
#             fullrange = F,  
#             size = 1,
#             alpha = 0.6) +        
#   scale_x_continuous(breaks = c(0, 500, 1000, 1500, 2000)) +
#   scale_y_continuous(limits = c(0, 3),
#                      breaks = c(0, 1, 2, 3)) +
#   labs(color = "Residue", x = "[pep213] (µM)", y = "CSP (ppm)") +
#   theme_bw() +
#   theme(axis.text = element_text(size = 14, color = "black"),
#         axis.title.x = element_text(size = 12),
#         axis.title.y = element_text(size = 12),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
#         axis.line.y = element_line(colour = 'black', size=1, linetype='solid')
#   )
# 
# ggsave("hFis1_pep213_Kd_CSP_select_fitting_plot_model2.pdf",
#        width = 20, height = 14, units = "cm")

# Using model 2: Calculate residuals for select residues
kd_csp_select_residuals_model2 <- kd_fitting %>%
  filter(., residue == 6 | residue == 8 | residue == 9 |
           residue == 36 | residue == 44 | residue == 45 |
           residue == 47 | residue == 51 | residue == 74 | residue == 76 |
           residue == 75 | residue == 79 | residue == 80 | residue == 83 |
           residue == 111 | residue == 119) %>%
  group_by(., residue) %>%
  do(augment(nlsLM(
    formula = csp ~ (dmax*conc) / (kd + conc),
    start = list(dmax = 0.5, kd = 5), 
    trace = TRUE,
    control = nls.control(maxiter = 100, tol = 1e-6),
    data = .))) 

write_csv(kd_csp_select_residuals_model2,
          "hFis1_pep213_Kd_csp_select_residuals_model2.csv")

kd_csp_select_residuals_model2 %>%
  ggplot(., aes(x = conc, y = .resid, color = as.factor(residue))) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linetype = 2) +
  facet_wrap(~ residue) +
  scale_x_continuous(limits = c(0, 2200),
                     breaks = c(0, 800, 1600)) +
  scale_y_continuous(limits = c(-1, 1),
                     breaks = c(-0.6, 0, 0.6)) +
  labs(x = "[pep213] (µM)", y = "Residuals") +
  scale_color_hue(h = c(180, 270)) +
  theme_bw() +
  theme(legend.position = "none")
ggsave("hFis1_pep213_Kd_csp_select_residuals_facet_model2.pdf", 
       width = 20, height = 18, units = "cm")

# Kd values from different fits   ----------------------------------------------

# Model1, Kd determined by fitting select residues (Change in Amplitude)
print(c(mean(select_kd_int$kd, na.rm = T), "±", sd(select_kd_int$kd, na.rm = T),
        "mean ± sd", "Kd from Intensity", "ligand depletion"))

# Model2, Kd determined by fitting select residues (Change in Amplitude)
print(c(mean(select_kd_int_model2$kd, na.rm = T), "±",
        sd(select_kd_int_model2$kd, na.rm = T),
        "mean ± sd", "Kd from Intensity", "simple Kd eq. / model 2"))

# Model1, Kd determined by fitting select residues (CSP against [pep213])
print(c(mean(kd_select_csp$kd, na.rm = T), "±", sd(kd_select_csp$kd, na.rm = T),
        "mean ± sd", "Kd from CSP, select residues", "ligand depletion"))

# Model2, Kd determined by fitting select residues (CSP against [pep213])
print(c(mean(kd_select_csp_model2$kd, na.rm = T), "±", 
        sd(kd_select_csp_model2$kd, na.rm = T),
        "mean ± sd", "Kd from CSP, select residues", "simple Kd eq. / model 2"))
