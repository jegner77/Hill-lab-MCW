#
# intrinsic TRP fluorescence of dN-hFis1 monitored
# pep13 titrated against ∆N-hFis1 residues 9 - 125
# [dN-hFis1] = 10 µM 
# [pep13] = 0 - 2 mM
# data collected on 20180209
# Buffer: 100 mM HEPES pH 7.4, 200 mM NaCl, 1 mM DTT, 0.02% Na-Azide
# TRP excitation = 295 nm
# emission = 300 - 400nm
# 

library(tidyverse)
library(broom)
library(readxl)
library(minpack.lm)
library(gridExtra)

# Import, tidy, and subtract peptide alone emission scan    --------------------

raw <- read_excel("20180209_dN_hFis1_pep13_TRP_fluor_data.xlsx", skip = 1)

# emission spectra without peptide alone subtraction
results <- raw %>%
  gather(., "tmp", "Fluorescence", 2:31) %>%
  separate(., "tmp", into = c("conc", "dN_hFis1", "ex_slit_width")) %>%
  separate(., "conc", into = c("conc", "del"), sep = "pep") %>%
  select(., -del) %>%
  mutate(., conc = parse_number(conc),
         dN_hFis1 = replace(dN_hFis1, dN_hFis1 == "dN", "+"),
         dN_hFis1 = replace(dN_hFis1, dN_hFis1 == "no", "-"),
         ex_slit_width = parse_number(ex_slit_width))

# difference emission spectra (ex slit width = 4): 
# difference emission spectra = dN-hFis1 em spectra - pep13 alone em spectra
corrected <- raw %>%
  mutate(., "0" = `0pep13_dN_4` - `0pep13_no_4`,
         "1" = `1pep13_dN_4` - `1pep13_no_4`,
         "3" = `3pep13_dN_4` - `3pep13_no_4`,
         "7" = `7pep13_dN_4` - `7pep13_no_4`,
         "10" = `10pep13_dN_4` - `10pep13_no_4`,
         "30" = `30pep13_dN_4` - `30pep13_no_4`,
         "70" = `70pep13_dN_4` - `70pep13_no_4`,
         "100" = `100pep13_dN_4` - `100pep13_no_4`,
         "300" = `300pep13_dN_4` - `300pep13_no_4`,
         "700" = `700pep13_dN_4` - `700pep13_no_4`,
         "1000" = `1000pep13_dN_4` - `1000pep13_no_4`,
         "1300" = `1300pep13_dN_4` - `1300pep13_no_4`,
         "2000" = `2000pep13_dN_4` - `2000pep13_no_4`) %>%
  select(., 1, 32:44) %>%
  gather(., "conc", "Fluorescence", 2:14) %>%
  mutate(., conc = parse_number(conc))

write_csv(corrected, "dN_hFis1_pep13_corrected.csv")

# Buffer subtraction from all emission spectra
buffer_corrected <- raw %>%
  mutate(., "0" = `0pep13_dN_4` - `0pep13_no_4`,
         "1" = `1pep13_dN_4` - `0pep13_no_4`,
         "3" = `3pep13_dN_4` - `0pep13_no_4`,
         "7" = `7pep13_dN_4` - `0pep13_no_4`,
         "10" = `10pep13_dN_4` - `0pep13_no_4`,
         "30" = `30pep13_dN_4` - `0pep13_no_4`,
         "70" = `70pep13_dN_4` - `0pep13_no_4`,
         "100" = `100pep13_dN_4` - `0pep13_no_4`,
         "300" = `300pep13_dN_4` - `0pep13_no_4`,
         "700" = `700pep13_dN_4` - `0pep13_no_4`,
         "1000" = `1000pep13_dN_4` - `0pep13_no_4`,
         "1300" = `1300pep13_dN_4` - `0pep13_no_4`,
         "2000" = `2000pep13_dN_4` - `0pep13_no_4`) %>%
  select(., 1, 32:44) %>%
  gather(., "conc", "Fluorescence", 2:14) %>%
  mutate(., conc = parse_number(conc))

# Visualize varied excitation slit width data  ---------------------------------

# first, plot 2 mM pep13 and dN-hFis1 emission spectra as function of slit width
results %>%
  filter(., conc == 0 & dN_hFis1 == "+" | conc == 2000 & dN_hFis1 == "-") %>% 
  ggplot(., aes(x = wavelength, y = Fluorescence, 
                shape = as.factor(ex_slit_width), color = as.factor(dN_hFis1))) +
  geom_point(size = 1) +
  geom_line(stat = "smooth", method = "loess", 
            span = 0.2, size = 1, alpha = 0.4) +
  geom_vline(xintercept = 333, color = "red", linetype = 2) +
  scale_x_continuous(limits = c(300, 400),
                     breaks = c(300, 320, 340, 360, 380, 400)) +
  #scale_y_continuous(limits = c(0, 320000),
  #                   breaks = c(0, 100000, 200000, 300000)) +
  labs(shape = "Excitation Slit Width:",
       #color = "dN-hFis1:",
       title = "TRP Emission Spectra as function of excitation slit width", 
       x = "Wavelength (nm)",
       y = "Fluorescence (AU)") +
  guides(shape = guide_legend(override.aes = list(size = 2))) +
  scale_color_manual("",
                     values = c("-" = "#51CC2C",
                                "+" = "#C22CCC"),
                     labels = c("-" = "2 mM pep13",
                                "+" = "10 µM dN-hFis1")) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        #axis.title.x = element_text(size = 18),
        #axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14, color = "black"),
        legend.position = "top"
  )

ggsave("20180209_dN_hFis1_pep13_TRP_em_spectra_slit_widths.pdf", 
       width = 20, height = 12, dpi = 300, units = "cm")

# zoomed-in view to more easily observe emission spectra profiles
results %>%
  filter(., conc == 0 & dN_hFis1 == "+" | conc == 2000 & dN_hFis1 == "-") %>% 
  ggplot(., aes(x = wavelength, y = Fluorescence, 
                shape = as.factor(ex_slit_width), color = as.factor(dN_hFis1))) +
  geom_point(size = 1) +
  geom_line(stat = "smooth", method = "loess", 
            span = 0.2, size = 1, alpha = 0.4) +
  geom_vline(xintercept = 333, color = "red", linetype = 2) +
  scale_x_continuous(limits = c(300, 400),
                     breaks = c(300, 320, 340, 360, 380, 400)) +
  scale_y_continuous(limits = c(0, 8000),
                     breaks = c(0, 2000, 4000, 6000, 8000)) +
  labs(shape = "Excitation Slit Width:",
       #color = "dN-hFis1:",
       title = "TRP Emission Spectra as function of excitation slit width", 
       x = "Wavelength (nm)",
       y = "Fluorescence (AU)") +
  guides(shape = guide_legend(override.aes = list(size = 2))) +
  scale_color_manual("",
                     values = c("-" = "#51CC2C",
                                "+" = "#C22CCC"),
                     labels = c("-" = "2 mM pep13",
                                "+" = "10 µM dN-hFis1")) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        #axis.title.x = element_text(size = 18),
        #axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14, color = "black"),
        legend.position = "top"
  )

ggsave("20180209_dN_hFis1_pep13_TRP_em_spectra_slit_widths_zoom.pdf", 
       width = 20, height = 12, dpi = 300, units = "cm")

# Does excitation slit width improve signal to noise between +/- pep13? --------
results %>%
  filter(., conc == 0 & dN_hFis1 == "+" & wavelength == 333 | 
           conc == 2000 & dN_hFis1 == "-" & wavelength == 333)

results %>%
  filter(., conc == 0 & dN_hFis1 == "+" & wavelength == 340 | 
           conc == 2000 & dN_hFis1 == "-" & wavelength == 340)

results %>%
  filter(., conc == 0 & dN_hFis1 == "+" & wavelength == 320 | 
           conc == 2000 & dN_hFis1 == "-" & wavelength == 320)
  
# the ratio (at a specific wavelength) between apo dN-hFis1 and 2000 uM pep13
# are equal regardless of excitation slit width

# Visualize emission spectra, prior to viewing difference emission spectra  ----

# plot full window (wavelength 300 - 400) of all emission spectra at ex slit = 4
results %>%
  filter(., ex_slit_width == 4) %>%
  ggplot(., aes(x = wavelength, y = Fluorescence,
                color = as.factor(conc), shape = as.factor(dN_hFis1))) +
  geom_point(size = 1) +
  geom_line(stat = "smooth", method = "loess", span = 0.2, alpha = 0.4) +
  scale_x_continuous(limits = c(300, 400),
                     breaks = c(300, 310, 320, 330, 340, 350, 
                                360, 370, 380, 390, 400)) +
  scale_y_continuous(limits = c(0, 45000),
                     breaks = c(0, 10000, 20000, 30000, 40000)) +
  labs(color = "[pep13] (µM)",
       shape = "dN-hFis1",
       title = "TRP Emission Spectra of pep13 and dN-hFis1 + pep13", 
       x = "Wavelength (nm)",
       y = "Fluorescence (AU)") +
  guides(shape = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        #axis.title.x = element_text(size = 18),
        #axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14, color = "black")
  )

ggsave("20180209_dN_hFis1_pep13_TRP_em_spectra_all.pdf", 
       width = 18, height = 12, dpi = 300, units = "cm")

# plot full window (wavelength 300 - 400) of pep13 only emission spectra 
  # at ex slit = 4
results %>%
  filter(., ex_slit_width == 4 & dN_hFis1 == "-") %>%
  ggplot(., aes(x = wavelength, y = Fluorescence,
                color = as.factor(conc), shape = as.factor(dN_hFis1))) +
  geom_point(size = 1) +
  geom_line(stat = "smooth", method = "loess", span = 0.2, alpha = 0.4) +
  scale_x_continuous(limits = c(300, 400),
                     breaks = c(300, 310, 320, 330, 340, 350, 
                                360, 370, 380, 390, 400)) +
  scale_y_continuous(limits = c(0, 45000),
                     breaks = c(0, 10000, 20000, 30000, 40000)) +
  labs(color = "[pep13] (µM)",
       shape = "dN-hFis1",
       title = "TRP Emission Spectra of pep13", 
       x = "Wavelength (nm)",
       y = "Fluorescence (AU)") +
  guides(shape = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        #axis.title.x = element_text(size = 18),
        #axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14, color = "black")
  )

ggsave("20180209_dN_hFis1_pep13_TRP_em_spectra_pep13_only.pdf", 
       width = 18, height = 12, dpi = 300, units = "cm")

# plot full window (wavelength 300 - 400) of dN-hFis1 only emission spectra 
  # at ex slit = 4
results %>%
  filter(., ex_slit_width == 4 & dN_hFis1 == "+") %>%
  ggplot(., aes(x = wavelength, y = Fluorescence,
                color = as.factor(conc), shape = as.factor(dN_hFis1))) +
  geom_point(size = 1) +
  geom_line(stat = "smooth", method = "loess", span = 0.2, alpha = 0.4) +
  scale_x_continuous(limits = c(300, 400),
                     breaks = c(300, 310, 320, 330, 340, 350, 
                                360, 370, 380, 390, 400)) +
  scale_y_continuous(limits = c(0, 45000),
                     breaks = c(0, 10000, 20000, 30000, 40000)) +
  labs(color = "[pep13] (µM)",
       shape = "dN-hFis1",
       title = "TRP Emission Spectra of dN-hFis1 + pep13", 
       x = "Wavelength (nm)",
       y = "Fluorescence (AU)") +
  guides(shape = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        #axis.title.x = element_text(size = 18),
        #axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14, color = "black")
  )

ggsave("20180209_dN_hFis1_pep13_TRP_em_spectra_dN_hFis1.pdf", 
       width = 18, height = 12, dpi = 300, units = "cm")

# Re-plot em spectra with narrower wavelength window   -------------------------
# plot narrow window (wavelength 310 - 370) of all emission spectra 
 # at ex slit = 4
results %>%
  filter(., ex_slit_width == 4) %>%
  ggplot(., aes(x = wavelength, y = Fluorescence,
                color = as.factor(conc), shape = as.factor(dN_hFis1))) +
  geom_point(size = 1) +
  geom_line(stat = "smooth", method = "loess", span = 0.2, alpha = 0.4) +
  scale_x_continuous(limits = c(310, 370),
                     breaks = c(300, 310, 320, 330, 340, 350, 
                                360, 370, 380, 390, 400)) +
  scale_y_continuous(limits = c(0, 20000),
                     breaks = c(0, 5000, 10000, 15000, 20000)) +
  labs(color = "[pep13] (µM)",
       shape = "dN-hFis1",
       title = "TRP Emission Spectra of dN-hFis1 + pep13", 
       x = "Wavelength (nm)",
       y = "Fluorescence (AU)") +
  guides(shape = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        #axis.title.x = element_text(size = 18),
        #axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14, color = "black")
  )

ggsave("20180209_dN_hFis1_pep13_TRP_em_spectra_all_zoom.pdf", 
       width = 18, height = 12, dpi = 300, units = "cm")

# plot narrow window (wavelength 310 - 370) of pep13 only emission spectra 
  # at ex slit = 4
results %>%
  filter(., ex_slit_width == 4 & dN_hFis1 == "-") %>%
  ggplot(., aes(x = wavelength, y = Fluorescence,
                color = as.factor(conc), shape = as.factor(dN_hFis1))) +
  geom_point(size = 1) +
  geom_line(stat = "smooth", method = "loess", span = 0.2, alpha = 0.4) +
  scale_x_continuous(limits = c(310, 370),
                     breaks = c(300, 310, 320, 330, 340, 350, 
                                360, 370, 380, 390, 400)) +
  scale_y_continuous(limits = c(0, 20000),
                     breaks = c(0, 5000, 10000, 15000, 20000)) +
  labs(color = "[pep13] (µM)",
       shape = "dN-hFis1",
       title = "TRP Emission Spectra of pep13", 
       x = "Wavelength (nm)",
       y = "Fluorescence (AU)") +
  guides(shape = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        #axis.title.x = element_text(size = 18),
        #axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14, color = "black")
  )

ggsave("20180209_dN_hFis1_pep13_TRP_em_spectra_pep13_only_zoom.pdf", 
       width = 18, height = 12, dpi = 300, units = "cm")

# plot narrow window (wavelength 310 - 370) of dN-hFis1 only emission spectra 
  # at ex slit = 4
results %>%
  filter(., ex_slit_width == 4 & dN_hFis1 == "+") %>%
  ggplot(., aes(x = wavelength, y = Fluorescence,
                color = as.factor(conc), shape = as.factor(dN_hFis1))) +
  geom_point(size = 1) +
  geom_line(stat = "smooth", method = "loess", span = 0.2, alpha = 0.4) +
  scale_x_continuous(limits = c(310, 370),
                     breaks = c(300, 310, 320, 330, 340, 350, 
                                360, 370, 380, 390, 400)) +
  scale_y_continuous(limits = c(0, 20000),
                     breaks = c(0, 5000, 10000, 15000, 20000)) +
  labs(color = "[pep13] (µM)",
       shape = "dN-hFis1",
       title = "TRP Emission Spectra of dN-hFis1 + pep13", 
       x = "Wavelength (nm)",
       y = "Fluorescence (AU)") +
  guides(shape = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        #axis.title.x = element_text(size = 18),
        #axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14, color = "black")
  )

ggsave("20180209_dN_hFis1_pep13_TRP_em_spectra_dN_hFis1_zoom.pdf", 
       width = 18, height = 12, dpi = 300, units = "cm")

# Visualize difference emission spectra   --------------------------------------

# plot full window (wavelength 310 - 400) of differene emission spectra
  # at ex slit = 4
  # some problem occured with 30 uM pep13 measurement so it was excluded
corrected %>%
  filter(., wavelength >= 310 & conc != 30) %>%
  ggplot(., aes(x = wavelength, y = Fluorescence,
                color = as.factor(conc))) +
  geom_point(size = 1) +
  geom_line(stat = "smooth", method = "loess", span = 0.2, alpha = 0.4) +
  scale_x_continuous(limits = c(310, 400),
                     breaks = c(310, 320, 330, 340, 350, 
                                360, 370, 380, 390, 400)) +
  scale_y_continuous(limits = c(0, 8000),
                     breaks = c(0, 2000, 4000, 6000, 8000)) +
  labs(color = "[pep13] (µM)",
       shape = "dN-hFis1",
       title = "Difference TRP Emission Spectra of dN-hFis1 + pep13", 
       x = "Wavelength (nm)",
       y = "Fluorescence (AU)") +
  guides(shape = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        #axis.title.x = element_text(size = 18),
        #axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14, color = "black")
  )

ggsave("20180209_dN_hFis1_pep13_TRP_em_spectra_all_corrected.pdf", 
       width = 18, height = 12, dpi = 300, units = "cm")

# plot full window (wavelength 310 - 400) of differene emission spectra
# at ex slit = 4
# some problem occured with 30 uM pep13 measurement so it was excluded
corrected %>%
  filter(., wavelength >= 310 & conc != 30) %>%
  ggplot(., aes(x = wavelength, y = Fluorescence,
                color = as.factor(conc))) +
  geom_line(stat = "smooth", method = "loess", span = 0.2, alpha = 0.8) +
  scale_x_continuous(limits = c(310, 400),
                     breaks = c(310, 320, 330, 340, 350, 
                                360, 370, 380, 390, 400)) +
  scale_y_continuous(limits = c(0, 8000),
                     breaks = c(0, 2000, 4000, 6000, 8000)) +
  labs(color = "[pep13] (µM)",
       shape = "dN-hFis1",
       title = "Difference TRP Emission Spectra of dN-hFis1 + pep13", 
       x = "Wavelength (nm)",
       y = "Fluorescence (AU)") +
  guides(shape = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        #axis.title.x = element_text(size = 18),
        #axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14, color = "black")
  )

ggsave("20180209_dN_hFis1_pep13_TRP_em_spectra_all_corrected_lines.pdf", 
       width = 18, height = 12, dpi = 300, units = "cm")

# plot window (wavelength 310 - 370) of differene emission spectra
  # at ex slit = 4, to be used for Kd determination
  # some problem occured with 30 uM pep13 measurement so it was excluded
corrected %>%
  filter(., wavelength >= 310 & wavelength <= 370 & conc != 30) %>%
  ggplot(., aes(x = wavelength, y = Fluorescence,
                color = as.factor(conc))) +
  geom_point(size = 1) +
  geom_line(stat = "smooth", method = "loess", span = 0.2, alpha = 0.4) +
  scale_x_continuous(limits = c(310, 370),
                     breaks = c(310, 320, 330, 340, 350, 360, 370)) +
  scale_y_continuous(limits = c(0, 8000),
                     breaks = c(0, 2000, 4000, 6000, 8000)) +
  labs(color = "[pep13] (µM)",
       shape = "dN-hFis1",
       title = "Difference TRP Emission Spectra of dN-hFis1 + pep13", 
       x = "Wavelength (nm)",
       y = "Fluorescence (AU)") +
  guides(shape = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        #axis.title.x = element_text(size = 18),
        #axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14, color = "black")
  )

ggsave("20180209_dN_hFis1_pep13_TRP_em_spectra_all_corrected_window.pdf", 
       width = 18, height = 12, dpi = 300, units = "cm")

# Determine Kd by average emission wavelength   --------------------------------

# plot window to determine Kd 
buffer_corrected %>%
  filter(., wavelength >= 310 & wavelength <= 370) %>%
  ggplot(., aes(x = wavelength, y = Fluorescence, color = as.factor(conc))) +
  geom_point(size = 1) +
  geom_smooth(method = "loess", span = 0.3, size = 0.5) +
  scale_x_continuous(limits = c(310, 370),
                     breaks = c(310, 320, 330, 340, 350, 360, 370)) +
  scale_y_continuous(limits = c(0, 17000),
                     breaks = c(0, 5000, 10000, 15000)) +
  labs(color = "[pep13] (µM)",
       title = "TRP Emission Spectra of dN-hFis1 + pep13", 
       x = "Wavelength (nm)",
       y = "Fluorescence (AU)") +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        #axis.title.x = element_text(size = 18),
        #axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=1, linetype='solid'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, color = "black")
  )

ggsave("20180209_dN_hFis1_pep13_TRP_em_spectra_avgEM_window.pdf", 
       width = 18, height = 12, dpi = 300, units = "cm")

# Determination of Kd by computing average emission wavelength against ln[conc]
  # Natural log (ln) is computed with log(x)
avg_em_lamda <- buffer_corrected %>%
  filter(., wavelength >= 310 & wavelength <= 370) %>%
  group_by(., conc) %>%
  mutate(., "IxF" = Fluorescence*wavelength) %>%
  summarise_at(., vars(IxF, Fluorescence), funs(sum)) %>%
  rename(., "sumF" = Fluorescence) %>%
  mutate(., "avgLAMDA" = IxF/sumF,
         "lnC" = log(conc)) #computes natural log of concentration

# visualize model fitting
plotCSM <- avg_em_lamda %>%
  filter(., lnC >= 0) %>%
  ggplot(., aes(x = lnC, y = avgLAMDA)) +
  geom_point(shape = 15) +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ a + ((b - a) / (1 + exp((kd - x)/c))),
            method.args = list(start = c(a = 330,
                                         b = 340,
                                         c = -1,
                                         kd = 10),
                               control = nls.control(maxiter = 100, tol = 1e-6)),
            se = FALSE,        
            fullrange = TRUE,  
            size = 1,
            alpha = 0.6, 
            color = "purple") +
  scale_x_continuous(limits = c(0, 8),
                     breaks = c(0, 2, 4, 6, 8)) +
  #scale_y_continuous(limits = c(330, 342),
   #                  breaks = c(330, 332, 334, 336, 338, 340, 342)) +
  labs(x = "ln([pep13]) (µM)", y = "Average Emission Wavelength (nm)") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )
plotCSM
ggsave("20180209_dN_hFis1_pep13_TRP_Kd_plot.pdf", 
       width = 18, height = 12, dpi = 300, units = "cm")

# calculate Kd from model fitting
Kd_fit <- avg_em_lamda %>%
  do(tidy(nlsLM(avgLAMDA ~ a + ((b - a) / (1 + exp((kd - lnC)/c))),
                start = list(a = 330, b = 340, c = -1, kd = 10), trace = TRUE, 
                control = nls.control(maxiter = 200, tol = 1e-6), .)))

write_csv(Kd_fit, "20180209_dN_hFis1_pep13_Kd_fit.csv")

Kd_fit

# Inspect fit (Avg emission wavelength)   --------------------------------------
Kd_residuals <- avg_em_lamda %>%
  do(augment(nlsLM(avgLAMDA ~ a + ((b - a) / (1 + exp((kd - lnC)/c))),
                start = list(a = 330, b = 340, c = -1, kd = 10), trace = TRUE, 
                control = nls.control(maxiter = 200, tol = 1e-6), .)))

write_csv(Kd_residuals, "20180209_dN_hFis1_pep13_Kd_residuals.csv")

Kd_residuals

res_csm <- Kd_residuals %>%
  filter(., lnC >= 0) %>%
  ggplot(., aes(x = lnC, y = .resid)) +
  geom_point(shape = 15) +
  geom_hline(yintercept = 0, color = "red", linetype = 2) +
  scale_x_continuous(limits = c(0, 8),
                     breaks = c(0, 2, 4, 6, 8)) +
  scale_y_continuous(limits = c(-0.17, 0.17),
                     breaks = c(-0.15, 0, 0.15)) +
  labs(x = "ln [pep13] (µM)", y = "Residuals") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )
res_csm
ggsave("20180209_dN_hFis1_pep13_TRP_Kd_residuals_plot.pdf", 
       width = 18, height = 6, dpi = 300, units = "cm")

# Final plot of central spectral mass -------------------------------------

# Generate layout matrix to plot residuals at 1/3 size of regression plot
lay <- cbind(c(1, 2, 2))

plot_csm <- grid.arrange(grobs = list(res_csm, plotCSM), layout_matrix = lay)

suppressWarnings(ggsave("final_plot_dN_pep13_center_spectral_mass.pdf", plot_csm,
                        width = 10, height = 10, units = "cm"))
ggsave("final_plot_dN_pep13_centrer_spectral_mass.png", plot_csm,
       width = 10, height = 10, units = "cm")

# Determine Kd by change in inverse Fluorescence intensity    --------------------------

# Determination of Kd by change in Fluorescence Intensity against ln[conc]
# After analysis, this is not an appropriate method to fit this data
# compute change in Fluorescence (from apo state)
tmp <- buffer_corrected %>%
  filter(., wavelength >= 320 & wavelength <= 370 & conc == 0) %>%
  mutate(., "zero" = Fluorescence)
  
df_changeF <- buffer_corrected %>%
  filter(., wavelength >= 320 & wavelength <= 370) %>%
  left_join(., tmp, by = c("wavelength")) %>%
  select(., wavelength, conc.x, Fluorescence.x, zero) %>%
  rename(., "conc" = conc.x,
         "Fluorescence" = Fluorescence.x) %>%
  mutate(., "changeF" = Fluorescence - zero,
         "inv_F" = 1/Fluorescence, 
         "inv_conc" = 1/conc,
         "lnC" = log(conc))

# visualize model fitting
df_changeF %>%
  filter(., wavelength >= 320 & wavelength <= 355 & lnC >= 0) %>%
  ggplot(., aes(x = lnC, y = inv_F)) +
  geom_point() +
  facet_wrap(~ wavelength) +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ a + ((b - a) / (1 + exp((kd - x)/c))),
            method.args = list(start = c(a = 0.0001,
                                         b = 0.0002,
                                         c = -2,
                                         kd = 10),
                               control = nls.control(maxiter = 100, tol = 1e-6)),
            se = FALSE,        
            fullrange = TRUE,  
            size = 1,
            alpha = 0.6, 
            color = "purple") +
  scale_x_continuous(limits = c(0, 10),
                     breaks = c(0, 2, 4, 6, 8, 10)) +
  #scale_y_continuous(limits = c(0.00008, 0.0003),
   #                  breaks = c(0.0001, 0.0002, 0.0003),
    #                 labels = scales::comma) +
  labs(title = "Kd Determination by nlsLM of 1/F vs. ln[pep13] at λ",
       x = "ln [pep13] (µM)", y = "1/F (AU)") +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )

ggsave("20180209_dN_hFis1_pep13_TRP_Kd_plot_intensity.png", 
       width = 28, height = 24, dpi = 300, units = "cm")

# calculate Kd from model fitting
Kd_intensity <-  df_changeF %>%
  filter(., wavelength >= 320 & wavelength <= 355 & lnC >= 0) %>%
  group_by(., wavelength) %>%
  do(tidy(nlsLM(inv_F ~ a + ((b - a) / (1 + exp((kd - lnC)/c))),
                start = list(a = 0.0001, b = 0.0002, c = -2, kd = 10), 
          trace = TRUE,
          control = nls.control(maxiter = 200, tol = 1e-6), .))) %>%
  spread(., key = term, value = estimate)

write_csv(Kd_intensity, "20180209_dN_hFis1_pep13_Kd_intensity.csv")

Kd_intensity

print(c(mean(Kd_intensity$kd, na.rm = T), "±", sd(Kd_intensity$kd, na.rm = T),
        "mean ± sd", "Kd from nlsLM() fit of 1/F vs. ln[pep13]"))

# Inspect fit (change in inverse Fluorescence intensity)    ----------------------------

Kd_intensity_residuals <-  df_changeF %>%
  filter(., wavelength >= 320 & wavelength <= 355 & lnC > 0) %>%
  group_by(., wavelength) %>%
  do(augment(nlsLM(inv_F ~ a + ((b - a) / (1 + exp((kd - lnC)/c))),
                start = list(a = 0.0001, b = 0.0002, c = -2, kd = 10), 
                trace = TRUE,
                control = nls.control(maxiter = 200, tol = 1e-6), .)))

write_csv(Kd_intensity_residuals, "20180209_dN_hFis1_pep13_Kd_intensity_residuals.csv")

Kd_intensity_residuals

Kd_intensity_residuals %>%
  filter(., lnC > 0) %>%
  ggplot(., aes(x = lnC, y = .resid)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linetype = 2) +
  facet_wrap(~ wavelength) +
  scale_x_continuous(limits = c(0, 10),
                     breaks = c(0, 2, 4, 6, 8, 10)) +
  #scale_y_continuous(limits = c(-12e-6, 12e-6),
   #                  breaks = c(-1e-5, 0, 1e-5)) +
  labs(title = "Residuals from Kd fit by nlsLM of 1/F vs. ln[pep13] at λ",
        x = "ln [pep13] (µM)", y = "Residuals") +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )

ggsave("20180209_dN_hFis1_pep13_TRP_Kd_intensity_residuals_plot.png", 
       width = 28, height = 24, dpi = 300, units = "cm")

# Determine Kd by change in Fluorescence Intensity of difference em spectra ----

# compute change in Fluorescence of difference em spectra (from each [titrant])
tmp2 <- corrected %>%
  filter(., wavelength >= 310 & wavelength <= 370 & conc == 0) %>%
  mutate(., "zero" = Fluorescence)

df_changeF_corrected <- corrected %>%
  filter(., wavelength >= 310 & wavelength <= 370 & conc != 30) %>%
  left_join(., tmp2, by = c("wavelength")) %>%
  select(., wavelength, conc.x, Fluorescence.x, zero) %>%
  rename(., "conc" = conc.x,
         "Fluorescence" = Fluorescence.x) %>%
  mutate(., "changeF" = Fluorescence - zero,
         "inv_F" = 1/Fluorescence, 
         "inv_conc" = 1/conc,
         "lnC" = log(conc))

# visualize model fitting
df_changeF_corrected %>%
  filter(., wavelength >= 310 & wavelength <= 370 & changeF >= 0) %>%
  ggplot(., aes(x = conc, y = changeF)) +
  geom_point() +
  facet_wrap(~ wavelength) +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ (Fmax*x)/(kd+x),
            method.args = list(start = c(Fmax = 2000,
                                         kd = 1),
                               control = nls.control(maxiter = 100, tol = 1e-6)),
            se = FALSE,        
            fullrange = TRUE, 
            size = 1,
            alpha = 0.6,
            color = "purple") +        
  scale_x_continuous(limits = c(0, 2000),
                     breaks = c(0, 800, 1600)) +
  scale_y_continuous(limits = c(0, 3030),
                     breaks = c(0, 1000, 2000, 3000)) +
  labs(title = "Kd Determination by nlsLM of ∆F vs. [pep13] at λ",
       x = "[pep13] (µM)", y = "∆F (AU)") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )

ggsave("20180209_dN_hFis1_pep13_TRP_Kd_plot_changeF_facet.png", 
        width = 30, height = 24, dpi = 300, units = "cm")

# visualize model fitting
df_changeF_corrected %>%
  filter(., wavelength >= 330 & wavelength <= 335 & changeF >= 0) %>%
  ggplot(., aes(x = conc, y = changeF, color = factor(wavelength))) +
  geom_point() +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ (Fmax*x)/(kd+x),
            method.args = list(start = c(Fmax = 2000,
                                         kd = 1),
                               control = nls.control(maxiter = 100, tol = 1e-6)),
            se = FALSE,        
            fullrange = TRUE, 
            size = 1,
            alpha = 0.6) +
  scale_x_continuous(limits = c(0, 2000),
                     breaks = c(0, 500, 1000, 1500, 2000)) +
  scale_y_continuous(limits = c(0, 3030),
                     breaks = c(0, 1000, 2000, 3000)) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Kd Determination by nlsLM of ∆F vs. [pep13] at λ",
       color = "Wavelength (nm)",
       x = "[pep13] (µM)", y = "∆F (AU)") +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )

ggsave("20180209_dN_hFis1_pep13_TRP_Kd_plot_changeF_330-335.png", 
       width = 20, height = 12, dpi = 300, units = "cm")
suppressWarnings(ggsave("20180209_dN_hFis1_pep13_TRP_Kd_plot_changeF_330-335.pdf", 
                        width = 20, height = 12, dpi = 300, units = "cm"))

# calculate Kd from model fitting
Kd_changeF <-  df_changeF_corrected %>%
  filter(., wavelength >= 330 & wavelength <= 335) %>%
  group_by(., wavelength) %>%
  do(tidy(nlsLM(changeF ~ (Fmax*conc)/(kd+conc),
                start = list(Fmax = 2000, kd = 1), 
                trace = TRUE,
                control = nls.control(maxiter = 100, tol = 1e-6), .))) %>%
  spread(., key = term, value = estimate)

write_csv(Kd_changeF, "20180209_dN_hFis1_pep13_Kd_changeF_330-335.csv")

Kd_changeF

print(c(mean(Kd_changeF$kd, na.rm = T), "±", sd(Kd_changeF$kd, na.rm = T),
        "mean ± sd", "Kd from nlsLM() of ∆F vs. [pep13]"))

# Inspect fit (change in Fluorescence intensity) of diff. em spectra   ---------

Kd_changeF_residuals <- df_changeF_corrected %>%
  filter(., wavelength >= 330 & wavelength <= 335) %>%
  group_by(., wavelength) %>%
  do(augment(nlsLM(changeF ~ (Fmax*conc)/(kd+conc),
                start = list(Fmax = 2000, kd = 1), 
                trace = TRUE,
                control = nls.control(maxiter = 100, tol = 1e-6), .)))

write_csv(Kd_changeF_residuals, "20180209_dN_hFis1_pep13_Kd_changeF_residuals_330-335.csv")

Kd_changeF_residuals

Kd_changeF_residuals %>%
  ggplot(., aes(x = conc, y = .resid, color = factor(wavelength))) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linetype = 2) +
  facet_wrap(~ wavelength) +
  scale_x_continuous(limits = c(0, 2000),
                     breaks = c(0, 800, 1600)) +
  scale_y_continuous(limits = c(-800, 800),
                     breaks = c(-600, 0, 600)) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Residuals from Kd fit by nlsLM of ∆F vs. [pep13] at λ",
       color = "Wavelength (nm)",
       x = "[pep13] (µM)", y = "Residuals") +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )

ggsave("20180209_dN_hFis1_pep13_TRP_Kd_changeF_residuals_facet_330-335.png", 
       width = 28, height = 24, dpi = 300, units = "cm")

Kd_changeF_residuals %>%
  ggplot(., aes(x = conc, y = .resid, color = factor(wavelength))) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linetype = 2) +
  scale_x_continuous(limits = c(0, 2000),
                     breaks = c(0, 500, 1000, 1500, 2000)) +
  scale_y_continuous(limits = c(-800, 800),
                     breaks = c(-800, -400, 0, 400, 800)) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Residuals from Kd fit by nlsLM of ∆F vs. [pep13] at λ",
       color = "Wavelength (nm)",
       x = "[pep13] (µM)", y = "Residuals") +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )

ggsave("20180209_dN_hFis1_pep13_TRP_Kd_changeF_residuals_330-335.png", 
       width = 20, height = 6, dpi = 300, units = "cm")
suppressWarnings(ggsave("20180209_dN_hFis1_pep13_TRP_Kd_changeF_residuals_330-335.pdf", 
       width = 20, height = 6, dpi = 300, units = "cm"))
