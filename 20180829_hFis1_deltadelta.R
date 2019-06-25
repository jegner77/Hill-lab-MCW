# 
# Delta-delta plots for fingerprinting and CSP mapping of
# FBS by NMR endpoint of HSQC titration experiments against 15N-hFis1
# Kd values previously determined by TREND fitting
#

# Load libraries, import data, and tidy data    --------------------------------

library(tidyverse)
library(broom)
library(minpack.lm)
library(readxl)

raw <- read_excel("20180829_hFis1_titration_endpoint_deltadelta.xlsx",
                  sheet = 1) %>%
  gather(., "tmp1", "cs", 2:25) %>%
  separate(., col = "tmp1", into = c("cmpd", "type"), sep = "-") %>%
  mutate(., residue = parse_number(label)) %>%
  select(., residue, everything()) %>%
  spread(., key = "type", value = "cs") %>%
  group_by(., cmpd) %>%
  nest()
  
# Calculate total chemical shift perturbation and plot against residue #  ------

deltadelta <- raw %>%
  unnest() %>%
  mutate(., CSP = sqrt(((5*(apo_H - end_H))^2) + (apo_N - end_N)^2))

mean_sd <- deltadelta %>%
  group_by(., cmpd) %>%
  summarise(., mean = mean(CSP),
            sd = sd(CSP)) %>%
  mutate(twoSD = 2*sd, sigma = mean + sd, twosigma = sigma*2)
mean_sd

write_csv(deltadelta, "hFis1_by_cmpd_CSP_sigmas.csv")

deltadelta <- deltadelta %>%
  left_join(., mean_sd, by = "cmpd")

# CSP by Residue # for apo and endpoint of cmpd concentration
deltadelta %>% 
  filter(., residue != "W40") %>%
  ggplot(., aes(x = residue, y = CSP)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ cmpd) +
  theme_bw()

# Plotting with purrr using split()
plots <- deltadelta %>%
  filter(., residue != "W40") %>%
  split(.$cmpd) %>%
  map(~ggplot(data = .x, aes(x = residue, y = CSP)) +
                       geom_bar(stat = "identity") + 
                       geom_hline(aes(yintercept = sd),
                                     color = "green") +
                       geom_hline(aes(yintercept = twoSD),
                                  color = "purple") +
                       scale_x_continuous(limits = c(0, 125),
                                          breaks = c(0, 20, 40, 60,
                                                     80, 100, 120)) +
                       scale_y_continuous(breaks = c(0, 0.5, 1, 1.5,
                                                     2, 2.5, 3)) +
                       labs(title = .$cmpd,
                            x = "Residue",
                            y = "CSP (ppm)") + 
                       theme_classic() +
                       theme(axis.text = element_text(color = "black",
                                                      size = 12),
                             axis.title = element_text(color = "black",
                                                       size = 12),
                             axis.line = element_line(color = "black",
                                                      size = 1, linetype = 1),
                             axis.ticks  = element_line(color = "black",
                                                        size = 1, linetype = 1),
                             axis.ticks.length = unit(0.2, "cm")
                             ))

# Print plots prior to saving; verify above code functions properly
walk(plots, print)

# Save plots
walk2(.x = plots, .y = names(plots), 
      ~ ggsave(filename = paste0(.y, "_CSP_plot.pdf"), plot = .x,
               width = 10, height = 8, units = "cm"))

# Plotting with purrr using nest()
plots2 <- deltadelta %>%
  filter(., residue != "W40") %>%
  group_by(., cmpd) %>%
  nest() %>%
  mutate(plot = map2(.x = data, .y = cmpd, ~ggplot(data = .x, aes(x = residue, y = CSP)) +
        geom_bar(stat = "identity") + 
        geom_hline(aes(yintercept = sd),
                   color = "green") +
        geom_hline(aes(yintercept = twoSD),
                   color = "purple") +
        scale_x_continuous(limits = c(0, 125),
                           breaks = c(0, 20, 40, 60,
                                      80, 100, 120)) +
        scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
        labs(title = .y,
             x = "Residue",
             y = "CSP (ppm)") + 
        theme_classic() +
        theme(axis.text = element_text(color = "black",
                                       size = 12),
              axis.title = element_text(color = "black",
                                        size = 12),
              axis.line = element_line(color = "black",
                                       size = 1, linetype = 1),
              axis.ticks  = element_line(color = "black",
                                         size = 1, linetype = 1),
              axis.ticks.length = unit(0.2, "cm")
        )))

# Print plots prior to saving; verify above code functions properly
walk(plots2, print)

# Save plots
map2(.x = paste0(plots2$cmpd, "_CSP_plot2.pdf"), .y = plots2$plot, ggsave)
    ## this works, but I am unable to control size of generated plot

# Write out CSV file for CSP mapping in PyMol   --------------------------------

CSP_mapping <- deltadelta %>%
  split(.$cmpd) %>%
  map(~mutate(.x, CSP_1sigma = CSP > sd,
      CSP_2sigma = CSP > twoSD)) %>%
  map(~select(.x, cmpd, residue, label, CSP, sd, twoSD, CSP_1sigma, CSP_2sigma))

walk2(.x = CSP_mapping, .y = names(CSP_mapping), 
     ~write_csv(.x, paste0(.y, "_CSP_mapping.csv")))
