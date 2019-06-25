# 
# Data analysis of 20171130 DLS data of 50 µM hFis1 +/- 1 mM fragment
# fragments are one of the 6 bad-acting fragments that induced peak broadening
# from zenobia 1 library during FBS by NMR screen
# Micanazole (MIC) and Nicardipine (NIC) were included as positive controls
# for colloidal aggregators. However, 0.45 um syringe filtering may have removed 
# aggregation formation of MIC and NIC. Prior to filtering MIC and NIC, 
# samples were cloudly and clearly aggregated, where solution was clear and 
# colorless after filtration.
#

library(tidyverse)
library(readxl)
library(scales)
library(stringr)
library(directlabels)

# Import and tidy-up DLS data   ------------------------------------------------

raw_data <- read_xls("20171130_hFis1_fragments_DLS_Intensity_by_size.xls")

tmp <- raw_data %>%
  select(., `Sample Name`, Sizes, Intensities) %>%
  slice(., 2:19) %>%
  rename(., sample = `Sample Name`)

DLS_data <- tmp %>%
  mutate(., V1 = strsplit(as.character(Sizes), " "),
            V2 = strsplit(as.character(Intensities), " ")) %>%
  unnest(V1, V2) %>%
  mutate(., V1 = as.numeric(V1),
         V2 = as.numeric(V2)) %>%
  select(., sample, V1, V2) %>%
  rename(., size = V1, intensity = V2)


# Plotting all DLS data: Intensity Particle Size Distribution plots   ----------

DLS_data %>%
  ggplot(., aes(x = size, y = intensity, color = sample)) +
  geom_line(aes(size = ifelse(sample == "apo_Fis1", 1, 0.5))) +
  scale_x_log10(limits = c(0.01, 10000),
                breaks = c(0.01, 0.1, 1, 10, 100, 1000, 10000),
                labels = comma) +
  scale_y_continuous(limits = c(0, 10),
                     breaks = c(0, 2.5, 5.0, 7.5, 10.0)) +
  guides(size = F) +
  scale_size_continuous(range = c(0.5, 2)) +
  labs(title = "Intensity Particle Size Distribution",
       x = "Size (radius, nm)",
       y = "Intensity (%)") +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=1, linetype='solid'),
        legend.title = element_text(size = 16),
        legend.justification = c(1, 1), legend.position = c(1, 1),
        legend.text = element_text(size = 14, color = "black")
  )

ggsave("20171130_hFis1_fragments_DLS_all.png",
       width = 40, height = 20, dpi = 300, units = "cm")

# Plotting all data containing Fis1   ------------------------------------------  
DLS_data %>%
  filter(., str_detect(sample, "Fis1")) %>% # only use samples containing Fis1
  ggplot(., aes(x = size, y = intensity, color = sample)) +
  geom_line(aes(size = ifelse(sample == "apo_Fis1", 1, 0.5))) +
  scale_x_log10(limits = c(0.01, 10000),
                breaks = c(0.01, 0.1, 1, 10, 100, 1000, 10000),
                labels = comma) +
  scale_y_continuous(limits = c(0, 10),
                     breaks = c(0, 2.5, 5.0, 7.5, 10.0)) +
  guides(size = F) +
  scale_size_continuous(range = c(0.5, 2)) +
  labs(title = "Intensity Particle Size Distribution",
       x = "Size (radius, nm)",
       y = "Intensity (%)") +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=1, linetype='solid'),
        legend.title = element_text(size = 16),
        legend.justification = c(1, 1), legend.position = c(1, 1),
        legend.text = element_text(size = 14, color = "black")
  )

ggsave("20171130_hFis1_fragments_DLS_Fis1_only.png",
       width = 40, height = 20, dpi = 300, units = "cm")

# Plotting Fis1 +/- badacting fragments     ------------------------------------
DLS_data %>% 
  filter(., str_detect(sample, c("buffer", "Fis1"))) %>%
  filter(., sample != "Fis1_ZT0327_2",
         sample != "Fis1_ZT0327_0.01percentTRITONX100") %>%
  ggplot(., aes(x = size, y = intensity, color = sample)) +
  geom_line(aes(size = ifelse(sample == "apo_Fis1", 1, 0.5))) +
  scale_x_log10(limits = c(0.01, 10000),
                breaks = c(0.01, 0.1, 1, 10, 100, 1000, 10000),
                labels = comma) +
  scale_y_continuous(limits = c(0, 10),
                     breaks = c(0, 2.5, 5.0, 7.5, 10.0)) +
  guides(size = F) +
  scale_size_continuous(range = c(0.5, 2)) +
  labs(title = "Intensity Particle Size Distribution",
       x = "Size (radius, nm)",
       y = "Intensity (%)") +
  scale_color_manual("Fis1 + Fragment",
                     values = c("buffer" = "#999999",
                                "apo_Fis1" = "black",
                                "Fis1_ZT0327" = "#e41a1c",
                                "Fis1_ZT0403" = "#377eb8",
                                "Fis1_ZT0405" = "#4daf4a",
                                "Fis1_ZT0406" = "#984ea3",
                                "Fis1_ZT0808" = "#ff7f00",
                                "Fis1_ZT0818" = "#f781bf"),
                     labels = c("Buffer", "Apo", "ZT0327", "ZT0403", "ZT0405", 
                                "ZT0406", "ZT0808", "ZT0818")) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=1, linetype='solid'),
        legend.title = element_text(size = 16),
        legend.justification = c(1, 1), legend.position = c(1, 1),
        legend.text = element_text(size = 14, color = "black")
  )

ggsave("20171130_hFis1_fragments_DLS_Fis1_fragment_final.png",
       width = 40, height = 20, dpi = 300, units = "cm")

# Plotting only fragments, including postive control MIC and NIC   -------------
DLS_data %>% 
  filter(., str_detect(sample, c("ZT", "buffer", "MIC", "NIC"))) %>%
  filter(., !grepl("Fis1", sample)) %>%
  ggplot(., aes(x = size, y = intensity, color = sample)) +
  geom_line(aes(size = ifelse(sample == "buffer", 1, 0.5))) +
  scale_x_log10(limits = c(0.01, 10000),
                breaks = c(0.01, 0.1, 1, 10, 100, 1000, 10000),
                labels = comma) +
  scale_y_continuous(limits = c(0, 10),
                     breaks = c(0, 2.5, 5.0, 7.5, 10.0)) +
  guides(size = F) +
  scale_size_continuous(range = c(0.5, 2)) +
  scale_color_brewer("Fragment", palette = "Set1") +
  labs(title = "Intensity Particle Size Distribution",
       x = "Size (radius, nm)",
       y = "Intensity (%)") +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=1, linetype='solid'),
        legend.title = element_text(size = 16),
        legend.justification = c(1, 1), legend.position = c(1, 1),
        legend.text = element_text(size = 14, color = "black")
  )

ggsave("20171130_hFis1_fragments_DLS_fragment_only.png",
       width = 40, height = 20, dpi = 300, units = "cm")

# Stacked DLS plot of Fis1 +/- badacting fragments for ChemBioChem paper  ------

DLS_data %>% 
  filter(., str_detect(sample, c("buffer", "Fis1"))) %>%
  filter(., sample != "Fis1_ZT0327_2",
         sample != "Fis1_ZT0327_0.01percentTRITONX100") %>%
  ggplot(., aes(x = size, y = intensity, height = intensity, color = sample)) +
  geom_line(aes(size = ifelse(sample == "apo_Fis1", 1, 0.5))) +
  facet_grid(sample ~ .) +
  scale_x_log10(limits = c(0.01, 10000),
                breaks = c(0.01, 0.1, 1, 10, 100, 1000, 10000),
                labels = comma) +
  scale_y_continuous(limits = c(0, 10),
                     breaks = c(0, 4.0, 8.0)) +
  guides(size = F) +
  scale_size_continuous(range = c(0.5, 2)) +
  labs(title = "Intensity Particle Size Distribution",
       x = "Size (radius, nm)",
       y = "Intensity (%)") +
  scale_color_manual("Fis1 + Fragment",
                     values = c("buffer" = "#999999",
                                "apo_Fis1" = "black",
                                "Fis1_ZT0327" = "#e41a1c",
                                "Fis1_ZT0403" = "#377eb8",
                                "Fis1_ZT0405" = "#4daf4a",
                                "Fis1_ZT0406" = "#984ea3",
                                "Fis1_ZT0808" = "#ff7f00",
                                "Fis1_ZT0818" = "#f781bf"),
                     labels = c("Buffer", "Apo", "ZT0327", "ZT0403", "ZT0405", 
                                "ZT0406", "ZT0808", "ZT0818")) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=1, linetype='solid'),
        legend.title = element_text(size = 16),
        #legend.justification = c(1, 1), legend.position = c(1, 1),
        legend.text = element_text(size = 14, color = "black"),
        strip.background = element_blank(),
        strip.text = element_blank()
  )

ggsave("20171130_hFis1_fragments_DLS_Fis1_fragment_stacked.png",
       width = 40, height = 20, dpi = 300, units = "cm")


# Stacked DLS plot of all data for ChemBioChem paper  ------

DLS_data %>%
  filter(., sample != "MIC", sample != "NIC",
         sample != "Fis1_ZT0327_2",
         sample != "Fis1_ZT0327_0.01percentTRITONX100") %>%
  ggplot(., aes(x = size, y = intensity, height = intensity, color = sample)) +
  geom_line() +
  facet_grid(sample ~ .) +
  geom_dl(aes(label = sample), 
          method = list(dl.trans(x = x - 2.5, y = y + 0.5), dl.combine("last.bumpup"))) +
  scale_x_log10(limits = c(0.01, 10000),
                breaks = c(0.1, 10, 1000),
                labels = comma) +
  scale_y_continuous(limits = c(0, 10),
                     breaks = c()) +
  guides(size = F, color = F) +
  scale_size_continuous(range = c(0.5, 2)) +
  labs(title = "",
       x = "Size (radius, nm)",
       y = "Intensity") +
  scale_color_manual("Fis1 + Fragment",
                     values = c("buffer" = "#999999",
                                "apo_Fis1" = "black",
                                "Fis1_ZT0327" = "#e41a1c",
                                "Fis1_ZT0403" = "#377eb8",
                                "Fis1_ZT0405" = "#4daf4a",
                                "Fis1_ZT0406" = "#984ea3",
                                "Fis1_ZT0808" = "#ff7f00",
                                "Fis1_ZT0818" = "#f781bf",
                                "ZT0327" = "#e41a1c",
                                "ZT0403" = "#377eb8",
                                "ZT0405" = "#4daf4a",
                                "ZT0406" = "#984ea3",
                                "ZT0808" = "#ff7f00",
                                "ZT0818" = "#f781bf")) +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        legend.title = element_text(size = 16),
        strip.background = element_blank(),
        strip.text = element_blank()
  )

ggsave("20171130_hFis1_fragments_DLS_all_stacked.png",
       width = 12.5, height = 18, dpi = 300, units = "cm")

