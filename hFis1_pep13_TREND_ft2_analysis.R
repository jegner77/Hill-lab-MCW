## # A tibble: 1 x 8
##    sigma isConv       finTol logLik   AIC   BIC deviance df.residual
##    <dbl> <lgl>         <dbl>  <dbl> <dbl> <dbl>    <dbl>       <int>
## 1 0.0299 TRUE   0.0000000149   17.9 -29.8 -29.5  0.00536           6# 
# Kd determination of pep13 against hFis1 (aa 1-125) using TREND on ft2 files
#

library(tidyverse)
library(broom)
library(minpack.lm)
library(gridExtra)

# Import normalized PC1 values from using TREND   ------------------------------

norm_pc1 <- as_tibble(list(conc = c(0, 25, 50, 150, 400, 800, 1600, 2000),
                           pc1 = c(0, 0.33, 0.44, 0.78, 0.90, 0.99, 0.99, 1.00)))

# Plot and fit normalized PC1 data for Kd determination   ----------------------

# using ligand depletion model (model 1):
a <- norm_pc1 %>%
  ggplot(., aes(x = conc, y = pc1)) +
  geom_point(shape = 15) +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ dmax*((kd+50+x) - ((kd+50+x)^2 - 4*50*x)^0.5)/(2*50),
            method.args = list(start = c(dmax = 1,
                                         kd = 700),
                               control = nls.control(maxiter = 100, tol = 1e-6)),
            se = FALSE,        # [protein] held constant at 50 µM
            fullrange = TRUE,  
            size = 1,
            alpha = 0.6) +
  theme_bw() +
  scale_x_continuous(limits = c(0, 2000),
                     breaks = c(0, 400, 800, 1200, 1600, 2000)) +
  scale_y_continuous(limits = c(0, 1.1),
                     breaks = c(0, 0.25, 0.50, 0.75, 1.0)) +
  labs(x = "[pep13] (µM)",
       y = "Normalized PC1") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        axis.line = element_line(color = "black", size = 1, linetype = 1),
        axis.ticks  = element_line(color = "black", size = 1, linetype = 1),
        axis.ticks.length = unit(0.2, "cm")
  )
a

# calculate Kd value from model 1
kd_model1 <- norm_pc1 %>%
  do(tidy(nlsLM(
    formula = pc1 ~ dmax*((kd+50+conc) - ((kd+50+conc)^2 - 4*50*conc)^0.5)/(2*50),
    start = list(dmax = 1, kd = 700), 
    trace = TRUE,
    control = nls.control(maxiter = 100, tol = 1e-6),
    data = .))) %>%
  spread(., key = term, value = estimate)

kd_model1

kd_model1_resid <- norm_pc1 %>%
  do(augment(nlsLM(
    formula = pc1 ~ dmax*((kd+50+conc) - ((kd+50+conc)^2 - 4*50*conc)^0.5)/(2*50),
    start = list(dmax = 1, kd = 200), 
    trace = TRUE,
    control = nls.control(maxiter = 100, tol = 1e-6),
    data = .)))

kd_model1_glance <- norm_pc1 %>%
  do(glance(nlsLM(
    formula = pc1 ~ dmax*((kd+50+conc) - ((kd+50+conc)^2 - 4*50*conc)^0.5)/(2*50),
    start = list(dmax = 1, kd = 200), 
    trace = TRUE,
    control = nls.control(maxiter = 100, tol = 1e-6),
    data = .)))
kd_model1_glance

b <- kd_model1_resid %>%
  ggplot(., aes(x = conc, y = .resid)) +
  geom_point(shape = 15) +
  geom_hline(yintercept = 0, color = "red", linetype = 2) +
  scale_x_continuous(limits = c(0, 2000),
                     breaks = c(0, 400, 800, 1200, 1600, 2000)) +
  scale_y_continuous(limits = c(-0.11, 0.11),
                     breaks = c(-0.1, 0, 0.1)) +
  labs(x = "",
       y = "Residuals") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        axis.line = element_line(color = "black", size = 1, linetype = 1),
        axis.ticks  = element_line(color = "black", size = 1, linetype = 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text.x = element_blank()
  )
b

# Plotting regression with residuals, model 1 ----------------------------------

# Generate layout matrix to plot residuals at 1/3 size of regression plot
lay <- cbind(c(1, 2, 2))

plot1 <- grid.arrange(grobs = list(b, a), layout_matrix = lay,
                      top = "Model 1: hFis1 + pep13")

ggsave("TREND_model1_hFis1_pep13.pdf", plot1,
       width = 10, height = 10, units = "cm")
ggsave("TREND_model1_hFis1_pep13.png", plot1,
       width = 10, height = 10, units = "cm")

# Using simplified Kd model (model 2): -----------------------------------------
c <- norm_pc1 %>%
  ggplot(., aes(x = conc, y = pc1)) +
  geom_point(shape = 15) +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ (dmax*x) / (kd + x),
            method.args = list(start = c(dmax = 1,
                                         kd = 700),
                               control = nls.control(maxiter = 100, tol = 1e-6)),
            se = FALSE,        
            fullrange = TRUE,  
            size = 1,
            alpha = 0.6) +
  scale_x_continuous(limits = c(0, 2000),
                     breaks = c(0, 400, 800, 1200, 1600, 2000)) +
  scale_y_continuous(limits = c(0, 1.1),
                     breaks = c(0, 0.25, 0.50, 0.75, 1.0)) +
  labs(x = "[pep13] (µM)",
       y = "Normalized PC1") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        axis.line = element_line(color = "black", size = 1, linetype = 1),
        axis.ticks  = element_line(color = "black", size = 1, linetype = 1),
        axis.ticks.length = unit(0.2, "cm")
  )
c

# calculate Kd value from model 2
kd_model2 <- norm_pc1 %>%
  do(tidy(nlsLM(
    formula = pc1 ~ (dmax*conc) / (kd + conc),
    start = list(dmax = 1, kd = 700), 
    trace = TRUE,
    control = nls.control(maxiter = 100, tol = 1e-6),
    data = .))) %>%
  spread(., key = term, value = estimate)

kd_model2

kd_model2_resid <- norm_pc1 %>%
  do(augment(nlsLM(
    formula = pc1 ~ (dmax*conc) / (kd + conc),
    start = list(dmax = 1, kd = 700), 
    trace = TRUE,
    control = nls.control(maxiter = 100, tol = 1e-6),
    data = .)))

kd_model2_glance <- norm_pc1 %>%
  do(glance(nlsLM(
    formula = pc1 ~ (dmax*conc) / (kd + conc),
    start = list(dmax = 1, kd = 700), 
    trace = TRUE,
    control = nls.control(maxiter = 100, tol = 1e-6),
    data = .)))

kd_model2_glance

d <- kd_model2_resid %>%
  ggplot(., aes(x = conc, y = .resid)) +
  geom_point(shape = 15) +
  geom_hline(yintercept = 0, color = "red", linetype = 2) +
  scale_x_continuous(limits = c(0, 2000),
                     breaks = c(0, 400, 800, 1200, 1600, 2000)) +
  scale_y_continuous(limits = c(-0.11, 0.11),
                     breaks = c(-0.1, 0, 0.1)) +
  labs(x = "",
       y = "Residuals") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        axis.line = element_line(color = "black", size = 1, linetype = 1),
        axis.ticks  = element_line(color = "black", size = 1, linetype = 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text.x = element_blank()
  )
d

# Plotting regression with residuals, model 2 ----------------------------------
plot2 <- grid.arrange(grobs = list(d, c), layout_matrix = lay,
                      top = "Model 2: hFis1 + pep13")

ggsave("TREND_model2_hFis1_pep13.pdf", plot2,
       width = 10, height = 10, units = "cm")
ggsave("TREND_model2_hFis1_pep13.png", plot2,
       width = 10, height = 10, units = "cm")


# Model discrimination using AIC, BIC, log likelihood, & deviance ---------

model_dis <- kd_model1_glance %>%
  union(kd_model2_glance) %>%
  mutate(., model = c("lig depletion", "Langmuir")) %>%
  select(model, everything())
print(model_dis, width = Inf)
