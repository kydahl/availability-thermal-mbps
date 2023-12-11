# Anopheles gambiae lifespan justification
library(tidyverse)

# Get appropriate data
data <- read_rds("data/clean/data_for_TPC_fitting.rds") %>% 
  filter(system_ID %in% c(
    "Anopheles gambiae / Plasmodium falciparum",
    "Anopheles gambiae / Plasmodium spp.",
    "Anopheles gambiae / none",
    "Anopheles spp. / Plasmodium falciparum",
    "Anopheles spp. / Plasmodium spp.",
    "Anopheles spp. / none"),
    trait.name %in% c("lf")
  ) %>% 
  arrange(T)

# Set up nonlinear fits
y = data$trait
x = data$T

# Get fits
Briere_fit <- nls(y ~ c * x * (x - T0) * sqrt(Tm - x),
                  start = list(c = 1, T0 = 0, Tm = 50))
Quad_fit <- nls(y ~ -c * (x - T0) * (x - Tm),
                start = list(c = 1, T0 = 0, Tm = 50))
z = 1/y
mu_Quad_fit <- nls(z ~ c * x^2 - r * x + t,
                   start = list(c = 1, r = 1, t = 1))

# Collect fits
fit_table <- tibble(T = seq(0,50, by = 0.1), 
                    Briere = predict(Briere_fit, newdata = list(x = seq(0,50, by = 0.1))), 
                    Quadratic = predict(Quad_fit, newdata = list(x = seq(0,50, by = 0.1))), 
                    mu_Quadratic = 1/predict(mu_Quad_fit, newdata = list(x = seq(0,50, by = 0.1)))
                    ) %>% 
  pivot_longer(cols = c("Briere", "Quadratic", "mu_Quadratic"), names_to = "model") %>% 
  group_by(T) %>% 
  mutate(abs_change = (value[model == "mu_Quadratic"] - value),
         perc_change = abs_change / value[model == "mu_Quadratic"]
         )

# Determine where V0 is greater than zero
temp_lims <- read_rds("results/medianVec_vals.rds") %>% 
  filter(system_ID ==  "Anopheles gambiae / Plasmodium falciparum",
         V0 > 0
         ) %>% 
  summarise(min_Temp = min(Temperature),
            max_Temp = max(Temperature))
  

trait_plots <- fit_table %>% 
  ggplot(aes(x = T)) +
  geom_point(data = data, aes(y = trait)) +
  geom_line(aes(y = value, color = model), lwd = 1) +
  geom_vline(xintercept = temp_lims$min_Temp) +
  geom_vline(xintercept = temp_lims$max_Temp) +
  scale_y_continuous("Lifespan", expand = c(0,0)) +
  scale_x_continuous("Temperature") +
  coord_cartesian(ylim = c(0,20)) +
  scale_color_discrete(name = "Nonlinear model") +
  theme_minimal(16)

# Change from largest value
perc_change_plots <- fit_table %>% 
  filter(model != "mu_Quadratic") %>% 
  ggplot(aes(x = T)) +
  geom_line(aes(y = perc_change, color = model), lwd = 1) +
  scale_y_continuous("Percent change from mu_Quadratic model",
                     breaks = seq(-0.3, 0.3, by = 0.1),
                     labels = unique(c(seq(-0.3, 0, by = 0.1), seq(0, 0.3, by = 0.1))),
                     limits = c(-0.25, 0.25),
                     expand = c(0,0)) +
  scale_x_continuous("Temperature") +
  coord_cartesian(xlim = c(temp_lims$min_Temp, temp_lims$max_Temp)) +
  scale_color_discrete(name = "Nonlinear model") +
  theme_minimal(16)

abs_change_plots <- fit_table %>% 
  filter(model != "mu_Quadratic") %>% 
  ggplot(aes(x = T)) +
  geom_line(aes(y = abs_change, color = model), lwd = 1) +
  scale_y_continuous("Absolute change from mu_Quadratic model",
                     # breaks = seq(-0.3, 0.3, by = 0.1),
                     # labels = unique(c(seq(-0.3, 0, by = 0.1), seq(0, 0.3, by = 0.1))),
                     limits = c(-2.5, 2.5),
                     expand = c(0,0)) +
  scale_x_continuous("Temperature") +
  coord_cartesian(xlim = c(temp_lims$min_Temp, temp_lims$max_Temp)) +
  scale_color_discrete(name = "Nonlinear model")