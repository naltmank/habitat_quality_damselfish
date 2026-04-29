rm(list=ls())
# install.packages("librarian")
librarian::shelf(here, tidyverse, effects, ggplot2, ggpubr, glmmTMB, DHARMa, performance, emmeans, car)

#### NUTRIENTS ####
nutrients <- read.csv(here::here("data", "lter_nutrient_data.csv"))

nutrients_sub <- nutrients %>%
  filter(site == "LTER_1") %>%
  filter(habitat != "Reefcrest") %>%
  filter (year <= 2022) %>%
  mutate(habitat_id = paste0(year, "-", habitat),
         year = as.factor(year)) 



# sargassum wasn't sampled consistently the way turb has been - focus on turbinaria for model
turbinaria <- subset(nutrients_sub, genus == "Turbinaria")
# missing data in 2006 - remove
turbinaria <- na.omit(turbinaria)

turb_cn_mod <- glmmTMB(CN_ratio ~ habitat + (1|year) , data = turbinaria)
plot(simulateResiduals(turb_cn_mod))
summary(turb_cn_mod)
Anova(turb_cn_mod)
plot(allEffects(turb_cn_mod))
performance::r2(turb_cn_mod)

##### PLOT #####
turbinaria_plot_df <- turbinaria %>%
  group_by(year, habitat) %>%
  summarize(cn_mean = mean(CN_ratio, na.rm = T),
            cn_se = sd(CN_ratio, na.rm = T)/sqrt(n()),
            group = paste0(habitat, "-", year))


(turbinaria_plot <- ggplot() +
    geom_point(data = turbinaria_plot_df, aes(x = year, y = cn_mean, colour = habitat),
               position = position_dodge(0.4), size = 3) +
    geom_errorbar(data = turbinaria_plot_df, aes(x = year, ymin = cn_mean - cn_se,
                                           ymax = cn_mean + cn_se,
                                           colour = habitat),
                  position = position_dodge(0.4)) +
    geom_line(data = turbinaria_plot_df, aes(x = year, y = cn_mean, colour = habitat, group = habitat),
              position = position_dodge(0.4)) +
    annotate("text", x = 1, y = 52,
             label = "Habitat - P < 0.0001",
             hjust = 0, size = 8) +
    labs(title = "a.\n", x = "Year", y = "C:N\n") +
    scale_colour_discrete(name = "Habitat") +
    theme_classic() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_text(colour = "black", size = 24),
          legend.text = element_text(colour = "black", size = 20),
          legend.position = "bottom",
          axis.line = element_line(colour = "black", size = 1.2),
          plot.title = element_text(color = "black", size = 25, hjust = 0, vjust = 0, face = "plain"),
          axis.text.x = element_text(color = "black", size = 24, angle = 45, hjust = .5, vjust = .5, face = "plain"),
          axis.title.x = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain"),
          axis.text.y = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain"),
          axis.title.y = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain", margin = margin(r = 15))
    )
)

#### SEDIMENTATION ####
sediment <- read.csv(here::here("data", "sedimentation_rates.csv"))

# some territories missing one of the pairs
missing_territories <- c("T3", "T31", "T32", "T40", "T45", "T5")

sediment_sub <- sediment[-which(sediment$territory_no %in% missing_territories),]

hist(sediment$sedimentation_rate) # gamma distribution
hist(sediment_sub$sedimentation_rate)

# alt - just the outside territories
outside_sed <- subset(sediment, territory_position == "Outside")

sed_mod <- glmmTMB(sedimentation_rate ~ territory_position * reef_location + (1|pair_id),
                   family = Gamma("log"),
                   data = sediment_sub)
summary(sed_mod)
Anova(sed_mod)
plot(allEffects(sed_mod))
plot(simulateResiduals(sed_mod)) # no issues
performance::r2(sed_mod)

# pairwise contrasts
sed.emm <- emmeans(sed_mod, ~ territory_position * reef_location)
contrast(sed.emm, "consec", simple = "each", combine = TRUE, adjust = "mvt")
sed_letters <- data.frame(
  reef_location = c("Backreef", "Backreef", "Fringing", "Fringing"),
  territory_position = c("Outside", "Inside", "Outside", "Inside"),
  letters = c("b", "a", "c", "c"),
  heights = rep(175, 4)
)


# just outside
outside_mod <- glmmTMB(sedimentation_rate ~ reef_location + (1|pair_id),
                       family = Gamma("log"),
                       data = outside_sed)
summary(outside_mod)
Anova(outside_mod)
plot(simulateResiduals(outside_mod))
performance::r2(outside_mod)

##### PLOT #####
(sediment_plot <- ggplot() +
  geom_boxplot(data = sediment_sub, aes(x = territory_position, y = sedimentation_rate,
                                        colour = reef_location), outlier.shape = NA) +
  geom_point(data = sediment_sub, aes(x = territory_position, y = sedimentation_rate,
                                      colour = reef_location), position = position_jitterdodge(0.4)) +
  geom_text(data = sed_letters, aes(x = territory_position, y = heights, group = reef_location,
                                    label = letters),
            position = position_dodge(0.75), size = 8) +
  scale_y_log10() +
  theme_classic() +
  labs(title = "b.\n", x = "Territory position",
       y = expression(atop("Sedimentation rate", (g ~ m^-2 ~ d^-1)))) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black", size = 1.2),
        plot.title = element_text(color = "black", size = 25, hjust = 0, vjust = 0, face = "plain"),
        axis.text.x = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain"),
        axis.title.x = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain", margin = margin(r = 15))
        )
)

##### PANEL GRAPH #####
(water_quality <- ggarrange(turbinaria_plot, sediment_plot, nrow = 2, common.legend = T, legend = "bottom",
                            heights = c(1,1)))

ggsave(filename = here::here("output", "water_quality.png"), water_quality, width = 10, height = 15,
       dpi = "retina")

#### FISHES ####
fish <- read.csv(here::here("data", "lter_fish_data_2022.csv")) %>%
  # negative values used to indicate no observations - remove
  filter(biomass > 0)

# subset piscivores
piscivore <- fish %>%
  filter(fine_trophic %in% "Piscivore") %>%
  # biomass bombs significantly affect biomass distributions - remove
  filter(!species %in% c("Carcharhinus melanopterus", "Triaenodon obesus")) %>%
  filter(lter_region == 1) %>%
  filter(habitat != "FO") %>%
  group_by(year, habitat, transect_no) %>%
  summarise(biomass_sum = sum(biomass, na.rm = T))

# subset scrapers and grazers - primary competitors with damselfishes
herbivore <- fish %>%
  filter(fine_trophic %in% c("Cropper", "Brusher", "Scraper", "Excavator", "Herbivore/Detritivore")) %>%
  
  # some damselfish caught on surveys (but not well) - remove
  filter(!species %in% c("Stegastes fasciolatus", "Stegastes punctatus")) %>%
  filter(lter_region == 1) %>%
  filter(habitat != "FO") %>%
  group_by(year, habitat, transect_no) %>%
  summarise(biomass_sum = sum(biomass, na.rm = T))

# add transect ID columns
herbivore$transect_id <- paste0(herbivore$habitat, "-", herbivore$transect_no)
piscivore$transect_id <- paste0(piscivore$habitat, "-", piscivore$transect_no)

herbivore$unique_id <- paste0(herbivore$year, "-", herbivore$habitat, "-", herbivore$transect_no)
piscivore$unique_id <- paste0(piscivore$year, "-", piscivore$habitat, "-", piscivore$transect_no)

# identify missing piscivore transects
missing_piscivores_df <- herbivore %>%
  anti_join(piscivore, by = c("year", "habitat", "transect_no", "transect_id")) %>%
  mutate(biomass_sum = 0)

# combine with existing piscivore data
piscivore_complete <- piscivore %>%
  bind_rows(missing_piscivores_df)

# add guild columns
piscivore_complete$guild <- rep("piscivore", nrow(piscivore_complete))
herbivore$guild <- rep("herbivore", nrow(herbivore))



##### PISCIVORE MODEL #####
hist(log(piscivore_complete$biomass_sum + 1)) # likely zero-inflated gamme works best
piscivore_complete$year <- as.factor(piscivore_complete$year)
piscivore_mod <- glmmTMB(log(biomass_sum + 1) ~ habitat + (1 | year) + (1|transect_id),
                         ziformula = ~ habitat, family = ziGamma("log"),
                         data = piscivore_complete)
summary(piscivore_mod) # more likely not to see piscivores in backreef
Anova(piscivore_mod)
emmeans(piscivore_mod, pairwise ~ habitat)
plot(simulateResiduals(piscivore_mod))
performance::r2(piscivore_mod)
##### HERBIVORE MODEL #####
hist(log(herbivore$biomass_sum + 1))
herbivore$year <- as.factor(herbivore$year)
herbivore_mod <- glmmTMB(log(biomass_sum + 1) ~ habitat + (1|year) + (1|transect_id),
                         family = Gamma("log"),
                         data = herbivore)
summary(herbivore_mod) # more likely not to see piscivores in fringing
Anova(herbivore_mod)
plot(simulateResiduals(herbivore_mod))
plot(allEffects(herbivore_mod))

emmeans(herbivore_mod, pairwise ~ habitat)

##### PLOT #####
relevant_fish <- rbind(herbivore, piscivore_complete)
fish_plot_df <- relevant_fish %>%
  group_by(year, habitat, guild) %>%
  summarise(biomass_mean = mean(biomass_sum, na.rm = T),
            biomass_se = sd(biomass_sum, na.rm = T)/sqrt(n())) %>%
  mutate(year = as.factor(year),
         group = paste0(habitat, "-", guild))

(biomass_plot <- ggplot() +
  geom_point(data = fish_plot_df, aes(x = year, y = biomass_mean, colour = habitat, shape = guild),
             position = position_dodge(0.4), size = 3) +
  geom_errorbar(data = fish_plot_df, aes(x = year, ymin = biomass_mean - biomass_se,
                                         ymax = biomass_mean + biomass_se,
                                         colour = habitat, group = group),
                position = position_dodge(0.4)) +
  geom_line(data = fish_plot_df, aes(x = year, y = biomass_mean, colour = habitat, group = group),
            position = position_dodge(0.4)) +
    annotate("text", x = 16.5, y = 1.5,
             label = "Herbivores - P = 0.051\nPiscivores - P < 0.001",
             hjust = 1, size = 7) +
  scale_y_log10(limits = c(1,75000)) +
  scale_colour_discrete(name = "Habitat",
                        labels = c("BA" = "Backreef",
                                   "FR" = "Fringing")) +
  scale_shape_manual(name = "Fish guild",
                     values = c("herbivore" = 19,
                                "piscivore" = 1),
                     labels = c("herbivore" = "Herbivore",
                                "piscivore" = "Piscivore")) +
  labs(title = "", x = "Year", y = expression(Biomass ~ (g ~ 250 ~ m^-2))) +
  theme_classic() +
  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(colour = "black", size = 24),
        legend.text = element_text(colour = "black", size = 20),
        axis.line = element_line(colour = "black", size = 1.2),
        plot.title = element_text(color = "black", size = 25, hjust = 0, vjust = 0, face = "plain"),
        axis.text.x = element_text(color = "black", size = 24, angle = 45, hjust = .5, vjust = .5, face = "plain"),
        axis.title.x = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain", margin = margin(r = 15))
       )
)

ggsave(filename = here::here("output", "herb_piscivore_plot.png"), biomass_plot, width = 12, height = 7.5,
       dpi = "retina")

