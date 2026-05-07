rm(list=ls())
# install.packages("librarian")
librarian::shelf(here, tidyverse, dplyr, effects, ggplot2, ggpubr, ggrepel,
                 glmmTMB, DHARMa, performance, emmeans, vegan, car)

#### READ AND FORMAT DATA ####
physiology <- read.csv(here::here("data", "damselfish_physiology.csv")) %>%
  # fish_id has formatting issues - rewrite them
  mutate(fish_id = paste0(territory_id, "-", fish_no)) %>%
  rename(habitat = reef_location)

# all relevant columns for first subset
relevant_cols <- c("territory_id", "pair_id", "habitat", "fish_id", "fertilizer_trt",
                   "sex", "intestinal_surface_area", "body_condition_factor",  "gonadosomatic_index",
                   "hepatosomatic_index")

# metadata cols for later subset
meta_cols <- c("territory_id", "pair_id", "habitat", "fish_id", "fertilizer_trt", "sex")

# subset relevant columns
physiology_sub <- physiology[ ,which(names(physiology) %in% relevant_cols)]

# split by sex for analyses
male <- subset(physiology_sub, sex == "M")
female <- subset(physiology_sub, sex == "F")

# prep for PCA on body condition indicators
male_cond <- as.matrix(male[ , which(names(male) %in% c("body_condition_factor",  "gonadosomatic_index",
                                                        "hepatosomatic_index"))])
rownames(male_cond) <- male$fish_id


female_cond <- as.matrix(female[ ,which(names(female) %in% c("body_condition_factor",  "gonadosomatic_index",
                                                             "hepatosomatic_index"))])
rownames(female_cond) <- female$fish_id

#### PCA ####
##### MALES #####
male_pca <- prcomp(male_cond, center = T, scale. = TRUE)
male_pca_points <- as.data.frame(male_pca$x)[,1:2]
# add fish_id back in as col for later merge
male_pca_points$fish_id <- rownames(male_pca_points)
male_pca_loadings <- as.data.frame(male_pca$rotation[, 1:2]) %>%
  rownames_to_column("variable") %>%
  mutate(variable = case_when(variable == "body_condition_factor" ~ "BCF" ,
                              variable == "gonadosomatic_index"   ~ "GSI" ,
                              variable == "hepatosomatic_index"   ~ "HSI" ,
                              T ~ variable))

male_pca_plotting <- merge(male, male_pca_points, by = "fish_id")

male_pca_hull <- 
  male_pca_plotting %>% 
  group_by(habitat, fertilizer_trt) %>% 
  slice(chull(PC1, PC2))

###### PLOT ######
(male_pca_fig <- ggplot() +
   geom_point(data = male_pca_plotting, aes(PC1, PC2, color = habitat, shape = fertilizer_trt), size = 3) +
   geom_polygon(data = male_pca_hull,
                aes(PC1, PC2, linetype = fertilizer_trt,
                    colour = habitat),
                fill = NA) +
   # linetype doesn't show up with just hull - cheese it by adding a path object
   geom_path(data = male_pca_hull,
             aes(PC1, PC2, linetype = fertilizer_trt,
                 colour = habitat)) +
   geom_segment(data = male_pca_loadings,
                aes(x = 0, y = 0, xend = PC1, yend = PC2),
                arrow = arrow(length = unit(0.2, "cm")), color = "black") +
   geom_text_repel(data = male_pca_loadings,
                   aes(x = PC1*1.2 , y = PC2*1.2, label = variable),
                   size = 7, hjust = 0.5, max.overlaps = 30) +

   scale_colour_discrete(name = "Habitat") +
   scale_shape_manual(name = "Treatment",
                      values = c("Y" = 19,
                                 "N" = 1),
                      labels = c("Y" = "Fertilized",
                                 "N" = "Unfertilized")) +
   scale_linetype_manual(name = "Treatment",
                         values = c("Y" = 1,
                                    "N" = 2),
                         labels = c("Y" = "Fertilized",
                                    "N" = "Unfertilized")) +
   labs(title = "a.", x = "PC1 (39.8%)", y = "PC2 (33.5%)") +
   theme_classic() +
   theme(axis.title.x = element_text(size=25), 
         axis.title.y = element_text(size=25), 
         axis.text.x = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain", margin = margin(r = 15)),
         axis.text.y = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain", margin = margin(r = 15)),
         strip.text.x = element_text(size = 25),
         plot.title = element_text(color = "black", size = 25, hjust = 0, vjust = 0, face = "plain"),
         panel.grid.major = element_blank(),  #remove major-grid labels
         panel.grid.minor = element_blank(),  #remove minor-grid labels
         legend.title = element_text(colour = "black", size = 24),
         legend.text = element_text(colour = "black", size = 20),
         legend.position = 'bottom', # move legend to bottom
         legend.box.background = element_rect(colour = "black")) +
   annotate("text", x = -3.15, y = -2.1,
            label = "Habitat - P = 0.001\nTreatment - P = 0.42\nInteraction - P = 0.79",
            hjust = 0, size = 7) 
)

##### FEMALES #####
female_pca <- prcomp(female_cond, center = T, scale. = TRUE)
female_pca_points <- as.data.frame(female_pca$x)[,1:2]
# add fish_id back in as col for later merge
female_pca_points$fish_id <- rownames(female_pca_points)
female_pca_loadings <- as.data.frame(female_pca$rotation[, 1:2]) %>%
  rownames_to_column("variable") %>%
  mutate(variable = case_when(variable == "body_condition_factor" ~ "BCF" ,
                              variable == "gonadosomatic_index"   ~ "GSI" ,
                              variable == "hepatosomatic_index"   ~ "HSI" ,
                              T ~ variable))

female_pca_plotting <- merge(female, female_pca_points, by = "fish_id")

female_pca_hull <- 
  female_pca_plotting %>% 
  group_by(habitat, fertilizer_trt) %>% 
  slice(chull(PC1, PC2))

###### PLOT ######
(female_pca_fig <- ggplot() +
   geom_point(data = female_pca_plotting, aes(PC1, PC2, color = habitat, shape = fertilizer_trt), size = 3) +
   geom_polygon(data = female_pca_hull,
                aes(PC1, PC2, linetype = fertilizer_trt,
                    colour = habitat),
                alpha = 0) +
   geom_path(data = female_pca_hull,
             aes(PC1, PC2, linetype = fertilizer_trt,
                 colour = habitat)) +
   geom_segment(data = female_pca_loadings,
                aes(x = 0, y = 0, xend = PC1, yend = PC2),
                arrow = arrow(length = unit(0.2, "cm")), color = "black") +
   geom_text_repel(data = female_pca_loadings,
                   aes(x = PC1*1.2 , y = PC2*1.2, label = variable),
                   size = 6, hjust = 0.5, max.overlaps = 30) +
   scale_colour_discrete(name = "Habitat") +
   scale_shape_manual(name = "Treatment",
                      values = c("Y" = 19,
                                 "N" = 1),
                      labels = c("Y" = "Fertilized",
                                 "N" = "Unfertilized")) +
   scale_linetype_manual(name = "Treatment",
                         values = c("Y" = 1,
                                    "N" = 2),
                         labels = c("Y" = "Fertilized",
                                    "N" = "Unfertilized")) +
   labs(title = "b.", x = "PC1 (41.0%)", y = "PC2 (32.3%)") +
   theme_classic() +
   theme(axis.title.x = element_text(size=25), 
         axis.title.y = element_text(size=25), 
         axis.text.x = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain", margin = margin(r = 15)),
         axis.text.y = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain", margin = margin(r = 15)),
         strip.text.x = element_text(size = 25),
         plot.title = element_text(color = "black", size = 25, hjust = 0, vjust = 0, face = "plain"),
         panel.grid.major = element_blank(),  #remove major-grid labels
         panel.grid.minor = element_blank(),  #remove minor-grid labels
         legend.title = element_text(colour = "black", size = 24),
         legend.text = element_text(colour = "black", size = 20),
         legend.position = 'bottom', # move legend to bottom
         legend.box.background = element_rect(colour = "black")) +
   annotate("text", x = -2.75, y = -2.85,
            label = "Habitat - P = 0.001\nTreatment - P = 0.22\nInteraction - P = 0.68",
            hjust = 0, size = 7) 
)

(panel_pca <- ggarrange(male_pca_fig, female_pca_fig, nrow = 2, common.legend = T, legend = "bottom",
                       heights = c(1,1)))

ggsave(filename = here::here("output", "body_quality_pca.png"), panel_pca, width = 10, height = 16.5,
       dpi = "retina")
##### PERMANOVA #####
adonis2(male_cond ~ habitat*fertilizer_trt, by = "terms", data = male)

adonis2(female_cond ~ habitat*fertilizer_trt, by = "terms", data = female)


#### INTESTINAL SURFACE AREA ####
# not a measure of wellness but can potentially help explain trends
# a measure of assimilation efficiency

# combined
comb_isa_mod <- glmmTMB(intestinal_surface_area ~ habitat*fertilizer_trt + (1|pair_id/territory_id),
                        family = Gamma("log"), data = physiology_sub)
summary(comb_isa_mod)
plot(simulateResiduals(comb_isa_mod))
Anova(comb_isa_mod)
comb_isa_emm <- emmeans(comb_isa_mod, ~ habitat*fertilizer_trt)
contrast(comb_isa_emm, "consec", simple = "each", combine = TRUE, adjust = "mvt")

# males only
male_isa_mod <- glmmTMB(intestinal_surface_area ~ habitat*fertilizer_trt + (1|pair_id/territory_id),
                   family = Gamma("log"), data = male)
summary(male_isa_mod)
plot(simulateResiduals(male_isa_mod))
Anova(male_isa_mod)
male_isa_emm <- emmeans(male_isa_mod, ~ habitat*fertilizer_trt)
contrast(male_isa_emm, "consec", simple = "each", maleine = TRUE, adjust = "mvt")

# females only
female_isa_mod <- glmmTMB(intestinal_surface_area ~ habitat*fertilizer_trt + (1|pair_id/territory_id),
                        family = Gamma("log"), data = female)
summary(female_isa_mod)
plot(simulateResiduals(female_isa_mod))
female_isa_emm <- emmeans(female_isa_mod, ~ habitat*fertilizer_trt)
contrast(female_isa_emm, "consec", simple = "each", femaleine = TRUE, adjust = "mvt")

##### PLOT #####
male_isa_letters <- data.frame(
  habitat = c("Fringing", "Fringing", "Backreef", "Backreef"),
  fertilizer_trt = c("N", "Y", "N", "Y"),
  letters = c("a", "a", "a", "b"),
  heights = rep(150, 4)
)

(male_isa_plot <- ggplot() +
  geom_boxplot(data = male, aes(x = fertilizer_trt, y = intestinal_surface_area,
                                        colour = habitat), outlier.shape = NA) +
  geom_point(data = male, aes(x = fertilizer_trt, y = intestinal_surface_area,
                                      colour = habitat), position = position_jitterdodge(0.4)) +
  geom_text(data = male_isa_letters, aes(x = fertilizer_trt, y = heights, group = habitat,
                                    label = letters),
            position = position_dodge(0.75), size = 8) +
  scale_y_log10() +
  theme_classic() +
  labs(title = "a.\n", x = "Treatment",
       y = expression(Intestinal ~ surface ~ area ~ (cm^2))) + 
  scale_x_discrete(labels = c("Unfertilized", "Fertilized")) +
  scale_colour_discrete(name = "Habitat") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 1.2),
        legend.title = element_text(colour = "black", size = 24),
        legend.text = element_text(colour = "black", size = 20),
        plot.title = element_text(color = "black", size = 25, hjust = 0, vjust = 0, face = "plain"),
        axis.text.x = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain"),
        axis.title.x = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain", margin = margin(r = 15))
  )
)

female_isa_letters <- data.frame(
  habitat = c("Fringing", "Fringing", "Backreef", "Backreef"),
  fertilizer_trt = c("N", "Y", "N", "Y"),
  letters = c("a", "a", "a", "b"),
  heights = rep(175, 4)
)

(female_isa_plot <- ggplot() +
    geom_boxplot(data = female, aes(x = fertilizer_trt, y = intestinal_surface_area,
                                  colour = habitat), outlier.shape = NA) +
    geom_point(data = female, aes(x = fertilizer_trt, y = intestinal_surface_area,
                                colour = habitat), position = position_jitterdodge(0.4)) +
    geom_text(data = female_isa_letters, aes(x = fertilizer_trt, y = heights, group = habitat,
                                           label = letters),
              position = position_dodge(0.75), size = 8) +
    scale_y_log10() +
    scale_x_discrete(labels = c("Unfertilized", "Fertilized")) +
    scale_colour_discrete(name = "Habitat") +
    theme_classic() +
    labs(title = "b.\n", x = "Treatment",
         y = expression(Intestinal ~ surface ~ area ~ (cm^2))) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", size = 1.2),
          legend.title = element_text(colour = "black", size = 24),
          legend.text = element_text(colour = "black", size = 20),
          plot.title = element_text(color = "black", size = 25, hjust = 0, vjust = 0, face = "plain"),
          axis.text.x = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain"),
          axis.title.x = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain"),
          axis.text.y = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain"),
          axis.title.y = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain", margin = margin(r = 15))
    )
)


(isa_panel <- ggarrange(male_isa_plot, female_isa_plot, nrow = 2, common.legend = T, legend = "bottom",
                        heights = c(1,1)))

ggsave(filename = here::here("output", "intestinal_surface_area_panel.png"), isa_panel, width = 10, height = 16.5,
       dpi = "retina")
