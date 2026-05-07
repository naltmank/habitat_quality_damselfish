rm(list=ls())
# install.packages("librarian")
librarian::shelf(here, reshape, tidyverse, dplyr, effects, ggplot2, ggpubr, ggrepel,
                 glmmTMB, DHARMa, performance, emmeans, car)

#### READ AND FORMAT DATA ####
fish_and_territory <- read.csv(here::here("data","territory_and_size.csv")) %>%
  mutate(week_no = as.ordered(week_no))

###### TERRITORY SIZE ######
# calculate territory surface area
# territory surface area fxn (Blanchette et al. 2019)
territory_surface_area = function(x, y, z){
  2*pi*((x^1.6*y^1.6 + x^1.6*z^1.6 + y^1.6*z^1.6)/3)^(1/1.6)/(100^2)
}
fish_and_territory$surface_area <- territory_surface_area(fish_and_territory$territory_length,
                                                          fish_and_territory$territory_width,
                                                          fish_and_territory$territory_height)

hist(fish_and_territory$surface_area)

# condense to just territory
# multiple fish sizes and counts causes duplicated territory rows
# take mean of territory size within each territory/timepoint to address this - should match the previous values
territory <- fish_and_territory %>%
  select(territory_no, territory_id, pair_id, week_no, habitat, fertilizer_treatment, 
         pre_post_treatment, pre_post_removal, surface_area) %>%
  group_by(territory_no, territory_id, pair_id, week_no, habitat, fertilizer_treatment, 
           pre_post_treatment, pre_post_removal, surface_area) %>%
  summarise(surface_area = mean(surface_area, na.rm = T))

# pre removal data
territory_pre_rmv <- subset(territory, pre_post_removal == "Pre-removal")

# post removal data
territory_post_rmv <- subset(territory, pre_post_removal == "Post-removal")

###### FISH COUNTS #####
# sum the counts per territory per time period
counts <- fish_and_territory %>%
  select(territory_no, territory_id, pair_id, week_no, habitat, fertilizer_treatment, 
         pre_post_treatment, pre_post_removal, surface_area, count) %>%
  group_by(territory_no, territory_id, pair_id, week_no, habitat, fertilizer_treatment, 
           pre_post_treatment, pre_post_removal, surface_area) %>%
  summarise(surface_area = mean(surface_area, na.rm = T), # for density calculation below
            fish_total = sum(count, na.rm = T),
            fish_density = fish_total/surface_area
            )

# some territories dropped to 0 fish and disappeared, making density == NaN - replace with 0
is.nan.data.frame <- function(x){
  do.call(cbind, lapply(x, is.nan))
}

counts[is.nan(counts)] <- 0

###### FISH SIZES ######
sizes_exp <- uncount(fish_and_territory, count) # expand table
zeroes <- subset(fish_and_territory, count == 0) %>% # retain zero counts so as not to lose transects
  # remove count column for merge
  select(-count)

sizes_cond <- rbind(sizes_exp, zeroes) %>%
  select(territory_no, territory_id, pair_id, week_no, habitat, fertilizer_treatment, 
         pre_post_treatment, pre_post_removal, surface_area, fish_size) %>%
  # models don't converge because of the repeated values - take the average fish size per territory
  # will still allow me to assess differences in size pre and post treatment
  group_by(territory_no, territory_id, pair_id, week_no, habitat, fertilizer_treatment, 
           pre_post_treatment, pre_post_removal) %>%
  summarise(size = mean(fish_size, na.rm = T))

 

#### PRE-REMOVAL ####
###### TERRITORY SIZE #####
# effect of treatment and habitat on surface area PRE-REMOVAL
# simplify model interpretation by doing it as relative change vs baseline
territory_post_trt <- territory_pre_rmv %>%
  group_by(territory_id) %>%
  mutate(delta_area = ((surface_area - surface_area[week_no == 0])/surface_area[week_no == 0])*100 ) %>%
  # filter out baseline - 0% change vs itself and not instructive
  filter(week_no != 0)

sa_trt_mod <- glmmTMB(
  delta_area ~ fertilizer_treatment * week_no * habitat + (1 | pair_id) + (1 | territory_id),
  data = territory_post_trt)
summary(sa_trt_mod)
Anova(sa_trt_mod)
plot(simulateResiduals(sa_trt_mod))
sa_trt_emm <- emmeans(sa_trt_mod, ~ fertilizer_treatment * week_no * habitat)
contrast(sa_trt_emm, "consec", simple = "each", combine = TRUE, adjust = "bonferroni")

# plot v1 - boxplots
(sa_treatment_plot_v1 <- ggplot() +
    geom_boxplot(data = territory_post_trt, aes(x = habitat, y = delta_area,
                                  colour = fertilizer_treatment), outlier.shape = NA) +
    geom_point(data = territory_post_trt, aes(x = habitat, y = delta_area,
                                colour = fertilizer_treatment), position = position_jitterdodge(0.4)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "black") +
    facet_wrap(~week_no) +
#    geom_text(data = male_isa_letters, aes(x = fertilizer_trt, y = heights, group = habitat,
#                                           label = letters),
#              position = position_dodge(0.75), size = 8) +
    theme_classic() +
    labs(title = "a.\n", x = "",
         y = "% Change") + 
    scale_colour_discrete(name = "Treatment") +
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

# plot v2 - mean +/- standard error
territory_post_trt_graph_df <- territory_post_trt %>%
  group_by(habitat, fertilizer_treatment, week_no) %>%
  summarise(change_mean = mean(delta_area, na.rm = T),
            change_se = sd(delta_area, na.rm = T)/sqrt(n())) %>%
  mutate(group = paste(habitat, fertilizer_treatment))

(sa_treatment_plot_v2 <- ggplot() +
    geom_point(data = territory_post_trt_graph_df, aes(x = week_no, y = change_mean,
                                              colour = habitat, shape = fertilizer_treatment),
               position = position_dodge(0.4), size = 4) +
    geom_errorbar(data = territory_post_trt_graph_df, aes(x = week_no, ymin = change_mean - change_se,
                                                          ymax = change_mean + change_se,
                                                       colour = habitat, group = group),
               position = position_dodge(0.4), width = 0.25) +
    geom_line(data = territory_post_trt_graph_df, aes(x = week_no, y = change_mean,
                                                      colour = habitat, group = group),
              position = position_dodge(0.4)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "black") +
    #    geom_text(data = male_isa_letters, aes(x = fertilizer_trt, y = heights, group = habitat,
    #                                           label = letters),
    #              position = position_dodge(0.75), size = 8) +
    scale_shape_manual(name = "Treatment",
                       values = c("Y" = 19,
                                  "N" = 1),
                       labels = c("Y" = "Fertilized",
                                  "N" = "Unfertilized")) +
    theme_classic() +
    labs(title = "a.\n", x = "",
         y = "% Change") + 
    scale_colour_discrete(name = "Habitat") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", size = 1.2),
          legend.title = element_text(colour = "black", size = 24),
          legend.text = element_text(colour = "black", size = 20),
          legend.position = "bottom",
          plot.title = element_text(color = "black", size = 25, hjust = 0, vjust = 0, face = "plain"),
          axis.text.x = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain"),
          axis.title.x = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain"),
          axis.text.y = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain"),
          axis.title.y = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain", margin = margin(r = 15))
    )
)


###### FISH DENSITY ######
counts_pre_rmv <- subset(counts, pre_post_removal == "Pre-removal")
# probably log transform with a gamma
hist(log(counts$fish_density + 1))

counts_trt_mod <- glmmTMB(
  fish_density ~ fertilizer_treatment * week_no * habitat + (1 | pair_id) + (1 | territory_id),
  family = Gamma("log"), data = counts_pre_rmv)
plot(simulateResiduals(counts_trt_mod))
summary(counts_trt_mod)
Anova(counts_trt_mod)
count_trt_emm <- emmeans(counts_trt_mod, ~ fertilizer_treatment * week_no * habitat)
contrast(count_trt_emm, "consec", simple = "each", combine = TRUE, adjust = "bonferroni")


(count_treatment_plot_v1 <- ggplot() +
    geom_boxplot(data = counts_pre_rmv, aes(x = habitat, y = fish_density,
                                                colour = fertilizer_treatment), outlier.shape = NA) +
    geom_point(data = counts_pre_rmv, aes(x = habitat, y = fish_density,
                                              colour = fertilizer_treatment), position = position_jitterdodge(0.4)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "black") +
    facet_wrap(~week_no) +
    scale_y_log10() +
    #    geom_text(data = male_isa_letters, aes(x = fertilizer_trt, y = heights, group = habitat,
    #                                           label = letters),
    #              position = position_dodge(0.75), size = 8) +
    theme_classic() +
    labs(title = "a.\n", x = "",
         y = "") + 
    scale_colour_discrete(name = "Treatment") +
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

# plot v2 - mean +/- standard error
counts_pre_rmv_graph_df <- counts_pre_rmv %>%
  group_by(habitat, fertilizer_treatment, week_no) %>%
  summarise(density_mean = mean(fish_density, na.rm = T),
            density_se = sd(fish_density, na.rm = T)/sqrt(n())) %>%
  mutate(group = paste(habitat, fertilizer_treatment))

(count_treatment_plot_v2 <- ggplot() +
    geom_point(data = counts_pre_rmv_graph_df, aes(x = week_no, y = density_mean,
                                                       colour = habitat, shape = fertilizer_treatment),
               position = position_dodge(0.4), size = 4) +
    geom_errorbar(data = counts_pre_rmv_graph_df, aes(x = week_no, ymin = density_mean - density_se,
                                                          ymax = density_mean + density_se,
                                                          colour = habitat, group = group),
                  position = position_dodge(0.4), width = 0.25) +
    geom_line(data = counts_pre_rmv_graph_df, aes(x = week_no, y = density_mean, linetype = fertilizer_treatment,
                                                      colour = habitat, group = group),
              position = position_dodge(0.4)) +
    #    geom_text(data = male_isa_letters, aes(x = fertilizer_trt, y = heights, group = habitat,
    #                                           label = letters),
    #              position = position_dodge(0.75), size = 8) +
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
    theme_classic() +
    labs(title = "b.\n", x = "",
         y = "Density") + 
    scale_colour_discrete(name = "Habitat") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", size = 1.2),
          legend.title = element_text(colour = "black", size = 24),
          legend.text = element_text(colour = "black", size = 20),
          legend.position = "bottom",
          plot.title = element_text(color = "black", size = 25, hjust = 0, vjust = 0, face = "plain"),
          axis.text.x = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain"),
          axis.title.x = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain"),
          axis.text.y = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain"),
          axis.title.y = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain", margin = margin(r = 15))
    )
)


###### SIZES ######
sizes_pre_rmv <- subset(sizes_cond, pre_post_removal == "Pre-removal")
# probably log transform with a gamma
hist( log(sizes_pre_rmv$size + 1) )

size_trt_mod <- glmmTMB(
  size ~ fertilizer_treatment * week_no * habitat + (1 | pair_id) + (1 | territory_id),
  family = Gamma("log"), data = sizes_pre_rmv)
plot(simulateResiduals(size_trt_mod))
summary(size_trt_mod)
Anova(size_trt_mod)
size_trt_emm <- emmeans(size_trt_mod, ~ fertilizer_treatment * week_no * habitat)
contrast(size_trt_emm, "consec", simple = "each", combine = TRUE, adjust = "bonferroni")


size_pre_rmv_graph_df <- sizes_pre_rmv %>%
  group_by(habitat, fertilizer_treatment, week_no) %>%
  summarise(size_mean = mean(size, na.rm = T),
            size_se = sd(size, na.rm = T)/sqrt(n())) %>%
  mutate(group = paste(habitat, fertilizer_treatment))

(size_treatment_plot <- ggplot() +
    geom_point(data = size_pre_rmv_graph_df, aes(x = week_no, y = size_mean,
                                                 colour = habitat, shape = fertilizer_treatment),
               position = position_dodge(0.4), size = 4) +
    geom_errorbar(data = size_pre_rmv_graph_df, aes(x = week_no, ymin = size_mean - size_se,
                                                    ymax = size_mean + size_se,
                                                    colour = habitat, group = group),
                  position = position_dodge(0.4), width = 0.25) +
    geom_line(data = size_pre_rmv_graph_df, aes(x = week_no, y = size_mean, linetype = fertilizer_treatment,
                                                colour = habitat, group = group),
              position = position_dodge(0.4)) +
    #    geom_text(data = male_isa_letters, aes(x = fertilizer_trt, y = heights, group = habitat,
    #                                           label = letters),
    #              position = position_dodge(0.75), size = 8) +
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
    theme_classic() +
    labs(title = "c.\n", x = "",
         y = "Size") + 
    scale_colour_discrete(name = "Habitat") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", size = 1.2),
          legend.title = element_text(colour = "black", size = 24),
          legend.text = element_text(colour = "black", size = 20),
          legend.position = "bottom",
          plot.title = element_text(color = "black", size = 25, hjust = 0, vjust = 0, face = "plain"),
          axis.text.x = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain"),
          axis.title.x = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain"),
          axis.text.y = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain"),
          axis.title.y = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain", margin = margin(r = 15))
    )
)

#### REMOVAL IMPACT ####
territory_removal <- territory %>%
  filter(week_no %in% c(6,9))

territory_removal_mod <- glmmTMB(
  surface_area ~ fertilizer_treatment * week_no * habitat  + (1 | pair_id) + (1 | territory_id),
  family = ziGamma("log"), ziformula = ~ habitat * fertilizer_treatment,
  data = territory_removal
)
plot(simulateResiduals(territory_removal_mod))
summary(territory_removal_mod)
Anova(territory_removal_mod)

counts_removal <- counts %>%
  filter(week_no %in% c(6,9))

counts_removal_mod <- glmmTMB(
  log(fish_density + 1) ~ fertilizer_treatment * week_no * habitat  + (1 | pair_id) + (1 | territory_id),
  data = counts_removal
)
plot(simulateResiduals(counts_removal_mod))
summary(counts_removal_mod)
Anova(counts_removal_mod)


sizes_removal <- sizes_cond %>%
  filter(week_no %in% c(6,9))

size_removal_mod <- glmmTMB(
  size ~ fertilizer_treatment * week_no * habitat  + (1 | pair_id) + (1 | territory_id),
  family = ziGamma("log"), ziformula = ~ habitat * fertilizer_treatment,
  data = sizes_removal
)
plot(simulateResiduals(size_removal_mod))
summary(size_removal_mod)
Anova(size_removal_mod)



#### POST-REMOVAL ####
###### TERRITORY SIZE ######
# treatments were still active in week 9 - add a phase of the experiment row
territory_post_rmv <- rbind(subset(territory_pre_rmv, week_no == 6), territory_post_rmv) 

sa_rmv_mod <- glmmTMB(
  surface_area ~ fertilizer_treatment * week_no * habitat  + (1 | pair_id) + (1 | territory_id),
  family = ziGamma("log"), ziformula = ~ habitat * fertilizer_treatment,
  data = territory_post_rmv)
plot(simulateResiduals(sa_rmv_mod))
summary(sa_rmv_mod)
Anova(sa_rmv_mod)
# plot(effects::allEffects(sa_rmv_mod))
sa_rmv_emm <- emmeans(sa_rmv_mod, ~ fertilizer_treatment * week_no * habitat)
contrast(sa_rmv_emm, "consec", simple = "each", combine = TRUE, adjust = "bonferroni")


###### COUNTS ######
counts_post_rmv <- rbind(subset(counts, week_no == 6), subset(counts, pre_post_removal == "Post-removal")) 

counts_rmv_mod <- glmmTMB(
  fish_density ~ fertilizer_treatment * week_no * habitat  + (1 | pair_id) + (1 | territory_id),
  family = ziGamma("log"), ziformula = ~ habitat * fertilizer_treatment,
  data = counts_post_rmv
)
plot(simulateResiduals(counts_rmv_mod))
summary(counts_rmv_mod)
Anova(counts_rmv_mod)

counts_rmv_emm <- emmeans(counts_rmv_mod, ~ fertilizer_treatment * week_no * habitat)
contrast(counts_rmv_emm, "consec", simple = "each", combine = TRUE, adjust = "bonferroni")

###### FISH SIZE ######
sizes_post_rmv <- rbind(subset(sizes_cond, week_no == 6), subset(sizes_cond, pre_post_removal == "Post-removal")) 


size_rmv_mod <- glmmTMB(
  size ~ fertilizer_treatment * week_no * habitat  + (1 | pair_id) + (1 | territory_id),
  family = ziGamma("log"), ziformula = ~ habitat * fertilizer_treatment,
  data = sizes_post_rmv
)
plot(simulateResiduals(size_rmv_mod))
summary(size_rmv_mod)
Anova(size_rmv_mod)

size_rmv_emm <- emmeans(size_rmv_mod, ~ fertilizer_treatment * week_no * habitat)
contrast(size_rmv_emm, "consec", simple = "each", combine = TRUE, adjust = "bonferroni")


