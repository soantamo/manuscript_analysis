#### Make quick plots for conceptual figure


library(plotly)

# read in ffd data
ffd <- readRDS('data/final_ffd.RDS')

# filter to k3
ffd_k3 <- ffd %>% 
  filter(k == 3, 
         including_electro == T)


p <- ggplot(data = ffd_k3) + 
  geom_line(aes(x = mean_last_7days, y = value, group = name))


# Assuming `p` is your ggplot2 map
ggplotly(p)

# Lota_lota
# Coregonus_sp_benthic_profundal
# Coregonus_sp_balchen
# Perca_fluviatilis
# Rutilus_rutilus

dir.create('figures/conceptual_figure/')


#### Community 1 ----

png('figures/conceptual_figure/Community_1_b_c.png', res = 300, width = 2500*0.9, height = 800*0.9)

ggplot(data = ffd_k3 %>% 
         filter(name %in% c('Lota_lota', 
                            'Coregonus_sp_benthic_profundal', 
                            #'Coregonus_sp_balchen', 
                            'Alosa_agone', 
                            'Rutilus_rutilus')) %>% 
         group_by(name) %>% 
         slice_sample(prop = 0.1)) + 
  geom_line(aes(x = mean_last_7days, y = value, group = name, col = ffd_per_0.1degree)) + 
  theme_bw() + 
  theme(aspect.ratio = 0.5, 
        panel.grid = element_blank(), 
        legend.position = 'none') + 
  ylab('Abundance') + 
  xlab('Temperature') +
  scale_color_gradient2(name = 'response derivative', 
                        low = "blue", mid = "gray50", high = "red",
                        midpoint = 0, limits = c(-0.1, 0.1), oob = scales::squish)  +
  

ggplot(data = ffd_k3 %>% 
         filter(name %in% c('Lota_lota', 
                            'Coregonus_sp_benthic_profundal', 
                           # 'Coregonus_sp_balchen', 
                            'Alosa_agone', 
                            'Rutilus_rutilus')) %>% 
         group_by(name) %>% 
         slice_sample(prop = 0.2)) + 
  geom_rect(aes(xmin = 4,
                xmax = 6,
                ymin = -0.17,
                ymax = 0.1),
             col = 'transparent', fill = '#ACD39E',
              ) +
  geom_rect(aes( xmin = 13,
                 xmax = 15,
                 ymin = -0.17,
                 ymax = 0.1),
            col = 'transparent', fill = '#ACD39E',
  )+ 
  geom_rect(aes( xmin = 23,
                 xmax = 25,
                 ymin = -0.17,
                 ymax = 0.1),
            col = 'transparent', fill = '#ACD39E',
  )+
  geom_line(aes(x = mean_last_7days, y = ffd_per_0.1degree, group = name, col = ffd_per_0.1degree)) + 
  geom_hline(aes(yintercept = 0)) + 
  theme_bw() + 
  theme(aspect.ratio = 0.5, 
        panel.grid = element_blank(), 
        legend.position = 'none') + 
  ylab('Response derivative') + 
  xlab('Temperature') +
  scale_color_gradient2(name = 'response derivative', 
                        low = "blue", mid = "gray50", high = "red",
                        midpoint = 0, limits = c(-0.1, 0.1), oob = scales::squish) 
dev.off()


#### Commmunity 2 ----

png('figures/conceptual_figure/Community_2_b_c.png', res = 300, width = 2500*0.9, height = 800*0.9)

ggplot(data = ffd_k3 %>% 
         filter(name %in% c('Lota_lota', 
                            'Gasterosteus_aculeatus',
                            'Coregonus_sp_benthic_profundal', 
                            'Coregonus_sp_balchen', 
                            'Gymnocephalus_cernua', 
                            'Squalius_cephalus',
                            'Rutilus_rutilus')) %>% 
         group_by(name) %>% 
         slice_sample(prop = 0.1)) + 
  geom_line(aes(x = mean_last_7days, y = value, group = name, col = ffd_per_0.1degree)) + 
  theme_bw() + 
  theme(aspect.ratio = 0.5, 
        panel.grid = element_blank(), 
        legend.position = 'none') + 
  ylab('Abundance') + 
  xlab('Temperature') +
  scale_color_gradient2(name = 'response derivative', 
                        low = "blue", mid = "gray50", high = "red",
                        midpoint = 0, limits = c(-0.1, 0.1), oob = scales::squish) +
  
  
  ggplot(data = ffd_k3 %>% 
           filter(name %in% c('Lota_lota', 
                              'Gasterosteus_aculeatus',
                              'Coregonus_sp_benthic_profundal', 
                              'Coregonus_sp_balchen', 
                              'Gymnocephalus_cernua', 
                              'Squalius_cephalus',
                              'Rutilus_rutilus')) %>% 
           group_by(name) %>% 
           slice_sample(prop = 0.2)) + 
  geom_rect(aes(xmin = 4,
                xmax = 6,
                ymin = -0.17,
                ymax = 0.2),
            col = 'transparent', fill = '#C2A5CF'
  ) +
  geom_rect(aes( xmin = 13,
                 xmax = 15,
                 ymin = -0.17,
                 ymax = 0.2),
            col = 'transparent', fill = '#C2A5CF',
  )+ 
  geom_rect(aes( xmin = 23,
                 xmax = 25,
                 ymin = -0.17,
                 ymax = 0.2),
            col = 'transparent', fill = '#C2A5CF',
  )+
  geom_line(aes(x = mean_last_7days, y = ffd_per_0.1degree, group = name, col = ffd_per_0.1degree)) + 
  geom_hline(aes(yintercept = 0)) + 
  theme_bw() + 
  theme(aspect.ratio = 0.5, 
        panel.grid = element_blank(), 
        legend.position = 'none') + 
  ylab('Response derivative') + 
  xlab('Temperature') +
  scale_color_gradient2(name = 'response derivative', 
                        low = "blue", mid = "gray50", high = "red",
                        midpoint = 0, limits = c(-0.1, 0.1), oob = scales::squish) 
dev.off()




#### Get the three communities and calcuate the response metrics ----

# read in functions for response diversity
source('scripts/functions/functions.R')

# Define fish communities
fish_communities <- list(
  com1 = c('Lota_lota', 
           'Coregonus_sp_benthic_profundal', 
           'Alosa_agone', 
           'Rutilus_rutilus'),
  
  com2 = c('Lota_lota', 
           'Gasterosteus_aculeatus',
           'Coregonus_sp_benthic_profundal', 
           'Coregonus_sp_balchen', 
           'Gymnocephalus_cernua', 
           'Squalius_cephalus',
           'Rutilus_rutilus')
)

# Community classification logic as a function
classify_community <- function(x) {
  case_when(
    x > 4  & x < 6  ~ '4-6',
    x > 13 & x < 15 ~ '13-15',
    x > 23 & x < 25 ~ '23-25',
    TRUE            ~ 'remove'
  )
}

# Summary function for each fish set
summarise_community <- function(fish_vector) {
  ffd_k3 %>%
    filter(name %in% fish_vector) %>%
    mutate(community = classify_community(mean_last_7days)) %>%
    filter(community != 'remove') %>%
    group_by(community) %>%
    summarise(
      dissimilarity_ross = resp_div(unlist(ffd_per_0.1degree), sign_sens = FALSE),
      divergence_ross    = resp_div(unlist(ffd_per_0.1degree), sign_sens = TRUE),
      mean_ffd = mean(unlist(ffd_per_0.1degree)),
      sd_ffd   = sd(unlist(ffd_per_0.1degree)),
      .groups = "drop"
    )
}

# Apply to each community
com_k3_list <- lapply(fish_communities, summarise_community)

# If needed, combine into one data frame
com_k3_all <- bind_rows(com_k3_list, .id = "fish_set")

com_k3_all <- com_k3_all %>% 
  rename(dissimilarity = dissimilarity_ross, 
         divergence = divergence_ross, 
         direction = mean_ffd, 
         SD = sd_ffd)

# Make long for plotting
com_long_rd <- com_k3_all %>% pivot_longer(3:ncol(.))
  
# relevel community factor
com_long_rd$community <- factor(com_long_rd$community, levels = unique(com_long_rd$community)[c(3,1,2)])
com_long_rd$name <- factor(com_long_rd$name, levels = unique(com_long_rd$name)[c(1,4,2,3)])

png('figures/conceptual_figure/RD_metrics.png', res = 300, width = 1200, height = 1200)
ggplot(data = com_long_rd) + 
  geom_hline(data = data.frame(name = factor(levels(com_long_rd$name), levels = levels(com_long_rd$name)),
                               yint = c(NA, NA, NA, 0)), aes(yintercept = yint)) + 
  geom_point(aes(x = community, y = value, col = fish_set), size = 3) +
  geom_line(aes(x = community, y = value, col = fish_set, group = fish_set), lwd = 1.5) + 
  facet_wrap(~name, scales = 'free_y', nrow =  2) + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        panel.grid = element_blank(), 
        legend.position = 'none') + 
  xlab('example temperature bins (\u00B0C)') + 
  ylab(NULL) + 
  scale_color_manual(values = c("com1" = "#ACD39E",  # Light Purple
                                "com2" = "#C2A5CF")) # Light Green
dev.off()

