# statistics overview predictions_temp
stats_predictions_temp <- function(df){
  require(broom)
  require(tidyverse)
  require(mgcv)
  require(modelsummary)
  require(writexl)
  
  species_list <- df |> 
    distinct(Species) |> 
    pull(Species)
  
  species_list <- sort(species_list)
  
  gam_output <- list()
  temp_data <- list()
  summary <- list()
  model_summary <- list()
  total_summary <- list()
  
  df$fLake <- as.factor(df$Lake)
  
  df$fProtocol <- as.factor(df$Protocol)
  
  for (i in species_list) {
    data <- df |> 
      filter(Species == i) |> 
      mutate(n_lake = n_distinct(Lake))
    
    # species occurring in one lake -> net-type as random effect
    if(max(data$n_lake) == 1) {
      
      # special case, because model not running with ZIP
      if (i == "Coregonus_sp_benthic_profundal")  { 
        
        gam_output <- gam(data = data, Presence ~ s(mean_last_7days, k = 3) +
                            s(fProtocol, bs = 're'), family = binomial)
        
        # stats
        glance_summary <- glance(gam_output)
        tidy_summary <- tidy(gam_output)
        
        # Add species column to both summaries
        glance_summary <- glance_summary |> 
          mutate(species = i)
        
        tidy_summary <- tidy_summary |> 
          mutate(species = i)
        
        # abundance column
        abundances <- data |> 
          rename(species = Species) |> 
          distinct(species, tot_abu)
        
        # combine both
        pre_total_summary <- merge(tidy_summary, glance_summary)
        
        total_summary <- merge(pre_total_summary, abundances)
        
        # Store the summary in the list
        summ <- summary(gam_output)
        
        dev_explained <- summ$dev.expl
        total_summary$deviance_explained <- dev_explained
        
        chi_sq <- summ[["chi.sq"]][["s(mean_last_7days)"]]
        
        total_summary$chi_square <- chi_sq
        
        p_value <- summ[["s.pv"]][1]
        total_summary$p_value <- p_value
        
        edf <- summ[["edf"]][1]
        total_summary$edf <- edf
        
        r_sq <- summ[["r.sq"]]
        total_summary$r_sq <- r_sq
        
        model_summary[[i]] <- total_summary
        
        
        # all species with abundance data, except for special case
      } else if (max(data$Abundance) > 1)  {
        
        gam_output <- gam(data = data, Abundance ~ s(mean_last_7days, k = 3) +
                            s(fProtocol, bs = 're'), family = ziP())
        
        # stats
        glance_summary <- glance(gam_output)
        tidy_summary <- tidy(gam_output)
        
        # Add species column to both summaries
        glance_summary <- glance_summary |> 
          mutate(species = i)
        
        tidy_summary <- tidy_summary |> 
          mutate(species = i)
        
        # abundance column
        abundances <- data |> 
          rename(species = Species) |> 
          distinct(species, tot_abu)
        
        # combine both
        pre_total_summary <- merge(tidy_summary, glance_summary)
        
        total_summary <- merge(pre_total_summary, abundances)
        
        # Store the summary in the list
        summ <- summary(gam_output)
        
        dev_explained <- summ$dev.expl
        total_summary$deviance_explained <- dev_explained
        
        chi_sq <- summ[["chi.sq"]][["s(mean_last_7days)"]]
        
        total_summary$chi_square <- chi_sq
        
        p_value <- summ[["s.pv"]][1]
        total_summary$p_value <- p_value
        
        edf <- summ[["edf"]][1]
        total_summary$edf <- edf
        
        r_sq <- summ[["r.sq"]]
        total_summary$r_sq <- r_sq
        
        model_summary[[i]] <- total_summary
        
        
      }
      # binomial data species occuring in one lake
      else {   
        
        
        gam_output <- gam(data = data, Abundance ~ s(mean_last_7days, k = 3) +
                            s(fProtocol, bs = 're'), family = binomial)
        
        # stats
        glance_summary <- glance(gam_output)
        tidy_summary <- tidy(gam_output)
        
        # Add species column to both summaries
        glance_summary <- glance_summary |> 
          mutate(species = i)
        
        tidy_summary <- tidy_summary |> 
          mutate(species = i)
        
        # abundance column
        abundances <- data |> 
          rename(species = Species) |> 
          distinct(species, tot_abu)
        
        # combine both
        pre_total_summary <- merge(tidy_summary, glance_summary)
        
        total_summary <- merge(pre_total_summary, abundances)
        
        # Store the summary in the list
        summ <- summary(gam_output)
        
        dev_explained <- summ$dev.expl
        total_summary$deviance_explained <- dev_explained
        
        chi_sq <- summ[["chi.sq"]][["s(mean_last_7days)"]]
        
        total_summary$chi_square <- chi_sq
        
        p_value <- summ[["s.pv"]][1]
        total_summary$p_value <- p_value
        
        edf <- summ[["edf"]][1]
        total_summary$edf <- edf
        
        r_sq <- summ[["r.sq"]]
        total_summary$r_sq <- r_sq
        
        model_summary[[i]] <- total_summary
        
      }
    }
    # all species occurring in multiple lakes
    else {
      
      
      # simplified models
      if (i %in% c("Alburnus_arborella", "Barbatula_sp_Lineage_I",
                   "Cyprinus_carpio", "Phoxinus_csikii", "Salmo_trutta")){
        
        
        gam_output<- gam(data = data, Presence ~ s(mean_last_7days, k = 3) + s(fLake, bs = 're')
                         +  s(fProtocol, bs = 're'), family = binomial)
        
        # stats
        glance_summary <- glance(gam_output)
        tidy_summary <- tidy(gam_output)
        
        # Add species column to both summaries
        glance_summary <- glance_summary |> 
          mutate(species = i)
        
        tidy_summary <- tidy_summary |> 
          mutate(species = i)
        
        # abundance column
        abundances <- data |> 
          rename(species = Species) |> 
          distinct(species, tot_abu)
        
        # combine both
        pre_total_summary <- merge(tidy_summary, glance_summary)
        
        total_summary <- merge(pre_total_summary, abundances)
        
        # Store the summary in the list
        summ <- summary(gam_output)
        
        dev_explained <- summ$dev.expl
        total_summary$deviance_explained <- dev_explained
        
        chi_sq <- summ[["chi.sq"]][["s(mean_last_7days)"]]
        
        total_summary$chi_square <- chi_sq
        
        p_value <- summ[["s.pv"]][1]
        total_summary$p_value <- p_value
        
        edf <- summ[["edf"]][1]
        total_summary$edf <- edf
        
        r_sq <- summ[["r.sq"]]
        total_summary$r_sq <- r_sq
        
        model_summary[[i]] <- total_summary
        
        
      }
      # models with ZIP and occurring in multiple lakes
      else if (max(data$Abundance) > 1)  { 
        
        
        gam_output<- gam(data = data, Abundance ~ s(mean_last_7days, k = 3) + s(fLake, bs = 're')
                         +  s(fProtocol, bs = 're'), family = ziP())
        
        # stats
        glance_summary <- glance(gam_output)
        tidy_summary <- tidy(gam_output)
        
        # Add species column to both summaries
        glance_summary <- glance_summary |> 
          mutate(species = i)
        
        tidy_summary <- tidy_summary |> 
          mutate(species = i)
        
        # abundance column
        abundances <- data |> 
          rename(species = Species) |> 
          distinct(species, tot_abu)
        
        # combine both
        pre_total_summary <- merge(tidy_summary, glance_summary)
        
        total_summary <- merge(pre_total_summary, abundances)
        
        # Store the summary in the list
        summ <- summary(gam_output)
        
        dev_explained <- summ$dev.expl
        total_summary$deviance_explained <- dev_explained
        
        chi_sq <- summ[["chi.sq"]][["s(mean_last_7days)"]]
        
        total_summary$chi_square <- chi_sq
        
        p_value <- summ[["s.pv"]][1]
        total_summary$p_value <- p_value
        
        edf <- summ[["edf"]][1]
        total_summary$edf <- edf
        
        r_sq <- summ[["r.sq"]]
        total_summary$r_sq <- r_sq
        
        model_summary[[i]] <- total_summary
        
        
        # binomial in multiple lakes
      } else {
        
        gam_output <- gam(data = data, Abundance ~ s(mean_last_7days, k = 3) + s(fLake, bs = 're')
                          +  s(fProtocol, bs = 're'), family = binomial)
        
        # stats
        glance_summary <- glance(gam_output)
        tidy_summary <- tidy(gam_output)
        
        # Add species column to both summaries
        glance_summary <- glance_summary |> 
          mutate(species = i)
        
        tidy_summary <- tidy_summary |> 
          mutate(species = i)
        
        # abundance column
        abundances <- data |> 
          rename(species = Species) |> 
          distinct(species, tot_abu)
        
        # combine both
        pre_total_summary <- merge(tidy_summary, glance_summary)
        
        total_summary <- merge(pre_total_summary, abundances)
        
        # Store the summary in the list
        summ <- summary(gam_output)
        
        dev_explained <- summ$dev.expl
        total_summary$deviance_explained <- dev_explained
        
        chi_sq <- summ[["chi.sq"]][["s(mean_last_7days)"]]
        
        total_summary$chi_square <- chi_sq
        
        p_value <- summ[["s.pv"]][1]
        total_summary$p_value <- p_value
        
        edf <- summ[["edf"]][1]
        total_summary$edf <- edf
        
        r_sq <- summ[["r.sq"]]
        total_summary$r_sq <- r_sq
        
        model_summary[[i]] <- total_summary
        
      }
      
    }
  }
  # Combine all summaries into a single data frame
  summary_df <- bind_rows(model_summary)
  
  # # Export the summary data frame to a text file
  write_xlsx(summary_df, path = "summary_stats_temperature_models.xlsx")   
}
