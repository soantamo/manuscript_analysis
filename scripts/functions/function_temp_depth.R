
predictions_depth_temp <- function(df){
  require(broom)
  require(tidyverse)
  require(DHARMa)
  require(mgcv)
  require(gratia)
  require(mgcViz)
  require(modelsummary)
  require(writexl)
  
  species_list <- df |> 
    distinct(Species) |> 
    pull(Species)
  
  species_list <- sort(species_list)
  
  gam_output <- list()
  model_prediction <- list()
  derivatives <- list()
  grid <- list()
  pred_df <- list()
  tiff_filename <- list()
  tiff_file_2 <- list()
  tiff_file_3 <- list()
  temp_data <- list()
  depth_data <- list()
  summary <- list()
  model_summary <- list()
  
  df$fLake <- as.factor(df$Lake)
  
  df$fProtocol <- as.factor(df$Protocol)
  
  for (i in species_list) {
    data <- df |> 
      filter(Species == i) |> 
      mutate(n_lake = n_distinct(Lake))
    
    unique_lakes <- unique(data$fLake)
    
    random_lake <- sample(levels(unique_lakes), 1)
    
    # temp depth models
    if(max(data$tot_abu) > 20){
      # one lake
      if(max(data$n_lake) == 1) {
        
        grid <- expand.grid(mean_last_7days = seq(
          from = min(data$mean_last_7days, na.rm = TRUE),
          to = max(data$mean_last_7days, na.rm = TRUE), by = 0.02), 
          Depth_sample = mean(data$Depth_sample, na.rm = TRUE),
          fProtocol = factor("VERT"))
        
        if (i == "Coregonus_sp_benthic_profundal")  { 
          
          print(paste(data$tot_abu[1], i))
          
          gam_output <- gam(data = data, Presence ~ s(mean_last_7days, k = 3) + s(Depth_sample, k = 3)
                            +  s(fProtocol, bs = 're'), family = binomial)
          
          # prepare residuals
          simulationOutput <- simulateResiduals(fittedModel = gam_output, plot = F)
          tiff_filename <- paste("total_models/gam_check/gam_check_", i, ".tiff", sep = "")
          tiff(tiff_filename, width = 800, height = 600)
          print(plot(simulationOutput))
          dev.off()
          
          # Plotting standardized residuals against predictors
          tiff_file_2 <- paste("total_models/gam_check/predictor_", i, ".tiff", sep = "")
          tiff(tiff_file_2, width = 800, height = 600)
          print(plotResiduals(simulationOutput, temp_data$mean_last_7days, xlab = "temp", main=NULL))
          dev.off()
          
          # Plotting standardized residuals against predictors: depth
          tiff_file_3 <- paste("total_models/gam_check/depth_predictor_", i, ".tiff", sep = "")
          tiff(tiff_file_3, width = 800, height = 600)
          print(plotResiduals(simulationOutput, depth_data$Depth_sample, xlab = "depth", main=NULL))
          dev.off()
          
          model_prediction <- predict.gam(gam_output, newdata = grid,
                                          exclude = "s(fProtocol)",
                                          type = "response", se.fit = TRUE)
          model_bind <- cbind(grid, as.data.frame(model_prediction))
          pred_df <- model_bind |>
            rename(temp = mean_last_7days) |>
            mutate(species = factor(i))
          
          # saving the predictions
          saveRDS(pred_df, paste0("total_models/predictions/depth_temp/predictions_",i,".rds"))
          
          # # stats
          # glance_summary <- glance(gam_output)
          # tidy_summary <- tidy(gam_output)
          # 
          # # Add species column to both summaries
          # glance_summary <- glance_summary |> 
          #   mutate(species = i)
          # 
          # tidy_summary <- tidy_summary |> 
          #   mutate(species = i)
          # 
          # # combine both
          # total_summary <- merge(tidy_summary, glance_summary)
          # 
          # # Store the summary in the list
          # summ <- summary(gam_output)
          # 
          # dev_explained <- summ$dev.expl
          # total_summary$deviance_explained <- dev_explained
          # 
          # chi_sq <- summ[["chi.sq"]][["s(mean_last_7days)"]]
          # 
          # total_summary$chi_square <- chi_sq
          # 
          # p_value <- summ[["s.pv"]][1]
          # total_summary$p_value <- p_value
          # 
          # edf <- summ[["edf"]][1]
          # total_summary$edf <- edf
          # 
          # r_sq <- summ[["r.sq"]]
          # total_summary$r_sq <- r_sq
          # 
          # model_summary[[i]] <- total_summary
          
        }
        # abundance
        else if (max(data$Abundance) > 1)  {
          print(paste(data$tot_abu[1], i))
          
          gam_output <- gam(data = data, Abundance ~ s(mean_last_7days, k = 3) + s(Depth_sample, k = 3)
                            +  s(fProtocol, bs = 're'), family = ziP())
          
          # prepare residuals
          simulationOutput <- simulateResiduals(fittedModel = gam_output, plot = F)
          tiff_filename <- paste("total_models/gam_check/gam_check_", i, ".tiff", sep = "")
          tiff(tiff_filename, width = 800, height = 600)
          print(plot(simulationOutput))
          dev.off()
          
          # Plotting standardized residuals against predictors
          tiff_file_2 <- paste("total_models/gam_check/predictor_", i, ".tiff", sep = "")
          tiff(tiff_file_2, width = 800, height = 600)
          print(plotResiduals(simulationOutput, temp_data$mean_last_7days, xlab = "temp", main=NULL))
          dev.off()
          
          # Plotting standardized residuals against predictors: depth
          tiff_file_3 <- paste("total_models/gam_check/depth_predictor_", i, ".tiff", sep = "")
          tiff(tiff_file_3, width = 800, height = 600)
          print(plotResiduals(simulationOutput, depth_data$Depth_sample, xlab = "depth", main=NULL))
          dev.off()
          
          model_prediction <- predict.gam(gam_output, newdata = grid,
                                          exclude = "s(fProtocol)",
                                          type = "response", se.fit = TRUE)
          model_bind <- cbind(grid, as.data.frame(model_prediction))
          pred_df <- model_bind |>
            rename(temp = mean_last_7days) |>
            mutate(species = factor(i))
          
          # saving the predictions
          saveRDS(pred_df, paste0("total_models/predictions/depth_temp/predictions_",i,".rds"))
          
          # # stats
          # glance_summary <- glance(gam_output)
          # tidy_summary <- tidy(gam_output)
          # 
          # # Add species column to both summaries
          # glance_summary <- glance_summary |> 
          #   mutate(species = i)
          # 
          # tidy_summary <- tidy_summary |> 
          #   mutate(species = i)
          # 
          # # combine both
          # total_summary <- merge(tidy_summary, glance_summary)
          # 
          # # Store the summary in the list
          # summ <- summary(gam_output)
          # 
          # dev_explained <- summ$dev.expl
          # total_summary$deviance_explained <- dev_explained
          # 
          # chi_sq <- summ[["chi.sq"]][["s(mean_last_7days)"]]
          # 
          # total_summary$chi_square <- chi_sq
          # 
          # p_value <- summ[["s.pv"]][1]
          # total_summary$p_value <- p_value
          # 
          # edf <- summ[["edf"]][1]
          # total_summary$edf <- edf
          # 
          # r_sq <- summ[["r.sq"]]
          # total_summary$r_sq <- r_sq
          # 
          # model_summary[[i]] <- total_summary
        }
        # binomial
        else{
          print(paste(data$tot_abu[1], i))
          
          gam_output <- gam(data = data, Abundance ~ s(mean_last_7days, k = 3) + s(Depth_sample, k = 3)
                            +  s(fProtocol, bs = 're'), family = binomial)
          
          # prepare residuals
          simulationOutput <- simulateResiduals(fittedModel = gam_output, plot = F)
          tiff_filename <- paste("total_models/gam_check/gam_check_", i, ".tiff", sep = "")
          tiff(tiff_filename, width = 800, height = 600)
          print(plot(simulationOutput))
          dev.off()
          
          # Plotting standardized residuals against predictors
          tiff_file_2 <- paste("total_models/gam_check/predictor_", i, ".tiff", sep = "")
          tiff(tiff_file_2, width = 800, height = 600)
          print(plotResiduals(simulationOutput, temp_data$mean_last_7days, xlab = "temp", main=NULL))
          dev.off()
          
          # Plotting standardized residuals against predictors: depth
          tiff_file_3 <- paste("total_models/gam_check/depth_predictor_", i, ".tiff", sep = "")
          tiff(tiff_file_3, width = 800, height = 600)
          print(plotResiduals(simulationOutput, depth_data$Depth_sample, xlab = "depth", main=NULL))
          dev.off()
          
          model_prediction <- predict.gam(gam_output, newdata = grid,
                                          exclude = "s(fProtocol)",
                                          type = "response", se.fit = TRUE)
          model_bind <- cbind(grid, as.data.frame(model_prediction))
          pred_df <- model_bind |>
            rename(temp = mean_last_7days) |>
            mutate(species = factor(i))
          
          # saving the predictions
          saveRDS(pred_df, paste0("total_models/predictions/depth_temp/predictions_",i,".rds"))
          
          # # stats
          # glance_summary <- glance(gam_output)
          # tidy_summary <- tidy(gam_output)
          # 
          # # Add species column to both summaries
          # glance_summary <- glance_summary |> 
          #   mutate(species = i)
          # 
          # tidy_summary <- tidy_summary |> 
          #   mutate(species = i)
          # 
          # # combine both
          # total_summary <- merge(tidy_summary, glance_summary)
          # 
          # # Store the summary in the list
          # summ <- summary(gam_output)
          # 
          # dev_explained <- summ$dev.expl
          # total_summary$deviance_explained <- dev_explained
          # 
          # chi_sq <- summ[["chi.sq"]][["s(mean_last_7days)"]]
          # 
          # total_summary$chi_square <- chi_sq
          # 
          # p_value <- summ[["s.pv"]][1]
          # total_summary$p_value <- p_value
          # 
          # edf <- summ[["edf"]][1]
          # total_summary$edf <- edf
          # 
          # r_sq <- summ[["r.sq"]]
          # total_summary$r_sq <- r_sq
          # 
          # model_summary[[i]] <- total_summary
          
        }
      }
      # multiple lakes
      else{
        
        unique_lakes <- distinct(data, Lake) |>
          pull()
        
        random_lake <- sample(unique_lakes, 1)
        
        grid <- expand.grid(mean_last_7days = seq(
          from = min(data$mean_last_7days, na.rm = TRUE),
          to = max(data$mean_last_7days, na.rm = TRUE), by = 0.02),
          Depth_sample = mean(data$Depth_sample, na.rm = TRUE),
          fProtocol = factor("VERT"), fLake = factor(random_lake))
        
        # special cases
        if (i %in% c("Alosa_agone", "Alburnus_arborella", "Barbatula_sp_Lineage_I",
                     "Coregonus_sp_large_pelagic", 
                     "Cyprinus_carpio", "Lepomis_gibbosus", "Phoxinus_csikii", "Salmo_trutta")){
          
          print(paste(data$tot_abu[1], i))
          
          gam_output <- gam(data = data, Presence ~ s(mean_last_7days, k = 3) + s(Depth_sample, k = 3)
                            + s(fLake, bs = 're')
                            +  s(fProtocol, bs = 're'), family = binomial)
          
          # prepare residuals
          simulationOutput <- simulateResiduals(fittedModel = gam_output, plot = F)
          tiff_filename <- paste("total_models/gam_check/gam_check_", i, ".tiff", sep = "")
          tiff(tiff_filename, width = 800, height = 600)
          print(plot(simulationOutput))
          dev.off()
          
          # Plotting standardized residuals against predictors
          tiff_file_2 <- paste("total_models/gam_check/predictor_", i, ".tiff", sep = "")
          tiff(tiff_file_2, width = 800, height = 600)
          print(plotResiduals(simulationOutput, temp_data$mean_last_7days, xlab = "temp", main=NULL))
          dev.off()
          
          # Plotting standardized residuals against predictors: depth
          tiff_file_3 <- paste("total_models/gam_check/depth_predictor_", i, ".tiff", sep = "")
          tiff(tiff_file_3, width = 800, height = 600)
          print(plotResiduals(simulationOutput, depth_data$Depth_sample, xlab = "depth", main=NULL))
          dev.off()
          
          model_prediction <- predict.gam(gam_output, newdata = grid,
                                          exclude = c("s(fProtocol)", "s(fLake)"),
                                          type = "response", se.fit = TRUE)
          model_bind <- cbind(grid, as.data.frame(model_prediction))
          pred_df <- model_bind |>
            rename(temp = mean_last_7days) |>
            mutate(species = factor(i))
          
          saveRDS(pred_df, paste0("total_models/predictions/depth_temp/predictions_",i,".rds"))
          
          # # stats
          # glance_summary <- glance(gam_output)
          # tidy_summary <- tidy(gam_output)
          # 
          # # Add species column to both summaries
          # glance_summary <- glance_summary |> 
          #   mutate(species = i)
          # 
          # tidy_summary <- tidy_summary |> 
          #   mutate(species = i)
          # 
          # # combine both
          # total_summary <- merge(tidy_summary, glance_summary)
          # 
          # # Store the summary in the list
          # summ <- summary(gam_output)
          # 
          # dev_explained <- summ$dev.expl
          # total_summary$deviance_explained <- dev_explained
          # 
          # chi_sq <- summ[["chi.sq"]][["s(mean_last_7days)"]]
          # 
          # total_summary$chi_square <- chi_sq
          # 
          # p_value <- summ[["s.pv"]][1]
          # total_summary$p_value <- p_value
          # 
          # edf <- summ[["edf"]][1]
          # total_summary$edf <- edf
          # 
          # r_sq <- summ[["r.sq"]]
          # total_summary$r_sq <- r_sq
          # 
          # model_summary[[i]] <- total_summary
        }
        # abundance
        else if (max(data$Abundance) > 1)  {
          
          print(paste(data$tot_abu[1], i))
          
          gam_output <- gam(data = data, Abundance ~ s(mean_last_7days, k = 3) + s(Depth_sample, k = 3)
                            + s(fLake, bs = 're')
                            +  s(fProtocol, bs = 're'), family = ziP())
          
          # prepare residuals
          simulationOutput <- simulateResiduals(fittedModel = gam_output, plot = F)
          tiff_filename <- paste("total_models/gam_check/gam_check_", i, ".tiff", sep = "")
          tiff(tiff_filename, width = 800, height = 600)
          print(plot(simulationOutput))
          dev.off()
          
          # Plotting standardized residuals against predictors
          tiff_file_2 <- paste("total_models/gam_check/predictor_", i, ".tiff", sep = "")
          tiff(tiff_file_2, width = 800, height = 600)
          print(plotResiduals(simulationOutput, temp_data$mean_last_7days, xlab = "temp", main=NULL))
          dev.off()
          
          # Plotting standardized residuals against predictors: depth
          tiff_file_3 <- paste("total_models/gam_check/depth_predictor_", i, ".tiff", sep = "")
          tiff(tiff_file_3, width = 800, height = 600)
          print(plotResiduals(simulationOutput, depth_data$Depth_sample, xlab = "depth", main=NULL))
          dev.off()
          
          model_prediction <- predict.gam(gam_output, newdata = grid,
                                          exclude = c("s(fProtocol)", "s(fLake)"),
                                          type = "response", se.fit = TRUE)
          model_bind <- cbind(grid, as.data.frame(model_prediction))
          pred_df <- model_bind |>
            rename(temp = mean_last_7days) |>
            mutate(species = factor(i))
          
          saveRDS(pred_df, paste0("total_models/predictions/predictions_",i,".rds"))
          
          # # stats
          # glance_summary <- glance(gam_output)
          # tidy_summary <- tidy(gam_output)
          # 
          # # Add species column to both summaries
          # glance_summary <- glance_summary |> 
          #   mutate(species = i)
          # 
          # tidy_summary <- tidy_summary |> 
          #   mutate(species = i)
          # 
          # # combine both
          # total_summary <- merge(tidy_summary, glance_summary)
          # 
          # # Store the summary in the list
          # summ <- summary(gam_output)
          # 
          # dev_explained <- summ$dev.expl
          # total_summary$deviance_explained <- dev_explained
          # 
          # chi_sq <- summ[["chi.sq"]][["s(mean_last_7days)"]]
          # 
          # total_summary$chi_square <- chi_sq
          # 
          # p_value <- summ[["s.pv"]][1]
          # total_summary$p_value <- p_value
          # 
          # edf <- summ[["edf"]][1]
          # total_summary$edf <- edf
          # 
          # r_sq <- summ[["r.sq"]]
          # total_summary$r_sq <- r_sq
          # 
          # model_summary[[i]] <- total_summary
        }
        # binomial
        else {
          
          print(paste(data$tot_abu[1], i))
          
          gam_output <- gam(data = data, Abundance ~ s(mean_last_7days, k = 3) + s(Depth_sample, k = 3)
                            + s(fLake, bs = 're')
                            +  s(fProtocol, bs = 're'), family = binomial)
          # prepare residuals
          simulationOutput <- simulateResiduals(fittedModel = gam_output, plot = F)
          tiff_filename <- paste("total_models/gam_check/gam_check_", i, ".tiff", sep = "")
          tiff(tiff_filename, width = 800, height = 600)
          print(plot(simulationOutput))
          dev.off()
          
          # Plotting standardized residuals against predictors
          tiff_file_2 <- paste("total_models/gam_check/predictor_", i, ".tiff", sep = "")
          tiff(tiff_file_2, width = 800, height = 600)
          print(plotResiduals(simulationOutput, temp_data$mean_last_7days, xlab = "temp", main=NULL))
          dev.off()
          
          # Plotting standardized residuals against predictors: depth
          tiff_file_3 <- paste("total_models/gam_check/depth_predictor_", i, ".tiff", sep = "")
          tiff(tiff_file_3, width = 800, height = 600)
          print(plotResiduals(simulationOutput, depth_data$Depth_sample, xlab = "depth", main=NULL))
          dev.off()
          
          model_prediction <- predict.gam(gam_output, newdata = grid,
                                          exclude = c("s(fProtocol)", "s(fLake)"),
                                          type = "response", se.fit = TRUE)
          model_bind <- cbind(grid, as.data.frame(model_prediction))
          pred_df <- model_bind |>
            rename(temp = mean_last_7days) |>
            mutate(species = factor(i))
          
          saveRDS(pred_df, paste0("total_models/predictions/depth_temp/predictions_",i,".rds"))
          
          # # stats
          # glance_summary <- glance(gam_output)
          # tidy_summary <- tidy(gam_output)
          # 
          # # Add species column to both summaries
          # glance_summary <- glance_summary |> 
          #   mutate(species = i)
          # 
          # tidy_summary <- tidy_summary |> 
          #   mutate(species = i)
          # 
          # # combine both
          # total_summary <- merge(tidy_summary, glance_summary)
          # 
          # # Store the summary in the list
          # summ <- summary(gam_output)
          # 
          # dev_explained <- summ$dev.expl
          # total_summary$deviance_explained <- dev_explained
          # 
          # chi_sq <- summ[["chi.sq"]][["s(mean_last_7days)"]]
          # 
          # total_summary$chi_square <- chi_sq
          # 
          # p_value <- summ[["s.pv"]][1]
          # total_summary$p_value <- p_value
          # 
          # edf <- summ[["edf"]][1]
          # total_summary$edf <- edf
          # 
          # r_sq <- summ[["r.sq"]]
          # total_summary$r_sq <- r_sq
          # 
          # model_summary[[i]] <- total_summary
        }
        
      }
      
    }
    # TEMPERATURE ONLY 
    else{
      
      # one lake
      if(max(data$n_lake) == 1) {
        
        grid <- expand.grid(mean_last_7days = seq(
          from = min(data$mean_last_7days, na.rm = TRUE),
          to = max(data$mean_last_7days, na.rm = TRUE), by = 0.02),
          fProtocol = factor("VERT"))
        
        if (i == "Coregonus_sp_benthic_profundal")  { 
          
          print(paste(data$tot_abu[1], i))
          
          gam_output <- gam(data = data, Presence ~ s(mean_last_7days, k = 3) +
                              s(fProtocol, bs = 're'), family = binomial)
          # prepare residuals
          simulationOutput <- simulateResiduals(fittedModel = gam_output, plot = F)
          tiff_filename <- paste("total_models/gam_check/gam_check_", i, ".tiff", sep = "")
          tiff(tiff_filename, width = 800, height = 600)
          print(plot(simulationOutput))
          dev.off()
          
          # Plotting standardized residuals against predictors
          tiff_file_2 <- paste("total_models/gam_check/predictor_", i, ".tiff", sep = "")
          tiff(tiff_file_2, width = 800, height = 600)
          print(plotResiduals(simulationOutput, temp_data$mean_last_7days, xlab = "temp", main=NULL))
          dev.off()
          
          
          model_prediction <- predict.gam(gam_output, newdata = grid,
                                          exclude = "s(fProtocol)",
                                          type = "response", se.fit = TRUE)
          model_bind <- cbind(grid, as.data.frame(model_prediction))
          pred_df <- model_bind |>
            rename(temp = mean_last_7days) |>
            mutate(species = factor(i))
          
          # saving the predictions
          saveRDS(pred_df, paste0("total_models/predictions/depth_temp/predictions_",i,".rds"))
          
          # # stats
          # glance_summary <- glance(gam_output)
          # tidy_summary <- tidy(gam_output)
          # 
          # # Add species column to both summaries
          # glance_summary <- glance_summary |> 
          #   mutate(species = i)
          # 
          # tidy_summary <- tidy_summary |> 
          #   mutate(species = i)
          # 
          # # combine both
          # total_summary <- merge(tidy_summary, glance_summary)
          # 
          # # Store the summary in the list
          # summ <- summary(gam_output)
          # 
          # dev_explained <- summ$dev.expl
          # total_summary$deviance_explained <- dev_explained
          # 
          # chi_sq <- summ[["chi.sq"]][["s(mean_last_7days)"]]
          # 
          # total_summary$chi_square <- chi_sq
          # 
          # p_value <- summ[["s.pv"]][1]
          # total_summary$p_value <- p_value
          # 
          # edf <- summ[["edf"]][1]
          # total_summary$edf <- edf
          # 
          # r_sq <- summ[["r.sq"]]
          # total_summary$r_sq <- r_sq
          # 
          # model_summary[[i]] <- total_summary
          
        }
        # abundance
        else if (max(data$Abundance) > 1)  {
          print(paste(data$tot_abu[1], i))
          
          gam_output <- gam(data = data, Abundance ~ s(mean_last_7days, k = 3) +
                              s(fProtocol, bs = 're'), family = ziP())
          # prepare residuals
          simulationOutput <- simulateResiduals(fittedModel = gam_output, plot = F)
          tiff_filename <- paste("total_models/gam_check/gam_check_", i, ".tiff", sep = "")
          tiff(tiff_filename, width = 800, height = 600)
          print(plot(simulationOutput))
          dev.off()
          
          # Plotting standardized residuals against predictors
          tiff_file_2 <- paste("total_models/gam_check/predictor_", i, ".tiff", sep = "")
          tiff(tiff_file_2, width = 800, height = 600)
          print(plotResiduals(simulationOutput, temp_data$mean_last_7days, xlab = "temp", main=NULL))
          dev.off()
          
          
          model_prediction <- predict.gam(gam_output, newdata = grid,
                                          exclude= "s(fProtocol)",
                                          type = "response", 
                                          se.fit = TRUE)
          model_bind <- cbind(grid, as.data.frame(model_prediction))
          pred_df <- model_bind |>
            rename(temp = mean_last_7days) |>
            mutate(species = factor(i))
          
          saveRDS(pred_df, paste0("total_models/predictions/depth_temp/predictions_",i,".rds"))
          
          # # stats
          # glance_summary <- glance(gam_output)
          # tidy_summary <- tidy(gam_output)
          # 
          # # Add species column to both summaries
          # glance_summary <- glance_summary |> 
          #   mutate(species = i)
          # 
          # tidy_summary <- tidy_summary |> 
          #   mutate(species = i)
          # 
          # # combine both
          # total_summary <- merge(tidy_summary, glance_summary)
          # 
          # # Store the summary in the list
          # summ <- summary(gam_output)
          # 
          # dev_explained <- summ$dev.expl
          # total_summary$deviance_explained <- dev_explained
          # 
          # chi_sq <- summ[["chi.sq"]][["s(mean_last_7days)"]]
          # 
          # total_summary$chi_square <- chi_sq
          # 
          # p_value <- summ[["s.pv"]][1]
          # total_summary$p_value <- p_value
          # 
          # edf <- summ[["edf"]][1]
          # total_summary$edf <- edf
          # 
          # r_sq <- summ[["r.sq"]]
          # total_summary$r_sq <- r_sq
          # 
          # model_summary[[i]] <- total_summary
          
        }
        # binomial
        else{
          print(paste(data$tot_abu[1], i))
          
          gam_output <- gam(data = data, Abundance ~ s(mean_last_7days, k = 3) +
                              s(fProtocol, bs = 're'), family = binomial)
          # prepare residuals
          simulationOutput <- simulateResiduals(fittedModel = gam_output, plot = F)
          tiff_filename <- paste("total_models/gam_check/gam_check_", i, ".tiff", sep = "")
          tiff(tiff_filename, width = 800, height = 600)
          print(plot(simulationOutput))
          dev.off()
          
          # Plotting standardized residuals against predictors
          tiff_file_2 <- paste("total_models/gam_check/predictor_", i, ".tiff", sep = "")
          tiff(tiff_file_2, width = 800, height = 600)
          print(plotResiduals(simulationOutput, temp_data$mean_last_7days, xlab = "temp", main=NULL))
          dev.off()
          
          
          model_prediction <- predict.gam(gam_output, newdata = grid,
                                          exclude = "s(fProtocol)",
                                          type = "response",
                                          se.fit = TRUE)
          model_bind <- cbind(grid, as.data.frame(model_prediction))
          pred_df <- model_bind |>
            rename(temp = mean_last_7days) |>
            mutate(species = factor(i))
          saveRDS(pred_df, paste0("total_models/predictions/depth_temp/predictions_",i,".rds"))
          
          # # stats
          # glance_summary <- glance(gam_output)
          # tidy_summary <- tidy(gam_output)
          # 
          # # Add species column to both summaries
          # glance_summary <- glance_summary |> 
          #   mutate(species = i)
          # 
          # tidy_summary <- tidy_summary |> 
          #   mutate(species = i)
          # 
          # # combine both
          # total_summary <- merge(tidy_summary, glance_summary)
          # 
          # # Store the summary in the list
          # summ <- summary(gam_output)
          # 
          # dev_explained <- summ$dev.expl
          # total_summary$deviance_explained <- dev_explained
          # 
          # chi_sq <- summ[["chi.sq"]][["s(mean_last_7days)"]]
          # 
          # total_summary$chi_square <- chi_sq
          # 
          # p_value <- summ[["s.pv"]][1]
          # total_summary$p_value <- p_value
          # 
          # edf <- summ[["edf"]][1]
          # total_summary$edf <- edf
          # 
          # r_sq <- summ[["r.sq"]]
          # total_summary$r_sq <- r_sq
          # 
          # model_summary[[i]] <- total_summary
          # 
        }
      }
      # multiple lakes
      else {
        
        # select random lake to predict them upon, random effects are excluded from
        # the prediction anyway
        unique_lakes <- distinct(data, Lake) |>
          pull()
        
        random_lake <- sample(unique_lakes, 1)
        
        grid <- expand.grid(mean_last_7days = seq(
          from = min(data$mean_last_7days, na.rm = TRUE),
          to = max(data$mean_last_7days, na.rm = TRUE), by = 0.02),
          fProtocol = factor("VERT"), fLake = factor(random_lake))
        
        # special cases
        if (i %in% c("Alosa_agone", "Alburnus_arborella", "Barbatula_sp_Lineage_I",
                     "Coregonus_sp_large_pelagic", 
                     "Cyprinus_carpio", "Lepomis_gibbosus", "Phoxinus_csikii", "Salmo_trutta")){
          print(i)
          
          gam_output<- gam(data = data, Presence ~ s(mean_last_7days, k = 3) + s(fLake, bs = 're')
                           +  s(fProtocol, bs = 're'), family = binomial)
          
          # prepare residuals
          simulationOutput <- simulateResiduals(fittedModel = gam_output, plot = F)
          tiff_filename <- paste("total_models/gam_check/gam_check_", i, ".tiff", sep = "")
          tiff(tiff_filename, width = 800, height = 600)
          print(plot(simulationOutput))
          dev.off()
          
          # Plotting standardized residuals against predictors
          tiff_file_2 <- paste("total_models/gam_check/predictor_", i, ".tiff", sep = "")
          tiff(tiff_file_2, width = 800, height = 600)
          print(plotResiduals(simulationOutput, temp_data$mean_last_7days, xlab = "temp", main=NULL))
          dev.off()
          
          
          model_prediction <- predict.gam(gam_output, newdata = grid,
                                          exclude = c("s(fProtocol)", "s(fLake)"),
                                          type = "response", se.fit = TRUE)
          model_bind <- cbind(grid, as.data.frame(model_prediction))
          pred_df <- model_bind |>
            rename(temp = mean_last_7days) |>
            mutate(species = factor(i))
          
          saveRDS(pred_df, paste0("total_models/predictions/depth_temp/predictions_",i,".rds"))
          
          # # stats
          # glance_summary <- glance(gam_output)
          # tidy_summary <- tidy(gam_output)
          # 
          # # Add species column to both summaries
          # glance_summary <- glance_summary |> 
          #   mutate(species = i)
          # 
          # tidy_summary <- tidy_summary |> 
          #   mutate(species = i)
          # 
          # # combine both
          # total_summary <- merge(tidy_summary, glance_summary)
          # 
          # # Store the summary in the list
          # summ <- summary(gam_output)
          # 
          # dev_explained <- summ$dev.expl
          # total_summary$deviance_explained <- dev_explained
          # 
          # chi_sq <- summ[["chi.sq"]][["s(mean_last_7days)"]]
          # 
          # total_summary$chi_square <- chi_sq
          # 
          # p_value <- summ[["s.pv"]][1]
          # total_summary$p_value <- p_value
          # 
          # edf <- summ[["edf"]][1]
          # total_summary$edf <- edf
          # 
          # r_sq <- summ[["r.sq"]]
          # total_summary$r_sq <- r_sq
          # 
          # model_summary[[i]] <- total_summary
        }
        # abundance
        else if (max(data$Abundance) > 1)  {
          print(paste(data$tot_abu[1], i))
          # should be one
          gam_output<- gam(data = data, Abundance ~ s(mean_last_7days, k = 3) + s(fLake, bs = 're')
                           +  s(fProtocol, bs = 're'), family = ziP())
          # prepare residuals
          simulationOutput <- simulateResiduals(fittedModel = gam_output, plot = F)
          tiff_filename <- paste("total_models/gam_check/gam_check_", i, ".tiff", sep = "")
          tiff(tiff_filename, width = 800, height = 600)
          print(plot(simulationOutput))
          dev.off()
          
          # Plotting standardized residuals against predictors
          tiff_file_2 <- paste("total_models/gam_check/predictor_", i, ".tiff", sep = "")
          tiff(tiff_file_2, width = 800, height = 600)
          print(plotResiduals(simulationOutput, temp_data$mean_last_7days, xlab = "temp", main=NULL))
          dev.off()
          
          
          model_prediction <- predict.gam(gam_output, newdata = grid,
                                          exclude = c("s(fProtocol)", "s(fLake)"),
                                          type = "response", se.fit = TRUE)
          model_bind <- cbind(grid, as.data.frame(model_prediction))
          pred_df <- model_bind |>
            rename(temp = mean_last_7days) |>
            mutate(species = factor(i))
          
          saveRDS(pred_df, paste0("total_models/predictions/depth_temp/predictions_",i,".rds"))
          
          # # stats
          # glance_summary <- glance(gam_output)
          # tidy_summary <- tidy(gam_output)
          # 
          # # Add species column to both summaries
          # glance_summary <- glance_summary |> 
          #   mutate(species = i)
          # 
          # tidy_summary <- tidy_summary |> 
          #   mutate(species = i)
          # 
          # # combine both
          # total_summary <- merge(tidy_summary, glance_summary)
          # 
          # # Store the summary in the list
          # summ <- summary(gam_output)
          # 
          # dev_explained <- summ$dev.expl
          # total_summary$deviance_explained <- dev_explained
          # 
          # chi_sq <- summ[["chi.sq"]][["s(mean_last_7days)"]]
          # 
          # total_summary$chi_square <- chi_sq
          # 
          # p_value <- summ[["s.pv"]][1]
          # total_summary$p_value <- p_value
          # 
          # edf <- summ[["edf"]][1]
          # total_summary$edf <- edf
          # 
          # r_sq <- summ[["r.sq"]]
          # total_summary$r_sq <- r_sq
          # 
          # model_summary[[i]] <- total_summary
          # 
        }
        # binomial
        else {
          print(paste(data$tot_abu[1], i))
          # should be 4
          
          gam_output <- gam(data = data, Abundance ~ s(mean_last_7days, k = 3) + s(fLake, bs = 're')
                            +  s(fProtocol, bs = 're'), family = binomial)
          # prepare residuals
          simulationOutput <- simulateResiduals(fittedModel = gam_output, plot = F)
          tiff_filename <- paste("total_models/gam_check/gam_check_", i, ".tiff", sep = "")
          tiff(tiff_filename, width = 800, height = 600)
          print(plot(simulationOutput))
          dev.off()
          
          # Plotting standardized residuals against predictors
          tiff_file_2 <- paste("total_models/gam_check/predictor_", i, ".tiff", sep = "")
          tiff(tiff_file_2, width = 800, height = 600)
          print(plotResiduals(simulationOutput, temp_data$mean_last_7days, xlab = "temp", main=NULL))
          dev.off()
        
          
          model_prediction <- predict.gam(gam_output, newdata = grid,
                                          exclude = c("s(fProtocol)", "s(fLake)"),
                                          type = "response", se.fit = TRUE)
          model_bind <- cbind(grid, as.data.frame(model_prediction))
          pred_df <- model_bind |>
            rename(temp = mean_last_7days) |> 
            mutate(species = factor(i))
          saveRDS(pred_df, paste0("total_models/predictions/depth_temp/predictions_",i,".rds"))
          
          # # stats
          # glance_summary <- glance(gam_output)
          # tidy_summary <- tidy(gam_output)
          # 
          # # Add species column to both summaries
          # glance_summary <- glance_summary |> 
          #   mutate(species = i)
          # 
          # tidy_summary <- tidy_summary |> 
          #   mutate(species = i)
          # 
          # # combine both
          # total_summary <- merge(tidy_summary, glance_summary)
          # 
          # # Store the summary in the list
          # summ <- summary(gam_output)
          # 
          # dev_explained <- summ$dev.expl
          # total_summary$deviance_explained <- dev_explained
          # 
          # chi_sq <- summ[["chi.sq"]][["s(mean_last_7days)"]]
          # 
          # total_summary$chi_square <- chi_sq
          # 
          # p_value <- summ[["s.pv"]][1]
          # total_summary$p_value <- p_value
          # 
          # edf <- summ[["edf"]][1]
          # total_summary$edf <- edf
          # 
          # r_sq <- summ[["r.sq"]]
          # total_summary$r_sq <- r_sq
          # 
          # model_summary[[i]] <- total_summary
          
        }
        
      }
      
    }
    
    
  }
  # # Combine all summaries into a single data frame
  # summary_df <- bind_rows(model_summary)
  # 
  # 
  # # # Export the summary data frame to a text file
  # write_xlsx(summary_df, path = "summary_statistics.xlsx")   
  
}


stats_gamms <- function(df){
  require(broom)
  require(tidyverse)
  require(DHARMa)
  require(mgcv)
  require(gratia)
  require(mgcViz)
  require(modelsummary)
  require(writexl)
  
  species_list <- df |> 
    distinct(Species) |> 
    pull(Species)
  
  species_list <- sort(species_list)
  
  gam_output <- list()
  model_prediction <- list()
  derivatives <- list()
  grid <- list()
  pred_df <- list()
  tiff_filename <- list()
  tiff_file_2 <- list()
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
  
    
    # temp depth models
    if(max(data$tot_abu) > 20){
      # one lake
      if(max(data$n_lake) == 1) {
        
        if (i == "Coregonus_sp_benthic_profundal")  { 
          
          print(paste(data$tot_abu[1], i))
          
          gam_output <- gam(data = data, Presence ~ s(mean_last_7days, k = 3) + s(Depth_sample, k = 3)
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
        # abundance
        else if (max(data$Abundance) > 1)  {
          print(paste(data$tot_abu[1], i))
          
          gam_output <- gam(data = data, Abundance ~ s(mean_last_7days, k = 3) + s(Depth_sample, k = 3)
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
        }
        # binomial
        else{
          print(paste(data$tot_abu[1], i))
          
          gam_output <- gam(data = data, Abundance ~ s(mean_last_7days, k = 3) + s(Depth_sample, k = 3)
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
      # multiple lakes
      else{
        
      
        # special cases
        if (i %in% c("Alosa_agone", "Alburnus_arborella", "Barbatula_sp_Lineage_I",
                     "Coregonus_sp_large_pelagic", 
                     "Cyprinus_carpio", "Lepomis_gibbosus", "Phoxinus_csikii", "Salmo_trutta")){
          
          print(paste(data$tot_abu[1], i))
          
          gam_output <- gam(data = data, Presence ~ s(mean_last_7days, k = 3) + s(Depth_sample, k = 3)
                            + s(fLake, bs = 're')
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
        # abundance
        else if (max(data$Abundance) > 1)  {
          
          print(paste(data$tot_abu[1], i))
          
          gam_output <- gam(data = data, Abundance ~ s(mean_last_7days, k = 3) + s(Depth_sample, k = 3)
                            + s(fLake, bs = 're')
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
        }
        # binomial
        else {
          
          print(paste(data$tot_abu[1], i))
          
          gam_output <- gam(data = data, Abundance ~ s(mean_last_7days, k = 3) + s(Depth_sample, k = 3)
                            + s(fLake, bs = 're')
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
    # TEMPERATURE ONLY 
    else{
      
      # one lake
      if(max(data$n_lake) == 1) {
        
    
        if (i == "Coregonus_sp_benthic_profundal")  { 
          
          print(paste(data$tot_abu[1], i))
          
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
          
        }
        # abundance
        else if (max(data$Abundance) > 1)  {
          print(paste(data$tot_abu[1], i))
          
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
        # binomial
        else{
          print(paste(data$tot_abu[1], i))
          
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
      # multiple lakes
      else {
        
        # select random lake to predict them upon, random effects are excluded from
        # the prediction anyway
      
        # special cases
        if (i %in% c("Alosa_agone", "Alburnus_arborella", "Barbatula_sp_Lineage_I",
                     "Coregonus_sp_large_pelagic", 
                     "Cyprinus_carpio", "Lepomis_gibbosus", "Phoxinus_csikii", "Salmo_trutta")){
          print(i)
          
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
        # abundance
        else if (max(data$Abundance) > 1)  {
          print(paste(data$tot_abu[1], i))
          # should be one
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
          
        }
        # binomial
        else {
          print(paste(data$tot_abu[1], i))
          # should be 4
          
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
    
    
  }
  # Combine all summaries into a single data frame
  summary_df <- bind_rows(model_summary)
  
  # # Export the summary data frame to a text file
  write_xlsx(summary_df, path = "summary_stats.xlsx")   
  
}
