################ Function GAMs
# Function to model species' abundances along temperature with generalized additive models. 
# Prediction of GAMs with mgcv and validation of model assumptions with DHARMa. 
# predictions() function is a loop to separate the species
# into the type of GAM they are predicted with based on 1) the number of lakes 
# the species occurs in, 2) the distribution of their abundance data. 

# We use 4 general types of models:
# 1) binomial with net-type as random intercept
# 2) ZIP with net-type as random intercept
# 3) binomial with net-type and lake as random intercepts
# 4) ZIP with net-type and lake as random intercepts

# and had 6 exceptions: 
# Coregonus_sp_benthic_profundal: Model did not run with ZIP and was thus modelled
# with binomial

# 5 models were simplified:
# 1 Alburnus_arborella    
# 2 Barbatula_sp_Lineage_I
# 3 Cyprinus_carpio       
# 4 Phoxinus_csikii       
# 5 Salmo_trutta  
# those 5 species have derivatives above the 90th percentile
# and modelled as binomial with net-type and lake as random intercepts


# When predicting the GAMs with random effects we followed the suggestion here: 
# https://fromthebottomoftheheap.net/2021/02/02/random-effects-in-gams/
# (see comment section)
# In the new data for the prediction a random level of the random effect is selected
# and then the effect of the random effect smooth(s) is excluded with the exclude
# argument.
# Here, "VERT" is always used as the level for the random effect of the net-type. 
# In the random effect for lakes, a random lake is chosen where the species is present
# and used in the prediction.

predictions_temp <- function(df){
  require(broom)
  require(tidyverse)
  require(DHARMa)
  require(mgcv)
  require(gratia)
  require(mgcViz)

  
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
  
  df$fLake <- as.factor(df$Lake)
  
  df$fProtocol <- as.factor(df$Protocol)
  
  for (i in species_list) {
    data <- df |> 
      filter(Species == i) |> 
      mutate(n_lake = n_distinct(Lake))
    
    unique_lakes <- unique(data$fLake)
    
    random_lake <- sample(levels(unique_lakes), 1)
    
    # species occurring in one lake -> net-type as random effect
    if(max(data$n_lake) == 1) {
      
      # newdata for all species that fall into this category -> occuring
      # in one lake only
      grid <- expand.grid(mean_last_7days = seq(
        from = min(data$mean_last_7days, na.rm = TRUE),
        to = max(data$mean_last_7days, na.rm = TRUE), by = 0.02),
        fProtocol = factor("VERT"))
      
      # special case, because model not running with ZIP
      if (i == "Coregonus_sp_benthic_profundal")  { 
        
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
        saveRDS(pred_df, paste0("total_models/predictions/temp/predictions_",i,".rds"))
        
        # additional plot with deviance explained and plot of species predictions
        
        # summary <- summary(gam_output)
        # 
        # print(signif(summary[["dev.expl"]]))
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/", i , "_temp_predictions.tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        # 
        # 
        
        # all species with abundance data, except for special case
      } else if (max(data$Abundance) > 1)  {
        
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
        
        # print(glance(gam_output[[i]]))
        
        model_prediction <- predict.gam(gam_output, newdata = grid,
                                        exclude= "s(fProtocol)",
                                        type = "response", 
                                        se.fit = TRUE)
        model_bind <- cbind(grid, as.data.frame(model_prediction))
        pred_df <- model_bind |>
          rename(temp = mean_last_7days) |>
          mutate(species = factor(i))
        
        saveRDS(pred_df, paste0("total_models/predictions/temp/predictions_",i,".rds"))
        
        # additional plot
        # summary <- summary(gam_output)
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/", i , "_temp_predictions.tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        
        
      }
      # binomial data species occuring in one lake
      else {   
        
        
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
        
        print(glance(gam_output))
        
        
        model_prediction <- predict.gam(gam_output, newdata = grid,
                                        exclude = "s(fProtocol)",
                                        type = "response",
                                        se.fit = TRUE)
        model_bind <- cbind(grid, as.data.frame(model_prediction))
        pred_df <- model_bind |>
          rename(temp = mean_last_7days) |>
          mutate(species = factor(i))
        saveRDS(pred_df, paste0("total_models/predictions/temp/predictions_",i,".rds"))
        
        # additional plot
        # summary <- summary(gam_output)
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/", i , "temp_predictions.tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        # 
        
      }
    }
    # all species occurring in multiple lakes
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
      
      # simplified models
      if (i %in% c("Alburnus_arborella", "Barbatula_sp_Lineage_I",
                   "Cyprinus_carpio", "Phoxinus_csikii", "Salmo_trutta")){
        
        
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
        
        print(glance(gam_output))
        
        model_prediction <- predict.gam(gam_output, newdata = grid,
                                        exclude = c("s(fProtocol)", "s(fLake)"),
                                        type = "response", se.fit = TRUE)
        model_bind <- cbind(grid, as.data.frame(model_prediction))
        pred_df <- model_bind |>
          rename(temp = mean_last_7days) |>
          mutate(species = factor(i))
        
        saveRDS(pred_df, paste0("total_models/predictions/temp/predictions_",i,".rds"))
        
        # additional plot
        # summary <- summary(gam_output)
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/temp_predictions_", i ,".tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        
        
      }
      # models with ZIP and occurring in multiple lakes
      else if (max(data$Abundance) > 1)  { 
        
        
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
        
        print(glance(gam_output))
        
        model_prediction <- predict.gam(gam_output, newdata = grid,
                                        exclude = c("s(fProtocol)", "s(fLake)"),
                                        type = "response", se.fit = TRUE)
        model_bind <- cbind(grid, as.data.frame(model_prediction))
        pred_df <- model_bind |>
          rename(temp = mean_last_7days) |>
          mutate(species = factor(i))
        saveRDS(pred_df, paste0("total_models/predictions/temp/predictions_",i,".rds"))
        
        # additional plot
        # summary <- summary(gam_output)
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/temp_predictions_", i ,".tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        
        # binomial in multiple lakes
      } else {
        
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
        
        print(glance(gam_output))
        
        model_prediction <- predict.gam(gam_output, newdata = grid,
                                        exclude = c("s(fProtocol)", "s(fLake)"),
                                        type = "response", se.fit = TRUE)
        model_bind <- cbind(grid, as.data.frame(model_prediction))
        pred_df <- model_bind |>
          rename(temp = mean_last_7days) |> 
          mutate(species = factor(i))
        saveRDS(pred_df, paste0("total_models/predictions/temp/predictions_",i,".rds"))
        
        # additional plot
        # summary <- summary(gam_output)
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/temp_predictions_", i ,".tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        
        
      }
      
    }
  }
}

predictions_temp <- function(df){
  require(broom)
  require(tidyverse)
  require(DHARMa)
  require(mgcv)
  require(gratia)
  require(mgcViz)

  
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
  
  df$fLake <- as.factor(df$Lake)
  
  df$fProtocol <- as.factor(df$Protocol)
  
  for (i in species_list) {
    data <- df |> 
      filter(Species == i) |> 
      mutate(n_lake = n_distinct(Lake))
    
    unique_lakes <- unique(data$fLake)
    
    random_lake <- sample(levels(unique_lakes), 1)
    
    # species occurring in one lake -> net-type as random effect
    if(max(data$n_lake) == 1) {
      
      # newdata for all species that fall into this category -> occuring
      # in one lake only
      grid <- expand.grid(mean_last_7days = seq(
        from = min(data$mean_last_7days, na.rm = TRUE),
        to = max(data$mean_last_7days, na.rm = TRUE), by = 0.02),
        fProtocol = factor("VERT"))
      
      # special case, because model not running with ZIP
      if (i == "Coregonus_sp_benthic_profundal")  { 
        
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
        saveRDS(pred_df, paste0("total_models/predictions/temp/predictions_",i,".rds"))
        
        # additional plot with deviance explained and plot of species predictions
        
        # summary <- summary(gam_output)
        # 
        # print(signif(summary[["dev.expl"]]))
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/", i , "_temp_predictions.tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        # 
        # 
        
        # all species with abundance data, except for special case
      } else if (max(data$Abundance) > 1)  {
        
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
        
        # print(glance(gam_output[[i]]))
        
        model_prediction <- predict.gam(gam_output, newdata = grid,
                                        exclude= "s(fProtocol)",
                                        type = "response", 
                                        se.fit = TRUE)
        model_bind <- cbind(grid, as.data.frame(model_prediction))
        pred_df <- model_bind |>
          rename(temp = mean_last_7days) |>
          mutate(species = factor(i))
        
        saveRDS(pred_df, paste0("total_models/predictions/temp/predictions_",i,".rds"))
        
        # additional plot
        # summary <- summary(gam_output)
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/", i , "_temp_predictions.tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        
        
      }
      # binomial data species occuring in one lake
      else {   
        
        
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
        
        print(glance(gam_output))
        
        
        model_prediction <- predict.gam(gam_output, newdata = grid,
                                        exclude = "s(fProtocol)",
                                        type = "response",
                                        se.fit = TRUE)
        model_bind <- cbind(grid, as.data.frame(model_prediction))
        pred_df <- model_bind |>
          rename(temp = mean_last_7days) |>
          mutate(species = factor(i))
        saveRDS(pred_df, paste0("total_models/predictions/temp/predictions_",i,".rds"))
        
        # additional plot
        # summary <- summary(gam_output)
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/", i , "temp_predictions.tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        # 
        
      }
    }
    # all species occurring in multiple lakes
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
      
      # simplified models
      if (i %in% c("Alburnus_arborella", "Barbatula_sp_Lineage_I",
                   "Cyprinus_carpio", "Phoxinus_csikii", "Salmo_trutta")){
        
        
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
        
        print(glance(gam_output))
        
        model_prediction <- predict.gam(gam_output, newdata = grid,
                                        exclude = c("s(fProtocol)", "s(fLake)"),
                                        type = "response", se.fit = TRUE)
        model_bind <- cbind(grid, as.data.frame(model_prediction))
        pred_df <- model_bind |>
          rename(temp = mean_last_7days) |>
          mutate(species = factor(i))
        
        saveRDS(pred_df, paste0("total_models/predictions/temp/predictions_",i,".rds"))
        
        # additional plot
        # summary <- summary(gam_output)
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/temp_predictions_", i ,".tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        
        
      }
      # models with ZIP and occurring in multiple lakes
      else if (max(data$Abundance) > 1)  { 
        
        
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
        
        print(glance(gam_output))
        
        model_prediction <- predict.gam(gam_output, newdata = grid,
                                        exclude = c("s(fProtocol)", "s(fLake)"),
                                        type = "response", se.fit = TRUE)
        model_bind <- cbind(grid, as.data.frame(model_prediction))
        pred_df <- model_bind |>
          rename(temp = mean_last_7days) |>
          mutate(species = factor(i))
        saveRDS(pred_df, paste0("total_models/predictions/temp/predictions_",i,".rds"))
        
        # additional plot
        # summary <- summary(gam_output)
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/temp_predictions_", i ,".tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        
        # binomial in multiple lakes
      } else {
        
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
        
        print(glance(gam_output))
        
        model_prediction <- predict.gam(gam_output, newdata = grid,
                                        exclude = c("s(fProtocol)", "s(fLake)"),
                                        type = "response", se.fit = TRUE)
        model_bind <- cbind(grid, as.data.frame(model_prediction))
        pred_df <- model_bind |>
          rename(temp = mean_last_7days) |> 
          mutate(species = factor(i))
        saveRDS(pred_df, paste0("total_models/predictions/temp/predictions_",i,".rds"))
        
        # additional plot
        # summary <- summary(gam_output)
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/temp_predictions_", i ,".tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        
        
      }
      
    }
  }
}

predictions <- function(df){
  require(broom)
  require(tidyverse)
  require(DHARMa)
  require(mgcv)
  require(gratia)
  require(mgcViz)
  
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
  
  df$fLake <- as.factor(df$Lake)
  
  df$fProtocol <- as.factor(df$Protocol)
  
  for (i in species_list) {
    data <- df |> 
      filter(Species == i) |> 
      mutate(n_lake = n_distinct(Lake))
    
    unique_lakes <- unique(data$fLake)
    
    random_lake <- sample(levels(unique_lakes), 1)
    
    # species occurring in one lake -> net-type as random effect
    if(max(data$n_lake) == 1) {
      
      # newdata for all species that fall into this category -> occuring
      # in one lake only
      grid <- expand.grid(mean_last_7days = seq(
        from = min(data$mean_last_7days, na.rm = TRUE),
        to = max(data$mean_last_7days, na.rm = TRUE), by = 0.02),
        fProtocol = factor("VERT"))
      
      # special case, because model not running with ZIP
      if (i == "Coregonus_sp_benthic_profundal")  { 
        
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
        saveRDS(pred_df, paste0("total_models/predictions/temp_test/predictions_",i,".rds"))
        
        # additional plot with deviance explained and plot of species predictions
        
        # summary <- summary(gam_output)
        # 
        # print(signif(summary[["dev.expl"]]))
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/", i , "_temp_predictions.tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        # 
        # 
        
        # all species with abundance data, except for special case
      } else if (max(data$Abundance) > 1)  {
        
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
        
        # print(glance(gam_output[[i]]))
        
        model_prediction <- predict.gam(gam_output, newdata = grid,
                                        exclude= "s(fProtocol)",
                                        type = "response", 
                                        se.fit = TRUE)
        model_bind <- cbind(grid, as.data.frame(model_prediction))
        pred_df <- model_bind |>
          rename(temp = mean_last_7days) |>
          mutate(species = factor(i))
        
        saveRDS(pred_df, paste0("total_models/predictions/temp_test/predictions_",i,".rds"))
        
        # additional plot
        # summary <- summary(gam_output)
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/", i , "_temp_predictions.tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        
        
      }
      # binomial data species occuring in one lake
      else {   
        
        
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
        
        print(glance(gam_output))
        
        
        model_prediction <- predict.gam(gam_output, newdata = grid,
                                        exclude = "s(fProtocol)",
                                        type = "response",
                                        se.fit = TRUE)
        model_bind <- cbind(grid, as.data.frame(model_prediction))
        pred_df <- model_bind |>
          rename(temp = mean_last_7days) |>
          mutate(species = factor(i))
        saveRDS(pred_df, paste0("total_models/predictions/temp_test/predictions_",i,".rds"))
        
        # additional plot
        # summary <- summary(gam_output)
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/", i , "temp_predictions.tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        # 
        
      }
    }
    # all species occurring in multiple lakes
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
      
      # # models with ZIP and occurring in multiple lakes
      if (max(data$Abundance) > 1)  { 
        
        
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
        
        print(glance(gam_output))
        
        model_prediction <- predict.gam(gam_output, newdata = grid,
                                        exclude = c("s(fProtocol)", "s(fLake)"),
                                        type = "response", se.fit = TRUE)
        model_bind <- cbind(grid, as.data.frame(model_prediction))
        pred_df <- model_bind |>
          rename(temp = mean_last_7days) |>
          mutate(species = factor(i))
        saveRDS(pred_df, paste0("total_models/predictions/temp_test/predictions_",i,".rds"))
        
        # additional plot
        # summary <- summary(gam_output)
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/temp_predictions_", i ,".tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        
        # binomial in multiple lakes
      } else {
        
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
        
        print(glance(gam_output))
        
        model_prediction <- predict.gam(gam_output, newdata = grid,
                                        exclude = c("s(fProtocol)", "s(fLake)"),
                                        type = "response", se.fit = TRUE)
        model_bind <- cbind(grid, as.data.frame(model_prediction))
        pred_df <- model_bind |>
          rename(temp = mean_last_7days) |> 
          mutate(species = factor(i))
        saveRDS(pred_df, paste0("total_models/predictions/temp_test/predictions_",i,".rds"))
        
        # additional plot
        # summary <- summary(gam_output)
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/temp_predictions_", i ,".tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        
        
      }
      
    }
  }
}

predictions_temp <- function(df){
  require(broom)
  require(tidyverse)
  require(DHARMa)
  require(mgcv)
  require(gratia)
  require(mgcViz)
  
  
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
  
  df$fLake <- as.factor(df$Lake)
  
  df$fProtocol <- as.factor(df$Protocol)
  
  for (i in species_list) {
    data <- df |> 
      filter(Species == i) |> 
      mutate(n_lake = n_distinct(Lake))
    
    unique_lakes <- unique(data$fLake)
    
    random_lake <- sample(levels(unique_lakes), 1)
    
    # species occurring in one lake -> net-type as random effect
    if(max(data$n_lake) == 1) {
      
      # newdata for all species that fall into this category -> occuring
      # in one lake only
      grid <- expand.grid(mean_last_7days = seq(
        from = min(data$mean_last_7days, na.rm = TRUE),
        to = max(data$mean_last_7days, na.rm = TRUE), by = 0.02),
        fProtocol = factor("VERT"))
      
      # special case, because model not running with ZIP
      if (i == "Coregonus_sp_benthic_profundal")  { 
        
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
        saveRDS(pred_df, paste0("total_models/predictions/temp/predictions_",i,".rds"))
        
        # additional plot with deviance explained and plot of species predictions
        
        # summary <- summary(gam_output)
        # 
        # print(signif(summary[["dev.expl"]]))
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/", i , "_temp_predictions.tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        # 
        # 
        
        # all species with abundance data, except for special case
      } else if (max(data$Abundance) > 1)  {
        
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
        
        # print(glance(gam_output[[i]]))
        
        model_prediction <- predict.gam(gam_output, newdata = grid,
                                        exclude= "s(fProtocol)",
                                        type = "response", 
                                        se.fit = TRUE)
        model_bind <- cbind(grid, as.data.frame(model_prediction))
        pred_df <- model_bind |>
          rename(temp = mean_last_7days) |>
          mutate(species = factor(i))
        
        saveRDS(pred_df, paste0("total_models/predictions/temp/predictions_",i,".rds"))
        
        # additional plot
        # summary <- summary(gam_output)
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/", i , "_temp_predictions.tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        
        
      }
      # binomial data species occuring in one lake
      else {   
        
        
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
        
        print(glance(gam_output))
        
        
        model_prediction <- predict.gam(gam_output, newdata = grid,
                                        exclude = "s(fProtocol)",
                                        type = "response",
                                        se.fit = TRUE)
        model_bind <- cbind(grid, as.data.frame(model_prediction))
        pred_df <- model_bind |>
          rename(temp = mean_last_7days) |>
          mutate(species = factor(i))
        saveRDS(pred_df, paste0("total_models/predictions/temp/predictions_",i,".rds"))
        
        # additional plot
        # summary <- summary(gam_output)
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/", i , "temp_predictions.tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        # 
        
      }
    }
    # all species occurring in multiple lakes
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
      
      # simplified models
      if (i %in% c("Alburnus_arborella", "Barbatula_sp_Lineage_I",
                   "Cyprinus_carpio", "Phoxinus_csikii", "Salmo_trutta")){
        
        
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
        
        print(glance(gam_output))
        
        model_prediction <- predict.gam(gam_output, newdata = grid,
                                        exclude = c("s(fProtocol)", "s(fLake)"),
                                        type = "response", se.fit = TRUE)
        model_bind <- cbind(grid, as.data.frame(model_prediction))
        pred_df <- model_bind |>
          rename(temp = mean_last_7days) |>
          mutate(species = factor(i))
        
        saveRDS(pred_df, paste0("total_models/predictions/temp/predictions_",i,".rds"))
        
        # additional plot
        # summary <- summary(gam_output)
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/temp_predictions_", i ,".tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        
        
      }
      # models with ZIP and occurring in multiple lakes
      else if (max(data$Abundance) > 1)  { 
        
        
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
        
        print(glance(gam_output))
        
        model_prediction <- predict.gam(gam_output, newdata = grid,
                                        exclude = c("s(fProtocol)", "s(fLake)"),
                                        type = "response", se.fit = TRUE)
        model_bind <- cbind(grid, as.data.frame(model_prediction))
        pred_df <- model_bind |>
          rename(temp = mean_last_7days) |>
          mutate(species = factor(i))
        saveRDS(pred_df, paste0("total_models/predictions/temp/predictions_",i,".rds"))
        
        # additional plot
        # summary <- summary(gam_output)
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/temp_predictions_", i ,".tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        
        # binomial in multiple lakes
      } else {
        
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
        
        print(glance(gam_output))
        
        model_prediction <- predict.gam(gam_output, newdata = grid,
                                        exclude = c("s(fProtocol)", "s(fLake)"),
                                        type = "response", se.fit = TRUE)
        model_bind <- cbind(grid, as.data.frame(model_prediction))
        pred_df <- model_bind |>
          rename(temp = mean_last_7days) |> 
          mutate(species = factor(i))
        saveRDS(pred_df, paste0("total_models/predictions/temp/predictions_",i,".rds"))
        
        # additional plot
        # summary <- summary(gam_output)
        # 
        # plot_pred <- pred_df |>
        #   ggplot(aes(temp, fit)) +
        #   geom_line() +
        #   geom_ribbon(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), alpha = 0.3) +
        #   theme_bw() +
        #   # facet_wrap(~fLake, scale = "free") +
        #   theme(strip.background = element_rect(fill="lightgrey")) +
        #   labs(title = paste("Species = ", i,
        #                      "deviance explained = ", signif(summary[["dev.expl"]])))
        # 
        # tiff(paste("total_models/plot_predictions/temp_predictions_", i ,".tiff", sep = ""), units="in", width=8, height=6, res=300)
        # 
        # print(plot(plot_pred))
        # 
        # dev.off()
        
        
      }
      
    }
  }
}

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
