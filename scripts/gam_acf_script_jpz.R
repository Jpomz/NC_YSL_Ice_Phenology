# GAM auto-correlation plots
# Justin Pomeranz, Fall 2023

# make sure to first source the `gam_script_jpz.R`

# libraries
library(tidyverse)

# read in full, combined data
full_data <- readRDS(file = "scripts/all_lake_data.rds")

# read in model fits
ice_response_list <- readRDS(
  file = "scripts/ice_response_list.RDS")

# custom function to make single ACF plot
plot_acf <- function(model, title, lake_name = NULL){
  ACF <- acf(resid(model,
                   type = "response"),
             plot = FALSE)
  ACF <- setNames(
    data.frame(
      unclass(ACF)[c("acf", "lag")]),
    c("ACF","Lag"))
  
  plot <- ggplot(ACF, aes(x = Lag, y = ACF)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = Lag, yend = 0)) +
    labs(title = title,
         subtitle = lake_name)
  
  print(plot)
}


# function to "scale-up" `plot_acf()` to list of model fits
plot_acf_list <- function(ice_list, response_name){
  for(lake in 1:length(ice_list)){
    lake_name = names(ice_days_list)[lake]
    plot_acf(ice_list[[lake]]$mod, response_name, lake_name = lake_name)
  }
}

# duration plots ####
ice_days_list <- ice_response_list$duration
plot_acf_list(ice_days_list, "Duration")

# ice-off plots ####
ice_off_list <- ice_response_list$off_date
plot_acf_list(ice_off_list, "Ice-off-date")

# ice-on plots ####
ice_on_list <- ice_response_list$on_date
plot_acf_list(ice_on_list, "ice-on-date")