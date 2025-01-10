## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----get data, fig.show='hold'------------------------------------------------


# Load libraries

library(BIEN)
library(geodata)
library(ggplot2)
library(S4DM)
library(sf)
library(terra)
library(tidyterra)

# Make a temporary directory to store climate data

  temp <- tempdir()

# Get some occurrence data

  #tv <- BIEN_occurrence_species(species = "Trillium vaseyi")
  data("sample_points")
  
  
# Get environmental data
# To make things a bit faster and easier, we'll limit ourselves to the 2 variables (mean temperature and annual precipitation)

  
  # env <- worldclim_global(var = "bio",
  #                          res = 10,
  #                          path = temp)
  # env <- env[[c(1,12)]]

  env <- rast(system.file('ex/sample_env.tif', package="S4DM"))  
  
# And we'll rescale the variables as well

  env <- scale(env)

# Just to take a look to make sure we didn't mess anything up

plot(env)


## ----rangebagging, echo=FALSE, results='asis'---------------------------------

data("sample_points")

sample_points_rangebagged <- 
make_range_map(occurrences = sample_points[c("longitude","latitude")],
               env = env,
               presence_method = "rangebagging",
               background_method = "none",
               background_buffer_width = 100000)


#Lets see what it looks like

# convert to polygon for easier visualization

sample_points_rangebagged_polygon <-
  sample_points_rangebagged |>
  as.polygons() |>
  st_as_sf()

# get a bbox for plotting

sample_points_bbox <-
sample_points_rangebagged_polygon |>
  st_buffer(dist = 500000) |>
  st_bbox()

#Now, we'll plot the standardized temperature raster, along with the occurrence records and the range map 

ggplot(env)+
  geom_raster(mapping = aes(x=x,y=y,fill=wc2.1_10m_bio_1))+
  scale_fill_viridis_c(name="Temp. C", na.value = "transparent")+
  scale_x_continuous(expand=c(0,0),
                     limits = c(sample_points_bbox[1],sample_points_bbox[3]))+
    scale_y_continuous(limits = c(sample_points_bbox[2],sample_points_bbox[4]),
                     expand=c(0,0))+
  theme_bw()+
  geom_sf(data = sample_points_rangebagged_polygon,
          fill = "grey",
          size=2,
          alpha=0.5)+
    geom_point(data = sample_points,
             mapping = aes(x=longitude,y=latitude))




## -----------------------------------------------------------------------------


# Here, we'll use the same data as before for Trillium vaseyi.

#First, we'll select the background data

sample_points_bg <- get_env_bg(coords = sample_points[c("longitude","latitude")],
                    env = env,
                    width = 50000,
                    standardize = TRUE) #note that we used a small set of background points to expedite model fitting

# The returned object 'xs_bg' contains two objects:
  # 1) sample_points_bg$env a matrix of environmental covariates. This is what we need for modeling.
  # 2) sample_points_bg$bg_cells a vector containing the environmental raster cell IDs that are present in tv_bg$env. This is useful for mapping the results.

# Next, we get the presence data:

  sample_points_presence <- get_env_pres(coords = sample_points[c("longitude","latitude")],
                              env = env,
                              env_bg = sample_points_bg)

#The returned object 'tv_presence' contains two objects:

  # 1) tv_presence$env a matrix of environmental covariates. This is what we need for modeling.
  # 2) tv_presence$occurrence_sf a sf object containing the coordinate data. This is useful for conducting spatially stratified cross-validation.




# Now, we can fit the model.  Previously we used rangebagging, this time we'll use a simple KDE estimation

  sample_points_kde_kde <- fit_plug_and_play(presence = sample_points_presence$env,
                                  background = sample_points_bg$env,
                                  method = "kde")


# The object that was returned is of the class "pnp_model", which is essentially a list of model fits and associated metadata.

# To view this data on a map, we can project it to the background data we used in fitting (or we could project to a new location entirely). In either case, we use the function `project_plu_and_play`.


  sample_points_kde_kde_predictions <- project_plug_and_play(pnp_model = sample_points_kde_kde,
                                                  data = sample_points_bg$env)



#Now we can make a blank raster to store our predictions
  
  sample_points_kde_kde_raster <- env[[1]]

  values(sample_points_kde_kde_raster) <- NA

#Add our predictions to the raster

  sample_points_kde_kde_raster[sample_points_bg$bg_cells] <-
    sample_points_kde_kde_predictions

#Now, we can plot our raster

  plot(sample_points_kde_kde_raster,
       xlim = c(sample_points_bbox[1],sample_points_bbox[3]),
       ylim = c(sample_points_bbox[2],sample_points_bbox[4]))
  points(sample_points[c("longitude","latitude")])
  

## ----thresholding-------------------------------------------------------------

#To threshold this continuous raster to yield a binary raster

sample_points_kde_kde_raster <- sdm_threshold(prediction_raster = sample_points_kde_kde_raster,
                                   occurrence_sf = sample_points_presence$occurrence_sf,
                                   quantile = 0.05,
                                   return_binary = T)


# As before, we'll plot this on top of temperature and occurrence records to see how well we did

  # convert to polygon for easier visualization
  
    sample_points_kde_kde_polygon <-
      sample_points_kde_kde_raster |>
      as.polygons()|>
      st_as_sf()

# Now, we'll plot the standardized temperature raster, along with the occurrence records and the range map 

  ggplot(env)+
    geom_raster(mapping = aes(x=x,y=y,fill=wc2.1_10m_bio_1))+
    scale_fill_viridis_c(name="Temp. C", na.value = "transparent")+
    scale_x_continuous(expand=c(0,0),
                       limits = c(sample_points_bbox[1],
                                  sample_points_bbox[3]))+
      scale_y_continuous(limits = c(sample_points_bbox[2],
                                    sample_points_bbox[4]),
                       expand=c(0,0))+
    theme_bw()+
    geom_sf(data = sample_points_kde_kde_polygon,
            fill = "grey",
            size=2,
            alpha=0.5)+
      geom_point(data = sample_points,
               mapping = aes(x=longitude,y=latitude))




## -----------------------------------------------------------------------------

# We'll rely on the same data as last time for simplicity.
# Since we're using different methods for estimating the presence and background distributions, we need to specify these separately:

sample_points_gaussian_kde <- fit_plug_and_play(presence = sample_points_presence$env,
                                     background = sample_points_bg$env,
                                     presence_method = "gaussian",
                                     background_method = "kde")

sample_points_gaussian_kde_predictions <- project_plug_and_play(pnp_model = sample_points_gaussian_kde,
                                                data = sample_points_bg$env)

# Now, we again convert everything to a raster and then to a polygon

  sample_points_gaussian_kde_raster <- env[[1]]

  values(sample_points_gaussian_kde_raster) <- NA

  sample_points_gaussian_kde_raster[sample_points_bg$bg_cells] <-  sample_points_gaussian_kde_predictions

# Now, we can plot our raster
  
  plot(sample_points_gaussian_kde_raster,
       xlim = c(sample_points_bbox[1],sample_points_bbox[3]),
       ylim = c(sample_points_bbox[2],sample_points_bbox[4]))
  points(sample_points[c("longitude","latitude")])
  

## ----thresholding gk----------------------------------------------------------
# To threshold this continuous raster to yield a binary raster

  sample_points_gaussian_kde_raster <- sdm_threshold(prediction_raster =
                                                       sample_points_gaussian_kde_raster,
                                     occurrence_sf = sample_points_presence$occurrence_sf,
                                     quantile = 0.05,
                                     return_binary = T)


# As before, we'll plot this on top of temperature and occurrence records to see how well we did


# Convert the raster to a polygon for visualization

  # convert to polygon for easier visualization
  
    sample_points_gaussian_kde_polygon <-
      sample_points_gaussian_kde_raster |>
      as.polygons()|>
      st_as_sf()

# Now, we'll plot the standardized temperature raster, along with the occurrence records and the range map 

  ggplot(env)+
    geom_raster(mapping = aes(x = x, y = y, fill = wc2.1_10m_bio_1))+
    scale_fill_viridis_c(name="Temp. C", na.value = "transparent")+
    scale_x_continuous(expand=c(0,0),
                       limits = c(sample_points_bbox[1],
                                  sample_points_bbox[3]))+
      scale_y_continuous(limits = c(sample_points_bbox[2],
                                    sample_points_bbox[4]),
                       expand=c(0,0))+
    theme_bw()+
    geom_sf(data = sample_points_gaussian_kde_polygon,
            fill = "grey",
            size=2,
            alpha=0.5)+
      geom_point(data = sample_points,
               mapping = aes(x=longitude,y=latitude))



## ----maxnet-------------------------------------------------------------------


  sample_points_maxnet <-
  fit_density_ratio(presence = sample_points_presence$env,
                    background = sample_points_bg$env,
                    method = "maxnet")
  
  
  sample_points_maxnet_predictions <-
    project_density_ratio(dr_model = sample_points_maxnet,
                          data = sample_points_bg$env)


#Now, we again convert everything to a raster and then to a polygon
  
  sample_points_maxnet_raster <- env[[1]]

  values(sample_points_maxnet_raster) <- NA

  sample_points_maxnet_raster[sample_points_bg$bg_cells] <-
    sample_points_maxnet_predictions


#Now, we can plot our raster
  
    plot(sample_points_maxnet_raster,
       xlim = c(sample_points_bbox[1],
                sample_points_bbox[3]),
       ylim = c(sample_points_bbox[2],
                sample_points_bbox[4]))
  points(sample_points[c("longitude","latitude")])





## ----ulsif--------------------------------------------------------------------


  sample_points_ulsif <-
  fit_density_ratio(presence = sample_points_presence$env,
                    background = sample_points_bg$env,
                    method = "ulsif")
  
  
  sample_points_ulsif_predictions <-
    project_density_ratio(dr_model = sample_points_ulsif,
                          data = sample_points_bg$env)


#Now, we again convert everything to a raster and then to a polygon
  
  sample_points_ulsif_raster <- env[[1]]

  values(sample_points_ulsif_raster) <- NA

  sample_points_ulsif_raster[sample_points_bg$bg_cells] <-
    sample_points_ulsif_predictions


#Now, we can plot our raster
  
    plot(sample_points_ulsif_raster,
       xlim = c(sample_points_bbox[1],
                sample_points_bbox[3]),
       ylim = c(sample_points_bbox[2],
                sample_points_bbox[4]))
  points(sample_points[c("longitude","latitude")])


## ----model evaluation---------------------------------------------------------

sample_points_gaussian_gaussian_fit <-
  evaluate_range_map(occurrences = sample_points[c("longitude","latitude")],
                     env = env,
                     presence_method = "gaussian",
                     background_method = "gaussian")


#Rather than looking at all of the results, we'll focus on just a few:

sample_points_gaussian_gaussian_fit$fold_results[c('testing_AUC','testing_sensitivity','testing_specificity')]

#The AUC gives us an overall idea of the discriminatory ability of the model, while the sensitivity and specificity tell us how well it discriminates presence vs. background points (respectively).



