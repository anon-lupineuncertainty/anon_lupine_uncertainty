# ==============================================================================
# Script: 02_VR_models.R
#
# Fitting vital rate models and performing model selection
#
# Purpose:
# This script takes the formatted lupine demographic data from the previous 
  # script combined with experimentally-collected recruitment data to fit vital
  # rate models and perform model selection for downstream construction of the
  # IPMs. Functions are written to extract parameter values to be used during
  # bootstrap resampling.
#
# Inputs:
  # - data/surv.csv: Full demographic dataset, subsetted for survival model
  # - data/grow.csv: Full demographic dataset, subsetted for growth model
  # - data/flow.csv: Full demographic dataset, subsetted for flowering model
  # - data/fert.csv: Full demographic dataset, subsetted for fertility model
  # - data/abort.csv: Full demographic dataset, subsetted for raceme abortion
  #     model
  # - data/cons.csv: Full demographic dataset, subsetted for raceme consumption
  #     model
  # - data/other_pars.csv: Dataframe of constant parameters for IPM construction
  # - data/fruits_per_raceme.csv: Experimentally-collected reproduction data
  # - data/seeds_per_fruit.csv: Experimentally-collected reproduction data
  # - data/seedl_size.csv: Recruitment size summary statistics
  # - data/germ_adj_factor.csv: Germination adjustment factors
#
# Outputs:
  # - data/pars_mean2.csv: Parameters needed to construct mean IPM with
  #     quadratic survival model
  # - data/pars_mean3.csv: Parameters needed to construct mean IPM with cubic
  #     survival model
  # - data/mf_mean2.csv: Model formulae of vital rates with quadratic survival
  #     model
  # - data/mf_mean3.csv: Model formulae of vital rates with cubic survival model
  # - data/other_pars.csv: Updated dataframe of constant parameters for IPM
  #     construction
  # - results/predictions/: Plots of model predictions over raw data
#
# Notes:
# This script requires access to the raw demographic dataset and is therefore
  # not fully reproducible from this repository alone.
# ==============================================================================

options( stringsAsFactors = F )
library( tidyverse )
library( bbmle )
library( patchwork )
library( mgcv )


# Data -------------------------------------------------------------------------

load_vr_dat <- function( ){
    surv       <<- read.csv( "data/surv.csv" )
    grow       <<- read.csv( "data/grow.csv" )
    flow       <<- read.csv( "data/flow.csv" )
    fert       <<- read.csv( "data/fert.csv" )
    abort      <<- read.csv( "data/abort.csv" )
    cons       <<- read.csv( "data/cons.csv" )
    other_pars <<- read.csv( "data/other_pars.csv" )
}

load_vr_dat()

fruit_rac   <- read.csv( "data/fruits_per_raceme.csv" )
seed_x_fr   <- read.csv( "data/seeds_per_fruit.csv" )
sl_size     <- read.csv( "data/seedl_size.csv" )
germ_adj    <- read.csv( "data/germ_adj_factor.csv" )

# replace 0 in seeds per fruit with very low value 0.01
seed_x_fr <- seed_x_fr %>% mutate( SEEDSPERFRUIT = replace( SEEDSPERFRUIT,
                                                            SEEDSPERFRUIT == 0,
                                                            0.01 ) )


# Functions to fit vital rate models -------------------------------------------

## Survival --------------------------------------------------------------------

fit_surv <- function( df ){
  
  surv_mod_0 <- glm( surv_t1 ~ 1, 
                     data = df, family = "binomial" )
  surv_mod_1 <- glm( surv_t1 ~ log_area_t0, 
                     data = df, family = "binomial" )
  surv_mod_2 <- glm( surv_t1 ~ log_area_t0 + log_area_t02, 
                     data = df, family = "binomial" )
  surv_mod_3 <- glm( surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03, 
                     data = df, family = "binomial" )
  
  mods <- list( surv_mod_0, surv_mod_1, surv_mod_2, surv_mod_3 )
  surv_dAIC <- AICtab( mods, weights = T, sort = F, mnames = paste0('model ', 1:length( mods ) ) )$dAIC
  
  surv_mods <- list( surv_mod_0, surv_mod_1, surv_mod_2, surv_mod_3, surv_dAIC )
  
  return( surv_mods )
}

m_surv    <- fit_surv( surv )


## Growth ----------------------------------------------------------------------

fit_grow <- function( df ){
  
  grow_mod_0 <- lm( log_area_t1 ~ 1, 
                    data = df )
  grow_mod_1 <- lm( log_area_t1 ~ log_area_t0, 
                    data = df )
  grow_mod_2 <- lm( log_area_t1 ~ log_area_t0 + log_area_t02, 
                    data = df )
  
  mods <- list( grow_mod_0, grow_mod_1, grow_mod_2 )
  grow_dAIC <- AICtab( mods, weights = T, sort = F, mnames = paste0('model ', 1:length( mods ) ) )$dAIC
  
  grow_mods <- list( grow_mod_0, grow_mod_1, grow_mod_2, grow_dAIC ) 
  
  return( grow_mods )
}

m_grow    <- fit_grow( grow )


## Reproduction ----------------------------------------------------------------

### Flowering ------------------------------------------------------------------
# Fecundity: probability of flowering

fit_flow <- function( df ){
  
  flow_mod_0 <- glm( flow_t0 ~ 1,
                     data = df, family = "binomial" )
  flow_mod_1 <- glm( flow_t0 ~ log_area_t0,
                     data = df, family = "binomial" )
  flow_mod_2 <- glm( flow_t0 ~ log_area_t0 + log_area_t02,
                     data = df, family = "binomial" )
  
  mods <- list( flow_mod_0, flow_mod_1, flow_mod_2 )
  flow_dAIC <- AICtab( mods, weights = T, sort = F, mnames = paste0('model ', 1:length( mods ) ) )$dAIC
  
  flow_mods <- list( flow_mod_0, flow_mod_1, flow_mod_2, flow_dAIC )
  
  return( flow_mods )
}

m_flow <- fit_flow( flow )


### Fertility ------------------------------------------------------------------
# Fecundity: number of racemes

fit_fert <- function( df ){
  
  fert_mod_0 <- glm( numrac_t0 ~ 1,
                     data = df, family = "poisson" )
  fert_mod_1 <- glm( numrac_t0 ~ log_area_t0,
                     data = df, family = "poisson" )
  fert_mod_2 <- glm( numrac_t0 ~ log_area_t0 + log_area_t02,
                     data = df, family = "poisson" )
  
  mods <- list( fert_mod_0, fert_mod_1, fert_mod_2 )
  fert_dAIC <- AICtab( mods, weights = T, sort = F, mnames = paste0('model ', 1:length( mods ) ) )$dAIC
  
  fert_mods <- list( fert_mod_0, fert_mod_1, fert_mod_2, fert_dAIC )
  
  return( fert_mods )
}

m_fert <- fit_fert( fert )


### Other reproductive processes -----------------------------------------------

# Abortion:

mod_abort <- glm( cbind( numab_t0, numrac_t0 - numab_t0 ) ~ 1, data = abort, 
                  family = "binomial" )

# Consumption:

mod_cons <- glm( cbind( numcl_t0, numint_t0 ) ~ 1, data = cons, 
                 family = "binomial" )

# Fruits per raceme:

m_fr_rac <- glm( NumFruits ~ 1, data = fruit_rac, family = "poisson" )

# Seeds per fruit:

m_seed_fr <- glm( SEEDSPERFRUIT ~ 1, data = seed_x_fr, 
                  family = Gamma( link = "log" ) )


# Plot vital rate models to check fit ------------------------------------------

# Function to predict from a model given a sequence of input values 

predictor_fun <- function( x, ranef ) {
  prediction <- ranef[1]
  if ( length( ranef ) >= 2 ) {
    prediction <- prediction + ranef[2] * x
  }
  if ( length( ranef ) >= 3 ) {
    prediction <- prediction + ranef[3] * x^2
  }
  if ( length( ranef ) >= 4 ) {
    prediction <- prediction + ranef[4] * x^3
  }
  
  return(prediction)
}

# Function to plot survival model predictions

plot_surv_mod <- function( df, mod, name, index ){
  
  seq_x <- seq( min( df$log_area_t0, na.rm = T ),
                max( df$log_area_t0, na.rm = T ),
                length.out = 100 )
  
  surv_pred <- predictor_fun( seq_x, coef( mod[[index]] ) ) %>%
    boot::inv.logit() %>%
    data.frame( log_area_t0 = seq_x,
                surv_t1 = . )
  
  p_surv <- plot_surv_raw( df )
  
  
  p_surv[[1]] <- p_surv[[1]] + geom_line( data = surv_pred, 
                                          aes( x = log_area_t0,
                                               y = surv_t1 ),
                                          color = "red", lwd = 2 ) +
    labs( title = name )
  
  p_surv[[2]] <- p_surv[[2]] + geom_line( data = surv_pred, 
                                          aes( x = log_area_t0,
                                               y = surv_t1 ),
                                          color = "red", lwd = 2 ) +
    labs( title = paste0( "dAIC = ", round( mod[[5]][index], 2 ) ) )
  
  ggsave( paste0( "results/predictions/surv_", 
                  index - 1, "_pred.png" ),
          p_surv, width = 8, height = 3, units = "in", dpi = 150 )
  
  return( p_surv )
}


# Function to plot growth model predictions

plot_grow_mod <- function( df, mod, name, index ){
  
  seq_x <- seq( min( df$log_area_t0, na.rm = T ),
                max( df$log_area_t0, na.rm = T ),
                length.out = 100 )
  
  grow_pred <- predictor_fun( seq_x, coef( mod[[index]] ) ) %>%
    data.frame( log_area_t0 = seq_x,
                log_area_t1 = . )
  
  p_grow <- plot_grow_raw( df )
  
  
  p_grow <- p_grow + geom_line( data = grow_pred, 
                                aes( x = log_area_t0,
                                     y = log_area_t1 ),
                                color = "red", lwd = 2 ) +
    labs( title = name,
          subtitle = paste0( "dAIC = ", round( mod[[4]][index], 2 ) ) )
  
  ggsave( paste0( "results/predictions/grow_", 
                  index - 1, "_pred.png" ),
          p_grow, width = 4, height = 3, units = "in", dpi = 150 )
  
  return( p_grow )
}


# Function to plot flowering model predictions

plot_flow_mod <- function( df, mod, name, index ){
  
  seq_x <- seq( min( df$log_area_t0, na.rm = T ),
                max( df$log_area_t0, na.rm = T ),
                length.out = 100 )
  
  flow_pred <- predictor_fun( seq_x, coef( mod[[index]] ) ) %>%
    boot::inv.logit() %>%
    data.frame( log_area_t0 = seq_x,
                flow_t0 = . )
  
  p_flow <- plot_flow_raw( df )
  
  
  p_flow[[1]] <- p_flow[[1]] + geom_line( data = flow_pred, 
                                          aes( x = log_area_t0,
                                               y = flow_t0 ),
                                          color = "red", lwd = 2 ) +
    labs( title = name )
  
  p_flow[[2]] <- p_flow[[2]] + geom_line( data = flow_pred, 
                                          aes( x = log_area_t0,
                                               y = flow_t0 ),
                                          color = "red", lwd = 2 ) +
    labs( title = paste0( "dAIC = ", round( mod[[4]][index], 2 ) ) )
  
  ggsave( paste0( "results/predictions/flow_", 
                  index - 1, "_pred.png" ),
          p_flow, width = 8, height = 3, units = "in", dpi = 150 )
  
  return( p_flow )
}


# Function to plot fertility model predictions

plot_fert_mod <- function( df, mod, name, index ){
  
  seq_x <- seq( min( df$log_area_t0, na.rm = T ),
                max( df$log_area_t0, na.rm = T ),
                length.out = 100 )
  
  fert_pred <- predictor_fun( seq_x, coef( mod[[index]] ) ) %>%
    exp() %>%
    data.frame( log_area_t0 = seq_x,
                numrac_t0 = . )
  
  p_fert <- plot_fert_raw( df )
  
  
  p_fert <- p_fert + geom_line( data = fert_pred, 
                                aes( x = log_area_t0,
                                     y = numrac_t0 ),
                                color = "red", lwd = 2 ) +
    labs( title = name,
          subtitle = paste0( "dAIC = ", round( mod[[4]][index], 2 ) ) )
  
  ggsave( paste0( "results/predictions/fert_", 
                  index - 1, "_pred.png" ),
          p_fert, width = 4, height = 3, units = "in", dpi = 150 )
  
  return( p_fert )
}

# Function to plot all model predictions

plot_mod_pred <- function( data, mod ){
  
  plot_out <- vector( mode = "list", length = ( length( mod ) - 1 ) )
  
  if( grepl( "surv", deparse( substitute( data ) ), fixed = TRUE ) ) {
    fx <- plot_surv_mod
  } else {
    if( grepl( "grow", deparse( substitute( data ) ), fixed = TRUE ) ){
      fx <- plot_grow_mod
    } else {
      if( grepl( "flow", deparse( substitute( data ) ), fixed = TRUE) ){
        fx <- plot_flow_mod
      } else 
        fx <- plot_fert_mod 
    }
  }
  mod_names <- c( "Intercept model", "Linear model", "Quadratic model", 
                  "Cubic model" )
  
  for( i in 1:( length( mod ) - 1 ) ){
    p_temp <- fx( data, mod, mod_names[i], i )
    
    plot_out[[i]] <- p_temp
    
  }
  
  return( plot_out )
}

p_surv_pred    <- plot_mod_pred( surv, m_surv )

p_grow_pred    <- plot_mod_pred( grow, m_grow )

p_flow_pred    <- plot_mod_pred( flow, m_flow )

p_fert_pred    <- plot_mod_pred( fert, m_fert )


# Model selection --------------------------------------------------------------

# Comparing dAIC values with visual inspection of model fit

## Survival --------------------------------------------------------------------

p_surv_pred

# Comparing the results of the uncertainty analysis between the quadratic and
  # cubic models will occur downstream

mod_surv2 <- m_surv[[3]]
mod_surv3 <- m_surv[[4]]


## Growth ----------------------------------------------------------------------

p_grow_pred 

# Selecting linear model

mod_grow <- m_grow[[2]]


## Reproduction ----------------------------------------------------------------

# Flowering

p_flow_pred

# Selecting linear model

mod_flow <- m_flow[[2]]


# Fertility

p_fert_pred

# Selecting linear model

mod_fert <- m_fert[[2]]


# Assembling parameters for IPMs -----------------------------------------------

other_pars$fruit_rac <- coef( m_fr_rac )[[1]] %>% exp
other_pars$seed_fruit <- coef( m_seed_fr )[[1]] %>% exp


# Two versions of the parameters and IPMs diverge at this point:
  # One version using the quadratic survival model
  # One version using the cubic survival model

# These models are differentiated by appending 2 (quadratic) or 3 (cubic) to the
  # end of object names

extr_pars <- function( type ){
  
  pars_shared <- list(
    grow_b0     = coef( mod_grow )[[1]],
    grow_b1     = coef( mod_grow )[[2]],
    grow_sig    = summary( mod_grow )$sigma,
    
    flow_b0     = coef( mod_flow )[[1]],
    flow_b1     = coef( mod_flow )[[2]],
    
    fert_b0     = coef( mod_fert )[[1]],
    fert_b1     = coef( mod_fert )[[2]],
    
    fruit_rac   = other_pars$fruit_rac,
    seed_fruit  = other_pars$seed_fruit,
    
    g0          = other_pars$g0,
    g1          = other_pars$g1,
    g2          = other_pars$g2,
    
    abort       = coef( mod_abort )[[1]] %>% boot::inv.logit(),
    clip        = coef( mod_cons )[[1]] %>% boot::inv.logit(),
    
    recr_sz     = other_pars$recr_sz,
    recr_sd     = other_pars$recr_sd,
    
    L           = range( c( grow$log_area_t0, grow$log_area_t1 ) )[1],
    U           = range( c( grow$log_area_t0, grow$log_area_t1 ) )[2],
    mat_siz     = 100
  )
  
  if( type == "3" ){
    
    pars_3 <- list(
      surv_b0     = coef( mod_surv )[[1]],
      surv_b1     = coef( mod_surv )[[2]],
      surv_b2     = coef( mod_surv )[[3]],
      surv_b3     = coef( mod_surv )[[4]]
    )
    
    pars_mean <- append( pars_3, pars_shared )
  } else {
    
    pars_2 <- list(
      surv_b0     = coef( mod_surv )[[1]],
      surv_b1     = coef( mod_surv )[[2]],
      surv_b2     = coef( mod_surv )[[3]]
    )
    
    pars_mean <- append( pars_2, pars_shared )
  }
  
  return( pars_mean )
}



# Exporting model formulae -----------------------------------------------------

# Again, we have two versions for the different survival models

extr_mod <- function( type ){
  if( type == "2" ){
    mf_mean <- list(
      mf_surv  = sub( "=*,", "", mod_surv2$call )[2],
      mf_grow  = sub( "=*,", "", mod_grow$call )[2],
      mf_flow  = sub( "=*,", "", mod_flow$call )[2],
      mf_fert  = sub( "=*,", "", mod_fert$call )[2],
      mf_abort = paste0( deparse( mod_abort$formula ) ),
      mf_cons  = paste0( deparse( mod_cons$formula ) )
    )
  } else mf_mean <- list(
    mf_surv  = sub( "=*,", "", mod_surv3$call )[2],
    mf_grow  = sub( "=*,", "", mod_grow$call )[2],
    mf_flow  = sub( "=*,", "", mod_flow$call )[2],
    mf_fert  = sub( "=*,", "", mod_fert$call )[2],
    mf_abort = paste0( deparse( mod_abort$formula ) ),
    mf_cons  = paste0( deparse( mod_cons$formula ) )
  )
  
  return( mf_mean )
}

mf_mean2 <- extr_mod( type = "2" )
mf_mean3 <- extr_mod( type = "3" )


# Fit all vital rate models and extract parameters -----------------------------

# Inputs: 
  # data = list of all dataframes to fit vital rate models on
    # (output from setup_vr_list)
  # mf = list of all model formulae to fit
  # type = which survival model to use, quadratic ("2") or cubic ("3")

# Must have models and data loaded to input shared parameters in extr_pars

model_vr <- function( datlist, mf, type ){
 
  mod_surv <- glm( formula = as.formula( mf$mf_surv ),
                   data = datlist$surv,
                   family = "binomial" )
  mod_grow <- lm( formula = as.formula( mf$mf_grow ),
                  data = datlist$grow )
  mod_flow <- glm( formula = as.formula( mf$mf_flow ),
                   data = datlist$flow,
                   family = "binomial" )
  mod_fert <- glm( formula = as.formula( mf$mf_fert ),
                   data = datlist$fert,
                   family = "poisson" )
  mod_abort <- glm( formula = as.formula( mf$mf_abort ),
                    data = datlist$abort,
                    family = "binomial" )
  mod_cons <- glm( formula = as.formula( mf$mf_cons ),
                   data = datlist$cons,
                   family = "binomial" )
  mod_list <- list( mod_surv, mod_grow, mod_flow, mod_fert, mod_abort, mod_cons )
  
  # Force 'extr_pars' to use the current local environment:
  environment( extr_pars ) <- environment()
    
  pars_out <- extr_pars( type = type )
  
  
  return( pars_out )
}


# Generating lists of parameter values

# Can use setup_vr_list and full dataframe to generate list input needed

# vr_list <- setup_vr_list( lupine_df )

# Or just put the dataframes back in a list

vr_list <- list( surv  = surv,
                 grow  = grow,
                 flow  = flow,
                 fert  = fert,
                 abort = abort,
                 cons  = cons )

pars_mean2 <- model_vr( vr_list, mf_mean2, "2" )
pars_mean3 <- model_vr( vr_list, mf_mean3, "3" )


# Save outputs -----------------------------------------------------------------

# Dataframes

write.csv( other_pars, "data/other_pars.csv", row.names = FALSE )

write.csv( pars_mean2, "data/pars_mean2.csv", row.names = FALSE )
write.csv( pars_mean3, "data/pars_mean3.csv", row.names = FALSE )

write.csv( mf_mean2, "data/mf_mean2.csv", row.names = FALSE )
write.csv( mf_mean3, "data/mf_mean3.csv", row.names = FALSE )
