# ==============================================================================
# Script: 01_VR_formatting.R
#
# Data formatting for vital rate modelling
#
# Purpose:
# This script takes the cleaned lupine demographic data combined with
  # experimentally-collected recruitment data to construct dataframes on which
  # the vital rate models will be fit. Functions are written to perform the
  # subsetting necessary to integrate into the downstream bootstrap resampling.
#
# Inputs:
  # - data/lupine_df.csv: Full demographic dataset
  # - data/seedbaskets.csv: Experimentally-collected recruitment dataset
  # - data/seedl_size.csv: Recruitment size summary statistics
#
# Outputs:
  # - data/surv.csv: Full demographic dataset, subsetted for survival model
  # - data/grow.csv: Full demographic dataset, subsetted for growth model
  # - data/flow.csv: Full demographic dataset, subsetted for flowering model
  # - data/fert.csv: Full demographic dataset, subsetted for fertility model
  # - data/abort.csv: Full demographic dataset, subsetted for raceme abortion
  #     model
  # - data/cons.csv: Full demographic dataset, subsetted for raceme consumption
  #     model
  # - data/other_pars.csv: Dataframe of constant parameters for IPM construction
  # - results/surv_point.png: Plot of raw survival data
  # - results/grow_hist.png: Histogram of raw size data
  # - results/grow_point.png: Plot of raw growth data
  # - results/flow_point.png: Plot of raw flowering data
  # - results/fert_point.png: Plot of raw fertility data
#
# Notes:
# This script requires access to the raw demographic dataset and is therefore
  # not fully reproducible from this repository alone.
# ==============================================================================

options( stringsAsFactors = F )
library( tidyverse )
library( patchwork )
library( binom )

# Data -------------------------------------------------------------------------

lupine_df   <- read.csv( "data/lupine_df.csv" )
germ        <- read.csv( "data/seedbaskets.csv" )
sl_size     <- read.csv( "data/seedl_size.csv" )


# Functions to setup dataframes for vital rate modelling -----------------------

## Survival --------------------------------------------------------------------

# Survival of individuals to year t1 based on size in year t0

setup_vr_surv <- function( df ){
  
  surv        <- subset( df, !is.na( surv_t1 ) ) %>%
    
    subset( area_t0     != 0 ) %>%
    mutate( year         = year + 1 )
  
  return( surv )
  
}

## Growth ----------------------------------------------------------------------

# Size in year t1 based on size in year t0

setup_vr_grow <- function( df ){
  
  grow        <- df %>% 
    
    subset( !( stage_t0 %in% c( "DORM", "NF" ) ) & 
              !( stage_t1 %in% c( "D", "NF", "DORM" ) ) ) %>% 
    subset( area_t0     != 0 ) %>%
    subset( area_t1     != 0 )
  
  return( grow )
  
}

## Reproduction ----------------------------------------------------------------

# There are four processes related to reproduction/fecundity:

### Flowering ------------------------------------------------------------------

# Probability of flowering in year t1 based on size in year t0

setup_vr_flow    <- function( df ){
  
  flow           <- subset( df, !is.na( flow_t0 ) ) %>% 
    
    subset( area_t0     != 0 ) %>% 
    mutate( year         = year + 1 )
  
  return( flow )
  
}

### Fertility ------------------------------------------------------------------

# Number of reproductive racemes produced in year t1 based on size in year t0

setup_vr_fert    <- function( df ){
  
  fert           <- subset( df, flow_t0 == 1 ) %>% 
    
    subset( area_t0 != 0) %>% 
    subset( !is.na( numrac_t0 ) ) %>% 
    subset( !( flow_t0 %in% 0 ) ) %>% 
    mutate( year = year + 1 ) %>%
    subset( !( numrac_t0 %in% 0 ) )
  
  return( fert )
  
}

### Abortion -------------------------------------------------------------------

# Probability of raceme abortion
  # Note: abortion was only recorded in 2010, 2011, and from 2013 onwards

setup_vr_abort    <- function( df ){
  
  abort           <- subset( df, !is.na( flow_t0 ) & flow_t0 == 1 ) %>% 
    
    subset( area_t0 != 0 ) %>% 
    subset( !is.na( numrac_t0 ) ) %>%
    subset( !is.na( numab_t0 ) ) %>%
    subset( !( flow_t0 %in% 0 ) ) %>% 
    subset( !( numrac_t0 %in% 0 ) ) %>%
    subset( year %in% c( 2010, 2011, 2013:2018 ) )
  
  return( abort )
  
}

### Consumption ----------------------------------------------------------------

# Probability of raceme consumption ("clipped" by mice)

setup_vr_cons     <- function( df ){
  
  cons            <- subset( df, !is.na( flow_t1 ) & flow_t1 == 1 ) %>%
    
    subset( area_t1 != 0 ) %>% 
    subset( !is.na( numrac_t1 ) ) %>%
    subset( !( flow_t1 %in% 0 ) ) %>% 
    mutate( year = year + 1 ) %>% 
    subset( !( numrac_t1 %in% 0 ) )
  
  return( cons )
  
}


## All together ----------------------------------------------------------------

# Stringing all formatting functions together

setup_vr_list <- function( df ){
  
  surv_temp  <- setup_vr_surv( df )
  grow_temp  <- setup_vr_grow( df )
  flow_temp  <- setup_vr_flow( df )
  fert_temp  <- setup_vr_fert( df )
  abort_temp <- setup_vr_abort( df )
  cons_temp  <- setup_vr_cons( df )
  
  vr_out <- list( surv  = surv_temp,
                  grow  = grow_temp,
                  flow  = flow_temp,
                  fert  = fert_temp,
                  abort = abort_temp,
                  cons  = cons_temp )
  
  return( vr_out )
  
}

vr_list <- setup_vr_list( lupine_df )


# Plot dataframes to check for anomalies in the raw data -----------------------

# Function to plot binned proportions (for survival data)

plot_binned_prop <- function( df, n_bins, siz_var, rsp_var ){
  
  size_var <- deparse( substitute( siz_var ) )
  resp_var <- deparse( substitute( rsp_var ) )
  
  # remove all NAs
  na_ids   <- c( which( is.na(df[,size_var] ) ),
                 which( is.na(df[,resp_var] ) )
  ) %>% unique
  
  # remove NAs only if there is at least one NA
  if( length( na_ids ) > 0 ) df <- df[-na_ids,]
  
  # binned survival probabilities
  h    <- ( max(df[,size_var], na.rm = T ) - min( df[,size_var], na.rm = T ) ) / n_bins
  lwr  <-  min( df[,size_var], na.rm = T ) + ( h*c( 0:( n_bins - 1 ) ) )
  upr  <- lwr + h
  mid  <- lwr + ( 1/2*h )
  
  # standard error of a bernoulli process
  # https://stats.stackexchange.com/questions/82720/confidence-interval-around-binomial-estimate-of-0-or-1
  se_bern <- function( x, lwr_upr ){
    
    # do not suppress potential convergence warnings
    surv_n <- sum( x )
    tot_n  <- length( x )
    binom.confint( surv_n, tot_n, methods = c( "wilson" ) )[,lwr_upr]
    
  }
  
  binned_prop <- function( lwr_x, upr_x, response, lwr_upr ){
    
    id  <- which(df[,size_var] > lwr_x & df[,size_var] < upr_x ) 
    tmp <- df[id,]
    
    if( response == 'prob' ){   return( sum( tmp[,resp_var], na.rm = T ) / nrow( tmp ) ) }
    if( response == 'n_size' ){ return( nrow( tmp ) ) }
    if( response == 'se' ){     return( se_bern( tmp[,resp_var], lwr_upr ) ) }
    
  }
  
  y_binned <- Map( binned_prop, lwr, upr, 'prob' )        %>% unlist
  x_binned <- mid
  y_n_size <- Map( binned_prop, lwr, upr, 'n_size' )      %>% unlist
  y_se_lwr <- Map( binned_prop, lwr, upr, 'se', 'lower' ) %>% unlist
  y_se_upr <- Map( binned_prop, lwr, upr, 'se', 'upper' ) %>% unlist
  
  data.frame( x_binned, 
              y_binned,
              n_s  = y_n_size,
              lwr  = y_se_lwr,
              upr  = y_se_upr ) %>% 
    mutate( n_prob = y_n_size/sum( y_n_size ) ) %>% 
    setNames( c( size_var, resp_var, 'n_s', 'lwr', 'upr', 'n_prob' ) )
  
}


### Survival plot --------------------------------------------------------------

# Function to plot survival data

plot_surv_raw <- function( df ){
  
  p_surv_raw <- ggplot( df, aes( x = log_area_t0, y = surv_t1 ) ) +
    geom_jitter( width = 0, height = 0.2, alpha = 0.3 ) + theme_bw()
  
  p_surv_bin <- ggplot( plot_binned_prop( df, 10, log_area_t0, surv_t1 ),
                        aes( x = log_area_t0, y = surv_t1 ) ) +
    geom_point( color = "red" ) +
    geom_errorbar( aes( x = log_area_t0, ymin = lwr, ymax = upr ),
                   size = 0.5, width = 0.5 ) +
    scale_y_continuous( breaks = c( 0.1, 0.5, 0.9 ) ) +
    ylim( 0, 1 ) + theme_bw() 
  
  p_surv <- p_surv_raw + p_surv_bin
  
  return( p_surv )
}

p_surv <- plot_surv_raw( vr_list$surv )


### Growth plots ---------------------------------------------------------------

# Function to plot growth data:
  # Histograms of the log area at time t0 and time t1

plot_grow_hist <- function( df ){
  p_grow_hist_t0 <- ggplot( df, aes ( x = log_area_t0 ) ) +
    geom_histogram( binwidth = 0.25, fill = "grey", color = "black" ) +
    labs( title = "log area at t0" ) +
    theme_bw() +
    theme( plot.title = element_text( hjust = 0.5 ) )
  
  p_grow_hist_t1 <- ggplot( df, aes ( x = log_area_t1 ) ) +
    geom_histogram( binwidth = 0.25, fill = "grey", color = "black" ) +
    labs( title = "log area at t1" ) +
    theme_bw() +
    theme( plot.title = element_text( hjust = 0.5 ) )
  
  p_grow_hist <- p_grow_hist_t0 + p_grow_hist_t1
  
  return( p_grow_hist )
}

# Growth: histograms of log area at time t0 and time t1 

p_grow_hist <- plot_grow_hist( vr_list$grow )

# Function to plot growth data:
  # Scatterplots of points

plot_grow_raw <- function( df ){
  p_grow_point <- ggplot( df, aes( x = log_area_t0, y = log_area_t1 ) ) +
    geom_point( alpha = 0.3 ) + theme_bw( )
  
  return( p_grow_point )
}

# Growth: points of log area at time t0 and time t1

p_grow_point <- plot_grow_raw( vr_list$grow )


### Reproduction plots ---------------------------------------------------------

# Function to plot flowering probability data

plot_flow_raw <- function( df ){
  p_flow_raw <- ggplot( df, aes( x = log_area_t0, y = flow_t0 ) ) +
    geom_jitter( width = 0, height = 0.2, alpha = 0.3 ) + theme_bw()
  
  p_flow_bin <- ggplot( plot_binned_prop( df, 10, log_area_t0, flow_t0 ),
                        aes( x = log_area_t0, y = flow_t0 ) ) +
    geom_point( color = "red" ) +
    geom_errorbar( aes( x = log_area_t0, ymin = lwr, ymax = upr ),
                   size = 0.5, width = 0.5 ) +
    scale_y_continuous( breaks = c( 0.1, 0.5, 0.9 ) ) +
    ylim( 0, 1 ) + theme_bw() 
  
  p_flow <- p_flow_raw + p_flow_bin
  
  return( p_flow )
}

p_flow <- plot_flow_raw( vr_list$flow )


# Function to plot fertility data (number of racemes)

plot_fert_raw <- function( df ){
  p_fert_raw <- ggplot( df, aes( x = log_area_t0, y = numrac_t0 ) ) +
    geom_point( alpha = 0.3 ) + theme_bw()
  
  return( p_fert_raw )
}

p_fert <- plot_fert_raw( vr_list$fert )


# Construct dataframe of other IPM parameters ----------------------------------

# These parameters are used in the construction of the "mean" IPM

germ_coef <- select( germ, g0:g2 ) %>% colMeans

other_pars <- data.frame( mat_siz = 100,
                          g0 = germ_coef[1],
                          g1 = germ_coef[2],
                          g2 = germ_coef[3],
                          recr_sz = sl_size$mean_sl_size,
                          recr_sd = sl_size$sd_sl_size )


# Save output ------------------------------------------------------------------

# Dataframes

write.csv( vr_list$surv, "data/surv.csv", row.names = F )
write.csv( vr_list$grow, "data/grow.csv", row.names = F )
write.csv( vr_list$flow, "data/flow.csv", row.names = F )
write.csv( vr_list$fert, "data/fert.csv", row.names = F )
write.csv( vr_list$abort, "data/abort.csv", row.names = F )
write.csv( vr_list$cons, "data/cons.csv", row.names = F )

write.csv( other_pars, "data/other_pars.csv", row.names = F )


# Plots

ggsave( "results/surv_point.png", p_surv, 
        width = 8, height = 3, units = "in", dpi = 150 )

ggsave( "results/grow_hist.png", p_grow_hist, 
        width = 8, height = 3, units = "in", dpi = 150 )

ggsave( "results/grow_point.png", p_grow_point, 
        width = 4, height = 3, units = "in", dpi = 150 )

ggsave( "results/flow_point.png", p_flow, 
        width = 8, height = 3, units = "in", dpi = 150 )

ggsave( "results/fert_point.png", p_fert, 
        width = 4, height = 3, units = "in", dpi = 150 )

