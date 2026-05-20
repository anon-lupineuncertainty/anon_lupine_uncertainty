# ==============================================================================
# Script: 04_Sampling.R
#
# Bootstrap resampling demographic data and refitting vital rate models
#
# Purpose:
# This script samples from the raw lupine demographic data and refits all vital
  # rate models to the subsetted data, outputting a dataframe of sampled vital
  # rate parameter values. 
#
# Inputs:
  # - data/lupine_df.csv: Full demographic dataset
  # - data/seedbaskets.xlsx: Experimentally-collected recruitment dataset
  # - data/pars_mean2.csv: Parameters needed to construct mean IPM with quadratic survival model
  # - data/pars_mean3.csv: Parameters needed to construct mean IPM with cubic survival model
  # - data/mf_mean2.csv: Model formulae of vital rates with quadratic survival model
  # - data/mf_mean3.csv: Model formulae of vital rates with cubic survival model
#
# Outputs:
  # - data/pars_sample2.csv: Sampled parameter values with quadratic survival model
  # - data/pars_sample3.csv: Sampled parameter values with cubic survival model
#
# Notes:
# This script requires access to the raw demographic dataset and is therefore
  # not fully reproducible from this repository alone.
# ==============================================================================

options( stringsAsFactors = F )
library( tidyverse )


# Data -------------------------------------------------------------------------

lupine_df    <- read.csv( "data/lupine_df.csv" )
germ         <- read_xlsx( "data/seedbaskets.xlsx" ) %>% 
  select( g0:g2 )
pars_mean2    <- read.csv( "data/pars_mean2.csv" )
pars_mean3    <- read.csv( "data/pars_mean3.csv" )
mf_mean2      <- read.csv( "data/mf_mean2.csv" )
mf_mean3      <- read.csv( "data/mf_mean3.csv" )


# Sample germination data ------------------------------------------------------

# A separate pipeline samples the germination data independently of the rest of
  # the demographic data, as these data were collected independently

# Define the number of samples to take

n <- 5000

# Calculate the average number of seeds produced by each reproductive raceme
  # (a constant value in this model)

seeds_per_raceme <- pars_mean2$fruit_rac * pars_mean2$seed_fruit


# Function to recalculate germination adjustment factor based on sampled data

calc_germ_adj <- function( df, g0 ){
  
  seed_raceme <- seeds_per_raceme
  
  # calculate number of seeds produced at t0
  
  rac_yr <- df %>% 
    subset( year %in% c( 2010, 2011, 2013:2017 ) ) %>%
    summarise( racemes = sum( numrac_t0, na.rm = T ),
               clipped = sum( numcl_t0, na.rm = T ),
               aborted = sum( numab_t0, na.rm = T ) ) %>%
    mutate( num_int = racemes - ( aborted + clipped ) ) %>%
    mutate( seeds = round( num_int * seed_raceme, 0 ) )
  
  # calculate number of new seedlings at t1
  
  sdl_yr <- df %>% 
    subset( year %in% c( 2011, 2012, 2014:2018 ) ) %>%
    subset( !is.na( stage_t0 ) ) %>% 
    subset( stage_t0 == 'SL' ) %>% 
    summarise( seedlings = n() )
  
  # calculate ratio of observed seedlings / estimated seeds
  
  ratio <- sdl_yr$seedlings / rac_yr$seeds
  
  # calculate germination adjustment factor based on germination field data
  
  germ_adj <- ( g0 - ratio ) / g0
  
  return( germ_adj )
  
}

# Function to draw a sample from the germination adjustment factors

sample_germ <- function( df, n ){
  
  germ_ii     <- germ[sample( 1:6, n, replace = T ),] %>% 
    summarise( across( where( is.numeric ), mean ) ) 
  
  germ_adj    <- calc_germ_adj( df, g0 = germ_ii$g0 )
  
  germ_coef   <- data.frame( g0 = germ_ii['g0'] * ( 1 - germ_adj ),
                             g1 = germ_ii['g1'] * ( 1 - germ_adj ),
                             g2 = germ_ii['g2'] * ( 1 - germ_adj ),
                             g_adj = germ_ii['g0'] - ( germ_ii['g0'] * germ_adj ) )
  
  return( germ_coef )
  
}


# Sample demographic data ------------------------------------------------------

# Subsample, with replacement, stratifying by year, specifying the proportion of
  # samples to select per year; refit vital rate models, output parameters, 
  # sensitivities, and lambdas for each sample
  # Inputs:
    # i: number of times to resample
    # seed: whether a seed should be set while sampling (default FALSE)

sample_params <- function( i, seed = F ){
  
  if( seed == T ){
    set.seed( i )
  }
  
  # create dataframe of number of samples to draw per year
  sample_yr <- lupine_df %>% 
    group_by( year ) %>%
    summarise( n = n() )
  
  # sample stratified by year with replacement
  sample_list <- vector( mode = "list", length = nrow( sample_yr ) )
  
  for( j in 1:nrow( sample_yr ) ){
    df_temp <- lupine_df[ which( lupine_df$year == sample_yr$year[j] ),]
    df_yr   <- df_temp[sample( 1:nrow( df_temp ),
                               sample_yr$n[j],
                               replace = T ),]
    sample_list[[j]] <- df_yr
  }
  
  sample_df <- bind_rows( sample_list )
  
  # set up vital rate dataframes
  vr_list          <- setup_vr_list( df = sample_df )
  
  # calculate parameter estimates for sampled dataset
  pars_temp2       <- model_vr( datlist = vr_list, mf = mf_mean2, type = "2" )
  pars_temp3       <- model_vr( datlist = vr_list, mf = mf_mean3, type = "3" )
  
  # sample from germination data and calculate new germination coefficients
  germ_coef        <- sample_germ( df = sample_df, n = 6 )
  pars_temp2$g0    <- germ_coef$g0
  pars_temp2$g1    <- germ_coef$g1
  pars_temp2$g2    <- germ_coef$g2
  pars_temp2$g_adj <- germ_coef$g_adj
  pars_temp3$g0    <- germ_coef$g0
  pars_temp3$g1    <- germ_coef$g1
  pars_temp3$g2    <- germ_coef$g2
  pars_temp3$g_adj <- germ_coef$g_adj
  
  # add suffixes to differentiate parameter values between IPMs
  names( pars_temp2 ) <- paste0( names( pars_temp2 ), "_2" )
  names( pars_temp3 ) <- paste0( names( pars_temp3 ), "_3" )
  
  # append both sets of parameters together to output a single list
  pars_temp <- append( pars_temp2, pars_temp3 )
  
  return( pars_temp )
}

pars_s <- lapply( 1:n, sample_params, seed = T ) %>% bind_rows()

# Split the output and remove the temporary suffixes

pars_s2 <- pars_s[,1:23]
pars_s3 <- pars_s[,24:46]

names( pars_s2 ) <- gsub( "_2$", "", names( pars_s2 ) )
names( pars_s3 ) <- gsub( "_3$", "", names( pars_s3 ) )


# Save output ------------------------------------------------------------------

# Dataframes

write.csv( pars_s2, "data/pars_sample2.csv", row.names = F )
write.csv( pars_s3, "data/pars_sample3.csv", row.names = F )

