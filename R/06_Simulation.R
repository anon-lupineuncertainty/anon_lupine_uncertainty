# ==============================================================================
# Script: 06_Simulation.R
#
# Simulate germination data to test uncertainty contribution
#
# Purpose:
# This script simulates germination data and compares uncertainty contributions
  # of the recruitment parameters as the sample size changes, in order to assess
  # the optimal number of germination trials to conduct in the field.
#
# Inputs:
  # - data/seedbaskets.xlsx: Experimentally-collected recruitment dataset
  # - data/pars_sample2.csv: Sampled parameter values with quadratic survival model
#
# Outputs:
# 
# 
#
# Notes:
# This script only requires access to the sampled parameter values and the 
  # recruitment dataset, and is therefore fully reproducible from this
  # repository and its associated data archive.
# ==============================================================================

options( stringsAsFactors = F )
library( readxl )
library( tidyverse )
library( extraDistr )


# Data -------------------------------------------------------------------------

germ_r <- read_xlsx( "data/seedbaskets.xlsx" )
s_pars <- read.csv( "data/pars_sample2.csv" )


# Simulating new germination data ----------------------------------------------

# Drawing from the beta binomial distribution

moment_match_beta <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return( data.frame(alpha = alpha, beta = beta) )
}

# i: number of samples to simulate
# df: raw germination data

germ_sim <- function( i, df ){
  
  moment_g0 <- moment_match_beta( mu = mean( df$g0 ), var = var( df$g0 ) )
  moment_g1 <- moment_match_beta( mu = mean( df$g1 ), var = var( df$g1 ) )
  moment_g2 <- moment_match_beta( mu = mean( df$g2 ), var = var( df$g2 ) )
  
  sim_out <- df[1,2:8]
  
  sim_out[1,2] <- rbbinom( 1, 50, alpha = moment_g0$alpha, beta = moment_g0$beta )
  sim_out[1,3] <- rbbinom( 1, 50, alpha = moment_g1$alpha, beta = moment_g1$beta )
  sim_out[1,4] <- rbbinom( 1, 50, alpha = moment_g2$alpha, beta = moment_g2$beta )
  
  sim_out[1,5] <- sim_out[1,2] / 50
  sim_out[1,6] <- sim_out[1,3] / 50
  sim_out[1,7] <- sim_out[1,4] / 50
  
  return( sim_out )
}


# Function to sample from the simulated dataframe and calculate new recruitment
  # parameters

resample_germ <- function( n, g_sim, g_adj ){
  
  germ_ii     <- g_sim[sample( 1:length( g_sim$g0 ), n, replace = T ),] %>% 
    summarise( across( where( is.numeric ), mean ) ) 
  
  germ_adj    <- ( germ_ii['g0'] - g_adj ) * germ_ii['g0']
  
  germ_coef   <- data.frame( g0 = germ_ii['g0'] * ( 1 - germ_adj ),
                             g1 = germ_ii['g1'] * ( 1 - germ_adj ),
                             g2 = germ_ii['g2'] * ( 1 - germ_adj ) )
  
  return( germ_coef )
  
}


# Function to place new recruitment parameters into the sampled parameters
  # dataframe
  # i: row to resample
  # n: number of germination trials to sample from the simulated trials
  # pars: dataframe of sampled parameter values
  # g_sim: dataframe of simulated germination trials


add_germ <- function( i, n, pars, g_sim ){
  
  pars_temp <- pars[i,]
  
  germ_adj_temp <- pars_temp$g_adj
  
  germ_temp <- resample_germ( n, g_sim, g_adj = germ_adj_temp )
  
  pars_temp$g0 <- germ_temp$g0
  pars_temp$g1 <- germ_temp$g1
  pars_temp$g2 <- germ_temp$g2
  
  return( pars_temp )
  
}


# Function to replace all germination parameters in the sampled parameter values
  # dataframe, based on simulated germination data of a specified sample size
  # i: number of times to resample
  # n: number of germination trials to simulate
  # pars: dataframe of sampled parameter values
  # germ: raw germination data to sample
  # seed: should a seed be set (default TRUE)

replace_germ <- function( i, n, pars, germ, seed = T ){
  
  if( seed == T ){
    set.seed( i )
  }
  
  # simulate germination data
  g_sim <- germ_sim( n, germ )
  
  # replace recruitment coefficients with sampled values drawn from simulated data
  pars_out <- lapply( 1:nrow( pars ), add_germ, 
                      n = n, pars = pars, g_sim = g_sim ) %>% bind_rows
  
  pars_out$rep <- i
  
  return( pars_out )

}

# Perform for varying sample sizes

germ6   <- lapply( 1:100, replace_germ, n = 6, pars = s_pars[1:1000,], germ = germ_r, seed = T )
germ10  <- lapply( 1:100, replace_germ, n = 10, pars = s_pars[1:1000,], germ = germ_r, seed = T )
germ20  <- lapply( 1:100, replace_germ, n = 20, pars = s_pars[1:1000,], germ = germ_r, seed = T )
germ30  <- lapply( 1:100, replace_germ, n = 30, pars = s_pars[1:1000,], germ = germ_r, seed = T )
germ40  <- lapply( 1:100, replace_germ, n = 40, pars = s_pars[1:1000,], germ = germ_r, seed = T )
germ50  <- lapply( 1:100, replace_germ, n = 50, pars = s_pars[1:1000,], germ = germ_r, seed = T )
germ75  <- lapply( 1:100, replace_germ, n = 75, pars = s_pars[1:1000,], germ = germ_r, seed = T )
germ100 <- lapply( 1:100, replace_germ, n = 100, pars = s_pars[1:1000,], germ = germ_r, seed = T )
germ200 <- lapply( 1:100, replace_germ, n = 200, pars = s_pars[1:1000,], germ = germ_r, seed = T )
germ500 <- lapply( 1:100, replace_germ, n = 500, pars = s_pars[1:1000,], germ = germ_r, seed = T )

# Perform uncertainty analyses on each set of parameters, reformat for plotting, export


# Save output ------------------------------------------------------------------

write.csv( germ6, "data/pars_gsim6.csv", row.names = F )
write.csv( germ10, "data/pars_gsim10.csv", row.names = F )
write.csv( germ20, "data/pars_gsim20.csv", row.names = F )
write.csv( germ30, "data/pars_gsim30.csv", row.names = F )
write.csv( germ40, "data/pars_gsim40.csv", row.names = F )
write.csv( germ50, "data/pars_gsim50.csv", row.names = F )
write.csv( germ75, "data/pars_gsim75.csv", row.names = F )
write.csv( germ100, "data/pars_gsim100.csv", row.names = F )
write.csv( germ200, "data/pars_gsim200.csv", row.names = F )
write.csv( germ500, "data/pars_gsim500.csv", row.names = F )

