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
  # - data/seedbaskets.csv: Experimentally-collected recruitment dataset
  # - data/pars_sample2.csv: Sampled parameter values with quadratic survival model
#
# Outputs:
  # - data/pars_gsim###.csv: Sampled model parameters with the specified sample
  #     size of simulated germination trials (10 files)
  # - data/g_uncert###.csv: Output of the iterative uncertainty analysis with
  #     the specified sample size of simulated germination trials (10 files)
  # - data/uncert_g_all.csv: Output of the iterative uncertainty analysis for
  #     all sample sizes, merged into one dataframe for plotting
  # - data/uncert_g_summary.csv: Summary of the iterative uncertainty analysis
  #     output for all sample sizes, for plotting
  # - data/corr_plot###.csv: Summarized parameter correlations for the specified
  #     sample size of simulated germination trials, for plotting (10 files)
#
# Notes:
# This script only requires access to the sampled parameter values and the 
  # recruitment dataset, and is therefore fully reproducible from this
  # repository and its associated data archive.
# ==============================================================================

options( stringsAsFactors = F )
library( tidyverse )
library( extraDistr )


# Data -------------------------------------------------------------------------

germ_r <- read.csv( "data/seedbaskets.csv" )
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
  
  germ_adj    <- ( germ_ii['g0'] - g_adj ) / germ_ii['g0']
  
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
  g_sim <- lapply( 1:n, germ_sim, df = germ ) %>% bind_rows
  
  # replace recruitment coefficients with sampled values drawn from simulated data
  pars_out <- lapply( 1:nrow( pars ), add_germ, 
                      n = n, pars = pars, g_sim = g_sim ) %>% bind_rows
  
  pars_out$rep <- i
  
  return( pars_out )

}

# Perform for varying sample sizes

germ6   <- lapply( 1:100, replace_germ, n = 6, pars = s_pars[1:1000,],
                   germ = germ_r, seed = T ) %>% bind_rows()
germ10  <- lapply( 1:100, replace_germ, n = 10, pars = s_pars[1:1000,],
                   germ = germ_r, seed = T ) %>% bind_rows()
germ20  <- lapply( 1:100, replace_germ, n = 20, pars = s_pars[1:1000,],
                   germ = germ_r, seed = T ) %>% bind_rows()
germ30  <- lapply( 1:100, replace_germ, n = 30, pars = s_pars[1:1000,],
                   germ = germ_r, seed = T ) %>% bind_rows()
germ40  <- lapply( 1:100, replace_germ, n = 40, pars = s_pars[1:1000,],
                   germ = germ_r, seed = T ) %>% bind_rows()
germ50  <- lapply( 1:100, replace_germ, n = 50, pars = s_pars[1:1000,],
                   germ = germ_r, seed = T ) %>% bind_rows()
germ75  <- lapply( 1:100, replace_germ, n = 75, pars = s_pars[1:1000,],
                   germ = germ_r, seed = T ) %>% bind_rows()
germ100 <- lapply( 1:100, replace_germ, n = 100, pars = s_pars[1:1000,],
                   germ = germ_r, seed = T ) %>% bind_rows()
germ200 <- lapply( 1:100, replace_germ, n = 200, pars = s_pars[1:1000,],
                   germ = germ_r, seed = T ) %>% bind_rows()
germ500 <- lapply( 1:100, replace_germ, n = 500, pars = s_pars[1:1000,],
                   germ = germ_r, seed = T ) %>% bind_rows()


# Uncertainty analysis ---------------------------------------------------------

# Initialize the ipmr IPM object (using quadratic survival model)

lupinus_ipm2 <- init_ipm( sim_gen = "general",
                          di_dd = "di",
                          det_stoch = "det",
                          kern_param = NULL ) %>% 
  define_kernel(
    name          = "P",
    formula       = s * g * d_sz,
    family        = "CC",
    s             = inv_logit( surv_b0 + ( surv_b1 * sz_1 ) + ( surv_b2 * ( sz_1 ^ 2 ) ) ),
    g             = dnorm( sz_2, g_mu, grow_sig ),
    g_mu          = grow_b0 + ( grow_b1 * sz_1 ),
    data_list     = pars_mean2,
    states        = list( c( 'sz' ) ),
    evict_cor     = FALSE
  ) %>%
  define_kernel(
    name          = "repr",
    formula       = r_r * r_s * fruit_rac * seed_fruit * ( 1 - ( abort + clip ) ) * g0 * recs * d_sz,
    family        = "CC",
    r_r           = inv_logit( flow_b0 + ( flow_b1 * sz_1 ) ),
    r_s           = exp( fert_b0 + ( fert_b1 * sz_1 ) ),
    recs          = dnorm( sz_2, recr_sz, recr_sd ),
    data_list     = pars_mean2,
    states        = list( c( 'sz' ) ),
    evict_cor     = FALSE
  ) %>%
  define_kernel(
    name          = "enter_SB1",
    formula       = f * g1 * d_sz,
    family        = "CD",
    f             = v_rac * fruit_rac * seed_fruit,
    v_rac         = tot_rac * ( 1 - ( abort + clip ) ),
    tot_rac       = inv_logit( flow_b0 + ( flow_b1 * sz_1 ) ) * exp( fert_b0 + ( fert_b1 * sz_1 ) ),
    data_list     = pars_mean2,
    states        = list( c( 'sz', 'sb1' ) ),
    evict_cor     = FALSE
  ) %>%
  define_kernel(
    name          = "enter_SB2",
    formula       = f * g2 * d_sz,
    family        = "CD",
    f             = v_rac * fruit_rac * seed_fruit,
    v_rac         = tot_rac * ( 1 - ( abort + clip ) ),
    tot_rac       = inv_logit( flow_b0 + ( flow_b1 * sz_1 ) ) * exp( fert_b0 + ( fert_b1 * sz_1 ) ),
    data_list     = pars_mean2,
    states        = list( c( 'sz', 'sb2' ) ),
    evict_cor     = FALSE
  ) %>%
  define_kernel(
    name      = "SB1_SB1",
    formula   = 0,
    family    = "DD",
    states    = list( c( 'sb1' ) ),
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name      = "SB2_SB2",
    formula   = 0,
    family    = "DD",
    states    = list( c( 'sb2' ) ),
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name          = "SB2_SB1",
    formula       = 1,
    family        = "DD",
    states        = list( c( 'sb2', 'sb1' ) ),
    evict_cor     = FALSE,
    uses_par_sets = FALSE
  ) %>%
  define_kernel(
    name      = "SB1_SB2",
    formula   = 0,
    family    = "DD",
    states    = list( c( 'sb1', 'sb2' ) ),
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name          = "SB1_germ",
    formula       = recs,
    family        = "DC",
    recs          = dnorm( sz_2, recr_sz, recr_sd ),
    data_list     = pars_mean2,
    states        = list( c( 'sb1', 'sz' ) ),
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions( 'norm', 'recs' )
  ) %>%
  define_kernel(
    name          = "SB2_germ",
    formula       = 0,
    family        = "DC",
    states        = list( c( 'sb2', 'sz' ) ),
    evict_cor     = FALSE
  ) %>%
  define_impl(
    list(
      P         = list( int_rule    = "midpoint",
                        state_start = "sz",
                        state_end   = "sz" ),
      repr      = list( int_rule    = "midpoint",
                        state_start = "sz",
                        state_end   = "sz" ),
      enter_SB1 = list( int_rule    = "midpoint",
                        state_start = "sz",
                        state_end   = "sb1" ),
      enter_SB2 = list( int_rule    = "midpoint",
                        state_start = "sz",
                        state_end   = "sb2" ),
      SB1_SB1   = list( int_rule    = "midpoint",
                        state_start = "sb1",
                        state_end   = "sb1" ),
      SB1_SB2   = list( int_rule    = "midpoint",
                        state_start = "sb1",
                        state_end   = "sb2" ),
      SB2_SB1   = list( int_rule    = "midpoint",
                        state_start = "sb2",
                        state_end   = "sb1" ),
      SB2_SB2   = list( int_rule    = "midpoint",
                        state_start = "sb2",
                        state_end   = "sb2" ),
      SB1_germ  = list( int_rule    = "midpoint",
                        state_start = "sb1",
                        state_end   = "sz" ),
      SB2_germ  = list( int_rule    = "midpoint",
                        state_start = "sb2",
                        state_end   = "sz" )
    )
  ) %>%
  define_domains(
    sz = c( pars_mean2$L, pars_mean2$U, pars_mean2$mat_siz )
  ) %>%
  define_pop_state(
    pop_vectors = list(
      n_sz  = rep( 1/100, 100 ),
      n_sb1 = 20,
      n_sb2 = 20
    )
  )

lupinus_ipm2 <- lupinus_ipm2 %>% make_ipm( iterations = 100,
                                           usr_funs = list( inv_logit = inv_logit ) )

# Setup the other necessary arguments

pars_var2   <- c( "surv_b0", "surv_b1", "surv_b2", 
                  "grow_b0", "grow_b1", "grow_sig",
                  "abort", "clip",
                  "flow_b0", "flow_b1",
                  "fert_b0", "fert_b1",
                  "g0", "g1", "g2" )

vr_tab2 <- data.frame( parameter = pars_var2,
                       vital_rate = c( rep( "survival", 3 ),
                                       rep( "growth", 3 ),
                                       rep( "reproduction", 6 ),
                                       rep( "recruitment", 3 ) ) )

ker <- c( "SB1_SB1", "SB2_SB1", "enter_SB1",
          "SB1_SB2", "SB2_SB2", "enter_SB2",
          "SB1_germ", "SB2_germ", "P", "repr" )


# Function to iteratively perform uncertainty analysis on subsets of sampled
  # parameter values

uncert_it <- function( i, df ){
  
  pars_temp <- df[which(df$rep == i),]
  
  uncert_temp <- uncertainty( ipm = lupinus_ipm2, pars = pars_var2,
                              samples = pars_temp, kernels = ker,
                              vr_table = vr_tab2, cores = 3 )
  
  uncert_out <- uncert_temp$vr_uncert
  uncert_out[5,1] <- "total"
  uncert_out[5,2] <- uncert_temp$mod_uncert
  uncert_out$rep <- i
  
  return( uncert_out )
}

g_uncert6   <- lapply( 1:100, uncert_it, germ6_c ) %>% bind_rows()
g_uncert10  <- lapply( 1:100, uncert_it, germ10_c ) %>% bind_rows()
g_uncert20  <- lapply( 1:100, uncert_it, germ20_c ) %>% bind_rows()
g_uncert30  <- lapply( 1:100, uncert_it, germ30_c ) %>% bind_rows()
g_uncert40  <- lapply( 1:100, uncert_it, germ40_c ) %>% bind_rows()
g_uncert50  <- lapply( 1:100, uncert_it, germ50_c ) %>% bind_rows()
g_uncert75  <- lapply( 1:100, uncert_it, germ75_c ) %>% bind_rows()
g_uncert100 <- lapply( 1:100, uncert_it, germ100_c ) %>% bind_rows()
g_uncert200 <- lapply( 1:100, uncert_it, germ200_c ) %>% bind_rows()
g_uncert500 <- lapply( 1:100, uncert_it, germ500_c ) %>% bind_rows()


# Wrangling for plotting

g_uncert6$sample_size <- 6
g_uncert10$sample_size <- 10
g_uncert20$sample_size <- 20
g_uncert30$sample_size <- 30
g_uncert40$sample_size <- 40
g_uncert50$sample_size <- 50
g_uncert75$sample_size <- 75
g_uncert100$sample_size <- 100
g_uncert200$sample_size <- 200
g_uncert500$sample_size <- 500

g_uncert_all <- bind_rows( g_uncert6, g_uncert10, g_uncert20, g_uncert30,
                           g_uncert40, g_uncert50, g_uncert75, g_uncert100,
                           g_uncert200, g_uncert500 )

g_uncert_summary <- g_uncert_all %>%
  group_by( sample_size, vital_rate ) %>%
  summarise( 
    mean = mean( variance_sum ),
    sd = sd( variance_sum ),
    se = sd / sqrt( n() ),
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se,
    .groups = "drop" )


# Calculate parameter correlations for plotting --------------------------------

pars_varg   <- c( "surv_b0", "surv_b1", "surv_b2", 
                  "grow_b0", "grow_b1", "grow_sig",
                  "abort", "clip",
                  "flow_b0", "flow_b1",
                  "fert_b0", "fert_b1",
                  "g0", "g1", "g2", "g_adj" )

# Function to iteratively calculate parameter correlations on subsets of sampled
  # parameter values

corr_it <- function( i, df ){
  
  pars_temp <- df[which(df$rep == i),]
  
  corr_temp <- cor( df[,c(pars_varg)] )
  
  corr_out <- pivot_longer(
    tibble::rownames_to_column(
      as.data.frame(corr_temp),
      "Var1"
    ),
    -Var1,
    names_to = "Var2",
    values_to = "correlation"
  )
  
  corr_out$rep <- i
  
  return( corr_out )
}

g_corr6   <- lapply( 1:100, corr_it, germ6_c ) %>% bind_rows()
g_corr10  <- lapply( 1:100, corr_it, germ10_c ) %>% bind_rows()
g_corr20  <- lapply( 1:100, corr_it, germ20_c ) %>% bind_rows()
g_corr30  <- lapply( 1:100, corr_it, germ30_c ) %>% bind_rows()
g_corr40  <- lapply( 1:100, corr_it, germ40_c ) %>% bind_rows()
g_corr50  <- lapply( 1:100, corr_it, germ50_c ) %>% bind_rows()
g_corr75  <- lapply( 1:100, corr_it, germ75_c ) %>% bind_rows()
g_corr100 <- lapply( 1:100, corr_it, germ100_c ) %>% bind_rows()
g_corr200 <- lapply( 1:100, corr_it, germ200_c ) %>% bind_rows()
g_corr500 <- lapply( 1:100, corr_it, germ500_c ) %>% bind_rows()


# Calculating means for plotting

corr_plot6 <- g_corr6 %>% group_by( Var1, Var2 ) %>% 
  summarise( mean_corr = mean( correlation ), .groups = "drop" )
corr_plot10 <- g_corr10 %>% group_by( Var1, Var2 ) %>% 
  summarise( mean_corr = mean( correlation ), .groups = "drop" )
corr_plot20 <- g_corr20 %>% group_by( Var1, Var2 ) %>% 
  summarise( mean_corr = mean( correlation ), .groups = "drop" )
corr_plot30 <- g_corr30 %>% group_by( Var1, Var2 ) %>% 
  summarise( mean_corr = mean( correlation ), .groups = "drop" )
corr_plot40 <- g_corr40 %>% group_by( Var1, Var2 ) %>% 
  summarise( mean_corr = mean( correlation ), .groups = "drop" )
corr_plot50 <- g_corr50 %>% group_by( Var1, Var2 ) %>% 
  summarise( mean_corr = mean( correlation ), .groups = "drop" )
corr_plot75 <- g_corr75 %>% group_by( Var1, Var2 ) %>% 
  summarise( mean_corr = mean( correlation ), .groups = "drop" )
corr_plot100 <- g_corr100 %>% group_by( Var1, Var2 ) %>% 
  summarise( mean_corr = mean( correlation ), .groups = "drop" )
corr_plot200 <- g_corr200 %>% group_by( Var1, Var2 ) %>% 
  summarise( mean_corr = mean( correlation ), .groups = "drop" )
corr_plot500 <- g_corr500 %>% group_by( Var1, Var2 ) %>% 
  summarise( mean_corr = mean( correlation ), .groups = "drop" )



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

write.csv( g_uncert6, "data/g_uncert6.csv", row.names = F )
write.csv( g_uncert10, "data/g_uncert10.csv", row.names = F )
write.csv( g_uncert20, "data/g_uncert20.csv", row.names = F )
write.csv( g_uncert30, "data/g_uncert30.csv", row.names = F )
write.csv( g_uncert40, "data/g_uncert40.csv", row.names = F )
write.csv( g_uncert50, "data/g_uncert50.csv", row.names = F )
write.csv( g_uncert75, "data/g_uncert75.csv", row.names = F )
write.csv( g_uncert100, "data/g_uncert100.csv", row.names = F )
write.csv( g_uncert200, "data/g_uncert200.csv", row.names = F )
write.csv( g_uncert500, "data/g_uncert500.csv", row.names = F )

write.csv( g_uncert_all, "data/uncert_g_all.csv", row.names = F )
write.csv( g_uncert_summary, "data/uncert_g_summary.csv", row.names = F )

write.csv( corr_plot6, "data/corr_plot6.csv", row.names = F )
write.csv( corr_plot10, "data/corr_plot10.csv", row.names = F )
write.csv( corr_plot20, "data/corr_plot20.csv", row.names = F )
write.csv( corr_plot30, "data/corr_plot30.csv", row.names = F )
write.csv( corr_plot40, "data/corr_plot40.csv", row.names = F )
write.csv( corr_plot50, "data/corr_plot50.csv", row.names = F )
write.csv( corr_plot75, "data/corr_plot75.csv", row.names = F )
write.csv( corr_plot100, "data/corr_plot100.csv", row.names = F )
write.csv( corr_plot200, "data/corr_plot200.csv", row.names = F )
write.csv( corr_plot500, "data/corr_plot500.csv", row.names = F )

