# ==============================================================================
# Script: 03_IPMs.R
#
# Constructing integral projection models
#
# Purpose:
# This script takes the parameter estimates from the mean lupine model,
  # assembled in the previous script, and uses them to construct integral
  # projection models. Two IPMs are constructed: one utilizing the quadratic
  # survival model, and one utilizing the cubic survival model. Each IPM is
  # constructed as a "raw," hand-built IPM, and also using `ipmr` functionality.
  # This script also includes sensitivity and elasticity functions used in the
  # bootstrap resampling.
#
# Inputs:
  # - data/pars_mean2.csv: Parameters needed to construct mean IPM with quadratic survival model
  # - data/pars_mean3.csv: Parameters needed to construct mean IPM with cubic survival model 
#

# Notes:
# This script only requires access to the mean parameter values and is therefore
  # fully reproducible from this repository and its associated data archive.
# ==============================================================================

options( stringsAsFactors = F )
library( tidyverse )
library( ipmr )

# Data -------------------------------------------------------------------------

pars_mean2 <- read.csv( "data/pars_mean2.csv" )
pars_mean3 <- read.csv( "data/pars_mean3.csv" )


# Setup functions for hand-built IPMs ------------------------------------------

# Inverse logit

inv_logit <- function( x ) exp( x )/( 1 + exp( x ) )


## Survival --------------------------------------------------------------------

# Survival of individual of size x at time t0 to time t1
# Two versions: one for quadratic model, one for cubic model

sx2 <- function( x, pars ){
  inv_logit(   pars$surv_b0 + 
                 pars$surv_b1 * x + 
                 pars$surv_b2 * ( x^2 )
  ) 
}

sx3 <- function( x, pars ){
  inv_logit(   pars$surv_b0 + 
                 pars$surv_b1 * x + 
                 pars$surv_b2 * ( x^2 ) +
                 pars$surv_b3 * ( x^3 )
  ) 
}


## Growth ----------------------------------------------------------------------

# Growth from size x at time t0 to size y at time t1
  # Returns a probability density distribution for each x value

gxy <- function( y, x, pars ){
  dnorm( y,  mean = pars$grow_b0 + pars$grow_b1*x, 
         sd   = pars$grow_sig )
}


## Transition: survival * growth -----------------------------------------------

# Two versions, each relying on one of the specified survival functions

pxy2 <- function( y, x, pars ){
  return( sx2( x, pars ) * gxy( y, x, pars ) )
}

pxy3 <- function( y, x, pars ){
  return( sx3( x, pars ) * gxy( y, x, pars ) )
}


## Reproduction ----------------------------------------------------------------

# Fecundity: production of seeds from x-sized mothers

fx <- function( x, pars ){
  # total racemes prod
  tot_rac  <- inv_logit( pars$flow_b0 + pars$flow_b1*x ) * 
    exp(       pars$fert_b0 + pars$fert_b1*x )
  
  # viable racs
  viab_rac <- tot_rac * ( 1 - ( pars$abort + pars$clip ) )
  
  # viable seeds
  viab_sd  <- viab_rac * pars$fruit_rac * pars$seed_fruit
  
  return( viab_sd )
}

# Size distribution of recruits

recs <- function( y, pars ){
  dnorm( y, mean = pars$recr_sz, sd = pars$recr_sd )
}

# Full fecundity function

fxy <- function( y, x, pars ){
  fx( x, pars ) * recs( y, pars )
}


# IPM kernel/matrix ------------------------------------------------------------

kernel <- function( pars ){
  
  # set up IPM domain ----------------------------------------------------------
  
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U
  h <- ( U - L ) / n
  b <- L + c( 0:n ) * h
  y <- 0.5 * ( b[1:n] + b[2:( n + 1 )] )
  
  # populate kernel ------------------------------------------------------------
  
  # seeds mini matrix
  s_mat <- matrix( 0, 2, 2 )
  
  # seeds that enter 1 yr-old seed bank
  plant_s1 <- fx( y, pars ) * pars$g1
  
  # seeds that enter 2 yr-old seed bank
  plant_s2 <- fx( y, pars ) * pars$g2
  
  # seeds that go directly to seedlings germinate right away
  Fmat <- ( outer( y, y, fxy, pars ) * pars$g0 * h )
  
  # recruits from the 1 yr-old seedbank
  s1_rec <- h * recs( y, pars )
  
  # seeds that enter the 2 yr-old seed bank
  s_mat[2,1] <- 1
  
  # recruits from the 2 yr-old seedbank
  s2_rec <- numeric( n )
  
  # survival and growth of adult plants
  if( is.null( pars$surv_b3 ) ){
    Tmat <- ( outer( y, y, pxy2, pars ) * h )
  } else {
    Tmat <- ( outer( y, y, pxy3, pars ) * h )
  }

  small_K <- Tmat + Fmat
  
  # assemble the kernel --------------------------------------------------------
  
  # top two vectors
  from_plant <- rbind( rbind( plant_s2, plant_s1 ),
                       small_K )
  
  # leftmost vectors
  from_seed <- rbind( s_mat,
                      cbind( s2_rec, s1_rec ) )
  
  k_yx <- cbind( from_seed, from_plant )
  
  return( k_yx )
  
}


# Function to calculate lambda from parameters ---------------------------------

lambda_mean <- function( pars ){
  ker <- kernel( pars )
  eig <- eigen( ker )
  
  return( Re( eig$values[which.max( Mod( eig$values ) )] ) )
}

lambda_mean( pars_mean2 )
lambda_mean( pars_mean3 )


# ipmr objects -----------------------------------------------------------------

# With quadratic survival function

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

lambda( lupinus_ipm2 )
lambda_mean( pars_mean2 )


# With cubic survival function

lupinus_ipm3 <- init_ipm( sim_gen = "general",
                          di_dd = "di",
                          det_stoch = "det",
                          kern_param = NULL ) %>% 
  define_kernel(
    name          = "P",
    formula       = s * g * d_sz,
    family        = "CC",
    s             = inv_logit( surv_b0 + ( surv_b1 * sz_1 ) + ( surv_b2 * ( sz_1 ^ 2 ) ) + ( surv_b3 * ( sz_1 ^ 3 ) ) ),
    g             = dnorm( sz_2, g_mu, grow_sig ),
    g_mu          = grow_b0 + ( grow_b1 * sz_1 ),
    data_list     = pars_mean3,
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
    data_list     = pars_mean3,
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
    data_list     = pars_mean3,
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
    data_list     = pars_mean3,
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
    data_list     = pars_mean3,
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
    sz = c( pars_mean3$L, pars_mean3$U, pars_mean3$mat_siz )
  ) %>%
  define_pop_state(
    pop_vectors = list(
      n_sz  = rep( 1/100, 100 ),
      n_sb1 = 20,
      n_sb2 = 20
    )
  )

lupinus_ipm3 <- lupinus_ipm3 %>% make_ipm( iterations = 100,
                                           usr_funs = list( inv_logit = inv_logit ) )

lambda( lupinus_ipm3 )
lambda_mean( pars_mean3 )


# Compare lambdas from different IPM constructions

( ( lambda( lupinus_ipm2 ) - lambda_mean( pars_mean2 ) ) / lambda_mean( pars_mean2 ) ) * 100
( ( lambda( lupinus_ipm3 ) - lambda_mean( pars_mean3 ) ) / lambda_mean( pars_mean3 ) ) * 100

# ipmr model predicts 0.004% greater than hand-built IPM, which is tolerable


# Sensitivity & elasticity analysis --------------------------------------------

# Function to perform sensitivity analysis on hand-built IPMs
  # Needs two parameter inputs:
    # "pars", which is all of the parameters and their values
    # "pars_list", which is a vector of the names of the parameters to perturb

sens <- function( pars, pars_list, dp = 0.01 ){
  nPar <- length( pars_list )
  sPar <- numeric( nPar )
  
  for( i in 1:nPar ){
    par.now <- pars_list[i]
    m.par <- pars
    m.par[[ which( names( m.par ) == par.now ) ]] <- m.par[[ which( names( m.par ) == par.now ) ]] - dp
    lambda.down <- lambda_mean( m.par )
    m.par[[ which( names( m.par ) == par.now ) ]] <- m.par[[ which( names( m.par ) == par.now ) ]] + 2*dp
    lambda.up <- lambda_mean( m.par )
    sj <- ( lambda.up - lambda.down ) / ( 2*dp )
    sPar[i] <- sj
  }
  
  df_out <- data.frame( parameter = pars_list,
                        sensitivity = sPar,
                        elasticity = sPar * abs( as.numeric( pars[ which( names( pars ) %in% pars_list ) ] ) ) / lambda_mean( pars ) )
  
  return( df_out )
}


# Perform the analysis

pars_var2 <- c( "surv_b0", "surv_b1", "surv_b2",
                "grow_b0", "grow_b1", "grow_sig",
                "abort", "clip",
                "flow_b0", "flow_b1",
                "fert_b0", "fert_b1",
                "g0", "g1", "g2" )

pars_var3 <- c( "surv_b0", "surv_b1", "surv_b2", "surv_b3",
               "grow_b0", "grow_b1", "grow_sig",
               "abort", "clip",
               "flow_b0", "flow_b1",
               "fert_b0", "fert_b1",
               "g0", "g1", "g2" )

sens_elas2 <- sens( pars_mean2, pars_list = pars_var2 )
sens_elas3 <- sens( pars_mean3, pars_list = pars_var3 )

