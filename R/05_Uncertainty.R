# ==============================================================================
# Script: 05_Uncertainty.R
#
# Perform uncertainty analysis
#
# Purpose:
# This script performs the uncertainty analysis on the two sets of sampled
  # parameter values, using the `uncertainty()` function in `ipmr`.
#
# Inputs:
  # - data/pars_sample2.csv: Sampled parameter values with quadratic survival model
  # - data/pars_sample3.csv: Sampled parameter values with cubic survival model
#
# Outputs:
  # 
  # 
#
# Notes:
# This script only requires access to the sampled parameter values and is
  # therefore fully reproducible from this repository and its associated data
  # archive.
# ==============================================================================

options( stringsAsFactors = F )
library( ipmr )


# Data -------------------------------------------------------------------------

s_pars2 <- read.csv( "data/pars_sample2.csv" )
s_pars3 <- read.csv( "data/pars_sample3.csv" )


# Initializing ipmr objects ----------------------------------------------------

# Repeated exactly as in script 03_IPMs.R

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


# Other arguments --------------------------------------------------------------

# Setting up other arguments needed to perform the uncertainty analysis

# pars: Parameter names on which to perform the uncertainty analysis

pars_var2   <- c( "surv_b0", "surv_b1", "surv_b2", 
                  "grow_b0", "grow_b1", "grow_sig",
                  "abort", "clip",
                  "flow_b0", "flow_b1",
                  "fert_b0", "fert_b1",
                  "g0", "g1", "g2" )

pars_var3   <- c( "surv_b0", "surv_b1", "surv_b2", "surv_b3", 
                  "grow_b0", "grow_b1", "grow_sig",
                  "abort", "clip",
                  "flow_b0", "flow_b1",
                  "fert_b0", "fert_b1",
                  "g0", "g1", "g2" )

# vr_table: Defining which vital rates each parameter corresponds to

vr_tab2 <- data.frame( parameter = pars_var2,
                       vital_rate = c( rep( "survival", 3 ),
                                       rep( "growth", 3 ),
                                       rep( "reproduction", 6 ),
                                       rep( "recruitment", 3 ) ) )

vr_tab3 <- data.frame( parameter = pars_var3,
                       vital_rate = c( rep( "survival", 4 ),
                                       rep( "growth", 3 ),
                                       rep( "reproduction", 6 ),
                                       rep( "recruitment", 3 ) ) )

# kernels: Kernel structure in row major order

ker <- c( "SB1_SB1", "SB2_SB1", "enter_SB1",
          "SB1_SB2", "SB2_SB2", "enter_SB2",
          "SB1_germ", "SB2_germ", "P", "repr" )


# Values of zero for g2 seem to break the IPM, replace with very small number

s_pars2[which(s_pars2$g2 == 0),"g2"] <- 0.00000001
s_pars3[which(s_pars3$g2 == 0),"g2"] <- 0.00000001

# And set delta to a smaller number

del <- 0.000000001

# Uncertainty analysis ---------------------------------------------------------


uncert2 <- uncertainty( ipm = lupinus_ipm2, pars = pars_var2, samples = s_pars2, 
                        kernels = ker, vr_table = vr_tab2, delta = del, cores = 3 )
uncert3 <- uncertainty( ipm = lupinus_ipm3, pars = pars_var3, samples = s_pars3, 
                        kernels = ker, vr_table = vr_tab3, delta = del, cores = 3 )

