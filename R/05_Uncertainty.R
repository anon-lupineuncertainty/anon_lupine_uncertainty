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
  # - data/pars_sample2.csv: Sampled parameter values with quadratic survival
  #     model
  # - data/pars_sample3.csv: Sampled parameter values with cubic survival model
  # - data/mean_lambdas.csv: Lambda values for model comparison between
  #     alternative survival model formulations
  # - data/pars_mean2.csv: Parameters needed to construct mean IPM with
  #     quadratic survival model
  # - data/pars_mean3.csv: Parameters needed to construct mean IPM with cubic
  #     survival model
#
# Outputs:
  # - data/uncert2.rds: Uncertainty analysis output, with quadratic survival
  #     model
  # - data/uncert3.rds: Uncertainty analysis output, with cubic survival model
  # - data/uncert2_ng.rds: Uncertainty analysis output, with quadratic survival
  #     model and germination coefficients held constant
  # - data/uncert3_ng.rds: Uncertainty analysis output, with cubic survival
  #     model and germination coefficients held constant
  # - data/uncert_comp_plot.csv: Table of uncertainty contributions by vital
  #     rate, for plotting comparisons between analysis outputs when germination
  #     coefficients were resampled vs held constant
  # - data/sampled_lambdas_comp.csv: Table of lambda values, for plotting
  #     comparisons between analysis outputs when germination coefficients were
  #     resampled vs held constant
  # - data/mean_lambdas_all.csv: Updated table of mean lambda values for
  #     plotting model comparisons between alternative survival model
  #     formulations
  # - data/sampled_lambdas_mf.csv: Table of sampled lambda values, for plotting
  #     comparisons between analysis outputs of alternative survival model
  #     formulations
  # - data/var_cont.csv: Vital rate uncertainty contributions for plotting
  #     comparisons between analysis outputs of alternative survival model
  #     formulations
  # - data/mean_lambdas_ng.csv: Table of mean lambda values for plotting model
  #     comparisons between alternative survival model formulations, with
  #     germination coefficients held constant
  # - data/sampled_lambdas_ng.csv: Table of sampled lambda values, for plotting
  #     comparisons between analysis outputs of alternative survival model
  #     formulations, with germination coefficients held constant
  # - data/var_cont_ng.csv: Vital rate uncertainty contributions for plotting
  #     comparisons between analysis outputs of alternative survival model
  #     formulations, with germination coefficients held constant
  # - data/cov_plot2.csv: Sampled parameter covariances with quadratic survival
  #     model, for plotting
  # - data/cov_plot3.csv: Sampled parameter covariances with cubic survival
  #     model, for plotting
  # - data/corr_plot2_ng.csv: Sampled parameter correlations with quadratic
  #     survival model and germination coefficients held constant, for plotting
  # - data/corr_plot3_ng.csv: Sampled parameter correlations with cubic survival
  #     model and germination coefficients held constant, for plotting
#
# Notes:
# This script only requires access to the sampled and mean parameter values and
  # is therefore fully reproducible from this repository and its associated data
  # archive.
# ==============================================================================

options( stringsAsFactors = F )
library( ipmr )


# Data -------------------------------------------------------------------------

s_pars2 <- read.csv( "data/pars_sample2.csv" )
s_pars3 <- read.csv( "data/pars_sample3.csv" )

mean_lambdas <- read.csv( "data/mean_lambdas.csv" )

pars_mean2 <- read.csv( "data/pars_mean2.csv" )
pars_mean3 <- read.csv( "data/pars_mean3.csv" )


# Initializing ipmr objects ----------------------------------------------------

# Repeated exactly as in script 03_IPMs.R

inv_logit <- function( x ) exp( x )/( 1 + exp( x ) )

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

# mega_mat: Subkernel structure in row major order
  # Both forms are equivalent

ker <- "c(P + repr, SB1_germ, SB2_germ,
           enter_SB1, SB1_SB1, SB2_SB1,
           enter_SB2, SB1_SB2, SB2_SB2)"

ker <- "c( SB2_SB2, SB1_SB2, enter_SB2,
           SB2_SB1, SB1_SB1, enter_SB1,
           SB2_germ, SB1_germ, P + repr )"


bnds <- list( abort = c(0, Inf),
              clip  = c(0, Inf),
              g0    = c(0, Inf),
              g1    = c(0, Inf),
              g2    = c(0, Inf) )


# Uncertainty analysis ---------------------------------------------------------


uncert2 <- uncertainty( ipm = lupinus_ipm2, pars = pars_var2, samples = s_pars2, 
                        mega_mat = ker, vr_table = vr_tab2, bounds = bnds,
                        cores = 3 )

uncert3 <- uncertainty( ipm = lupinus_ipm3, pars = pars_var3, samples = s_pars3,
                        mega_mat = ker, vr_table = vr_tab3, bounds = bnds,
                        cores = 3 )


# Holding germination constant -------------------------------------------------

# Uncertainty analysis, holding germination coefficients constant

s_pars2_ng <- s_pars2
s_pars3_ng <- s_pars3

s_pars2_ng$g0 <- pars_mean2$g0
s_pars2_ng$g1 <- pars_mean2$g1
s_pars2_ng$g2 <- pars_mean2$g2

s_pars3_ng$g0 <- pars_mean3$g0
s_pars3_ng$g1 <- pars_mean3$g1
s_pars3_ng$g2 <- pars_mean3$g2

pars_var2_ng <- c( "surv_b0", "surv_b1", "surv_b2",
                   "grow_b0", "grow_b1", "grow_sig",
                   "abort", "clip",
                   "flow_b0", "flow_b1",
                   "fert_b0", "fert_b1" )

pars_var3_ng <- c( "surv_b0", "surv_b1", "surv_b2", "surv_b3",
                   "grow_b0", "grow_b1", "grow_sig",
                   "abort", "clip",
                   "flow_b0", "flow_b1",
                   "fert_b0", "fert_b1" )

vr_tab2_ng <- data.frame( parameter = pars_var2_ng,
                       vital_rate = c( rep( "survival", 3 ),
                                       rep( "growth", 3 ),
                                       rep( "reproduction", 6 ) ) )

vr_tab3_ng <- data.frame( parameter = pars_var3_ng,
                          vital_rate = c( rep( "survival", 4 ),
                                          rep( "growth", 3 ),
                                          rep( "reproduction", 6 ) ) )

uncert2_ng <- uncertainty( ipm = lupinus_ipm2, pars = pars_var2_ng, 
                           samples = s_pars2_ng, mega_mat = ker,
                           vr_table = vr_tab2_ng, cores = 3 )

uncert3_ng <- uncertainty( ipm = lupinus_ipm3, pars = pars_var3_ng, 
                           samples = s_pars3_ng, mega_mat = ker,
                           vr_table = vr_tab3_ng, cores = 3 )


# Formatting for plotting ------------------------------------------------------

# Comparing explained uncertainty with recruitment parameters varying and constant

uncert_comp <- uncert2$vr_uncert
uncert_comp[5,1] <- "total"
uncert_comp[5,2] <- uncert2$mod_uncert
uncert_comp$mod <- "With varying recruitment"

uncert_comp_ng <- uncert2_ng$vr_uncert
uncert_comp_ng[4,1] <- "total"
uncert_comp_ng[4,2] <- uncert2_ng$mod_uncert
uncert_comp_ng$mod <- "With constant recruitment"

uncert_comp_plot <- bind_rows( uncert_comp, uncert_comp_ng )

# Lambdas for plotting - varying vs constant recruitment parameters

lam_tab2 <- uncert2$lambdas
lam_tab2$mod <- "Varying"

lam_tab_ng <- uncert2_ng$lambdas
lam_tab_ng$mod <- "Constant"

lam_tab_2ng <- bind_rows( lam_tab2, lam_tab_ng )


# Mean of sampled lambdas for plotting 

lam_tab2 <- uncert2$lambdas

lam_tab3 <- uncert3$lambdas

mean_lam <- data.frame( type  = c( 2, 2, 3, 3 ),
                        model = c( "Mean model", "Mean of sampled models",
                                   "Mean model", "Mean of sampled models" ),
                        value = c( mean_lambdas[1,2], mean( lam_tab2$lambda ),
                                   mean_lambdas[2,2], mean( lam_tab3$lambda ) ) )

# Lambdas for plotting

lam_tab2 <- uncert2$lambdas
lam_tab2$type <- 2

lam_tab3 <- uncert3$lambdas
lam_tab3$type <- 3

lam_tab_23 <- bind_rows( lam_tab2, lam_tab3 )

# Vital rate uncertainty contributions for plotting 

var_tab2 <- uncert2$vr_uncert
var_tab2[5,] <- list( "total", uncert2$mod_uncert )
var_tab2$type <- 2

var_tab3 <- uncert3$vr_uncert
var_tab3[5,] <- list( "total", uncert3$mod_uncert )
var_tab3$type <- 3

var_tab <- bind_rows( var_tab2, var_tab3 )


# Mean of sampled lambdas for plotting - germination coefs constant

lam_tab2_ng <- uncert2_ng$lambdas
lam_tab3_ng <- uncert3_ng$lambdas

mean_lam_ng <- data.frame( type  = c( 2, 2, 3, 3 ),
                           model = c( "Mean model", "Mean of sampled models",
                                      "Mean model", "Mean of sampled models" ),
                           value = c( mean_lambdas[1,2], mean( lam_tab2_ng$lambda ),
                                      mean_lambdas[2,2], mean( lam_tab3_ng$lambda ) ) )

# Lambdas for plotting - germination coefs constant

lam_tab2_ng <- uncert2_ng$lambdas
lam_tab2_ng$type <- 2

lam_tab3_ng <- uncert3_ng$lambdas
lam_tab3_ng$type <- 3

lam_tab_ng <- bind_rows( lam_tab2_ng, lam_tab3_ng )

# Vital rate uncertainty contributions for plotting - germination coefs constant

var_tab2_ng <- uncert2_ng$vr_uncert
var_tab2_ng[4,] <- list( "total", uncert2_ng$mod_uncert )
var_tab2_ng$type <- 2

var_tab3_ng <- uncert3_ng$vr_uncert
var_tab3_ng[4,] <- list( "total", uncert3_ng$mod_uncert )
var_tab3_ng$type <- 3

var_tab_ng <- bind_rows( var_tab2_ng, var_tab3_ng )


# Parameter covariances

cov2 <- cov( s_pars2[,pars_var2] )
cov3 <- cov( s_pars3[,pars_var3] )

cov2_plot <- pivot_longer(
  tibble::rownames_to_column(
    as.data.frame(cov2),
    "Var1"
  ),
  -Var1,
  names_to = "Var2",
  values_to = "covariance"
)

cov3_plot <- pivot_longer(
  tibble::rownames_to_column(
    as.data.frame(cov3),
    "Var1"
  ),
  -Var1,
  names_to = "Var2",
  values_to = "covariance"
)


# Parameter correlation matrices - germ coefs constant

corr2_ng <- cor( s_pars2_ng[,pars_var2_ng] )
corr3_ng <- cor( s_pars3_ng[,pars_var3_ng] )

corr2_ng_plot <- pivot_longer(
  tibble::rownames_to_column(
    as.data.frame(corr2_ng),
    "Var1"
  ),
  -Var1,
  names_to = "Var2",
  values_to = "correlation"
)

corr3_ng_plot <- pivot_longer(
  tibble::rownames_to_column(
    as.data.frame(corr3_ng),
    "Var1"
  ),
  -Var1,
  names_to = "Var2",
  values_to = "correlation"
)

sig2_tbl <- cor_pmat( s_pars2_ng[,pars_var2_ng] )
sig2 <- as.matrix(sig2_tbl[, -1])
rownames(sig2) <- sig2_tbl$rowname
sig3_tbl <- cor_pmat( s_pars3_ng[,pars_var3_ng] )
sig3 <- as.matrix(sig3_tbl[, -1])
rownames(sig3) <- sig3_tbl$rowname

sig2_star <- sig2
sig3_star <- sig3

for( i in 1:length(sig2)){
  if(sig2[i] < 0.001){
    sig2_star[i] <- "***"
  } else {
    if(sig2[i] < 0.01){
      sig2_star[i] <- "**"
    } else {
      if(sig2[i] < 0.05){
        sig2_star[i] <- "*"
      } else 
        sig2_star[i] <- ""
    }
  }
}

for( i in 1:length(sig3)){
  if(sig3[i] < 0.001){
    sig3_star[i] <- "***"
  } else {
    if(sig3[i] < 0.01){
      sig3_star[i] <- "**"
    } else {
      if(sig3[i] < 0.05){
        sig3_star[i] <- "*"
      } else 
        sig3_star[i] <- ""
    }
  }
}

star2_plot <- pivot_longer(
  tibble::rownames_to_column(
    as.data.frame(sig2_star),
    "Var1"
  ),
  -Var1,
  names_to = "Var2",
  values_to = "correlation"
)

star3_plot <- pivot_longer(
  tibble::rownames_to_column(
    as.data.frame(sig3_star),
    "Var1"
  ),
  -Var1,
  names_to = "Var2",
  values_to = "correlation"
)

corr2_ng_plot$text <- star2_plot$correlation
corr3_ng_plot$text <- star3_plot$correlation

# Save output ------------------------------------------------------------------

saveRDS( uncert2, "data/uncert2.rds" )
saveRDS( uncert3, "data/uncert3.rds" )
saveRDS( uncert2_ng, "data/uncert2_ng.rds" )
saveRDS( uncert3_ng, "data/uncert3_ng.rds" )


write.csv( uncert_comp_plot, "data/uncert_comp_plot.csv", row.names = FALSE )
write.csv( lam_tab_2ng, "data/sampled_lambdas_comp.csv", row.names = FALSE )

write.csv( mean_lam, "data/mean_lambdas_all.csv", row.names = FALSE )
write.csv( lam_tab_23, "data/sampled_lambdas_mf.csv", row.names = FALSE )
write.csv( var_tab, "data/var_cont.csv", row.names = FALSE )

write.csv( mean_lam_ng, "data/mean_lambdas_ng.csv", row.names = FALSE )
write.csv( lam_tab_ng, "data/sampled_lambdas_ng.csv", row.names = FALSE )
write.csv( var_tab_ng, "data/var_cont_ng.csv", row.names = FALSE )

write.csv( cov2_plot, "data/cov_plot2.csv", row.names = FALSE )
write.csv( cov3_plot, "data/cov_plot3.csv", row.names = FALSE )

write.csv( corr2_ng_plot, "data/corr_plot2_ng.csv", row.names = F )
write.csv( corr3_ng_plot, "data/corr_plot3_ng.csv", row.names = F )

