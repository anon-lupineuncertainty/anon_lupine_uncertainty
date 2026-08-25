# ==============================================================================
# Script: 08_Benchmarking.R
#
# Purpose: 
# This script runs a single iteration of the uncertainty analysis to benchmark
  # computation time, from bootstrap resampling through the complete uncertainty
  # analysis framework.
#
# Inputs:
# Just run this script directly after running the rest while the relevant
  # objects and functions are loaded in the environment.
#
# Notes: 
# This script requires access to the raw demographic dataset and is therefore
  # not fully reproducible from this repository alone.
# ==============================================================================


# Bootstrap resampling and model fitting ---------------------------------------

n <- 5000

sample_params_bm <- function( i, seed = F ){
  
  if( seed == T ){
    set.seed( i )
  }
  
  sample_yr <- lupine_df %>% 
    group_by( year ) %>%
    summarise( n = n() )
  
  sample_list <- vector( mode = "list", length = nrow( sample_yr ) )
  
  for( j in 1:nrow( sample_yr ) ){
    df_temp <- lupine_df[ which( lupine_df$year == sample_yr$year[j] ),]
    df_yr   <- df_temp[sample( 1:nrow( df_temp ),
                               sample_yr$n[j],
                               replace = T ),]
    sample_list[[j]] <- df_yr
  }
  
  sample_df <- bind_rows( sample_list )
  
  vr_list          <- setup_vr_list( df = sample_df )
  
  pars_temp2       <- model_vr( datlist = vr_list, mf = mf_mean2, type = "2" )

  germ_coef        <- sample_germ( df = sample_df, n = 6 )
  pars_temp2$g0    <- germ_coef$g0
  pars_temp2$g1    <- germ_coef$g1
  pars_temp2$g2    <- germ_coef$g2
  pars_temp2$g_adj <- germ_coef$g_adj

  pars_temp <- pars_temp2
  
  return( pars_temp )
}


start_boot <- Sys.time()

pars_s <- lapply( 1:n, sample_params_bm, seed = T ) %>% bind_rows()

end_boot <- Sys.time()


# Uncertainty analysis ---------------------------------------------------------

start_uncert <- Sys.time()

uncert2 <- uncertainty( ipm = lupinus_ipm2, pars = pars_var2, samples = s_pars, 
                        mega_mat = ker, vr_table = vr_tab2, bounds = bnds,
                        cores = 3 )

end_uncert <- Sys.time()


# Timing summary ---------------------------------------------------------------

boot_time <- difftime( end_boot, start_boot, units = "hours" )
uncert_time <- difftime( end_uncert, start_uncert, units = "hours" )
total_time <- boot_time + uncert_time

cat("Bootstrap/model-fitting time:", boot_time, "\n")
cat("Uncertainty-analysis time:", uncert_time, "\n")
cat("Total time:", total_time, "\n")

timing <- data.frame(
  bootstrap_hours = as.numeric(boot_time),
  uncertainty_hours = as.numeric(uncert_time),
  total_hours = as.numeric(total_time)
)

write.csv(timing, "data/benchmark_timing.csv", row.names = FALSE)