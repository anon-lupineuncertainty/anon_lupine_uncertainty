# ==============================================================================
# Script: 07_Figures.R
#
# Prepare figures
#
# Purpose:
# This script prepares figures for inclusion in the manuscript.
#
# Inputs:
  # - data/surv.csv: Full demographic dataset, subsetted for survival model
  # - data/uncert2.rds: Output of uncertainty analysis using quadratic survival model
  # - data/uncert3.rds: Output of uncertainty analysis using cubic survival model
#
# Outputs:
# 
# 
#
# Notes:
# This script primarily requires access to the objects produced by the previous
  # scripts and is therefore mostly reproducible from this repository and its
  # associated data archive. Components of Figure 2 are not reproducible, as
  # they include plots of the raw demographic data.
# ==============================================================================

options( stringsAsFactors = F )
library( ipmr )
library( ggplot2 )
library( patchwork )


# Data -------------------------------------------------------------------------

surv    <- read.csv( "data/surv.csv" )

uncert2 <- readRDS( "data/uncert2.rds" )
uncert3 <- readRDS( "data/uncert3.rds" )




# Figure 1 ---------------------------------------------------------------------

fig1a <- plot( uncert2, type = "param" )

fig1b <- plot( uncert2, type = "vr" )


# Figure 2 ---------------------------------------------------------------------

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

surv_binned <- plot_binned_prop( surv, 10, log_area_t0, surv_t1 )

fig2_base <- ggplot() + 
  geom_jitter( data = surv, aes( x = log_area_t0, y = surv_t1 ),
               width = 0.05, height = 0.1, alpha = 0.1 ) +
  geom_errorbar( data = surv_binned, 
                 aes( x = log_area_t0, ymin = lwr, ymax = upr ),
                 size = 0.5, width = 0.2 ) +
  geom_point( data = surv_binned,
              aes( x = log_area_t0, y = surv_t1 ), color = "red", cex = 2 ) +
  ylim( -0.1, 1.1 ) +
  labs( x = "log( Size at time t0)",
        y = "Probability of survival to time t1" ) +
  theme_bw( )

surv_mod2 <- glm( surv_t1 ~ log_area_t0 + log_area_t02, 
                   data = surv, family = "binomial" )
surv_mod3 <- glm( surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03, 
                   data = surv, family = "binomial" )

model_AIC  <- data.frame( model = c( 2, 3 ), size = c( 0.5, 0.5 ), 
                          survival = c( 0.75, 0.75 ), 
                          AIC = c( paste0( "AIC = ", round( AIC( surv_mod2 ), 2 ) ), 
                                   paste0( "AIC = ", round( AIC( surv_mod3 ), 2 ) ) ) )

seq_x <- seq( min( surv$log_area_t0, na.rm = T ),
              max( surv$log_area_t0, na.rm = T ),
              length.out = 100 )

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

surv_pred2 <- predictor_fun( seq_x, coef( surv_mod2 ) ) %>%
  boot::inv.logit() %>%
  data.frame( log_area_t0 = seq_x,
              surv_t1 = . )

surv_pred3 <- predictor_fun( seq_x, coef( surv_mod3 ) ) %>%
  boot::inv.logit() %>%
  data.frame( log_area_t0 = seq_x,
              surv_t1 = . )

# Figure 2a

fig2a <- fig2_base +
  geom_text( data = model_AIC[which( model_AIC$model == 2 ),], 
             aes( x = size, y = survival, label = AIC ) ) +
  geom_line( data = surv_pred2, 
             aes( x = log_area_t0,
                  y = surv_t1 ),
             color = "red", lwd = 1.5, alpha = 0.8 ) +
  labs( title = "(a) Quadratic model" )

# Figure 2b

fig2b <- fig2_base +
  geom_text( data = model_AIC[which( model_AIC$model == 3 ),], 
             aes( x = size, y = survival, label = AIC ) ) +
  geom_line( data = surv_pred3, 
             aes( x = log_area_t0,
                  y = surv_t1 ),
             color = "red", lwd = 1.5, alpha = 0.8 ) +
  labs( title = "(b) Cubic model" )

# Figure 2c

model_labs <- c( "Quadratic", "Cubic" )

var_tab2 <- uncert2$vr_uncert
var_tab2[5,] <- c( "total", uncert2$mod_uncert )
var_tab2$type <- "quadratic"

var_tab3 <- uncert3$vr_uncert
var_tab3[5,] <- c( "total", uncert3$mod_uncert )
var_tab3$type <- "cubic"

var_tab <- bind_rows( var_tab2, var_tab3 )


fig2c <- ggplot( ) + 
  geom_bar( data = var_tab[which(var_tab$vital_rate != "total" ),], 
            aes( x = type, y = variance_sum, fill = vital_rate ), 
            stat = "identity", position = "stack" ) +
  geom_bar( data = var_tab[which(var_tab$vital_rate == "total" ),], 
            aes( x = type, y = variance_sum ), 
            stat = "identity", position = "stack", color = "black", fill = NA ) +
  scale_fill_manual( values = c( "#D0873C", "#8FB339", "#7A5195", "#88CCEE" ),
                     labels = c( "Growth", "Recruitment", "Reproduction", "Survival" ),
                     guide = "legend" ) +
  scale_x_continuous( breaks = c( 2, 3 ), labels = model_labs ) +
  guides( fill = guide_legend( "Vital rate contribution" ) ) +
  labs( x = "Survival model", 
        y = "Uncertainty",
        title = "(c)",
        fill = "Vital rate contribution" ) +
  theme_bw( )

# Figure 2, all together

fig2ab <- fig2a + fig2b + plot_layout( ncol = 1, axes = "collect" )

fig2 <- wrap_plots( fig2ab ) + fig2c + plot_layout( ncol = 2 )


# Figure 3 ---------------------------------------------------------------------


