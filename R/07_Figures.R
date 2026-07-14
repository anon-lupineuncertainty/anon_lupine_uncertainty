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
  # - data/mean_lambdas.csv: Lambda values from the mean models
  # - data/pars_sample2.csv: Sampled parameter values with quadratic survival model
  # - data/pars_sample3.csv: Sampled parameter values with cubic survival model
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
library( ggsignif )


# Data -------------------------------------------------------------------------

uncert2    <- readRDS( "data/uncert2.rds" )
uncert2_ng <- readRDS( "data/uncert2_ng.rds" )
uncert3    <- readRDS( "data/uncert3.rds" )

uncert_comp_plot <- read.csv( "data/uncert_comp_plot.csv" )

mean_lam <- read.csv( "data/mean_lambdas_all.csv" )

surv    <- read.csv( "data/surv.csv" )

s_pars2 <- read.csv( "data/pars_sample2.csv" )
s_pars3 <- read.csv( "data/pars_sample3.csv" )

corr2_plot <- read.csv( "data/corr_plot2.csv" )
corr3_plot <- read.csv( "data/corr_plot3.csv" )


# Figure 1 ---------------------------------------------------------------------

fig1 <- plot( uncert2, type = "param" )


# Figure 2 ---------------------------------------------------------------------

fig2 <- plot( uncert2, type = "vr" )


# Figure 3 --------------------------------------------------------------------

# Figure 3a: comparison of explained uncertainty with recruitment parameters
  # varying and constant

fig3a <- ggplot( ) + 
  geom_bar( data = uncert_comp_plot[which(uncert_comp_plot$vital_rate != "total" ),], 
            aes( x = mod, y = variance_sum, fill = vital_rate ), 
            stat = "identity", position = "stack", width = 0.7 ) +
  geom_bar( data = uncert_comp_plot[which(uncert_comp_plot$vital_rate == "total" ),], 
            aes( x = mod, y = variance_sum, alpha = "Model uncertainty" ), 
            stat = "identity", position = "stack", color = "black", fill = NA,
            linewidth = 0.8, width = 0.7 ) +
  scale_y_sqrt() +
  scale_fill_manual( values = c( "#8FB339", "#88CCEE", "#7A5195", "#D0873C" ),
                     labels = c( "Growth", "Recruitment", "Reproduction", "Survival" ),
                     guide = "legend" ) +
  scale_alpha_manual( name = NULL,
                      values = 1,
                      breaks = "Model uncertainty",
                      guide = guide_legend( override.aes = list( color = "black" ) )
  ) +
  scale_x_discrete( labels = c( "Constant", "Varying (resampled)" ) ) +
  guides( fill = guide_legend( "Vital rate contribution" ) ) +
  labs( x = "Recruitment parameters", 
        y = "Uncertainty",
        fill = "Vital rate contribution" ) +
  theme_bw( )


# Figure 3b: comparison of sampled and mean lambda values, with recruitment
  # parameters varying and constant

lam_tab2 <- uncert2$lambdas
lam_tab2$mod <- "With varying recruitment"

lam_tab_ng <- uncert2_ng$lambdas
lam_tab_ng$mod <- "With constant recruitment"

lam_tab <- bind_rows( lam_tab2, lam_tab_ng )

lam_t <- t.test( lambda ~ type, data = lam_tab )

fig3b <- ggplot() +
  geom_violin( data = lam_tab, aes( x = factor( mod ), y = lambda,
                                    fill = factor( mod ) ),
               color = NA, alpha = 0.4 ) +
  geom_signif( data = lam_tab, aes( x = factor( mod ), y = lambda ),
               comparisons = list( c("With varying recruitment","With constant recruitment")),
               map_signif_level = TRUE ) +
  scale_fill_manual( values = c( "#f5a351", "#a85603" ), guide = "none" ) +
  geom_segment( aes( x = 0.55, y = mean( lam_tab2$lambda ), 
                     xend = 1.45, yend = mean( lam_tab2$lambda ),
                     alpha = "Mean of samples" ),
                size = 0.9, color = "#f5a351" ) +
  geom_segment( aes( x = 1.55, y = mean( lam_tab_ng$lambda ), 
                     xend = 2.45, yend = mean( lam_tab_ng$lambda ) ),
                size = 0.9, color = "#a85603"  ) +
  geom_point( data = mean_lam[c(1,3),], aes( x = factor( type ), y = value,
                                             color = factor( type ), 
                                             alpha = "Mean model" ),
              cex = 2, pch = 16 ) +
  scale_color_manual( values = c( "#f5a351", "#a85603" ), guide = "none" ) +
  scale_x_discrete( labels = model_labs ) +
  scale_alpha_manual( name = NULL,
                      values = c( 1, 1 ),
                      breaks = c( "Mean model", "Mean of samples" ),
                      guide = guide_legend( override.aes = list( linetype = c(0,1),
                                                                 shape = c(16,NA),
                                                                 color = "#D0873C" ) )
  ) +
  ylim( min( lam_tab$lambda ) - 0.03, max( lam_tab$lambda ) + 0.3 ) +
  labs( x = "Survival model", 
        y = "Lambda",
        title = "(a)" ) +
  theme_bw( )


# Figure 3 ---------------------------------------------------------------------

model_labs <- c( "Quadratic", "Cubic" )

lam_tab2 <- uncert2$lambdas
lam_tab2$type <- 2

lam_tab3 <- uncert3$lambdas
lam_tab3$type <- 3

lam_tab <- bind_rows( lam_tab2, lam_tab3 )

lam_t <- t.test( lambda ~ type, data = lam_tab )

# need to update this to rely on mean_lambdas instead of mean_lam

# Figure 3a

fig3a <- ggplot() +
  geom_violin( data = lam_tab, aes( x = factor( type ), y = lambda,
                                    fill = factor( type ) ),
               color = NA, alpha = 0.4 ) +
  geom_signif( data = lam_tab, aes( x = factor( type ), y = lambda ),
               comparisons = list( c("2","3")),
               map_signif_level = TRUE ) +
  scale_fill_manual( values = c( "#f5a351", "#a85603" ), guide = "none" ) +
  geom_segment( aes( x = 0.55, y = mean_lam[2,3], 
                     xend = 1.45, yend = mean_lam[2,3],
                     alpha = "Mean of samples" ),
                size = 0.9, color = "#f5a351" ) +
  geom_segment( aes( x = 1.55, y = mean_lam[4,3], 
                     xend = 2.45, yend = mean_lam[4,3] ),
                size = 0.9, color = "#a85603"  ) +
  geom_point( data = mean_lam[c(1,3),], aes( x = factor( type ), y = value,
                                             color = factor( type ), 
                                             alpha = "Mean model" ),
              cex = 2, pch = 16 ) +
  scale_color_manual( values = c( "#f5a351", "#a85603" ), guide = "none" ) +
  scale_x_discrete( labels = model_labs ) +
  scale_alpha_manual( name = NULL,
                      values = c( 1, 1 ),
                      breaks = c( "Mean model", "Mean of samples" ),
                      guide = guide_legend( override.aes = list( linetype = c(0,1),
                                                                shape = c(16,NA),
                                                                color = "#D0873C" ) )
                      ) +
  ylim( min( lam_tab$lambda ) - 0.03, max( lam_tab$lambda ) + 0.3 ) +
  labs( x = "Survival model", 
        y = "Lambda",
        title = "(a)" ) +
  theme_bw( )

# Figure 3b

var_tab2 <- uncert2$vr_uncert
var_tab2[5,] <- list( "total", uncert2$mod_uncert )
var_tab2$type <- 2

var_tab3 <- uncert3$vr_uncert
var_tab3[5,] <- list( "total", uncert3$mod_uncert )
var_tab3$type <- 3

var_tab <- bind_rows( var_tab2, var_tab3 )


fig3b <- ggplot( ) + 
  geom_bar( data = var_tab[which(var_tab$vital_rate != "total" ),], 
            aes( x = type, y = variance_sum, fill = vital_rate ), 
            stat = "identity", position = "stack", width = 0.7 ) +
  geom_bar( data = var_tab[which(var_tab$vital_rate == "total" ),], 
            aes( x = type, y = variance_sum, alpha = "Model uncertainty" ), 
            stat = "identity", position = "stack", color = "black", fill = NA,
            linewidth = 0.8, width = 0.7 ) +
  scale_fill_manual( values = c( "#8FB339", "#88CCEE", "#7A5195", "#D0873C" ),
                     labels = c( "Growth", "Recruitment", "Reproduction", "Survival" ),
                     guide = "legend" ) +
  scale_x_continuous( breaks = c( 2, 3 ), labels = model_labs ) +
  scale_alpha_manual( name = NULL,
                      values = 1,
                      breaks = "Model uncertainty",
                      guide = guide_legend( override.aes = list( color = "black" ) )
  ) +
  guides( fill = guide_legend( "Vital rate contribution" ) ) +
  labs( x = "Survival model", 
        y = "Uncertainty",
        title = "(b)",
        fill = "Vital rate contribution" ) +
  theme_bw( )

# Figure 3c

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

fig3_base <- ggplot() + 
  geom_jitter( data = surv, aes( x = log_area_t0, y = surv_t1 ),
               width = 0.05, height = 0.1, alpha = 0.1 ) +
  geom_errorbar( data = surv_binned, 
                 aes( x = log_area_t0, ymin = lwr, ymax = upr ),
                 size = 0.5, width = 0.2 ) +
  geom_point( data = surv_binned,
              aes( x = log_area_t0, y = surv_t1 ), color = "#D0873C", cex = 2 ) +
  ylim( -0.1, 1.1 ) +
  labs( x = expression(log(Size~at~time~italic(t[0]))),
        y = expression(Probability~of~survival~to~time~italic(t[1])) ) +
  theme_bw( )

surv_mod2 <- glm( surv_t1 ~ log_area_t0 + log_area_t02, 
                   data = surv, family = "binomial" )
surv_mod3 <- glm( surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03, 
                   data = surv, family = "binomial" )

model_AIC  <- data.frame( model = c( 2, 3 ), size = c( 5, 5 ), 
                          survival = c( 0.25, 0.25 ), 
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

# Figure 3c

fig3c <- fig3_base +
  geom_text( data = model_AIC[which( model_AIC$model == 2 ),], 
             aes( x = size, y = survival, label = AIC ) ) +
  geom_line( data = surv_pred2, 
             aes( x = log_area_t0,
                  y = surv_t1 ),
             color = "#f5a351", lwd = 1.5, alpha = 0.8 ) +
  labs( title = "(c) Quadratic model" )

# Figure 3d

fig3d <- fig3_base +
  geom_text( data = model_AIC[which( model_AIC$model == 3 ),], 
             aes( x = size, y = survival, label = AIC ) ) +
  geom_line( data = surv_pred3, 
             aes( x = log_area_t0,
                  y = surv_t1 ),
             color = "#a85603", lwd = 1.5, alpha = 0.8 ) +
  labs( title = "(d) Cubic model" )


base_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    plot.margin = margin(5.5, 10, 5.5, 5.5)
  )

# Figure 3, all together

fig3ab <- fig3a + fig3b + plot_layout( ncol = 1, axes = "collect_x" )

fig3cd <- fig3c + fig3d + plot_layout( ncol = 1, axes = "collect" )

fig3 <- wrap_plots( fig3ab ) + wrap_plots( fig3cd ) + 
  plot_layout( ncol = 2, widths = c(1,2) ) & base_theme


# Figure 3 ---------------------------------------------------------------------

# Revised version of Figure 3, visualizing the uncertainty analyses holding the
  # germination coefficients constant
# Clean up or remove later

# Figure 3a

lam_tab2_ng <- uncert2_ng$lambdas
lam_tab2_ng$type <- 2

lam_tab3_ng <- uncert3_ng$lambdas
lam_tab3_ng$type <- 3

lam_tab_ng <- bind_rows( lam_tab2_ng, lam_tab3_ng )

lam_t_ng <- t.test( lambda ~ type, data = lam_tab_ng )

mean_lam_ng <- data.frame( type  = c( 2, 2, 3, 3 ),
                        model = c( "Mean model", "Mean of sampled models",
                                   "Mean model", "Mean of sampled models" ),
                        value = c( mean_lambdas[1,2], mean( lam_tab2_ng$lambda ),
                                   mean_lambdas[2,2], mean( lam_tab3_ng$lambda ) ) )

fig3a_ng <- ggplot() +
  geom_violin( data = lam_tab_ng, aes( x = factor( type ), y = lambda,
                                    fill = factor( type ) ),
               color = NA, alpha = 0.4 ) +
  geom_signif( data = lam_tab_ng, aes( x = factor( type ), y = lambda ),
               comparisons = list( c("2","3")),
               map_signif_level = TRUE ) +
  scale_fill_manual( values = c( "#f5a351", "#a85603" ), guide = "none" ) +
  geom_segment( aes( x = 0.55, y = mean_lam_ng[2,3], 
                     xend = 1.45, yend = mean_lam_ng[2,3],
                     alpha = "Mean of samples" ),
                size = 0.9, color = "#f5a351" ) +
  geom_segment( aes( x = 1.55, y = mean_lam_ng[4,3], 
                     xend = 2.45, yend = mean_lam_ng[4,3] ),
                size = 0.9, color = "#a85603"  ) +
  geom_point( data = mean_lam_ng[c(1,3),], aes( x = factor( type ), y = value,
                                             color = factor( type ), 
                                             alpha = "Mean model" ),
              cex = 2, pch = 16 ) +
  scale_color_manual( values = c( "#f5a351", "#a85603" ), guide = "none" ) +
  scale_x_discrete( labels = model_labs ) +
  scale_alpha_manual( name = NULL,
                      values = c( 1, 1 ),
                      breaks = c( "Mean model", "Mean of samples" ),
                      guide = guide_legend( override.aes = list( linetype = c(0,1),
                                                                 shape = c(16,NA),
                                                                 color = "#D0873C" ) )
  ) +
  ylim( min( lam_tab_ng$lambda ) - 0.03, max( mean_lam_ng$value ) + 0.3 ) +
  labs( x = "Survival model", 
        y = "Lambda",
        title = "(a)" ) +
  theme_bw( )

# Figure 3b

var_tab2_ng <- uncert2_ng$vr_uncert
var_tab2_ng[4,] <- list( "total", uncert2_ng$mod_uncert )
var_tab2_ng$type <- 2

var_tab3_ng <- uncert3_ng$vr_uncert
var_tab3_ng[4,] <- list( "total", uncert3_ng$mod_uncert )
var_tab3_ng$type <- 3

var_tab_ng <- bind_rows( var_tab2_ng, var_tab3_ng )


fig3b_ng <- ggplot( ) + 
  geom_bar( data = var_tab_ng[which(var_tab_ng$vital_rate != "total" ),], 
            aes( x = type, y = variance_sum, fill = vital_rate ), 
            stat = "identity", position = "stack", width = 0.7 ) +
  geom_bar( data = var_tab_ng[which(var_tab_ng$vital_rate == "total" ),], 
            aes( x = type, y = variance_sum, alpha = "Model uncertainty" ), 
            stat = "identity", position = "stack", color = "black", fill = NA,
            linewidth = 0.8, width = 0.7 ) +
  scale_fill_manual( values = c( "#8FB339", "#7A5195", "#D0873C" ),
                     labels = c( "Growth", "Reproduction", "Survival" ),
                     guide = "legend" ) +
  scale_x_continuous( breaks = c( 2, 3 ), labels = model_labs ) +
  scale_alpha_manual( name = NULL,
                      values = 1,
                      breaks = "Model uncertainty",
                      guide = guide_legend( override.aes = list( color = "black" ) )
  ) +
  guides( fill = guide_legend( "Vital rate contribution" ) ) +
  labs( x = "Survival model", 
        y = "Uncertainty",
        title = "(b)",
        fill = "Vital rate contribution" ) +
  theme_bw( )

# Figure 3, all together

fig3ab_ng <- fig3a_ng + fig3b_ng + plot_layout( ncol = 1, axes = "collect_x" )

fig3cd <- fig3c + fig3d + plot_layout( ncol = 1, axes = "collect" )

fig3_ng <- wrap_plots( fig3ab_ng ) + wrap_plots( fig3cd ) + 
  plot_layout( ncol = 2, widths = c(1,2) ) & base_theme





# Figure 4 ---------------------------------------------------------------------




# Figure S1 --------------------------------------------------------------------

figs1a <- ggplot( corr2_plot, aes( x = Var1, y = Var2, fill = correlation ) ) + 
  geom_tile() +
  geom_text( aes( label = text ) ) +
  base_theme +
  labs( title = "(a) With quadratic survival model" ) + 
  scale_fill_gradient2( low = "#D55E00", mid = "white", high = "#0072B2",
                        midpoint = 0, limits = c(-1,1), name = "Correlation" ) +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1 ) ) +
  annotate( "rect", xmin = 0.5, xmax = 6.5, ymin = 16.7, ymax = 17.0,
            fill = "#7A5195" ) +
  annotate( "rect", xmin = 6.5, xmax = 10.5, ymin = 16.7, ymax = 17.0,
            fill = "#88CCEE" ) +
  annotate( "rect", xmin = 10.5, xmax = 13.5, ymin = 16.7, ymax = 17.0,
            fill = "#8FB339" ) +
  annotate( "rect", xmin = 13.5, xmax = 16.5, ymin = 16.7, ymax = 17.0,
            fill = "#D0873C" ) +
  annotate( "rect", xmin = 16.7, xmax = 17.0, ymin = 0.5, ymax = 6.5,
            fill = "#7A5195" ) +
  annotate( "rect", xmin = 16.7, xmax = 17.0, ymin = 6.5, ymax = 10.5,
            fill = "#88CCEE" ) +
  annotate( "rect", xmin = 16.7, xmax = 17.0, ymin = 10.5, ymax = 13.5,
            fill = "#8FB339" ) +
  annotate( "rect", xmin = 16.7, xmax = 17.0, ymin = 13.5, ymax = 16.5, 
            fill = "#D0873C" ) +
  geom_vline( xintercept = c( 6.5, 10.5, 13.5 ), color = "white", linewidth = 0.7 ) +
  geom_hline( yintercept = c( 6.5, 10.5, 13.5 ), color = "white", linewidth = 0.7 )


figs1b <- ggplot( corr3_plot, aes( x = Var1, y = Var2, fill = correlation ) ) + 
  geom_tile() +
  geom_text( aes( label = text ) ) +
  labs( title = "(b) With cubic survival model" ) + 
  base_theme + 
  scale_fill_gradient2( low = "#D55E00", mid = "white", high = "#0072B2",
                        midpoint = 0, limits = c(-1,1), name = "Correlation" ) +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1 ) ) +
  annotate( "rect", xmin = 0.5, xmax = 6.5, ymin = 17.7, ymax = 18.0,
            fill = "#7A5195" ) +
  annotate( "rect", xmin = 6.5, xmax = 10.5, ymin = 17.7, ymax = 18.0,
            fill = "#88CCEE" ) +
  annotate( "rect", xmin = 10.5, xmax = 13.5, ymin = 17.7, ymax = 18.0,
            fill = "#8FB339" ) +
  annotate( "rect", xmin = 13.5, xmax = 17.5, ymin = 17.7, ymax = 18.0,
            fill = "#D0873C" ) +
  annotate( "rect", xmin = 17.7, xmax = 18.0, ymin = 0.5, ymax = 6.5,
            fill = "#7A5195" ) +
  annotate( "rect", xmin = 17.7, xmax = 18.0, ymin = 6.5, ymax = 10.5,
            fill = "#88CCEE" ) +
  annotate( "rect", xmin = 17.7, xmax = 18.0, ymin = 10.5, ymax = 13.5,
            fill = "#8FB339" ) +
  annotate( "rect", xmin = 17.7, xmax = 18.0, ymin = 13.5, ymax = 17.5, 
            fill = "#D0873C" ) +
  geom_vline( xintercept = c( 6.5, 10.5, 13.5 ), color = "white", linewidth = 0.7 ) +
  geom_hline( yintercept = c( 6.5, 10.5, 13.5 ), color = "white", linewidth = 0.7 )

figs1b <- figs1b +
  theme(
    legend.position = "none"
  )

figs1 <- figs1a / figs1b


# Save output ------------------------------------------------------------------

ggsave( "results/Figure1.pdf", fig1, width = 180, height = 100, units = "mm",
        device = cairo_pdf )
ggsave( "results/Figure2.pdf", fig2, width = 85, height = 100, units = "mm",
        device = cairo_pdf )
ggsave( "results/Figure3.pdf", fig3, width = 180, height = 130, units = "mm",
        device = cairo_pdf )


ggsave( "results/FigureS1.pdf", figs1, width = 180, height = 220, units = "mm",
        device = cairo_pdf )
