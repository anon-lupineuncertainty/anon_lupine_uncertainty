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

g_uncert_all <- read.csv( "data/uncert_g_all.csv" )
g_uncert_summary <- read.csv( "data/uncert_g_summary.csv" )

surv    <- read.csv( "data/surv.csv" )
grow    <- read.csv( "data/grow.csv" )
flow    <- read.csv( "data/flow.csv" )
fert    <- read.csv( "data/fert.csv" )

s_pars2 <- read.csv( "data/pars_sample2.csv" )
s_pars3 <- read.csv( "data/pars_sample3.csv" )

corr2_plot <- read.csv( "data/corr_plot2.csv" )
corr3_plot <- read.csv( "data/corr_plot3.csv" )

corr_plot6   <- read.csv( "data/corr_plot6.csv" )
corr_plot10  <- read.csv( "data/corr_plot10.csv" )
corr_plot20  <- read.csv( "data/corr_plot20.csv" )
corr_plot30  <- read.csv( "data/corr_plot30.csv" )
corr_plot40  <- read.csv( "data/corr_plot40.csv" )
corr_plot50  <- read.csv( "data/corr_plot50.csv" )
corr_plot75  <- read.csv( "data/corr_plot75.csv" )
corr_plot100 <- read.csv( "data/corr_plot100.csv" )
corr_plot200 <- read.csv( "data/corr_plot200.csv" )
corr_plot500 <- read.csv( "data/corr_plot500.csv" )


# Setup ------------------------------------------------------------------------

base_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    plot.margin = margin(5.5, 10, 5.5, 5.5)
  )


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
        fill = "Vital rate contribution",
        title = "(a)" ) +
  base_theme


# Figure 3b: comparison of sampled and mean lambda values, with recruitment
  # parameters varying and constant

lam_tab2 <- uncert2$lambdas
lam_tab2$mod <- "Varying"

lam_tab_ng <- uncert2_ng$lambdas
lam_tab_ng$mod <- "Constant"

lam_tab <- bind_rows( lam_tab2, lam_tab_ng )

lam_t <- t.test( lambda ~ mod, data = lam_tab )

mean_lam_plot <- mean_lam[c(1,1),]
colnames( mean_lam_plot ) <- c( "mod", "type", "value" )
mean_lam_plot[1,1] <- "Varying"
mean_lam_plot[2,1] <- "Constant"

fig3b <- ggplot() +
  geom_violin( data = lam_tab, aes( x = factor( mod ), y = lambda ),
               color = NA, alpha = 0.2, fill = "black" ) +
  geom_signif( data = lam_tab, aes( x = factor( mod ), y = lambda ),
               comparisons = list( c("Varying","Constant")),
               map_signif_level = TRUE ) +
  geom_segment( aes( x = 2.45, y = mean( lam_tab2$lambda ), 
                     xend = 1.55, yend = mean( lam_tab2$lambda ),
                     alpha = "Mean of samples" ),
                size = 0.9 ) +
  geom_segment( aes( x = 0.55, y = mean( lam_tab_ng$lambda ), 
                     xend = 1.45, yend = mean( lam_tab_ng$lambda ) ),
                size = 0.9 ) +
  geom_point( data = mean_lam_plot, aes( x = factor( mod ), y = value, 
                                             alpha = "Mean model" ),
              cex = 3, pch = 16 ) +
  scale_x_discrete( labels = c( "Constant", "Varying (resampled)" ) ) +
  ylim( 0.7, 1.9 ) +
  scale_alpha_manual( name = NULL,
                      values = c( 1, 1 ),
                      breaks = c( "Mean model", "Mean of samples" ),
                      guide = guide_legend( override.aes = list( linetype = c(0,1),
                                                                 shape = c(16,NA),
                                                                 color = "black" ) )
  ) +
  labs( x = "Recruitment parameters", 
        y = "Lambda",
        title = "(b)" ) +
  base_theme


# Patching them together

fig3 <- fig3a + fig3b + plot_layout( ncol = 1, axes = "collect_x" )



# Figure 4 ---------------------------------------------------------------------

# g_uncert_summary <- g_uncert_summary %>%
#   mutate( vital_rate = factor( vital_rate, levels = c( "survival",
#                                                        "growth",
#                                                        "reproduction",
#                                                        "recruitment",
#                                                        "total" ) ) )

fig4a <- ggplot( filter( g_uncert_summary, vital_rate == "total" ),
                 aes( x = sample_size, y = mean ) ) +
  geom_ribbon( aes( ymin = lower, ymax = upper ), alpha = 0.2 ) +
  geom_line( linewidth = 1 ) +
  geom_point( data = filter( g_uncert_summary, vital_rate == "total" ),
              aes( x = sample_size, y = mean ),
              size = 2.2, 
              color = "black" ) +
  scale_x_log10( breaks = c( 6, 10, 20, 30, 40, 50, 75, 100, 200, 500 ) ) +
  guides( alpha = "none" ) +
  labs( x = "Simulated germination trial sample size",
        y = "Variance in λ") +
  base_theme

fig4b <- ggplot( filter( g_uncert_summary, vital_rate != "total" ),
                 aes( x = sample_size, y = mean, color = vital_rate ) ) +
  geom_line( linewidth = 1 ) +
  scale_color_manual( values = c( "#8FB339", "#88CCEE", "#7A5195", "#D0873C" ),
                      labels = c( "Growth", "Recruitment",
                                  "Reproduction", "Survival" ),
                      guide = "legend",
                      name = "Vital rate" ) +
  scale_x_log10( breaks = c( 6, 10, 20, 30, 40, 50, 75, 100, 200, 500 ) ) +
  labs( x = "Simulated germination trial sample size",
        y = "Uncertainty contribution") +
  base_theme

fig4_test <- ggplot() +
  geom_ribbon( data = filter( g_uncert_summary, vital_rate == "total" ),
               aes( x = sample_size, ymin = lower, ymax = upper ),
               fill = "black",
               alpha = 0.2 ) +
  geom_area( data = filter( g_uncert_summary, vital_rate != "total" ),
             aes( x = sample_size, y = mean, fill = vital_rate ),
             alpha = 0.4, 
             color = NA ) +
  geom_line( data = filter( g_uncert_summary, vital_rate == "total" ),
             aes( x = sample_size, y = mean, linetype = "Model uncertainty" ),
             linewidth = 1.2,
             color = "black" ) +
  geom_point( data = filter( g_uncert_summary, vital_rate == "total" ),
              aes( x = sample_size, y = mean, shape = "Model uncertainty" ),
              size = 2.2, 
              color = "black" ) +
  scale_x_log10( breaks = c( 6, 10, 20, 30, 40, 50, 75, 100, 200, 500 ) ) +
  scale_fill_manual( values = c( growth      = "#8FB339",
                                 recruitment = "#88CCEE",
                                 reproduction = "#7A5195",
                                 survival    = "#D0873C" ),
                     labels = c( "Growth",
                                 "Recruitment", 
                                 "Reproduction",
                                 "Survival" ),
                     name = "Vital rate" ) +
  scale_linetype_manual( values = 1, guide = "none" ) +
  scale_shape_manual( values = 16,
                      guide = guide_legend( title = NULL, override.aes = list( color = "black",
                                                                               linewidth = 1.2,
                                                                               linetype = 1,
                                                                               shape = 16 ) ) ) +
  labs( x = "Simulated germination trial sample size",
        y = "Uncertainty" ) +
  base_theme


fig4_test2 <- ggplot( filter( g_uncert_summary, vital_rate == "total" ),
                      aes( x = sample_size, y = mean ) ) +
  geom_bar( data = filter( g_uncert_summary, vital_rate != "total" ), 
            aes( x = sample_size, y = mean, fill = vital_rate ), 
            stat = "identity", position = "stack",  width = 0.7 ) +
  geom_ribbon( aes( ymin = lower, ymax = upper ), alpha = 0.2 ) +
  geom_line( linewidth = 1 ) +
  geom_point( data = filter( g_uncert_summary, vital_rate == "total" ),
              aes( x = sample_size, y = mean ),
              size = 2.2, 
              color = "black" ) +
  scale_x_log10( breaks = c( 6, 10, 20, 30, 40, 50, 75, 100, 200, 500 ) ) +
  scale_fill_manual( values = c( growth      = "#8FB339",
                                 recruitment = "#88CCEE",
                                 reproduction = "#7A5195",
                                 survival    = "#D0873C" ),
                     labels = c( "Growth",
                                 "Recruitment", 
                                 "Reproduction",
                                 "Survival" ),
                     name = "Vital rate" ) +
  guides( alpha = "none" ) +
  labs( x = "Simulated germination trial sample size",
        y = "Uncertainty" ) +
  base_theme


g_bar <- g_uncert_summary %>%
  filter(vital_rate != "total") %>%
  arrange(sample_size, vital_rate) %>%
  group_by(sample_size) %>%
  mutate(
    ymax = cumsum(mean),
    ymin = lag(ymax, default = 0),
    xmin = sample_size * 0.98,
    xmax = sample_size * 1.02
  )

fig4_test3 <- ggplot( filter( g_uncert_summary, vital_rate == "total" ),
                      aes( x = sample_size, y = mean ) ) +
  geom_rect( data = g_bar, aes( xmin = xmin,
                                xmax = xmax,
                                ymin = ymin,
                                ymax = ymax,
                                fill = vital_rate ) ) +
  geom_ribbon( aes( ymin = lower, ymax = upper ), alpha = 0.2 ) +
  geom_line( linewidth = 1 ) +
  geom_point( data = filter( g_uncert_summary, vital_rate == "total" ),
              aes( x = sample_size, y = mean ),
              size = 2.2, 
              color = "black" ) +
  scale_x_log10( breaks = c( 6, 10, 20, 30, 40, 50, 75, 100, 200, 500 ) ) +
  scale_fill_manual( values = c( growth      = "#8FB339",
                                 recruitment = "#88CCEE",
                                 reproduction = "#7A5195",
                                 survival    = "#D0873C" ),
                     labels = c( "Growth",
                                 "Recruitment", 
                                 "Reproduction",
                                 "Survival" ),
                     name = "Vital rate" ) +
  guides( alpha = "none" ) +
  labs( x = "Simulated germination trial sample size",
        y = "Uncertainty" ) +
  base_theme

# Supplemental figures =========================================================

# Function to plot model predictions

plot_mod_pred <- function( data, mod ){
  
  plot_out <- vector( mode = "list", length = ( length( mod ) - 1 ) )
  
  if( grepl( "surv", deparse( substitute( data ) ), fixed = TRUE ) ) {
    fx <- plot_surv_mod
  } else {
    if( grepl( "grow", deparse( substitute( data ) ), fixed = TRUE ) ){
      fx <- plot_grow_mod
    } else {
      if( grepl( "flow", deparse( substitute( data ) ), fixed = TRUE) ){
        fx <- plot_flow_mod
      } else 
        fx <- plot_fert_mod 
    }
  }
  mod_names <- c( "(a) Intercept model", "(b) Linear model", "(c) Quadratic model", 
                  "(d) Cubic model" )
  
  for( i in 1:( length( mod ) - 1 ) ){
    p_temp <- fx( data, mod, mod_names[i], i )
    
    plot_out[[i]] <- p_temp
    
  }
  
  return( plot_out )
}


# Figure S2 --------------------------------------------------------------------

# Fit of growth models

plot_grow_mod <- function( df, mod, name, index ){
  
  seq_x <- seq( min( df$log_area_t0, na.rm = T ),
                max( df$log_area_t0, na.rm = T ),
                length.out = 100 )
  
  grow_pred <- predictor_fun( seq_x, coef( mod[[index]] ) ) %>%
    data.frame( log_area_t0 = seq_x,
                log_area_t1 = . )
  
  p_grow <- ggplot( df, aes( x = log_area_t0, y = log_area_t1 ) ) +
    geom_point( alpha = 0.1 ) +
    geom_line( data = grow_pred,
               aes( x = log_area_t0,
                    y = log_area_t1 ),
               color = "#8FB339", lwd = 2 ) +
    labs( title = name,
          x = expression(log(Size~at~time~italic(t[0]))),
          y = expression(log(Size~at~time~italic(t[0]))) ) +
    base_theme
  
  return( p_grow )
}


p_grow_pred <- plot_mod_pred( grow, m_grow )

figs2 <- p_grow_pred[[1]] + p_grow_pred[[2]] + p_grow_pred[[3]] + plot_layout( axes = "collect" )


# Figure S3 --------------------------------------------------------------------

# Fit of survival models

plot_surv_mod <- function( df, mod, name, index ){
  
  seq_x <- seq( min( df$log_area_t0, na.rm = T ),
                max( df$log_area_t0, na.rm = T ),
                length.out = 100 )
  
  surv_pred <- predictor_fun( seq_x, coef( mod[[index]] ) ) %>%
    boot::inv.logit() %>%
    data.frame( log_area_t0 = seq_x,
                surv_t1 = . )
  
  p_surv <- ggplot() + 
    geom_jitter( data = df, aes( x = log_area_t0, y = surv_t1 ),
                 width = 0.05, height = 0.1, alpha = 0.1 ) +
    geom_errorbar( data = plot_binned_prop( df, 10, log_area_t0, surv_t1 ), 
                   aes( x = log_area_t0, ymin = lwr, ymax = upr ),
                   size = 0.5, width = 0.2 ) +
    geom_point( data = plot_binned_prop( df, 10, log_area_t0, surv_t1 ),
                aes( x = log_area_t0, y = surv_t1 ), color = "#D0873C", cex = 2 ) +
    geom_line( data = surv_pred, 
               aes( x = log_area_t0,
                    y = surv_t1 ),
               color = "#D55E00", lwd = 2 ) +
    ylim( -0.1, 1.1 ) +
    labs( x = expression(log(Size~at~time~italic(t[0]))),
          y = expression(Probability~of~survival~to~time~italic(t[1])),
          title = name ) +
    base_theme
  
  return( p_surv )
}

p_surv_pred <- plot_mod_pred( surv, m_surv )

figs3 <- p_surv_pred[[1]] + p_surv_pred[[2]] + 
  p_surv_pred[[3]] + p_surv_pred[[4]] + plot_layout( ncol = 2, axes = "collect" )


# Figure S4 --------------------------------------------------------------------

# Fit of flowering models

plot_flow_mod <- function( df, mod, name, index ){
  
  seq_x <- seq( min( df$log_area_t0, na.rm = T ),
                max( df$log_area_t0, na.rm = T ),
                length.out = 100 )
  
  flow_pred <- predictor_fun( seq_x, coef( mod[[index]] ) ) %>%
    boot::inv.logit() %>%
    data.frame( log_area_t0 = seq_x,
                flow_t0 = . )
  
  p_flow <- ggplot( data = df, aes( x = log_area_t0, y = flow_t0 ) ) +
    geom_jitter( width = 0, height = 0.1, alpha = 0.1 ) + 
    geom_point( data = plot_binned_prop( df, 10, log_area_t0, flow_t0 ),
                aes( x = log_area_t0, y = flow_t0 ), 
                color = "#7A5195",
                cex = 2 ) +
    geom_errorbar( data = plot_binned_prop( df, 10, log_area_t0, flow_t0 ),
                   aes( x = log_area_t0, ymin = lwr, ymax = upr ),
                   size = 0.5, width = 0.5 ) +
    geom_line( data = flow_pred, 
               aes( x = log_area_t0,
                    y = flow_t0 ),
               color = "#7A5195", lwd = 2 ) +
    scale_y_continuous( breaks = c( 0.1, 0.5, 0.9 ) ) +
    ylim( -0.1, 1.1 ) +
    labs( x = expression(log(Size~at~time~italic(t[0]))),
          y = expression(Probability~of~flowering~at~time~italic(t[0])),
          title = name ) +
    base_theme 

  return( p_flow )
}

p_flow_pred <- plot_mod_pred( flow, m_flow )

figs4 <- p_flow_pred[[1]] + p_flow_pred[[2]] + p_flow_pred[[3]] + plot_layout( axes = "collect" )


# Figure S5 --------------------------------------------------------------------

# Fit of fertility models

plot_fert_mod <- function( df, mod, name, index ){
  
  seq_x <- seq( min( df$log_area_t0, na.rm = T ),
                max( df$log_area_t0, na.rm = T ),
                length.out = 100 )
  
  fert_pred <- predictor_fun( seq_x, coef( mod[[index]] ) ) %>%
    exp() %>%
    data.frame( log_area_t0 = seq_x,
                numrac_t0 = . )
  
  p_fert <- ggplot( data = df, aes( x = log_area_t0, y = numrac_t0 ) ) +
    geom_point( alpha = 0.1 ) + 
    geom_line( data = fert_pred,
               aes( x = log_area_t0,
                    y = numrac_t0 ),
               color = "#7A5195", lwd = 2 ) +
    labs( title = name,
          x = expression(log(Size~at~time~italic(t[0]))),
          y = expression(Number~of~racemes~produced~at~time~italic(t[0]))) +
    base_theme
  
  return( p_fert )
}

p_fert_pred <- plot_mod_pred( fert, m_fert )

figs5 <- p_fert_pred[[1]] + p_fert_pred[[2]] + p_fert_pred[[3]] + plot_layout( axes = "collect" )


# Figure S7 --------------------------------------------------------------------

figs7 <- ggplot( corr2_plot, aes( x = Var1, y = Var2, fill = correlation ) ) + 
  geom_tile() +
  geom_text( aes( label = text ) ) +
  base_theme +
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


# Figure S8 --------------------------------------------------------------------

corr_plot_fxn <- function( df, caption ){
  
  temp_plot <- ggplot( df, aes( x = Var1, y = Var2, fill = mean_corr ) ) + 
    geom_tile() +
    base_theme +
    labs( title = caption ) + 
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
  
  return( temp_plot )
}


# Correlation plots for the simulated germination datasets

figs8a <- corr_plot_fxn( corr_plot6,
                         caption = "(a) 6 samples" ) +
  theme( legend.position = "none" )

figs8b <- corr_plot_fxn( corr_plot10,
                         caption = "(b) 10 samples" ) +
  theme( legend.position = "none" )

figs8c <- corr_plot_fxn( corr_plot20,
                         caption = "(c) 20 samples" ) +
  theme( legend.position = "none" )

figs8d <- corr_plot_fxn( corr_plot30,
                         caption = "(d) 30 samples" ) +
  theme( legend.position = "none" )

figs8e <- corr_plot_fxn( corr_plot40,
                         caption = "(e) 40 samples" ) +
  theme( legend.position = "none" )

figs8f <- corr_plot_fxn( corr_plot50,
                         caption = "(f) 50 samples" ) +
  theme( legend.position = "none" )

figs8g <- corr_plot_fxn( corr_plot75,
                         caption = "(g) 75 samples" ) +
  theme( legend.position = "none" )

figs8h <- corr_plot_fxn( corr_plot100, 
                         caption = "(h) 100 samples" ) +
  theme( legend.position = "none" )

figs8i <- corr_plot_fxn( corr_plot200,
                         caption = "(i) 200 samples" ) +
  theme( legend.position = "none" )

figs8j <- corr_plot_fxn( corr_plot500,
                         caption = "(j) 500 samples" ) +
  theme( legend.position = "none" )

figs8 <- figs8a + figs8b + figs8c + figs8d + figs8e + figs8f + figs8g + figs8h + 
  figs8i + figs8j + plot_layout( ncol = 4 )


# Figure S9 ---------------------------------------------------------------------

model_labs <- c( "Quadratic", "Cubic" )

lam_tab2 <- uncert2$lambdas
lam_tab2$type <- 2

lam_tab3 <- uncert3$lambdas
lam_tab3$type <- 3

lam_tab <- bind_rows( lam_tab2, lam_tab3 )

lam_t <- t.test( lambda ~ type, data = lam_tab )

mean_lam_plot2 <- mean_lam[c(1:4),]

# Figure S9a

figs9a <- ggplot() +
  geom_violin( data = lam_tab, aes( x = factor( type ), y = lambda,
                                    fill = factor( type ) ),
               color = NA, alpha = 0.4 ) +
  geom_signif( data = lam_tab, aes( x = factor( type ), y = lambda ),
               comparisons = list( c("2","3")),
               map_signif_level = TRUE ) +
  scale_fill_manual( values = c( "#f5a351", "#a85603" ), guide = "none" ) +
  geom_segment( aes( x = 0.55, y = mean_lam_plot2[2,3], 
                     xend = 1.45, yend = mean_lam_plot2[2,3],
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
  base_theme

# Figure S9b

var_tab2 <- uncert2$vr_uncert
var_tab2[5,] <- list( "total", uncert2$mod_uncert )
var_tab2$type <- 2

var_tab3 <- uncert3$vr_uncert
var_tab3[5,] <- list( "total", uncert3$mod_uncert )
var_tab3$type <- 3

var_tab <- bind_rows( var_tab2, var_tab3 )


figs9b <- ggplot( ) + 
  geom_bar( data = var_tab[which(var_tab$vital_rate != "total" ),], 
            aes( x = type, y = variance_sum, fill = vital_rate ), 
            stat = "identity", position = "stack", width = 0.7 ) +
  geom_bar( data = var_tab[which(var_tab$vital_rate == "total" ),], 
            aes( x = type, y = variance_sum, alpha = "Model uncertainty" ), 
            stat = "identity", position = "stack", color = "black", fill = NA,
            linewidth = 0.8, width = 0.7 ) +
  scale_fill_manual( values = c( "#8FB339", "#88CCEE", "#D0873C" ),
                     labels = c( "Growth", "Recruitment", "Survival" ),
                     guide = "legend" ) +
  # scale_fill_manual( values = c( "#8FB339", "#88CCEE", "#7A5195", "#D0873C" ),
  #                    labels = c( "Growth", "Recruitment", "Reproduction", "Survival" ),
  #                    guide = "legend" ) +
  scale_x_continuous( breaks = c( 2, 3 ), labels = model_labs ) +
  scale_y_sqrt() +
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
  base_theme

# Figure S9c

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

figs9_base <- ggplot() + 
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

# Figure S9c

figs9c <- figs9_base +
  geom_text( data = model_AIC[which( model_AIC$model == 2 ),], 
             aes( x = size, y = survival, label = AIC ) ) +
  geom_line( data = surv_pred2, 
             aes( x = log_area_t0,
                  y = surv_t1 ),
             color = "#f5a351", lwd = 1.5, alpha = 0.8 ) +
  labs( title = "(c) Quadratic survival model" )

# Figure S9d

figs9d <- figs9_base +
  geom_text( data = model_AIC[which( model_AIC$model == 3 ),], 
             aes( x = size, y = survival, label = AIC ) ) +
  geom_line( data = surv_pred3, 
             aes( x = log_area_t0,
                  y = surv_t1 ),
             color = "#a85603", lwd = 1.5, alpha = 0.8 ) +
  labs( title = "(d) Cubic survival model" )


# Figure S9, all together

figs9ab <- figs9a + figs9b + plot_layout( ncol = 1, axes = "collect_x" )

figs9cd <- figs9c + figs9d + plot_layout( ncol = 1, axes = "collect" )

figs9 <- wrap_plots( figs9ab ) + wrap_plots( figs9cd ) + 
  plot_layout( ncol = 2, widths = c(1,2) ) & base_theme


# Figure S10 --------------------------------------------------------------------

figs10a <- ggplot( corr2_plot, aes( x = Var1, y = Var2, fill = correlation ) ) + 
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


figs10b <- ggplot( corr3_plot, aes( x = Var1, y = Var2, fill = correlation ) ) + 
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

figs10b <- figs10b +
  theme(
    legend.position = "none"
  )

figs10 <- figs10a / figs10b


# Save output ------------------------------------------------------------------

ggsave( "results/Figure1.pdf", fig1, width = 180, height = 100, units = "mm",
        device = cairo_pdf )
ggsave( "results/Figure2.pdf", fig2, width = 85, height = 100, units = "mm",
        device = cairo_pdf )
ggsave( "results/Figure3.pdf", fig3, width = 130, height = 150, units = "mm",
        device = cairo_pdf )


ggsave( "results/FigureS2.pdf", figs2, width = 220, height = 110, units = "mm",
        device = cairo_pdf )
ggsave( "results/FigureS3.pdf", figs3, width = 200, height = 160, units = "mm",
        device = cairo_pdf )
ggsave( "results/FigureS4.pdf", figs4, width = 220, height = 110, units = "mm",
        device = cairo_pdf )
ggsave( "results/FigureS5.pdf", figs5, width = 220, height = 110, units = "mm",
        device = cairo_pdf )


ggsave( "results/FigureS7.pdf", figs7, width = 180, height = 110, units = "mm",
        device = cairo_pdf )
ggsave( "results/FigureS8.pdf", figs8, width = 280, height = 220, units = "mm",
        device = cairo_pdf )
ggsave( "results/FigureS9.pdf", figs9, width = 180, height = 110, units = "mm",
        device = cairo_pdf )
ggsave( "results/FigureS10.pdf", figs10, width = 180, height = 220, units = "mm",
        device = cairo_pdf )

