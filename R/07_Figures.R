# ==============================================================================
# Script: 07_Figures.R
#
# Prepare figures
#
# Purpose:
# This script prepares figures for inclusion in the manuscript.
#
# Inputs:
  # See previous scripts for descriptions of input files needed for plotting
#
# Outputs:
  # - results/Figure#.pdf: PDF of Figures 1-4 for the main text
  # - results/FigureS#.pdf: PDF of Figures S2-S11 for the appendices
# 
# Notes:
# This script primarily requires access to the objects produced by the previous
  # scripts and is therefore mostly reproducible from this repository and its
  # associated data archive. Figures S2-S5 and Figure S10e,f are not
  # reproducible, as they include plots of the raw demographic data.
# ==============================================================================

options( stringsAsFactors = F )
library( ipmr )
library( ggplot2 )
library( patchwork )
library( ggsignif )
library( ggh4x )


# Data -------------------------------------------------------------------------

uncert2 <- readRDS( "data/uncert2.rds" )

uncert_comp_plot <- read.csv( "data/uncert_comp_plot.csv" )
lam_tab_comp     <- read.csv( "data/sampled_lambdas_comp.csv" )

mean_lam   <- read.csv( "data/mean_lambdas_all.csv" )
lam_tab_23 <- read.csv( "data/sampled_lambdas_mf.csv" )
var_tab    <- read.csv( "data/var_cont.csv" )

mean_lam_ng <- read.csv( "data/mean_lambdas_ng.csv" )
lam_tab_ng  <- read.csv( "data/sampled_lambdas_ng.csv" )
var_tab_ng  <- read.csv( "data/var_cont_ng.csv" )

g_uncert_all <- read.csv( "data/uncert_g_all.csv" )
g_uncert_summary <- read.csv( "data/uncert_g_summary.csv" )

surv    <- read.csv( "data/surv.csv" )
grow    <- read.csv( "data/grow.csv" )

s_pars2 <- read.csv( "data/pars_sample2.csv" )
s_pars3 <- read.csv( "data/pars_sample3.csv" )

pars_mean2 <- read.csv( "data/pars_mean2.csv" )

germ6   <- read.csv( "data/pars_gsim6.csv" )
germ10  <- read.csv( "data/pars_gsim10.csv" )
germ20  <- read.csv( "data/pars_gsim20.csv" )
germ30  <- read.csv( "data/pars_gsim30.csv" )
germ40  <- read.csv( "data/pars_gsim40.csv" )
germ50  <- read.csv( "data/pars_gsim50.csv" )
germ75  <- read.csv( "data/pars_gsim75.csv" )
germ100 <- read.csv( "data/pars_gsim100.csv" )
germ200 <- read.csv( "data/pars_gsim200.csv" )
germ500 <- read.csv( "data/pars_gsim500.csv" )

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

cov2_plot <- read.csv( "data/cov_plot2.csv" )
cov3_plot <- read.csv( "data/cov_plot3.csv" )


# Setup ------------------------------------------------------------------------

base_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    plot.margin = margin(5.5, 10, 5.5, 5.5)
  )


# Figure 1 ---------------------------------------------------------------------

# Adapted from plot( uncert2, type = "param" ) to add lettering

# fig1 <- plot( uncert2, type = "param" ) +
#   theme( plot.subtitle = element_text( size = 8, face = "bold", 
#                                        hjust = 0.2, vjust = 0.8 ),
#          plot.margin = margin( 10,10,5.5,5.5))


default_cols <- c(
  "#D0873C", # soft orange
  "#8FB339", # olive green
  "#7A5195", # muted purple
  "#88CCEE" # light blue
)

vr_levels <- unique(uncert2$params_uncert$vital_rate)

vr_colors <- setNames( rep(default_cols, length.out = length(vr_levels)), 
                       vr_levels )

df <- uncert2$params_uncert

df$parameter <- factor( df$parameter, 
                        levels = rev(df$parameter[order(df$vital_rate)]) )

title_wrap <- function(x) str_wrap(x, width = 26)

p1 <- ggplot(df, aes(x = variance, y = parameter, color = vital_rate)) +
  geom_vline(xintercept = 0) +
  geom_point(size = 2.5) +
  geom_segment(aes(x = 0, xend = variance, yend = parameter), linewidth = 1) +
  scale_x_sqrt() +
  labs(x = "Variance",
       y = "Parameter",
       title = title_wrap("(a) Parameter uncertainty")) +
  scale_color_manual(values = vr_colors, guide = "none", 
                     labels = c( "Growth", "Recruitment", 
                                 "Reproduction", "Survival" )) +
  base_theme

p2 <- ggplot(df, aes(x = sensitivity, y = parameter, color = vital_rate)) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 1, color = "gray50", lty = 2, lwd = 0.3 ) +
  geom_vline(xintercept = -1, color = "gray50", lty = 2, lwd = 0.3) +
  geom_point(size = 2.5) +
  geom_segment(aes(x = 0, xend = sensitivity, yend = parameter), 
               linewidth = 1) +
  labs(x = "Sensitivity",
       y = NULL,
       title = title_wrap("(b) Model sensitivity to parameter")) +
  scale_color_manual(values = vr_colors, guide = "none",
                     labels = c( "Growth", "Recruitment", 
                                 "Reproduction", "Survival" )) +
  base_theme

p3 <- ggplot(df, aes(x = variance_prop, y = parameter, color = vital_rate)) +
  geom_vline(xintercept = 0) +
  geom_point(size = 2.5) +
  geom_segment(aes(x = 0, xend = variance_prop, yend = parameter), 
               linewidth = 1) +
  labs(x = "Uncertainty contribution",
       y = NULL,
       title = title_wrap("(c) Parameter contribution to uncertainty in \u03BB"),
       color = "Vital rate") +
  scale_color_manual(values = vr_colors, labels = c( "Growth", "Recruitment", 
                                                     "Reproduction", "Survival" )) +
  base_theme

fig1 <- p1 + p2 + p3 + plot_layout(ncol = 3, axes = 'collect_y') & 
  theme(plot.title = element_text(hjust = 0.5, size = 8),
        axis.text.x = element_text(size = 6))


# Figure 2 --------------------------------------------------------------------

# Figure 2a: comparison of explained uncertainty with recruitment parameters
  # varying and constant

fig2a <- ggplot( ) + 
  geom_bar( data = uncert_comp_plot[which(uncert_comp_plot$vital_rate != "total" ),], 
            aes( x = mod, y = variance_sum, fill = vital_rate ), 
            stat = "identity", position = "stack", width = 0.7 ) +
  geom_bar( data = uncert_comp_plot[which(uncert_comp_plot$vital_rate == "total" ),], 
            aes( x = mod, y = variance_sum, alpha = "Model uncertainty" ), 
            stat = "identity", position = "stack", color = "black", fill = NA,
            linewidth = 0.8, width = 0.7 ) +
  scale_y_sqrt() +
  scale_fill_manual( values = c( "#8FB339", "#88CCEE", "#7A5195", "#D0873C" ),
                     labels = c( "Growth", "Recruitment", 
                                 "Reproduction", "Survival" ),
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


# Figure 2b: comparison of sampled and mean lambda values, with recruitment
  # parameters varying and constant

lam_t <- t.test( lambda ~ mod, data = lam_tab_comp )

mean_lam_plot <- mean_lam[c(1,1),]
colnames( mean_lam_plot ) <- c( "mod", "type", "value" )
mean_lam_plot[1,1] <- "Varying"
mean_lam_plot[2,1] <- "Constant"

mean_lam_v <- mean( lam_tab_comp[which(lam_tab_comp$mod == "Varying"),'lambda'])
mean_lam_c <- mean( lam_tab_comp[which(lam_tab_comp$mod == "Constant"),'lambda'])

fig2b <- ggplot() +
  geom_violin( data = lam_tab_comp, aes( x = factor( mod ), y = lambda ),
               color = NA, alpha = 0.2, fill = "black" ) +
  geom_segment( aes( x = 2.45, y = mean_lam_v, 
                     xend = 1.55, yend = mean_lam_v,
                     alpha = "Mean of samples" ),
                size = 0.9 ) +
  geom_segment( aes( x = 0.55, y = mean_lam_c, 
                     xend = 1.45, yend = mean_lam_c ),
                size = 0.9 ) +
  geom_point( data = mean_lam_plot, aes( x = factor( mod ), y = value, 
                                             alpha = "Mean model" ),
              cex = 3, pch = 16 ) +
  geom_signif( data = lam_tab_comp, aes( x = factor( mod ), y = lambda ),
               comparisons = list( c("Varying","Constant")),
               map_signif_level = TRUE ) +
  scale_x_discrete( labels = c( "Constant", "Varying (resampled)" ) ) +
  ylim( 1.0, 1.8 ) +
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

fig2 <- fig2a + fig2b + plot_layout( ncol = 1, axes = "collect_x" )


# Figure 3 ---------------------------------------------------------------------

fig3 <- ggplot() +
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
  scale_y_sqrt() +
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
                aes( x = log_area_t0, y = surv_t1 ), 
                color = "#D0873C", cex = 2 ) +
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
  p_surv_pred[[3]] + p_surv_pred[[4]] + plot_layout( ncol = 2, 
                                                     axes = "collect" )


# Figure S4 --------------------------------------------------------------------

pars_var2   <- c( "surv_b0", "surv_b1", "surv_b2", 
                  "grow_b0", "grow_b1", "grow_sig",
                  "abort", "clip",
                  "flow_b0", "flow_b1",
                  "fert_b0", "fert_b1",
                  "g0", "g1", "g2" )

s_pars_plot <- s_pars2[,pars_var2]

s_pars_plot <- pivot_longer( s_pars_plot, cols = pars_var2,
                             names_to = "par")

s_pars_mean <- s_pars_plot %>% group_by( par ) %>% 
  summarise( mean = mean( value ), .groups = "drop" )

pars_mean <- pars_mean2[,pars_var2]
pars_mean <- pivot_longer( pars_mean, cols = pars_var2,
                           names_to = "par")

vr_key <- c(
  surv_b0 = "Survival",
  surv_b1 = "Survival",
  surv_b2 = "Survival",
  grow_b0 = "Growth",
  grow_b1 = "Growth",
  grow_sig = "Growth",
  abort = "Reproduction",
  clip = "Reproduction",
  flow_b0 = "Reproduction",
  flow_b1 = "Reproduction",
  fert_b0 = "Reproduction",
  fert_b1 = "Reproduction",
  g0 = "Recruitment",
  g1 = "Recruitment",
  g2 = "Recruitment"
)

figs4 <- ggplot() +
  geom_violin( data = s_pars_plot, aes( x = par, y = value, 
                                        fill = vr_key[par] ),
               color = NA, alpha = 0.6 ) +
  facet_wrap( ~ par, scales = "free" ) +
  geom_point( data = pars_mean, aes( x = par, y = value,
                                     alpha = "Mean model" ),
              cex = 2, pch = 16, color = "black" ) +
  geom_segment( data = s_pars_mean, aes( x = 0.5, xend = 1.5,
                                         y = mean, yend = mean,
                                         alpha = "Mean of samples") ) +
  scale_alpha_manual( name = NULL,
                      values = c( 1, 1 ),
                      breaks = c( "Mean model", "Mean of samples" ),
                      guide = guide_legend( override.aes = list( linetype = c(0,1),
                                                                 shape = c(16,NA),
                                                                 color = "black" ) ) ) +
  scale_fill_manual( name = "Vital rate", values = c( Survival = "#D0873C",
                                                      Growth = "#8FB339",
                                                      Reproduction = "#7A5195",
                                                      Recruitment = "#88CCEE" ) ) +
  labs( x = "Parameter",
        y = "Value" ) +
  base_theme +
  theme( strip.text.x = element_blank() )


# Figure S5 --------------------------------------------------------------------

figs5 <- ggplot( corr2_plot, aes( x = Var1, y = Var2, fill = correlation ) ) + 
  geom_tile() +
  geom_text( aes( label = text ) ) +
  base_theme +
  scale_fill_gradient2( low = "#D55E00", mid = "white", high = "#0072B2",
                        midpoint = 0, limits = c(-1,1), 
                        name = "Correlation" ) +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1 ) ) +
  annotate( "rect", xmin = 0.5, xmax = 6.5, ymin = 15.7, ymax = 16.0,
            fill = "#7A5195" ) +
  annotate( "rect", xmin = 6.5, xmax = 9.5, ymin = 15.7, ymax = 16.0,
            fill = "#88CCEE" ) +
  annotate( "rect", xmin = 9.5, xmax = 12.5, ymin = 15.7, ymax = 16.0,
            fill = "#8FB339" ) +
  annotate( "rect", xmin = 12.5, xmax = 15.5, ymin = 15.7, ymax = 16.0,
            fill = "#D0873C" ) +
  annotate( "rect", xmin = 15.7, xmax = 16.0, ymin = 0.5, ymax = 6.5,
            fill = "#7A5195" ) +
  annotate( "rect", xmin = 15.7, xmax = 16.0, ymin = 6.5, ymax = 9.5,
            fill = "#88CCEE" ) +
  annotate( "rect", xmin = 15.7, xmax = 16.0, ymin = 9.5, ymax = 12.5,
            fill = "#8FB339" ) +
  annotate( "rect", xmin = 15.7, xmax = 16.0, ymin = 12.5, ymax = 15.5, 
            fill = "#D0873C" ) +
  geom_vline( xintercept = c( 6.5, 9.5, 12.5 ), color = "white", 
              linewidth = 0.7 ) +
  geom_hline( yintercept = c( 6.5, 9.5, 12.5 ), color = "white", 
              linewidth = 0.7 )


# Figure S6 --------------------------------------------------------------------

corr_plot_fxn <- function( df, caption ){
  
  temp_plot <- ggplot( df, aes( x = Var1, y = Var2, fill = mean_corr ) ) + 
    geom_tile() +
    base_theme +
    labs( title = caption ) + 
    scale_fill_gradient2( low = "#D55E00", mid = "white", high = "#0072B2",
                          midpoint = 0, limits = c(-1,1), 
                          name = "Correlation" ) +
    theme( axis.title.x = element_blank(),
           axis.title.y = element_blank(),
           axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1 ) ) +
    annotate( "rect", xmin = 0.5, xmax = 6.5, ymin = 15.7, ymax = 16.0,
              fill = "#7A5195" ) +
    annotate( "rect", xmin = 6.5, xmax = 9.5, ymin = 15.7, ymax = 16.0,
              fill = "#88CCEE" ) +
    annotate( "rect", xmin = 9.5, xmax = 12.5, ymin = 15.7, ymax = 16.0,
              fill = "#8FB339" ) +
    annotate( "rect", xmin = 12.5, xmax = 15.5, ymin = 15.7, ymax = 16.0,
              fill = "#D0873C" ) +
    annotate( "rect", xmin = 15.7, xmax = 16.0, ymin = 0.5, ymax = 6.5,
              fill = "#7A5195" ) +
    annotate( "rect", xmin = 15.7, xmax = 16.0, ymin = 6.5, ymax = 9.5,
              fill = "#88CCEE" ) +
    annotate( "rect", xmin = 15.7, xmax = 16.0, ymin = 9.5, ymax = 12.5,
              fill = "#8FB339" ) +
    annotate( "rect", xmin = 15.7, xmax = 16.0, ymin = 12.5, ymax = 15.5, 
              fill = "#D0873C" ) +
    geom_vline( xintercept = c( 6.5, 9.5, 12.5 ), color = "white", 
                linewidth = 0.7 ) +
    geom_hline( yintercept = c( 6.5, 9.5, 12.5 ), color = "white", 
                linewidth = 0.7 )
  
  return( temp_plot )
}


# Correlation plots for the simulated germination datasets

figs6a <- corr_plot_fxn( corr_plot6,
                         caption = "(a) 6 samples" ) +
  theme( legend.position = "none" )

figs6b <- corr_plot_fxn( corr_plot10,
                         caption = "(b) 10 samples" ) +
  theme( legend.position = "none" )

figs6c <- corr_plot_fxn( corr_plot20,
                         caption = "(c) 20 samples" ) +
  theme( legend.position = "none" )

figs6d <- corr_plot_fxn( corr_plot30,
                         caption = "(d) 30 samples" ) +
  theme( legend.position = "none" )

figs6e <- corr_plot_fxn( corr_plot40,
                         caption = "(e) 40 samples" ) +
  theme( legend.position = "none" )

figs6f <- corr_plot_fxn( corr_plot50,
                         caption = "(f) 50 samples" ) +
  theme( legend.position = "none" )

figs6g <- corr_plot_fxn( corr_plot75,
                         caption = "(g) 75 samples" ) +
  theme( legend.position = "none" )

figs6h <- corr_plot_fxn( corr_plot100, 
                         caption = "(h) 100 samples" ) +
  theme( legend.position = "none" )

figs6i <- corr_plot_fxn( corr_plot200,
                         caption = "(i) 200 samples" ) +
  theme( legend.position = "none" )

figs6j <- corr_plot_fxn( corr_plot500,
                         caption = "(j) 500 samples" ) +
  theme( legend.position = "none" )

figs6 <- figs6a + figs6b + figs6c + figs6d + figs6e + figs6f + figs6g + figs6h + 
  figs6i + figs6j + plot_layout( ncol = 4 )


# Figure S7 --------------------------------------------------------------------

# Comparing sampled germination parameter values to mean model values

g_format <- function( df, n ){
  df_temp <- df[,c("g0", "g1", "g2")]
  df_temp <- pivot_longer( df_temp, cols = c("g0", "g1", "g2"),
                           names_to = "par" )
  df_temp$type <- n
  
  return( df_temp )
}

gtab6   <- g_format( germ6, 6 )
gtab10  <- g_format( germ10, 10 )
gtab20  <- g_format( germ20, 20 )
gtab30  <- g_format( germ30, 30 )
gtab40  <- g_format( germ40, 40 )
gtab50  <- g_format( germ50, 50 )
gtab75  <- g_format( germ75, 75 )
gtab100 <- g_format( germ100, 100 )
gtab200 <- g_format( germ200, 200 )
gtab500 <- g_format( germ500, 500 )

gtabmean <- pars_mean2[,c("g0", "g1", "g2")]

gtabmean <- pivot_longer( gtabmean, cols = c("g0", "g1", "g2"),
                          names_to = "par")

gtab <- bind_rows( gtab6, gtab10, gtab20, gtab30, gtab40, gtab50, gtab75,
                   gtab100, gtab200, gtab500 )

gtab_sample_mean <- gtab %>% group_by( type, par ) %>% 
  summarise( mean = mean( value ), .groups = "drop" )

par_levels <- unique(gtab$par)

gtab <- gtab %>% mutate(par = factor(par, 
                                     levels = par_levels, 
                                     labels = c( expression('g'[0]),
                                                 expression('g'[1]),
                                                 expression('g'[2]) ) ) )

gtabmean <- gtabmean %>% mutate(par = factor(par, 
                                             levels = par_levels, 
                                             labels = c( expression('g'[0]),
                                                         expression('g'[1]),
                                                         expression('g'[2]) )))

gtab_sample_mean <- gtab_sample_mean %>% mutate( par = factor(par, 
                                                              levels = par_levels, 
                                                              labels = c( expression('g'[0]),
                                                                          expression('g'[1]),
                                                                          expression('g'[2]) )) )



nice_top <- function(x) {
  magnitude <- 10^floor(log10(x))
  scaled <- x / magnitude
  
  multiplier <- dplyr::case_when(
    scaled <= 1   ~ 1.5,
    scaled <= 1.5 ~ 2,
    scaled <= 2   ~ 2.5,
    scaled <= 2.5 ~ 3,
    scaled <= 3   ~ 4,
    scaled <= 4   ~ 5,
    scaled <= 5   ~ 6,
    scaled <= 6   ~ 7.5,
    scaled <= 7.5 ~ 8,
    scaled <= 8   ~ 9,
    TRUE          ~ 10
  )
  
  multiplier * magnitude
}


g0_scale <- scale_y_continuous( breaks = function(x){
  top <- nice_top(max(x))
  c(0,top/2, top)
},
labels = scales::label_comma(),
expand = expansion(mult=c(0,0.05)),
sec.axis = sec_axis(~ . , name = "Simulated sample size",
                    breaks = NULL, labels = NULL))

g1_scale <- scale_y_continuous( breaks = c(0, 30000),
                                limits = c(0,30000),
                                labels = scales::label_comma(),
                                sec.axis = sec_axis(~ . , name = "Simulated sample size",
                                                    breaks = NULL, labels = NULL))

g2_scale <- scale_y_continuous( breaks = c(0, 60000),
                                limits = c(0, 60000),
                                labels = scales::label_comma(),
                                sec.axis = sec_axis(~ . , name = "Simulated sample size",
                                                    breaks = NULL, labels = NULL))

g0_x <- scale_x_continuous( limits = c(0.01,0.03),
                            breaks = c(0.01, 0.02, 0.03),
                            oob = scales::oob_keep,
                            sec.axis = sec_axis(~ . , name = "Recruitment parameter", 
                                                breaks = NULL, labels = NULL ) )

g1_x <- scale_x_continuous( limits = c(0,0.04),
                            breaks = c(0,0.02,0.04),
                            oob = scales::oob_keep,
                            sec.axis = sec_axis(~ . , name = "Recruitment parameter", 
                                                breaks = NULL, labels = NULL ) )

g2_x <- scale_x_continuous( limits = c(0,0.0022),
                            n.breaks = 3,
                            oob = scales::oob_keep,
                            sec.axis = sec_axis(~ . , name = "Recruitment parameter", 
                                                breaks = NULL, labels = NULL ) )

y_scales <- rep(list(g0_scale, g1_scale, g2_scale), times = 10)

x_scales <- rep(list(g0_x, g1_x, g2_x), times = 10)


figs7 <- ggplot() +
  geom_histogram( data = gtab, aes( x = value ),
                  color = NA, fill = "#88CCEE", bins = 50 ) +
  facet_grid2( type ~ par, scales = "free", labeller = label_parsed, 
               independent = "y" ) +
  geom_vline( data = gtabmean, aes( xintercept = value,
                                    alpha = "Mean model" ),
              linewidth = 0.8, linetype = 3, color = "black" ) +
  geom_vline( data = gtab_sample_mean,
                aes( xintercept = mean,
                     alpha = "Mean of samples" ),
                linewidth = 0.8, color = "black" ) +
  scale_alpha_manual( name = NULL,
                      values = c( 1, 1 ),
                      breaks = c( "Mean model", "Mean of samples" ),
                      guide = guide_legend( override.aes = list( linetype = c(3,1),
                                                                 color = "black" ) ) ) +
  labs( x = "Value", 
        y = "Frequency" ) +
  facetted_pos_scales( x = x_scales, y = y_scales ) +
  scale_y_continuous( sec.axis = sec_axis(~ . , name = "Simulated sample size",
                                          breaks = NULL, labels = NULL) ) +
  scale_x_continuous( sec.axis = sec_axis(~ . , name = "Recruitment parameter", 
                                          breaks = NULL, labels = NULL ) ) +
  base_theme +
  theme( strip.text = element_text( face = "bold", hjust = 0, size = 10 ),
         strip.text.x = element_text( hjust = 0.5),
         strip.text.y = element_text(angle = 0),
         axis.title = element_text(size = 10, face = "bold"),
         panel.spacing.y = unit( 0.9, "lines" ) )


# Figure S8 -------------------------------------------------------------------

model_labs <- c( "Quadratic", "Cubic" )

mean_lam_plot2 <- mean_lam[c(1:4),]

# Figure S8a

figs8a <- ggplot() +
  geom_violin( data = lam_tab_23, aes( x = factor( type ), y = lambda,
                                    fill = factor( type ) ),
               color = NA, alpha = 0.4 ) +
  geom_signif( data = lam_tab_23, aes( x = factor( type ), y = lambda ),
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
  ylim( min( lam_tab_23$lambda ) - 0.03, max( lam_tab_23$lambda ) + 0.1 ) +
  labs( x = "Survival model", 
        y = "Lambda",
        title = "(a) Varying recruitment" ) +
  base_theme

# Figure S8b

figs8b <- ggplot() +
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
  ylim( min( lam_tab_ng$lambda ) - 0.03, max( lam_tab_ng$lambda ) + 0.05 ) +
  labs( x = "Survival model", 
        y = "Lambda",
        title = "(b) Constant recruitment" ) +
  base_theme

# Figure S8c

figs8c <- ggplot( ) + 
  geom_bar( data = var_tab[which(var_tab$vital_rate != "total" ),], 
            aes( x = type, y = variance_sum, fill = vital_rate ), 
            stat = "identity", position = "stack", width = 0.7 ) +
  geom_bar( data = var_tab[which(var_tab$vital_rate == "total" ),], 
            aes( x = type, y = variance_sum, alpha = "Model uncertainty" ), 
            stat = "identity", position = "stack", color = "black", fill = NA,
            linewidth = 0.8, width = 0.7 ) +
  scale_fill_manual( values = c( "#8FB339", "#88CCEE", "#7A5195", "#D0873C" ),
                     labels = c( "Growth", "Recruitment", 
                                 "Reproduction", "Survival" ),
                     guide = "legend" ) +
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
        title = "(c) Varying recruitment",
        fill = "Vital rate contribution" ) +
  base_theme

# Figure S8d

figs8d <- ggplot( ) + 
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
  scale_y_sqrt() +
  scale_alpha_manual( name = NULL,
                      values = 1,
                      breaks = "Model uncertainty",
                      guide = guide_legend( override.aes = list( color = "black" ) )
  ) +
  guides( fill = guide_legend( "Vital rate contribution" ) ) +
  labs( x = "Survival model", 
        y = "Uncertainty",
        title = "(d) Constant recruitment",
        fill = "Vital rate contribution" ) +
  base_theme +
  theme( legend.position = "none" )

# Figure S8e-f

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

figs8_base <- ggplot() + 
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
  base_theme

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

# Figure S8e

figs8e <- figs8_base +
  geom_text( data = model_AIC[which( model_AIC$model == 2 ),], 
             aes( x = size, y = survival, label = AIC ) ) +
  geom_line( data = surv_pred2, 
             aes( x = log_area_t0,
                  y = surv_t1 ),
             color = "#f5a351", lwd = 1.5, alpha = 0.8 ) +
  labs( title = "(e) Quadratic survival model" )

# Figure S8f

figs8f <- figs8_base +
  geom_text( data = model_AIC[which( model_AIC$model == 3 ),], 
             aes( x = size, y = survival, label = AIC ) ) +
  geom_line( data = surv_pred3, 
             aes( x = log_area_t0,
                  y = surv_t1 ),
             color = "#a85603", lwd = 1.5, alpha = 0.8 ) +
  labs( title = "(f) Cubic survival model" )


# Figure S8, all together

figs8ab <- figs8a + figs8b + plot_layout( ncol = 1, axes = "collect", 
                                               guides = "collect" )

figs8cd <- figs8c + figs8d + plot_layout( ncol = 1, axes = "collect", 
                                               guides = "collect" )

figs8ef <- figs8e + figs8f + plot_layout( ncol = 2, axes = "collect" )

figs8abcd <- wrap_plots( figs8ab ) + wrap_plots( figs8cd )

figs8 <- wrap_plots( figs8abcd ) / wrap_plots( figs8ef ) + 
  plot_layout( heights = c( 3, 2 ) )


# Figure S9 --------------------------------------------------------------------

figs9a <- ggplot( corr2_plot, aes( x = Var1, y = Var2, fill = correlation ) ) + 
  geom_tile() +
  geom_text( aes( label = text ) ) +
  base_theme +
  labs( title = "(a) With quadratic survival model" ) + 
  scale_fill_gradient2( low = "#D55E00", mid = "white", high = "#0072B2",
                        midpoint = 0, limits = c(-1,1), name = "Correlation" ) +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1 ) ) +
  annotate( "rect", xmin = 0.5, xmax = 6.5, ymin = 15.7, ymax = 16.0,
            fill = "#7A5195" ) +
  annotate( "rect", xmin = 6.5, xmax = 9.5, ymin = 15.7, ymax = 16.0,
            fill = "#88CCEE" ) +
  annotate( "rect", xmin = 9.5, xmax = 12.5, ymin = 15.7, ymax = 16.0,
            fill = "#8FB339" ) +
  annotate( "rect", xmin = 12.5, xmax = 15.5, ymin = 15.7, ymax = 16.0,
            fill = "#D0873C" ) +
  annotate( "rect", xmin = 15.7, xmax = 16.0, ymin = 0.5, ymax = 6.5,
            fill = "#7A5195" ) +
  annotate( "rect", xmin = 15.7, xmax = 16.0, ymin = 6.5, ymax = 9.5,
            fill = "#88CCEE" ) +
  annotate( "rect", xmin = 15.7, xmax = 16.0, ymin = 9.5, ymax = 12.5,
            fill = "#8FB339" ) +
  annotate( "rect", xmin = 15.7, xmax = 16.0, ymin = 12.5, ymax = 15.5, 
            fill = "#D0873C" ) +
  geom_vline( xintercept = c( 6.5, 9.5, 12.5 ), color = "white", 
              linewidth = 0.7 ) +
  geom_hline( yintercept = c( 6.5, 9.5, 12.5 ), color = "white", 
              linewidth = 0.7 )


figs9b <- ggplot( corr3_plot, aes( x = Var1, y = Var2, fill = correlation ) ) + 
  geom_tile() +
  geom_text( aes( label = text ) ) +
  labs( title = "(b) With cubic survival model" ) + 
  base_theme + 
  scale_fill_gradient2( low = "#D55E00", mid = "white", high = "#0072B2",
                        midpoint = 0, limits = c(-1,1), name = "Correlation" ) +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1 ) ) +
  annotate( "rect", xmin = 0.5, xmax = 6.5, ymin = 16.7, ymax = 17.0,
            fill = "#7A5195" ) +
  annotate( "rect", xmin = 6.5, xmax = 9.5, ymin = 16.7, ymax = 17.0,
            fill = "#88CCEE" ) +
  annotate( "rect", xmin = 9.5, xmax = 12.5, ymin = 16.7, ymax = 17.0,
            fill = "#8FB339" ) +
  annotate( "rect", xmin = 12.5, xmax = 16.5, ymin = 16.7, ymax = 17.0,
            fill = "#D0873C" ) +
  annotate( "rect", xmin = 16.7, xmax = 17.0, ymin = 0.5, ymax = 6.5,
            fill = "#7A5195" ) +
  annotate( "rect", xmin = 16.7, xmax = 17.0, ymin = 6.5, ymax = 9.5,
            fill = "#88CCEE" ) +
  annotate( "rect", xmin = 16.7, xmax = 17.0, ymin = 9.5, ymax = 12.5,
            fill = "#8FB339" ) +
  annotate( "rect", xmin = 16.7, xmax = 17.0, ymin = 12.5, ymax = 16.5, 
            fill = "#D0873C" ) +
  geom_vline( xintercept = c( 6.5, 9.5, 12.5 ), color = "white", 
              linewidth = 0.7 ) +
  geom_hline( yintercept = c( 6.5, 9.5, 12.5 ), color = "white", 
              linewidth = 0.7 )


figs9ab <- figs9a / figs9b + plot_layout( guides = "collect" )


# Covariances instead of correlations

figs9c <- ggplot( cov2_plot, aes( x = Var1, y = Var2, fill = covariance ) ) + 
  geom_tile() +
  base_theme +
  labs( title = "(c) With quadratic survival model" ) + 
  scale_fill_gradient2( low = "#D55E00", mid = "white", high = "#0072B2",
                        midpoint = 0, name = "Covariance",
                        limits = c( -0.09, 0.11 ) ) +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1 ) ) +
  annotate( "rect", xmin = 0.5, xmax = 6.5, ymin = 15.7, ymax = 16.0,
            fill = "#7A5195" ) +
  annotate( "rect", xmin = 6.5, xmax = 9.5, ymin = 15.7, ymax = 16.0,
            fill = "#88CCEE" ) +
  annotate( "rect", xmin = 9.5, xmax = 12.5, ymin = 15.7, ymax = 16.0,
            fill = "#8FB339" ) +
  annotate( "rect", xmin = 12.5, xmax = 15.5, ymin = 15.7, ymax = 16.0,
            fill = "#D0873C" ) +
  annotate( "rect", xmin = 15.7, xmax = 16.0, ymin = 0.5, ymax = 6.5,
            fill = "#7A5195" ) +
  annotate( "rect", xmin = 15.7, xmax = 16.0, ymin = 6.5, ymax = 9.5,
            fill = "#88CCEE" ) +
  annotate( "rect", xmin = 15.7, xmax = 16.0, ymin = 9.5, ymax = 12.5,
            fill = "#8FB339" ) +
  annotate( "rect", xmin = 15.7, xmax = 16.0, ymin = 12.5, ymax = 15.5, 
            fill = "#D0873C" ) +
  geom_vline( xintercept = c( 6.5, 9.5, 12.5 ), color = "white", 
              linewidth = 0.7 ) +
  geom_hline( yintercept = c( 6.5, 9.5, 12.5 ), color = "white", 
              linewidth = 0.7 )


figs9d <- ggplot( cov3_plot, aes( x = Var1, y = Var2, fill = covariance ) ) + 
  geom_tile() +
  labs( title = "(d) With cubic survival model" ) + 
  base_theme + 
  scale_fill_gradient2( low = "#D55E00", mid = "white", high = "#0072B2",
                        midpoint = 0, name = "Covariance",
                        limits = c( -0.09, 0.11 ) ) +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1 ) ) +
  annotate( "rect", xmin = 0.5, xmax = 6.5, ymin = 16.7, ymax = 17.0,
            fill = "#7A5195" ) +
  annotate( "rect", xmin = 6.5, xmax = 9.5, ymin = 16.7, ymax = 17.0,
            fill = "#88CCEE" ) +
  annotate( "rect", xmin = 9.5, xmax = 12.5, ymin = 16.7, ymax = 17.0,
            fill = "#8FB339" ) +
  annotate( "rect", xmin = 12.5, xmax = 16.5, ymin = 16.7, ymax = 17.0,
            fill = "#D0873C" ) +
  annotate( "rect", xmin = 16.7, xmax = 17.0, ymin = 0.5, ymax = 6.5,
            fill = "#7A5195" ) +
  annotate( "rect", xmin = 16.7, xmax = 17.0, ymin = 6.5, ymax = 9.5,
            fill = "#88CCEE" ) +
  annotate( "rect", xmin = 16.7, xmax = 17.0, ymin = 9.5, ymax = 12.5,
            fill = "#8FB339" ) +
  annotate( "rect", xmin = 16.7, xmax = 17.0, ymin = 12.5, ymax = 16.5, 
            fill = "#D0873C" ) +
  geom_vline( xintercept = c( 6.5, 9.5, 12.5 ), color = "white", 
              linewidth = 0.7 ) +
  geom_hline( yintercept = c( 6.5, 9.5, 12.5 ), color = "white", 
              linewidth = 0.7 )


figs9cd <- figs9c / figs9d + plot_layout( guides = "collect" )

figs9 <- wrap_plots( figs9ab ) + wrap_plots( figs9cd )


# Save output ------------------------------------------------------------------

ggsave( "results/Figure1.pdf", fig1, width = 180, height = 100, units = "mm",
        device = cairo_pdf )
ggsave( "results/Figure2.pdf", fig2, width = 130, height = 150, units = "mm",
        device = cairo_pdf )
ggsave( "results/Figure3.pdf", fig3, width = 180, height = 100, units = "mm",
        device = cairo_pdf )


ggsave( "results/FigureS2.pdf", figs2, width = 220, height = 110, units = "mm",
        device = cairo_pdf )
ggsave( "results/FigureS3.pdf", figs3, width = 200, height = 160, units = "mm",
        device = cairo_pdf )
ggsave( "results/FigureS4.pdf", figs4, width = 220, height = 160, units = "mm",
        device = cairo_pdf )
ggsave( "results/FigureS5.pdf", figs5, width = 160, height = 110, units = "mm",
        device = cairo_pdf )
ggsave( "results/FigureS6.pdf", figs6, width = 280, height = 240, units = "mm",
        device = cairo_pdf )
ggsave( "results/FigureS7.pdf", figs7, width = 260, height = 220, units = "mm",
        device = cairo_pdf )
ggsave( "results/FigureS8.pdf", figs8, width = 220, height = 180, units = "mm",
        device = cairo_pdf )
ggsave( "results/FigureS9.pdf", figs9, width = 280, height = 220, units = "mm",
        device = cairo_pdf )

