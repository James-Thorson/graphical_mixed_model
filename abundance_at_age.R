
# devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)', force = TRUE, dep = FALSE )
root_dir = R'(C:\Users\James.Thorson\Desktop\Git\GGMM)'
data_dir = file.path(root_dir, "data" )

library(tinyVAST)
library(cv)

dataset = c( "ages", "lengths" )[1]

Date = Sys.Date()
date_dir = file.path(root_dir, Date )
  dir.create( date_dir )

#######################
# LENGTHS
#######################

if( dataset == "lengths" ){
  out = scan( file.path(data_dir,"GOA_Rex_8_2021_run7.dat"), skip = 205, nlines = 17 )
  #length(out) / 64 # 29 length bins
  mat = matrix( out, ncol = 64, byrow = TRUE )
  #mat[,1]

  comps = mat[,6+1:29]
  comps = sweep( comps, MARGIN=1, STATS = rowSums(comps), FUN = "/" )
  dimnames(comps) = list( 'year' = mat[,1], 'bin' = paste0("bin_",seq_len(ncol(comps))) )
}

#######################
# AGES
#######################

if( dataset == "ages" ){
  out = scan( file.path(data_dir,"GOA_Rex_8_2021_run7.dat"), skip = 263, nlines = 18 )
  #length(out) / 49
  mat = matrix( out, ncol = 49, byrow = TRUE )
  #mat[,1]

  comps = mat[,9+1:20]
  comps = sweep( comps, MARGIN=1, STATS = rowSums(comps), FUN = "/" )
  dimnames(comps) = list( 'year' = mat[,1], 'bin' = paste0("bin_",1:20) )
}

#######################
# RUN
#######################

#
which_cols = which( colSums(comps) > 0 )
comps = comps[,which_cols]
bins = colnames(comps)

#
long = expand.grid( dimnames(comps) )
long$prop = as.vector(comps)
long$year = as.numeric(as.character(long$year))
long$bin = factor( long$bin, levels = bins )

#
cohorts_rw = paste0( bins[-length(bins)], " -> ", bins[-1], ", 1, NA, 1" )
sds = paste0( bins, " <-> ", bins, ", 0, sd" )
bins_est = paste0( bins[-length(bins)], " -> ", bins[-1], ", 0, age" )
years_est = paste0( bins, " -> ", bins, ", 1, year" )
cohorts_est = paste0( bins[-length(bins)], " -> ", bins[-1], ", 1, cohorts" )
mods = list(
  "null" = paste0( c(sds), collapse = "\n" ),
  #"cohort_RW" = paste0( c(cohorts_rw, sds), collapse = "\n" ),
  "bins" = paste0( c(bins_est, sds), collapse = "\n" ),
  "years" = paste0( c(years_est, sds), collapse = "\n" ),
  "cohorts" = paste0( c(cohorts_est, sds), collapse = "\n" ),
  "bin_year" = paste0( c(bins_est, years_est, sds), collapse = "\n" ),
  "cohort_year" = paste0( c(cohorts_est, years_est, sds), collapse = "\n" ),
  "bin_cohort" = paste0( c(bins_est, cohorts_est, sds), collapse = "\n" ),
  "bin_year_cohort" = paste0( c(bins_est, years_est, cohorts_est, sds), collapse = "\n" )
)

fit = CV = NULL
for( index in seq_along(mods) ){
  message( Sys.time(), ": running model ", names(mods)[index] )
  #
  fit[[index]] = tinyVAST(
    data = long,
    formula = prop ~ 0 + bin,
    family = tweedie(),
    time_column = "year",
    variable_column = "bin",
    time_term = mods[[index]],
    control = tinyVASTcontrol( trace = 0 )
  )

  #
  CV[[index]] = cv( fit[[index]] )
}

#
Table_S1 = data.frame(
  model = names(mods),
  num_par = sapply( fit, FUN = \(x) length(x$opt$par) ),
  RMSE = sqrt(sapply( CV, cvInfo )),
  Prop_var = 1 - sapply( CV, cvInfo ) / cvInfo( CV[[1]] ),
  cAIC = sapply( fit, cAIC ),
  AIC = sapply( fit, AIC )
)
Table_S1[,4] = sapply( Table_S1[,4], FUN = round, digits = 2 )
Table_S1[,3] = sapply( Table_S1[,3], FUN = round, digits = 3 )
Table_S1[,5:6] = sapply( Table_S1[,5:6], FUN = round, digits = 1 )
write.csv( Table_S1, file = file.path(date_dir,"Table_S1.csv") )

if(dataset=="ages") which_best = 6

#
make_matrix = function( dimnames ){matrix(NA, nrow=length(dimnames[[1]]), ncol=length(dimnames[[2]]), dimnames=dimnames)}
p_ct = make_matrix( list('year' = 1992:2022, 'bin' = levels(long$bin)) )
newdata = expand.grid( dimnames(p_ct) )
#newdata$binnum = paste0( "bin_", newdata$bin )
#newdata$binnum = factor( newdata$binnum, levels = paste0("bin_", unique(newdata$bin)) )

# Do prediction
#mu_gz = sample_variable( fit, newdata = newdata, n_samples = 1000, sample_fixed = FALSE, variable_name = "mu_g" )
#p_ct[] = rowMeans(mu_gz)
p_ct[] = predict( fit[[which_best]], newdata = newdata ) #, se.fit = TRUE, bias.correct = TRUE )
#p_ct = sweep( p_ct, MARGIN=1, STATS = rowSums(p_ct), FUN = "/" )
newdata$phat = as.vector(p_ct)

#
drop_years = min(long$year):max(long$year)
fit_z = NULL
newdata$p_cv = NULL
for( z in seq_along(drop_years) ){
  message( Sys.time(), ": running model ", drop_years[z] )
  long_z = subset( long, year != drop_years[z] )
  long_z = droplevels( long_z )

  fit_z = tinyVAST(
    data = long_z,
    formula = prop ~ 0 + bin,
    family = tweedie(),
    time_column = "year",
    times = min(long$year):max(long$year),
    variable_column = "bin",
    time_term = mods[[which_best]],
    control = tinyVASTcontrol( trace = 0 )
  )

  #
  which_rows = which( newdata$year == drop_years[z] )
  newdata[which_rows,'p_cv'] = predict( fit_z, newdata = newdata[which_rows,] )
}

library(viridisLite)

if( (names(mods)[which_best] == "cohort_year") & (dataset == "ages") ){
  library(igraph)
  parlist = fit[[which_best]]$obj$env$parList()
  DF = data.frame( from = c("(a,t)","(a,t)"),
                   to = c("(a+1,t+1)","(a,t+1)"),
                   label = parlist$nu_z[1:2] )
  DF$label = round(DF$label, digits=3)
  variables = unique(c(DF$from, DF$to))

  # Create and plotgraph
  pg <- graph_from_data_frame( d = DF,
                               directed = TRUE,
                               vertices = data.frame(variables) )

  coords = layout_(pg, with_sugiyama())
  coords$layout = cbind( x = c(0,1,1), y = c(0,1,0) )

  png( file = file.path(date_dir,"abundance_at_age_results.png"), width = 4, height = 6, res = 200, units = "in")
    par( mfrow=c(3,1), mar=c(0,0,2,1), mgp = c(1.5,0.25,0), tck=-0.02, oma = c(3,3,0,0) )

    plot( pg, layout = coords, vertex.color = "white", vertex.frame.color = NA,
          vertex.label = variables, vertex.label.cex = 1.5, edge.label.cex = 1.5 )
    title("Estimated graph")

    p_tc = make_matrix( list('bin' = unique(full$bin), 'year' = unique(full$year)) )
    p_tc[] = full$prop
    image( x = as.numeric(unique(full$year)),
           y = unique(full$bin),
           z = log(t(p_tc)),
           col = viridisLite::viridis(12),
           xlab = "", ylab = "", xaxt = "n" )
    title("Observed proportional abundance (log-scale)")

    p_tc[] = full$phat
    image( x = as.numeric(unique(full$year)),
           y = unique(full$bin),
           z = log(t(p_tc)),
           col = viridisLite::viridis(12),
           xlab = "", ylab = "" )
    title("Predicted proportional abundance (log-scale)")
    mtext( side = c(1,2), outer = TRUE, text = c("Year (t)", "Age (a)"), line = 1.75 )
  dev.off()
}


library(ggplot2)
full = merge( long, newdata, all = TRUE )
full$bin = as.numeric(gsub( as.character(full$bin), pattern = "bin_", replace = "" ))
ggplot( full ) +
  geom_point( aes(x = bin, y = prop) ) +
  facet_wrap( vars(year) ) +
  geom_line( aes(x = bin, y = phat, group = year) ) +
  geom_line( aes(x = bin, y = p_cv, group = year), col = "red" ) +
  scale_x_continuous( ) +
  theme(panel.border = element_rect(color = "black", fill = NA),
       panel.background = element_blank(),
       panel.grid = element_blank(),
       panel.spacing.x = unit(0,"line"))
ggsave( file.path(date_dir,paste0(dataset,"_results.png")), width = 6, height = 6)

