
library(RTMB)
library(Matrix)
library(mvtnorm)
library(ape)
library(sem)

root_dir = R'(C:\Users\James.Thorson\Desktop\Git\GGMM)'
data_dir = file.path( root_dir, "data" )

# Compress VertTree for github
if( FALSE ){
  local_dir = R'(C:\Users\James.Thorson\Desktop\Work files\AFSC\2025-04 -- Mammal tree)'
  tree = read.nexus( file.path( local_dir, "output.nex" ) )
  tree = tree$tree_1998
  write.tree( tree,
              file = file.path(data_dir,"VertTree_mammals.tre") )
}

#
Date = Sys.Date()
run_dir = file.path(root_dir, Date )
  dir.create(run_dir)

# Functions
source( file.path(root_dir, "trait_functions.R") )

dataset = c("mammals", "rhino", "small_rhino", "line", "tree", "bivar_tree")[1]

if( dataset == "mammals" ){
  tree = read.tree( file.path( data_dir, "VertTree_mammals.tre" ) )
  max_edge = max(tree$edge.length)
  tree$edge.length = tree$edge.length / max_edge
  s_root = Ntip(tree) + 1
  n_nodes = Nnode(tree)
  n_tips = Ntip(tree)

  traits = read.table( file.path(data_dir,"PanTHERIA_1-0_WR05_Aug2008.txt"),
              row.names = NULL )

  tree$tip.label = tolower( tree$tip.label )
  traits$traits_binom = paste0( tolower(traits$MSW05_Species), "_", traits$MSW05_Binomial )

  trait_to_tip = match( traits$traits_binom, tree$tip.label )
  traits = traits[ which(!is.na(trait_to_tip)), ]

  add_NA = function(vec) ifelse( as.numeric(vec) == -999, NA, as.numeric(vec) )
  data = data.frame(
    ln_metabolism = log(add_NA(traits[,'X18.1_BasalMetRate_mLO2hr'])),
    ln_range = log(add_NA(traits[,'X22.1_HomeRange_km2'])),
    ln_size = log(add_NA(traits[,'X5.1_AdultBodyMass_g']))
  )
  rownames(data) = traits$traits_binom

  #
  sem = "
    ln_size -> ln_metabolism, b1
    ln_size -> ln_range, b2
  "
  ou_j = c(TRUE, TRUE, TRUE)

  #
  plm = phylolm::phylolm(
    ln_metabolism ~ ln_size,
    data = data,
    phy = tree,
    model = "OUrandomRoot"
  )

  # CI relationship
  sem_dsep = "
    ln_metabolism -> ln_range, target
    ln_size -> ln_range, cond1
  "

  # Table of data
  has_data = ifelse( is.na(data), 0, 1 )
  count = t(has_data) %*% has_data
  write.csv( count, file = file.path(run_dir,"traits_data_available.csv") )
}
if( dataset == "rhino" ){
  # Load data set
  data(rhino, rhino_tree, package="phylopath")

  # Run phylosem
  sem = "
    DD -> RS, p1
    BM -> LS, p2
    BM -> NL, p3
    NL -> DD, p4
  "
  tree = rhino_tree
  n_nodes = Nnode(tree)
  n_tips = Ntip(tree)
  s_root = Ntip(tree) + 1

  # Settings
  data = rhino[,c("BM","NL","DD","RS","LS")]
  #data = rhino[,c("BM","NL"),drop = FALSE]
  ou_j = rep( FALSE, ncol(data) )
}
if( dataset == "small_rhino" ){
  # Load data set
  data(rhino, rhino_tree, package="phylopath")

  # Run phylosem
  sem = "
    DD -> RS, p1
  "
  # Settings
  data = rhino[1:4,c("DD","RS")]
  #data = rhino[,c("BM","NL"),drop = FALSE]
  ou_j = c( FALSE, FALSE )

  tree = keep.tip( rhino_tree, tip = rownames(data) )
  n_nodes = Nnode(tree)
  n_tips = Ntip(tree)
  s_root = Ntip(tree) + 1
}
if( dataset == "line" ){
  n_tips = 100
  n_nodes = 0
  x = sort(runif(n_tips, min=0, max=10))
  #sigma = 0.8
  #theta = 1.2
  sigma = 1
  theta = 0.001
  (margvar = sigma^2 / (2*theta)) # margvar
  Dist_ss = outer( x, x, FUN = \(a,b){abs(a-b)} )
  V_ss = margvar * exp(-theta * Dist_ss ) #
  Q_ss = solve(V_ss)
  #y = RTMB:::rgmrf0( n=1, Q = Q_ss)[,1]

  # Sample
  y = ytrue = rmvnorm( n = 1, sigma = V_ss )[1,]
  which_na = sample( seq_along(y), replace = FALSE, size = floor(length(y)/2) )
  y[which_na] = NA

  tree = list(
    edge = cbind( 1:(length(x)-1), 2:length(x) ),
    edge.length = diff(x),
    tip.label = paste0("t",seq_len(n_tips))
  )
  s_root = setdiff( seq_along(x), tree$edge[,2] )
  # c( mean(y^2), margvar )

  #
  sem = ""
  ou_j = TRUE
  data = data.frame(y=y)[seq_len(n_tips),,drop=FALSE]
  rownames(data) = tree$tip.label
}
if( dataset == "tree"){
  #x = sort(runif(100, min=0, max=10))
  n_tips = 4
  tree = rtree(n_tips)
  n_nodes = Nnode(tree)
  sigma = 0.8
  theta = 1
  (margvar = sigma^2 / (2*theta))
  Dist_ss = dist.nodes(tree)
  V_ss = margvar * exp(-theta * Dist_ss ) #
  Q_ss = solve(V_ss)
  #y = RTMB:::rgmrf0( n=1, Q = Q_ss)[,1]
  s_root = Ntip(tree) + 1

  # Sample
  y = ytrue = rmvnorm( n = 1, sigma = V_ss )[1,]
  y = y[seq_len(n_tips)]
  which_na = sample( seq_along(y), replace = FALSE, size = floor(length(y)/2) )
  y[which_na] = NA

  #
  sem = ""
  ou_j = FALSE
  data = data.frame(y=y)
  rownames(data) = tree$tip.label
}
if( dataset == "bivar_tree"){
  #x = sort(runif(100, min=0, max=10))
  n_tips = 4
  tree = rtree(n_tips)
  n_nodes = Nnode(tree)
  sigma = 0.8
  theta = 1
  (margvar = sigma^2 / (2*theta))
  rhoxy = 0.5
  Dist_ss = dist.nodes(tree)
  V_ss = margvar * exp(-theta * Dist_ss ) #
  Q_ss = solve(V_ss)
  #y = RTMB:::rgmrf0( n=1, Q = Q_ss)[,1]
  s_root = Ntip(tree) + 1

  # Sample
  x = xtrue = rmvnorm( n = 1, sigma = V_ss )[1,]
  y = ytrue = rhoxy*x + rmvnorm( n = 1, sigma = V_ss )[1,]

  x = x[seq_len(n_tips)]
  y = y[seq_len(n_tips)]
  #which_na = sample( seq_along(y), replace = FALSE, size = floor(length(y)/2) )
  #y[which_na] = NA

  #
  sem = "x -> y, rhoxy"
  ou_j = c( TRUE, TRUE )
  data = data.frame(x=x, y=y)
  rownames(data) = tree$tip.label
}


#Q_network( log_theta = 0, tree = tree, estimate_ou = TRUE)$Q
#Q_network( log_theta = 0, tree = tree, estimate_ou = FALSE)$Q


#############
# Fit in RTMB
#############

# Whether to drop root for ISAR
drop_bm_root = TRUE
which_drop = ((n_tips+n_nodes) * (seq_along(ou_j)-1) + (n_tips+1))[(ou_j==FALSE)]
assemble_version = 3

#
if( is.null(tree$node.label) & (n_nodes > 0) ){
  tree$node.label = paste0("node_",seq_len(n_nodes))
}

#
SEM_model = specifyModel( text=sem,
                          exog.variances=TRUE,
                          endog.variances=TRUE,
                          covs=colnames(data),
                          quiet=TRUE )
model = build_ram( model = SEM_model,
           vars = colnames(data) )

parlist = list(
  y_sj = as.matrix(data[match(c(tree$tip.label,tree$node.label),rownames(data)),,drop=FALSE]),
  beta_p = rep(1,max(model$parameter)),
  ln_theta = rep(0,ncol(data)),
  xbar = rep(0,ncol(data))
)
map = list(
  #ln_theta = factor( rep(1,ncol(data)) ),
  y_sj = ifelse( is.na(parlist$y_sj), seq_len(prod(dim(parlist$y_sj))), NA ),
  ln_theta = factor(ifelse(ou_j, seq_len(ncol(data)), NA))
)
if(isTRUE(drop_bm_root)) map$y_sj[which_drop] = NA
map$y_sj = factor(map$y_sj)

parlist$y_sj = ifelse( is.na(parlist$y_sj), 0, parlist$y_sj )

#method = "GMRF"
obj = MakeADFun( func = get_nll,
                  random = "y_sj",
                  map = map,
                  parameters = parlist,
                  silent = TRUE )
opt = nlminb( obj$par, obj$fn, obj$gr,
              control = list() )
rep = obj$report()
sdrep = sdreport(obj)
#Q1 = get_inner_H( obj1, opt_par = opt1$par )
#Q2 = rep1$Qhat_ss
#range( Q1 - Q2 ) # Should be 0

# Compare with phylosem
psem = phylosem::phylosem(
  data = data,
  sem = sem,
  estimate_ou = all(ou_j),
  tree = tree
)

if( FALSE ){
  parlist$beta_p <- psem$opt$par
}

if( dataset %in% c("tree") ){
  plm = phylolm::phylolm(
    data = data,
    formula = y ~ 1,
    phy = tree,
    model = "BM"
  )
}
if( dataset %in% c("small_rhino") ){
  plm = phylolm::phylolm(
    data = data,
    formula = RS ~ DD,
    phy = tree,
    model = "BM"
  )
}
if( dataset %in% c("bivar_tree") ){
  plm = phylolm::phylolm(
    data = data,
    formula = y ~ x,
    phy = tree,
    model = "OUrandomRoot"
  )
}

# Compare SDs
c( opt$par['beta_p'], psem$opt$par['beta_z'] )

#
if( dataset == "mammals" ){
  theta = exp( obj$env$parList()$ln_theta ) / max_edge
  ln_theta_se = as.list(sdrep, what = "Std. Error")$ln_theta
  theta_lo = theta * exp(-1.96 * ln_theta_se)
  theta_hi = theta * exp(1.96 * ln_theta_se)
  max_depth = max(node.depth.edgelength( tree ))
  xset = seq(0, max_depth * max_edge, length = 1000)
  FUN = function(p) exp(-p*xset)
  yset = sapply( theta, FUN = FUN )

  library(viridisLite)
  library(igraph)

  DF = data.frame( from=model$first,
                   to=model$second,
                   label=obj$env$parList()$beta_p )
  DF = subset( DF, DF$from != DF$to )
  DF$label = round(DF$label, digits=3)
  f = function(char_vec) gsub(x=char_vec, pattern="ln_", replace = "")
  variables = f(colnames(data))
  #variables = paste0( variables, "(t)" )
  DF$from = f(DF$from)
  DF$to = f(DF$to)

  #DF2 = data.frame(
  #  from = paste0( f(DF$from), "(t)" )
  #)

  # Create and plotgraph
  pg <- graph_from_data_frame( d = DF,
                               directed = TRUE,
                               vertices = data.frame(variables) )

  png( file = file.path(run_dir,"traits.png"), width=4, height = 5, res = 200, units = "in")
    par( mfrow = c(2,1), mar=c(0,2.5,1,0.5), mgp = c(1.5,0.25,0), tck = -0.02, oma = c(3,0,0,0) )

    # Panel-1
    coords = layout_(pg, with_sugiyama())
    plot( pg, layout = coords, vertex.color = "white", vertex.frame.color = NA,
          vertex.label = variables )
    title( "Trait interactions" )

    # Panel-2
    matplot(
      x = xset,
      y = yset,
      type = "l",
      lwd = 2,
      lty = "solid",
      col = viridisLite::viridis(3),
      xlab = "",
      ylab = "O-U correlation",
      xaxs = "i",
      yaxs = "i"
    )
    title( "Correlation over time" )
    mtext( side=1, text = "Evolutionary time (million years)", line=1.5 )
    legend( "topright", bty = "n", fill = viridis(3),
            legend = variables )
    for( i in seq_len(ncol(data)) ){
      polygon( x = c(xset, rev(xset)),
               y = c(FUN(theta_lo[i]), rev(FUN(theta_hi[i]))),
               col = viridis(3, alpha = 0.2)[i],
               border = NA )
    }
  dev.off()

}
