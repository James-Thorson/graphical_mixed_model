
# FROM:  SDSEM_2025-04-18.R

root_dir = R'(C:\Users\James.Thorson\Desktop\Git\GGMM)'

#
Date = Sys.Date()
run_dir = file.path(root_dir, Date )
  dir.create(run_dir)

library(Matrix)
library(sf)
library(ggplot2)
library(RTMB)
#library(fmesher)
library(viridisLite)
library(igraph)
source( file.path(root_dir,"make_sdsem_ram.R") )

ind = \(n,i){ x=rep(0,n); x[i]=1; return(x)}

# Settings
size = 0.1
extent = 1
n_t = 3
n_missing = 0

# make extrapolation grid
domain = st_sfc(st_polygon(list(extent*cbind(c(-1,1,1,-1,-1),c(-1,-1,1,1,-1)))))
domain = st_make_grid( domain, cellsize=c(size,size) )
domain_locs = st_coordinates(st_centroid(domain))

# make plotting grid
extrap = st_sfc(st_polygon(list(extent*cbind(c(-1,1,1,-1,-1),c(-1,-1,1,1,-1)))))
extrap = st_make_grid( extrap, cellsize=c(size,size) / 10 )
extrap_locs = st_coordinates(st_centroid(extrap))

#
mesh_locs = pracma::poisson2disk(n=length(domain)*2, a=extent*2, b=extent*2 )
mesh = fm_mesh_2d( mesh_locs-extent, refine=TRUE )

png( file=file.path(run_dir,"diffusion.png"), width=7, height=6, res=200, units="in" )
  par( mfrow=c(3,4), mar=c(0,0,0,0), mgp=c(2,0.5,0), tck=-0.02, oma=c(0,3,2,0) )
  for(row_index in 1:3){

    #
    method = "grid"
    pal = list( viridis(100), magma(100) )[[ifelse(method=="grid",1,2)]]

    # conditional stuff
    spatial_graph = domain
    which_mid = which.min( rowMeans(scale(domain_locs,center=TRUE)^2) )
    use_area = TRUE
    n_s = length(domain)
    A_gs = st_within( st_centroid(extrap), domain )
    A_gs = as(A_gs,"sparseMatrix")
    if(row_index == 1){
      sdsem = "
        var <-> var, 0, 0, sd_var, 1
        var -> var, 0, 1, space_cor, 0.1
        var -> var, 1, 0, time_cor, 0.8
      "
    }
    if(row_index==2){
      sdsem = "
        var <-> var, 0, 0, sd_var, 1
        var -> var, 0, 1, space_cor, 0.1
        var -> var, 1, 0, time_cor, 0.8
        var -> var, 1, 1, diffusion, -0.04
      "
    }
    if(row_index==3){
      sdsem = "
        var <-> var, 0, 0, sd_var, 1
        var -> var, 0, 1, space_cor, 0.1
        var -> var, 1, 0, time_cor, 0.8
        var -> var, 1, 1, diffusion, -0.08
      "
    }
    mid_area = 1

    # Build RAM
    model = read_model(
              sem = sdsem,
              times = seq_len(n_t),
              variables = c("var"),
              covs = NULL,
              quiet = FALSE )
    ram = make_matrices(
              model = model,
              times = seq_len(n_t),
              variables = c("var"),
              spatial_graph = spatial_graph,
              use_area = use_area )

    # Generate indicator
    logD1 = ind(n_s, which_mid) / mid_area
    logD_k = c( logD1, rep(0,nrow(ram$P_kk)-n_s) )
    # Reformat either one
    IminusP_kk = ram$IminusP_kk
    logD_st = matrix( solve(IminusP_kk, logD_k), ncol=n_t )
    logD_st = as.matrix(A_gs %*% logD_st)

    # Edit data
    colSums(logD_st)
    exp(diff(log(colSums(logD_st))))
    N_st = array( rpois(prod(dim(logD_st)), lambda=100*exp(logD_st)), dim=dim(logD_st) )
    Z_stc = N_st %o% 1
    dimnames(Z_stc) = list( seq_len(dim(Z_stc)[1]), seq_len(dim(Z_stc)[2]), c("var") )

    # Calculate diffusion ~propto~ MSD
    Mean_tz = matrix(NA, nrow=n_t, ncol=2)
    MSD_t = rep(NA, n_t)
    for(t in seq_len(n_t) ){
      Mean_tz[t,] = apply(extrap_locs, MARGIN=2, FUN=weighted.mean, w=logD_st[,t])
      MSD_t[t] = weighted.mean( rowSums( (extrap_locs-rep(1,nrow(extrap_locs))%o%Mean_tz[t,])^2), w=logD_st[,t] )
    }
    #diff( MSD_t )

    if(row_index == 1){
      DF = data.frame( from = rep("(s,t)", 2),
                       to = c("(s+1,t)","(s,t+1)"),
                       label = model$start[-1] )
    }
    if(row_index %in% 2:3){
      DF = data.frame( from = rep("(s,t)", 3),
                       to = c("(s+1,t)","(s,t+1)","(s+1,t+1)"),
                       label = model$start[-1] )
    }
    variables = rbind( "(s,t)" = c(0,0), "(s+1,t)" = c(1,0), "(s,t+1)" = c(0,1), "(s+1,t+1)" = c(1,1) )

    # Create and plotgraph
    pg <- graph_from_data_frame( d = DF,
                                 directed = TRUE,
                                 vertices = data.frame(rownames(variables)) )
    # Panel-1
    coords = layout_(pg, with_sugiyama())
    coords$layout = variables
    plot( pg, layout = coords, vertex.color = "white", vertex.frame.color = NA,
          vertex.label = rownames(variables),
          vertex.label.cex = 1.5, edge.label.cex = 1.5,
          xlim=1.4*c(-1,1), ylim=1.0*c(-1,1) )
    mtext(side=2, line=1, text=c("diffusing","intermediate","separable")[row_index] )
    #mtext(side=3, line=1.5, text="Graph" )

    # Plot
    for( col_index in seq_len(n_t) ){
      plot( st_sf(extrap, logD_st[,col_index]), border=NA,    # [,column_index]
            key.pos=NULL, reset=FALSE, pal=pal,
            breaks=seq(0,max(logD_st[,col_index]),length=101),
            main = "" )
      if( row_index==1 ) mtext(side=3, line=1.5, text=paste0("time ", col_index) )
      mtext( side=3, line = 0, paste0("MSD = ", round(MSD_t[col_index],3)) )
    }
    H = ram$Q_kk
    Title = paste0( "Nonzero=", length(H@x)," Sparsity=", round(length(H@x)/prod(dim(H)),3) )
    print(Title)
    H = drop0(H)
    Title = paste0( "Nonzero=", length(H@x)," Sparsity=", round(length(H@x)/prod(dim(H)),3) )
    print(Title)
  }
dev.off()
