
Q_network <-
function( log_theta,
          tree,
          estimate_ou = TRUE,
          bm_distance_ratio = 1000 ){
  # As bm_distance_ratio -> Inf, the BM approaches rank-deficient ICAR model

  # Necessary in packages
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")

  # Locals
  parent_s = tree$edge[,1]  # max(table(tree$edge[,1])) = 2
  child_s = tree$edge[,2]   # max(table(tree$edge[,2])) = 1
  dist_s = tree$edge.length
  n_s = max( tree$edge )
  n_edge = length(parent_s)
  if(any(table(child_s))>1) stop("Not a tree or stream")

  # Assemble pieces
  #  NOTE:  can't use M$x and assume it's in the same order as vectors i and j in sparseMatrix(i,j)
  if( isTRUE(estimate_ou) ){
    theta = exp( log_theta );
    margvar_hat = 1 / ( 2 * theta ) #
    v1_s = exp(-theta * dist_s) / (1 - exp(-2 * theta * dist_s))
    v2_s = exp(-2 * theta * dist_s) / (1 - exp(-2 * theta * dist_s))

    # Make components
    P_ss = AD(sparseMatrix( i = child_s, j = parent_s, x = 0, dims = c(n_s,n_s) ))
    P_ss[cbind(child_s, parent_s)] = v1_s

    #
    D0_ss = AD(sparseMatrix( i = child_s, j = child_s, x = 0, dims = c(n_s,n_s) ))
    D0_ss[cbind(child_s, child_s)] = v2_s
    D_ss = D0_ss + AD(Diagonal( n_s ))
  }else{
    margvar_hat = 1
    # As theta -> 0:
    #  P_ss = 1 / ( dist * 2*theta)
    # D_ss = 1 / ( dist * 2*theta )  OR 1 for the root
    P_ss = AD(sparseMatrix( i = child_s,
                            j = parent_s,
                            x = 1 / dist_s,
                            dims = c(n_s,n_s) ))
    D0_ss = AD(sparseMatrix( i = child_s,
                            j = child_s,
                            x = 1 / dist_s, # dist_s,
                            dims=c(n_s,n_s) ))
    #D1_ss = AD(Diagonal( n_s ))

    # if isTRUE(drop_bm_root), then bm_distance_ratio doesn't matter
    roots = setdiff( seq_len(n_s), child_s )
    if(isTRUE(drop_bm_root)){
      D1_ss = AD(sparseMatrix( i = roots, j = roots, x = 0, dims = c(n_s,n_s) ))
    }else{
      D1_ss = AD(sparseMatrix( i = roots, j = roots, x = min(dist_s) / bm_distance_ratio, dims = c(n_s,n_s) ))
    }
    D_ss = D0_ss + D1_ss
    v1_s = v2_s = NULL
  }

  # Rescale
  D_ss = D_ss / margvar_hat
  P_ss = P_ss / margvar_hat

  # Transformations
  #invD_ss = AD(sparseMatrix( i = seq_len(n_s), j = seq_len(n_s), x = 0, dims=c(n_s,n_s) ))
  #invD_ss@x = 1 / D_ss@x

  # Assemble
  #DminusP_ss = D_ss - P_ss
  #Q = (t(DminusP_ss) %*% invD_ss %*% DminusP_ss) / margvar_hat

  # Return stuff
  network = list(
    #Q = Q,
    #margvar_hat = margvar_hat
    P_ss = P_ss,
    D_ss = D_ss,
    v1_s = v1_s,
    v2_s = v2_s
  )
  return(network)
}

parse_path <-
function( path ){
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  path.1 <- gsub("-", "", gsub(" ", "", path))
  direction <- if(regexpr("<>", path.1) > 0){
    2
  }else if(regexpr("<", path.1) > 0){
    -1
  }else if(regexpr(">", path.1) > 0){
    1
  }else{
    stop(paste("ill-formed path:", path))
  }
  path.1 <- strsplit(path.1, "[<>]")[[1]]
  out = list(first = path.1[1], second = path.1[length(path.1)], direction = direction)
  return(out)
}

# Function that converts SEM model to a RAM, see `?sem` for more context
build_ram = function( model, vars ){
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")

  #
  vars = sapply( vars, FUN=function(char){gsub("-", "", gsub(" ", "", char))} )
  n.paths = nrow(model)
  par.names = model[, 2]
  startvalues = model[,3]

  # EXCERPT FROM `getAnywhere("sem.semmod")`
  heads = from = to = rep(0, n.paths)
  for (p in 1:n.paths) {
    #path = sem:::parse.path(model[p, 1])
    path = parse_path(model[p, 1])
    heads[p] = abs(path$direction)
    to[p] = path$second
    from[p] = path$first
    if (path$direction == -1) {
      to[p] = path$first
      from[p] = path$second
    }
  }
  missing_vars = setdiff( c(from,to), vars )
  if( length(missing_vars) > 0 ) stop( "Check `build_ram`:", paste0(missing_vars,sep=", ") )

  ram = data.frame(matrix(0, nrow=p, ncol=5))
  pars = na.omit(unique(par.names))
  ram[, 1] = heads
  ram[, 2] = vars[apply(outer(vars, to, "=="), 2, which)]
  ram[, 3] = vars[apply(outer(vars, from, "=="), 2, which)]
  par.nos = apply(outer(pars, par.names, "=="), 2, which)
  if(length(par.nos) > 0){
    ram[, 4] = unlist(lapply(par.nos, function(x) if (length(x)==0){0}else{x}))
  }
  ram[, 5] = startvalues
  colnames(ram) = c("direction", "second", "first", "parameter", "start")
  return(ram)
}


get_nll <-
function( parlist ){

  model <- as.data.frame(model)
  model[,'parameter'] = as.integer(model[,'parameter'])
  if( any((model$direction==2) & (model$first!=model$second)) ){
    stop("Exogenous covariance not allowed")
  }

  # Combine fixed, estimated, and mapped parameters into vector
  beta_i = rep(0, nrow(model))
  off = which(model[,'parameter'] == 0)
  if( length(off) > 0 ){
    beta_i[off] = as.numeric(model[off,'start'])
  }
  not_off = which(model[,'parameter'] > 0)
  if( length(not_off) > 0 ){
    beta_i[not_off] = parlist$beta_p[model[not_off,'parameter']]
  }

  #Q_kk = drop0(sparseMatrix( i=1, j=1, x=0, dims=c(1,1) ))
  variables = colnames(data)
  n_s = n_nodes + n_tips
  #n_s = max( tree$edge )

  P_kk = invV_kk = drop0(sparseMatrix( i=1, j=1, x=0, dims=rep(length(variables)*n_s,2) ))   # Make with a zero
  I_ss = Diagonal( n = n_s )
  #I_jj = Diagonal( n = length(variables) )
  for( i in seq_len(nrow(model)) ){

    # Interaction matrix
    P_jj = sparseMatrix( i = match(model[i,'second'],variables),
                         j = match(model[i,'first'],variables),
                         x = 1,
                         dims = rep(length(variables),2) )

    # Should it be first or second ???
    Ind_jj = sparseMatrix( i = match(model[i,'second'],variables),
                           j = match(model[i,'second'],variables),
                           x = 1,
                           dims = rep(length(variables),2) )
    which_theta = match(model[i,'second'], variables)
    #theta = exp(parlist$ln_theta[which_theta])
    network = Q_network( log_theta = parlist$ln_theta[which_theta],
                         estimate_ou = ou_j[which_theta],
                         tree = tree )

    # Transformations
    D_ss = network$D_ss
    invD_ss = AD(sparseMatrix( i = seq_len(n_s), j = seq_len(n_s), x = 0, dims=c(n_s,n_s) ))
    invD_ss@x = 1 / D_ss@x

    # Assemble P
    if(abs(as.numeric(model[i,'direction']))==1){
      if( assemble_version %in% c(1,2) ){
        tmp_kk = kronecker(P_jj, I_ss)
      }else if( assemble_version == 3 ){
        Pstar_ss = invD_ss %*% network$P_ss
        #tmp_kk = kronecker(P_jj, Pstar_ss)
        IminusPstar_ss = I_ss - Pstar_ss
        #Istar_ss = invD_ss %*% I_ss
        tmp_kk = kronecker(P_jj, IminusPstar_ss)
      }
      P_kk = P_kk + beta_i[i] * tmp_kk
    }

    # Assemble D
    if(abs(as.numeric(model[i,'direction']))==2){
      # Previous
      if( assemble_version %in% c(1,2) ){
        if( assemble_version == 1 ){
          # Assemble version-1
          DminusP_ss = D_ss - network$P_ss
          Qhat_ss = (t(DminusP_ss) %*% invD_ss %*% DminusP_ss) # / network$margvar_hat
        }else if( assemble_version == 2 ){
          # Assemble version-2 (equivalent to version-1)
          Pstar_ss = invD_ss %*% network$P_ss
          IminusPstar_ss = Diagonal( n=nrow(invD_ss) ) - Pstar_ss
          Qhat_ss = t(IminusPstar_ss) %*% D_ss %*% IminusPstar_ss
        }
        # Combine elements for version-1 or version-2
        tmp_kk = kronecker(P_jj, Qhat_ss)
        invV_kk = invV_kk + (1 / beta_i[i]^2) * tmp_kk # AD(tmp_kk)
      }else if( assemble_version == 3 ){
        # Assemble version-3
        Pstar_ss = invD_ss %*% network$P_ss
        tmp1_kk = kronecker(Ind_jj, Pstar_ss)
        P_kk = P_kk + tmp1_kk
        tmp2_kk = (1 / beta_i[i]^2) * kronecker(Ind_jj, D_ss)
        invV_kk = invV_kk + tmp2_kk
      }
    }
  }
  #Q_kk = invV_kk

  # Assemble
  I_kk = Diagonal(nrow(P_kk))
  IminusP_kk = AD(I_kk - P_kk)
  Qprime_kk = Q_kk = t(IminusP_kk) %*% invV_kk %*% IminusP_kk

  #P_jj = sparseMatrix( i = 1,
  #                     j = 1,
  #                     x = 1,
  #                     dims = rep(1,2) )
  #Q_jj = Q_jj + P_jj / exp(parlist$ln_sigma)^2
  #Q_kk = kronecker( Q_jj, Qhat_ss )

  #
  yprime_k = y_k = as.vector( parlist$y_sj )
  muprime_k = mu_k = rep( parlist$xbar, each = nrow(parlist$y_sj) )

  #
  if( length(which_drop)>0 & isTRUE(drop_bm_root) ){
    Qprime_kk = Qprime_kk[-which_drop,-which_drop]
    yprime_k = yprime_k[-which_drop]
    muprime_k = muprime_k[-which_drop]
  }

  loglik = dgmrf( yprime_k,
                  mu = muprime_k,
                  Q = Qprime_kk,
                  #Q = Qhat_ss / exp(2 * parlist$ln_sigma),
                  log = TRUE )

  REPORT( invV_kk )
  REPORT( P_kk )
  REPORT( Qprime_kk )
  REPORT( loglik )
  return( -1* loglik )
}

get_inner_H <-
function( obj,
          opt_par = obj$env$last.par.best ){

  par_null = obj$env$parList( x = opt_par )
  obj_null = MakeADFun( func = get_nll,
                   random = "y_sj",
                   parameters = par_null )
  H_null = obj_null$env$spHess( random=TRUE )
  return( H_null )
}
