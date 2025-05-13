
#' @title Classify variables path
#'
#' @description \code{classify_variables} is copied from \code{sem:::classifyVariables}
#'
#' @details
#' Copied from package `sem` under licence GPL (>= 2) with permission from John Fox
#'
#' @param model SEM model
#'
#' @return Tagged-list defining exogenous and endogenous variables
#' @export
classify_variables <-
function( model ){

    variables <- logical(0)
    for (paths in model[, 1]) {
        vars <- gsub(pattern=" ", replacement="", x=paths)
        vars <- sub("-*>", "->", sub("<-*", "<-", vars))
        if (grepl("<->", vars)) {
            vars <- strsplit(vars, "<->")[[1]]
            if (is.na(variables[vars[1]]))
                variables[vars[1]] <- FALSE
            if (is.na(variables[vars[2]]))
                variables[vars[2]] <- FALSE
        }
        else if (grepl("->", vars)) {
            vars <- strsplit(vars, "->")[[1]]
            if (is.na(variables[vars[1]]))
                variables[vars[1]] <- FALSE
            variables[vars[2]] <- TRUE
        }
        else if (grepl("<-", vars)) {
            vars <- strsplit(vars, "<-")[[1]]
            if (is.na(variables[vars[2]]))
                variables[vars[2]] <- FALSE
            variables[vars[1]] <- TRUE
        }
        else stop("incorrectly specified model")
    }
    list(endogenous = names(variables[variables]), exogenous = names(variables[!variables]))
}


#' @title Parse path
#'
#' @description \code{parse_path} is copied from \code{sem::parse.path}
#'
#' @details
#' Copied from package `sem` under licence GPL (>= 2) with permission from John Fox
#'
#' @return Tagged-list defining variables and direction for a specified path coefficient
#'
#' @param path text to parse
#' @export
parse_path <-
function( path ){
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



#' @title Make a RAM (Reticular Action Model)
#'
#' @description \code{make_dsem_ram} converts SEM arrow notation to \code{ram} describing SEM parameters
#'
#' @inheritParams dsem
#' @param times A character vector listing the set of times in order
#' @param variables A character vector listing the set of variables
#' @param covs A character vector listing variables for which to estimate a standard deviation
#' @param quiet Boolean indicating whether to print messages to terminal
#' @param remove_na Boolean indicating whether to remove NA values from RAM (default) or not.
#'            \code{remove_NA=FALSE} might be useful for exploration and diagnostics for
#'            advanced users
#'
#' @details
#' \strong{RAM specification using arrow-and-lag notation}
#'
#'
#' @export
read_model <-
function( sem,
          times,
          variables,
          covs = NULL,
          quiet = FALSE ){

  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")

  ####### Error checks
  if( !is.numeric(times) ) stop("`times` must be numeric in `make_dsem_ram`")

  ####### Define local functions
  # helper function
  match_row = function( df, x ) which( df[1]==x[1] & df[2]==x[2] )
  #
  add.variances <- function() {
      variables <- need.variance()
      nvars <- length(variables)
      if (nvars == 0)
          return(model)
      message("NOTE: adding ", nvars, " variances to the model")
      paths <- character(nvars)
      par.names <- character(nvars)
      for (i in 1:nvars) {
          paths[i] <- paste(variables[i], "<->", variables[i])
          par.names[i] <- paste("V[", variables[i], "]", sep = "")
      }
      model.2 <- cbind(
        'path' = c(model[, 1], paths),
        'lag' = c(model[,2], rep(0,nvars)),
        'name' = c(model[, 3], par.names),
        'start' = c(model[, 4], rep(NA, length(paths))) )
      model.2
  }
  need.variance <- function() {
      all.vars <- classify_variables(model)
      exo.vars <- all.vars$exogenous
      end.vars <- all.vars$endogenous
      variables <- logical(0)
      for (i in seq_len(nrow(model))) {
          paths = model[i,1]
          lag = model[i,2]
          vars <- gsub(pattern=" ", replacement="", x=paths)
          vars <- sub("-*>", "->", sub("<-*", "<-", vars))
          vars <- sub("<->|<-", "->", vars)
          vars <- strsplit(vars, "->")[[1]]
          if ((vars[1] != vars[2]) | (lag != 0)) {
              for (a.variable in vars) {
                if (is.na(variables[a.variable]))
                  variables[a.variable] <- TRUE
              }
          }
          else {
              variables[vars[1]] <- FALSE
          }
      }
      if (!exog.variances && length(exo.vars) > 0)
          variables[exo.vars] <- FALSE
      if (!endog.variances && length(end.vars) > 0)
          variables[end.vars] <- FALSE
      names(variables)[variables]
  }

  ####### Step 2 -- Make RAM
  # convert to data frame
  model = scan( text = sem,
                what = list(path = "", time_lag = 1, space_lag = 1, par = "", start = 1, dump = ""),
                sep = ",",
                strip.white = TRUE,
                comment.char = "#",
                fill = TRUE,
                quiet = quiet)
  model$path <- gsub("\\t", " ", model$path)
  model$par[model$par == ""] <- NA
  model <- data.frame( "path"=model$path, "time_lag"=model$time_lag, "space_lag"=model$space_lag,
                  "name"=model$par, "start"=model$start)

  # Adding a SD automatically
  if( !is.null(covs) ){
    for (cov in covs) {
      vars <- strsplit(cov, "[ ,]+")[[1]]
      nvar <- length(vars)
      for (i in 1:nvar) {
      for (j in i:nvar) {
        p1 = paste(vars[i], "<->", vars[j])
        p2 = if (i==j) paste("V[", vars[i], "]", sep = "") else paste("C[",vars[i], ",", vars[j], "]", sep = "")
        p3 = NA
        row <- c(p1, 0, 0, p2, p3)
        if( any((row[1]==model[,1]) & (row[2]==model[,2])) ){
          next
        }else{
          model <- rbind(model, row, deparse.level = 0)
        }
      }}
    }
  }

  exog.variances = endog.variances = TRUE
  model = add.variances()

  # Add parameter column
  par.names = model[, 4]
  pars = na.omit(unique(par.names))
  par.nos = apply(outer(pars, par.names, "=="), 2, which)
  par.nos = unlist(sapply( par.nos, FUN=\(x) ifelse(length(x)==0, 0, x) ))
  model = cbind( model, "parameter"=par.nos )

  # Add incidence to model
  model = cbind( model, first=NA, second=NA, direction=NA )
  for( i in seq_len(nrow(model)) ){
    path = parse_path(model[i,1])
    model[i,c('first','second','direction')] = unlist( path[c('first','second','direction')] )
  }

  ####### Step 2 -- Make RAM

  # Deal with fixed values
  #if( max(model$parameter) != length(beta_p) ) stop("Check beta_p")

  return(model)
}

make_matrices <-
function( beta_p,
          model,
          times,
          variables,
          spatial_graph,
          use_area = FALSE ){

  #
  if( is(spatial_graph,"sfc") ){
    # Adjacency
    Method = "grid"
    if( use_area==TRUE ){
      area_s = st_area( spatial_graph )
      units(area_s) = NULL
    }else{
      area_s = rep(1,length(spatial_graph))
    }
    st_rook = function(m, ...) st_relate(m, m, pattern="F***1****", ... )
    # M1_ss
    A_ss = st_rook( spatial_graph, sparse=TRUE )
    A_ss = as(A_ss, "sparseMatrix")
    A_ss = Diagonal(nrow(A_ss), x=1/rowSums(A_ss)) %*% A_ss   # rowsums = 1
    diag(A_ss) = -1                                           # rowsums = 0 to match G from SPDE method
    M1_ss = Diagonal(nrow(A_ss), x=1/area_s) %*% A_ss   # scale-invariance
    # M0
    M0_ss = Diagonal(nrow(A_ss))
    # Diagonal matrix for SPDE
    D_ss = Diagonal(nrow(A_ss))
  }else if( is(spatial_graph,"inla.mesh") ){
    spde = fm_fem( mesh )
    # Original parameterization
    #M0_ss = spde$c0
    #M1_ss = -1 * spde$g1
    # Transform parameterization
    #
    M0_ss = Diagonal( mesh$n )
    #
    diag_c0 = diag(spde$c0)
    M1_ss = -1 * Diagonal( n=nrow(spde$c0), x=1/diag_c0 ) %*% spde$g1
    # Only used for two-headed arrows
    D_ss = Diagonal( n=nrow(spde$c0), x=1/sqrt(diag_c0) )
  }else{
    stop("Invalid spatial_graph")
  }

  if(missing(beta_p)){
    model_unique = model[match(unique(model$parameter),model$parameter),]
    beta_p =  model_unique$start
  }

  # Loop through paths
  G2_kk = P2_kk = G_kk = P_kk = drop0(sparseMatrix( i=1, j=1, x=0, dims=rep(nrow(M0_ss)*length(variables)*length(times),2) ))   # Make with a zero
  for( i in seq_len(nrow(model)) ){
    time_lag = as.numeric(model[i,2])
    space_lag = as.numeric(model[i,3])
    L_tt = sparseMatrix( i = seq(time_lag+1,length(times)),
                         j = seq(1,length(times)-time_lag),
                         x = 1,
                         dims = rep(length(times),2) )

    P_jj = sparseMatrix( i = match(model[i,'second'],variables),
                         j = match(model[i,'first'],variables),
                         x = 1,
                         dims = rep(length(variables),2) )

    # Grid
    if(space_lag==0){
      M_ss = M0_ss
    }else{
      M_ss = M1_ss
    }

    # Assemble
    if(abs(as.numeric(model[i,'direction']))==1){
      tmp_kk = kronecker(kronecker(P_jj, L_tt), M_ss)
      P_kk = P_kk + beta_p[model$parameter[i]] * tmp_kk
    }else{
      tmp_kk = kronecker(kronecker(P_jj, L_tt), D_ss)
      G_kk = G_kk + beta_p[model$parameter[i]] * tmp_kk
    }
  }

  # Diagonal component
  I_kk = Diagonal(nrow(P_kk))

  # Assemble
  IminusP_kk = I_kk - P_kk
  invV_kk = AD(G_kk)
  invV_kk@x = 1 / G_kk@x^2
  Q_kk = t(IminusP_kk) %*% invV_kk %*% IminusP_kk

  out = list(
    "P_kk" = P_kk,
    "G_kk" = G_kk,
    "invV_kk" = invV_kk,
    "IminusP_kk" = IminusP_kk,
    "Q_kk" = Q_kk # ,
  )
  return(out)
}
