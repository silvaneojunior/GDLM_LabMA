#' generic_smoother
#'
#' Generic smoother for all models.
#'
#' @param mt Matrix: A matrix containing the filtered mean of the latent variables at each time. Each row should represent one variable.
#' @param Ct Array: A 3D-array representing the filtered covariance matrix of the latent variables at each time. The third dimension should represent the time index.
#' @param at Matrix: A matrix containing the one-step-ahead mean of the latent variables at each time based upon the filtered mean. Each row should represent one variable.
#' @param Rt Array: A 3D-array representing the one-step-ahead covariance matrix of the latent variables at each time based upon the filtered covariance matrix. The third dimension should represent the time index.
#' @param G  Matrix: A matrix representing the transition matrix of the model.
#'
#' @return List: The smoothed mean (mts) and covariance (Cts) of the latent variables at each time. Their dimension follows, respectivelly, the dimensions of mt and Ct.
#' @export
#'
#' @examples
#' T=20
#'
#' mt=matrix(c(cumsum(rnorm(T)+1),rep(1,T)),2,T,byrow=TRUE)
#' Ct=array(diag(c(1,1)),c(2,2,T))
#' G=matrix(c(1,0,1,1),2,2)
#' at=G%*%mt
#' Rt=array(G%*%t(G)+diag(c(0.1,0.1)),c(2,2,T))
#'
#' smoothed_values=generic_smoother(mt,Ct,at,Rt,G)
generic_smoother= function(mt,Ct,at,Rt,G){
  T=dim(mt)[2]
  n=dim(mt)[1]
  mts <- mt
  Cts <- Ct

  var_index=matrix(apply(Ct,3,diag),n,T)!=0

  for(t in (T-1):1){
    var_ref=var_index[,t]
    restricted_Rt=Rt[var_ref,var_ref,t+1]
    restricted_Ct=Ct[var_ref,var_ref,t]
    simple_Rt_inv=restricted_Ct%*%t(G[var_ref,var_ref])%*%solve(restricted_Rt)

    mts[var_ref,t] <- mt[var_ref,t] + simple_Rt_inv%*%(mts[var_ref,t+1] - at[var_ref,t+1])
    Cts[var_ref,var_ref,t] <- restricted_Ct - simple_Rt_inv%*%(restricted_Rt - Cts[var_ref,var_ref,t+1])%*%t(simple_Rt_inv)
  }
  return(list('mts'=mts,'Cts'=Cts))
}

#' poisson_filter
#'
#' Filtering function for the Poisson model.
#'
#' @param y Vector: The observed value at time t.
#' @param m0 Vector: The prior mean for the latent vector at time t.
#' @param C0 Matrix: The prior covariance matrix for the latent vector at time t.
#' @param FF Matrix: The regression matrix at time t. Should be compatible with the dimensions of m0 and y, i.e., if m0 has dimension n and y has dimension m, FF should be n x m.
#' @param G Matrix: The state evolution matrix.
#' @param D Matrix: The discount factor matrix at time t.
#' @param W Matrix: The noise matrix at time t.
#' @param offset Vector: Same dimension as y. A vector contaning the offset at time t.
#' @param parms list: a list contaning extra arguments. In this model, extra parameters are not used.
#'
#' @return A list containing:
#' \itemize{
#'  \item at: One-step-ahead mean for the latent vectors.
#'  \item Rt: One-step-ahead covariance matrix for the latent vectors.
#'  \item ft: One-step-ahead linear predictor.
#'  \item qt: One-step-ahead covariance matrix for the linear predictior.
#'  \item a: The alpha parameter of the compatibilized gamma prior (see Ref. Raíra).
#'  \item b: The beta parameter of the compatibilized gamma prior (see Ref. Raíra).
#'  \item a.post: The alpha parameter of the gamma posteirior (see Ref. Raíra).
#'  \item b.post: The beta parameter of the gamma posteirior (see Ref. Raíra).
#'  \item mt: The filtered mean of the latent vector at time t.
#'  \item Ct: The filtered covariance matrix of the latent vector at time t.
#'  \item y: The observed value at time t.
#'  \item parms: The same as the argument.
#'  }
#' @export
#'
#' @examples
#' y=5
#' m0=c(4,1)
#' C0=diag(c(1,1))
#' G=matrix(c(1,0,1,1),2,2)
#' FF=matrix(c(1,0),2,1)
#' D=diag(c(0.1,0.1),2,2)+1
#' W=diag(c(0,0),2,2)
#' offset=1
#'
#' filtered_data=poisson_filter(y,m0,C0,FF,G,D,W,offset)
poisson_filter = function(y,m0,C0,FF,G,D,W,offset=1,parms=list()){
  at <- G%*%m0
  Rt <-G%*%C0%*%(t(G))*D+W

  reduc_RFF=Rt%*%FF

  # One-step-ahead prediction
  ft <- t(FF)%*%at + log(offset)
  qt <- t(FF)%*%reduc_RFF

  # Compatibilizing priors

    a <- (1/qt)
    b <- (exp(-ft -0.5*qt)/(qt))

  # Calculating posterior

  a.post <- a + y
  b.post <- b + 1

  # Compatibilizing posterior

    gt <- log(a.post/b.post) + 1/(2*a.post)
    pt <- (2*a.post-1)/(2*a.post^2)

  mt <- at+reduc_RFF*as.vector((gt-ft)*(1/(qt)))
  if(length(qt)>1){
    Ct <- Rt - (reduc_RFF%*%t(reduc_RFF))%*%diag((1 - pt/qt)*(1/qt))
  }else{
    Ct <- Rt - (reduc_RFF%*%t(reduc_RFF))*as.vector((1 - pt/qt)*(1/qt))
  }

  return(list('at'=at,   'Rt'=Rt,
              'ft'=ft,   'Qt'=qt,
              'a'=a,     'b'=b,
              'a.post'=a.post,'b.post'=b.post,
              'mt'=mt,   'Ct'=Ct,
              'y'=y,'parms'=parms))
}

#' poisson_pred
#'
#' Calculate the values for the predictive distribuition given the values of the parameter of the conjugated distribuition of the linear predictor.
#' It's worth noting that, since the conjugated distribuition of the linear predictior is Gamma (see Ref. Raíra), then the predictive distribuition is Negative Binomial.
#'
#' @param model List: contaning the parameters of the conjugated distribuition for the linear predictor (may be based on the prior, filtered or smoothed distribuition).
#' @param IC_prob Numeric: the desired credibility for the credibility interval
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item pred vector/matrix: the mean of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item var.pred vector/matrix: the variance of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item icl.pred vector/matrix: the percentile of 100*((1-IC_prob)/2)% of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item icu.pred vector/matrix: the percentile of 100*(1-(1-IC_prob)/2)% of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#' }
#' @export
#'
#' @examples
#' # A fitted model shoulb be used as argument, but you can also pass only the parameter themselves.
#'
#' model=list(
#' 'a'=c(1:3),
#' 'b'=c(3:1))
#' )
#'
#' poisson_pred(model)
poisson_pred=function(model,IC_prob=0.95){
  a=model$a
  b=model$b
  list(
    'pred'     = a/ b,
    'var.pred' = a*(b+1)/(b)^2,
    'icl.pred' = qnbinom((1-IC_prob)/2, a, (b/(b +1))),
    'icu.pred' = qnbinom(1-(1-IC_prob)/2, a, (b/(b +1)))
  )
}

#' poisson_fit
#'
#' Fit the  poisson model giver the observed value and the model parameters.
#'
#' @param y Matrix: The observed data. It's dimension shoulb be T x m, where T is the length of the time series and m is the number of outcomes at each time.
#' @param m0 Vector: The prior mean for the latent vector.
#' @param C0 Matrix: The prior covariance matrix for the latent vector.
#' @param FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#' @param G Matrix: The state evolution matrix.
#' @param D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#' @param W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#' @param offset Matrix: The offset of the model. It's dimension should be the same as y.
#' @param IC_prob Numeric: the desired credibility for the credibility interval.
#' @param parms list: a list contaning extra arguments. In this model, extra parameters are not used.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item mt Matrix: The filtered mean of the latent variables for each time. Dimensions are n x T.
#'    \item Ct Array: A 3D-array containing the filtered covariance matrix of the latent variable for each time. Dimensions are n x n x T.
#'    \item ft Matrix: The one-step-ahead linear predictor for each time. Dimensions are m x T.
#'    \item qt Array: A 3D-array containing the one-step-ahead covariance matrix for the linear predictor for each time. Dimensions are m x T.
#'    \item a Matrix: The alpha parameter for the gamma prior for each time. Dimensions are m x T.
#'    \item b Matrix: The beta parameter for the gamma prior for each time. Dimensions are m x T.
#'    \item a.post Matrix: The alpha parameter for the gamma posteior for each time. Dimensions are m x T.
#'    \item b.post Matrix: The beta parameter for the gamma posteior for each time. Dimensions are m x T.
#'    \item FF Array: The same as the argument (same values).
#'    \item G Matrix: The same as the argument (same values).
#'    \item D Array: The same as the argument (same values).
#'    \item W Array: The same as the argument (same values).
#'    \item pred Matrix: The one-step-ahead predictions for each time. Dimensions are m x T.
#'    \item var.pred Matrix: The variance for the one-step-ahead predictions for each time. Dimensions are m x T. Note that, in the multivariate Poisson case, the series are supossed independent, so, in particular, they are uncorrelated.
#'    \item icl.pred Matrix: The lower credibility interval for the prediction at each time. Dimensions are m x T.
#'    \item icu.pred Matrix: The upper credibility interval for the prediction at each time. Dimensions are m x T.
#'    \item mts Matrix: The smoothed mean of the latent variables for each time. Dimensions are n x T.
#'    \item Cts Array: A 3D-array containing the smoothed covariance matrix of the latent variable for each time. Dimensions are n x n x T.
#'    \item IC_prob Numeric: Deprecated
#'    \item offset Vector: The same as the argument (same values).
#'    \item data_out Matrix: The same as the argument y (same values).
#'    \item parms: The same as the argument.
#' }
#' @export
#'
#' @examples
#' # Not ideal way: should use fit_model function.
#' T=200
#' w=(200/40)*2*pi
#' y= matrix(rpois(T,20*(sin(w*1:T/T)+2)),T,1)
#' m0=c(0,0,0)
#' C0=diag(c(1,1,1))
#' G=as.matrix(Matrix::bdiag(1,matrix(c(cos(w),sin(w),-sin(w),cos(w)),2,2)))
#' FF=array(matrix(c(1,1,0),3,1),c(3,1,T))
#' D=array(diag(c(0.1,0,0)),c(3,3,T))+1
#' W=array(diag(c(0,0,0)),c(3,3,T))
#' offset=matrix(1,T,1)
#'
#' fitted_data=GDLM::poisson_fit(y=y,m0 = m0, C0 = C0, FF=FF,G=G,D=D,W=W, offset=offset, IC_prob=0.95)
#'
#' plot(y)
#' lines(fitted_data$pred[1,])
poisson_fit <- function(y,m0 = 0, C0 = 1, FF,G,D,W, offset, IC_prob=0.95,parms=list()){
  T <- dim(y)[1]
  n <- dim(FF)[1]
  r <- dim(FF)[2]

  D.aux <- D
  D <- ifelse(D.aux == 0, 1, D.aux)

  m0 <- matrix(m0,n,1)
  C0 <- C0
  mt <- matrix(0,nrow=n,ncol=T)
  Ct <- array(rep(diag(n),T),dim=c(n,n,T))
  Rt <- array(rep(0,T),dim=c(n,n,T))
  ft <- matrix(0,nrow=T,ncol=r)
  at <- matrix(0,nrow=n,ncol=T)
  qt <-  array(0,dim=c(r,r,T))

  #r=1
  at <- matrix(0, nrow=n, ncol=T)
  mt <- matrix(0, nrow=n, ncol=T)
  ft <- matrix(0, nrow=r, ncol=T)
  qt <- matrix(0, nrow=r, ncol=T)
  Ct <- array(rep(diag(n),T),dim=c(n,n,T))
  Rt <- array(rep(diag(n),T),dim=c(n,n,T))
  pred = var.pred = icl.pred = icu.pred = matrix(0, nrow=r, ncol=T)
  media.post = var.post = icu.post = icl.post = matrix(0, nrow=r, ncol=T)
  a = b = matrix(0,r,T)
  a.post = b.post = matrix(0,r,T)

  norm_ic=qnorm(1-(1-IC_prob)/2)

  # Prior

  last_m=m0
  last_C=C0

  start<- proc.time()
  for(t in 1:T){
    filter=poisson_filter(y[t,],last_m,last_C,FF[,,t],G,D[,,t],W[,,t],offset[t,])

    at[,t]  <- filter$at
    Rt[,,t] <- filter$Rt
    ft[,t]  <- filter$ft
    qt[,t]  <- filter$Qt
    a[,t]  <- filter$a
    b[,t] <- filter$b
    a.post[,t]  <- filter$a.post
    b.post[,t] <- filter$b.post
    mt[,t]  <- filter$mt
    Ct[,,t] <- filter$Ct

    last_m=mt[,t]
    last_C=Ct[,,t]

    media.post[t] <- a.post[t]/(b.post[t])
    var.post[t] <- a.post[t]/(b.post[t]^2)
    icu.post[t] <- media.post[t] + norm_ic*sqrt(var.post[t])
    icl.post[t] <- media.post[t] - norm_ic*sqrt(var.post[t])

    prediction=poisson_pred(filter,IC_prob)

    pred[,t]     <- prediction$pred
    var.pred[,t] <- prediction$var.pred
    icl.pred[,t] <- prediction$icl.pred
    icu.pred[,t] <- prediction$icu.pred
  }

  smoothed=generic_smoother(mt,Ct,at,Rt,G)
  mts <- smoothed$mts
  Cts <- smoothed$Cts


  result <- list('mt'=mt,'Ct'=Ct,
                 'ft'=ft,'qt'=qt,
                 'a'=a,'b'=b,
                 'a.post'=a.post, 'b.post'=b.post,
                 'FF'=FF, 'G'=G, 'D'=D,'W'=W,
                 'pred'=pred, 'var.pred'=var.pred, 'icl.pred'=icl.pred, 'icu.pred'=icu.pred,
                 'mts'=mts, 'Cts'=Cts ,
                 'IC_prob'=IC_prob,'offset'=offset,
                 'data_out'=y,'parms'=parms)
  return(result)

}


#' @export
poisson_LB_kernel=list('fit'=poisson_fit,
                    'filter'=poisson_filter,
                    'smoother'=generic_smoother,
                    'pred'=poisson_pred,
                    'multi_var'=FALSE)


#' @export
kernel_list=list('poisson_lb'=poisson_LB_kernel)
