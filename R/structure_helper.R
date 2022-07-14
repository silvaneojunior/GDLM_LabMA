#' polynomial_block
#'
#' Creates the structure for a polynomial block with desired order.
#'
#' @param order Positive integer: The order of the polimial structure.
#' @param values Matrix, vector or scalar: The values to be used on the first row of the regression matrix. If values is a matrix, it's dimensions should be k x t, where k is the number of outcomes of the model and t is the length of the outcome. If values is a vector and it's dimesions is equal to k (or k is Null), then it's values will be repeated for each times.  If values is a vector and it's dimesions is not equal to k (and k is not Null), then it's values will be repeated for each outcome (it's length will be used as time length). If values is a scalar, it's value will be used for all series and all times.
#' @param name String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.
#' @param D Array, Matrix, vector or  scalar: The values for the discount factors at each time. If D is a array, it's dimensions should be nxnxt, where n is the order of the polynomial block and t is the length of the outcomes. If D is a matrix, it's dimesions should be nxn and it's values will be used for each time. If D is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D in the diagonal.
#' @param W Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor at each time. If W is a array, it's dimensions should be nxnxt, where n is the order of the polynomial block and t is the length of the outcomes. If W is a matrix, it's dimesions should be nxn and it's values will be used for each time. If W is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of W in the diagonal.
#' @param m0 Vector or scalar: The prior mean for the latent variables associated with this block. If m0 is a vector, it's dimesion should be equal to the order of the polynomial block. If m0 is a scalar, it's value will be used for all latent variables.
#' @param C0 Matrix, vector or scalar: The prior covariance matrix for the latent variables associated with this block. If C0 is a matrix, it's dimesions should be nxn. If W is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0 in the diagonal.
#' @param k Positive integer: The number of outcomes in the model. Must be consistent with the dimension of values.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#'    \item G Matrix: The state evolution matrix.
#'    \item D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#'    \item W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#'    \item m0 Vector: The prior mean for the latent vector.
#'    \item C0 Matrix: The prior covariance matrix for the latent vector.
#'    \item names list: A list containg the variables indexes by their name.
#'    \item order Positive integer: Same as argument.
#'    \item n Positive integer: The number of latent variables associated with this block (same value as order).
#'    \item t Positive integer: The number of time steps associated with this block. If 1, the block is compatible with blocks of any time length, but if t is greater than 1, this block can only be used with blocks of the same time length.
#'    \item k Positive integer: The number of outcomes associated with this block. This block can only be used with blocks with the same outcome length.
#' }
#'
#' @export
#' @examples
#' # Creating a first order structure for a model with 2 outcomes.
#' # One block is created for each outcome, with each block being associated with only one of the outcomes.
#' level_1=polynomial_block(order=1,values=c(1,0))
#' level_2=polynomial_block(order=1,values=c(1,0))
#'
#' # Creating a block with shared effect between the oucomes
#' level_3=polynomial_block(order=2,values=c(1,1))
polynomial_block <- function(order,values=1,name='Var_Poly',D=1,W=0,m0=0,C0=1,k=NULL){
  G=diag(order)
  if(is.null(k)){
    byrow_flag=FALSE
    multiple_block=FALSE
    t=ifelse(is.null(dim(values)),1,dim(values)[2])
    k=ifelse(is.null(dim(values)),length(values),dim(values)[1])
    if(k==1 & length(D)>1){
      if(length(dim(D))>1){
        k=dim(D)[1]
      }else{
        k=length(D)
      }
    }
    if(k==1 & length(W)>1){
      if(length(dim(W))>1){
        k=dim(W)[1]
      }else{
        k=length(W)
      }
    }
  }else{
    byrow_flag=TRUE
    multiple_block=length(dim(values))<2
    t=ifelse(is.null(dim(values)),length(values),dim(values)[2])
  }
  if(t==1 & length(dim(D))==3){
    t=dim(D)[3]
  }
  if(t==1 & length(dim(W))==3){
    t=dim(W)[3]
  }
  #FF=matrix(c(ifelse(is.na(values),0,values),rep(0,(order-1)*t)),order,t,byrow = TRUE)
  FF=array(0,c(order,k,t))
  FF[1,,]=matrix(values,k,t,byrow=byrow_flag)

  if(order==2){
    G[1,2]=1
  }else{if(order>2){
    diag(G[1:(order-1),2:order])=1
  }}
  if(length(m0)<order){
    m0=rep(m0,order)
  }
  if(length(C0)==1){
    C0=diag(order)*C0
  }
  if(length(dim(C0))==1){
    C0=diag(C0)
  }
  if(length(dim(C0))>2){
    stop(paste0('ERROR: C0 must be a matrix, but it has ',length(dim(C0)),' dimensions.'))
  }
  if(any(dim(C0)!=c(order,order))){
    stop(paste0('ERROR: C0 must have dimensions ',order,'x',order,'. Got ',dim(C0)[1],'x',dim(C0)[2],'.'))
  }

  if(length(D)==1){
    D=array(1,c(order,order,t))*D
    D[,,apply(is.na(FF),3,any)]=1
  }else{if(length(dim(D))==0 & length(D)>1){
    D=array(diag(D),c(order,order,t))
    D[,,apply(is.na(FF),3,any)]=1
  }}
  if(length(W)==1){
    W=array(diag(order),c(order,order,t))*W
    W[,,apply(is.na(FF),3,any)]=0
  }else{if(length(dim(W))==0 & length(W)>1){
    W=array(diag(W),c(order,order,t))
    W[,,apply(is.na(FF),3,any)]=0
  }}

  FF=ifelse(is.na(FF),0,FF)

  names=list()
  names[[name]]=c(1:order)
  block=list('FF'=FF,
             'G'=G,
             'D'=D,
             'W'=W,
             'm0'=m0,
             'C0'=diag(order)*C0,
             'names'=names,
             'order'=order,
             'n'=order,
             't'=t,
             'k'=k)

  if(multiple_block){
    block=multiple_block(block,k,values=values)
  }


  return(block)
}


#' harmonic_block
#'
#' Creates the structure for a harmonic block with desired periodicity.
#'
#' @param period Positive integer: The size of the harmonic cycle.
#' @param values Matrix, vector or scalar: The values to be used on the first row of the regression matrix. If values is a matrix, it's dimensions should be k x t, where k is the number of outcomes of the model and t is the length of the outcome. If values is a vector and it's dimesions is equal to k (or k is Null), then it's values will be repeated for each times.  If values is a vector and it's dimesions is not equal to k (and k is not Null), then it's values will be repeated for each outcome (it's length will be used as time length). If values is a scalar, it's value will be used for all series and all times.
#' @param name String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.
#' @param D Array, Matrix, vector or  scalar: The values for the discount factors at each time. If D is a array, it's dimensions should be nxnxt, where n is the order of the polynomial block and t is the length of the outcomes. If D is a matrix, it's dimesions should be nxn and it's values will be used for each time. If D is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D in the diagonal.
#' @param W Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor at each time. If W is a array, it's dimensions should be nxnxt, where n is the order of the polynomial block and t is the length of the outcomes. If W is a matrix, it's dimesions should be nxn and it's values will be used for each time. If W is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of W in the diagonal.
#' @param m0 Vector or scalar: The prior mean for the latent variables associated with this block. If m0 is a vector, it's dimesion should be equal to the order of the polynomial block. If m0 is a scalar, it's value will be used for all latent variables.
#' @param C0 Matrix, vector or scalar: The prior covariance matrix for the latent variables associated with this block. If C0 is a matrix, it's dimesions should be nxn. If W is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0 in the diagonal.
#' @param k Positive integer: The number of outcomes in the model. Must be consistent with the dimension of values.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#'    \item G Matrix: The state evolution matrix.
#'    \item D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#'    \item W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#'    \item m0 Vector: The prior mean for the latent vector.
#'    \item C0 Matrix: The prior covariance matrix for the latent vector.
#'    \item names list: A list containg the variables indexes by their name.
#'    \item period Positive integer: Same as argument.
#'    \item n Positive integer: The number of latent variables associated with this block (2).
#'    \item t Positive integer: The number of time steps associated with this block. If 1, the block is compatible with blocks of any time length, but if t is greater than 1, this block can only be used with blocks of the same time length.
#'    \item k Positive integer: The number of outcomes associated with this block. This block can only be used with blocks with the same outcome length.
#' }
#'
#' @export
#' @examples
#' # Creating seasonal structure for a model with 2 outcomes.
#' # One block is created for each outcome, with each block being associated with only one of the outcomes.
#' season_1=harmonic_block(period=3,values=c(1,0))
#' season_2=harmonic_block(order=6,values=c(1,0))
#'
#' # Creating a block with shared effect between the oucomes
#' season_3=harmonic_block(order=12,values=c(1,1))
harmonic_block <- function(period,values=1,name='Var_Sazo',D=1,W=0,m0=0,C0=1,k=NULL){
  w=2*pi/period
  order=2
  if(is.null(k)){
    byrow_flag=FALSE
    multiple_block=FALSE
    t=ifelse(is.null(dim(values)),1,dim(values)[2])
    k=ifelse(is.null(dim(values)),length(values),dim(values)[1])
    if(k==1 & length(D)>1){
      if(length(dim(D))>1){
        k=dim(D)[1]
      }else{
        k=length(D)
      }
    }
    if(k==1 & length(W)>1){
      if(length(dim(W))>1){
        k=dim(W)[1]
      }else{
        k=length(W)
      }
    }
  }else{
    byrow_flag=TRUE
    multiple_block=length(dim(values))<2
    t=ifelse(is.null(dim(values)),length(values),dim(values)[2])
  }
  if(t==1 & length(dim(D))==3){
    t=dim(D)[3]
  }
  if(t==1 & length(dim(W))==3){
    t=dim(W)[3]
  }
  #FF=matrix(c(ifelse(is.na(values),0,values),rep(0,(order-1)*t)),order,t,byrow = TRUE)
  FF=array(0,c(order,k,t))
  FF[1,,]=matrix(values,k,t)

  G <- matrix(c(cos(w),-sin(w),sin(w),cos(w)),order,order)

  if(length(m0)<order){
    m0=rep(m0,order)
  }
  if(length(C0)==1){
    C0=diag(order)*C0
  }
  if(length(dim(C0))==1){
    C0=diag(C0)
  }
  if(length(dim(C0))>2){
    stop(paste0('ERROR: C0 must be a matrix, but it has ',length(dim(C0)),' dimensions.'))
  }
  if(any(dim(C0)!=c(order,order))){
    stop(paste0('ERROR: C0 must have dimensions ',order,'x',order,'. Got ',dim(C0)[1],'x',dim(C0)[2],'.'))
  }

  if(length(D)==1){
    D=array(1,c(order,order,t))*D
    D[,,apply(is.na(FF),3,any)]=1
  }else{if(length(dim(D))==0 & length(D)>1){
    D=array(diag(D),c(order,order,t))
    D[,,apply(is.na(FF),3,any)]=1
  }}
  if(length(W)==1){
    W=array(diag(order),c(order,order,t))*W
    W[,,apply(is.na(FF),3,any)]=0
  }else{if(length(dim(W))==0 & length(W)>1){
    W=array(diag(W),c(order,order,t))
    W[,,apply(is.na(FF),3,any)]=0
  }}
  FF=ifelse(is.na(FF),0,FF)
  names=list()
  names[[name]]=c(1:order)

  block=list('FF'=FF,
             'G'=G,
             'D'=D,
             'W'=W,
             'm0'=m0,
             'C0'=diag(order)*C0,
             'names'=names,
             'period'=period,
             'n'=order,
             't'=t,
             'k'=k)

  if(multiple_block){
    block=multiple_block(block,k,values=values)
  }
  return(block)
}

transf_block <- function(lag,values,name='Var_Poly_transf',D=1,m0=0,C0=1,W=0){
  G=diag(order)
  t=dim(values)[2]
  k=dim(values)[1]

  x=c(1:lag)
  mat=c()
  for(j in c(1:lag)){
    mat=c(mat,x**j)
  }

  M=matrix(mat,lag,lag,byrow=TRUE)
  pre_time=matrix(0,out_var,lag)
  extended_values=cbind(pre_time,values)
  pre_FF=matrix(0,(lag+1)*out_var,T_final)
  for(t in c(1:T_final)){
    for(out in c(1:out_var)){
      pre_FF[out_var*(k-1)+1:out_var,t]=M%*%t(extended_values[,(t):(t+k-1)])
    }
  }

  for(i in c(1:out_var)){
    placeholder=matrix(0,out_var,T_final)
    placeholder[i,]=vac_flag
    W=array(0,c(1,1,T_final))
    W[,,true_indice_inter]=1
    bloc_final=concat_bloco(bloc_final,
                            gera_bloco_poly(order=1,
                                            values=placeholder,
                                            name='vac_serie_' %>% paste0(i,'_',0),
                                            D=1/1,
                                            m0=0,
                                            C0=0,
                                            W=W))
  }


  FF=array(0,c(order,k,t))
  FF[1,,]=values
  if(order==2){
    G[1,2]=1
  }else{if(order>2){
    diag(G[1:(order-1),2:order])=1
  }}
  if(length(m0)<order){
    m0=rep(m0,order)
  }

  if(length(D)==1){
    D=array(1,c(order,order,t))*D
    D[,,apply(is.na(values),2,any)]=1
  }
  if(length(W)==1){
    W=array(diag(order),c(order,order,t))*W
    W[,,apply(is.na(values),2,any)]=0
  }

  names=list()
  names[[name]]=c(1:order)
  return(list('FF'=FF,
              'G'=G,
              'D'=D,
              'W'=W,
              'm0'=m0,
              'C0'=diag(order)*C0,
              'names'=names,
              'order'=order,
              'n'=order,
              't'=t,
              'k'=k))
}


#' block_merge
#'
#' An auxiliary function to merge blocks.
#'
#' @param ... <undefined class object> or list: A sequence of block to be merged.
#'
#' @return The merged block as a <undefined class object>.
#' @export
#'
#' @examples
#' level_1=polynomial_block(order=1,values=matrix(c(rep(1,T),rep(0,T)),2,T,byrow=TRUE))
#' level_2=polynomial_block(order=2,values=matrix(c(rep(0,T),rep(1,T)),2,T,byrow=TRUE))
#' season_2=harmonic_block(period=20,values=matrix(c(rep(0,T),rep(1,T)),2,T,byrow=TRUE))
#'
#' final_block=block_merge(level_1,level_2,season_2)
block_merge <- function(...){
  blocks=list(...)

  n=0
  t=1
  k=1
  names=list()
  for(block in blocks){
    ref_names=block$names
    for(name in names(ref_names)){
      ref_names[[name]]=ref_names[[name]]+n
    }
    names=c(names,ref_names)
    if(block$t>1){
      if(block$t!=t & t>1){
        stop(paste('Error: Blocks should have same length or length equal 1. Got',block$t,'and',t))
      }
      t=block$t
    }
    n=n+block$n
    k=max(block$k,k)
  }
  for(name in names(names)){
    ref_idx=which(names(names)==name)
    n_names=length(ref_idx)
    if(n_names>1){
      names(names)[ref_idx]=paste0(names(names)[ref_idx],'_',c(1:n_names))
    }
  }

  FF=array(0,c(n,k,t))
  G=matrix(0,n,n)
  D=array(0,c(n,n,t))
  W=array(0,c(n,n,t))
  m0=c()
  C0=matrix(0,n,n)
  position=1
  for(block in blocks){
    current_range=position:(position+block$n-1)
    FF[current_range,,]=block$FF
    G[current_range,current_range]=block$G
    D[current_range,current_range,]=block$D
    W[current_range,current_range,]=block$W
    m0=c(m0,block$m0)
    C0[current_range,current_range]=block$C0
    position=position+block$n
  }
  return(list('FF'=FF,'G'=G,'D'=D,'W'=W,'m0'=m0,'C0'=C0,'n'=n,'t'=t,'k'=k,'names'=names))
}

#' multiple_block
#'
#' An auxiliary function to repeat the same block for multiple outcomes. Useful if a regressor will be used in all outcomes, but with different effects in each one.
#'
#' @param ref_block <undefined class object> or list: The block to be copied.
#' @param k positive integer: The number of outcomes.
#' @param values vector: The values of the regressor at each time. If not provided, the first row of the FF array from the ref_block will be used.
#'
#' @return The merged block as a <undefined class object>.
#' @export
#'
#' @examples
#' level_i=polynomial_block(order=1,values=0)
#' level=multiple_block(level_i)
multiple_block <- function(ref_block,k,values=ref_block$FF[1,1,]){
  if(ref_block$t>1 & ref_block$t!=length(values)){
    stop('ERROR: ref_block have time length greater than 1, but not equal to the length of values')
  }
  n=ref_block$n
  t=length(values)
  if(ref_block$t==1 | ref_block$k==1){
    ref_block$FF=array(ref_block$FF,c(n,k,t))
  }
  aux_func=function(i){
    ref_block_i=ref_block
    ref_block_i$FF[1,-i,]=0
    ref_block_i$FF[1,i,]=values
    ref_block_i$k=k
    ref_block_i$t=t
    ref_block_i$name=paste(ref_block_i$name,i,sep='_')
    return(ref_block_i)
  }
  block_list=lapply(1:k,aux_func)

  return(do.call(block_merge,block_list))
}
