% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/structure_helper.R
\name{polynomial_block}
\alias{polynomial_block}
\title{polynomial_block}
\usage{
polynomial_block(
  order,
  values = 1,
  name = "Var_Poly",
  D = 1,
  W = 0,
  m0 = 0,
  C0 = 1,
  k = NULL
)
}
\arguments{
\item{order}{Positive integer: The order of the polimial structure.}

\item{values}{Matrix, vector or scalar: The values to be used on the first row of the regression matrix. If values is a matrix, it's dimensions should be k x t, where k is the number of outcomes of the model and t is the length of the outcome. If values is a vector and it's dimesions is equal to k (or k is Null), then it's values will be repeated for each times.  If values is a vector and it's dimesions is not equal to k (and k is not Null), then it's values will be repeated for each outcome (it's length will be used as time length). If values is a scalar, it's value will be used for all series and all times.}

\item{name}{String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.}

\item{D}{Array, Matrix, vector or  scalar: The values for the discount factors at each time. If D is a array, it's dimensions should be nxnxt, where n is the order of the polynomial block and t is the length of the outcomes. If D is a matrix, it's dimesions should be nxn and it's values will be used for each time. If D is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D in the diagonal.}

\item{W}{Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor at each time. If W is a array, it's dimensions should be nxnxt, where n is the order of the polynomial block and t is the length of the outcomes. If W is a matrix, it's dimesions should be nxn and it's values will be used for each time. If W is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of W in the diagonal.}

\item{m0}{Vector or scalar: The prior mean for the latent variables associated with this block. If m0 is a vector, it's dimesion should be equal to the order of the polynomial block. If m0 is a scalar, it's value will be used for all latent variables.}

\item{C0}{Matrix, vector or scalar: The prior covariance matrix for the latent variables associated with this block. If C0 is a matrix, it's dimesions should be nxn. If W is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0 in the diagonal.}

\item{k}{Positive integer: The number of outcomes in the model. Must be consistent with the dimension of values.}
}
\value{
A list containing the following values:
\itemize{
   \item FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
   \item G Matrix: The state evolution matrix.
   \item D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
   \item W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
   \item m0 Vector: The prior mean for the latent vector.
   \item C0 Matrix: The prior covariance matrix for the latent vector.
   \item names list: A list containg the variables indexes by their name.
   \item order Positive integer: Same as argument.
   \item n Positive integer: The number of latent variables associated with this block (same value as order).
   \item t Positive integer: The number of time steps associated with this block. If 1, the block is compatible with blocks of any time length, but if t is greater than 1, this block can only be used with blocks of the same time length.
   \item k Positive integer: The number of outcomes associated with this block. This block can only be used with blocks with the same outcome length.
}
}
\description{
Creates the structure for a polynomial block with desired order.
}
\examples{
# Creating a first order structure for a model with 2 outcomes.
# One block is created for each outcome, with each block being associated with only one of the outcomes.
level_1=polynomial_block(order=1,values=c(1,0))
level_2=polynomial_block(order=1,values=c(1,0))

# Creating a block with shared effect between the oucomes
level_3=polynomial_block(order=2,values=c(1,1))
}
