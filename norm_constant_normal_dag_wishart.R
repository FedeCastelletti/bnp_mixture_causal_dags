###################################################################################
## This function computes the log normalizing constant of a (Normal) DAG Wishart ##
###################################################################################

## It can be used to compute marginal likelihoods if applied to prior and posterior DAG Wishart distributions
## Notice that the normalizing constants of the Normal density (prior/posterior on eta) and of the likelihood
## are omitted, because they both cancel out in the M.H. ratio (indeed does not depend on the DAG)

## INPUT:

# j : a node
# A : adjacency matrix of DAG D
# . : hyper-parameters of the (prior/posterior) DAG Wishart (a, U, m_0, a_mu)

## OUTPUT:

# log normalizing constant of a Normal DAG Wishart (up to those terms not depending on the DAG structure)


## Preliminary functions

library(gRbase)

# to find parents of nodes 'set' in graph 'object'

pa = function (set, object){
  amat <- as(object,"matrix")
  rownames(amat) = colnames(amat) = 1:q
  if (is_ugMAT(amat)) 
    return(NULL)
  pa <- names(which(amat[, set] > 0))
  pa <- setdiff(pa, set)
  if (length(pa)) 
    as.numeric(pa)
  else NULL
}

# to find family of nodes 'set' in graph 'object'

fa = function(set, object){
  as.numeric(c(set, pa(set, object)))
}


norm_const_j = function(j, A, a, U, m, a_mu){
  
  q = ncol(A)
  
  # j:  a node (j = 1, ..., q)
  
  pa_j  = pa(j, A)
  
  # prior/posterior hyper-parameters of the prior induced on node j
  
  p_j  = length(pa_j)
  aD_j = a - q + p_j + 1
  
  if(length(pa_j) == 0){
    
    U_jj = U[j,j]
    
    const_j = - lgamma(aD_j/2) + (aD_j/2)*log(U_jj/2)
    
    
  }else{
    
    U_jj = U[j,j] - U[j,pa_j]%*%solve(U[pa_j,pa_j])%*%U[pa_j,j]
    
    const_j = - lgamma(aD_j/2) + (aD_j/2)*log(U_jj/2) + 0.5*log(det(as.matrix(U[pa_j,pa_j])))
    
  }
  
  return(const_j)
  
}


# #############
# ## Example ##
# #############
# 
# library(pcalg)
# library(mvtnorm)
# 
# ## Generate DAG
# 
# q = 10
# 
# D = randomDAG(q, prob = 0.2)
# A = t(as(D, "matrix"))
# A[A != 0] = 1
# 
# ################################
# ## Prior normalizing constant ##
# ################################
# 
# ## Prior hyper-parameters
# 
# a    = q
# U    = diag(rep(1,q))
# m    = rep(0,q)
# a_mu = 1
# 
# norm_const_j(j = 1, A, a, U, m, a_mu)
# norm_const_j(j = 2, A, a, U, m, a_mu)
# norm_const_j(j = 3, A, a, U, m, a_mu)
# norm_const_j(j = 4, A, a, U, m, a_mu)
# norm_const_j(j = 5, A, a, U, m, a_mu)
# 
# ## Now generates true parameters (mu, D, L) and data
# 
# mu = runif(q, -1, 1)
# 
# L = A*matrix(runif(q*q, 0.1, 1), q, q)*sample(c(-1, 1), size = q*q, replace = TRUE); diag(L) = 1
# D = diag(rep(1, q))
# S = solve(t(L))%*%D%*%solve(L)
# 
# library(mvtnorm)
# 
# n = 500
# X = rmvnorm(n, mu, S)
# 
# ####################################
# ## Posterior normalizing constant ##
# ####################################
# 
# ## Posterior hyper-parameters
# 
# x_bar  = colMeans(X)
# X_zero = t((t(X) - x_bar))
# 
# a_tilde    = a + n
# U_tilde    = U + t(X_zero)%*%X_zero + (a_mu*n)/(a_mu + n)*(x_bar - m)%*%t(x_bar - m)
# a_mu_tilde = a_mu + n
# m_tilde  = a_mu/(a_mu + n)*m + n/(a_mu + n)*x_bar
# 
# norm_const_j(j = 1, A, a = a_tilde, U = U_tilde, m = m_tilde, a_mu = a_mu_tilde)
# norm_const_j(j = 2, A, a = a_tilde, U = U_tilde, m = m_tilde, a_mu = a_mu_tilde)
# norm_const_j(j = 3, A, a = a_tilde, U = U_tilde, m = m_tilde, a_mu = a_mu_tilde)
# norm_const_j(j = 4, A, a = a_tilde, U = U_tilde, m = m_tilde, a_mu = a_mu_tilde)
# 
# ## Log marginal likelihood relative to node j can be computed from
# 
# norm_const_j(j = 1, A, a = a, U = U, m = m, a_mu = a_mu) -
#   norm_const_j(j = 1, A, a = a_tilde, U = U_tilde, m = m_tilde, a_mu = a_mu_tilde)
# 
# norm_const_j(j = 2, A, a = a, U = U, m = m, a_mu = a_mu) -
#   norm_const_j(j = 2, A, a = a_tilde, U = U_tilde, m = m_tilde, a_mu = a_mu_tilde)
