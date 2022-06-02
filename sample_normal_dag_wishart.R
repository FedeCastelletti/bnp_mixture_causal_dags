#################################################################################
## This is a function to sample from a Normal DAG Wishart (prior or posterior) ##
#################################################################################

## INPUT:

# A       : (q,q) adjacency matrix of the input DAG
# m, a_mu : hyper-parameters of the Normal prior (posterior)
# a, U    : hyper-parameters of the Wishart prior (posterior)

## OUTPUT:

# one draw from the prior (posterior) of
# mu    : a (q,1) vector
# Sigma : a (q,q) matrix

# see also our Supplementary Material for details


library(gRbase)

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

fa = function(set, object){
  as.numeric(c(set, pa(set, object)))
}

sample_omega_mu = function(A, a, U, m, a_mu){
  
  q = ncol(A)
  
  # A:  a DAG, represented by its (q,q) adjacency matrix
  # S:  number of MCMC iterations
  
  sample_chol_eta_j = function(j){
    
    # j:  a node (j = 1, ..., q)
    
    pa_j  = pa(j, A)
    
    # set node hyper-parameters
    
    p_j  = length(pa_j)
    aD_j = a - q + p_j + 1
    
    L_j = matrix(0, 1, ncol = p_j)
    
    if(p_j == 0){
      
      U_jj  = U[j,j]
      D_jj  = (rgamma(1, aD_j/2, U_jj/2))^(-1)
      eta_j = rnorm(1, m[j], sqrt(D_jj/a_mu))
      
    } else{
      
      U_jj  = U[j,j] - U[j,pa_j]%*%solve(U[pa_j,pa_j])%*%U[pa_j,j]
      D_jj  = (rgamma(1, aD_j/2, U_jj/2))^(-1)
      L_j   = rmvnorm(1, -solve(U[pa_j,pa_j])%*%U[pa_j,j], D_jj*solve(U[pa_j,pa_j]))
      eta_j = rnorm(1, m[j] + L_j%*%m[pa_j], sqrt(D_jj/a_mu))
      
    }
    
    return(list(L_j = L_j, D_jj = D_jj, eta_j = eta_j))
    
  }
  
  
  chol_nodes = lapply(X = 1:q, FUN = sample_chol_eta_j)
  j_parents  = sapply(X = 1:q, FUN = pa, object = A)
  
  L   = matrix(0, q, q); diag(L) = 1
  D   = matrix(0, q, q)
  eta = c()
  
    for(j in 1:q){
        
        L[unlist(j_parents[j]),j] = chol_nodes[[j]]$L_j
        D[j,j] = chol_nodes[[j]]$D_jj
        eta[j] = chol_nodes[[j]]$eta_j
        
    }
  
  
  Sigma  = t(solve(L))%*%D%*%solve(L)
  mu     = t(solve(L))%*%eta
  
  
  return(list(Sigma = Sigma, mu = mu))
  
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
# ##########################################
# ## Sample from Normal DAG Wishart prior ##
# ##########################################
# 
# ## Prior hyper-parameters
# 
# a    = q
# U    = diag(rep(1,q))
# m    = rep(0,q)
# a_mu = 1
# 
# sample_omega_mu(A = A, a = a, U = U, m = m, a_mu = a_mu)
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
# ##############################################
# ## Sample from Normal DAG Wishart posterior ##
# ##############################################
# 
# ## Posterior hyper-parameters
# 
# x_bar  = colMeans(X)
# X_zero = t((t(X) - x_bar))
# 
# a_tilde    = a + n
# U_tilde    = U + t(X_zero)%*%X_zero + (a_mu*n)/(a_mu + n)*(x_bar - m)%*%t(x_bar - m)
# a_mu_tilde = a_mu + n
# m_tilde    = a_mu/(a_mu + n)*m + n/(a_mu + n)*x_bar
# 
# sample_omega_mu(A = A, a = a_tilde, U = U_tilde, m = m_tilde, a_mu = a_mu_tilde)
