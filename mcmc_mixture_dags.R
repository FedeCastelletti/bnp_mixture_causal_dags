#####################################################################
## MCMC scheme for posterior inference of DP mixture of DAG models ##
#####################################################################

## Required packages:

library(pcalg)
library(gRbase)
library(mvtnorm)
library(abind)

## INPUT:

# X    :  (n,q) data matrix
# S    :  number of MCMC iterations
# burn :  number of iterations to be discarded (burn-in period)

# a_alpha    : shape hyper-parameter of the Gamma prior on precision parameter alpha_0
# b_alpha    : rate hyper-parameter of the Gamma prior on precision parameter alpha_0
# b_pi, b_pi : hyper-parameters of the Beta prior on probability of edge inclusion pi

# y_set : response variable for causal-effect evaluation (by default node 1); it can be any subset of {1,...,q}


## OUTPUT:

# simil_probs : an (n,n) posterior similarity matrix between individuals
# graph_probs : a (q,q,n) array with n (q,q) matrices collecting subject-specific posterior probabilities of edge inclusion
# est_causal  : a list of length |y_set|; each element is an (n,q) matrix with subject-specific causal effects computed across nodes {1,...,q} relative to one response variable in y_set


mcmc_mixture_dags = function(X, S, burn, a_alpha, b_alpha, a_pi, b_pi, y_set = 1){
  
  # Required auxiliary functions
  
  source("sample_normal_dag_wishart.r")
  source("norm_constant_normal_dag_wishart.r")
  source("move_dag.r")
  source("sample_from_baseline.r")
  
  n = nrow(X)
  q = ncol(X)
  
  n_base    = S
  burn_base = burn
  
  out_baseline = sample_baseline_dags(S = n_base, burn = burn_base, q = q, a_pi, b_pi)$DAG_chain
  
  # Store space for parameters
  
  Xi_chain = matrix(NA, n, S)
  
  # this is a (n,S) matrix collecting the cluster indicators for each sample unit i
  
  A_chain = vector(mode = "list", length = S)
  
  # this is a list collecting the adjacency matrices of the DAGs visited by the chain
  # the length of the list is S and each element t is a collection (array) of K(t) adjacency matrices
  # clearly, the number of mixture components (clusters) K(t) is random itself
  
  Sigma_chain = vector(mode = "list", length = S)
  mu_chain    = vector(mode = "list", length = S)
  
  alpha_0_chain = rep(NA, S)
  
  
  ## Prior hyper-parameters
  
  U    = diag(1,q)
  a    = q
  m_0  = rep(0, q)
  a_mu = 1
  
  
  ########################
  ## Set initial values ##
  ########################
  
  set.seed(123)
  
  K_inits = 2 ## number of clusters
  
  alpha_0 = 1 ## precision parameter
  
  ## DAGs D_1, ..., D_K
  
  A_0 = array(0, c(q, q, K_inits)); colnames(A_0) = rownames(A_0) = 1:q; A_chain[[1]] = A_0
  
  ## Cluster allocators
  
  set.seed(123)
  
  xi = sample(K_inits, n, replace = TRUE)
  
  while(length(table(xi)) < K_inits){
    
    xi = sample(K_inits, n, replace = TRUE)
    
  }
  
  Xi_chain[,1] = xi
  
  ## DAG parameters
  
  Sigma_0 = array(NA, c(q, q, K_inits))
  L_0     = array(NA, c(q, q, K_inits))
  D_0     = array(NA, c(q, q, K_inits))
  mu_0    = array(NA, c(q, 1, K_inits))
  
  for(k in 1:K_inits){
    
    X_k = matrix(X[xi == k,], ncol = q)
    n_k = nrow(X_k)
    
    x_bar  = colMeans(X_k)
    X_zero = t((t(X_k) - x_bar))
    
    # Posterior hyperparameters for DAG k
    
    a_tilde    = a + n_k
    U_tilde    = U + t(X_zero)%*%X_zero + (a_mu*n_k)/(a_mu + n_k)*(x_bar - m_0)%*%t(x_bar - m_0)
    a_mu_tilde = a_mu + n_k
    m_0_tilde  = a_mu/(a_mu + n_k)*m_0 + n_k/(a_mu + n_k)*x_bar
    
    Post_0 = sample_omega_mu(A_0[,,k], a_tilde, U_tilde, m_0_tilde, a_mu_tilde)
    
    Sigma_0[,,k] = Post_0$Sigma
    mu_0[,,k]    = Post_0$mu
    
  }
  
  Sigma_chain[[1]] = Sigma_0
  mu_chain[[1]]    = mu_0
  
  Dags  = A_0
  Sigma = Sigma_0
  mu    = mu_0
  
  r = table(xi)
  
  
  #####################
  ## MCMC iterations ##
  #####################
  
  for(t in 2:S){
    
    ###################################################
    ## Slice sampler to update cluster indicators Xi ##
    ###################################################
    
    maxCl <- length(r) # maximum number of clusters
    ind <- which(r!=0) # indexes of non empty clusters
    
    r = r[ind]
    
    ## [1] ## Blocked Gibbs to update u and v
    
    if(length(r) == 1){
      
      v = rbeta(length(r), 1 + r, alpha_0)
      
    }else{
      
      v = rbeta(length(r), 1 + r, alpha_0 + c(rev(cumsum(rev(r)))[-1]))
      
    }
    
    v = c(v, rbeta(1, 1, alpha_0))
    
    omega_tmp = c(); omega_tmp = v[1]
    
    
    for(k in 2:(length(r)+1)){
      
      omega_tmp[k] = v[k]*prod(1 - v[1:(k-1)])
      
    }
    
    R = omega_tmp[length(omega_tmp)] # R is the rest, i.e. the weight for potential new clusters
    
    omega <- numeric(maxCl)
    
    omega[ind] = omega_tmp[-length(omega_tmp)]
    
    u = stats::runif(n)*omega[xi]
    u_star = min(u)
    
    h <- 0 # the number of new non empty clusters
    
    while(R > u_star){
      
      h <- h+1
      beta_temp <- stats::rbeta(n = 1, shape1 = 1, shape2 = alpha_0)
      
      omega = c(omega, R*beta_temp) # weight of the new cluster
      
      R <- R * (1 - beta_temp) # remaining weight
      
      Dag_star   = out_baseline[,,sample(n_base - burn_base, 1)]
      prior_star = sample_omega_mu(A = Dag_star, a = a, U = U, m = m_0, a_mu = a_mu)
      
      Sigma_star  = prior_star$Sigma
      mu_star     = prior_star$mu
      
      Dags  = abind(Dags, Dag_star)
      Sigma = abind(Sigma, Sigma_star)
      mu    = abind(mu, mu_star)
      
    }
    
    ## [2] ## Update of indicator variables xi
    
    K_star = dim(Dags)[3]
    
    probs = sapply(1:(K_star), function(k) dmvnorm(X, mean = mu[,,k], sigma = Sigma[,,k])*(u <= omega[k]))
    
    xi_star = sapply(1:n, function(i) sample(1:(K_star), size = 1, prob = probs[i,]))
    
    labs = as.integer(names(table(xi_star)))
    
    K_star = length(labs)
    
    Dags  = array(Dags[,,labs], c(q, q, K_star))
    Sigma = array(Sigma[,,labs], c(q, q, K_star))
    mu    = array(mu[,,labs], c(q, 1, K_star))
    
    xi_star = as.factor(xi_star); levels(xi_star) = 1:K_star   # update labels
    
    r = table(xi_star)
    
    omega = omega[labs]
    
    xi = c(xi_star)
    
    K = dim(Dags)[3]
    
    
    ###############################
    ## Update of alpha_0 given K ##
    ###############################
    
    eta = rbeta(1, alpha_0 + 1, n)
    
    alpha_0 = c(rgamma(1, shape = a_alpha + K, rate = b_alpha - log(eta)), rgamma(1, shape = a_alpha + K - 1, rate = b_alpha - log(eta)))[sample(c(1,2), 1, prob = c(a_alpha + K - 1, n*(b_alpha - log(eta))))]
    
    alpha_0_chain[t] = alpha_0
    
    
    ##################################################
    ## Update DAGs D_1, ..., D_K and DAG parameters ##
    ##################################################
    
    set = 1:K
    
    for(k in set){
      
      DAG = Dags[,,k]
      
      move_star = move(A = DAG, q = q)
      
      DAG_star   = move_star$A_new
      type.op    = move_star$type.operator
      nodes_star = move_star$nodes
      
      X_k = matrix(X[xi == k,], ncol = q)
      n_k = nrow(X_k)
      
      x_bar  = colMeans(X_k)
      X_zero = t((t(X_k) - x_bar))
      
      # Posterior hyper-parameters for DAG k
      
      a_tilde    = a + n_k
      U_tilde    = U + t(X_zero)%*%X_zero + (a_mu*n_k)/(a_mu + n_k)*(x_bar - m_0)%*%t(x_bar - m_0)
      a_mu_tilde = a_mu + n_k
      m_0_tilde  = a_mu/(a_mu + n_k)*m_0 + n_k/(a_mu + n_k)*x_bar
      
      
      # Multiplicity correction (log)prior
      
      logprior.new = lgamma(n.edge(DAG_star) + a_pi) + 
        lgamma(q*(q-1)/2 - n.edge(DAG_star) + b_pi - 1)
      
      logprior.old = lgamma(n.edge(DAG) + a_pi) + 
        lgamma(q*(q-1)/2 - n.edge(DAG) + b_pi - 1)
      
      logprior = logprior.new - logprior.old
      
      
      # Distinguish 3 cases:
      
      if(type.op == 1){
        
        # (1) Insert a directed edge
        
        marg_star = norm_const_j(j = nodes_star[2], A = DAG_star, U = U, a = a, m = m_0, a_mu = a_mu) -
                    norm_const_j(j = nodes_star[2], A = DAG_star, U = U_tilde, a = a_tilde, m = m_0_tilde, a_mu = a_mu_tilde)
        
        marg      = norm_const_j(j = nodes_star[2], A = DAG, U = U, a = a, m = m_0, a_mu = a_mu) -
                    norm_const_j(j = nodes_star[2], A = DAG, U = U_tilde, a = a_tilde, m = m_0_tilde, a_mu = a_mu_tilde)
        
      }else{
        
        if(type.op == 2){
          
          # (2) Delete a directed edge
          
          marg_star = norm_const_j(j = nodes_star[2], A = DAG_star, U = U, a = a, m = m_0, a_mu = a_mu) -
                      norm_const_j(j = nodes_star[2], A = DAG_star, U = U_tilde, a = a_tilde, m = m_0_tilde, a_mu = a_mu_tilde)
          
          marg      = norm_const_j(j = nodes_star[2], A = DAG, U = U, a = a, m = m_0, a_mu = a_mu) -
                      norm_const_j(j = nodes_star[2], A = DAG, U = U_tilde, a = a_tilde, m = m_0_tilde, a_mu = a_mu_tilde)
          
        }else{
          
          # (3) Reverse a directed edge
          
          marg_star = norm_const_j(j = nodes_star[1], A = DAG_star, U = U, a = a, m = m_0, a_mu = a_mu) -
                      norm_const_j(j = nodes_star[1], A = DAG_star, U = U_tilde, a = a_tilde, m = m_0_tilde, a_mu = a_mu_tilde) +
                      norm_const_j(j = nodes_star[2], A = DAG_star, U = U, a = a, m = m_0, a_mu = a_mu) -
                      norm_const_j(j = nodes_star[2], A = DAG_star, U = U_tilde, a = a_tilde, m = m_0_tilde, a_mu = a_mu_tilde)
          
          marg      = norm_const_j(j = nodes_star[1], A = DAG, U = U, a = a, m = m_0, a_mu = a_mu) -
                      norm_const_j(j = nodes_star[1], A = DAG, U = U_tilde, a = a_tilde, m = m_0_tilde, a_mu = a_mu_tilde) +
                      norm_const_j(j = nodes_star[2], A = DAG, U = U, a = a, m = m_0, a_mu = a_mu) -
                      norm_const_j(j = nodes_star[2], A = DAG, U = U_tilde, a = a_tilde, m = m_0_tilde, a_mu = a_mu_tilde)
          
        }
        
      }
      
      # acceptance ratio
      
      ratio_D = min(0, marg_star - marg + logprior)
      
      # accept DAG
      
      if(log(runif(1)) < ratio_D){
        
        Dags[,,k] = DAG_star
        
      }
      
      ############################################################################
      ## Sample from the posterior of Sigma_k and mu_k conditionally on DAG D_k ##
      ############################################################################
      
      Post_Sigma_mu = sample_omega_mu(A = Dags[,,k], a_tilde, U_tilde, m_0_tilde, a_mu_tilde)
      
      Sigma_post = Post_Sigma_mu$Sigma
      mu_post    = Post_Sigma_mu$mu
      
      Sigma[,,k] = Sigma_post
      mu[,,k]    = mu_post
      
      
    }
    
    A_chain[[t]]     = Dags
    Sigma_chain[[t]] = Sigma
    mu_chain[[t]]    = mu
    Xi_chain[,t]     = xi
    
    if(t%%100 == 0) print(paste0("Iteration ", t))
    
    
  }
  
  #####################
  ## End of the MCMC ##
  #####################
  
  #########################
  ## Posterior summaries ##
  #########################
  
  ## Construct subject-specific matrices with posterior probabilities of edge inclusion
  
  graph_probs = array(NA, c(q, q, n))
  
  for(i in 1:n){
    
    probs_i = sapply((burn + 1): S, function(t) c(A_chain[[t]][,,Xi_chain[i,t]]))
    graph_probs[,,i] = matrix(rowMeans(probs_i), q, q)
    
  }
  
  
  ## Compute subject-specific causal effects (for each response variable in y_set)
  
  # function to compute causal effects from Sigma, Dag and y as in Equation ()
  
  all_causal_effects = function(Sigma, Dag, y){
    
    x_set = setdiff(1:q, y)
    
    return(sapply(x_set, function(x) (Sigma[y,fa(set = x, Dag)]%*%solve(Sigma[fa(set = x, Dag),fa(set = x, Dag)]))[1]))
    
  }
  
  # function to compute subject-specific causal effect estimates for response node y
  
  est_causal_y = function(y){
    
    est_causal_list = vector(mode = "list", length = S)
    
    for(s in 1:S){
      
      K = dim(Sigma_chain[[s]])[3]
      
      est_causal_list[[s]] = sapply(1:K, function(k) all_causal_effects(Sigma_chain[[s]][,,k], A_chain[[s]][,,k], y))
      
    }
    
    est_causal = matrix(0, n, q)
    
    for(i in 1:n){
      
      out_ind = sapply((burn + 1):S, function(s) est_causal_list[[s]][,Xi_chain[i,s]])
      
      est_causal[i,setdiff(1:q, y)] = rowMeans(out_ind)
      
    }
    
    return(est_causal)
    
  }
  
  est_causal_y_set = vector(mode = "list", length = length(y_set))
  
  for(j in 1:length(y_set)){
    
    est_causal_y_set[[j]] = est_causal_y(y_set[j])
    
  }
  
  names_list = paste0("response y = ", y_set)
  
  names(est_causal_y_set) = names_list
  
  A_chain = NULL
  Sigma_chain = NULL
  
  
  ## Construct the (n,n) similarity matrix
  
  simil_mat = matrix(0, nrow = n, ncol = n)
  
  for(t in (burn + 1):S){
    
    simil_mat = simil_mat + (matrix(Xi_chain[,t], nrow = n, ncol = n) == t(matrix(Xi_chain[,t], nrow = n, ncol = n)))*1
    
  }
  
  simil_probs = simil_mat/(S-burn)
  
  
  return(list(simil_probs = simil_probs,
              graph_probs = graph_probs,
              est_causal  = est_causal_y_set
              ))
}


# #############
# ## Example ##
# #############

# See "example_for_mcmc.r"