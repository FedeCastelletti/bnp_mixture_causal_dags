#######################################################################
## This is a function to sample from the baseline over the DAG space ##
#######################################################################

## INPUT:

# S    : number of draws
# burn : burn in period
# q    : number of nodes in the DAGs

# b_pi, b_pi : hyper-parameters of the Beta prior on probability of edge inclusion pi

## OUTPUT:

# (S - burn) DAGs sampled from the space of DAGs with q nodes


library(pcalg)
library(gRbase)

source("move_dag.r")

sample_baseline_dags = function(S, burn, q, a_pi, b_pi){
  
  DAG_chain = array(NA, c(q, q, S))
  
  # Set intial value
  
  DAG = matrix(0, q, q)
  
  DAG_chain[,,1] = DAG
  
  for(s in 2:S){
    
    move_star = move(A = DAG, q = q)
    
    DAG_star  = move_star$A_new
    
    
    # Multiplicity correction (log)prior
    
    logprior.new = lgamma(n.edge(DAG_star) + a_pi) + 
      lgamma(q*(q-1)/2 - n.edge(DAG_star) + b_pi - 1)
    
    logprior.old = lgamma(n.edge(DAG) + a_pi) + 
      lgamma(q*(q-1)/2 - n.edge(DAG) + b_pi - 1)
    
    logprior = logprior.new - logprior.old
    
    
    # acceptance ratio
    
    ratio_D = min(0, logprior)
    
    # accept DAG
    
    if(log(runif(1)) < ratio_D){
      
      DAG = DAG_star
      
    }
    
    DAG_chain[,,s] = DAG
    
  }
  
  return(list(DAG_chain = DAG_chain[,,(burn + 1):S]))
  
}
