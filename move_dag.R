#################################################################
## This function performs a move from a DAG to an adjacent DAG ##
#################################################################

## It can be used as a proposal distribution within the MCMC

## INPUT:

# A : adjacency graph of the initial DAG D_0
# q : number of nodes in the initial DAG

## OUTPUT:

# A_new         : adjacency matrix of a direct successor (adjacent) DAG of D_0
# type.operator : the type of operator applied (see details below)
# nodes         : nodes involved in the local move


n.edge = function(A){
  length(which(A[lower.tri(A)] == 1 | t(A)[lower.tri(A)] == 1))
}

names   = c("action","test","x","y")
actions = c("id","dd","rd")

library(gRbase)

# types are then indexed by (1,2,3) (respectively insert, delete and reverse a directed edge)

move = function(A, q = q){

  A_na = A
  diag(A_na) = NA
  
  id_set = c()
  dd_set = c()
  rd_set = c()
  
  # set of nodes for id (insert directed edge)
  
  set_id = which(A_na == 0, TRUE)
  
  if(length(set_id) != 0){
    id_set = cbind(1, rbind(set_id))
  }
  
  # set of nodes for dd (delete directed edge)
  
  set_dd = which(A_na == 1, TRUE)
  
  if(length(set_dd != 0)){
    dd_set = cbind(2, set_dd)
  }
  
  # set of nodes for rd (reverse directed edge)
  
  set_rd = which(A_na == 1, TRUE)
  
  if(length(set_rd != 0)){
    rd_set = cbind(3, set_rd)
  }
  
  O = rbind(id_set, dd_set, rd_set)

  repeat {
    
    i = sample(dim(O)[1],1)
    
      act_to_exe  = paste0(actions[O[i,1]],"(A=A,c(",as.vector(O[i,2]),",",as.vector(O[i,3]),"))")
      A_succ      = eval(parse(text = act_to_exe))
      act_to_eval = paste0("is.DAG(A_succ)")
      val = eval(parse(text = act_to_eval))
    
    if (val != 0){
      break
    }
  }
  
  A_new = A_succ
  
  return(list(A_new = A_new, type.operator = O[i,1], nodes = O[i,2:3]))
  
}

# 3 actions

# A     adjacency matrix of DAG D
# x,y   nodes involved in the action

id = function(A, nodes){
  x = nodes[1]
  y = nodes[2]
  A[x,y] = 1
  return(A)
}

dd = function(A, nodes){
  x = nodes[1]
  y = nodes[2]
  A[x,y] = 0
  return(A)
}

rd = function(A, nodes){ # reverse D x -> y
  x = nodes[1]
  y = nodes[2]
  A[x,y] = 0 
  A[y,x] = 1
  return(A)
}
