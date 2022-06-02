# Generate two different DAGs

q   = 10   # number of variables
n_1 = 100  # sample size for group 1,2
n_2 = 100

p.star = 2/q # probability of edge inclusion (equal for all DAGs)


library(pcalg)

set.seed(11)
D_1 = randomDAG(q, prob = p.star)
A_1 = t(as(D_1, "matrix"))
A_1[A_1 != 0] = 1

set.seed(12)
D_2 = randomDAG(q, prob = p.star)
A_2 = t(as(D_2, "matrix"))
A_2[A_2 != 0] = 1

library(network)

D_1 = network(A_1, label = 1:q)
D_2 = network(A_2, label = 1:q)

vertex_col = "gray90"

par(mar = c(0,0,0,0), oma = c(0,0,0,0), mfrow = c(1,2))

out = plot.network(D_1, displaylabels = TRUE, vertex.col = vertex_col,
                   mode = "circle",
                   label.pos = 5,
                   usecurve = TRUE, edge.curve = 0, vertex.cex = 3.5,
                   label.cex = 0.8, edge.lwd = 0.1, arrowhead.cex = 1.5)

out = plot.network(D_2, displaylabels = TRUE, vertex.col = vertex_col,
                   mode = "circle",
                   label.pos = 5,
                   usecurve = TRUE, edge.curve = 0, vertex.cex = 3.5,
                   label.cex = 0.8, edge.lwd = 0.1, arrowhead.cex = 1.5)


# Generate parameters and data

B_1 = A_1*matrix(runif(q*q, .1, 1), q, q)*sample(c(-1, 1), size = q*q, replace = TRUE); diag(B_1) = 1
B_2 = A_2*matrix(runif(q*q, .1, 1), q, q)*sample(c(-1, 1), size = q*q, replace = TRUE); diag(B_2) = 1

Sigma_cond_1 = diag(rep(1, q))
Sigma_cond_2 = diag(rep(1, q))

S_1 = solve(t(B_1))%*%Sigma_cond_1%*%solve(B_1)
S_2 = solve(t(B_2))%*%Sigma_cond_2%*%solve(B_2)

mu_1 = runif(q,1,2)
mu_2 = runif(q,-2,-1)

mu_1 = rnorm(q,1,1)
mu_2 = rnorm(q,-1,1)


Sigma_true = list(S_1, S_2)
Dag_true   = list(A_1, A_2)

library(mvtnorm)

X_1 = rmvnorm(n_1, mu_1, S_1)
X_2 = rmvnorm(n_2, mu_2, S_2)


Xi = c(rep(1, n_1), rep(2, n_2)) # (true) cluster indicators

X = rbind(X_1, X_2) # dataset

n = n_1 + n_2 # total sample size

par(mfrow = c(1,1), mar = c(4,4,1,2))
plot(X[,1], xlab = "i", ylab = "y", col = c(rep("blue", n_1), rep("hotpink", n_2), rep("black", n_2)), lwd = 1) # see the heterogeneity in the data

par(mfrow = c(1,1), mar = c(4,4,1,2))
plot(X[,5], xlab = "i", ylab = "y", col = c(rep("blue", n_1), rep("hotpink", n_2), rep("black", n_2)), lwd = 1) # see the heterogeneity in the data


## Run the MCMC

S = 2500
burn = 500

a_pi = 1
b_pi = 2*(q-2)/3

source("mcmc_mixture_dags.r")

a_alpha = 1
b_alpha = 3


t0 = proc.time()
out_mcmc = mcmc_mixture_dags(X = X, S = S, burn = burn, a_alpha = a_alpha, b_alpha = b_alpha, a_pi = a_pi, b_pi = b_pi, y_set = c(1,2))
t1 = proc.time() - t0


###############################################################
## Output (1) : Subject-specific BMA causal effect estimates ##
###############################################################

## True causal effects

all_causal_effects = function(Sigma, Dag, y){
  
  x_set = setdiff(1:q, y)
  
  return(sapply(x_set, function(x) (Sigma[y,fa(set = x, Dag)]%*%solve(Sigma[fa(set = x, Dag),fa(set = x, Dag)]))[1]))
  
}

y = 1

true_causal = matrix(0, n, q)

for(i in 1:n){
  
  true_causal[i,] = c(0, all_causal_effects(Sigma = Sigma_true[[Xi[i]]], Dag = Dag_true[[Xi[i]]], y = y))
  
}

## Extract estimated causal effects and compare with true causal effects

est_causal = out_mcmc$est_causal$`response y = 1`

seq_pos = c(0,0.11,0.22,0.33,0.44,0.55,0.66,0.77,0.88,1)

par(mfrow = c(1,2))

library(fields)
colori = colorRampPalette(c('white','black'))
image.plot(abs(t(est_causal)), col = colori(100), zlim = c(0,max(abs(est_causal), abs(true_causal))), cex.sub = 1, xlab = "j", ylab = "i", axes = F,
           horizontal = F, legend.shrink = 1, axis.args = list(cex.axis = 1), cex.lab = 1)

axis(1, at = seq_pos, lab = 1:q, las = 1)
axis(2, at = seq(1/n*10,1,by=1/n*10), lab = seq(10,n, by = 10), las = 1)

colori = colorRampPalette(c('white','black'))
image.plot(abs(t(true_causal)), col = colori(100), zlim = c(0,max(abs(est_causal), abs(true_causal))), cex.sub = 1, xlab = "j", ylab = "i", axes = F,
           horizontal = F, legend.shrink = 1, axis.args = list(cex.axis = 1), cex.lab = 1)

axis(1, at = seq_pos, lab = 1:q, las = 1)
axis(2, at = seq(1/n*10,1,by=1/n*10), lab = seq(10,n, by = 10), las = 1)


##########################################################
## Output (2) : Subject-specific matrices with edge PPI ##
##########################################################

# PPI for four selected subjects

PPIs = out_mcmc$graph_probs

par(mfrow = c(1,2))

library(fields)

colori = colorRampPalette(c('white','black'))
image.plot(abs(t(PPIs[,,1])), col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "v", ylab = "u", axes = F,
           horizontal = F, legend.shrink = 1, axis.args = list(cex.axis = 1), cex.lab = 1)
axis(1, at = seq_pos, lab = 1:q, las = 1)
axis(2, at = seq_pos, lab = 1:q, las = 1)

colori = colorRampPalette(c('white','black'))
image.plot(abs(t(PPIs[,,101])), col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "v", ylab = "u", axes = F,
           horizontal = F, legend.shrink = 1, axis.args = list(cex.axis = 1), cex.lab = 1)
axis(1, at = seq_pos, lab = 1:q, las = 1)
axis(2, at = seq_pos, lab = 1:q, las = 1)


## Evaluate SHD between subject-specific (estimated and true) DAGs

A_true = abind(A_1, A_2, along = 3)

library(pcalg)

SHDs = c()

for(i in 1:n){
  
  med_i = round(out_mcmc$graph_probs[,,i])
  
  SHDs[i] = shd(as(med_i, "graphNEL"), as(A_true[,,Xi[i]], "graphNEL"))
  
}

SHDs


############################################
## Output (3) Posterior similarity matrix ##
############################################

simil_probs = out_mcmc$simil_probs

par(mfrow = c(1,1))

library(fields)
colori = colorRampPalette(c('white','black'))
par(mar = c(4,4,1,2), oma = c(0.5,0.5,0.5,0.5), cex = 1, mgp = c(2,0.5,0))
image.plot(t(simil_probs), col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "i'", ylab = "i", axes = F,
           horizontal = F, legend.shrink = 1, axis.args = list(cex.axis = 1), cex.lab = 1)

axis(1, at = seq(1/n*10,1,by=1/n*10), lab = seq(10,n, by = 10), las = 1)
axis(2, at = seq(1/n*10,1,by=1/n*10), lab = seq(10,n, by = 10), las = 1)


#######################
## Cluster estimates ##
#######################

library(mcclust)

# Estimated clustering (i and i' in the same cluster iff simil.probs > 0.5)

from_simil_to_clust = function(simil_probs){
  
  simil_mat = round(simil_probs)
  
  clust_ind = c()
  
  for(i in n:1){
    
    clust_ind[simil_mat[i,] == 1] = i
    
  }
  
  clust_ind = as.factor(clust_ind)
  
  levels(clust_ind) = 1:(length(levels(clust_ind)))
  
  return(clust_ind)
  
}

from_simil_to_clust(simil_probs)

# Variation of Information (VI) between true and estimated clustering

vi.dist(Xi, from_simil_to_clust(simil_probs))

