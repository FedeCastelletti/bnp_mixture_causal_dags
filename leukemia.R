# Leukemia Data

library(data.table)
library(pcalg)

data = read.csv("leukemia.csv", sep = ";")

X = data[,-1]

# Standardize the data

q = ncol(X)
X = as.matrix(X)

m  = colMeans(X)
sd = apply(X = X, FUN = sd, MARGIN = 2)

X = t((t(X) - m)/sd)

# Set hyperparameters of the gamma prior on precision parameter alpha0

a_alpha = 4
b_alpha = 1

# Set hyperparameters of the beta prior on prior probabilities of edge inclusion pi

a_pi = 1
b_pi = 10

# Set the number of MCMC iterations and burn in period

S    = 120000
burn = 20000

source("mcmc_mixture_dags.r")

# Response variables of interest are:

y_set = c(1,2,3)


t0 = proc.time()
out_mcmc = mcmc_mixture_dags(X = X, S = S, burn = burn, a_alpha = a_alpha, b_alpha = b_alpha, a_pi = a_pi, b_pi = b_pi, y_set = y_set)
t1 = proc.time() - t0


str(out_mcmc)

simil_probs = out_mcmc$simil_probs
graph_probs = out_mcmc$graph_probs


est_causal_1 = out_mcmc$est_causal$`response y = 1`
est_causal_2 = out_mcmc$est_causal$`response y = 2`
est_causal_3 = out_mcmc$est_causal$`response y = 3`


from_simil_to_clust = function(simil_probs){
  
  simil_mat = round(simil_probs)
  
  clust_ind = c()
  
  for(i in ncol(simil_mat):1){
    
    clust_ind[simil_mat[i,] == 1] = i
    
  }
  
  clust_ind = as.factor(clust_ind)
  
  levels(clust_ind) = 1:(length(levels(clust_ind)))
  
  return(clust_ind)
  
}

clus_hat = from_simil_to_clust(simil_probs)

n = nrow(X)

ordered_individuals = unlist(sapply(1:length(table(clus_hat)), function(k) c(which(clus_hat == k))))

simil_ordered = simil_probs[ordered_individuals,ordered_individuals]



#############
## Results ##
#############

#######################
## Similarity matrix ##
#######################

grid_labs = c(1, round(n/4), n/2, round(3*n/4), n)
pos_labs  = grid_labs/n; pos_labs[1] = 0

library(fields)

#pdf("simil_map.pdf", width = 7, height = 6.3)

colori = colorRampPalette(c('white','black'))
par(mar = c(5,5,1,1), oma = c(0.5,0.5,0.5,0.5), cex = 1, mgp = c(3.5,1,0), mfrow = c(1,1), mai = c(1,1,0.5,0.5))
image.plot(simil_ordered, col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "subjects (i)", ylab = "subjects (i')", axes = F, horizontal = F, legend.shrink = 1, border = "black", lwd = 2)
axis(1, at = pos_labs, lab = grid_labs, las = 2)
axis(2, at = pos_labs, lab = grid_labs, las = 2)
abline(h = 1.0025)
abline(v = 1.0025)
abline(h = -0.002)
abline(v = -0.002)

#dev.off()


######################
## Graph estimation ##
######################

## Subject-specific graphs

prob_1 = graph_probs[,,1]
prob_6 = graph_probs[,,6]

rownames(prob_1) = rownames(prob_6) = colnames(X)
colnames(prob_1) = colnames(prob_6) = colnames(X)

#pdf("probs_two.pdf", width = 9.8, height = 4.8)

par(oma = c(3,3,1,1.5), mar = c(5,4,0.3,1), mai = c(0.5,0.5,0.1,0.5))
set.panel(1,2)

colori = colorRampPalette(c('white','black'))

image.plot(t(prob_1), col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "", axes = F, horizontal = F, legend.shrink = 1)
axis(1, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2)
axis(2, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2)
abline(h = 1.030)
abline(v = 1.030)
abline(h = -0.03)
abline(v = -0.03)

image.plot(t(prob_6), col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "", axes = F, horizontal = F, legend.shrink = 1)
axis(1, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2)
axis(2, at = seq(0, 1, l = q), lab = NA, las = 2, srt = 35, cex = 1.2)
abline(h = 1.030)
abline(v = 1.030)
abline(h = -0.03)
abline(v = -0.03)

#dev.off()


####################
## Causal effects ##
####################

## Sort individuals according to the "estimated" cluster allocations

est_causal_1_sort = est_causal_1[ordered_individuals,]
est_causal_2_sort = est_causal_2[ordered_individuals,]
est_causal_3_sort = est_causal_3[ordered_individuals,]

colnames(X)[1:3]

# to set a range (equal for all heat maps)

min(c(est_causal_1, est_causal_2, est_causal_3))
max(c(est_causal_1, est_causal_2, est_causal_3))

#pdf("causal_1.pdf", width = 10, height = 5.5)

par(mar = c(4.5,5.5,3.5,2), oma = c(0,0,0,0), mai = c(1.02,1,0.82,0.4))
set.panel(1,1)

colori = colorRampPalette(c('cyan4','white','brown'))

image.plot(est_causal_1_sort, col = colori(100), zlim = c(-0.5,0.5), cex.sub = 1, xlab = "subjects (i)", axes = F, horizontal = F, legend.shrink = 1)
axis(1, at = pos_labs, lab = grid_labs, las = 2)
axis(2, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2)
mtext(line = 1, side = 3, expression(AKT), outer = F)
abline(h = 1.030)
abline(v = 1.002)
abline(h = -0.03)
abline(v = -0.002)

#dev.off()


#pdf("causal_2.pdf", width = 10, height = 5.5)

par(mar = c(4.5,5.5,3.5,2), oma = c(0,0,0,0), mai = c(1.02,1,0.82,0.4))
set.panel(1,1)

colori = colorRampPalette(c('cyan4','white','brown'))

image.plot(est_causal_2_sort, col = colori(100), zlim = c(-0.5,0.5), cex.sub = 1, xlab = "subjects (i)", axes = F, horizontal = F, legend.shrink = 1)
axis(1, at = pos_labs, lab = grid_labs, las = 2)
axis(2, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2)
mtext(line = 1, side = 3, expression(AKT.p308), outer = F)
abline(h = 1.030)
abline(v = 1.002)
abline(h = -0.03)
abline(v = -0.002)

#dev.off()


#pdf("causal_3.pdf", width = 10, height = 5.5)

par(mar = c(4.5,5.5,3.5,2), oma = c(0,0,0,0), mai = c(1.02,1,0.82,0.4))
set.panel(1,1)

colori = colorRampPalette(c('cyan4','white','brown'))

image.plot(est_causal_3_sort, col = colori(100), zlim = c(-0.5,0.5), cex.sub = 1, xlab = "subjects (i)", axes = F, horizontal = F, legend.shrink = 1)
axis(1, at = pos_labs, lab = grid_labs, las = 2)
axis(2, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2)
mtext(line = 1, side = 3, expression(AKT.p473), outer = F)
abline(h = 1.030)
abline(v = 1.002)
abline(h = -0.03)
abline(v = -0.002)

#dev.off()
