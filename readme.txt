These R codes implement our DP mixture of DAG model.
Specifically:

mcmc_mixture_dags.R : contains the main MCMC algorithm for posterior inference
move_dag.R          : performs one move from a DAG to an adjacent DAG (implements the proposal distribution over the space of DAGs)

sample_from_baseline.R 		   : sample from the baseline over the space of DAGs

norm_constant_normal_dag_wishart.R : computes the prior/posterior normalizing constant of a Normal-DAG-Wishart distribution
sample_normal_dag_wishart.R        : samples from a (prior/posterior) Normal-DAG-Wishart distribution

example_for_mcmc.R : implements mcmc_mixture_dags.R on a simulated dataset
leukemia.R  : implements mcmc_mixture_dags.R on the AML data

leukemia.csv    : the AML dataset