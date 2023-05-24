## ---------------------------
## Project name: TIGER
## Purpose of script: source code
## Author: Chen Chen
## Date Created: 2022-10-05
##
## Copyright (c) Chen Chen, 2022
## Email: cchen22@arizona.edu



# # install cmdstanr and cmdstan
# ## we recommend running this is a fresh R session or restarting your current session
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# library(cmdstanr)
# ## to check C++ toolchain
# check_cmdstan_toolchain()
# ## install cmdstan
# install_cmdstan(cores = 4)


#' TIGER main function
#'
#' @param expr A normalized log-transformed gene expressison matrix. Rows are genes
#' and columns are sampeles (cells).
#' @param prior A prior regulatory network in adjacency matrix format. Rows are TFs
#' and columns target genes.
#' @param method Method used for Bayesian inference. "VB" or "MCMC". Defaults to "VB".
#' @param signed Prior network is signed or not. Defaults to TRUE.
#' @param baseline Include baseline or not. Defaults to TRUE.
#' @param psis_loo Use pareto smoothed importance sampling leave-one-out cross
#' validation to check model fitting or not. Defaults to FALSE.
#' @param seed Seed for reproducible results. Defaults to 123.
#' @param out_path (Optional) output path for CmdStanVB or CmdStanMCMC object. Defaults to NULL.
#' @param out_size Posterior sampling size. Default = 300.
#' @param a_sigma Hyperparameter of error term. Default = 1.
#' @param b_sigma Hyperparameter of error term. Default = 1.
#' @param a_alpha Hyperparameter of edge weight W. Default = 1.
#' @param b_alpha Hyperparameter of edge weight W. Default = 1.
#' @param sigmaZ Standard deviation of TF activity Z. Default = 10.
#' @param sigmaB Standard deviation of baseline term. Default = 1.
#'
#' @return A TIGER list object.
#' * W is the estimated regulatory network, but different from prior network,
#' rows are genes and columns are TFs.
#' * Z is the estimated TF activities, rows are TFs and columns are samples.
#' * TF.name, TG.name, and sample.name are the used TFs, target genes and samples.
#' * If psis_loo is TRUE, loocv is a table of psis_loo result for model checking.
#' * If psis_loo is TRUE, elpd_loo is the Bayesian LOO estimate of the expected log pointwise predictive
#' density, which can be used for Bayesian stacking to handle multi-modality later.
#' @export
#'
#' @examples
TIGER = function(expr,prior,method="VB",
                     signed=TRUE,baseline=TRUE,psis_loo = FALSE,
                     seed=123,out_path=NULL,out_size = 300,
                     a_sigma=1,b_sigma=1,a_alpha=1,b_alpha=1,sigmaZ=10,sigmaB=1){
  # check data
  sample.name = colnames(expr)
  TF.name = intersect(rownames(prior),rownames(expr)) # TF needs to express
  TG.name = intersect(rownames(expr),colnames(prior))
  if (length(TG.name)==0 | length(TF.name)==0){
    stop("No matched gene names in the two inputs...")
  }

  #0. prepare stan input
  if (signed){
    prior = prior.pp(prior[TF.name,TG.name],expr)
    P = prior
  }else{
    P = prior[TF.name,TG.name]
  }
  X = expr[TG.name,]
  n_genes = dim(X)[1]
  n_samples = dim(X)[2]
  n_TFs = dim(P)[1]

  P = as.vector(t(P)) ## row=TG, col=TF
  P_zero = as.array(which(P==0))
  P_ones = as.array(which(P!=0))
  P_negs = as.array(which(P==-1))
  P_poss = as.array(which(P==1))
  P_blur = as.array(which(P==1e-6))
  n_zero = length(P_zero)
  n_ones = length(P_ones)
  n_negs = length(P_negs)
  n_poss = length(P_poss)
  n_blur = length(P_blur)
  n_all = length(P)
  sign = as.integer(signed)
  baseline = as.integer(baseline)
  psis_loo = as.integer(psis_loo)

  data_to_model = list(n_genes = n_genes, n_samples = n_samples, n_TFs = n_TFs,X = as.matrix(X), P = P,
                       P_zero = P_zero, P_ones = P_ones,P_negs = P_negs, P_poss = P_poss, P_blur = P_blur,
                       n_zero = n_zero,n_ones = n_ones,n_negs = n_negs, n_poss = n_poss, n_blur = n_blur,
                       n_all = n_all,sign = sign,baseline = baseline,psis_loo = psis_loo,
                       sigmaZ = sigmaZ, sigmaB = sigmaB,
                       a_sigma = a_sigma,b_sigma = b_sigma,a_alpha = a_alpha,b_alpha = b_alpha)

  #1. compile stan model, only once
  f = cmdstanr::write_stan_file(TIGER_C) # save to .stan file in root folder
  #mod = cmdstanr::cmdstan_model(f,cpp_options = list(stan_threads = TRUE)) # compile stan program, allow within-chain parallel
  mod = cmdstanr::cmdstan_model(f)

  #2. run VB or MCMC
  if (method=="VB"){
    fit <- mod$variational(data = data_to_model, algorithm = "meanfield",seed = seed,
                           iter = 50000, tol_rel_obj = 0.005,output_samples = out_size)
  }else if (method=="MCMC"){
    fit <- mod$sample(data = data_to_model,chains=1,seed = seed,max_treedepth=10,
                      iter_warmup = 1000,iter_sampling=out_size,adapt_delta=0.99)
  }

  ## optional: save stan object
  if (!is.null(out_path)){
    fit$save_object(paste0(out_path,"fit_",seed,".rds"))
  }

  #3. posterior distributions

  ## point summary of W
  W_sample = fit$draws("W",format = "draws_matrix") ## matrix, each row is a sample of vectorized matrix
  one_ind = which(colSums(W_sample==0)<nrow(W_sample)) ## index for all-zero columns
  W_pos = colMeans(W_sample[,one_ind]) ## average W
  W_pos2 = rep(0,ncol(W_sample))
  W_pos2[one_ind] = W_pos
  W_pos = matrix(W_pos2,nrow = n_genes,ncol = n_TFs) ## convert to matrix genes*TFs

  ## point summary of Z
  Z_sample = fit$draws("Z",format = "draws_matrix")
  Z_pos = colMeans(Z_sample) ## average Z
  Z_pos = matrix(Z_pos,nrow = n_TFs,ncol = n_samples) ## convert to matrix TFs*samples

  ## rescale
  IZ = Z_pos*(apply(abs(W_pos),2,sum)/apply(W_pos!=0,2,sum))
  IW = t(t(W_pos)*apply(Z_pos,1,sum)/n_samples)

  ## output
  rownames(IW) = TG.name
  colnames(IW) = TF.name
  rownames(IZ) = TF.name
  colnames(IZ) = sample.name

  # check model fitting
  if (psis_loo){
    message("Pareto Smooth Importance Sampling...")
    loocv = loo::loo(fit$draws("log_lik",format = "draws_array"),
                     r_eff=loo::relative_eff(fit$draws("log_lik",format = "draws_array")),
                     moment_match=TRUE)
    print(loocv)
    elpd_loo = loocv$pointwise[,"elpd_loo"]
  }else{
    loocv = NA
    elpd_loo = NA
  }

  # output
  tiger_fit = list(W = IW, Z = IZ,
                   TF.name = TF.name, TG.name = TG.name, sample.name = sample.name,
                   loocv = loocv, elpd_loo=elpd_loo)

  return(tiger_fit)
}



#' Convert bipartite edge list to adjacency mat
#'
#' @param el An edge list dataframe with three columns. First column is TF name,
#' second column is gene name, and third column is edge weight.
#'
#' @return An adjacency matrix with rows as TFs and columns as genes.
#' @export
#'
#' @examples
el2adj = function(el){
  el = as.data.frame(el)
  all.A = unique(el[,1])
  all.B = unique(el[,2])
  adj =  array(0,dim=c(length(all.A),length(all.B)))
  rownames(adj) = all.A
  colnames(adj) = all.B
  adj[as.matrix(el[,1:2])] = as.numeric(el[,3])
  return(adj)
}


#' Convert a bipartite adjacency matrix to an edgelist
#'
#' @param adj An adjacency matrix, with rows as TFs and columns as genes.
#'
#' @return An edge list dataframe with three columns. First column is TF name,
#' second column is gene name, and third column is edge weight.
#' @export
#'
#' @examples
adj2el = function(adj){
  el = matrix(NA,nrow(adj)*ncol(adj),3)
  el[,1]=rep(row.names(adj), ncol(adj))
  el[,2]=rep(colnames(adj), each=nrow(adj))
  el=as.data.frame(el)
  el[,3]=as.numeric(unlist(c(adj),use.names = F))
  el[,1]=as.character(el[,1])
  el[,2]=as.character(el[,2])
  colnames(el)=c("from","to","weight")
  return(el)
}


#' Convert a bipartite edgelist to regulon
#'
#' @param el An edge list dataframe with three columns. First column is TF name,
#' second column is gene name, and third column is edge weight.
#'
#' @return A VIPER required regulon object
#' @export
#'
#' @examples
el2regulon = function(el) {
  regulon_list = split(el, el$from)
  viper_regulons = lapply(regulon_list, function(regulon) {
    tfmode = stats::setNames(regulon$weight, regulon$to)
    list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
  })

  return(viper_regulons)
}

#' Convert bipartite adjacency to regulon
#'
#' @param adj An adjacency matrix, with rows as TFs and columns as genes.
#'
#' @return A VIPER required regulon object.
#' @export
#'
#' @examples
adj2regulon = function(adj){
  el = adj2el(adj)
  el = el[el[,3]!=0,]
  regulon = el2regulon(el)
  return(regulon)
}


#' Filter low confident edge signs in the prior network using GeneNet
#'
#' @param prior A prior network (adjacency matrix) with rows as TFs and columns as genes.
#' @param expr A normalized log-transformed gene expression matrix.
#'
#' @return A filtered prior network (adjacency matrix).
#' @export
#'
#' @examples
prior.pp = function(prior,expr){

  # filter tfs and tgs
  tf = intersect(rownames(prior),rownames(expr)) ## TF needs to express
  tg = intersect(colnames(prior),rownames(expr))
  all.gene = unique(c(tf,tg))

  # create coexp net
  coexp = GeneNet::ggm.estimate.pcor(t(expr[all.gene,]), method = "static")
  diag(coexp)= 0

  # prior and coexp nets
  P_ij = prior[tf,tg] ## prior ij
  C_ij = coexp[tf,tg]*abs(P_ij) ## coexpression ij

  # signs
  sign_P = sign(P_ij) ## signs in prior
  sign_C = sign(C_ij) ## signs in coexp

  # blurred edge index
  blurs = which((sign_P*sign_C)<0,arr.ind = T) ## inconsistent edges
  P_ij[blurs] = 1e-6

  # remove all zero TFs (in case prior has all zero TFs)
  A_ij = P_ij
  A_ij = A_ij[rowSums(A_ij!=0)>0,]
  A_ij = A_ij[,colSums(A_ij!=0)>0]

  return(A_ij)
}



# stan model, conditional likelihood
TIGER_C <-
'
data {
  int<lower=0> n_genes;                       // Number of genes
  int<lower=0> n_samples;                     // Number of samples
  int<lower=0> n_TFs;                         // Number of TFs
  int<lower=0> n_zero;                        // length of zero elements in P
  int<lower=0> n_ones;                        // length of non-zero elements in P
  int<lower=0> n_negs;                        // length of repression elements in P
  int<lower=0> n_poss;                        // length of activation elements in P
  int<lower=0> n_blur;                        // length of blurred elements in P
  int<lower=0> n_all;                         // length of all elements in P
  matrix[n_genes,n_samples] X;                // Gene expression matrix X
  vector[n_all] P;                            // Prior connection probability
  array[n_zero] int P_zero;                   // index of zero probablity edges
  array[n_ones] int P_ones;                   // index of non-zero prob edges
  array[n_negs] int P_negs;                   // index of repression prob edges
  array[n_poss] int P_poss;                   // index of activation prob edges
  array[n_blur] int P_blur;                   // index of blurred prob edges
  int sign;                                   // use signed prior network or not
  int baseline;                               // inclue baseline term or not
  int psis_loo;                               // use loo to check model or not
  real sigmaZ;                                // prior sd of Z
  real sigmaB;                                // prior sd of baseline
  real a_alpha;                               // hyparameter for inv_gamma
  real b_alpha;
  real a_sigma;                               // hyparameter for inv_gamma
  real b_sigma;
}

transformed data {
  vector[n_genes*n_samples] X_vec;            // gene expression X
  X_vec = to_vector(X);
}

parameters {
  matrix<lower=0>[n_TFs,n_samples] Z;         // TF activity matrix Z
  vector<lower=0>[n_genes] sigma2;            // variances of noise term
  vector[baseline ? n_genes : 0] b0;          // baseline expression for each gene
  vector<lower=0>[sign ? n_blur : 0] alpha0;  // Noise precision of W_blur
  vector<lower=0>[sign ? n_poss : 0] alpha2;  // Noise precision of W_poss
  vector<lower=0>[sign ? n_negs : 0] alpha3;  // Noise precision of W_negs
  vector[sign ? n_blur : 0] beta0;            // Regulatory network blurred edge weight
  vector<upper=0>[sign ? n_negs : 0] beta3;   // Regulatory network negative edge weight
  vector<lower=0>[sign ? n_poss : 0] beta2;   // Regulatory network positive edge weight
  vector<lower=0>[sign ? 0 : n_ones] alpha1;  // Noise precision of W_ones
  vector[sign ? 0 : n_ones] beta1;            // Regulatory network non-zero edge weight
}

transformed parameters {
  matrix[n_genes, n_TFs] W;                   // Regulatory netwrok W
  vector[n_all] W_vec;                        // Regulatory vector W_vec

  W_vec[P_zero]=rep_vector(0,n_zero);
  if (sign) {
    vector[n_negs] W_negs = beta3.*sqrt(alpha3);  // Regulatory network negative edge weight
    vector[n_poss] W_poss = beta2.*sqrt(alpha2);  // Regulatory network positive edge weight
    vector[n_blur] W_blur = beta0.*sqrt(alpha0);  // Regulatory network blurred edge weight
    W_vec[P_negs]=W_negs;
    W_vec[P_poss]=W_poss;
    W_vec[P_blur]=W_blur;
  }else{
    vector[n_ones] W_ones = beta1.*sqrt(alpha1);  // Regulatory network non-zero edge weight
    W_vec[P_ones]=W_ones;
  }
  W = to_matrix(W_vec,n_genes,n_TFs); // by column

  matrix[n_genes,n_samples] mu = W*Z; // mu for gene expression X
  if (baseline){
    matrix[n_genes,n_samples] mu0;    // baseline
    mu0 = rep_matrix(b0,n_samples);
    mu =  mu + mu0;
  }

  vector[n_genes*n_samples] X_mu = to_vector(mu);
  vector[n_genes*n_samples] X_sigma = to_vector(rep_matrix(sqrt(sigma2),n_samples));
}

model {
  // priors
  sigma2 ~ inv_gamma(a_sigma,b_sigma);

  if (baseline){
    b0 ~ normal(0,sigmaB);
  }

  if (sign) {
    // student-t
    alpha2 ~ inv_gamma(a_alpha,b_alpha);
    beta2 ~ normal(0,1);

    alpha3 ~ inv_gamma(a_alpha,b_alpha);
    beta3 ~ normal(0,1);

    alpha0 ~ inv_gamma(a_alpha,b_alpha);
    beta0 ~ normal(0,1);

  }else{
    alpha1 ~ inv_gamma(a_alpha,b_alpha);
    beta1 ~ normal(0,1);
  }

  to_vector(Z) ~ normal(0,sigmaZ);

  // likelihood
  X_vec ~ normal(X_mu, X_sigma);

}


generated quantities {
  vector[n_genes*n_samples] log_lik;
  if (psis_loo){
    // leave one element out
    for (i in 1:n_genes*n_samples){
      log_lik[i] = normal_lpdf(X_vec[i]|X_mu[i],X_sigma[i]);
    }
  }
}

'
