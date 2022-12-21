# helper functions as in the netzooPy impl
# def Scale(X):
#   X_temp = X
#   X_std = np.std(X_temp, axis=0)
#   X_mean = np.mean(X_temp, axis=0)
#   return (X_temp - X_mean) / X_std

scale = function(x,bias=F)
{
  # sd does 1/(n-1), python does 1/n
  # use the bias option for exact match with python
  n = length(x)
  if(bias)
  {
    n = length(x)
    numer = (x-mean(x))
    denom = sqrt((n-1)/n)*sd(x)
    return(numer/denom)
  }
  return((x-mean(x))/sd(x))
}

# def VarS(X):
#   xbar = np.mean(X, 0)
# n = X.shape[0]
# x_minus_xbar = X - xbar
# a = x_minus_xbar*x_minus_xbar
# wbar = np.cov(X.T, bias=True)#x_minus_xbar.T@x_minus_xbar/n
# varS = a.T@a
# varS += - n*wbar**2
# varS *= n/(n-1)**3
# return(varS)

# VarS calculates the unbiased estimate of the entries of 
# S according to the formula in Appendix A of Schafer and Strimmer 
# 2005

VarS = function(x)
{ 
  # x is an n x p matrix of data
  # xbar = np.mean(X, 0)
  # n = X.shape[0]
  n = nrow(x)
  p = ncol(x)
  # x_minus_xbar = X - xbar
  x_minus_xbar = apply(x,2,function(x){x-mean(x)}) # center
  # sanity check
  apply(x_minus_xbar,2,mean)
  a = x_minus_xbar*x_minus_xbar
  varHat_s = n/((n-1)^3)* apply(a,2,sum)
  wbar = 1/n*t(x_minus_xbar) %*% x_minus_xbar # this should be an outer product
  
    # varS = a.T@a
  # varS += - n*wbar**2
  # varS *= n/(n-1)**3
  
  # try doing this with an array instead of matrix multiplication
  # xbar = apply(x,2,mean)
  # p = ncol(x)
  # w = array(dim=c(n,p,p))
  # for(k in 1:n)
  # {
  #   for(i in 1:p)
  #   {
  #     for(j in 1:p) # this will be symmetric, but leaving it now for clarity
  #     {
  #       w[k,i,j] = (x[k,i] - xbar[i])*(x[k,j]-xbar[j])
  #     }
  #   }
  # }
  # 
  # wbar = 1/n*apply(w,2:3,sum)
  # summand = 0
  # for(k in 1:n)
  #   summand = summand + (w[k,,]-wbar)^2
  # 
  # varhat_s = n/((n-1)^3)*summand
  varS = t(a) %*% a + -n*wbar^2
  varHat_s = n/((n-1)^3)*varS
  return(varHat_s)
  # this code matches with the python output 
}

# def EsqS(X):
#   #xbar = np.mean(X, 0)
#   n = X.shape[0]
# #x_minus_xbar = X - xbar
# wbar = np.cov(X.T, bias=True) #x_minus_xbar.T@x_minus_xbar/n
# ES2 = wbar**2*n**2/(n-1)**2 #
# return(ES2)

EsqS = function(x)
{
  n = nrow(x)
  wbar = (n-1)/n*cov(x)
  ES2 = wbar^2*n^2/((n-1)^2)
  return(ES2)
  # this code matches with the python output
}

# def risk_orig(lam):
#   R = const + (lam[0]*T1_1 + lam[1]*T1_2
#                + lam[0]**2*T2_1 + lam[1]**2*T2_2
#                + lam[0]*lam[1]*T3 + np.sqrt(1-lam[0])*np.sqrt(1-lam[1])*T4)     #reparametrize lamx = 1-lamx**2 for better optimization
# return(R)

risk_orig = function(lambda1, lambda2, t11, t12, t21, t22, t3, t4)
{
  print("[dragonR] risk_orig(): This is not the reparameterized version.")
  
  return(lambda1*t11 + lambda2*t12 + 
           lambda1^2*t21 + lambda2^2*t22 + 
           lambda1*lambda2*t3 + 
           sqrt(1-lambda1)*sqrt(1-lambda2)*t4)
}

# def risk(lam):
#   R = const + ((1.-lam[0]**2)*T1_1 + (1.-lam[1]**2)*T1_2
#                + (1.-lam[0]**2)**2*T2_1 + (1.-lam[1]**2)**2*T2_2
#                + (1.-lam[0]**2)*(1.-lam[1]**2)*T3 + lam[0]*lam[1]*T4)     #reparametrize lamx = 1-lamx**2 for better optimization
# return(R)

# reparameterize gamma1 = (1-lambda1^2), gamma2 = (1-lambda2^2)
risk = function(gamma, const, t11, t12, t21, t22, t3, t4)
{
  gamma1 = gamma[1]
  gamma2 = gamma[2]
  
  R = const + (1-gamma1^2)*t11 + (1-gamma2^2)*t12 +
    (1-gamma1^2)^2*t21 + (1-gamma2^2)^2*t22 +
     (1-gamma1^2)*(1-gamma2^2)*t3 + gamma1*gamma2*t4 
  return(R)
}

# def estimate_penalty_parameters_dragon(X1, X2):
estimatePenaltyParameters = function(X1,X2)
{
  # X1 = matrix(c(1,2,3,1,5,12),nrow=3,byrow=T)
  # X2 = matrix(c(9,7,8),nrow=3,byrow=T)
  # X1 is omics matrix 1, dimensions n x p1 
  # X2 is omics matrix 2, dimensions n x p2
  # The matrices should have the same ordering by samples
  n = nrow(X1)
  p1 = ncol(X1)
  p2 = ncol(X2)
  # n = X1.shape[0]
  # p1 = X1.shape[1]
  # p2 = X2.shape[1]
  X = cbind.data.frame(X1, X2)
  # X = np.append(X1, X2, axis=1)
  
  # varS = VarS(X)
  # eSqs = EsqS(X)
  
  varS = VarS(X)
  esqS = EsqS(X)
  # IDs = np.cumsum([p1,p2])
  IDs = cumsum(c(p1,p2))
  
  # varS1 = varS[0:IDs[0],0:IDs[0]]
  varS1 = varS[1:IDs[1],1:IDs[1]]
  # varS12 = varS[0:IDs[0],IDs[0]:IDs[1]]
  varS12 = varS[1:IDs[1],(IDs[1]+1):IDs[2]]
  # varS2 = varS[IDs[0]:IDs[1],IDs[0]:IDs[1]]
  varS2 = varS[(IDs[1]+1):IDs[2],(IDs[1]+1):IDs[2]]
  # eSqs1 = eSqs[0:IDs[0],0:IDs[0]]
  esqS1 = esqS[1:IDs[1],1:IDs[1]]
  # eSqs12 = eSqs[0:IDs[0],IDs[0]:IDs[1]]
  esqS12 = esqS[1:IDs[1],(IDs[1]+1):IDs[2]]
  # eSqs2 = eSqs[IDs[0]:IDs[1],IDs[0]:IDs[1]]
  esqS2 =  esqS[(IDs[1]+1):IDs[2],(IDs[1]+1):IDs[2]]
  # 
  # const = (np.sum(varS1) + np.sum(varS2) - 2.*np.sum(varS12)
  #          + 4.*np.sum(eSqs12))

  # T1_1 = -2.*(np.sum(varS1) - np.trace(varS1) + np.sum(eSqs12))
  T1_1 = -2*(sum(varS1) - matrix.trace(as.matrix(varS1)) + sum(esqS12))
  # T1_2 = -2.*(np.sum(varS2) - np.trace(varS2) + np.sum(eSqs12))
  T1_2 = -2*(sum(varS2) - matrix.trace(as.matrix(varS2)) + sum(esqS12))
  # T2_1 = np.sum(eSqs1) - np.trace(eSqs1)
  T2_1 = sum(esqS1) - matrix.trace(as.matrix(esqS1))
  # T2_2 = np.sum(eSqs2) - np.trace(eSqs2)
  T2_2 = sum(esqS2) - matrix.trace(as.matrix(esqS2))
  # T3 = 2.*np.sum(eSqs12)
  T3 = 2*sum(esqS12)
  # T4 = 4.*(np.sum(varS12)-np.sum(eSqs12))
  T4 = 4*(sum(varS12)-sum(esqS12))
  
  const = (sum(varS1) + sum(varS2) - 2*sum(varS12)
           + 4*sum(esqS12))
  
  # x = np.arange(0., 1.01, 0.01)
  x = seq(0,1,by=0.01)
  riskgrid = matrix(nrow = length(x),ncol=length(x))
  for(i in 1:length(x))
  {
    for(j in 1:length(x))
    {
      riskgrid[i,j] = risk(gamma=c(x[i],x[j]),
                           const = const,
                           t11=T1_1,
                           t12=T1_2,
                           t21=T2_1,
                           t22=T2_2,
                           t3=T3,
                           t4=T4)
    }
  }
  
  dim(riskgrid)
  # lamgrid = meshgrid(x, x)
  # risk_grid = risk(lamgrid)
  # indices = np.unravel_index(np.argmin(risk_grid.T, axis=None), risk_grid.shape)
  # lams = [x[indices[0]],x[indices[1]]]
  # 
  lams = x[arrayInd(which(riskgrid == min(riskgrid)),.dim=c(101,101))]
  # this is seeding with grid search and then we use optimization
  print(lams)
  
  # res = minimize(risk, lams, method='L-BFGS-B',#'TNC',#'SLSQP',
  #                tol=1e-12,
  #                bounds = [[0.,1.],[0.,1.]])
  
  # python result:  array([0.79372539, 0.        ])
  res = optim(lams,risk, 
              const=const,
              t11=T1_1,
              t12=T1_2,
              t21=T2_1,
              t22=T2_2,
              t3=T3,
              t4=T4,
              method="L-BFGS-B",
              lower=c(0,0),
              upper=c(1,1),
              control = list(trace=T,pgtol = 1e-15))
  
  # reparameterize
  lambdas = c(1-res$par[1]^2, 1-res$par[2]^2)
  return(list("lambdas"=lambdas,"gammas"=res$par,"optim_result"=res,"risk_grid" = riskgrid))
  # penalty_parameters = (1.-res.x[0]**2), (1.-res.x[1]**2)
  # 
  # def risk_orig(lam):
  #   R = const + (lam[0]*T1_1 + lam[1]*T1_2
  #                + lam[0]**2*T2_1 + lam[1]**2*T2_2
  #                + lam[0]*lam[1]*T3 + np.sqrt(1-lam[0])*np.sqrt(1-lam[1])*T4)     #reparametrize lamx = 1-lamx**2 for better optimization
  # return(R)
  # risk_grid_orig = risk_orig(lamgrid)
  # return(penalty_parameters, risk_grid_orig)
  # 
}

get_shrunken_covariance_dragon = function(X1,X2, lambdas)
{
  n = nrow(X1)
  p1 = ncol(X1)
  p2 = ncol(X2)
  p = p1 + p2
  X = cbind.data.frame(X1,X2)
  S = cov(X) # the R implementation of cov() uses the unbiased (1/(n-1)); we need the unbiased version for the lemma of Ledoit and Wolf
  
  # target matrix
  Targ = diag(diag(S))
  Sigma = matrix(nrow=p, ncol=p)
  
  # Sigma = np.zeros((p,p))
  # IDs = np.cumsum([0,p1,p2])
  IDs = c(cumsum(c(p1,p2)))
  
  idx1 = 1:IDs[1]
  idx2 = (IDs[1]+1):IDs[2]
  
  # Fill in Sigma_11
  Sigma[idx1,idx1] = (1-lambdas[1])*S[idx1,idx1] + lambdas[1]*Targ[idx1,idx1]
  
  # Fill in Sigma_22
  Sigma[idx2,idx2] = (1-lambdas[2])*S[idx2,idx2] + lambdas[2]*Targ[idx2,idx2]
  
  # Fill in Sigma_12 
  Sigma[idx1,idx2] = sqrt((1-lambdas[1])*(1-lambdas[2]))*S[idx1,idx2] + sqrt(lambdas[1]*lambdas[2])*Targ[idx1,idx2]
  
  # Fill in Sigma_21
  Sigma[idx2,idx1] = sqrt((1-lambdas[1])*(1-lambdas[2]))*S[idx2,idx1] + sqrt(lambdas[1]*lambdas[2])*Targ[idx2,idx1]
  
  return(Sigma)
}

get_precision_matrix_dragon = function(X1, X2, lambdas)
{
   Sigma = get_shrunken_covariance_dragon(X1, X2, lambdas)
   Theta = solve(Sigma)
   return(Theta)
   
   # in the python implementation, mean is also returned. Omitting here
  #  X = np.hstack((X1, X2))
  #  mu = np.mean(X, axis=0)
}

get_partial_correlation_from_precision = function(Theta,selfEdges=F)
{
  # by default, does not return self edges (diagonal is set to zero)
  ggm = -cov2cor(Theta)
  if(!selfEdges)
    ggm[diag(ggm)] = 0
  return(ggm)
}

get_partial_correlation_dragon = function(X1,X2,lambdas)
{
  Theta = get_precision_matrix_dragon(X1, X2, lambdas)
  ggm = get_partial_correlation_from_precision(Theta)
  return(ggm)
}

# The functions below are for p-value estimation on the DRAGON results

# The functions logli, estimate_kappa, and estimate_p_values are for benchmarking
# with standard GGM; omitting here
# logli = function(X, Theta, mu)
# estimate_kappa = function(n, p, lambda0, seed)
# estimate_p_values(r, n, lambda0, kappa='estimate', seed=1): 

log_lik_shrunken = function(kappa, p, lambda, rhos)
{
  # kappa is to be optimized, so comes first in the arguments
  # p is fixed (number of predictors)
  # lambda is fixed (as estimated by DRAGON)
  # rhos is fixed (observed partial correlations from the data)
  mysum = 0
  
  for(i in 1:p)
  {
    first_term = (kappa-3)/2*log((1-lambda)^2-rhos[i]^2)
    second_term = lbeta(1/2, (kappa-1)/2)
    third_term = (kappa-2)*log(1-lambda)
    mysum = mysum + first_term - second_term - third_term
  }
  
  return(mysum)
}
                   
# def estimate_kappa_dragon(n, p1, p2, lambdas, seed, simultaneous = False):
estimate_kappa_dragon = function(n, p1, n2, lambdas, seed, simultaneous = F)
{
  
}

# def estimate_p_values_dragon(r, n, p1, p2, lambdas, kappa='estimate', seed=1, simultaneous = False):
estimate_p_values_dragon = function(r, n, p1, p2, lambdas, kappa="estimate",seed=1, simultaneous = F)
{
  
}

#' Run DRAGON in R.
#' 
#' Description: Estimates a multi-omic Gaussian graphical model for two input layers of paired omic data.
#'
#' @param layer1 : first layer of omics data; rows: samples (order must match layer2), columns: variables
#' @param layer2 : second layer of omics data; rows: samples (order must match layer1), columns: variables.
#' @param pval : calculate p-values for network edges. Not yet implemented in R; available in netZooPy.
#' @param gradient : method for estimating parameters of p-value distribution, applies only if p-val == T. default = "finite_difference"; other option = "exact"
#' @return A list of model results. cov : the shrunken covariance matrix
#' \itemize{
#'  \item{\code{cov}}{  the shrunken covariance matrix}
#'  \item{\code{prec}}{  the shrunken precision matrix}
#'  \item{\code{ggm}}{ the shrunken Gaussian graphical model; matrix of partial correlations. Self-edges (diagonal elements) are set to zero.}
#'  \item{\code{lambdas}}{  Vector of omics-specific tuning parameters (lambda1, lambda2) for \code{layer1} and \code{layer2}}
#'  \item{\code{gammas}}{  Reparameterized tuning parameters; gamma = 1 - lambda^2}
#'  \item{\code{risk_grid}}{  Risk grid, for assessing optimization. Grid boundaries are in terms of gamma.}
#' }
#' 
#' @export
dragon = function(layer1,layer2,pval = F,gradient = "finite_difference", verbose = F)
{
  if(verbose)
    print("[netZooR::dragon] Estimating penalty parameters...")
  # estimate penalty parameters
  myres = estimatePenaltyParameters(layer1, layer2)
  lambdas = myres$lambdas
  
  if(verbose)
    print(paste(c("[netZooR::dragon] Estimated parameters:",lambdas),collapse=" "))
  
  if(verbose)
    print("[netZooR::dragon] Calculating shrunken matrices...")
  # apply penalty parameters to return shrunken covariance and ggm
  shrunken_cov = get_shrunken_covariance_dragon(layer1, layer2,lambdas)
  precmat = get_precision_matrix_dragon(layer1, layer2, lambdas)
  ggm = get_partial_correlation_dragon(layer1, layer2, lambdas)
  
  # if pval, return pval approx with finite difference
  if(pval)
  {
    print("[netZooR::dragon] p-value calculation not yet implemented in R; to estimate p-values, use netZooPy.")
  }
  
  return(list("cov"=shrunken_cov,
              "prec"=precmat,
              "ggm"=ggm,
              "lambdas"=lambdas,
              "gammas"=myres$gammas,
              "risk_grid"=myres$risk_grid))
}
