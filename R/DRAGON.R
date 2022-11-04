# Function to implement DRAGON: https://arxiv.org/abs/2104.01690

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
  x = matrix(c(1,2,3,1,5,12),nrow=3,byrow=T)
  # x is an n x p matrix of data
  # xbar = np.mean(X, 0)
  # n = X.shape[0]
  n = nrow(x)
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
  xbar = apply(x,2,mean)
  p = ncol(x)
  w = array(dim=c(n,p,p))
  for(k in 1:n)
  {
    for(i in 1:p)
    {
      for(j in 1:p) # this will be symmetric, but leaving it now for clarity
      {
        w[k,i,j] = (x[k,i] - xbar[i])*(x[k,j]-xbar[j])
      }
    }
  }
  
  wbar = 1/n*apply(w,2:3,sum)
  summand = 0
  for(k in 1:n)
    summand = summand + (w[k,,]-wbar)^2

  varhat_s = n/((n-1)^3)*summand
  return(varhat_s)
  # this code matches with the python output 
}

# def EsqS(X):
#   #xbar = np.mean(X, 0)
#   n = X.shape[0]
# #x_minus_xbar = X - xbar
# wbar = np.cov(X.T, bias=True) #x_minus_xbar.T@x_minus_xbar/n
# ES2 = wbar**2*n**2/(n-1)**2
# return(ES2)

EsqS = function(x)
{
  # following E[x^2] = Var(x) - E[x]^2
  return(VarS(x) - mean(x)^2)
}

#     def risk(lam):
#R = const + ((1.-lam[0]**2)*T1_1 + (1.-lam[1]**2)*T1_2
#             + (1.-lam[0]**2)**2*T2_1 + (1.-lam[1]**2)**2*T2_2
#             + (1.-lam[0]**2)*(1.-lam[1]**2)*T3 + lam[0]*lam[1]*T4)     #reparametrize lamx = 1-lamx**2 for better optimization
# return(R)

risk = function(lambda1, lambda2, t11, t12, t21, t22, t3, t4)
{
  return(lambda1*t11 + lambda2*t12 + 
           lambda1^2*t21 + lambda2^2*t22 + 
           lambda1*lambda2*t3 + 
           sqrt(1-lambda1)*sqrt(1-lambda2)*t4)
}

# def estimate_penalty_parameters_dragon(X1, X2):
estimatePenaltyParameters = function(X1,X2)
{
  
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
  
  varX = VarS(X)
  esqS = EsqS(X)
  
  # IDs = np.cumsum([p1,p2])
  IDs = cumsum(c(1:p1,1:p2))
  
  # varS1 = varS[0:IDs[0],0:IDs[0]]
  # varS12 = varS[0:IDs[0],IDs[0]:IDs[1]]
  # varS2 = varS[IDs[0]:IDs[1],IDs[0]:IDs[1]]
  # 
  # eSqs1 = eSqs[0:IDs[0],0:IDs[0]]
  # eSqs12 = eSqs[0:IDs[0],IDs[0]:IDs[1]]
  # eSqs2 = eSqs[IDs[0]:IDs[1],IDs[0]:IDs[1]]
  # 
  # const = (np.sum(varS1) + np.sum(varS2) - 2.*np.sum(varS12)
  #          + 4.*np.sum(eSqs12))
  # T1_1 = -2.*(np.sum(varS1) - np.trace(varS1) + np.sum(eSqs12))
  # T1_2 = -2.*(np.sum(varS2) - np.trace(varS2) + np.sum(eSqs12))
  # T2_1 = np.sum(eSqs1) - np.trace(eSqs1)
  # T2_2 = np.sum(eSqs2) - np.trace(eSqs2)
  # T3 = 2.*np.sum(eSqs12)
  # T4 = 4.*(np.sum(varS12)-np.sum(eSqs12))
  # 
  # def risk(lam):
  #   R = const + ((1.-lam[0]**2)*T1_1 + (1.-lam[1]**2)*T1_2
  #                + (1.-lam[0]**2)**2*T2_1 + (1.-lam[1]**2)**2*T2_2
  #                + (1.-lam[0]**2)*(1.-lam[1]**2)*T3 + lam[0]*lam[1]*T4)     #reparametrize lamx = 1-lamx**2 for better optimization
  # return(R)
  # 
  # x = np.arange(0., 1.01, 0.01)
  # lamgrid = meshgrid(x, x)
  # risk_grid = risk(lamgrid)
  # indices = np.unravel_index(np.argmin(risk_grid.T, axis=None), risk_grid.shape)
  # lams = [x[indices[0]],x[indices[1]]]
  # 
  # res = minimize(risk, lams, method='L-BFGS-B',#'TNC',#'SLSQP',
  #                tol=1e-12,
  #                bounds = [[0.,1.],[0.,1.]])
  # 
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

dragon = function()
{

}
