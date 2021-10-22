import numpy as np
import ot
import ot.plot

def cal_bary_wass( x, M=None,w=1, reg=2e-2 ): 
  A = np.vstack(x)                    
  n_distributions = A.shape[1]  
  if (M is None):
    M = ot.utils.dist0(A.shape[0])                      
    M /= M.max() 
  if w==1:
    alpha = w/n_distributions
    weights = np.repeat(alpha, n_distributions)
  if w!=1:
    weights = np.array(w)
    weights = weights / sum(weights) #require sum(weights)=1
  
  bary_wass = ot.bregman.barycenter(A, M, reg, weights)
  if(bary_wass[0]!=bary_wass[0] or max(bary_wass)<np.quantile(A, 0.3)):
    reg = 2e-3
    bary_wass = ot.bregman.barycenter(A, M, reg, weights)
  if(bary_wass[0]!=bary_wass[0] or max(bary_wass)<np.quantile(A, 0.3)):
    reg = 2e-4
    bary_wass = ot.bregman.barycenter(A, M, reg, weights)
  return bary_wass
