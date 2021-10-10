import numpy as np

from numpy.linalg import inv

import math
#import cmath
import scipy
#from scipy.linalg import null_space

from scipy.optimize import minimize
from scipy.linalg import pinvh

# ----------------------------------------------------------------------
def solve_Ax_b(A, b):
  sz = len(b)
  x = np.zeros((sz))

  sz_b = 0
  b_use = [False] * sz
  for i in range( sz ):
    if( abs(b[i]) > 1.e-15 ):
      b_use[i] = True
      sz_b += 1
  
  sz_A = 0
  A_use = [False] * sz
  for r in range( sz ):
    if not b_use[r]:
      continue
    for c in range( sz ):
      if( abs(A[r,c]) > 1.e-15 ):
        A_use[c] = True
        sz_A += 1
  
  b_new = np.zeros((sz_b))
  A_new = np.zeros((sz_b,sz_A))
  r = -1
  for i in range( sz ):
    if not b_use[i]:
      continue
    r += 1
    b_new[r] = b[i]
    
    c = -1
    for j in range( sz ):
      if not A_use[j]:
        continue
      c += 1
      A_new[r,c] = A[i,j]

  zct = np.dot(b_new, A_new)

  def fun(vctr):
    return np.linalg.norm(np.dot(A_new,vctr)-b_new)**2.0

  def dfun(vctr):
    wct = np.dot(A_new,vctr)
    wct = np.dot(A_new.T,wct)
    return 2.0*(wct-zct)

  x_new = np.zeros((sz_A))
  x_new = scipy.optimize.minimize(fun, x_new, method='Newton-CG', \
                                  jac=dfun, tol=1.e-9).x
  for i in range(len(x_new)):
    if abs(x_new[i]) <= 1.e-8:
      x_new[i] = 0.0

  j = -1
  for i in range(sz):
    if not A_use[i]:
      continue
    j += 1
    x[i] = x_new[j]

  return x
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
def solve_Ax_b_naive(A, b):
  zct = np.dot(b, A)

  def fun(vctr):
    return np.linalg.norm(np.dot(A,vctr)-b)**2.0

  def dfun(vctr):
    wct = np.dot(A,vctr)
    wct = np.dot(A.T,wct)
    return 2.0*(wct-zct)

  #x = np.zeros(( len(b) ))
  x = np.random.rand(len(b))
  x = scipy.optimize.minimize(fun, x, method='Newton-CG', \
                                  jac=dfun, tol=1.e-12).x
  for i in range(len(b)):
    if abs(x[i]) <= 1.e-8:
      x[i] = 0.0

  return x
# ----------------------------------------------------------------------

  
# ----------------------------------------------------------------------
def find_a_vec_c(psi, indx_h, HSgm, SgmSgm, nq):
# psi is a vector
  Nbasis = SgmSgm.shape[2]

  a_vec = np.zeros((Nbasis))
  b_vec = np.zeros((Nbasis))
  A_mtrx = np.zeros((Nbasis,Nbasis))

  psi_l = psi.reshape((1,2**nq))
  psi_r = psi.reshape((2**nq,1))
  
# b vector calculation  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
  for i in range(Nbasis):
    val = np.dot(np.conj(psi_l), np.dot(HSgm[:,:,indx_h,i],psi_r)).item()
    if abs(val.imag) > 1.e-15:
      print("Strange in b!!! Imaginary part is not zero!!!", val)
    b_vec[i] = val.real
    
    #if abs(b_vec[i]) > 1.e-15:
      #print(i, b_vec[i])
  if np.amax(abs(b_vec)) == 0.0:
    return a_vec


# A matrix calculation  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
  for i in range(Nbasis):
    ##if i%100 == 0:
      ##print(i, Nbasis)
    for j in range(i,Nbasis):
      val = np.dot(np.conj(psi_l), np.dot(SgmSgm[:,:,i,j],psi_r)).item()
      if abs(val.imag) > 1.e-15:
        print("Strange in A_mtrx!!! Imaginary part is not zero!!!", val)

      A_mtrx[i,j] = val.real
      A_mtrx[j,i] = val.real
      #if abs(val.real) > 1.e-15 and i != j:
        #print(i,j,A_mtrx[i,j])

# Solve the system of equation  -  -  -  -  -  -  -  -  -  -  -  -  -  - 

  a_vec = solve_Ax_b(A_mtrx, b_vec)

  return a_vec
# ----------------------------------------------------------------------
