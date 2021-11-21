import numpy as np

from numpy.linalg import inv

import math
#import cmath
import scipy
#from scipy.linalg import null_space

from scipy.optimize import minimize
from scipy.linalg import pinvh

# arxiv:1608.04571

# -------------------------------------------------------------------
def solve_Ax_b_L_curve(A, b, lmbd_min=1.e-4, lmbd_max=1.e-2):
  sz = len(b)

  phi = 0.5*(1 + math.sqrt(5))

  lmbd = np.zeros(4)
  x = np.zeros(4)
  
  def menger(Pj, Pk, Pl):
    den = np.linalg.norm(Pj - Pk) \
        * np.linalg.norm(Pk - Pl) \
        * np.linalg.norm(Pl - Pj)

    num = 2 * (  Pj[0]*Pk[1] + Pk[0]*Pl[1] + Pl[0]*Pj[1] 
               - Pj[0]*Pl[1] - Pk[0]*Pj[1] - Pl[0]*Pk[1])

    #print("den = ", den)
    return num/den 

  lmbd[0] = lmbd_min # extremum value
  x[0] = math.log10(lmbd[0])
  
  lmbd[3] = lmbd_max # extremum value
  x[3] = math.log10(lmbd[3])

  x[1] = (x[3] + phi*x[0]) / (1 + phi)
  lmbd[1] = 10**x[1]

  x[2] = x[0] + (x[3] - x[1])
  lmbd[2] = 10**x[2]
  
  P = []
  solution = []
  for l in lmbd:
    obj = solve_Ax_b_Tikhonov(A, b, l)
    solution.append( obj[0] )
    P.append( obj[1] )


  while (lmbd[3] - lmbd[0])/lmbd[3] > 1.e-3:
    C1 = menger(P[0], P[1], P[2])
    C2 = menger(P[1], P[2], P[3])

    while C2 < 0:
      lmbd[3] = lmbd[2]
      x[3] = x[2]
      P[3] = P[2]

      lmbd[2] = lmbd[1]
      x[2] = x[1]
      P[2] = P[1]
      
      x[1] = (x[3] + phi*x[0]) / (1 + phi)
      lmbd[1] = 10**x[1]
      solution[1], P[1] = solve_Ax_b_Tikhonov(A, b, lmbd[1])
      #print(P[1])

      C2 = menger(P[1], P[2], P[3])
      
    if C1 > C2:
      l_save = lmbd[1]
      x_save = x[1]
      P_save = P[1]
      solution_save = solution[1]

      lmbd[3] = lmbd[2]
      x[3] = x[2]
      P[3] = P[2]

      lmbd[2] = lmbd[1]
      x[2] = x[1]
      P[2] = P[1]

      x[1] = (x[3] + phi*x[0]) / (1 + phi)
      lmbd[1] = 10**x[1]
      solution[1], P[1] = solve_Ax_b_Tikhonov(A, b, lmbd[1])
      #print(P[1])
    else:
      l_save = lmbd[2]
      x_save = x[2]
      P_save = P[2]
      solution_save = solution[2]

      lmbd[0] = lmbd[1]
      x[0] = x[1]
      P[0] = P[1]

      lmbd[1] = lmbd[2]
      x[1] = x[2]
      P[1] = P[2]

      x[2] = x[0] + (x[3] - x[1])
      lmbd[2] = 10**x[2]
      solution[2], P[2] = solve_Ax_b_Tikhonov(A, b, lmbd[2])
      #print(P[2])

  return solution_save
# -------------------------------------------------------------------


# -------------------------------------------------------------------
def solve_Ax_b_Tikhonov(A, b, lmbd):
  sz = len(b)
  
  bA = np.dot(b, A)

  def fun(x):
    return np.linalg.norm(np.dot(A,x)-b)**2.0 \
        + lmbd * np.linalg.norm(x)**2.0

  def dfun(x):
    return 2.0*(np.dot(A.T, np.dot(A,x) ) - bA + lmbd*x)

  x = np.zeros((sz))
  x = scipy.optimize.minimize(fun, x, method='Newton-CG',
                              jac=dfun, tol=1.e-9).x

  xi = math.log( np.linalg.norm(np.dot(A,x)-b)**2.0 )
  eta = math.log( np.linalg.norm(x)**2.0 )
  return x, np.array([xi,eta])
# -------------------------------------------------------------------


# -------------------------------------------------------------------
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
  obj = scipy.optimize.minimize(fun, x_new, 
                                  method='Newton-CG', \
                                  jac=dfun, tol=1.e-9)
  x_new = obj.x
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
# -------------------------------------------------------------------


# -------------------------------------------------------------------
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
# -------------------------------------------------------------------

  
# -------------------------------------------------------------------
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
