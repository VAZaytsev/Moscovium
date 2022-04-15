import numpy as np

kp_to_name = {
  -1: "s",
   1: "p1", 
  -2: "p3",
   2: "d3",
  -3: "d5"
  }


# ===================================================================
def bin_to_vec(bn):
  wf_bn = 1

  for b in bn:
    if b == "0":
      wf_bn = np.kron( wf_bn, np.array([1,0]) )
    else:
      wf_bn = np.kron( wf_bn, np.array([0,1]) )

  return wf_bn
# ===================================================================


# ===================================================================
def get_Nsteps_for_exp_tauH(psi, H_mtrx, exp_tau_H):
  psi_in = psi.copy()
  Niter_max = 10000

  converged = False
  E_old = (psi.conj().T).dot( H_mtrx.dot(psi) ).item().real

  for it in range(Niter_max):
    psi_in = exp_tau_H.dot(psi_in)
    psi_in /= np.linalg.norm(psi_in)


    E_new = (psi_in.conj().T).dot( H_mtrx.dot(psi_in) ).item().real

    if abs(E_new - E_old) < 1.e-6:
      converged = True
      return it, converged, E_new

    E_old = E_new

  return Niter_max, converged, E_new
# ===================================================================


# ===================================================================
def average_value_for_bn_mtrx(bn, mtrx):
  wf_bn = bin_to_vec(bn)
  return np.conj( wf_bn.T ).dot( mtrx.dot( wf_bn )).item(0)
# ===================================================================


# ===================================================================
def bin_to_conf( orb_arr, bn ):
  sz = len(orb_arr)
  
  n_arr = []
  k_arr = []
  q_arr = []

# Construct an array of radial orbitals  
  for x in orb_arr:
    Found = False
    for i in range(len(n_arr)):
      if x.n == n_arr[i] and x.k == k_arr[i]:
        Found = True
        break

    if not Found:
      n_arr.append(x.n)
      k_arr.append(x.k)
      q_arr.append(0)


# Assign the ocupation numbers to the radial orbitals
  for i in range(sz):
    x = bn[sz-1-i]
    if x == "0":
      continue

    for ii in range(len(n_arr)):
      if orb_arr[i].n == n_arr[ii] and orb_arr[i].k == k_arr[ii]:
        q_arr[ii] += 1
        break

  conf = ""
  for i in range(len(n_arr)):
    if q_arr[i] == 0:
      continue
    conf += str(n_arr[i]) \
    + kp_to_name[k_arr[i]] \
    + "^" + str(q_arr[i]) \
    + " "

  return conf
# ===================================================================
