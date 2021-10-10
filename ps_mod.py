import numpy as np
import cmath

# Dictionaries =========================================================

M_P_MP_PM = {"M": [[1,2], [0.5,-0.5*1j]], \
             "P": [[1,2],[0.5,0.5*1j]], \
             "MP": [[0,3],[0.5,-0.5]], \
             "PM": [[0,3],[0.5,0.5]]
           }

SgmSgm = {
  "12": [ 1j, "3"],
  "13": [-1j, "2"],
  "23": [ 1j, "1"],
  "21": [-1j, "3"],
  "31": [ 1j, "2"],
  "32": [-1j, "1"],
  }

PP_P = np.matrix( [[0,1,2,3],
                   [1,0,3,2],
                   [2,3,0,1],
                   [3,2,1,0]], dtype=int)
PP_C = np.matrix( [[1, 1,  1,  1],
                   [1, 1, 1j,-1j],
                   [1,-1j, 1, 1j],
                   [1, 1j,-1j, 1]], dtype=complex)

IXYZ = {
  0: np.matrix([[1, 0],[0, 1]]),
  1: np.matrix([[0, 1],[1, 0]]),
  2: np.matrix([[0,-1j],[1j, 0]],dtype=complex),
  3: np.matrix([[1, 0],[0,-1]]),
  }
# ===================================================================


# ===================================================================
def ps2mtrx(ps):
  nq = len(ps)
  for i in range(nq):
    if i == 0:
      mtrx = IXYZ[int(ps[i])]
    else:
      mtrx = np.kron(mtrx, IXYZ[int(ps[i])])
  return mtrx
# ===================================================================


# ===================================================================
def ps_2_vec(ps, nq):
  if int(ps[0]) == 1:
    vec = np.array([0, 1])
  else:
    vec = np.array([1, 0])

  for i in range(1,nq):
    if int(ps[i]) == 1:
      vec = np.kron(vec, np.array([0,1]) )
    else:
      vec = np.kron(vec, np.array([1,0]) )

  return vec
# ===================================================================


# ===================================================================
def ps_ps(ps_in1, ps_in2):
  nq = len(ps_in1)
  coef = 1.0
  arr = [""]*nq
  for i in range(nq):
    arr[i] = str(PP_P[int(ps_in1[i]), int(ps_in2[i])])
    coef *= PP_C[int(ps_in1[i]), int(ps_in2[i])]

  #for i in range(nq):
    #char = str(int(ps_in1[i])) + str(int(ps_in2[i]))

    #for k in range(4):
      #char = char.replace("0" + str(k), str(k))
      #char = char.replace(str(k) + "0", str(k))
      #char = char.replace(str(k) + str(k), "0")

      #if char in SgmSgm.keys():
        #coef *= SgmSgm[char][0]
        #char = SgmSgm[char][1]

    #arr[i] = char

  return coef, arr
# ======================================================================


# ======================================================================
def dict_2_vec(dct, nq):
  vec_tot = np.zeros((2**nq))

  nrm = 0.0
  for k in dct.keys():
    if int(k[0]) == 1:
      vec = np.array([0, 1])
    else:
      vec = np.array([1, 0])

    for i in range(1,nq):
      if int(k[i]) == 1:
        vec = np.kron(vec, np.array([0,1]) )
      else:
        vec = np.kron(vec, np.array([1,0]) )
    
    nrm += dct[k]
    vec_tot += vec * dct[k]

  vec_tot = vec_tot / nrm
  return vec_tot
# ======================================================================


# ======================================================================
def pq_2_ps(p, q, nq):
  arr = [""]*nq
  arr[nq-1-p] = arr[nq-1-p] + "M"
  for j in range(p+1,nq):
    arr[nq-1-j] = arr[nq-1-j] + "Z"
    
  arr[nq-1-q] = arr[nq-1-q] + "P"
  for j in range(q+1,nq):
    arr[nq-1-j] = arr[nq-1-j] + "Z"

  coef = 1.0
  for j in range(nq):
    arr[j] = arr[j].replace("ZZ","")
    arr[j] = arr[j].replace("ZP","P")
    arr[j] = arr[j].replace("MZ","M")
    arr[j] = arr[j].replace("PZ","-P")
    if "-" in arr[j]:
      coef = -coef
      arr[j] = arr[j].replace("-","")
    arr[j] = arr[j].replace("ZM","-M")
    if "-" in arr[j]:
      coef = -coef
      arr[j] = arr[j].replace("-","")

  nsplits = 0
  for j in range(nq):
    if arr[j] == "":
      arr[j] = 0
      continue
      
    if arr[j] == "Z":
      arr[j] = 3
      continue
    
    nsplits += 1

  coef_out = [coef] * (2**nsplits)
  arr_out = np.zeros((2**nsplits,nq))
  for j in range(2**nsplits):
    j_bin = bin(j)[2:].zfill(nsplits)
    icntr = 0
    for i in range(nq):
      if arr[i] == 0 or arr[i] == 3:
        arr_out[j,i] = arr[i]
        continue
      arr_out[j,i] = M_P_MP_PM[arr[i]][0][int(j_bin[icntr])]
      coef_out[j] = coef_out[j] * M_P_MP_PM[arr[i]][1][int(j_bin[icntr])]
      icntr += 1

  return coef_out, arr_out
# ======================================================================


# ======================================================================
def pqrs_2_ps(p, q, r, s, nq):
  arr = [""]*nq
  arr[nq-1-p] = arr[nq-1-p] + "M"
  for j in range(p+1,nq):
    arr[nq-1-j] = arr[nq-1-j] + "Z"
    
  arr[nq-1-q] = arr[nq-1-q] + "M"
  for j in range(q+1,nq):
    arr[nq-1-j] = arr[nq-1-j] + "Z"

  arr[nq-1-r] = arr[nq-1-r] + "P"
  for j in range(r+1,nq):
    arr[nq-1-j] = arr[nq-1-j] + "Z"
    
  arr[nq-1-s] = arr[nq-1-s] + "P"
  for j in range(s+1,nq):
    arr[nq-1-j] = arr[nq-1-j] + "Z"

  coef = 1.0
  for j in range(nq):
    arr[j] = arr[j].replace("ZZ","")
    arr[j] = arr[j].replace("ZP","P")
    arr[j] = arr[j].replace("MZ","M")
    arr[j] = arr[j].replace("PZ","-P")
    if "-" in arr[j]:
      coef = -coef
      arr[j] = arr[j].replace("-","")
    arr[j] = arr[j].replace("ZM","-M")
    if "-" in arr[j]:
      coef = -coef
      arr[j] = arr[j].replace("-","")

  nsplits = 0
  for j in range(nq):
    if arr[j] == "":
      arr[j] = 0
      continue
      
    if arr[j] == "Z":
      arr[j] = 3
      continue
    
    nsplits += 1

  coef_out = [coef] * (2**nsplits)
  arr_out = np.zeros((2**nsplits,nq))
  for j in range(2**nsplits):
    j_bin = bin(j)[2:].zfill(nsplits)
    icntr = 0
    for i in range(nq):
      if arr[i] == 0 or arr[i] == 3:
        arr_out[j,i] = arr[i]
        continue
      arr_out[j,i] = M_P_MP_PM[arr[i]][0][int(j_bin[icntr])]
      coef_out[j] = coef_out[j] * M_P_MP_PM[arr[i]][1][int(j_bin[icntr])]
      icntr += 1
   
  return coef_out, arr_out
# ======================================================================


# ======================================================================
def find_ps_in_arr(ps, ps_arr):
  found = False
  for i in range( len(ps_arr) ):
    if all(ps[k] == ps_arr[i][k] for k in range( len(ps) ) ):
      found = True
      break

  if found:
    return i
  else:
    return -1
# ======================================================================


# ======================================================================
def exp_alpha_PS_test(alpha, ps):
  nq = len(ps)
  for i in range(nq):
    if i == 0:
      mtrx_pwr = IXYZ[ps[i]]
    else:
      mtrx_pwr = np.kron(mtrx_pwr, IXYZ[ps[i]])

  eigen_val, eigen_vec = np.linalg.eig(mtrx_pwr)

  mtrx_out = np.zeros((2**nq,2**nq))
  for row in range(2**nq):
    P = np.outer(eigen_vec[:,row], np.conj(eigen_vec[:,row]) )
    mtrx_out = np.add(mtrx_out,cmath.exp(alpha * eigen_val[row]) * P)

  return mtrx_out
# ======================================================================


# ======================================================================
def exp_mtrx(mtrx):
# Check whether matrix is antihermitian
  if np.amax(abs(mtrx + mtrx.conj().T)) <= 1.e-15:
    un_tets = True
    fctr = 1j
    H = -1j*mtrx # Hermitian matrix is constructed
    eigen_val, eigen_vec = np.linalg.eigh(H)
  else:
    un_tets = False
    fctr = 1.0
    eigen_val, eigen_vec = np.linalg.eig(mtrx)

  sz = len(eigen_val)

  mtrx_out = np.zeros((sz,sz),dtype=complex)
  for row in range(len(eigen_val)):
    P = np.outer(eigen_vec[:,row], np.conj(eigen_vec[:,row]) )
    mtrx_out = np.add(mtrx_out,cmath.exp(fctr*eigen_val[row]) * P)

  if un_tets:
    non_Un = np.amax( \
      abs( \
        (mtrx_out.conj().T).dot(mtrx_out) - np.identity(sz) \
      ) \
    )
    if non_Un > 1.e-14:
      print("non_un = ", non_Un)
  else:
    non_Un = None

  return mtrx_out
# ======================================================================
