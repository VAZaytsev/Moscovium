import numpy as np
import math
import cmath


IXYZ = {
  0: np.matrix([[1, 0],[0, 1]]),
  1: np.matrix([[0, 1],[1, 0]]),
  2: np.matrix([[0,-1j],[1j, 0]],dtype=complex),
  3: np.matrix([[1, 0],[0,-1]]),
  }


# ===================================================================
def one_qubit_rot_mtrx(tp, theta):
  if tp == 0:
    return np.identity((2))

  if tp == 1:
    c = math.cos(0.5 * theta)
    s = math.sin(0.5 * theta)
    return np.array([[c,-1j*s],[-1j*s, c]])

  if tp == 2:
    c = math.cos(0.5 * theta)
    s = math.sin(0.5 * theta)
    return np.array([[c,-s],[s, c]])

  if tp == 3:
    return np.array([[cmath.exp(-1j*0.5*theta),0],[0, cmath.exp(1j*0.5*theta)]])

  print("Wrong input for rot_mtrx_sub! Abort")
  exit()
# ===================================================================




# ===================================================================
def two_qubit_rot_mtrx(theta):
  s = math.sin(theta)
  c = math.cos(theta)

  return np.array([[ 1,  0,  0,  0],
                   [ 0, -s,  c,  0],
                   [ 0,  c,  s,  0],
                   [ 0,  0,  0,  1]],dtype = float)
# ===================================================================




# ===================================================================
def deriv_mtrx(anzatz_tp, rot, nq):
  sz = 2**nq
  
  if anzatz_tp == 0:
    mtrx = np.zeros((nq, 4, sz, sz),dtype = complex)

    for r in range(4):
      for iq in range(nq):
        mtrx[iq,r,:,:] = -0.5 * 1j * np.kron( 
          np.kron( np.identity(2**(nq-1-iq) ),IXYZ[r] ),
          np.identity( 2**iq )
          ) # The order is important!!!

    return mtrx
  
  if anzatz_tp == 1:
    mtrx = np.zeros(( nq-1, sz, sz ),dtype=int)

    d_left = np.zeros((4,4), dtype=int)
    d_left[1,2] = -1
    d_left[2,1] =  1

    # Cycle over even numbers, the derivative goes to left
    for i in range(nq-1):
      mtrx[i,:,:] = np.kron(
        np.kron( np.identity( 2**(nq-2-i) ), d_left ), 
        np.identity( 2**i ) 
        ) # The order is important!!!

    return mtrx
# ===================================================================




# ===================================================================
def entangling_mtrx(nq):
  # ---.-------
  #    | 
  # ---x-.-----
  #      |
  # -----x-.---
  #        |
  # -------x---
  c0x1 = np.zeros((4,4))
  c0x1[0,0] = 1
  c0x1[2,2] = 1
  c0x1[1,3] = 1
  c0x1[3,1] = 1
  
  mtrx = np.identity(2**nq)
  for i in range(nq-1):
    tmp = np.kron( 
      np.kron( np.identity(2**(nq-2-i)), c0x1), 
      np.identity(2**i) 
      )

    mtrx = tmp.dot(mtrx) # The order is important!!!

  return mtrx
# ===================================================================




# ===================================================================
def apply_anzatz(anzatz_tp, rot, nlayers, ang, nq, psi_i, Uent=None):
  psi_f = psi_i.copy()


#  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
  if anzatz_tp == 0: 
    ang_per_layer = len(rot) * nq

    for ilayer in range(nlayers+1):
      for ir, r in enumerate(rot):

        mtrx = 1
        for iq in range(nq):
          i_ang = ilayer * ang_per_layer + ir*nq + iq

          mtrx = np.kron( 
            one_qubit_rot_mtrx(r, ang[i_ang]), mtrx 
            ) # The order is important!!!
        psi_f = mtrx.dot(psi_f)

      if ilayer != nlayers:
        psi_f = Uent.dot(psi_f)

    return psi_f


#  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
  if anzatz_tp == 1:
    ang_per_layer = nq - 1

    for ilayer in range(nlayers):

      mtrx = 1
      for iq in range(0,nq,2):
        i_ang = ilayer * ang_per_layer + int(iq/2)

        mtrx = np.kron( 
          two_qubit_rot_mtrx(ang[i_ang]), mtrx 
          ) # The order is important!!!
      if nq%2 == 1:
        mtrx = np.kron(np.identity(2), mtrx) # The order is important!!!

      psi_f = mtrx.dot(psi_f)


      mtrx = np.identity(2)
      for i in range(1,nq-1,2):
        i_ang = ilayer * ang_per_layer + math.floor(0.5*nq) + int((iq-1)/2)

        mtrx = np.kron( 
          two_qubit_rot_mtrx(ang[i_ang]), mtrx
          ) # The order is important!!!
      if nq%2 == 0:
        mtrx = np.kron(np.identity(2), mtrx) # The order is important!!!

      psi_f = mtrx.dot(psi_f)

    return psi_f
# ===================================================================




# ===================================================================
def anzatz_matrices(anzatz_tp, rot, nlayers, ang, nq):
  sz = 2**nq

#  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
  if anzatz_tp == 0: 
    ang_per_layer = len(rot) * nq

    mtrx_out = np.zeros((nlayers+1,len(rot),sz,sz),dtype=complex)

    for ilayer in range(nlayers+1):
      for ir, r in enumerate(rot):

        mtrx = 1
        for iq in range(nq):
          i_ang = ilayer * ang_per_layer + ir*nq + iq

          mtrx = np.kron( 
            one_qubit_rot_mtrx(r, ang[i_ang]), mtrx 
            ) # The order is important!!!
        mtrx_out[ilayer,ir,:,:] = mtrx.copy()

    return mtrx_out


#  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
  if anzatz_tp == 1:
    ang_per_layer = nq - 1

    mtrx_out = np.zeros((nlayers,2,sz,sz),dtype=complex)

    for ilayer in range(nlayers):

      mtrx = 1
      for iq in range(0,nq,2):
        i_ang = ilayer * ang_per_layer + int(iq/2)

        mtrx = np.kron( 
          two_qubit_rot_mtrx(ang[i_ang]), mtrx 
          ) # The order is important!!!
      if nq%2 == 1:
        mtrx = np.kron(np.identity(2), mtrx) # The order is important!!!

      mtrx_out[ilayer,0,:,:] = mtrx.copy()


      mtrx = np.identity(2)
      for i in range(1,nq-1,2):
        i_ang = ilayer * ang_per_layer + math.floor(0.5*nq) + int((iq-1)/2)

        mtrx = np.kron( 
          two_qubit_rot_mtrx(ang[i_ang]), mtrx
          ) # The order is important!!!
      if nq%2 == 0:
        mtrx = np.kron(np.identity(2), mtrx) # The order is important!!!

      mtrx_out[ilayer,1,:,:] = mtrx.copy()

    return mtrx_out
# ===================================================================




# ===================================================================
def apply_danzatz(anzatz_tp, rot, nlayers, ang, nq,
                  U, dU, a, psi_i, Uent=None):
  dpsi_f = psi_i.copy()

#  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
  if anzatz_tp == 0: 
    ang_per_layer = len(rot) * nq

    a_layer = math.floor(a / ang_per_layer)
    a_r = math.floor( (a - a_layer*ang_per_layer) / nq )
    a_q = a - a_layer*ang_per_layer - a_r*nq

    for ilayer in range(nlayers+1):
      for ir, r in enumerate(rot):
        dpsi_f = U[ilayer,ir,:,:].dot(dpsi_f)

        if a_layer == ilayer and a_r == ir:
          dpsi_f = dU[a_q,r,:,:].dot(dpsi_f)

      if ilayer != nlayers:
        dpsi_f = Uent.dot(dpsi_f)

    return dpsi_f


#  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
  if anzatz_tp == 1:
    ang_per_layer = nq - 1

    if nq%2 == 0:
      ang_per_sublayer = int(nq/2)
    else:
      ang_per_sublayer = int((nq-1)/2)

    a_layer = math.floor(a / ang_per_layer)
    a_r = math.floor( (a - a_layer*ang_per_layer) / ang_per_sublayer )
    a_q = ( a  - a_layer*ang_per_layer - a_r*ang_per_sublayer ) * 2 + a_r

    for ilayer in range(nlayers):
      for ir in range(2):
        dpsi_f = U[ilayer,ir,:,:].dot(dpsi_f)

        if a_layer == ilayer and a_r == ir:
          dpsi_f = dU[a_q,:,:].dot(dpsi_f)
      
  return dpsi_f
# ===================================================================
