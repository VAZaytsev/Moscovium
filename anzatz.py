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
def XYZ_vec(tp, q, vec):
  if tp == 1:
    mask = [i ^ (1 << q) for i in range(len(vec))]

    return vec[mask]

  if tp == 2: 
    mask = [i ^ (1 << q) for i in range(len(vec))]
    sgn = np.array([ 1j if (i & (1<<q)) == 0 else 
                    -1j for i in range(len(vec))])

    return vec[mask] * sgn[mask]

  if tp == 3:
    return vec * np.array([ 1 if (i & (1<<q)) == 0 else 
                           -1 for i in range(len(vec))])

  exit("Wrong input for XYZ_vec! Abort")
# ===================================================================




# ===================================================================
def CX_vec(ctrl, targ, vec):
  mask = [i ^ (1 << targ) if (i & (1<<ctrl)) != 0 else 
          i for i in range(len(vec))]

  return vec[mask]
# ===================================================================




# ===================================================================
def oCX_vec(ctrl, targ, vec):
  mask = [i ^ (1 << targ) if (i & (1<<ctrl)) == 0 else 
          i for i in range(len(vec))]

  return vec[mask]
# ===================================================================




# ===================================================================
def Rq_vec(tp, theta, q, vec):
  if tp == 2: 
    mask = [i ^ (1 << q) for i in range(len(vec))]
    sgn = np.array([ 1 if (i & (1<<q)) == 0 else 
                    -1 for i in range(len(vec))])

    return math.cos(0.5*theta) * vec \
         + math.sin(0.5*theta) * vec[mask] * sgn[mask]

  if tp == 3:
    z1 = cmath.exp(-1j*0.5*theta)
    z2 = cmath.exp( 1j*0.5*theta)

    return vec * np.array([z1 if (i & (1<<q)) == 0 else 
                           z2 for i in range(len(vec))])
  
  if tp == 1: 
    mask = [i ^ (1 << q) for i in range(len(vec))]

    return math.cos(0.5*theta) * vec \
         - 1j * math.sin(0.5*theta) * vec[mask]

  if tp == 0:
    return vec

  exit("Wrong input for Rq_vec! Abort")
# ===================================================================




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

  exit("Wrong input for rot_mtrx_sub! Abort")
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

  if anzatz_tp == 2:
    return None

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
  
#  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
  if anzatz_tp == 2:
    psi_f = Rq_vec(2, ang[0], 0, psi_f)
    
    for ctrl in range(6):
      psi_f = CX_vec(ctrl, ctrl+1, psi_f)

    psi_f = Rq_vec(2, 0.5*ang[1], 1, psi_f)
    psi_f = CX_vec(0, 1, psi_f)
    psi_f = Rq_vec(2,-0.5*ang[1], 1, psi_f)
    psi_f = CX_vec(0, 1, psi_f)

    psi_f = CX_vec(2, 3, psi_f)

    psi_f = Rq_vec(2, 0.5*ang[2], 5, psi_f)
    psi_f = CX_vec(4, 5, psi_f)
    psi_f = Rq_vec(2,-0.5*ang[2], 5, psi_f)
    psi_f = CX_vec(4, 5, psi_f)

    psi_f = XYZ_vec(1, 3, psi_f)
    
    psi_f = oCX_vec(1, 2, psi_f)
    psi_f = oCX_vec(2, 3, psi_f)
    psi_f = oCX_vec(3, 4, psi_f)

    psi_f = Rq_vec(2, 0.5*ang[3], 5, psi_f)
    psi_f = CX_vec(4, 5, psi_f)
    psi_f = Rq_vec(2,-0.5*ang[3], 5, psi_f)
    psi_f = CX_vec(4, 5, psi_f)

    psi_f = oCX_vec(3, 4, psi_f)
    psi_f = oCX_vec(2, 3, psi_f)
    
    psi_f = oCX_vec(5, 6, psi_f)

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

#  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
  if anzatz_tp == 2:
    return None
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


#  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
  if anzatz_tp == 2:
    for i in range(2):
      if i == 1 and a == 0:
        break

      psi = Rq_vec(2, ang[0], 0, psi_i)
      if a == 0:
        psi = XYZ_vec(2, 0, psi)
      
      for ctrl in range(6):
        psi = CX_vec(ctrl, ctrl+1, psi)

      psi = Rq_vec(2, 0.5*ang[1], 1, psi)
      if a == 1 and i == 0:
        psi = XYZ_vec(2, 1, psi)
      psi = CX_vec(0, 1, psi)
      psi = Rq_vec(2,-0.5*ang[1], 1, psi)
      if a == 1 and i == 1:
        psi = XYZ_vec(2, 1, psi)
      psi = CX_vec(0, 1, psi)

      psi = CX_vec(2, 3, psi)

      psi = Rq_vec(2, 0.5*ang[2], 5, psi)
      if a == 2 and i == 0:
        psi = XYZ_vec(2, 5, psi)
      psi = CX_vec(4, 5, psi)
      psi = Rq_vec(2,-0.5*ang[2], 5, psi)
      if a == 2 and i == 1:
        psi = XYZ_vec(2, 5, psi)
      psi = CX_vec(4, 5, psi)

      psi = XYZ_vec(1, 3, psi)
      
      psi = oCX_vec(1, 2, psi)
      psi = oCX_vec(2, 3, psi)
      psi = oCX_vec(3, 4, psi)

      psi = Rq_vec(2, 0.5*ang[3], 5, psi)
      if a == 3 and i == 0:
        psi = XYZ_vec(2, 5, psi)
      psi = CX_vec(4, 5, psi)
      psi = Rq_vec(2,-0.5*ang[3], 5, psi)
      if a == 3 and i == 1:
        psi = XYZ_vec(2, 5, psi)
      psi = CX_vec(4, 5, psi)

      psi = oCX_vec(3, 4, psi)
      psi = oCX_vec(2, 3, psi)
      
      psi = oCX_vec(5, 6, psi)
      
      if i == 0:
        dpsi_f = psi
      else:
        dpsi_f -= psi
        dpsi_f *= 0.5

    return -0.5*1j*dpsi_f
# ===================================================================
