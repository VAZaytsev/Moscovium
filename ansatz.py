import numpy as np
import math
import cmath

import inspect

IXYZ = {
  0: np.matrix([[1, 0],[0, 1]]),
  1: np.matrix([[0, 1],[1, 0]]),
  2: np.matrix([[0,-1j],[1j, 0]],dtype=complex),
  3: np.matrix([[1, 0],[0,-1]]),
  }


# matrices of different gates =======================================
def Rx_mtrx(theta):
  c = math.cos(0.5 * theta)
  s = math.sin(0.5 * theta)
  return np.array([[c,-1j*s],[-1j*s, c]])
def dRx_mtrx(theta):
  c = math.cos(0.5 * theta)
  s = math.sin(0.5 * theta)
  return -0.5*np.array([[s,1j*c],[1j*c, s]])


def Ry_mtrx(theta):
  c = math.cos(0.5 * theta)
  s = math.sin(0.5 * theta)
  return np.array([[c,-s],[s, c]])
def dRy_mtrx(theta):
  c = math.cos(0.5 * theta)
  s = math.sin(0.5 * theta)
  return -0.5*np.array([[s,c],[-c,s]])


def Rz_mtrx(theta):
  return np.array([[cmath.exp(-1j*0.5*theta),0],
                   [0, cmath.exp(1j*0.5*theta)]])
def dRz_mtrx(theta):
  return 0.5*1j*np.array([[cmath.exp(-1j*0.5*theta),0],
                          [0, cmath.exp(1j*0.5*theta)]])


C0X1_mtrx = np.matrix([[1,0,0,0],
                       [0,0,0,1],
                       [0,0,1,0],
                       [0,1,0,0]])


def fSim_3cx_mtrx(theta):
  s = math.sin(theta)
  c = math.cos(theta)

  return np.array([[ 1,  0,  0,  0],
                   [ 0, -s,  c,  0],
                   [ 0,  c,  s,  0],
                   [ 0,  0,  0,  1]],dtype = float)
def dfSim_3cx_mtrx(theta):
  s = math.sin(theta)
  c = math.cos(theta)

  return np.array([[ 0,  0,  0,  0],
                   [ 0, -c, -s,  0],
                   [ 0, -s,  c,  0],
                   [ 0,  0,  0,  0]],dtype = float)
# ===================================================================



# ===================================================================
class ansatz_cls:
  def __init__(self):
    self.layers = []

  def clc_nparams(self):
    self.nparams = sum([layer.get_num_params() 
                        for layer in self.layers])

  def act_on_vctr(self, args, vctr_in):
    vctr_out = vctr_in.copy()

    indx = 0 
    for layer in self.layers:
      nargs = layer.get_num_params()
      vctr_out = layer.act_on_vctr(args[indx:indx+nargs], vctr_out)
      indx += nargs
    return vctr_out


  def dact_on_vctr(self, args, vctr_in):
    res = [vctr_in.copy()]*len(args)

    indx = 0 
    for layer in self.layers:
      nargs = layer.get_num_params()

      for i in range(indx):
        res[i] = layer.act_on_vctr(args[indx:indx+nargs], res[i])

      if nargs != 0:
        res[indx:indx+nargs] = layer.dact_on_vctr(args[indx:indx+nargs], res[indx])

      tmp = layer.act_on_vctr(args[indx:indx+nargs], res[-1])
      for i in range(indx+nargs,len(args)):
        res[i] = tmp.copy()

      indx += nargs
    return res
# ===================================================================


# ===================================================================
def hardware_efficient_ansatz(Nlayers, Nq, rot):
  # in rot can be any number of rotations
  # values can be only 1,2, or 3 related to Rx, Ry, and Rz, respectively
  #
  # --- rot[1](ang_0) --- rot[2](ang_4) ---.-----
  #                                        |
  # --- rot[1](ang_1) --- rot[2](ang_5) ---x-.---
  #                                          |
  # --- rot[1](ang_2) --- rot[2](ang_6) ---.-x---
  #                                        |
  # --- rot[1](ang_3) --- rot[2](ang_7) ---x-----
  #
  # + final layer of rotations
  ansatz = ansatz_cls()
  layer_rot = []
  for r in rot:
    layer_rot.append( layer_cls() )
    for iq in range(Nq):
      if int(r) == 1:
        layer_rot[-1].add_gate(iq, iq, Rx_mtrx, dRx_mtrx)
      if int(r) == 2:
        layer_rot[-1].add_gate(iq, iq, Ry_mtrx, dRy_mtrx)
      if int(r) == 3:
        layer_rot[-1].add_gate(iq, iq, Rz_mtrx, dRz_mtrx)
    layer_rot[-1].sort()


  layer_ent1 = layer_cls()
  layer_ent2 = layer_cls()
  for iq in range(0,Nq,2):
    if iq+1 < Nq:
      layer_ent1.add_gate(iq+1,iq,C0X1_mtrx)
    if iq+1 < Nq and iq+2 < Nq:
      layer_ent2.add_gate(iq+2,iq+1,C0X1_mtrx)

  for ilayer in range(Nlayers):
    for layer in layer_rot:
      ansatz.layers.append(layer)

    ansatz.layers.append(layer_ent1)
    ansatz.layers.append(layer_ent2)

  # last layer of rotations
  for layer in layer_rot:
    ansatz.layers.append(layer)

  ansatz.clc_nparams()
  return ansatz
# ===================================================================


# ===================================================================
def two_qubit_rot_ansatz(Nlayers, Nq):
  #      _______
  # --- |       | -------------
  #     | ang_0 |      _______
  # --- |_______| --- |       |
  #      _______      | ang_2 |
  # --- |       | --- |_______|
  #     | ang_1 |
  # --- |_______| -------------
  #
  ansatz = ansatz_cls()

  layer_ent1 = layer_cls()
  layer_ent2 = layer_cls()
  for iq in range(0,Nq,2):
    if iq+1 < Nq:
      layer_ent1.add_gate(iq+1, iq, fSim_3cx_mtrx, dfSim_3cx_mtrx)

    if iq+1 < Nq and iq+2 < Nq:
      layer_ent2.add_gate(iq+2, iq+1, fSim_3cx_mtrx, dfSim_3cx_mtrx)
      
  for ilayer in range(Nlayers):
    ansatz.layers.append(layer_ent1)
    ansatz.layers.append(layer_ent2)
  
  ansatz.clc_nparams()
  return ansatz
# ===================================================================


# ===================================================================
def UCC_ansatz(ClusterOps, exc):
  ansatz = ansatz_cls()

  def mtrx_fun(ClOp):
    def f(x):
      res = np.zeros( ClOp.P[0].shape, dtype=complex )
      for i, e in enumerate(ClOp.e_val):
        res += cmath.exp(e * x) * ClOp.P[i]
      return res
    return f

  def dmtrx_fun(ClOp):
    def f(x):
      res = np.zeros( ClOp.P[0].shape, dtype=complex )
      for i, e in enumerate(ClOp.e_val):
        res += e * cmath.exp(e * x) * ClOp.P[i]
      return res
    return f

  nq = ClusterOps[0].Tq.num_qubits

  for ClusterOp in ClusterOps:
    if ClusterOp.nex in exc:
      layer = layer_cls()
      layer.add_gate(nq-1, 0, 
                     mtrx_fun(ClusterOp), dmtrx_fun(ClusterOp))
      ansatz.layers.append(layer)

  ansatz.clc_nparams()
  return ansatz
# ===================================================================


# ===================================================================
class layer_cls:
  def __init__(self):
    self.rules = []
    self.ql = []
    self.qr = []
    self.drules = []

  
  def add_gate(self, ql, qr, mtrx_fun, dmtrx_fun=None):
    if ql < qr:
      exit("Wrong values for ql and qr")
      
    if any([qr in range(self.qr[i],self.ql[i]+1) 
            for i in range(len(self.qr))] ):
      exit("ql and qr overlap with previous matrices")

    self.rules.append( mtrx_fun )
    self.ql.append( ql )
    self.qr.append( qr )
    self.drules.append( dmtrx_fun )


  def get_num_params(self):
    nparam = 0
    for rule in self.rules:
      if inspect.isfunction(rule):
        nparam += len(inspect.signature(rule).parameters)

    return nparam

  def sort(self):
    indxs = sorted(range(len(self.qr)), key=lambda i: self.qr[i])
    self.rules = [self.rules[i] for i in indxs]
    self.ql = [self.ql[i] for i in indxs]
    self.qr = [self.qr[i] for i in indxs]
    self.drules = [self.drules[i] for i in indxs]


  def act_on_vctr(self, args, vctr):
    nq = (len(vctr)-1).bit_length()

    tmp = 1
    indx = 0
    l = -1
    for i,rule in enumerate(self.rules):
      if inspect.isfunction(rule):
        nargs = len(inspect.signature(rule).parameters)
        mtrx = rule(*tuple(args[indx:indx+nargs]))
        indx += nargs
      else:
        mtrx = rule

      tmp = np.kron( np.kron( mtrx, np.identity( 2**(self.qr[i]-1-l) ) ), tmp)
      l = self.ql[i]

    tmp = np.kron( np.identity( 2**(nq-1-l ) ), tmp )

    return tmp.dot(vctr).reshape(2**nq,1)


  def dact_on_vctr(self, args, vctr):
    nq = (len(vctr)-1).bit_length()
    
    res = [vctr.copy()] * len(args)

    indx = 0
    l = -1
    tmp_arr = [1] * len(args)
    for i,rule in enumerate(self.rules):
      nargs = 0
      if inspect.isfunction(rule):
        nargs = len(inspect.signature(rule).parameters)
        mtrx = rule(*tuple(args[indx:indx+nargs]))

        if self.drules[i] != None:
          dmtrx = self.drules[i](*tuple(args[indx:indx+nargs]))
          
          if nargs == 1:
            tmp_arr[indx] = np.kron( np.kron( dmtrx, np.identity( 2**(self.qr[i]-1-l) ) ), tmp_arr[indx])
          else:
            for ii,d in enumerate(dmtrx):
              j = ii + indx
              tmp_arr[j] = np.kron( np.kron( d, np.identity( 2**(self.qr[i]-1-l) ) ), tmp_arr[j])
        else:
          for ii in range(nargs):
            j = ii + indx
            tmp_arr[j] = np.zeros(( 2**(self.ql[i]+1), 2**(self.ql[i]+1) ))
      else:
        mtrx = rule


      for j,tmp in enumerate(tmp_arr):
        if j in range(indx,indx+nargs):
          continue
        tmp_arr[j] = np.kron( np.kron( mtrx, np.identity( 2**(self.qr[i]-1-l) ) ), tmp)

      indx += nargs
      l = self.ql[i]

    for i,tmp in enumerate(tmp_arr):
      tmp_arr[i] = np.kron( np.identity( 2**(nq-1-l ) ), tmp )
      res[i] = tmp_arr[i].dot( res[i] ).reshape(2**nq,1)
    
    return res
# ===================================================================
