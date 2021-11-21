import os
os.environ["OMP_NUM_THREADS"] = "8" # export OMP_NUM_THREADS=4


 
import numpy as np
from numpy import pi

import math
import cmath
import scipy

import sympy

import time

import sys

from qiskit_nature.problems.second_quantization.electronic.builders.fermionic_op_builder import build_ferm_op_from_ints
from qiskit_nature.converters.second_quantization import QubitConverter
from qiskit_nature.mappers.second_quantization import JordanWignerMapper, ParityMapper

from qiskit_nature.operators.second_quantization import FermionicOp

from find_a import *
from second_q_op_vaz import *
from ps_mod import *
from ansatz import *
from routines import *

import tappering
import optimizers

from quantum_mod import *

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# a- = q+ = 0.5(X + iY)
# a+ = q- = 0.5(X - iY)
#
# N = a+a- = q-q+ = 0.5*(I - Z)
#
# NII = mtrx(I) (x) mtrx(I) (x) mtrx(N)
# I+ = mtrx(+) (x) mtrx(Z)

# q3 q2 q1 q0 - numeration of qubits
# o3 o2 o1 o0 - numeration of orbitals
# f0 f1 f2 f3 - numeration of fermionic operators
#
# a+_1 = I+II = mtrx(I) (x) mtrx(I) (x) mtrx(a+) (x) mtrx(Z)

name2l = {
  "s": 0,
  "p": 1,
  "d": 2
  }

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def jl2k(j,l):
  return (-1)**int(l+j+0.5) * int(j + 0.5)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

class rad_orb:
  def __init__(self,n,l,j):
    self.n = n
    self.j = j
    self.l = l

    self.k = (-1)**int(l+j+0.5) * int(j + 0.5)

    self.i = -1
  def __eq__(self, other):
    return self.k == other.k and self.n == other.n
  
  def __str__(self):
    return str(self.n) + " " + str(self.k)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

class orb_cls:
  def __init__(self,n,k,m):
    self.n = n
    self.k = k
    self.m = m

    self.i = -1

    self.j = abs(k) - 0.5
    self.l = int( abs(k + 0.5) - 0.5 )

  def __eq__(self, other):
    return self.k == other.k and self.n == other.n and self.m == other.m

  def __str__(self):
    return str(self.n) + " " + str(self.k) + " " + str(self.m)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def find_orb(arr, indx, m):
  for i,x in enumerate(arr):
    if indx == x.i and m == x.m:
      return i
  return -1  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def yellow_highlighting(word):
  return '\033[30;103m' + word + '\033[0m'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Read files --------------------------------------------------------
#fl_1b_int = open(sys.argv[1],"r")
#fl_2b_int = open(sys.argv[2],"r")
#norb_max = int(sys.argv[3])
fl_inp = open(sys.argv[1],"r")

for ln in fl_inp.readlines():
  if ln == "\n":
    continue

  key, val = ln.strip().replace(" ","").split("=")
  if key == "OBI":
    fl_1b_int = open(val,"r")

  if key == "TBI":
    fl_2b_int = open(val,"r")

  if key == "Nelec":
    Nelec = int(val)

  if key == "norb":
    norb_max = int(val)

  if key == "ansatz_tp":
    ansatz_tp = int(val)

  if key == "Nlayers":
    Nlayers = int(val)

  if key == "rot":
    rot = [int(x) for x in val.replace("[","").replace("]","").split(",")]

  if key == "Q_on":
    Q_on = (val == "True")

  if key == "noise_off":
    noise_off = (val == "True")

  if key == "Nrep":
    Nrep = int(val)

  if key == "Parity":
    Parity_tot = int(val)

#exit()


# reading one-electron orbitals - - - - - - - - - - - - - - - - - - -
read_states = False
rad_orb_arr = []
for ln in fl_1b_int.readlines():
  if ln == "\n":
    continue

  arr = ln.split()

  if arr[0] == "Ni":
    read_states = True
    continue
  
  if arr[0] == "n1":
    break

  if read_states:
    sz = len(arr[1])
    n = int( arr[1][:sz-1] )
    l = name2l[ arr[1][sz-1] ]
    j = 0.5*float( arr[3].split("/")[0] )
    indx = int(arr[0])
    
    rad_orb_arr.append( rad_orb(n,l,j) )
    rad_orb_arr[-1].i = indx
# End reading one-electron orbitals - - - - - - - - - - - - - - - - -


# creating spin-orbitals - - - - - - - - - - - - - - - - - - - - - - 
orb_arr = []

# There can be not enough space for all orbitals
# some orbitals are deleted
norb = 0
delete_from = -1
for i,x in enumerate(rad_orb_arr):
  if norb + 2*x.j + 1 > norb_max:
    delete_from = i
    break
  norb += int(2*x.j + 1)


sz = len(rad_orb_arr) - delete_from
for i in range(sz):
  if delete_from == -1:
    continue
  rad_orb_arr.pop(delete_from)
  
# sort orbitals with respect to parity
rad_orb_arr = sorted( rad_orb_arr, key = lambda x: (-1)**x.l, reverse=False )


for x in rad_orb_arr:
  for m in range( 1, int(2*x.j+1)+1, 2 ):
    for s in range(-1,2,2):
      orb_arr.append( orb_cls(x.n, x.k, 0.5*m*s) )
      orb_arr[-1].i = x.i


print("\nFollowing orbitals are used")
iq_even_first = -1
print("Odd")
for i,x in enumerate(orb_arr):
  if x.l%2 == 0 and iq_even_first == -1:
    iq_even_first = i
    print("Even")
    
  print("q"+str(i), x, x.l)
print()

one_b_int = np.zeros([norb,norb])
two_b_int = np.zeros([norb,norb,norb,norb])

#exit()
fl_1b_int.seek(0)
# end creating spin-orbitals - - - - - - - - - - - - - - - - - - - - 


# Read one-body integrals - - - - - - - - - - - - - - - - - - - - - -
read_int = False
for ln in fl_1b_int.readlines():
  if ln == "\n":
    continue

  arr = ln.split()
  if arr[0] == "n1":
    read_int = True
    continue

  if read_int:
    i = int(arr[0])
    j = int(arr[1])
    h_ij = float( arr[2] )
    
    for iorb,x in enumerate(orb_arr):
      for jorb,y in enumerate(orb_arr):
        if x.m == y.m and i == x.i and j == y.i:
          one_b_int[iorb,jorb] = h_ij
          one_b_int[jorb,iorb] = h_ij
# End reading one-body integrals - - - - - - - - - - - - - - - - - - 


# Read two-body integrals - - - - - - - - - - - - - - - - - - - - - -
read_int = False
for ln in fl_2b_int.readlines():
  if ln == "\n":
    continue

  arr = ln.split()

  if arr[0] == "n1":
    read_int = True
    continue

  if read_int:
    i = np.zeros((4),dtype=int)
    for ii in range(4):
      i[ii] = int( arr[ii] )

    m = np.zeros((4),dtype=float)
    for ii in range(4):
      m[ii] = 0.5*float( arr[4+ii] )
      
    u_coul_dir = float( arr[8] )
    u_coul_ex = float( arr[9] )
    u_br_dir = float( arr[10] )
    u_br_ex = float( arr[11] )

    p = find_orb(orb_arr, i[0], m[0])
    if p == -1:
      continue
    q = find_orb(orb_arr, i[1], m[1])
    if q == -1:
      continue
    r = find_orb(orb_arr, i[2], m[2])
    if r == -1:
      continue    
    s = find_orb(orb_arr, i[3], m[3])
    if s == -1:
      continue

    # direct
    two_b_int[p,q,r,s] = u_coul_dir + u_br_dir
    two_b_int[q,p,s,r] = u_coul_dir + u_br_dir

    two_b_int[s,r,q,p] = u_coul_dir + u_br_dir
    two_b_int[r,s,p,q] = u_coul_dir + u_br_dir

    # exchange
    two_b_int[p,q,s,r] = u_coul_ex + u_br_ex
    two_b_int[q,p,r,s] = u_coul_ex + u_br_ex

    two_b_int[r,s,q,p] = u_coul_ex + u_br_ex
    two_b_int[s,r,p,q] = u_coul_ex + u_br_ex

# ===================================================================


# Number of particles -----------------------------------------------
N_part_op = build_ferm_op_from_ints(one_body_integrals=np.identity(norb))
# -------------------------------------------------------------------


# Jz ----------------------------------------------------------------
tmp = np.identity(norb)
for i,orb in enumerate(orb_arr):
  tmp[i,i] = orb.m

Jz_op = build_ferm_op_from_ints(one_body_integrals=tmp)
# -------------------------------------------------------------------


# Parity ------------------------------------------------------------
tmp = np.identity(norb)
for i,orb in enumerate(orb_arr):
  tmp[i,i] = orb.l

Parity_op = build_ferm_op_from_ints(one_body_integrals=tmp)
# -------------------------------------------------------------------


# Hamiltonian -------------------------------------------------------
#qubit_converter = QubitConverter(mapper=JordanWignerMapper())
qubit_converter = QubitConverter(mapper=ParityMapper())

H_op = 0*FermionicOp("I"*norb).reduce()


# This term will add some energy to the states with wrong particle number
#if anzatz_tp == 0:
  #N_part_op = build_ferm_op_from_ints(one_body_integrals=np.identity(norb))
  #H_op += 10*(N_part_op - Nelec*FermionicOp("I"*norb))**2


#Add one-body integrals
for p in range(norb):
  for q in range(norb):
    if one_b_int[p,q] == 0:
      continue
    
    coef, lbl = second_q_op_pq(p, q, norb)
    H_op += coef * FermionicOp(lbl) * one_b_int[p,q]


#Add two-body integrals
for p in range(norb):
  for q in range(norb):
    if p == q:
      continue
    for r in range(norb):
      for s in range(norb):
        if r == s:
          continue

        if two_b_int[p,q,r,s] != 0.0:
          coef, lbl = second_q_op_pqrs(p, q, s, r, norb)

          H_op += 0.5 * coef * FermionicOp(lbl) * two_b_int[p,q,r,s]

H_op = H_op.reduce()
#print("\n",H_op)
#exit()


# Rewrite Hamiltonian as a sum of Pauli strings ---------------------
H_q = qubit_converter.convert(H_op)
#for h in H_q.primitive.to_list():
  #print(h)
Parity_q = qubit_converter.convert(Parity_op)
N_part_q = qubit_converter.convert(N_part_op)
Jz_q = qubit_converter.convert(Jz_op)



Nq = H_q.num_qubits
print( "Nq before tappering = ", Nq, flush=True )



# Tappering qubits --------------------------------------------------
Tappering_on = True
if Tappering_on:
  # tapper the qubit which is responsible for the total number of 
  # particles

  H_q = tappering.tapper(H_q, Nq-1, Nelec%2)
  N_part_q = tappering.tapper(N_part_q,Nq-1,Nelec%2)
  Jz_q = tappering.tapper(Jz_q,Nq-1,Nelec%2)
  Parity_q = tappering.tapper(Parity_q, Nq-1, Nelec%2)


  # tapper the qubit which is responsible for the total number of 
  # particles
  H_q = tappering.tapper(H_q, iq_even_first-1, Parity_tot)
  N_part_q = tappering.tapper(N_part_q,iq_even_first-1,Parity_tot)
  Jz_q = tappering.tapper(Jz_q,iq_even_first-1,Parity_tot)
  Parity_q = tappering.tapper(Parity_q, iq_even_first-1, Parity_tot)



Nq = H_q.num_qubits
print( "Nq after tappering = ", Nq, flush=True )

#exit()
# -------------------------------------------------------------------



# Matrices ----------------------------------------------------------
H_mtrx = H_q.to_matrix().real
N_part_mtrx = N_part_q.to_matrix().real
Jz_mtrx = Jz_q.to_matrix().real

Parity_mtrx = Parity_q.to_matrix().real
for i in range(len(Parity_mtrx[:,0])):
  Parity_mtrx[i,i] %= 2


# Find eigenvalues and eigenvectors of the Hamiltonian
energy, wf = np.linalg.eigh(H_mtrx)
#for e in energy:
  #print(e)
#exit()


# Diagonalize Jz matrix ---------------------------------------------
indx_arr = [0]
for i in range(1,2**Nq):
  diff = energy[i]-energy[i-1]

  if diff < 5.e-12:
    indx_arr.append(i)
    continue

  # Cluster of the states with the same energies, 
  # number of particles, and parities is created
  # Calculate Jz matrix for this cluster
  sz = len(indx_arr)
  Jz_mtrx_new = np.zeros((sz,sz))
  for r in range(sz):
    for c in range(r, sz):
      jz = np.conj(wf[:,indx_arr[r]].T).dot( 
        Jz_mtrx.dot( wf[:,indx_arr[c]]) 
        ).item(0)

      Jz_mtrx_new[r,c] = jz
      Jz_mtrx_new[c,r] = jz


  val, vec = np.linalg.eigh(Jz_mtrx_new)

  wf[:,indx_arr] = np.matmul(wf[:,indx_arr],vec)
  
  # Preparation for a new cluster
  indx_arr = [i]
#exit()
# -------------------------------------------------------------------



# -------------------------------------------------------------------
class sl_dets_with_symmetries:
  def __init__(self, bn, n, p, j):
    self.bn_arr = [bn]
    self.n = n
    self.p = p
    self.j = j


SD_sym = []
for i in range(2**Nq):
  bn = bin(i)[2:].zfill(Nq)

  wf_bn = bin_to_vec(bn)

  jz = np.conj( wf_bn.T ).dot( Jz_mtrx.dot( wf_bn )).item(0)
  parity = np.conj( wf_bn.T ).dot( Parity_mtrx.dot( wf_bn )).item(0)
  npart = np.conj(wf_bn.T).dot( N_part_mtrx.dot( wf_bn) ).item(0)

  found = False
  for el in SD_sym:
    if el.n == npart and el.p == parity and el.j == jz:
      el.bn_arr.append(bn)
      found = True
      break

  if not found:
    SD_sym.append( sl_dets_with_symmetries(bn, npart, parity, jz) )



counter = 0
for el in SD_sym:
  if int(el.n) != Nelec:
    continue

  if el.p != Parity_tot:
    continue

  print("J_z = ", el.j, "P = ", int(el.p), "N_SlDet = ", len(el.bn_arr))
  for b in el.bn_arr:
    print(b)
  print()
#exit()
# -------------------------------------------------------------------




# Find the energy of the ground and first excited states
# for each possible number of electrons
ne_arr = range(norb+1)
indx_g = [-1]*(norb+1)
indx_e = [-1]*(norb+1)

for i in range(2**Nq):
  npart = np.conj(wf[:,i].T).dot( N_part_mtrx.dot( wf[:,i]) ).item(0)  
  jz = np.conj(wf[:,i].T).dot( Jz_mtrx.dot( wf[:,i]) ).item(0)

  #if int( round(npart)) == Nelec:
    #print( energy[i], jz )

  for ii, x in enumerate(ne_arr):
    if int(round(npart)) == x:
      if indx_g[ii] == -1:
        indx_g[ii] = i
        break
      # condition abs(...) > 1.e-13 is needed to get rid of degenerate states
      if indx_e[ii] == -1 and abs(energy[i] - energy[indx_g[ii]]) > 1.e-13:
        indx_e[ii] = i
        #break


    
print("\n","Ne", 
      f'{"Eg": >8}', 
      f'{"Ee": >12}', 
      f'{"Jz": >4}')
for i, ne in enumerate(ne_arr):
  #if ne != Nelec:
    #continue

  if indx_g[i] == -1:
    continue

  #indx = np.argsort( abs(wf[:,indx_g[i]]), axis=0 )

  jz_g = np.conj( wf[:,indx_g[i]].T ).dot( 
    Jz_mtrx.dot( wf[:,indx_g[i]] )).item(0)
  
  parity_g = int(np.conj( wf[:,indx_g[i]].T ).dot( 
    Parity_mtrx.dot( wf[:,indx_g[i]] )).item(0))


  if indx_e[i] != -1:
    print(f'{ne: >3}', 
          f'{energy[indx_g[i]]: 10.5f}', 
          f'{energy[indx_e[i]]: 10.5f}', 
          f'{jz_g: 4.1f}')
  else:
    print(f'{ne: >3}', 
          f'{energy[indx_g[i]]: 10.5f}', 
          " "*10,
          f'{jz_g: 4.1f}')

#exit()
# -------------------------------------------------------------------




# Imaginary time propagation ========================================
#Q_on = True
#noise_off = True
create_noise_model(Nq, Nrep)


# Not all Pauli matrices have to be measured - - - - - - - - - - - - 
if Q_on:
  ps_meas, indx_ps2meas = extract_ps_for_measurement(H_q)
  print("Hamiltonian consists of",len(indx_ps2meas),"PS", flush=True )
  print("It is sufficient to measure", len(ps_meas), "of them", flush=True)
#exit()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


dtau = 1.e-2
# In principle the value of dtau and the difference in energy between 
# ground and excited states defines the number of necessary steps
# for propagation in imaginary time


# Construct exp(-tau H) matrix ======================================
exp_tau_H = np.zeros((2**Nq,2**Nq))
for i in range(2**Nq):
  P = np.outer( wf[:,i], np.conj(wf[:,i]) )
  exp_tau_H += math.exp(-dtau*energy[i]) * P
# ===================================================================



# Preparation -------------------------------------------------------
ansatz = ansatz_cls(ansatz_tp, Nlayers, Nq, rot)

Nparam = ansatz.nparams 
print("Total number of parameters in anzatz = ", Nparam, flush=True)

rot_angs = math.pi*(2*np.random.rand(Nparam)-1)
#rot_angs = np.zeros(Nparam)

# Exact ground-state wave function for given number of electrons


psi_exact = wf[:,indx_g[Nelec]].copy()
indx = np.argsort( abs(psi_exact), axis=0 )

# Dominant contribution comes from the Slater determinant
# with the binary representation 
psi_i_bn = bin( 
  int( np.argsort( abs(psi_exact), axis=0 )[-1] ) 
  )[2:].zfill(Nq)
print("\n","Initial guess = ", psi_i_bn, flush=True)


jz_g = average_value_for_bn_mtrx(psi_i_bn, Jz_mtrx)
print("For this state: jz = ", jz_g,
      "Energy = ", average_value_for_bn_mtrx(psi_i_bn, H_mtrx))



print("\n","Exact wave function")
print(f'{"weight": >18}', " "*2, f'{"contribution to E": >17}')
for i in range(2**Nq):
  ii = indx[2**Nq-1-i]
  bn_ii = bin(ii)[2:].zfill(Nq)

  npart = average_value_for_bn_mtrx(bn_ii, N_part_mtrx)
  if npart != Nelec:
    continue

  jz = average_value_for_bn_mtrx(bn_ii, Jz_mtrx)
  if jz != jz_g:
    continue

  parity = average_value_for_bn_mtrx(bn_ii, Parity_mtrx)
  if parity != Parity_tot:
    continue

  print( f'{i: >2}', bn_ii, 
        f'{psi_exact[ ii ]**2: .6f}',
        average_value_for_bn_mtrx(bn_ii, H_mtrx)*psi_exact[ ii ]**2)


     
psi_in = bin_to_vec(psi_i_bn)


def calculate_E(angs, psi_in):
    psi_out = ansatz.act_on_vctr(angs, psi_in)

    E = (psi_out.conj().T).dot( H_mtrx.dot(psi_out) ).item().real
    return E


adam = optimizers.Adam_cls(Nparam, eta=0.05)
NG = optimizers.NatGrad_cls(Nparam, eta=0.05)
ITE = optimizers.ITE_cls(Nparam, eta=0.05)

optims = [adam, NG, ITE]
#optims = [adam, NG]
#optims = [adam]

angs = [rot_angs.copy()] * len(optims)
E_old = [None] * len(optims)

for i,opt in enumerate(optims):
  opt.f = calculate_E(angs[i], psi_in)


Niters_max = 2000
eps = 1.e-8

txt = ""
for opt in optims:
  txt += "   " + opt.name + "   "
print("\n"," "*5, txt)  

iters = 1
while any([not x.converged for x in optims]) and \
  all([x.t <= Niters_max for x in optims]):

  txt_out = ""

  for i,opt in enumerate(optims):
    if opt.converged:
      txt_out += " "*10
      continue
    
    angs[i] = optimizers.update_angles(ansatz, angs[i], 
                                       opt, psi_in, H_mtrx)
    E_old[i] = opt.f
    opt.f = calculate_E(angs[i], psi_in)

    if abs(opt.f - E_old[i]) < eps:
      opt.converged = True

    txt_out += f' {opt.f: .6f}'

  if iters%50 == 0:
    print(f'{iters:>5}', txt_out)
    
  iters += 1
      
      

psi_out = []
txt_E = f'{energy[indx_g[Nelec]]: .5f}'
txt_name = "  Exact   "
for i,opt in enumerate(optims):
  psi_out.append( ansatz.act_on_vctr(angs[i], psi_in) )

  txt_name += "   " + opt.name + "   "
  txt_E += f' {opt.f: .6f}'

print("\n",txt_name, "\n", txt_E, "\n")



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#for i in range(2**Nq):
  #if any([abs(psi[i])**2 > 1.e-5 for psi in psi_out]):
    #bn = bin(i)[2:].zfill(Nq)

    #wf_bn = bin_to_vec(bn)

    #jz = np.conj( wf_bn.T ).dot( Jz_mtrx.dot( wf_bn )).item(0)
    #npart = np.conj(wf_bn.T).dot( N_part_mtrx.dot( wf_bn) ).item(0)

    #txt_out = str(bn)
    #if abs(psi_exact[i])**2 > 1.e-5:
      #txt_out = yellow_highlighting(txt_out)

    #txt_out += f'{abs(psi_exact[i])**2: .5f}'

    #for psi in psi_out:
      #txt_out += f'{abs(psi[i].item())**2: .5f}'
      
    #print( txt_out, jz, npart )
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


## Newton-CG - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#obj = scipy.optimize.minimize(fun, rot_ang_c, 
                              ##method='Nelder-Mead',
                              ##method='COBYLA',
                              #method='Newton-CG', jac=dfun, 
                              ##method='CG', jac=dfun,
                              ##method='BFGS', jac=dfun,
                              #tol=1.e-12,
                              #options={"maxiter": 10000, "disp": False })

#rot_ang_bst = obj.x
#psi_f_vec = apply_anzatz(anzatz_tp, rot, Nlayers, rot_ang_bst, Nq,
                         #psi_i_vec, Uent)
#E_Newton = (psi_f_vec.conj().T).dot( H_mtrx.dot(psi_f_vec) ).item().real
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -





## VQE ===============================================================
#VQE_on = False
#if VQE_on:
  #E_q = measure_Ham(anzatz_tp, psi_i_circ, Uent_circ,
                    #Nlayers, rot, rot_ang_bst, 
                    #ps_meas, indx_ps2meas, H_q, 
                    #noise_off=noise_off)
  #print("E_VQE = ", E_q)
## ===================================================================


## Test new imaginary time propagation ===============================
#psi = apply_anzatz(anzatz_tp, rot, Nlayers, rot_ang_c, Nq, psi_i_vec, Uent)
#Nstps, converged, H_av = get_Nsteps_for_exp_tauH(psi, H_mtrx, exp_tau_H)
#if converged:
  #print("Converged in ", Nstps, "got <H> = ", H_av)
#else:
  #print("NOT converged in ", Nstps, "got <H> = ", H_av)

#exit()
## ===================================================================


## Propagation in accordance to McArdle ==============================
#psi_tauH = apply_anzatz(anzatz_tp, rot, Nlayers, rot_ang_c, Nq,
                        #psi_i_vec, Uent)

#U = anzatz_matrices(anzatz_tp, rot, Nlayers, rot_ang_c, Nq)
#for it in range(10000):
  #A_mtrx_c = np.zeros((Nparam,Nparam),dtype=float)
  #C_vec_c = np.zeros((Nparam),dtype=float)

  #if Q_on:
    #A_mtrx_q = np.zeros((Nparam,Nparam),dtype=float)
    #C_vec_q = np.zeros((Nparam),dtype=float)

  #psi_c = apply_anzatz(anzatz_tp, rot, Nlayers, rot_ang_c, Nq,
                       #psi_i_vec, Uent)

  ## apply exact exp(-tau H) - - - - - - - - - - - - - - - - - - - - -
  #psi_tauH = exp_tau_H.dot(psi_c)
  #psi_tauH /= np.linalg.norm(psi_tauH)


## Calculate C -------------------------------------------------------
  #for a in range(Nparam):
    #dpsi_c_a = apply_danzatz(anzatz_tp, rot, Nlayers, rot_ang_c, Nq, 
                               #U, dU, a, psi_i_vec, Uent)

    #C_vec_c[a] = -2 * (dpsi_c_a.conj().T).dot( H_mtrx.dot(psi_c) ).item(0).real


    #if Q_on:
      #C_vec_q[a] = -clc_c_vec(anzatz_tp, psi_i_circ, Uent_circ, 
                             #Nlayers, rot, rot_ang_q, a,
                             #ps_meas, indx_ps2meas, H_q,
                             #noise_off=noise_off)
      ##print("C_vec_q[a] = ", C_vec_q[a])

      
## Calculated A matrix -----------------------------------------------
    #for b in range(a, Nparam):
      #dpsi_c_b = apply_danzatz(anzatz_tp, rot, Nlayers, rot_ang_c, Nq, 
                                #U, dU, b, psi_i_vec, Uent)

      #A_mtrx_c[a,b] = 2*(dpsi_c_a.conj().T).dot( dpsi_c_b ).item(0).real
      #A_mtrx_c[b,a] = A_mtrx_c[a,b]

      #if Q_on:
        #A_mtrx_q[a,b] = clc_a_mtrx(anzatz_tp, psi_i_circ, Uent_circ,
                                   #Nlayers, rot, rot_ang_q, a, b,
                                   #noise_off=noise_off)
        #A_mtrx_q[b,a] = A_mtrx_q[a,b]


## Calculate new values of angles ------------------------------------
  ##dtheta_c = solve_Ax_b(A_mtrx_c, C_vec_c)
  #dtheta_c = solve_Ax_b_L_curve(A_mtrx_c, C_vec_c)
  ##print("dtheta_c = ", dtheta_c - dtheta_tst)
  #rot_ang_c += dtau * dtheta_c
  


## Average value of Hamiltonian and wf measurements ------------------
  #U = anzatz_matrices(anzatz_tp, rot, Nlayers, rot_ang_c, Nq)
  #psi_c = apply_anzatz(anzatz_tp, rot, Nlayers, rot_ang_c, Nq,
                       #psi_i_vec, Uent)

  #E = (psi_c.conj().T).dot( H_mtrx.dot(psi_c) ).item().real

  #stps_left, converged, H_av = get_Nsteps_for_exp_tauH(psi_c, H_mtrx, exp_tau_H)
  ##print(np.amax(abs(psi_tauH - psi_c)))

## It will be nice to reduce the step (dtau), when 
## the local minima is found
  #if converged:
    #print("<E> = ", f'{E: .5f}', 
          #"|psi - psi_tauH| = ", 
          #f'{ np.linalg.norm(psi_c - psi_tauH)**2.0: .1e}',
          #"Steps left: ", stps_left, 
          #flush=True)
  #else:
    #print("<E> = ", f'{E: .5f}', 
          #"|psi - psi_tauH| = ", 
          #f'{ np.linalg.norm(psi_c - psi_tauH)**2.0: .1e}',
          #"Steps left: >", stps_left, 
          #flush=True)

  ##if it%100 == 0:
    ##print("<E> = ", f'{E: .5f}', 
          ##"|psi - psi_exct| = ", 
          ##f'{ np.linalg.norm(psi_c - psi_exact)**2.0: .1e}',
          ##np.amax(abs(dtheta_c)), flush=True)

  ##if np.amax(abs(dtheta_c)) < 1.e-5:
    ##break
  ##exit()


##print("Best rotation angles \n", rot_ang_bst, flush=True )
#print("Angles from imaginary time\n", rot_ang_c, flush=True )


## On output we have
#print( "\n", " "*Nq, f'{"exact": >7}', f'{"approx": >8}', flush=True )
#for i in range(2**Nq):
  #if abs(psi_c[i])**2 > 1.e-5 or abs(psi_exact[i])**2 > 1.e-5:
    #print( bin(i)[2:].zfill(Nq), 
          #f'{abs(psi_exact[i])**2: .5f}', 
          #f'{abs(psi_c[i])**2: .5f}', 
          #flush=True)
