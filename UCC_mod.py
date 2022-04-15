import numpy as np

import itertools

from qiskit_nature.operators.second_quantization import FermionicOp
from qiskit_nature.converters.second_quantization import QubitConverter
from qiskit.opflow.primitive_ops.pauli_sum_op import PauliSumOp

# ===================================================================
class ClusterOperator():
  def __init__(self, ann, crt):
    self.ann = ann
    self.crt = crt
    self.nex = len(crt)
    self.q = []


  def fer_op(self, nq):
    T = FermionicOp("", register_length=nq)
    for o in self.crt:
      T @= FermionicOp("+_"+str(o), register_length=nq)
    for o in self.ann:
      T @= FermionicOp("-_"+str(o), register_length=nq)
    T -= ~T
    self.T = T.reduce()


  def q_op(self, q_converter, num_part):
    self.Tq = q_converter.convert(self.T, num_particles=num_part)


  def mtrx(self):
    # the matrix is real and antihermitian (T^+ = -T)
    self.mtrx = self.Tq.to_matrix().real

    e_val, e_vec = np.linalg.eig(self.mtrx)

    self.e_val = [-1j, 0, 1j]
    self.P = [ np.zeros(e_vec.shape, dtype=complex) ] * 3
    for i,e in enumerate(e_val):
      for j,val in enumerate(self.e_val):
        if abs(e - val) < 1.e-13:
          out = np.outer( e_vec[:,i], np.conj(e_vec[:,i]) )
          self.P[j] = np.add(self.P[j],
                             np.outer( e_vec[:,i], 
                                      np.conj(e_vec[:,i]) 
                                      )
                             )
          break


  def q_act_on(self):
    nq = self.Tq.num_qubits

    for x in self.Tq.primitive.to_list():
      for i, ps in enumerate(x[0]):
        if ps != "I" and not nq-1-i in self.q:
          self.q.append(nq-1-i)

    self.q = sorted(self.q, reverse=True)
          
# ===================================================================


# ===================================================================
def create_cluster_operators(psi_bn_p_tappered, 
             Nelec, 
             Parity_tot, 
             iq_even_first, 
             orbs):
  nq = len(psi_bn_p_tappered)
  # restore tappered qubits
  #print(psi_bn_p_tappered,"tappered")
  
  psi_bn_p = str(Nelec%2) \
                + psi_bn_p_tappered[:nq-iq_even_first+1] \
                + str(Parity_tot) \
                + psi_bn_p_tappered[nq-iq_even_first+1:]
  nq += 2
  #print("nq = ", nq)
        
  #print(psi_bn_p,"parity")
  
  # convert to Jordan-Wigner mapping
  psi_bn_jw = [str((int(psi_bn_p[i])+int(psi_bn_p[i+1]))%2) 
               for i in range(len(psi_bn_p)-1)]
  psi_bn_jw.append( psi_bn_p[-1] )
  psi_bn_jw = "".join(psi_bn_jw)
  #print(psi_bn_jw,"jw")

  # occupied orbitals
  occ_in = [nq-1-i for i,x in enumerate(psi_bn_jw) if x == "1"]
  #print("occ = ", occ_in)
  
  # vacant orbitals
  vac = [orb for orb in range(nq) if not orb in occ_in]
  #print("vac = ", vac)

  cluster_ops = []
  # Evangelista order is used to excite from the reference state
  # The parity and momentum projection are conserved
  # nes if the deepest active electron
  for nes in range(1,Nelec+1):
    for nex in range(nes,0,-1):
      if nex > len(vac):
        continue

      for add in list(itertools.combinations(occ_in[:nes-1],nex-1)):
        ann = sorted([occ_in[nes-1]] + list(add))
        
        for crt in list(itertools.combinations(vac,nex)):
          m_old = sum( [orbs[a].m for a in ann] )
          m_new = sum( [orbs[a].m for a in crt] )
          if m_old != m_new:
            continue

          p_old = sum( [orbs[a].l for a in ann] )%2
          p_new = sum( [orbs[a].l for a in crt] )%2
          if p_old != p_new:
            continue

          cluster_ops.append( ClusterOperator(list(crt), ann) )
          cluster_ops[-1].fer_op(nq)

  return cluster_ops
# ===================================================================


# ===================================================================
def s_cluster_operators(nq):
  cluster_ops = []

  # Singles of only closest neighbours
  for i in range(nq-1):
    ann = [i]
    crt = [i+1]

    cluster_ops.append( ClusterOperator(crt, ann) )
    cluster_ops[-1].fer_op(nq)
    
  # Doubles 1100-0011
  #for i in range(nq):
    #if i+3 >= nq:
      #continue
    #ann = [i,i+1]
    #crt = [i+2,i+3]

    #cluster_ops.append( ClusterOperator(crt, ann) )
    #cluster_ops[-1].fer_op(nq)

  # Doubles 1010-0101
  for i in range(nq):
    if i+3 >= nq:
      continue
    ann = [i,i+2]
    crt = [i+1,i+3]

    cluster_ops.append( ClusterOperator(crt, ann) )
    cluster_ops[-1].fer_op(nq)

  # Doubles 0110-1001
  #for i in range(nq):
    #if i+3 >= nq:
      #continue
    #ann = [i,i+3]
    #crt = [i+1,i+2]

    #cluster_ops.append( ClusterOperator(crt, ann) )
    #cluster_ops[-1].fer_op(nq)

  return cluster_ops
# ===================================================================
