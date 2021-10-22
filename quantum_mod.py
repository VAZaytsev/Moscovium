import math

import numpy as np
from numpy import pi
np.version.version

from qiskit import *
#from qiskit.providers.aer import QasmSimulator

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

#provider = IBMQ.load_account()

Nshots = 8192
Nrep = 1

err_stat = 2/math.sqrt(Nshots*Nrep)

backend = Aer.get_backend('qasm_simulator')

# For noise model
from qiskit.providers.aer.noise import NoiseModel, pauli_error
from qiskit.providers.aer.noise import depolarizing_error
from qiskit.transpiler.coupling import CouplingMap

coupling_map = None
basis_gates = ['cx', 'id', 'reset', 'rz', 'sx', 'x']
noise_model = NoiseModel(basis_gates)



# ===================================================================
def create_noise_model(nq, nrep):
  global noise_model
  global coupling_map
  global Nrep
  global err_stat
  
  Nrep = nrep
  err_stat = 2/math.sqrt(Nshots*Nrep)
  
# linear connectivity
  coupling_map = [[q,q+1] for q in range(nq-1)]

# future ancilla qubits
  for q in range(nq):
    coupling_map.append( [q,nq+q] )
    

# Measurement error
  p_meas = 1.126e-2
  error_meas = pauli_error([('X',p_meas), 
                            ('I', 1 - p_meas)])
  noise_model.add_all_qubit_quantum_error(error_meas, "measure")


# Depolarizing quantum errors
  prob_1 = 2.935e-4
  prob_2 = 7.925e-3

  error_1 = depolarizing_error(prob_1, 1)
  error_2 = depolarizing_error(prob_2, 2)

  noise_model.add_all_qubit_quantum_error(error_1, ['u1', 'u2', 'u3'])
  noise_model.add_all_qubit_quantum_error(error_2, ['cx'])

  return 
# ===================================================================




# ===================================================================
def entangling_circ(nq):
  qc = QuantumCircuit(nq)

  for q in range(nq-1):
    qc.cx(q,q+1)

  return qc
# ===================================================================




# ===================================================================
def ps_2_circ(ps):
  qc = QuantumCircuit( len(ps) )

  for i, s in enumerate(ps):
    if s == "1" or s == 1:
      qc.x(len(ps)-1-i)

  return qc
# ===================================================================




# ===================================================================
def ps2meas_basis(ps):
  nq = len(ps)
  circ = QuantumCircuit(nq,nq)

  for i in range(nq):
    if ps[i] == "1" or ps[i] == "X":
      circ.h(nq-1-i)
      continue

    if ps[i] == "2" or ps[i] == "Y":
      circ.rx(-0.5*pi,nq-1-i)
      continue

  return circ
# ===================================================================




# ===================================================================
def extract_ps_for_measurement(H):
  ps_meas = []
  indx_ps2meas = np.zeros( len(H.primitive.to_list()), dtype=int )

  for i,x in enumerate(H.primitive.to_list()):
    ps_tmp = x[0].replace("Z","I")

    found = False
    for ii,ps in enumerate(ps_meas):
      if ps_tmp == ps:
        found = True
        break

    if found:
      indx_ps2meas[i] = ii
    else:
      ps_meas.append(ps_tmp)
      indx_ps2meas[i] = len(ps_meas)-1

  return ps_meas, indx_ps2meas
# ===================================================================




# ===================================================================
def anzatz_circ(anzatz_tp, nq, ent_circ, nlayers, rot, angs):

#  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
  if anzatz_tp == 0:
    qc = QuantumCircuit(nq)

    ang_per_layer = len(rot) * nq

    for ilayer in range( nlayers+1 ):
      for ir,r in enumerate(rot):
        for iq in range( nq ):
          i_ang = ilayer * ang_per_layer + ir*nq + iq

          if r == 1:
            qc.rx( angs[i_ang], iq )
          if r == 2:
            qc.ry( angs[i_ang], iq )
          if r == 3:
            qc.rz( angs[i_ang], iq )

      if ilayer != nlayers:
        qc = qc.compose(ent_circ, [i for i in range(nq)])

    return qc
  

#  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
  if anzatz_tp == 2:
    qc = QuantumCircuit(nq)

    qc.ry(angs[0], 0)

    for ctrl in range(6):
      qc.cx(ctrl, ctrl+1)

    qc.ry(0.5*angs[1], 1)
    qc.cx(0, 1)
    qc.ry(-0.5*angs[1], 1)
    qc.cx(0, 1)

    qc.cx(2, 3)

    qc.ry( 0.5*angs[2], 5)
    qc.cx(4, 5)
    qc.ry(-0.5*angs[2], 5)
    qc.cx(4, 5)

    qc.x(3)
      
    qc.x(1)
    qc.cx(1, 2)
    qc.x(1)
    qc.x(2)
    qc.cx(2, 3)
    qc.x(3)
    qc.cx(3, 4)

    qc.ry( 0.5*angs[3], 5)
    qc.cx(4, 5)
    qc.ry(-0.5*angs[3], 5)
    qc.cx(4, 5)

    qc.cx(3, 4)
    qc.x(3)
    qc.cx(2, 3)
    qc.x(2)

    qc.x(5)
    qc.cx(5, 6)
    qc.x(5)
    
    return qc
# ===================================================================




# ===================================================================
def c_vec_circ(anzatz_tp, psi_i_circ, ent_circ, 
               nlayers, rot, angs, a):
  nq = psi_i_circ.num_qubits

  qc_arr = []
  coef_arr = []
  targ_arr = []
  coef = 1
  targ = 0

#  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
  if anzatz_tp == 0:
    qc = QuantumCircuit(nq+1, nq+1).compose(psi_i_circ, 
                                            [i for i in range(nq)])
    qc.h(nq)
    qc.s(nq)

    ang_per_layer = len(rot) * nq

    for ilayer in range( nlayers+1 ):
      for ir,r in enumerate(rot):
        for iq in range( nq ):
          i_ang = ilayer * ang_per_layer + ir*nq + iq

          if i_ang == a:
            targ = iq
            if r == 1:
              qc.cx(nq, iq)
            if r == 2:
              qc.cy(nq, iq)
            if r == 3:
              qc.cz(nq, iq)

            qc.x(nq)
            qc.h(nq)

          if r == 1:
            qc.rx( angs[i_ang], iq )
          if r == 2:
            qc.ry( angs[i_ang], iq )
          if r == 3:
            qc.rz( angs[i_ang], iq )

      if ilayer != nlayers:
        qc = qc.compose(ent_circ, [i for i in range(nq)])

    return [targ], [1.0], [qc]
  

#  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
  if anzatz_tp == 2:
    for icounter in range(2):
      qc = QuantumCircuit(nq+1, nq+1).compose(psi_i_circ, 
                                              [i for i in range(nq)])
      qc.h(nq)
      qc.s(nq)

      if a == 0:
        coef = 1
        targ = 0
        qc.cy(nq, 0)
      qc.ry(angs[0], 0)


      for ctrl in range(6):
        qc.cx(ctrl, ctrl+1)

      if a == 1 and icounter == 0:
        coef = 0.5
        targ = 1
        qc.cy(nq, 1)
      qc.ry(0.5*angs[1], 1)
      qc.cx(0, 1)
      if a == 1 and icounter == 1:
        coef = -0.5
        targ = 1
        qc.cy(nq, 1)
      qc.ry(-0.5*angs[1], 1)
      qc.cx(0, 1)

      qc.cx(2, 3)

      if a == 2 and icounter == 0:
        coef = 0.5
        targ = 5
        qc.cy(nq, 5)
      qc.ry( 0.5*angs[2], 5)
      qc.cx(4, 5)
      if a == 2 and icounter == 1:
        coef = -0.5
        targ = 5
        qc.cy(nq, 5)
      qc.ry(-0.5*angs[2], 5)
      qc.cx(4, 5)

      qc.x(3)
      
      qc.x(1)
      qc.cx(1, 2)
      qc.x(1)
      qc.x(2)
      qc.cx(2, 3)
      qc.x(3)
      qc.cx(3, 4)

      if a == 3 and icounter == 0:
        coef = 0.5
        targ = 5
        qc.cy(nq, 5)
      qc.ry( 0.5*angs[3], 5)
      qc.cx(4, 5)
      if a == 3 and icounter == 1:
        coef = -0.5
        targ = 5
        qc.cy(nq, 5)
      qc.ry(-0.5*angs[3], 5)
      qc.cx(4, 5)

      qc.cx(3, 4)
      qc.x(3)
      qc.cx(2, 3)
      qc.x(2)
      
      qc.x(5)
      qc.cx(5, 6)
      qc.x(5)
      
      qc.x(nq)
      qc.h(nq)

      qc_arr.append(qc)
      coef_arr.append(coef)
      targ_arr.append(targ)

      if a == 0:
        break
    
    return targ_arr, coef_arr, qc_arr
# ===================================================================




# ===================================================================
def a_mtrx_circ(anzatz_tp, psi_i_circ, ent_circ, 
                nlayers, rot, angs, a, b):
# it is assumed that a <= b 
  nq = psi_i_circ.num_qubits

  qc_arr = []
  coef_arr = []
  targ_arr = []

#  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
  if anzatz_tp == 0:
    ang_per_layer = len(rot) * nq

    layer_a = math.floor(a / ang_per_layer)
    ir_a = math.floor((a  - layer_a*ang_per_layer) / float(nq))
    targ_a = a - layer_a * ang_per_layer - ir_a * nq

    layer_b = math.floor(b / ang_per_layer)
    ir_b = math.floor((b  - layer_b*ang_per_layer) / float(nq))
    targ_b = b - layer_b * ang_per_layer - ir_b * nq

    if targ_a == targ_b:
      qc = QuantumCircuit(nq+1, nq+1).compose(psi_i_circ, 
                                              [i for i in range(nq)])
    else:
      qc = QuantumCircuit(nq+2, nq+2).compose(psi_i_circ, 
                                              [i for i in range(nq)])

    qc.h(nq)

    global_break = False
    for ilayer in range( nlayers+1 ):
      for ir,r in enumerate(rot):
        for iq in range( nq ):
          i_ang = ilayer * ang_per_layer + ir*nq + iq

          if i_ang == a:
            if r == 1:
              qc.cx(nq, iq)
            if r == 2:
              qc.cy(nq, iq)
            if r == 3:
              qc.cz(nq, iq)

          if i_ang == b:
            if targ_a == targ_b:
              qc.x(nq)
              ctrl = nq
            else:
              qc.h(nq)
              ctrl = nq+1
              qc.h(ctrl)

            if r == 1:
              qc.cx(ctrl, iq)
            if r == 2:
              qc.cy(ctrl, iq)
            if r == 3:
              qc.cz(ctrl, iq)
            qc.h(ctrl)

            global_break = True
            break
              
          cond1 = (ilayer == layer_b) and iq == targ_a and ir < ir_a
          cond2 = (ilayer == layer_b) and iq == targ_b
          cond3 = (ilayer != layer_b)
          if cond1 or cond2 or cond3:
            if r == 1:
              qc.rx( angs[i_ang], iq )
            if r == 2:
              qc.ry( angs[i_ang], iq )
            if r == 3:
              qc.rz( angs[i_ang], iq )

        if global_break:
          break
        
      if global_break:
        break
      
      if ilayer != nlayers:
        qc = qc.compose(ent_circ, [i for i in range(nq)])

    return [[targ_a, targ_b]], [1.0], [qc]


#  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
  if anzatz_tp == 2:
    targ_a = {0:0, 1:1, 2:5, 3:5}[a]
    targ_b = {0:0, 1:1, 2:5, 3:5}[b]

    ncounter_a = {0:1, 1:2, 2:2, 3:2}[a]
    ncounter_b = {0:1, 1:2, 2:2, 3:2}[b]

    ctrl_a = nq
    if a == b or (a == 2 and b == 3):
      na = 1
      ctrl_b = nq
    else:
      na = 2
      ctrl_b = nq+1

    for icounter_a in range(ncounter_a):
      for icounter_b in range(ncounter_b):

        if na == 2:
          qc = QuantumCircuit(nq+2, nq+2).compose(psi_i_circ, 
                                                  [i for i in range(nq)])
          qc.h(nq)
          qc.h(nq+1)
        else:
          qc = QuantumCircuit(nq+1, nq+1).compose(psi_i_circ, 
                                                  [i for i in range(nq)])
          qc_empty = QuantumCircuit(nq+1, nq+1).compose(psi_i_circ, 
                                                  [i for i in range(nq)])
          qc.h(nq)

        derivs = 0

        if a == 0:
          coef_a = 1
          if na == 1 and derivs == 1:
            qc.x(nq)
          qc.cy(ctrl_a, 0)
          derivs += 1
        if b == 0:
          coef_b = 1
          if na == 1 and derivs == 1:
            qc.x(nq)
          qc.cy(ctrl_b, 0)
          derivs += 1

        if derivs == 2:
          qc.h(nq)
          if na == 2:
            qc.h(nq+1)
      
          if a == b and icounter_a == icounter_b:
            qc_arr.append(qc_empty)
          else:
            qc_arr.append(qc)
          coef_arr.append(coef_a*coef_b)
          targ_arr.append([targ_a,targ_b])
          continue

        qc.ry(angs[0], 0)
        qc.cx(0,1)
        qc.cx(1,2)
        
        # Rotation on theta_1
        if a == 1 and icounter_a == 0:
          coef_a = 0.5
          if na == 1 and derivs == 1:
            qc.x(nq)
          qc.cy(ctrl_a, 1)
          derivs += 1
        if b == 1 and icounter_b == 0:
          coef_b = 0.5
          if na == 1 and derivs == 1:
            qc.x(nq)
          qc.cy(ctrl_b, 1)
          derivs += 1

        if derivs == 2:
          qc.h(nq)
          if na == 2:
            qc.h(nq+1)
          if a == b and icounter_a == icounter_b:
            qc_arr.append(qc_empty)
          else:
            qc_arr.append(qc)
          coef_arr.append(coef_a*coef_b)
          targ_arr.append([targ_a,targ_b])
          continue

        qc.ry(0.5*angs[1], 1)
        qc.cx(0, 1)
        if a == 1 and icounter_a == 1:
          coef_a = -0.5
          if na == 1 and derivs == 1:
            qc.x(nq)
          qc.cy(ctrl_a, 1)
          derivs += 1
        if b == 1 and icounter_b == 1:
          coef_b = -0.5
          if na == 1 and derivs == 1:
            qc.x(nq)
          qc.cy(ctrl_b, 1)
          derivs += 1

        if derivs == 2:
          qc.h(nq)
          if na == 2:
            qc.h(nq+1)
          if a == b and icounter_a == icounter_b:
            qc_arr.append(qc_empty)
          else:
            qc_arr.append(qc)
          coef_arr.append(coef_a*coef_b)
          targ_arr.append([targ_a,targ_b])
          continue
        qc.ry(-0.5*angs[1], 1)


        qc.cx(2,3)
        qc.cx(3,4)
        qc.cx(4,5)
        qc.cx(5,6)


        # Rotation on theta_2
        if a == 2 and icounter_a == 0:
          coef_a = 0.5
          if na == 1 and derivs == 1:
            qc.x(nq)
          qc.cy(ctrl_a, 5)
          derivs += 1
        if b == 2 and icounter_b == 0:
          coef_b = 0.5
          if na == 1 and derivs == 1:
            qc.x(nq)
          qc.cy(ctrl_b, 5)
          derivs += 1

        if derivs == 2:
          qc.h(nq)
          if na == 2:
            qc.h(nq+1)
          if a == b and icounter_a == icounter_b:
            qc_arr.append(qc_empty)
          else:
            qc_arr.append(qc)
          coef_arr.append(coef_a*coef_b)
          targ_arr.append([targ_a,targ_b])
          continue

        qc.ry( 0.5*angs[2], 5)
        qc.cx(4, 5)
        if a == 2 and icounter_a == 1:
          coef_a = -0.5
          if na == 1 and derivs == 1:
            qc.x(nq)
          qc.cy(ctrl_a, 5)
          derivs += 1
        if b == 2 and icounter_b == 1:
          coef_b = -0.5
          if na == 1 and derivs == 1:
            qc.x(nq)
          qc.cy(ctrl_b, 5)
          derivs += 1

        if derivs == 2:
          qc.h(nq)
          if na == 2:
            qc.h(nq+1)
          if a == b and icounter_a == icounter_b:
            qc_arr.append(qc_empty)
          else:
            qc_arr.append(qc)
          coef_arr.append(coef_a*coef_b)
          targ_arr.append([targ_a,targ_b])
          continue

        qc.ry(-0.5*angs[2], 5)
        qc.cx(4, 5)

        # Rotation on theta_3
        if a == 3 and icounter_a == 0:
          coef_a = 0.5
          if na == 1 and derivs == 1:
            qc.x(nq)
          qc.cy(ctrl_a, 5)
          derivs += 1
        if b == 3 and icounter_b == 0:
          coef_b = 0.5
          if na == 1 and derivs == 1:
            qc.x(nq)
          qc.cy(ctrl_b, 5)
          derivs += 1

        if derivs == 2:
          qc.h(nq)
          if na == 2:
            qc.h(nq+1)
          if a == b and icounter_a == icounter_b:
            qc_arr.append(qc_empty)
          else:
            qc_arr.append(qc)
          coef_arr.append(coef_a*coef_b)
          targ_arr.append([targ_a,targ_b])
          continue

        qc.ry( 0.5*angs[3], 5)
        
        qc.cx(0, 1)
        qc.x(1)
        qc.cx(2, 3)
        qc.x(3)
        qc.cx(1, 2)
        qc.x(2)
        qc.cx(2, 3)
        qc.x(3)
        qc.cx(3, 4)
        qc.cx(4, 5)

        if a == 3 and icounter_a == 1:
          coef_a = -0.5
          if na == 1 and derivs == 1:
            qc.x(nq)
          qc.cy(ctrl_a, 5)
          derivs += 1
        if b == 3 and icounter_b == 1:
          coef_b = -0.5
          if na == 1 and derivs == 1:
            qc.x(nq)
          qc.cy(ctrl_b, 5)
          derivs += 1

        if derivs == 2:
          qc.h(nq)
          if na == 2:
            qc.h(nq+1)
          if a == b and icounter_a == icounter_b:
            qc_arr.append(qc_empty)
          else:
            qc_arr.append(qc)
          coef_arr.append(coef_a*coef_b)
          targ_arr.append([targ_a,targ_b])
          continue

    return targ_arr, coef_arr, qc_arr
# ===================================================================




# ===================================================================
def measure_ps(circ_in, ps, 
               noise_off=True, layout=None, real_QC=False):
  if circ_in.num_clbits == 0:
    cr = ClassicalRegister(len(ps), 'c')
    circ_in.add_register(cr)

  circ = circ_in.compose( ps2meas_basis(ps), 
                         [i for i in range(len(ps))] )

  circ.measure_all()
  #for q in range(circ_in.num_qubits):
    #circ.measure(q,q)

  # Simulator without noise
  if noise_off:
    result = execute(circ, backend, shots=Nshots*Nrep).result()
    counts = result.get_counts(circ)
  else:

    # Simulator with noise
    if not real_QC:
      tcirc = transpile(circ,
                        backend=backend,
                        coupling_map=coupling_map,
                        initial_layout=layout,
                        basis_gates=basis_gates)
      #print(tcirc.depth())
      #print(circ.num_nonlocal_gates(), tcirc.num_nonlocal_gates())

      result = execute(tcirc, backend,
                       coupling_map=coupling_map,
                       basis_gates=basis_gates,
                       noise_model=noise_model,
                       shots=Nshots*Nrep).result()

      counts = result.get_counts(circ)

# Real Quantum Computer
    else:
      counts = {}
      for irep in range(Nrep):
        result = execute(circ, backend_QC, shots=Nshots).result()
        counts_new = result.get_counts()

        for k in counts_new.keys():
          if k in counts.keys():
            counts[k] += counts_new[k]
          else:
            counts[k] = counts_new[k]

  for k in counts.keys():
    counts[k] /= Nshots*Nrep

  return counts
# ===================================================================




# ===================================================================
def measure_Ham(anzatz_tp, psi_i_circ, ent_circ,
                nlayers, rot, angs, 
                ps_meas, indx_ps2meas, H, 
                noise_off=True):
  res = 0.0

  nq = psi_i_circ.num_qubits

  circ = anzatz_circ(anzatz_tp, nq, ent_circ, nlayers, rot, angs)
  qc = psi_i_circ.compose(circ)

  ps_meas_data = [None]*len(ps_meas)


  # Measure all necessary Pauli strings
  for i,x in enumerate(ps_meas):
    ps_meas_data[i] = measure_ps(qc, x, noise_off=noise_off)


  # Calculate Hamiltonian
  for i,x in enumerate(H.primitive.to_list()):
    dct = ps_meas_data[ indx_ps2meas[i] ]
    for k in dct.keys():
      pwr = sum( [int(k[ii]) if x[0][ii] != "I" else 
                  0 for ii in range(nq)] )

      res += (-1)**pwr * x[1].real * dct[k]
        
  return res
# ===================================================================




# ===================================================================
def clc_c_vec(anzatz_tp, psi_i_circ, ent_circ,
              nlayers, rot, angs, a, 
              ps_meas, indx_ps2meas, H, 
              noise_off=True):
  res = 0.0

  nq = psi_i_circ.num_qubits

  targ_a, coef, c_circ = c_vec_circ(anzatz_tp, psi_i_circ, ent_circ,
                                    nlayers, rot, angs, a)

  for indx, circ in enumerate(c_circ):
    ps_meas_data = [None]*len(ps_meas)

    layout = [i for i in range(nq)]
    layout.append(targ_a[indx] + nq)

    # Measure all necessary Pauli strings
    for i,x in enumerate(ps_meas):
      ps_meas_data[i] = measure_ps(circ, x,
                                   noise_off=noise_off,
                                   layout=layout)

    # Calculate C_vec 
    for i,x in enumerate(H.primitive.to_list()):
      dct = ps_meas_data[ indx_ps2meas[i] ]
      for k in dct.keys():
        if k[0] == "1":
          pwr = 0
        else:
          pwr = 1
        pwr += sum( [int(k[1+ii]) if x[0][ii] != "I" else 
                     0 for ii in range(nq)] )

        res += coef[indx] * (-1)**pwr * x[1].real * dct[k]

  return res
# ===================================================================




# ===================================================================
def clc_a_mtrx(anzatz_tp, psi_i_circ, ent_circ,
               nlayers, rot, angs, a, b, 
               noise_off=True):
  res = 0.0

  nq = psi_i_circ.num_qubits

  targ, coef, a_circ = a_mtrx_circ(anzatz_tp, psi_i_circ,
                                               ent_circ, 
                                               nlayers, rot, angs, a, b)

  for indx, circ in enumerate(a_circ):
    layout = [i for i in range(nq)]
    layout.append(targ[indx][0] + nq)
    

    na = 1
    if targ[indx][0] != targ[indx][1]:
      layout.append(targ[indx][1] + nq)
      na = 2

    # Measure
    counts = measure_ps(circ, "I"*(nq+na),
                        noise_off=True, 
                        layout=layout)

    # Calculate A_mtrx
    if na == 1:
      p0 = sum( [counts[k] if k[0] == "0" else 
                 0 for k in counts.keys()] )

      res += coef[indx]*(p0 - 0.5)

    else:
      p00 = sum([counts[k] if k[:2] == "00" else 
                 0 for k in counts.keys()])
      p11 = sum([counts[k] if k[:2] == "11" else 
                 0 for k in counts.keys()])

      res += coef[indx]*(p00 + p11 - 0.5)

  return res
# ===================================================================
