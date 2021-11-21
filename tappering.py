from qiskit.opflow.primitive_ops.pauli_sum_op import PauliSumOp

def tapper(H, q, state):
  nq = H.num_qubits

  H_new = PauliSumOp.from_list([("I"*(nq-1),0.0)])
  for h in H.primitive.to_list():
    if h[0][nq-1-q] == "X" or h[0][nq-1-q] == "Y":
      continue
    
    lbl = h[0][:nq-1-q] + h[0][nq-q:]

    if h[0][nq-1-q] == "I":
      coef = 1
    else:
      coef = (-1)**state

    H_new += PauliSumOp.from_list([(lbl,h[1]*coef)])

  return H_new.reduce()
