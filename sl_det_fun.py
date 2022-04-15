from qiskit_nature.mappers.second_quantization import ParityMapper
import itertools

# ===================================================================
class sl_det_cls:
  def __init__(self, bn, n, p, j):
    self.bn_arr = [bn]
    self.n = n # number of particles
    self.p = p # parity
    self.j = j # total angular momentum
# ===================================================================


# ===================================================================
def get_jz( orb_arr, bn ):
  jz = 0
  for i,x in enumerate(bn):
    jz += orb_arr[len(bn)-1-i].m * int(x)
    
  return jz
# ===================================================================


# ===================================================================
def get_parity( orb_arr, bn ):
  parity = 0
  for i,x in enumerate(bn):
    parity += orb_arr[len(bn)-1-i].l * int(x)
    
  return parity%2
# ===================================================================


# ===================================================================
def sl_dets_with_n_j_p(orb_arr, n, j, p, qubit_conv, tapper):
  SD_sym = []

  if tapper:
    iq_even_first = -1
    for i,x in enumerate(orb_arr):
      if x.l%2 == 0 and iq_even_first == -1:
        iq_even_first = i
        break
    

# Jordan Wigner mapping
  nq = len(orb_arr)
  for pos in list(itertools.combinations(range(nq),n)):
    bn = ''.join(['1' if i in pos else '0' for i in range(nq)])

    if get_jz( orb_arr, bn ) != 2*j:
      continue
  
    if get_parity( orb_arr, bn ) != p:
      continue
  
    SD_sym.append(bn)

# Convert to Parity mapping if needed
  if isinstance(qubit_conv.mapper,type(ParityMapper())):
    for i in range(len(SD_sym)):
      SD_sym[i] = ''.join([str(SD_sym[i][ii:].count('1')%2) 
                           for ii in range(nq)])
      if tapper:
        SD_sym[i] = SD_sym[i][1:nq-iq_even_first] + SD_sym[i][nq-iq_even_first+1:]

  return SD_sym
