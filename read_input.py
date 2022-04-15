import numpy as np

ansatzes = {
  "he":0,
  "ducc":2
  }

name2l = {
  "s": 0,
  "p": 1,
  "d": 2,
  "f": 3,
  "g": 4,
  "h": 5,
  "i": 6
  }

# ===================================================================
class ansatz_info_cls():
  def __init__(self, tp, L, rot_exc, split):
    self.tp = tp
    self.L = L
    self.rot_exc = rot_exc
    self.split = split
# ===================================================================


# ===================================================================
def state_to_nlj(state):
  for key in name2l.keys():
    if key in state:
      l = name2l[key]
      n,j = [int(x) for x in state.replace(key," ").split()]
      break

  return n, l, j
# ===================================================================


# ===================================================================
class rad_orb_cls:
  def __init__(self, n, l, j, indx=0):
    self.n = n
    self.l = l
    self.j = j
    self.indx = indx

    self.k = (-1)**int(l+j*0.5+0.5) * int(j*0.5 + 0.5)

    self.i = -1

  def __eq__(self, other):
    return self.k == other.k and self.n == other.n

  def __str__(self):
    return str(self.n) + " " + str(self.l) + " " + str(self.j)
# ===================================================================


# ===================================================================
def find_rad_orb(arr, n, l, j):
  for i,x in enumerate(arr):
    if x.n != n:
      continue

    if x.l != l:
      continue
    
    if x.j != j:
      continue
    
    return i
  return -1  
# ===================================================================


# ===================================================================
class orb_cls:
  def __init__(self, n, k, m):
    self.n = n
    self.k = k
    self.m = m

    self.i = -1

    self.j = int(2*abs(k) - 1)
    self.l = int( abs(k + 0.5) - 0.5 )

  def __eq__(self, other):
    return self.k == other.k and self.n == other.n and self.m == other.m

  def __str__(self):
    return str(self.n) + " " + str(self.k) + " " + str(self.m)
# ===================================================================


# ===================================================================
def find_orb(arr, n, l, j, m):
  for i,x in enumerate(arr):
    if x.n != n:
      continue

    if x.l != l:
      continue
    
    if x.j != j:
      continue
    
    if x.m != m:
      continue
    
    return i
  return -1  
# ===================================================================


# ===================================================================
def read(fl_name):
  #print("entered")
  fl_inp = open(fl_name,"r")

  basis = []
  read_basis = False

  for ln in fl_inp.readlines():
    if ln == "\n":
      continue
    
    if "#" in ln:
      continue

    if read_basis:
      if ("basis" in ln) and ("end" in ln):
        read_basis = False
        continue

      basis.append( ln.split()[0] )
      continue

    key, val = ln.strip().replace(" ","").split("=")

    if key == "OBI":
      fl_1b_int = open(val,"r")

    if key == "TBI":
      fl_2b_int = open(val,"r")

    if key == "Nelec":
      Nelec = int(val)

    if key == "ansatz_tp":
      ansatz_tp = ansatzes[val]

    if key == "Nlayers":
      Nlayers = int(val)

    if key == "rot_exc":
      rot_exc = [int(x) for x in val.replace("[","").replace("]","").split(",")]

    if key == "Q_on":
      Q_on = (val == "True")

    if key == "noise_off":
      noise_off = (val == "True")

    if key == "Nrep":
      Nrep = int(val)

    if key == "Parity":
      Parity_tot = int(val)
      
    if key == "Jz":
      Jz_tot = float(val)

    if key == "tappering":
      tappering_on = (val == "True")

    if key == "split":
      split_ent_layers = (val == "True")

    if key == "basis":
      if val == "begin":
        read_basis = True


# Fill in radial orbitals from basis - - - - - - - - - - - - - - - -
  #print("Basis")
  rad_orb_arr = []
  for b in basis:
    n, l, j = state_to_nlj(b)
    rad_orb_arr.append( rad_orb_cls(n, l, j) )
    #print(rad_orb_arr[-1])
  #exit()


# Sort orbitals with respect to parity
  rad_orb_arr = sorted(rad_orb_arr, 
                       key = lambda x: (-1)**x.l, 
                       reverse=False )


# Creating spin-orbitals - - - - - - - - - - - - - - - - - - - - - - 
  orb_arr = []

  norb = 0
  for x in rad_orb_arr:
    for m in range( 1, x.j+1, 2 ):
      for s in range(-1,2,2):
        orb_arr.append( orb_cls(x.n, x.k, m*s) )
        orb_arr[-1].i = x.i
  norb = len(orb_arr)


  #print("\nFollowing orbitals are used")
  iq_even_first = -1
  #print("Odd")
  for i,x in enumerate(orb_arr):
    if x.l%2 == 0 and iq_even_first == -1:
      iq_even_first = i
      #print("Even")
      
    #print("q"+str(i), x, x.l)
  #print()

  one_b_int = np.zeros([norb,norb])
  two_b_int = np.zeros([norb,norb,norb,norb])

  #exit()
# end creating spin-orbitals - - - - - - - - - - - - - - - - - - - - 


# Read one-body integrals - - - - - - - - - - - - - - - - - - - - - -
  fl_1b_int.seek(0)
  read_int = False
  for ln in fl_1b_int.readlines():
    if ln == "\n":
      continue
    
    if "na" in ln:
      read_int = True
      continue

    arr = ln.split()
    if read_int:
      na,la,ja = [int(arr[i]) for i in range(3)]
      nb,lb,jb = [int(arr[3+i]) for i in range(3)]

      h_ab = float( arr[6] )
      
      for a,x in enumerate(orb_arr):
        if x.n != na or x.l != la or x.j != ja:
          continue

        for b,y in enumerate(orb_arr):
          if y.m != x.m:
            continue

          if y.n != nb or y.l != lb or y.j != jb:
            continue

          one_b_int[a,b] = h_ab
          one_b_int[b,a] = h_ab
  #print("One-body integrals read")


# Read two-body integrals - - - - - - - - - - - - - - - - - - - - - -
  read_int = False
  for ln in fl_2b_int.readlines():
    if ln == "\n":
      continue
    
    if "na" in ln:
      read_int = True
      continue

    arr = ln.split()
    if read_int:
      n = [int( arr[ii] ) for ii in range(4)]
      l = [int( arr[4+ii] ) for ii in range(4)]
      j = [int( arr[8+ii] ) for ii in range(4)]
      m = [int( arr[12+ii] ) for ii in range(4)]

      u_coul_dir = float( arr[16] )
      u_coul_ex = float( arr[17] )
      u_br_dir = float( arr[18] )
      u_br_ex = float( arr[19] )

      p = find_orb(orb_arr, n[0], l[0], j[0], m[0])
      if p == -1:
        continue
      q = find_orb(orb_arr, n[1], l[1], j[1], m[1])
      if q == -1:
        continue
      r = find_orb(orb_arr, n[2], l[2], j[2], m[2])
      if r == -1:
        continue
      s = find_orb(orb_arr, n[3], l[3], j[3], m[3])
      if s == -1:
        continue

      # direct
      two_b_int[p,q,r,s] = u_coul_dir + u_br_dir
      two_b_int[q,p,s,r] = u_coul_dir + u_br_dir

      two_b_int[s,r,q,p] = u_coul_dir + u_br_dir
      two_b_int[r,s,p,q] = u_coul_dir + u_br_dir

      # exchange
      two_b_int[p,q,s,r] = -(u_coul_ex + u_br_ex)
      two_b_int[q,p,r,s] = -(u_coul_ex + u_br_ex)

      two_b_int[r,s,q,p] = -(u_coul_ex + u_br_ex)
      two_b_int[s,r,p,q] = -(u_coul_ex + u_br_ex)
  #print("Two-body integrals read")

  ansatz = ansatz_info_cls(ansatz_tp, Nlayers, rot_exc, split_ent_layers)

  return Nelec, Jz_tot, Parity_tot, tappering_on, orb_arr, one_b_int, two_b_int, ansatz
