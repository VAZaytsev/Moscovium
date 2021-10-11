import numpy as np

kp_to_name = {
  -1: "s",
   1: "p1", 
  -2: "p3",
   2: "d3",
  -3: "d5"
  }


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
def bin_to_conf( orb_arr, bn ):
  sz = len(orb_arr)
  
  n_arr = []
  k_arr = []
  q_arr = []

# Construct an array of radial orbitals  
  for x in orb_arr:
    Found = False
    for i in range(len(n_arr)):
      if x.n == n_arr[i] and x.k == k_arr[i]:
        Found = True
        break

    if not Found:
      n_arr.append(x.n)
      k_arr.append(x.k)
      q_arr.append(0)


# Assign the ocupation numbers to the radial orbitals
  for i in range(sz):
    x = bn[sz-1-i]
    if x == "0":
      continue

    for ii in range(len(n_arr)):
      if orb_arr[i].n == n_arr[ii] and orb_arr[i].k == k_arr[ii]:
        q_arr[ii] += 1
        break

  conf = ""
  for i in range(len(n_arr)):
    if q_arr[i] == 0:
      continue
    conf += str(n_arr[i]) \
    + kp_to_name[k_arr[i]] \
    + "^" + str(q_arr[i]) \
    + " "

  return conf
# ===================================================================
