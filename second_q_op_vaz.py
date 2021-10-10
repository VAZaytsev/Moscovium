# indexes of the creation and annihilation operators must go in the
# decreasing order, i.e.
# a_3 a^+_1 a_0
# permutation to such order of the initial form give a rise to a sign
# which is recorded in coef 

# ===================================================================
def second_q_op_pq(p, q, n):
  coef = 1

# p = q
  if p == q:
    tmp = "I"*p + "N" + "I"*(n - p - 1)
    return coef, tmp

# no equal pairs
  if p < q:
    coef *= -1

  tmp = ""
  for i in range(n):
    if i == p:
      tmp += "-"
      continue
    if i == q:
      tmp += "+"
      continue
    tmp += "I"

  return coef, tmp
# ===================================================================


# ===================================================================
def second_q_op_pqrs(p, q, r, s, n):
  coef = 1

# two equal pairs
  if p == r and q == s:
    coef = -1
    tmp = ""
    for i in range(n):
      if i == p or i == q:
        tmp += "N"
      else:
        tmp += "I"
    return coef, tmp

  if p == s and q == r:
    tmp = ""
    for i in range(n):
      if i == p or i == q:
        tmp += "N"
      else:
        tmp += "I"
    return coef, tmp

# one equal pair
  if p == r:
    coef *= -1
    if q < s:
      coef *= -1
    tmp = ""
    for i in range(n):
      if i == p:
        tmp += "N"
        continue
      if i == q:
        tmp += "-"
        continue
      if i == s:
        tmp += "+"
        continue
      tmp += "I"
    return coef, tmp

  if p == s:
    if q < r:
      coef *= -1
    tmp = ""
    for i in range(n):
      if i == p:
        tmp += "N"
        continue
      if i == q:
        tmp += "-"
        continue
      if i == r:
        tmp += "+"
        continue
      tmp += "I"
    return coef, tmp

  if q == r:
    if p < s:
      coef *= -1
    tmp = ""
    for i in range(n):
      if i == q:
        tmp += "N"
        continue
      if i == p:
        tmp += "-"
        continue
      if i == s:
        tmp += "+"
        continue
      tmp += "I"
    return coef, tmp

  if q == s:
    coef *= -1
    if p < r:
      coef *= -1
    tmp = ""
    for i in range(n):
      if i == q:
        tmp += "N"
        continue
      if i == p:
        tmp += "-"
        continue
      if i == r:
        tmp += "+"
        continue
      tmp += "I"
    return coef, tmp
  
# no equal pairs
  arr = [p, q, r, s]
  for i in range(3):
    swapped = False
    for j in range(0,3-i):
      if arr[j] > arr[j + 1]:
        arr[j], arr[j + 1] = arr[j + 1], arr[j]
        swapped = True
        coef *= -1
    if not swapped:
      break

  tmp = ""
  for i in range(n):
    if i == p or i == q:
      tmp += "-"
      continue
    if i == r or i == s:
      tmp += "+"
      continue
    tmp += "I"
  return coef, tmp
# ======================================================================
