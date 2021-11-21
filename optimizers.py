import numpy as np
from find_a import *

# -------------------------------------------------------------------
def update_angles(ansatz, args, optimizer, psi_in, H_mtrx):
  def fun(angs):
    psi = ansatz.act_on_vctr(angs, psi_in)

    E = (psi.conj().T).dot( H_mtrx.dot(psi) ).item().real
    return E


  def dfun(angs):
    Nparam = len(angs)
    res = np.zeros((Nparam),dtype=float)

    psi = ansatz.act_on_vctr(angs, psi_in)
    dpsi = ansatz.dact_on_vctr(angs, psi_in)

    for a in range(Nparam):
      res[a] = 2 * (dpsi[a].conj().T).dot( H_mtrx.dot(psi) ).item().real
    return res


  def QuantumGemetricTensor(angs):
    Nparam = len(angs)
    psi = ansatz.act_on_vctr(angs, psi_in)
    dpsi = ansatz.dact_on_vctr(angs, psi_in)

    QGT = np.zeros((Nparam,Nparam), dtype = complex)
    for row in range(Nparam):
      for col in range(row,Nparam):
        QGT[row,col] = (dpsi[row].conj().T).dot( dpsi[col] ).item() \
          - (dpsi[row].conj().T).dot( psi ) * (psi.conj().T).dot( dpsi[col] )

      QGT[col,row] = QGT[row,col]
    return QGT.real


  def QuantumGemetricTensor_ITE(angs):
    Nparam = len(angs)
    psi = ansatz.act_on_vctr(angs, psi_in)
    dpsi = ansatz.dact_on_vctr(angs, psi_in)

    QGT_ITE = np.zeros((Nparam,Nparam))
    for row in range(Nparam):
      for col in range(row,Nparam):
        QGT_ITE[row,col] = (dpsi[row].conj().T).dot( dpsi[col] ).item().real
        QGT_ITE[col,row] = QGT_ITE[row,col]
    return QGT_ITE.real


  if isinstance(optimizer, Adam_cls):
    dx = dfun(args)
    return optimizer.update(args, dx)


  if isinstance(optimizer, NatGrad_cls):
    QGT_mtrx = QuantumGemetricTensor(args)
    RHS_vec = -optimizer.eta * dfun(args)

    return optimizer.update(args, QGT_mtrx, RHS_vec)


  if isinstance(optimizer, ITE_cls):
    QGT_ITE_mtrx = QuantumGemetricTensor(args)
    RHS_vec = -optimizer.eta * dfun(args)

    return optimizer.update(args, QGT_ITE_mtrx, RHS_vec)
# -------------------------------------------------------------------


# -------------------------------------------------------------------
class Adam_cls():
  def __init__(self, n, 
               eta=0.01, 
               beta1=0.9, 
               beta2=0.999, 
               epsilon=1.e-8):
    self.name = "Adam"

    self.m = np.zeros(n)
    self.v = np.zeros(n)

    self.eta = eta # learning rate
    self.beta1 = beta1
    self.beta2 = beta2
    self.epsilon = epsilon
    
    self.t = 1
    self.f = None
    self.converged = False
    

  def update(self, x, dx):
    # biased first momentum estimate
    self.m = self.beta1*self.m + (1-self.beta1)*dx

    # biased second raw momentum estimate
    self.v = self.beta2*self.v + (1-self.beta2)*(dx**2)

    # bias-corrected first momentum estimate
    m_corr = self.m/(1-self.beta1**self.t)
    
    # bias-corrected second raw momentum estimate
    v_corr = self.v/(1-self.beta2**self.t)

    # update
    self.t += 1
    
    dx = self.eta*(m_corr/(np.sqrt(v_corr) + self.epsilon))

    return x - dx
# -------------------------------------------------------------------


# -------------------------------------------------------------------
class NatGrad_cls():
  def __init__(self, n, eta=0.01):
    self.eta = eta # learning rate
    self.name = " NG "

    self.converged = False
    self.t = 1
    self.f = None

  def update(self, x, A_mtrx, b_vec):
    self.t += 1

    dx = solve_Ax_b_L_curve(A_mtrx, b_vec, lmbd_min=1.e-3, lmbd_max=1.e-1)
    #dx = solve_Ax_b(A_mtrx, b_vec)

    return x + dx
# -------------------------------------------------------------------


# -------------------------------------------------------------------
class ITE_cls():
  def __init__(self, n, eta=0.01):
    self.eta = eta # learning rate
    self.name = " ITE"
        
    self.converged = False
    self.t = 1
    self.f = None

  def update(self, x, A_mtrx, b_vec):
    self.t += 1

    dx = solve_Ax_b_L_curve(A_mtrx, b_vec, lmbd_min=1.e-3, lmbd_max=1.e-1)
    #dx = solve_Ax_b(A_mtrx, b_vec)

    return x + dx
# -------------------------------------------------------------------
