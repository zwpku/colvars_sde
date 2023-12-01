#!/usr/bin/env python
import time
from tqdm import tqdm
import random
import numpy as np

def set_all_seeds(seed):
    torch.manual_seed(seed)
    random.seed(seed)

class StiffPot : 
    def __init__(self):
        self.eps = 0.5
        self.x0 = [-1, 0]

    def init_path(self, N):
        return np.array(np.linspace(self.min_A, self.min_B, N))

    def V(self, X):
      return (X[0]**2 - 1)**2 + 1.0 / self.eps * (X[0]**2 + X[1] - 1)**2

    def gradV(self, X): 
      return np.array(( 4.0 * X[0] * (X[0]**2 - 1.0 + 1.0 / self.eps * (X[0]**2 + X[1] - 1)), 2.0 / self.eps * (X[0]**2 + X[1] - 1)) )

def UnbiasedTraj(pot, X_0, beta, delta_t=1e-3, N=1000, seed=0):

    r = np.random.RandomState(seed)
    X = X_0
    dim = X.shape[0]
    for i in range(N):
        #b = np.array([0.1, 0.1])
        b = r.normal(size=(dim,))
        X = X - pot.gradV(X) * delta_t + np.sqrt(2 * delta_t/beta) * b
    #print (X, pot.V(X))

pot = StiffPot()
x_0 = np.array(pot.x0)
beta = 1.0
N = 1000000
delta_t = 0.01

t0 = time.time()
UnbiasedTraj(pot, x_0, beta, delta_t, N)
t1 = time.time()

print ('Runtime: %.2f sec' % (t1 - t0))

