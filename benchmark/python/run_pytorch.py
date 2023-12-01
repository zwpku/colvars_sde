#!/usr/bin/env python
import time
from tqdm import tqdm
import random
import numpy as np
import torch

def set_all_seeds(seed):
    torch.manual_seed(seed)
    random.seed(seed)

class StiffPot(torch.nn.Module):
    def __init__(self):
        super().__init__()
        
    def forward(self, X):
        return (X[0]**2 - 1)**2 + 2.0 * (X[0]**2 + X[1] - 1)**2

def UnbiasedTraj(pot, X_0, beta, delta_t=1e-3, N=1000, seed=0):

    r = np.random.RandomState(seed)
    X = X_0
    dim = X.shape[0]
    for i in range(N):
        #b = np.array([0.1, 0.1])
        b = r.normal(size=(dim,))
        input_tensor = torch.tensor(X, requires_grad=True)
        V = pot(input_tensor)
        grad = torch.autograd.grad(V, input_tensor)[0].numpy()
        X = X - grad * delta_t + np.sqrt(2 * delta_t/beta) * b
    #print (X, pot.V(X))

pot = StiffPot()
x_0 = np.array([-1.0,0])
beta = 1.0
N = 1000000
delta_t = 0.01

t0 = time.time()
UnbiasedTraj(pot, x_0, beta, delta_t, N)
t1 = time.time()

print ('Runtime: %.2f sec' % (t1 - t0))
