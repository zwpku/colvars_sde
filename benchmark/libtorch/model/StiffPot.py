import torch

class MyModel(torch.nn.Module):
    def __init__(self, eps):
        super().__init__()
        self.eps = eps
    def forward(self, X):
        return ((X[:,0]**2 - 1)**2 + 1.0 / self.eps * (X[:,0]**2 + X[:,1] - 1)**2).reshape((-1,1))

model = MyModel(0.5)
scripted_cv_filename = f'./stiff_pot.pt'
torch.jit.script(model).save(scripted_cv_filename)

