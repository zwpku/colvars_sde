import torch

class MyModel(torch.nn.Module):
    def __init__(self):
        super().__init__()
    def forward(self, X):
        return (X[0]**2 - 1)**2 + 2.0 * (X[0]**2 + X[1] - 1)**2

model = MyModel()
scripted_cv_filename = f'./stiff_pot.pt'
torch.jit.script(model).save(scripted_cv_filename)

