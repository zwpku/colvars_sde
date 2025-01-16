### Example 1:

Sampling a two-dimensional system with a stiff potential using adaptive biasing force (ABF) method. The CV is defined as the $x$ coordinate, i.e. $\xi(x,y)=x$.

#### Files:

- `params.cfg` contains general simulation parameters. 
- `colvars-abf-coordinate.in` contains the parameters in ABF and the definition of CV. The format is similar to the config file in colvars package. 

The config file is specified by the following line in params.cfg:

> colvars colvars-abf-coordinate.in

#### Simulation: 

```
colvars_sde params.cfg
```

#### Examine the result:

A Jupyter notebook is provided, which contains code to view the results, such as the trajectory and the calculated free energy. 


```
jupyter notebook analysis.ipynb

```

### Example 2:

Sampling the same system. The CV is defined by a TorchScript model, stored in the file `scripted_cv_cpu.pt`. 

In this case, edit the parameter file `params.cfg` to specify the config file as follows:

> #colvars colvars-abf-coordinate.in

>  colvars colvars-abf.in

The steps for simulation and examination are the same as above.

