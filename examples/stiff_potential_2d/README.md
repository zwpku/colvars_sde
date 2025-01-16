### Example:


Sampling a two-dimensional stiff potential using adaptive biasing force (ABF) method. The CV is defined as the $x$ coordinate, i.e. $\xi(x,y)=x$.

#### Files:

- `params.cfg` contains general simulation parameters. 
- `colvars-abf-coordinate.in` contains the parameters in ABF and the definition of CV. The format is similar to the config file in colvars package. 

#### Simulation:

```
colvars_sde params.cfg
```

#### Examine the result:

A Jupyter notebook is provided, which contains code to view the results, such as the trajectory and the calculated free energy. 


```
jupyter notebook analysis.ipynb

```

