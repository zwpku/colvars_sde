## Enhanced sampling of general stochastic differential equations (SDEs) 

The implementation of enhanced sampling algorithms in this package is largely based on modification of the code in [colvars](https://github.com/Colvars/colvars).

#### Build 

1. (Optional) Download [libtorch C++](https://pytorch.org/cppdocs/installing.html) library, in order to allow the use of collective variables (CVs) defined by a TorchScript model. 

2. Configure

```
   mkdir build 
   cd build
   cmake ../cmake/CMakeLists.txt
   
```

   To enable libtorch support, edit the generated `CMakeCache.txt` file, and set 1) `COLVARS_SDE_TORCH` to `ON` and 2) `Torch_DIR` to the path to libtorch (e.g. path/to/libtorch/share/cmake/Torch). Otherwise, set `COLVARS_SDE_TORCH` to `OFF`. Run 

```
   cmake ../cmake/CMakeLists.txt
```

once again. 

3. Compile

```
   make 
   make install
```   

