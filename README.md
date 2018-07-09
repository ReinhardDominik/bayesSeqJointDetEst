Bayesian Sequential Joint Detection and Estimation
==================================================

This repository contains the code to reproduce the results of the paper "Bayesian Sequential Joint Detection and Estimation" (submitted to Sequential Analysis). If you want to use this code please cite
```
@article{reinhard2018,
author =    {Reinhard, Dominik and Fauss{}, Michael and Zoubir, Abdelhak M.},
  title     = {{Bayesian Sequential Joint Detection and Estimation}},
  journal   = {{Sequential Analysis (submitted)}},
  year      = {2018},
}
```

Authors
----------
Dominik Reinhard, Michael Fauss and Abdelhak M. Zoubir are with the Signal Processing Group at TU Darmstadt.

### Contact ####
Dominik Reinhard

reinhard@spg.tu-darmstadt.de

Requirements
------------
* MATLAB (at least 2016b)
* MATLAB Signal Processing Toolbox
* CVX
* Java
* Gurobi (may be replaced by another solver)

Check dependencies by calling
```matlab
>> checkDependencies();
```

Tested with
-----------
* MATLAB version 2016b
* MATLAB Signal Processing Toolbox version 7.3
* CVX version 2.1
* Gurobi version 7.02
* Java version 1.7.0_60

Usage
-----
Check all dependencies by running
```matlab
checkDependencies()
```
Please make sure, that you have a working cvx installation and the path variable is set correctly.
If you do not have a gurobi license have a look at getPrefStruct.m and change the cvx_solver accordingly.

The folder examples containts the scripts to reproduce the results presented in the paper.

If you have any troubles running the code, please contact the author.

