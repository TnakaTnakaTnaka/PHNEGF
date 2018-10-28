PHNEGF
====

### Version 0.9.1

# Overview
The PHNEGF are the scripts interfaced with ALAMODE software to investigate the ballistic phonon transport in the bulk system on the basis of Nonequilibrium Green’s function (NEGF) method. ALAMODE is a software package designed for analyzing lattice anharmonicity and lattice thermal conductivity of solids by using an external DFT package. By using PHNEGF, you can calculate the phonon transmittance and thermal conductance from the results of harmonic interatomic force constants (IFCs) in ALAMODE.


* mod_dymat.py is a python module for generating dynamical matrix by using hessian file in ALAMODE. 

* NEGF(-mulp).py is a script to calculate q-resolved phonon transmittance by using NEGF method. This script needs to be combined with mod_dymat.py.

* tran.py is a script to calculate transmittance averaged in 1st Brillouin zone from the data obtained by NEGF.py.

* kappa.py is a script to calculate phonon thermal ocnductance from transmittance data.

# Download
You can download the latest versions at https://github.com/TnakaTnakaTnaka/PHNEGF. If you download the github version, please use the 'master' branch.

```
$ git clone http://github.com/TnakaTnakaTnaka/PHNEGF.git
```

# Documentation
For more details about PHNEGF, see the manual (manual.pdf).

# Reference
ALAMODE software package

http://alamode.readthedocs.io/

# Author
Yuto Tanaka 
