# scg-symmetry-search

Search the symmetry-adapted tensors based on the spin space/Shubnikov groups for a given magnetic material.

## Preparation

The programs consisting of Julia(`.jl`) files are on top of `spglib` and `spinspg` packages developed in Python language.

1. prepare Julia and Python3 environment in your machine.

Note that `spinspg` requires the Python interpreter, `Python>=3.8.0`.  
REF : https://pypi.org/project/spinspg/


Install the required packages for Python.
```
# recommend to use pip in installation of packages
# > sudo apt install python3-pip
> python3 -m pip install sympy
> python3 -m pip install pymatgen
> python3 -m pip install spinspg
```

Julia installation should be `>=1.8.5`  
To compute the program, install the follwing Julia packages
```
> julia 
julia> ]
Pkg> add LinearAlgebra
Pkg> add PyCall
Pkg> add Optim
Pkg> add SymPy
```

2. switch default interpreter for Python called in Julia 

```
> julia
julia > using PyCall
julia > PyCall.python # check which python interpreter is called
julia> ENV["PYTHON"] = raw"(path for the Python executer installing spglib, spinspg )"
julia > ]
Pkg > build PyCall
julia > PyCall.python # check the python interpreter changed
```
Note that PyCall.jl currently does not call the python interpreter prepared in virtual environment through `conda create` (as of 2023/07).


## Introduction 

Samples for searching symmetry-adapted tensors are prepared in `sample.ipynb`. 





## Citation

If you use scg-symmetry-search in your research, please cite `spglib`, `spinspg`, and the following bibliography

```
@misc{Watanabe2023_scg,
      title={Symmetry Analysis with Spin Crystallographic Groups: Disentangling Spin-Orbit-Free Effects in Emergent Electromagnetism}, 
      author={Hikaru Watanabe and Kohei Shinohara and Takuya Nomoto and Atsushi Togo and Ryotaro Arita},
      year={2023},
      eprint={2307.11560},
      archivePrefix={arXiv},
      primaryClass={cond-mat.mtrl-sci}
}
```

## Licence 

The program 

