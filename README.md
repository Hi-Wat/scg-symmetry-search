# scg-symmetry-search

Search the symmetry-adapted tensors based on the spin space/Shubnikov groups for a given magnetic material.

## Preparation

The program consisting of Julia(`.jl`) files depends on `spglib` and `spinspg` packages of Python.

1. prepare Julia and Python3 environment in your machine.

`scg-symmetry-search` also depends on the `SymPy` of Python.

```
> pip install sympy
> pip install spinspg # which automatically install spglib
```

2. 


3. switch default interpreter for Python called in Julia 

```
> julia
julia > using PyCall
julia > PyCall.python # check which python interpreter is called
julia> ENV["PYTHON"] = raw"(path for the Python executer installing spglib, spinspg )"
julia > ]
Pkg > build PyCall
julia > PyCall.python # check the python interpreter changed
```
Note that PyCall.jl currently does not call the python executer prepared in virtual environment through `conda create`.


## Introduction 

Sample cases for searching symmetry-adapted tensors are prepared by `sample.ipynb`. 





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


