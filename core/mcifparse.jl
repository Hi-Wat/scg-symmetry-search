module Mcifparse

include("./comp.jl")
using .Compute
using PyCall

function spindim_search(x::String)
    if x == "COLLINEAR"
        spindim = 1
    elseif x == "COPLANAR"
        spindim = 2
    elseif x == "NONCOPLANAR"
        spindim = 3
    else 
        spindim = 3
    end
    return spindim
end


##### from https://github.com/spglib/magndata-bench/blob/main/tests/conftest.py
py"""
from __future__ import annotations
import numpy as np
from pathlib import Path
import os,copy,itertools

from spglib import get_magnetic_symmetry
from spinspg import get_spin_symmetry
from pymatgen.core import Structure



def mgs2spglib(structure, magmom, symprec=1e-4, angle_tolerance=1e-1, mag_symprec=1e-3):
    lattice     = structure.lattice.matrix
    frac_coords = structure.frac_coords
    unique_species: list[Element|Species] = []
    zs = []
    for species, g in itertools.groupby(structure, key=lambda s: s.species):
        if species in unique_species:
            ind = unique_species.index(species)
            zs.extend([ind + 1] * len(tuple(g)))
        else:
            unique_species.append(species)
            zs.extend([len(unique_species)] * len(tuple(g)))
    cell = lattice, frac_coords, zs, magmom
    return cell

def spinspacegroup_identify(name: str, symprec=1e-3,angle_tolerance=1e-2, mag_symprec=1e-3):
    magnname = name
    struct = Structure.from_file(str(magnname))
    magmom = [list(site['properties']['magmom']) for site in struct.as_dict()['sites']]
    cell = mgs2spglib(struct, magmom, symprec=1e-4, angle_tolerance=angle_tolerance, mag_symprec=mag_symprec)
    cell = (cell[0], cell[1], cell[2], np.array(cell[3]))

    species,pgmat,trans,spmat = get_spin_symmetry(*cell,symprec=symprec)
    normalized_lat = (cell[0]/np.linalg.norm(cell[0],axis=1)).T
    return name,normalized_lat,str(species.spin_only_group_type),species.axis,list(pgmat),list(trans),list(spmat)


def shubnikovgroup_identify(name: str, symprec=1e-3):
    magnname = name
    struct = Structure.from_file(str(magnname))
    magmom = [list(site['properties']['magmom']) for site in struct.as_dict()['sites']]
    cell = mgs2spglib(struct, magmom, symprec=1e-4, angle_tolerance=1e-1, mag_symprec=1e-3)
    cell = (cell[0], cell[1], cell[2], np.array(cell[3]))

    magdata = get_magnetic_symmetry(cell, symprec=1e-3)
    normalized_lat = (cell[0]/np.linalg.norm(cell[0],axis=1)).T


    pgmat = np.array(magdata["rotations"])
    trans = np.array(magdata["translations"])
    spmat = copy.deepcopy(pgmat)

    improp_signs = list(map(lambda x : False if np.linalg.det(x)==1 else True, pgmat))
    spmat[improp_signs] *= -1
    spmat[magdata["time_reversals"]] *= -1

    return name,normalized_lat,"SHUBNIKOV",[0.0,0.0,0.0],list(pgmat),list(trans),list(spmat)


"""

function get_spinspacegroup(name::String,symprec=1.0e-4)
    magndat = py"spinspacegroup_identify"(name,symprec)
    spindim = spindim_search(magndat[3])
    if spindim == 3
        pure_spinaxis = [0.0,0.0,0.0]
    else
        pure_spinaxis = magndat[4]
    end

    opsets = Vector{SpinspaceOperation}([])
    for (pg,trans,sp) in zip(magndat[5],magndat[6],magndat[7])
        push!(opsets,SpinspaceOperation("",magndat[2]*pg*inv(magndat[2]),magndat[2]*trans,sp))
    end
    purespinaxis_args = Compute.purespin_axisanglesearch(spindim,pure_spinaxis)
    Compute.purespinops_add!(opsets,spindim,purespinaxis_args)
    
    return opsets, spindim, purespinaxis_args
end


function get_shubnikovgroup(name::String,symprec=1.0e-4)
    magndat = py"shubnikovgroup_identify"(name,symprec)
    spindim = spindim_search(magndat[3])
    opsets = Vector{Compute.SpinspaceOperation}([])
    for (pg,trans,sp) in zip(magndat[5],magndat[6],magndat[7])
        push!(opsets,Compute.SpinspaceOperation("",magndat[2]*pg*inv(magndat[2]),magndat[2]*trans,magndat[2]*sp*inv(magndat[2])))
    end

    return opsets, 3, [0.0,0.0]
end

function magnetization_sp(opsets::Vector{SpinspaceOperation},spindim::Int64,purespinaxis_args::Vector{Float64};message::Bool=false)
    
    physar = Array([
        phys_property("S",1,true,1,-1,false),
        ])
    select_mode = "eq"
    appind = 2
    ar,reshapear = tensor_search(
                    physar,opsets,spindim,purespinaxis_args
                    ,appind,mode=select_mode,interactive_bool=message,emrep=true
                    )
    display(reshapear)
end

function magnetization_orb(opsets::Vector{SpinspaceOperation},spindim::Int64,purespinaxis_args::Vector{Float64};message::Bool=false)
    
    physar = Array([
        phys_property("O",1,false,1,-1,false),
        ])
    select_mode = "eq"
    appind = 2
    ar,reshapear = tensor_search(
                    physar,opsets,spindim,purespinaxis_args
                    ,appind,mode=select_mode,interactive_bool=message,emrep=true
                    )
    display(reshapear)
end

function conductivity(opsets::Vector{SpinspaceOperation},spindim::Int64,purespinaxis_args::Vector{Float64};message::Bool=false)
    
    physar = Array([
        phys_property("J",1,false,-1,-1,false)
        ])
    select_mode = "linear_same"
    appind = 2
    ar,reshapear = tensor_search(
                    physar,opsets,spindim,purespinaxis_args
                    ,appind,mode=select_mode,interactive_bool=message,emrep=true
                    )
    display(reshapear)
end




function magnetoelectricity_sp(opsets::Vector{SpinspaceOperation},spindim::Int64,purespinaxis_args::Vector{Float64};message::Bool=false)
    
    physar = Array([
        phys_property("S",1,true,1,-1,false),
        phys_property("E",1,false,-1,1,false)
        ])
    select_mode = "linear"
    appind = 2
    ar,reshapear = tensor_search(
                    physar,opsets,spindim,purespinaxis_args
                    ,appind,mode=select_mode,interactive_bool=message,emrep=true
                    )
    reshapear[:,:,1] |>display
    reshapear[:,:,2] |>display
end

function magnetoelectricity_orb(opsets::Vector{SpinspaceOperation},spindim::Int64,purespinaxis_args::Vector{Float64};message::Bool=false)
    
    physar = Array([
        phys_property("O",1,false,1,-1,false),
        phys_property("E",1,false,-1,1,false)
        ])
    select_mode = "linear"
    appind = 2
    ar,reshapear = tensor_search(
                    physar,opsets,spindim,purespinaxis_args
                    ,appind,mode=select_mode,interactive_bool=message,emrep=true
                    )
    reshapear[:,:,1] |>display
    reshapear[:,:,2] |>display
end

function E_to_spincurrent(opsets::Vector{SpinspaceOperation},spindim::Int64,purespinaxis_args::Vector{Float64};message::Bool=false)
    
    physar = Array([
        phys_property("J",1,false,-1,-1,false),
        phys_property("s",1,true,1,-1,false),
        phys_property("E",1,false,-1,1,false)
        ])
    select_mode = "linear"
    appind = 3
    ar,reshapear = tensor_search(
                    physar,opsets,spindim,purespinaxis_args
                    ,appind,mode=select_mode,interactive_bool=message,emrep=true
                    )
    for i in 1:3
        display(reshapear[:,i,:,1])
        display(reshapear[:,i,:,2])
    end
end

function E_to_orbitalcurrent(opsets::Vector{SpinspaceOperation},spindim::Int64,purespinaxis_args::Vector{Float64};message::Bool=false)
    
    physar = Array([
        phys_property("J",1,false,-1,-1,false),
        phys_property("o",1,false,1,-1,false),
        phys_property("E",1,false,-1,1,false)
        ])
    select_mode = "linear"
    appind = 3
    ar,reshapear = tensor_search(
                    physar,opsets,spindim,purespinaxis_args
                    ,appind,mode=select_mode,interactive_bool=message,emrep=true
                    )
    for i in 1:3
        display(reshapear[:,i,:,1])
        display(reshapear[:,i,:,2])
    end
end



end ### end module 


