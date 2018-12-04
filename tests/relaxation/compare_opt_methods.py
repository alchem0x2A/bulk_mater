from __future__ import print_function
import sys
import os, os.path
import shutil
#Importing issue with gpaw-python?
#However cannot use this for module importing
sys.path.append(os.path.join(os.path.dirname(__file__), "../../"))
import unittest
from bulk.build import StructureBuilder, convert_name
from gpaw import GPAW, PW, FermiDirac
import ase.optimize
from ase.constraints import StrainFilter, UnitCellFilter, ExpCellFilter
from ase.io.trajectory import Trajectory
from ase.parallel import rank, size, parprint, world
import json
import numpy


cur_dir = os.path.abspath(os.path.dirname(__file__))

def main(begin=0, end=None):
    sb = StructureBuilder()
    assert begin >= 0
    entries = sb.entries
    if end is None:
        end = len(entries)
    candidates = entries[begin : end]  # End is not calculated!
    results_all = []
    for i, line in enumerate(candidates):
        name = line["formula"]
        prototype = line["prototype"]
        res_single = run_single(sb, name, prototype)
        for entry in res_single:
            results_all.append(entry)

    parprint(type(results_all), results_all)
    world.barrier()
    if rank == 0:
        with open(os.path.join(cur_dir,
                               "../../tmp/",
                               "relax_result_{}-{}.json".format(begin, end)), "w") as f:
            json.dump(results_all, f)
    return

# 
def run_single(sb, name="Si", prototype=None,
               root_dir=os.path.join(cur_dir, "../../tmp/")):
    """
    root_dir --> base_dir --> files
    """
    results = []
    mater = sb.get_structure(formula=name,
                             prototype=prototype)
    parprint("Before", mater)
    if len(mater) != 1:
        return
    mater = mater[0]
    base_dir = os.path.join(root_dir, "{}-{}/".format(name, prototype))
    if rank == 0:
        if not os.path.exists(base_dir):
            os.makedirs(base_dir)
    world.barrier()
    def relax_use_method(xc, method,
                         fmax=0.002, steps=100,
                         clear=False):  # allow only up to 100 steps
        atoms = mater.copy()                      # same start point!
        if xc.upper() not in  ("PBE", "LDA"):
            raise ValueError("XC method not known!")
        if method.upper() not in ("UCF", "ECF", "IT"):
            raise ValueError("Optimization method not known!")
        parprint("Running {}-{} for {}-{}".format(xc, method, name, prototype))
        
        calc = GPAW(mode=dict(name="pw",
                              ecut=800),
                    occupations=dict(name="fermi-dirac",
                                     width=0.01),
                    basis="dzp",
                    xc=xc.upper(),
                    kpts=dict(gamma=True,
                              density=4.0),  # Very rough k-density
                    txt=os.path.join(base_dir,
                                     "{}-{}.txt".format(xc.upper(),
                                                        method.upper())))
        atoms.set_calculator(calc)
        traj_filename = os.path.join(base_dir,
                                     "{}-{}.traj".format(xc.upper(),
                                                         method.upper()))
        
        log_filename = os.path.join(base_dir,
                                    "{}-{}.log".format(xc.upper(),
                                                       method.upper()))

        if clear is True:
            with open(traj_filename, "w") as _:
                pass
            with open(log_filename, "w") as _:
                pass

            res_fmax = None
            if method.upper() == "UCF":  # UnitCellFilter
                ff = UnitCellFilter(atoms)
                opt = ase.optimize.BFGS(ff, trajectory=traj_filename, logfile=log_filename)
                opt.run(fmax=fmax, steps=steps)
                f = ff.get_forces()
                smax = numpy.sqrt((f ** 2).sum(axis=1).max())
            elif method.upper() == "ECF":  # ExpCellFilter
                ff = ExpCellFilter(atoms)
                opt = ase.optimize.BFGS(ff, trajectory=traj_filename, logfile=log_filename)
                opt.run(fmax=fmax, steps=steps)
                f = ff.get_forces()
                smax = numpy.sqrt((f ** 2).sum(axis=1).max())
            elif method.upper() == "IT":  # Iterative 
                sf = StrainFilter(atoms)
                opt_strain = ase.optimize.BFGS(sf,
                                               trajectory=traj_filename,
                                               logfile=log_filename)
                opt_force = ase.optimize.BFGS(atoms,
                                              trajectory=traj_filename,
                                              logfile=log_filename)
                opt_strain.run(fmax=fmax, steps=steps)
                opt_force.run(fmax=fmax, steps=steps)
                f = atoms.get_forces()
                smax = numpy.sqrt((f ** 2).sum(axis=1).max())
            
        a, b, c, *_ = atoms.get_cell_lengths_and_angles()
        if "smax" not in dir():
            smax=None
        res_entry = dict(name=name,
                         prototype=prototype,
                         xc=xc.upper(),
                         method=method.upper(),
                         abc=(a, b, c),
                         max_force=smax)
        return res_entry


    for xc in ("PBE", "LDA"):
        for method in ("UCF", "ECF", "IT"):
    # for xc in ("PBE", ):
        # for method in ("UCF", "ECF", "IT"):
            res = relax_use_method(xc, method, clear=True)
            results.append(res)
    return results

if __name__ == "__main__":
    if len(sys.argv) != 3:
        parprint("Nedd begin and end!")
    else:
        begin = int(sys.argv[1])
        end = int(sys.argv[2])
        if end < 0:
            end = None
        main(begin, end)
    # unittest.main()
