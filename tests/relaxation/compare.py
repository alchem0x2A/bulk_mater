from __future__ import print_function
import sys
import os, os.path
import shutil
#Importing issue with gpaw-python?
#However cannot use this for module importing
sys.path.append(os.path.join(os.path.dirname(__file__), "../../"))
import pkgutil
import unittest
from bulk.build import StructureBuilder
from gpaw import GPAW, PW, FermiDirac
import ase.optimize
from ase.constraints import StrainFilter, UnitCellFilter, ExpCellFilter
from ase.io.trajectory import Trajectory
from ase.parallel import rank, size, parprint, world

# class TestPath(unittest.TestCase):
    # def test_module(self):
        # def sub_module(package):
            # submod = []
            # for importer, modname, ispkg in pkgutil.iter_modules(package.__path__):
                # print("Found submodule {} (is a package: {})".format(modname,
                                                                     # ispkg))
                # submod.append(modname)
            # return submod
        # submod = sub_module(bulk)
        # self.assertGreater(len(submod), 0)

def main():
    candidates = [
                  ("Si", "diamond"),
                  ("Ge", "diamond"),
                  ("GaAs", "zincblende"),
                  ("ZrN", "rocksalt"),
                  ("BN", "wurtzite"),
                  ("GaN", "wurtzite"),
                  ("TlCl", "cesiumchloride"),
    ]
    for name, proto in candidates:
        run_single(name=name,
                   prototype=proto)
    return


def run_single(name="Si", prototype=None):
    sb = StructureBuilder()
    mater = sb.get_structure(formula=name,
                             prototype=prototype)
    parprint("Before", mater)
    if len(mater) != 1:
        return
    mater = mater[0]
    base_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                            "../../tmp/", "{}-{}/".format(name, prototype))
    if rank == 0:
        if not os.path.exists(base_dir):
            os.makedirs(base_dir)
    world.barrier()
    calc = GPAW(mode=dict(name="pw",
                          ecut=800),
                occupations=dict(name="fermi-dirac",
                                 width=0.01),
                basis="dzp",
                xc="PBE",
                kpts=dict(gamma=True,
                          density=4.0),  # Very rough k-density
                txt=os.path.join(base_dir,
                                 "{}-{}.txt".format(name, prototype)))
    mater.set_calculator(calc)
    sf = StrainFilter(mater)
    traj_filename = os.path.join(base_dir,
                                 "{}-{}.traj".format(name,
                                                     prototype))
    log_filename = os.path.join(base_dir,
                                "{}-{}.log".format(name,
                                                   prototype))
    # Overwrite file
    with open(traj_filename, "w") as _:
        pass
    with open(log_filename, "w") as _:
        pass
    
    opt_strain = ase.optimize.BFGS(sf,
                                   trajectory=traj_filename,
                                   logfile=log_filename)
    opt_force = ase.optimize.BFGS(mater,
                                  trajectory=traj_filename,
                                  logfile=log_filename)
    opt_strain.run(fmax=0.02)
    # opt_force.run(fmax=0.02)
    
    parprint("After", mater)

if __name__ == "__main__":
    main()
    # unittest.main()
