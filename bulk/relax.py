import os, os.path
import numpy
import json
import ase.optimize
import ase.db
from ase.io import read
from ase.constraints import StrainFilter
from ase.parallel import world, rank, parprint
from ase.optimize import BFGS

from gpaw import GPAW, PW, FermiDirac

# Relax single atom
def relax(atoms, name="",
          base_dir="./",
          smax=2e-4):
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    param_file = os.path.join(curr_dir, "../parameters.json")
    gpw_file = os.path.join(base_dir, "gs.gpw")
    if os.path.exists(gpw_file):
        parprint("Relaxation already done, will use gpw directly!")
        return 0
    if os.path.exists(param_file):
        params = json.load(open(param_file, "r"))
    else:
        raise FileNotFoundError("no parameter file!")
    
    # calculation asign
    calc = GPAW(**params["relax"])
    atoms.set_calculator(calc)
    traj_filename = os.path.join(base_dir,
                                 "{}_relax.traj".format(name))
    log_filename = os.path.join(base_dir,
                                 "{}_relax.log".format(name))
    opt = QuasiNewton(atoms,
                      trajectory=traj_filename,
                      logfile=log_filename)
    mask = [1, 1, 1, 0, 0, 0]   # Relax for bulk
    opt.run(fmax=0.01, smax=smax, smask=mask)
    
    # Calculate the ground state 
    calc.set(**params["gs"])
    atoms.get_potential_energy()
    calc.write(gpw_file)
    
    
def optimize_method(atom,
                    method="IT",
                    fmax=0.005,  # max force or stess * V
                    steps=400,):      # max steps
    method = method.upper()
    
    
