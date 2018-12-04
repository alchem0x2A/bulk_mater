import os
import numpy
import json
import ase.db
from ase.constraints import StrainFilter
from ase.parallel import world, rank, parprint
from ase.optimize import BFGS
from ase.atoms import Atoms

from gpaw import GPAW, PW, FermiDirac

cur_dir = os.path.dirname(__file__)
default_json_file = os.path.abspath(os.path.join(cur_dir,
                                                 "../config/parameters.json"))

class MaterCalc(object):
    """Class wrapper for relaxation and calculation
    """
    def __init__(self,
                 atoms,          # A ase.atoms.Atoms object
                 base_dir,      
                 param_file=default_json_file):
        # Read the parameters
        if os.path.exists(param_file):
            with open(param_file, "r") as f:
                params = json.load(f)
            self.__params = params
        else:
            raise FileNotFoundError("No parameter file!")
        
        # Base_dir for all data
        if not os.path.exists(base_dir):
            if rank == 0:
                os.makedirs(base_dir)   # Recursively makedirs
        world.barrier()
        self.__base_dir = os.path.abspath(base_dir)

        if isinstance(atoms, Atoms):
            self.atoms = atoms
        else:
            raise TypeError("Atom must be an Atoms instance!")
        return
    
    @property
    def params(self):
        return self.__params

    @property
    def base_dir(self):
        return self.__base_dir

    def relax(self,
              method="IT",
              fmax=0.01,        # maximum force or stress * Volume
              steps=500):        # maximum steps
        atoms_copy = self.atoms.copy()  # makesure nothing happens
        method = method.upper()
        res_file = os.path.join(self.base_dir, "relaxed.traj")  # relaxed trajectory
        # Continue without calculation
        if os.path.exists(res_file):
            parprint("Relaxation Already Done!")
            return True
        
        calc = GPAW(**self.params["relax"],
                    txt=os.path.join(self.base_dir, "relax.txt"))
        atoms_copy.set_calculator(calc)
        traj_filename = os.path.join(self.base_dir,
                                     "relax-{}.traj".format(method))
        log_filename = os.path.join(self.base_dir,
                                     "relax-{}.log".format(method))

        # Now choose the method
        if method in ("UCF", "ECF"):  # UnitCellFilter
            filt = UnitCellFilter if method == "UCF" else ExpCellFilter
            ff = filt(atoms_copy)
            opt = BFGS(ff, trajectory=traj_filename, logfile=log_filename)
            opt.run(fmax=fmax, steps=steps)
            converged = opt.converged
        elif method.upper() == "IT":  # Iterative
            converged = False
            sf = StrainFilter(atoms_copy)
            loop = 0
            while (loop < 5) and (converged is False):
                opt_strain = BFGS(sf,
                                  trajectory=traj_filename,
                                  logfile=log_filename)
                opt_force = BFGS(atoms_copy,
                                 trajectory=traj_filename,
                                 logfile=log_filename)
                opt_strain.run(fmax=fmax, steps=steps)
                opt_force.run(fmax=fmax, steps=steps)
                parprint(opt_strain.nsteps, opt_force.nsteps)
                conv_step = (opt_strain.nsteps < 2) and (opt_force.nsteps < 2)
                conv_job = opt_strain.converged and opt_force.converged
                converged = conv_step and conv_job
                loop += 1
                
        if converged:
            self.atoms = atoms_copy  # copy back
            atoms_copy.write(res_file)  # parallel?
            parprint("Relaxation Done!")
            return True
        else:
            parprint("Something wrong with relaxation!")
            return False

    def ground_state(self):
        pass

    def bandgap(self):
        pass 

    def excited_state(self):
        pass

    def dielectric(self):
        pass


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
    
    
