import os
import numpy
import json
import ase.db
from ase.constraints import StrainFilter
from ase.parallel import world, rank, parprint
from ase.optimize import BFGS
from ase.atoms import Atoms
from ase.io import write, read
from ase.dft.bandgap import bandgap as bg
from ase.dft.kpoints import special_paths, get_cellinfo
from ase.units import Ha

from gpaw import GPAW, PW, FermiDirac
from gpaw.response.df import DielectricFunction

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
        # Handle file path-related issues
        self.__base_dir = os.path.abspath(base_dir)
        self.__relaxed_traj = os.path.join(self.__base_dir,
                                           "relaxed.traj")
        self.__gs_file = os.path.join(self.__base_dir,
                                     "gs.gpw")  # ground state in PBE
        self.__gs_gllb_file = os.path.join(self.__base_dir,
                                           "gs_gllb.gpw")  # ground state in PBE
        self.__bg_file_template = os.path.join(self.__base_dir,
                                               "bg_{}.gpw")  # to add later
        self.__es_file = os.path.join(self.__base_dir,
                                      "es.gpw")  # excited states in PBE
        self.__eps_file_template = os.path.join(self.__base_dir,
                                                "eps_{}.npz")

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
        # Continue without calculation
        if os.path.exists(self.__relaxed_traj):
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
            atoms_copy.write(self.__relaxed_traj)  # parallel?
            parprint("Relaxation Done!")
            return True
        else:
            parprint("Something wrong with relaxation!")
            return False

    def ground_state(self):
        """Calculate ground state and generate excited wavefunction
        """
        if not os.path.exists(self.__relaxed_traj):
            raise FileNotFoundError("Need to relax first!")
        
        if os.path.exists(self.__gs_file):
            parprint("Ground state already done!")
            return True
        # relaxed.traj must present
        traj = read(self.__relaxed_traj)
        calc = GPAW(**self.params["gs"],
                    txt=os.path.join(self.base_dir,
                                     "gs.txt"))
        traj.set_calculator(calc)
        try:
            traj.get_potential_energy()
            calc.write(self.__gs_file)  # write
        except Exception as e:
            parprint("Something wrong with gs calculation!")
            parprint("ErrMsg: {}".format(e))
            return False
        else:
            return True
            

    def bandgap(self, method="GLLB", save=True):
        # Check method input
        method = method.strip().lower()
        if method not in ("pbe", "gllb",
                          "hse", "gw"):
            raise ValueError("Method for bandgap calculation not known!")
        # Check ground state
        if not os.path.exists(self.__gs_file):
            raise FileNotFoundError("Ground state not calculated!")
        
        # Check bandgap?
        self.__bg_file = self.__bg_file_template.format(method)
        if os.path.exists(self.__bg_file):
            parprint("Bandgap for method {} is already calculated!".format(method))
            return True
        
        # Real calculation now
        lattice_type = get_cellinfo(self.atoms.cell).lattice
        use_bs = lattice_type in special_paths.keys()
        if method == "pbe":
            # Use band_structure or simple sampling kpts?
            if use_bs:
                kpts = dict(path=special_paths[lattice_type],
                            npoints=self.params["gap"][method]["npoints"])
            else:
                kpts = dict(density=12)  # denser
            calc = GPAW(restart=self.__gs_file)
            calc.set(kpts=kpts,
                     fixdensity=True,
                     symmetry="off",
                     txt=None)  # non-SC
            calc.get_potential_energy()  # Otherwise the wavefunction not known?
            bg_min, *_ = bg(calc, direct=False)
            bg_dir, *_ = bg(calc, direct=True)
            if use_bs:
                bs = calc.band_structure()
                xcoords, label_xcoords, labels = bs.get_labels()
                res_bs = bs.todict()
                res_bs.update(xcoords=xcoords,
                              label_xcoords=label_xcoords,
                              labels=labels)
            else:
                res_bs = None
            return bg_min, bg_dir, res_bs
        elif method == "gllb":
            # Need to restart the SC calculation
            if not os.path.exists(self.__gs_gllb_file):
                calc_ = GPAW(restart=self.__gs_file, txt=None)
                calc = GPAW(**self.params["gs"])
                calc.atoms = calc_.atoms.copy()
                calc.set(xc="GLLBSC")
                del calc_
                calc.get_potential_energy()  # Recalculate GS
                calc.write(self.__gs_gllb_file)
            #TODO: merge with PBE method
            if use_bs:
                kpts = dict(path=special_paths[lattice_type],
                            npoints=self.params["gap"][method]["npoints"])
            else:
                kpts = dict(density=12)
            calc_bs = GPAW(restart=self.__gs_gllb_file)
            calc_bs.set(kpts=kpts,
                        fixdensity=True,
                        symmetry="off")
            calc_bs.get_potential_energy()
            homolumo = calc_bs.get_homo_lumo()
            response = calc_bs.hamiltonian.xc.xcs["RESPONSE"]
            response.calculate_delta_xc( homolumo / Ha)
            EKs, Dxc = response.calculate_delta_xc_perturbation()
            ns = calc_bs.get_number_of_spins()
            nk = len(calc_bs.get_ibz_k_points())
            e_kn = numpy.array([[calc_bs.get_eigenvalues(kpt=k, spin=s) \
                                for k in range(nk)] for s in range(ns)])
            efermi = calc_bs.get_fermi_level()
            e_kn[e_kn > efermi] += Dxc
            bg_gllb_min, *_ = bg(eigenvalues=e_kn,
                                      efermi=efermi,
                                      direct=False)
            bg_gllb_dir, *_  = bg(eigenvalues=e_kn,
                                       efermi=efermi,
                                       direct=True)
            if use_bs:
                bs = calc_bs.band_structure()
                bs.energies[bs.energies > bs.reference] += Dxc  # Shift gllbsc discontinuous
                res_bs = bs.todict()
                xcoords, label_xcoords, labels = bs.get_labels()
                res_bs.update(xcoords=xcoords,
                              label_xcoords=label_xcoords,
                              labels=labels)
            else:
                res_bs = None
            if save:
                if rank == 0:
                    numpy.savez(self.__bg_file,
                                Eg_min=bg_gllb_min,
                                Eg_dir=bg_gllb_dir,
                                **res_bs)
                    print("Bandgap saved!")
            world.barrier()
            return bg_gllb_min, bg_gllb_dir, res_bs
        else:
            raise NotImplementedError("Method not implemented yet")

    def excited_state(self):
        # Only depend on the ground state
        if not os.path.exists(self.__gs_file):
            raise FileNotFoundError("Ground state not calculated!")
        
        if os.path.exists(self.__es_file):
            parprint("Excited state calculated!")
            return True

        calc = GPAW(restart=self.__gs_file,
                    **self.params["es"])
        calc.get_potential_energy()
        nv = calc.get_number_of_electrons() // 2
        nbands = max(70, nv * 3)  # include empty bands
        # calc.set(**self.params["es"])  # only parallel over 1 kpt
        calc.diagonalize_full_hamiltonian(nbands=nbands)  # diagonalize
        calc.write(self.__es_file, mode="all")
        parprint("Excited states calculated!")
        return True
        
    def dielectric(self, method="rpa"):
        method = method.lower()
        if method not in ("rpa", "gw"):
            raise ValueError("Dielectric Method not known!")
        if not os.path.exists(self.__es_file):
            raise FileNotFoundError("Ground state not calculated!")
        
        self.__eps_file = self.__eps_file_template.format(method)
        if os.path.exists(self.__eps_file):
            parprint(("Dielectricfunction using"
                      " method {} already calculated!").format(method))
            return True
        if os.path.exists(self.__es_file):
            parprint("Excited state done, will use directly!")
        if method == "rpa":
            df = DielectricFunction(calc=self.__es_file,
                                    **self.params[method])
            epsx0, epsx = df.get_dielectric_function(direction="x",
                                                     filename=None)
            epsy0, epsy = df.get_dielectric_function(direction="y",
                                                     filename=None)
            epsz0, epsz = df.get_dielectric_function(direction="z",
                                                     filename=None)
            freq = df.get_frequencies()
            data = dict(frequencies=freq,
                        eps_x=epsx, eps_x0=epsx0,
                        eps_y=epsy, eps_y0=epsy0,
                        eps_z=epsz, eps_z0=epsz)
            # write result
            if rank == 0:
                numpy.savez(self.__eps_file, **data)
            parprint("Dielectric function using {} calculated!".format(method))
            return True
        else:
            raise NotImplementedError("{} not implemented".format(method))


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
    
    
