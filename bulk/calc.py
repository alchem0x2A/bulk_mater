import os
import numpy
import json
import ase.db
from ase.constraints import StrainFilter, UnitCellFilter, ExpCellFilter
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
        self.__relaxed_gllb_traj = os.path.join(self.__base_dir,
                                                "relaxed_gllb.traj")
        self.__gs_file = os.path.join(self.__base_dir,
                                     "gs.gpw")  # ground state in PBE
        self.__gs_gllb_file = os.path.join(self.__base_dir,
                                           "gs_gllb.gpw")  # ground state in PBE
        self.__bg_file_template = os.path.join(self.__base_dir,
                                               "bg_{}.npz")  # to add later
        self.__es_file = os.path.join(self.__base_dir,
                                      "es.gpw")  # excited states in PBE
        self.__eps_file_template = os.path.join(self.__base_dir,
                                                "eps_{}.npz")

        if isinstance(atoms, Atoms):
            self.atoms = atoms
        elif atoms is None:     # Dummy instance for checking only
            self.atoms = None
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
              steps=500,
              xc="pbe",
              skip=False):        # maximum steps
        atoms_copy = self.atoms.copy()  # makesure nothing happens
        method = method.upper()
        # Continue without calculation
        xc = xc.lower()
        if xc == "pbe":
            relaxed_traj = self.__relaxed_traj
        elif xc == "gllb":
            relaxed_traj = self.__relaxed_gllb_traj
        else:
            raise NotImplementedError("xc method not implemented")
            
        if skip:                # Do not relax
            atoms_copy.write(relaxed_traj)
            return True
        
        if os.path.exists(relaxed_traj):
            parprint("Relaxation Already Done!")
            return True
        
        calc = GPAW(**self.params["relax"],
                    txt=os.path.join(self.base_dir, "relax_{}.txt".format(xc)))
        if xc == "gllb":
            calc.set(xc="PBEsol")
        atoms_copy.set_calculator(calc)
        if xc == "pbe":
            traj_filename = os.path.join(self.base_dir,
                                         "relax-{}-{}.traj".format(method, xc))
            log_filename = os.path.join(self.base_dir,
                                        "relax-{}-{}.log".format(method, xc))
        else:
            traj_filename = None
            log_filename = None

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
            max_loop = 10
            while (loop < max_loop) and (converged is False):
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
                
        if not converged:
            # But still use the "relaxed trajectory"
            parprint("Relaxation result might be unstable!")
        else:
            parprint("Relaxation Done!")
        if xc == "pbe":
            self.atoms = atoms_copy  # copy back only for pbe
        atoms_copy.write(relaxed_traj)  # parallel?
        return True

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
            calc.write(self.__gs_file, mode="all")  # write
        except Exception as e:
            parprint("Something wrong with gs calculation!")
            parprint("ErrMsg: {}".format(e))
            return False
        else:
            return True
            

    def bandgap(self, method="GLLB", save=True, skip=False):
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
            data = numpy.load(self.__bg_file)
            return data["Eg_min"], data["Eg_dir"], None
        
        # Real calculation now
            
        if method == "pbe":
            # Use band_structure or simple sampling kpts?
            calc = GPAW(restart=self.__gs_file)
            try:
                lattice_type = get_cellinfo(calc.atoms.cell).lattice
                use_bs = lattice_type in special_paths.keys()
            except (ValueError, AssertionError):      # Monolithic cell?
                use_bs = False
            parprint("Use bandpath?", use_bs)
            if use_bs:
                kpts = dict(path=special_paths[lattice_type],
                            npoints=self.params["gap"][method]["npoints"])
                calc.set(kpts=kpts,
                         fixdensity=True,
                         symmetry="off",
                         txt=None)  # non-SC
            calc.get_potential_energy()  # Otherwise the wavefunction not known?
            bg_min, *_ = bg(calc, direct=False)
            bg_dir, *_ = bg(calc, direct=True)
            # Possibly no gap
            if bg_min is None:
                bg_min = 0
            if bg_dir is None:
                bg_dir = 0
            if use_bs:
                bs = calc.band_structure()
                xcoords, label_xcoords, labels = bs.get_labels()
                res_bs = bs.todict()
                res_bs.update(xcoords=xcoords,
                              label_xcoords=label_xcoords,
                              labels=labels)
            else:
                res_bs = {}     # no result
            if save:
                if rank == 0:
                    numpy.savez(self.__bg_file,
                                Eg_min=bg_min,
                                Eg_dir=bg_dir,
                                **res_bs)
                    print("Bandgap saved!")
            world.barrier()
            return bg_min, bg_dir, res_bs
        elif method == "gllb":
            # Need to restart the SC calculation
            if not os.path.exists(self.__gs_gllb_file):
                if not os.path.exists(self.__relaxed_gllb_traj):
                    self.relax(fmax=0.002, steps=200, xc="gllb", skip=skip)
                world.barrier()
                atoms = read(self.__relaxed_gllb_traj)
                calc = GPAW(**self.params["gs"], txt=os.path.join(self.__base_dir,
                                                                  "gs_gllb.txt"))
                atoms.set_calculator(calc)
                calc.set(xc="GLLBSC")
                atoms.get_potential_energy()  # Recalculate GS
                calc.write(self.__gs_gllb_file, mode="all")
            #TODO: merge with PBE method
            calc_bs = GPAW(restart=self.__gs_gllb_file)
            try:
                lattice_type = get_cellinfo(calc_bs.atoms.cell).lattice
                use_bs = lattice_type in special_paths.keys()
            except (ValueError, AssertionError):      # Monolithic cell?
                use_bs = False
            if use_bs:
                kpts = dict(path=special_paths[lattice_type],
                            npoints=self.params["gap"][method]["npoints"])
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
            if bg_gllb_min is None:
                bg_gllb_min = 0
            if bg_gllb_dir is None:
                bg_gllb_dir = 0
            if use_bs:
                bs = calc_bs.band_structure()
                bs.energies[bs.energies > bs.reference] += Dxc  # Shift gllbsc discontinuous
                res_bs = bs.todict()
                xcoords, label_xcoords, labels = bs.get_labels()
                res_bs.update(xcoords=xcoords,
                              label_xcoords=label_xcoords,
                              labels=labels)
            else:
                res_bs = {}
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

        calc = GPAW(restart=self.__gs_file, parallel=dict(kpt=1))
        nv = calc.get_number_of_electrons() // 2
        nbands = max(70, nv * 3)  # include empty bands
        parprint(self.params["es"])
        calc.set(**self.params["es"])
        calc.set(nbands=int(nbands))
        calc.get_potential_energy()
        
        # calc.set(**self.params["es"])  # only parallel over 1 kpt
        calc.diagonalize_full_hamiltonian(nbands=int(nbands))  # diagonalize
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
    
    def check_status(self,
                     bg_method="gllb",
                     eps_method="rpa"):
        """Check the status of calculation.
        """
        exists = os.path.exists
        result = dict(relax=exists(self.__relaxed_traj),
                      gs=exists(self.__gs_file),
                      bandgap=exists(self.__bg_file_template.format(bg_method.lower())),
                      es=exists(self.__es_file),
                      dielectric=exists(self.__eps_file_template.format(eps_method.lower())))
        return result
