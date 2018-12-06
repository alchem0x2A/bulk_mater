import sys
import os
#Importing issue with gpaw-python?
#However cannot use this for module importing
sys.path.append(os.path.join(os.path.dirname(__file__), "../../"))
import numpy
from bulk.build import StructureBuilder
from bulk.calc import MaterCalc
import unittest
from bulk.build import StructureBuilder
from ase.parallel import parprint, world
from ase.dft.bandgap import bandgap
from ase.dft.kpoints import special_paths, get_cellinfo
from ase.units import Ha

from gpaw import GPAW, PW, FermiDirac
from gpaw.xc.exx import EXX
from gpaw.xc.tools import vxc

class Test(unittest.TestCase):
    def test_bs(self):
        sb = StructureBuilder()
        atoms, *_ = sb.get_structure("C", "diamond")
        # print(atoms)
        base_dir = os.path.join(os.path.dirname(__file__),
                                "../../tmp/C-class/")
        m_calc = MaterCalc(atoms=atoms,
                           base_dir=base_dir)
        self.assertTrue(m_calc.relax(fmax=0.002))  # Very tight limit!
        self.assertTrue(m_calc.ground_state())
        # get the PBE BS
        lattice_type = get_cellinfo(m_calc.atoms.cell).lattice
        self.assertTrue(lattice_type in special_paths.keys())
        kpts_bs = dict(path=special_paths[lattice_type],
                       npoints=120)
        # HSE06 base generate
        gs_file = os.path.join(base_dir, "gs.gpw")
        _calc = GPAW(restart=gs_file)
        atoms = _calc.atoms.copy()
        calc = GPAW(**_calc.parameters)
        calc.set(kpts=dict(gamma=True,
                           density=4))  # low density calculations
        calc.atoms = atoms
        del _calc
        calc.get_potential_energy()
        calc.write(os.path.join(base_dir, "hse.gpw"), mode="all")
        calc = GPAW(restart=os.path.join(base_dir, "hse.gpw"),
                    txt=None)
        ns = calc.get_number_of_spins()
        nk = len(calc.get_ibz_k_points())
        nbands = calc.get_number_of_bands()
        eigen_pbe = numpy.array([[calc.get_eigenvalues(spin=s,
                                                       kpt=k) \
                                  for k in range(nk)]\
                                 for s in range(ns)])
        parprint("EIGEN_PBE", eigen_pbe.shape)
        vxc_pbe = vxc(calc, "PBE")
        parprint("VXC_PBE", vxc_pbe.shape)
        # world.barrier()
        # HSE06 now
        calc_hse = EXX(os.path.join(base_dir, "hse.gpw"),
                       xc="HSE06",
                       bands=[0, nbands])
        calc_hse.calculate()
        vxc_hse = calc_hse.get_eigenvalue_contributions()
        parprint(vxc_hse.shape)
        parprint(vxc_hse)
        eigen_hse = eigen_pbe - vxc_pbe + vxc_hse
        
        # HSE bandgap from just kpts
        bg_hse_min, *_ = bandgap(eigenvalues=eigen_hse,
                                 efermi=calc.get_fermi_level(),
                                 direct=False)
        bg_hse_dir, *_ = bandgap(eigenvalues=eigen_hse,
                                 efermi=calc.get_fermi_level(),
                                 direct=True)
        parprint("HSE: E_min \t E_dir")
        parprint("{:.3f}\t{:.3f}".format(bg_hse_min, bg_hse_dir))
        """
        # get the gllbsc by steps
        calc_ = GPAW(restart=gpw_name)
        calc_gllb = GPAW(**calc_.parameters)
        calc_gllb.set(xc="GLLBSC",
                      txt=os.path.join(base_dir, "gllb-gs.txt"))
        calc_gllb.atoms = calc_.atoms
        del calc_
        calc_gllb.get_potential_energy()  # SC calculation
        calc_gllb.write("gllb-gs.gpw")
        calc_gllb_bs = GPAW(restart="gllb-gs.gpw",
                            kpts=kpts_bs,
                            fixdensity=True,
                            symmetry="off",
                            txt=os.path.join(base_dir,"gllb-bs.txt"))
        world.barrier()
        calc_gllb_bs.get_potential_energy()
        homolumo = calc_gllb_bs.get_homo_lumo()
        bg_gllb_ks = homolumo[1] - homolumo[0]
        response = calc_gllb.hamiltonian.xc.xcs["RESPONSE"]
        response.calculate_delta_xc( homolumo / Ha)
        EKs, Dxc = response.calculate_delta_xc_perturbation()
        bg_gllb_deltaxc = EKs + Dxc

        ibz_kpts = calc_gllb_bs.get_ibz_k_points()
        e_kn = numpy.array([calc_gllb_bs.get_eigenvalues(kpt=k) \
                            for k in range(len(ibz_kpts))])
        efermi = calc_gllb_bs.get_fermi_level()
        e_kn[e_kn > efermi] += Dxc
        bg_gllb_min, *_ = bandgap(eigenvalues=e_kn,
                               efermi=efermi,
                               direct=False)
        bg_gllb_dir, *_  = bandgap(eigenvalues=e_kn,
                               efermi=efermi,
                               direct=True)

        parprint("PBE: E_min \t E_dir")
        parprint("{:.3f}\t{:.3f}".format(bg_pbe_min, bg_pbe_dir))
        parprint("Gllb: EKS \t E_deltaxc")
        parprint("{:.3f}\t{:.3f}".format(bg_gllb_ks, bg_gllb_deltaxc))
        parprint("Gllb: E_min \t E_dir")
        parprint("{:.3f}\t{:.3f}".format(bg_gllb_min, bg_gllb_dir))
        bs_gllb = calc_gllb_bs.band_structure()
        bs_gllb.energies[bs_gllb.energies > bs_gllb.reference] += Dxc
        bs_gllb.plot(emin=-10, emax=10,
                    filename=os.path.join(base_dir, "gllb-bs.png"))
            
        calc_gllb_bs.write(os.path.join(base_dir,
                                        "gllb-bs.gpw"))
        """

    
        
if __name__ == "__main__":
    unittest.main()
