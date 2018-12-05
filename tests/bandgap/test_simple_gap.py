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
from ase.parallel import parprint
from gpaw import GPAW, PW, FermiDirac
from ase.dft.bandgap import bandgap

class Test(unittest.TestCase):
    def test_relax_single(self):
        sb = StructureBuilder()
        atoms, *_ = sb.get_structure("C", "diamond")
        # print(atoms)
        m_calc = MaterCalc(atoms=atoms,
                           base_dir="../../tmp/C-class/")
        m_calc.relax(fmax=0.002)  # Very tight limit!
        m_calc.ground_state()
        
        # self.assertTrue(res)
        base_dir = "../../tmp/C-class/"
        gpw_name = os.path.join(base_dir, "gs.gpw")
        self.assertTrue(os.path.exists(gpw_name))
        calc = GPAW(restart=gpw_name, txt="gp.txt")
        # PBE bandgap
        bg_pbe_min, *_ = bandgap(calc, direct=False)
        bg_pbe_dir, *_ = bandgap(calc, direct=True)
        # gllbsc
        calc_gllb = GPAW(**calc.parameters)
        calc_gllb.atoms = calc.atoms
        calc_gllb.set(xc="GLLBSC", txt="gllb.txt")
        calc_gllb.get_potential_energy()  # SC calculation
        response = calc_gllb.hamiltonian.xc.xcs["RESPONSE"]
        response.calculate_delta_xc()
        EKs, Dxc = response.calculate_delta_xc_perturbation()
        gllb_gap = EKs + Dxc

        parprint("Eg-PBE-min, Eg-PBE-dir, Eg-gllb")
        parprint(bg_pbe_min, bg_pbe_dir, gllb_gap)
        # Use kpts?
        ibz_kpts = calc.get_ibz_k_points()
        e_kn = numpy.array([calc.get_eigenvalues(kpt=k) \
                            for k in range(len(ibz_kpts))])
        efermi = calc.get_fermi_level()
        e_kn[e_kn > efermi] += Dxc
        gllb_gap_min = bandgap(eigenvalues=e_kn,
                               efermi=efermi,
                               direct=False)
        gllb_gap_dir = bandgap(eigenvalues=e_kn,
                               efermi=efermi,
                               direct=True)
        
        parprint("Efermi", efermi)
        parprint("Eg-gllb-min, Eg-gllb-dir")
        parprint(gllb_gap_min, gllb_gap_dir)


        
if __name__ == "__main__":
    unittest.main()
