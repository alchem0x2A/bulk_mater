import sys
import os
#Importing issue with gpaw-python?
#However cannot use this for module importing
sys.path.append(os.path.join(os.path.dirname(__file__), "../../"))
import numpy
from bulk.build import StructureBuilder
from bulk.calc import MaterCalc
import unittest
from ase.parallel import parprint

class Test(unittest.TestCase):
    def test_relax_single(self):
        sb = StructureBuilder()
        atoms, *_ = sb.get_structure("C", "diamond")
        base_dir = os.path.join(os.path.dirname(__file__),
                                "../../tmp/Si-module-test/")
        m_calc = MaterCalc(atoms=atoms, base_dir=base_dir)
        self.assertTrue(m_calc.relax(fmax=0.002))
        self.assertTrue(m_calc.ground_state())
        # PBE
        Eg_min, Eg_dir, bs = m_calc.bandgap("PBE")
        with open("PBE-gap.npz", "wb") as f:
            numpy.savez(f, Eg_min=Eg_min,
                        Eg_dir=Eg_dir,
                        **bs)
        parprint("PBE Gap: min \t dir")
        parprint("{:.3f}\t{:.3f}".format(Eg_min, Eg_dir))
        parprint("Bandstructure PBE:", bs)
        # GLLB
        Eg_min, Eg_dir, bs = m_calc.bandgap("GLLB")
        parprint("GLLB Gap: min \t dir")
        parprint("{:.3f}\t{:.3f}".format(Eg_min, Eg_dir))
        parprint("Bandstructure PBE:", bs)
        
        


        
if __name__ == "__main__":
    unittest.main()
