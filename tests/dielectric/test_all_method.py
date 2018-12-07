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
    def test_single(self):
        sb = StructureBuilder()
        atoms, *_ = sb.get_structure("Si", "diamond")
        base_dir = os.path.join(os.path.dirname(__file__),
                                "../../tmp/Si-eps-test/")
        m_calc = MaterCalc(atoms=atoms, base_dir=base_dir)
        self.assertTrue(m_calc.relax(fmax=0.002))
        self.assertTrue(m_calc.ground_state())
        # PBE
        Eg_min, Eg_dir, bs = m_calc.bandgap(method="PBE")
        parprint(Eg_min, Eg_dir)
        self.assertLess(abs(Eg_min - 0.61), 0.15)
        Eg_min, Eg_dir, bs = m_calc.bandgap(method="GLLB")
        parprint(Eg_min, Eg_dir)
        self.assertLess(abs(Eg_min - 1.1), 0.15)
        self.assertTrue(m_calc.excited_state())
        self.assertTrue(m_calc.dielectric())
        
        


        
if __name__ == "__main__":
    unittest.main()
