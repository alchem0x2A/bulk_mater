import sys
import os
#Importing issue with gpaw-python?
#However cannot use this for module importing
sys.path.append(os.path.join(os.path.dirname(__file__), "../../"))

from bulk.build import StructureBuilder
from bulk.calc import MaterCalc
import unittest
from bulk.build import StructureBuilder
from ase.parallel import parprint

class Test(unittest.TestCase):
    def test_gs_single(self):
        sb = StructureBuilder()
        atoms, *_ = sb.get_structure("Si", "diamond")
        print(atoms)
        m_calc = MaterCalc(atoms=atoms,
                           base_dir="../../tmp/Si-class/")
        self.assertTrue(m_calc.relax(fmax=0.002))
        self.assertTrue(m_calc.ground_state())
if __name__ == "__main__":
    unittest.main()
