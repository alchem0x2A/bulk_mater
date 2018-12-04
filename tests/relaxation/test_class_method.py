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
    def test_calc_single(self):
        sb = StructureBuilder()
        atoms, *_ = sb.get_structure("Si", "diamond")
        print(atoms)
        m_calc = MaterCalc(atoms=atoms,
                           base_dir="../../tmp/Si-class/")
        res = m_calc.relax(fmax=0.002)  # Very tight limit!
        self.assertTrue(res)
if __name__ == "__main__":
    unittest.main()
