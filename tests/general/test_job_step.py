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
        # PBE
        parprint(m_calc.check_status())
        
        


        
if __name__ == "__main__":
    unittest.main()
