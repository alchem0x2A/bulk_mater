from __future__ import print_function
import sys
import os
#Importing issue with gpaw-python?
#However cannot use this for module importing
sys.path.append(os.path.join(os.path.dirname(__file__), "../../"))
import pkgutil
import unittest
from bulk.build import StructureBuilder
from ase.parallel import parprint

class TestPath(unittest.TestCase):
    def test_build_entries(self):
        sb = StructureBuilder()
        entries = sb.entries
        # [print(l) for l in entries]
        self.assertGreater(len(entries), 0)
    def test_simple_build(self):
        sb = StructureBuilder()
        entries = sb.entries
        failed = []
        for e in entries:
            try:
                mater = sb._convert_struct(e)
            except Exception as err:
                failed.append((e["formula"],
                               e["prototype"],
                               err))
        if len(failed) > 0:
            parprint(len(failed), failed)
        else:
            parprint("All structure passed!")
        self.assertEqual(len(failed), 0)
    def test_build_by_name(self):
        sb = StructureBuilder()
        test_samples = [("AlAs", "Z"),
                        ("ZnO", "z"),
                        ("ZnO", "W"),
                        ("ZnO", "R"),
                        ("TlCl", "C"),
                        ("HgO", None),
                        ("Si3N4", "other"),
                        ("Si3N4", "unknown"),
                        ("C", "D")]
        failed = []
        for s in test_samples:
            formula, prototype = s
            try:
                mat = sb.get_structure(formula, prototype)
                if len(mat) == 0:
                    failed.append((formula, prototype, "zero-len"))
            except Exception as e:
                failed.append((formula, prototype, e))
        if len(failed) > 0:
            parprint(len(failed), failed)
        else:
            parprint("All structure passed!")
        self.assertEqual(len(failed), 0)
if __name__ == "__main__":
    unittest.main()
