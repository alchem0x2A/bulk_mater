from __future__ import print_function
import sys
import os
#Importing issue with gpaw-python?
#However cannot use this for module importing
sys.path.append(os.path.join(os.path.dirname(__file__), "../../"))
import pkgutil
import unittest
import bulk

class TestPath(unittest.TestCase):
    def test_module(self):
        def sub_module(package):
            submod = []
            for importer, modname, ispkg in pkgutil.iter_modules(package.__path__):
                print("Found submodule {} (is a package: {})".format(modname,
                                                                     ispkg))
                submod.append(modname)
            return submod
        submod = sub_module(bulk)
        self.assertGreater(len(submod), 0)

if __name__ == "__main__":
    unittest.main()
