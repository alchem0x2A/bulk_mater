import numpy
import os
import matplotlib.pyplot as plt
from ase.dft.band_structure import BandStructure, BandStructurePlot


data = numpy.load("./PBE-gap.npz")

bs = BandStructure(cell=data["cell"],
                   kpts=data["kpts"],
                   energies=data["energies"],
                   reference=data["reference"])
print(bs)
bsp = BandStructurePlot(bs)
bsp.plot(emin=-8, emax=8, filename="bs.png")
