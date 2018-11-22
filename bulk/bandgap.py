import os, os.path
import numpy
import json
import ase.db
from gpaw import GPAW, PW, FermiDirac
from ase.parallel import parprint, rank
import numpy as np

'''
Currently gllb-sc bandgap calculated
'''
def gap(base_dir="./", mode="gllb"):
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    param_file = os.path.join(curr_dir, "../parameters.json")
    gpw_file = os.path.join(base_dir, "gs.gpw")
    gap_res = os.path.join(base_dir, "gap.txt")
    # Check is gap calculated
    if os.path.exists(gap_res):
        parprint("Bandgap calculated, will use gpw directly!")
        return 0
    # Check if gpw file exists
    if not os.path.exists(gpw_file):
        raise FileNotFoundError("Structure relaxation not done yet!")
    if os.path.exists(param_file):
        params = json.load(open(param_file, "r"))
    else:
        raise FileNotFoundError("no parameter file!")
    # The true gllb-sc calculation
    calc = GPAW(gpw_file,
                **params["gap"])
    calc.get_potential_energy()

    # Now for bg
    response = calc.hamiltonian.xc.xcs["RESPONSE"]
    response.calculate_delta_xc()
    EKs, Dxc = response.calculate_delta_xc_perturbation()
    gap = EKs + Dxc
    if rank == 0:
        print("{:.4f}".format(gap),
              open(gap_res, "w"))  # write answer
    return 0
    
    
