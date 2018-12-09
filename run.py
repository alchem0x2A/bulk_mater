import sys
import os
from os.path import exists, join, dirname, abspath
# May need this for the path issue
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from bulk.build import StructureBuilder, convert_name
from bulk.calc import MaterCalc
from ase.parallel import paropen, parprint, world, rank, broadcast
import shutil

# run single job
def run_single(formula, prototype, 
               root="/cluster/scratch/ttian/bulk",
               clean=False):
    prototype = convert_name(prototype)[0]
    name = "{}-{}".format(formula, prototype)
    base_dir = join(root, name)
    # Directory manipulation
    if rank == 0:    
        if clean:
            shutil.rmtree(base_dir, ignore_errors=True)
                
        if not exists(base_dir):
            os.makedirs(base_dir)
    world.barrier()
    
    sb = StructureBuilder()
    atoms, *_ = sb.get_structure(formula, prototype)
    m_calc = MaterCalc(atoms=atoms, base_dir=base_dir)
    m_calc.relax(fmax=0.002)
    m_calc.ground_state()
    eg_min, eg_dir, *_ = m_calc.bandgap(method="pbe")
    parprint("PBE min/dir: {:.3f}\t{:.3f}".format(eg_min, eg_dir))
    eg_min, eg_dir, *_ = m_calc.bandgap(method="gllb")
    parprint("GLLB min/dir: {:.3f}\t{:.3f}".format(eg_min, eg_dir))
    m_calc.excited_state()
    m_calc.dielectric(method="rpa")

    return 0

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise ValueError("Please provide at least 1 parameters")

    formula = sys.argv[1]
    if len(sys.argv) == 2:
        proto = None
        run_single(formula, proto)
    elif len(sys.argv) >= 3:
        if sys.argv[2] == "clean":
            run_single(formula, clean=True)
        else:
            proto = sys.argv[2]
            if proto.lower() in ("none",
                                "other",
                                "unknown"):
                proto = "other"
            run_single(formula, proto)
    else:
        raise ValueError("Parameter ill defined!")
