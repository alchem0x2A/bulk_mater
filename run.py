import sys
import os, os.path
# May need this for the path issue
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from src.build import get_structure
from src.relax import relax
from src.bandgap import gap
from src.dielectric import excited, permittivity
import shutil
from ase.parallel import paropen, parprint, world, rank, broadcast

def run_single(formula, kind=None, 
         root="/cluster/scratch/ttian/bulk",
         clean=False):
    # candidates = {}
    # if rank == 0:
    mol = get_structure(formula, cs=kind)
    if mol is None:
        return False
    name = "{}-{}".format(formula, kind)

    # Directory manipulation
    if rank == 0:    
        base_dir = os.path.join(root, name)
        if clean:
            shutil.rmtree(base_dir, ignore_errors=True)
                
        if not os.path.exists(base_dir):
            os.makedirs(base_dir)
        
    world.barrier()
    if clean:
        return                  # on all ranks

    # On all ranks
    base_dir = os.path.join(root, name)
    # Relaxation and gs
    relax(mol, name=name,
          base_dir=base_dir)
    parprint("Relaxation for {} finished!".format(name))

    gap(base_dir=base_dir, mode="gllb")
    parprint("Bandgap for {} calculated!".format(name))
    
    excited(base_dir=base_dir)
    parprint("Excitation for {} finished!".format(name))
    
    permittivity(base_dir=base_dir, mode="df")
    parprint("Permittivity for {} finished!".format(name))
    
        # polarizability(base_dir=base_dir, mode="tetra")  # 
        # parprint("Polarizability using tetra {} finished!".format(name))
    return 0

if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise ValueError("Please provide 2 parameters")
    elif len(sys.argv) == 2:
        formula = sys.argv[1]
        main(formula)
    elif len(sys.argv) == 3:
        formula = sys.argv[1]
        if sys.argv[2] == "clean":
            main(formula, clean=True)
        else:
            kind = sys.argv[2]
            main(formula, kind)
    else:
        raise ValueError("Parameter ill defined!")
