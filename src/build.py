import ase.db
from ase.atoms import Atoms, string2symbols
from ase.data import atomic_numbers, covalent_radii
from ase.build import bulk
from ase.parallel import parprint, world, broadcast
import os, os.path

# Construct the bulk material using given crystalline geometry
def get_structure(formula, kind=None):
    # Use ase.build.bulk
    try:
        symbols = string2symbols(formula)
    except Exception:
        return None
    guess_a = sum([covalent_radii[atomic_numbers[s]] for s in symbols]) * (3 ** 0.5)
    try:
        mater = bulk(formula, kind, a=guess_a)
    except Exception:
        return None
    # Material ok?
    mater.set_pbc((True, True, True))
    return mater

    

if __name__ == "__main__":
    import os, sys
    if len(sys.argv) < 2:
        raise ValueError("not enough parameters")
    else:
        formula = sys.argv[1]
        print(get_structure(formula))
