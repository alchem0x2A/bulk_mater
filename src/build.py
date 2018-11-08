import ase.db
from ase.atoms import Atoms, string2symbols
from ase.io import read
from ase.data import atomic_numbers, covalent_radii
from ase.build import bulk
from ase.parallel import parprint, world, broadcast, rank
import json
import os, os.path

json_file = os.path.abspath("../structure/prototype.json")
data = json.load(open(json_file, "r"))

def conver_name(cs):
    cs = cs.lower()
    if cs[0] == "z":
        return "zincblende"
    elif cs[0] == "w":
        return "wurtzite"
    elif cs[0] in ("n", "r"):
        return "rocksalt"
    elif cs[0] == "c":
        return "cesiumchloride"
    else:
        return "other"
        
# Construct the bulk material using given crystalline geometry
def get_structure(formula, cs=None):
    # Use ase.build.bulk
    if cs is None:
        cs = "other"
    else:
        cs = conver_name(cs)
    
    item = data[cs][formula]
    if isinstance(item, dict):  # wurtzite
        mater = bulk(name=formula,
                     crystalstructure=cs,
                     a=item["a"],
                     c=item["c"])
    elif isinstance(item, float):  # only a
        mater = bulk(name=formula,
                     crystalstructure=cs,
                     a=item)
    elif isinstance(item, str):  # other structure
        cif_file = os.path.abspath(os.path.join("../structure",
                                                item))
        mater = read(cif_file)
    else:
        mater = None
    return mater

    

if __name__ == "__main__":
    import os, sys
    if len(sys.argv) < 3:
        raise ValueError("not enough parameters")
    else:
        formula = sys.argv[1]
        cs = sys.argv[2]
        print(get_structure(formula, cs))
