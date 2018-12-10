import numpy
from bulk.build import StructureBuilder
from bulk.calc import MaterCalc
import os
from os.path import dirname, abspath, join, exists

def get_single_data(root_dir, formula, prototype):
    name = "{}-{}".format(formula, prototype)
    base_dir = os.path.join(root_dir, name)
    results = dict(formula=formula,
                   prototype=prototype,
                   name=name,
                   source="gpaw")
    bg_pbe_file = join(base_dir, "bg_pbe.npz")
    bg_gllb_file = join(base_dir, "bg_gllb.npz")
    eps_file = join(base_dir, "eps_rpa.npz")
    if all(map(exists, [bg_pbe_file, bg_gllb_file, eps_file])):  # all exists
        with numpy.load(bg_pbe_file) as f_pbe:
            results.update(gap_dir_pbe=float(f_pbe["Eg_dir"]),
                           gap_min_pbe=float(f_pbe["Eg_min"]))
            
        with numpy.load(bg_gllb_file) as f_gllb:
            results.update(gap_dir_gllb=float(f_gllb["Eg_dir"]),
                           gap_min_gllb=float(f_gllb["Eg_min"]))

        with numpy.load(eps_file) as f_eps:
            results.update(freq_rpa=f_eps["frequencies"],
                           eps_x=f_eps["eps_x"],
                           eps_y=f_eps["eps_y"],
                           eps_z=f_eps["eps_z"])
        print(name, "OK")
        return results
    else:
        return None

def main(root_dir="/cluster/scratch/ttian/bulk/"):
    sb = StructureBuilder()
    entries = sb.entries
    all_data = []
    for entry in entries:
        formula = entry["formula"]
        prototype = entry["prototype"]
        single_data = get_single_data(root_dir,
                                      formula,
                                      prototype)
        if single_data is not None:
            all_data.append(single_data)
    f_save = join(root_dir, "bulk_eps.npz")
    numpy.savez_compressed(f_save, data=all_data)
    print("done!")

if __name__ == "__main__":
    main()
