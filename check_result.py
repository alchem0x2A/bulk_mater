from bulk.build import StructureBuilder
from bulk.calc import MaterCalc
# from run import run_single
from subprocess import run
import os

root = "/cluster/scratch/ttian/bulk/"
def sub_job(formula, prototype, ncores=24):
    sub_string = ("bsub -n {core} -W 24:00 -R "
                  "\"rusage[mem=2048]\" -J {jobname} "
                  "\"mpirun -n {core} gpaw-python run.py {form} {proto}\"")
    name = "{}-{}".format(formula, prototype)
    proc = run(sub_string.format(core=ncores,
                          jobname=name,
                          form=formula,
                          proto=prototype),
        shell=True,
        timeout=60)
    return proc.returncode

def main():
    sb = StructureBuilder()
    entries = sb.entries
    finished = 0
    for entry in entries:
        formula = entry["formula"]
        prototype = entry["prototype"]
        base_dir = os.path.join(root,
                                "{}-{}".format(formula, prototype))
        mc = MaterCalc(atoms=None, base_dir=base_dir)
        if all([value for key, value in mc.check_status().items()]):
            finished += 1
        else:
            print("{}-{}:".format(formula, prototype),
                  mc.check_status())
    print("Finished jobs {}".format(finished))
    return True
    
if __name__ == "__main__":
    main()
