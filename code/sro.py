from utilities.python.tools import *
from utilities.python.abaqus_inp import *
from os import listdir, remove, mkdir, popen, _exit
from os.path import exists, dirname, abspath, join, isfile
from shutil import rmtree, copy, move
from argparse import ArgumentParser


def pre():
    log("")
    log("")

    log("Cleaning project directory", 1)
    for file_dir in listdir():
        if file_dir.endswith(".exe"):
            remove(file_dir)
        # if file_dir.endswith(".inp"):
        #     remove(file_dir)
        # if file_dir.endswith(".mtx"):
        #     remove(file_dir)
        # if file_dir.endswith(".info"):
        #     remove(file_dir)
        if file_dir == "_out":
            rmtree(file_dir)
        if file_dir == "_log":
            rmtree(file_dir)
        if file_dir == "_job":
            rmtree(file_dir)
        if file_dir == "_vtk":
            rmtree(file_dir)

    log("Done", 1)

    log("Updating folder structure", 1)
    mkdir("_out")
    mkdir("_log")
    mkdir("_job")
    mkdir("_vtk")
    log("Done", 1)


def post(inputfile, user, disbond):
    log("")
    log("")

    job = inputfile.split("/")[1]
    extension = ".inp"
    outdir = "_job/"
    input = job + extension
    inputfnm = job + "_fnm" + extension
    info = job + "_fnm.info"
    user = user + ".f"
    path = dirname(abspath(__file__))

    if disbond:
        jobda = job + "_da"
        inputda = jobda + extension

        log("Copying input files to _job directory", 2)
        copy(input, outdir + input)
        copy(inputda, outdir + inputda)
        log("Done", 2)

        log("Running job 1\n", 2)
        abaquscmd = "abaqus interactive job=" + job + " input=" + input
        ost = popen("cd _job && " + abaquscmd + " && cd ..")
        out = ost.readlines()
        log("Done with the following output:", 2)
        for ol in out:
            log(ol.strip(), 4)
        log("")
        mkdir(outdir + job)
        for file in listdir(outdir):
            item = join(join(path, outdir[:-1]), file)
            if isfile(item):
                if not (".mtx" in file or "_da.inp" in file):
                    move(outdir + file, join(outdir + job, file))

        log("Running job 2\n", 2)
        abaquscmd = "abaqus interactive job=" + jobda + " input=" + inputda
        ost = popen("cd _job && " + abaquscmd + " && cd ..")
        out = ost.readlines()
        log("Done with the following output:", 2)
        for ol in out:
            log(ol.strip(), 4)
        log("")
        mkdir(outdir + jobda)
        for file in listdir(outdir):
            item = join(join(path, outdir[:-1]), file)
            if isfile(item):
                if not (".mtx" in file):
                    move(outdir + file, join(outdir + jobda, file))

        log("Copying matrix files to root directory", 2)
        for file in listdir(outdir):
            if ".mtx" in file:
                copy(outdir + file, file)
        log("Done", 2)
        log("")
        log("")

    # if exists(info) and exists(inputfnm):
    #     log("Copying input files to _job directory", 2)
    #     copy(inputfnm, outdir + inputfnm)
    #     copy(info, outdir + info)
    #     log("Done", 2)

    #     abaquscmd = "abaqus interactive job=" + job + " input=" + \
    #         inputfnm + " user=..\\" + user + " scratch=\"%cd%\""
    #     log("Running job : " + abaquscmd + "\n\n", 2)
    #     ost = popen("cd _job && " + abaquscmd + " && cd ..")
    #     out = ost.readlines()
    #     log("Done with the following output:", 2)
    #     for ol in out:
    #         log(ol.strip(), 4)
    #     log("")
    #     log("")


def sro(filename, lsinit, opt, steps):
    log("")
    log("")
    log("Input file specificied: " + filename, 1)

    log("Attempting to open the input file", 1)
    inp = AbaqusInp(filename)
    log("File read successfully", 1)

    log("Applying disbond and printing {0}_da.inp and {0}.inp".format(
        filename), 2)
    disb = inp.apply_disbond()
    log("Done", 2)

    log("Applying FNM user elements and printing {0}_fnm.inp {0}_fnm.info".format(
        filename), 2)
    fnm = inp.apply_fnm(lsinit[0], lsinit[1], opt[0], opt[1], steps)
    log("Done", 2)

    return fnm, disb


def parse_input():
    ap = ArgumentParser(usage="python sro.py [-h/--help] -i/--input filename [-u/--user filename] [-l/--lsinit nX,nY] [-o/--opt Alpha,VolFrac] [-s/--steps N]",
                        description="Front-end CLI for the stringerRunoutOptimisation fortran Abaqus subroutines")

    ap.add_argument("-i", "--input", required=False, type=str, default="sub",
                    help="Filename of the abaqus input file (*.inp) residing in the _in/ directory.")
    ap.add_argument("-u", "--user", required=False, type=str, default="user",
                    help="Filename of the abaqus fortran subroutines file (*.f) residing in the base directory.")
    ap.add_argument("-l", "--lsinit",  required=False, type=str, default="0,0",
                    help="Hole pattern for the levelset initialisation in the form \"HolesX,HolesY\". (default = 0,0)")
    ap.add_argument("-o", "--opt",  required=False, type=str, default="0.5,0.5",
                    help="Parameters for the optimisation problem in the form \"ObjWeight,VolFraction\". (default = 0.5,0.5)")
    ap.add_argument("-s", "--steps", required=False, type=int, default=1,
                    help="Number of optimisation steps needed. (default=1)")

    args = ap.parse_known_args()[0]

    input = "_in/" + args.input
    user = args.user
    lsinit = []
    opt = []
    steps = 1
    if args.input.endswith(".inp"):
        input = args.input[:-4]
    if args.user.endswith(".f"):
        user = args.user[:-1]
    aux = args.lsinit.split(",")
    if not len(aux) == 2:
        log("Input error - lsinit option should be in the form \"HolesX,HolesY\"")
        _exit(0)
    else:
        lsinit = [int(x) for x in args.lsinit.split(",")]
    aux = args.opt.split(",")
    if not len(aux) == 2:
        log("Input error - lsinit option should be in the form \"ObjWeight,VolFraction\"")
        _exit(0)
    else:
        opt = [float(x) for x in args.opt.split(",")]
    if args.steps > 0 and args.steps < 5000:
        steps = args.steps
    else:
        log("Input error - steps should be an integer between 1 and 4999")
        _exit(0)

    return input, user, lsinit, opt, steps


if __name__ == "__main__":
    log("")
    log("")
    log("Stringer Run-out Optimisation - START")

    input_file, code_file, lsinit, opt, steps = parse_input()

    pre()
    fnm, disb = sro(input_file, lsinit, opt, steps)
    if fnm:
        post(input_file, code_file, disb)
    else:
        log("Fnm generation did not complete. SRO will not go further.", 1)

    log("Stringer Run-out Optimisation - END")
    log("")
    log("")
