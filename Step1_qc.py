#!/usr/bin/env python3

# In the config_parameters.txt, you need add the number of parallel
# jobs. e.g. PARALLEL_JOB = 10


import sys
import argparse as ap
import itertools as it
import subprocess as sub
import multiprocessing as mul


class THREAD:
    # Number of parallel jobs
    GROUP = 2


def base(reads_1):
    if reads_1.endswith(".gz"):
        out = reads_1.split('.')
        out = '.'.join(out[:len(out) - 2])
    else:
        out = reads_1.split('.')
        out = '.'.join(out[:len(out) - 1])
    if out.endswith("_R1"):
        out = out.rstrip('1').rstrip('R').rstrip('_')
    elif out.endswith("_1"):
        out = out.rstrip('1').rstrip('_')
    elif out.endswith(".R1"):
        out = out.rstrip('1').rstrip('R').rstrip('.')
    elif out.endswith(".1"):
        out = out.rstrip('1').rstrip('.')
    return out


def to2(reads_1):
    if reads_1.endswith(".gz"):
        out = reads_1.split('.')
        prefix = '.'.join(out[:len(out) - 2])
        suffix = '.'.join(out[len(out) - 2:len(out)])
    else:
        out = reads_1.split('.')
        prefix = '.'.join(out[:len(out) - 1])
        suffix = '.'.join(out[len(out) - 1:len(out)])
    if prefix.endswith("R1"):
        prefix = prefix.replace("R1", "R2")
        out = '.'.join([prefix, suffix])
    elif prefix.endswith(".1"):
        prefix = prefix.replace(".1", ".2")
        out = '.'.join([prefix, suffix])
    elif prefix.endswith("_1"):
        prefix = prefix.replace("_1", "_2")
        out = '.'.join([prefix, suffix])
    return out


def fastqc(reads_1, reads_2, threads='2'):
    print("fastqc", "-t", threads, reads_1, reads_2)
    sub.call(["fastqc", "-t", threads, reads_1, reads_2])


def fastq_report(sample, species):
    d = "'" + base(sample) + "/1'"
    rmd = "-e \"rmarkdown::render('fastq_report.Rmd',output_dir=" + \
          d + ")\""
    command = ' '.join(["Rscript", rmd, sample, species])
    print(command)
    sub.call(command, shell=True)
    # Reads_2
    d = "'" + base(sample) + "/2'"
    rmd = "-e \"rmarkdown::render('fastq_report.Rmd',output_dir=" + \
          d + ")\""
    command = ' '.join(["Rscript", rmd, to2(sample), species])
    print(command)
    sub.call(command, shell=True)


# Help function for parallel running
def helper_fastqc(args):
    fastqc(args[0], args[1])


def helper_fastq_report(args):
    fastq_report(args[0], args[1])


def setting():
    parser = \
        ap.ArgumentParser(description=
                          "Automatic RNA-seq paired-end fastq QC pipeline.")
    parser.add_argument("-F", "--list_R1", required=True,
                        help="List of R1 Fastq files.")
    parser.add_argument("-R", "--list_R2",
                        help="List of R2 Fastq files. [Optional]")
    parser.add_argument("-S", "--species", default="Hsapiens",
                        help="Name of species from BSgenome packages. "
                             "e.g. Hsapiens, Celegans [Default: Hsapiens]")
    parser.add_argument("-P", "--parameter_file", type=str,
                        default="./config_parameters.txt",
                        help="Some constant parameters. See example. "
                             "[Default: ./config_parameters.txt]")
    parser.add_argument("-v", "--version", action="version",
                        version="2019-08.", help="Show version")
    a = parser.parse_args()

    if a.list_R2 is None:
        a.list_R2 = -1
    if a.species not in ["Btaurus", "Celegans", "Cfamiliaris", "Dmelanogaster",
                         "Drerio", "Hsapiens", "Mmulatta", "Mmusculus",
                         "Ptroglodytes", "Rnorvegicus", "Scerevisiae",
                         "Sscrofa"]:
        print("We temporarily don't support this species for GC distribution "
              "... use Hsapiens instead!")
        a.species = "Hsapiens"
    if a.parameter_file is None:
        print("Wrong parameter file!")
        sys.exit(1)

    return a


def main():
    load = setting()
    if load.list_R2 == -1:
        with open(load.list_R1, 'r') as fr:
            r_list = [to2(x.strip()) for x in fr.readlines()]
    else:
        with open(load.list_R2, 'r') as fr:
            r_list = [x.strip() for x in fr.readlines()]
    with open(load.list_R1, 'r') as fr:
        f_list = [x.strip() for x in fr.readlines()]
    if len(f_list) != len(r_list):
        print("Different length of R1/R2 fastq files!")
        sys.exit(1)

    # Processing parameter file
    with open(load.parameter_file, 'r') as fr:
        lines = [x.strip().replace(' ', '') for x in fr.readlines()]
        for i in lines:
            if not i.startswith('#') and i != '':
                opt = i.split('=')
                if opt[0] == "PARALLEL_JOB":
                    THREAD.GROUP = int(opt[1])

    # Parallel running
    pool = mul.Pool(THREAD.GROUP)
    print("Parallel running", THREAD.GROUP, "jobs ...")
    pool.map(helper_fastqc, it.zip_longest(f_list, r_list))
    pool.map(helper_fastq_report,
             it.zip_longest(f_list, it.repeat(load.species, len(f_list))))


if __name__ == "__main__":
    main()
