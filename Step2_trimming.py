#!/usr/bin/env python3

# In the config_parameters.txt, you need add the number of parallel
# jobs. e.g. PARALLEL_JOB = 10

# In the config_parameters.txt, you need add the path of Trimmomatic jar
# jobs. e.g. TRIMMOMATIC = /home/xxx/trimmomatic.jar


import os
import sys
import argparse as ap
import itertools as it
import subprocess as sub
import multiprocessing as mul


class THREAD:
    # Number of parallel jobs
    GROUP = 2

class PARAMETERS:
    OPTIONS = {}


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
    return(out)

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
    return(out)

def trimming(reads_1, reads_2, num_5, num_3, adapter, threads):
    if adapter != None:
        adapter = ':'.join(["ILLUMINACLIP", adapter, "2:30:10"])
    else:
        adapter = ''
    # Shell pipe has much better performance for file reading
    cmd = ' '.join(["zcat", reads_1, "| head -n 2 | tail -n 1"])
    # Alternative method:
    # read = sub.getoutput(cmd)
    # keep_3 = str(len(read.splitlines()[1]) - int(num_3))
    read = sub.run(cmd, shell = True, stdout = sub.PIPE, stderr = sub.STDOUT)
    keep_3 = str(len(read.stdout) - num_3)
    # Path
    options = ' '.join(["PE", "-threads", str(threads), reads_1, reads_2,
                        '_'.join([base(reads_1), "paired_R1.fastq.gz"]),
                        '_'.join([base(reads_1), "unpaired_R1.fastq.gz"]),
                        '_'.join([base(reads_1), "paired_R2.fastq.gz"]),
                        '_'.join([base(reads_1), "unpaired_R2.fastq.gz"]),
                        adapter, ':'.join(["HEADCROP", str(num_5)]),
                        ':'.join(["CROP", keep_3]), ':'.join(["LEADING", "30"]),
                        ':'.join(["TRAILING", "30"]),
                        ':'.join(["MINLEN", "36"])])
    command = ' '.join(["java", "-jar", PARAMETERS.OPTIONS["TRIMMOMATIC"],
                        options])
    print(command)
    sub.call(command, shell = True)

def helper_trimming(args):
    trimming(args[0], args[1], args[2], args[3], args[4], args[5])

def setting():
    parser = \
        ap.ArgumentParser(description =
                          "Automatic RNA-seq paired-end fastq trimming "
                          "pipeline.")
    parser.add_argument("-F", "--list_R1", required = True,
                        help = "List of R1 Fastq files.")
    parser.add_argument("-R", "--list_R2",
                        help = "List of R2 Fastq files. [Optional]")
    parser.add_argument("-5", "--cut_5", default = 0, type = int,
                        help = "Numbers of bases to cut from 5' head. "
                               "[Default: 0]")
    parser.add_argument("-3", "--cut_3", default = 0, type = int,
                        help = "Numbers of bases to cut from 3' head. "
                               "[Default: 0]")
    parser.add_argument("-A", "--adapter", default = None, type = str,
                        help = "Path of adapter file for Trimmomatic. "
                               "[Optinal]")
    parser.add_argument("-T", "--threads", default = 2, type = int,
                        help = "Number of threads for each job. [Default: 2]")
    parser.add_argument("-P", "--parameter_file", type=str,
                        default="./config_parameters.txt",
                        help="Some constant parameters. See example. "
                             "[Default: ./config_parameters.txt]")
    parser.add_argument("-v", "--version", action = "version",
                        version = "2019-08.", help = "Show version")
    a = parser.parse_args()

    if a.threads < 1:
        print("Positive number for threads!")
        sys.exit(1)
    if a.parameter_file == None:
        print("Wrong parameter file!")
        sys.exit(1)

    return(a)

def main():
    load = setting()
    if load.list_R2 == None:
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
        lines = [x.strip().replace(' ','') for x in fr.readlines()]
        for i in lines:
            if not i.startswith('#') and i != '':
                opt = i.split('=')
                if opt[0] != "PARALLEL_JOB":
                    PARAMETERS.OPTIONS[opt[0]] = opt[1]
                else:
                    THREAD.GROUP = int(opt[1])

    # Parallel running
    pool = mul.Pool(THREAD.GROUP)
    print("Parallel running", THREAD.GROUP, "jobs ...")
    pool.map(helper_trimming,
             it.zip_longest(f_list, r_list,
                            it.repeat(load.cut_5, len(f_list)),
                            it.repeat(load.cut_3, len(f_list)),
                            it.repeat(load.adapter, len(f_list)),
                            it.repeat(load.threads, len(f_list))))

if __name__ == "__main__":
    main()
