#/usr/bin/env python3

# Combine fastq files, different file groups should be separated by empty line.
# e.g.
# ---
# abc_1_R1.fastq.gz
# abc_2_R1.fastq.gz
#
# xyz_1_R2.fastq.gz
# xyz_2_R2.fastq.gz
# ---
# This file list will cause:
# abc_1_R1.fastq.gz + abc_2_R1.fastq.gz -> combined_abc_1_R1.fastq.gz
# xyz_1_R2.fastq.gz + xyz_2_R2.fastq.gz -> combined_xyz_1_R2.fastq.gz


import os
import argparse as ap
import subprocess as sub


def setting():
    parser = \
        ap.ArgumentParser(description =
                          "Automatic RNA-seq paired-end fastq QC pipeline.")
    parser.add_argument("-F", "--reads_list", required = True,
                        help = "List of Fastq files need to be combined.")
    parser.add_argument("-v", "--version", action = "version",
                        version = "2019-08.", help = "Show version")
    a = parser.parse_args()

    return(a)

def main():
    load = setting()
    with open(load.reads_list, 'r') as fr:
        temp = []
        for line in fr.readlines():
            line = line.strip()
            if line != '':
                temp.append(line)
            else:
                sub.call([' '.join(["cat",' '.join(temp), '>',
                                    os.path.dirname(temp[0]) + "/combined_" +
                                    os.path.basename(temp[0])])], shell = True)
                temp = []

    sub.call([' '.join(["cat",' '.join(temp), '>',
                        os.path.dirname(temp[0]) + "/combined_" +
                        os.path.basename(temp[0])])], shell = True)


if __name__ == "__main__":
    main()
