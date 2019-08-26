#!/usr/bin/env python3

# In the config_parameters.txt, you need add some path for tools. e.g.
# STAR_REF = /home/xxx/references/star_ref
# RSEM_REF = /home/xxx/references/rsem_ref
# HISAT2_REF = /home/xxx/references/hisat2_ref
# GTF = /home/xxx/references/ref.gtf
# KALLISTO_REF = /home/xxx/references/kallisto_ref
# BED = /home/xxx/references/ref.bed --- BED 12 file
# SHELL = /bin/bash

# In the config_parameters.txt, you need add the number of parallel
# jobs. e.g. PARALLEL_JOB = 10


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
    BAM_LIST = []


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


def do_rsem(bamfile, threads, strand, sample_name):
    options = ' '.join(["--alignments", "--paired-end", "-p", threads,
                        "--strandedness", strand])
    command = ' '.join(["rsem-calculate-expression", options, bamfile,
                        PARAMETERS.OPTIONS["RSEM_REF"], sample_name])
    print(command)
    sub.call(command, shell=True)


def do_featureCounts(bamfiles, threads, strand, sample_name):
    bamfiles = ' '.join(bamfiles)
    gtf = PARAMETERS.OPTIONS["GTF"]
    # Strandedness
    if strand == "none":
        sp = ''
    elif strand == "forward":
        sp = "-s 1"
    else:
        sp = "-s 2"
    options = ' '.join(["-T", threads, "-a", gtf, "-g", "gene_id", "-t", "exon",
                        "-p", "-o", sample_name, sp])
    command = ' '.join(["featureCounts", options, bamfiles])
    if not os.path.exists(os.path.dirname(sample_name)):
        os.makedirs(os.path.dirname(sample_name))
    print(command)
    sub.call(command, shell=True)


def star_quantify(reads_1, reads_2, strand, threads):
    # Compressed fastq option
    if reads_1.endswith(".gz"):
        com = " --readFilesCommand zcat"
    else:
        com = ""
    # ENCODE project options
    encode_option = " --outFilterType BySJout --outFilterMultimapNmax 20 " \
                    "--alignSJoverhangMin 8 --alignSJDBoverhangMin 12 " \
                    "--outFilterMismatchNmax 999 " \
                    "--outFilterMismatchNoverReadLmax 0.04 " \
                    "--alignIntronMin 20 --alignIntronMax 100000 " \
                    "--alignMatesGapMax 100000"
    # Acceleration for multiprocessing, NoSharedMemory is the only option
    # compatible with --twopassMode Basic
    shared_mem = " --genomeLoad NoSharedMemory"
    # output dir
    out = " --outFileNamePrefix " + base(reads_1) + "/star/"
    if not os.path.exists(base(reads_1) + "/star/"):
        os.makedirs(base(reads_1) + "/star/")
    # Cufflinks compatibility
    cufflinks = " --outSAMstrandField intronMotif" + \
                " --outFilterIntronMotifs RemoveNoncanonical"
    # File type of output
    ot = " --outSAMtype BAM SortedByCoordinate"
    # Chimeric options, from STAR-Fusion
    chim = " --chimOutType Junctions" + " --chimSegmentMin 12" + \
           " --chimMultimapScoreRange 10" + " --chimMultimapNmax 10" + \
           " --chimNonchimScoreDropMin 10" + " --chimJunctionOverhangMin 12" + \
           " --chimSegmentReadGapMax 3" + \
           " --alignSJstitchMismatchNmax 5 -1 5 5" + \
           " --outSAMunmapped Within" + " --outSAMattrRGline ID:GRPundef" + \
           " --chimOutJunctionFormat 1"
    # Quantification mode
    qm = " --quantMode TranscriptomeSAM GeneCounts"
    # 2-pass mapping
    tm = " --twopassMode Basic"
    # Merging and mapping of overlapping paired-end reads
    po = " --peOverlapNbasesMin 12" + " --peOverlapMMp 0.1"
    # Allocating more space for reads alignment
    atp = " --alignTranscriptsPerReadNmax 20000"

    com = "STAR --runThreadN " + threads + " --genomeDir " + \
          PARAMETERS.OPTIONS["STAR_REF"] + " --readFilesIn " + reads_1 + " " + \
          reads_2 + encode_option + shared_mem + out + cufflinks + ot + chim + \
          qm + tm + po + atp + com
    print(com)
    sub.call(com, shell=True)
    # Quality Checking
    sub.call([PARAMETERS.OPTIONS["SHELL"], "BamQC.sh",
              base(reads_1) + "/star/Aligned.sortedByCoord.out.bam",
              PARAMETERS.OPTIONS["BED"]])
    # Quantification
    if not os.path.exists(base(reads_1) + "/rsem/"):
        os.makedirs(base(reads_1) + "/rsem/")
    do_rsem(base(reads_1) + "/star/" + "Aligned.toTranscriptome.out.bam",
            threads, strand,
            base(reads_1) + "/rsem/" + os.path.basename(base(reads_1)))


def hisat2_quantify(reads_1, reads_2, strand, threads):
    # Choose reference
    ref = PARAMETERS.OPTIONS["HISAT2_REF"]
    # Specify strand-specific information. For Illumina, reads 1 is first
    # strand corresponding to R; reads 2 if second strand corresponding to F.
    # So, we can use 'unstrand' for normal library and 'RF' for paired-end and
    # Illumina strand-specific library.
    sp = {"none": "", "forward": "--rna-strandness FR",
          "reverse": "--rna-strandness RF"}
    sp = sp[strand]
    # Path
    options = ' '.join(["--phred33", sp, "--dta", "-t",
                        "--new-summary", "--summary-file",
                        '/'.join([base(reads_1), "hisat2", "summary.txt"]),
                        "--rg-id", "GRPundef", "-p", threads, "--reorder"])
    sam = '/'.join([base(reads_1), "hisat2/alignment.sam"])
    bam = '/'.join([base(reads_1), "hisat2/alignment.bam"])
    command = ' '.join(["hisat2", "-x", ref, options, "-1", reads_1, "-2",
                        reads_2, "-S", sam])
    if not os.path.exists(os.path.dirname(sam)):
        os.makedirs(os.path.dirname(sam))
    print(command)
    sub.call(command, shell=True)
    print("samtools", "view", "-b", "-o", bam, "-@", threads, sam)
    sub.call(["samtools", "view", "-b", "-o", bam, "-@", threads, sam])
    if os.path.exists(sam) and os.path.exists(bam):
        os.remove(sam)
        # Quality checking
        sub.call([PARAMETERS.OPTIONS["SHELL"], "BamQC.sh", bam,
                  PARAMETERS.OPTIONS["BED"]])
        # List of bam files
        return bam
    else:
        print("No alignment sam/bam file.")
        sys.exit(1)


def kallisto_quantify(reads_1, reads_2, strand, threads):
    # Choose species
    ref = PARAMETERS.OPTIONS["KALLISTO_REF"]
    # Strand-specific option
    sp = {"none": "", "forward": "--fr-stranded", "reverse": "--rf-stranded"}
    sp = sp[strand]
    # Path
    options = ' '.join(["-i", ref, "-o",
                        '/'.join([base(reads_1), "kallisto"]), "--bias",
                        "-b 1000", "--seed=1984", "--fusion", sp, "-t",
                        threads])
    if not os.path.exists('/'.join([base(reads_1), "kallisto"])):
        os.makedirs('/'.join([base(reads_1), "kallisto"]))
    command = ' '.join(["kallisto quant", options, reads_1, reads_2])
    print(command)
    sub.call(command, shell=True)


def quantify(reads_1, reads_2, strand, mode, threads):
    if mode == 1:
        star_quantify(reads_1, reads_2, strand, threads)
    elif mode == 2:
        bam = hisat2_quantify(reads_1, reads_2, strand, threads)
        return bam
    elif mode == 3:
        kallisto_quantify(reads_1, reads_2, strand, threads)


def helper_quant(args):
    bam = quantify(args[0], args[1], args[2], args[3], args[4], )
    return bam


def setting():
    parser = \
        ap.ArgumentParser(description=
                          "Automatic RNA-seq paired-end quantification.")
    parser.add_argument("-F", "--list_R1", required=True,
                        help="List file of reads 1 fastq.")
    parser.add_argument("-R", "--list_R2",
                        help="List file of reads 2 fastq. [Optional]")
    parser.add_argument("-S", "--strand", default="none",
                        help="Strandedness. <none|forward|reverse> "
                             "[Default: none].\n"
                             "For Illumina TruSeq Stranded protocols, "
                             "use 'reverse'.")
    parser.add_argument("-M", "--mode", default=3, type=int,
                        help="Different tool chain. "
                             "<1: STAR-RSEM|2: HISAT2-feature|3: kallisto> "
                             "[Default: 3]")
    parser.add_argument("-T", "--threads", default=2, type=int,
                        help="Number of threads. [Default: 2]")
    parser.add_argument("-P", "--parameter_file", type=str,
                        default="./config_parameters.txt",
                        help="Some constant parameters. See example. "
                             "[Default: ./config_parameters.txt]")
    parser.add_argument("-v", "--version", action="version",
                        version="2019-08.", help="Show version")
    a = parser.parse_args()

    if a.strand not in ["none", "forward", "reverse"]:
        print("Please select corrected library.")
        sys.exit(1)
    if a.mode not in range(1, 4):
        print("Please select corrected tool chain.")
        sys.exit(1)
    if a.threads < 1:
        print("Wrong number of threads!")
        sys.exit(1)
    if a.parameter_file is None:
        print("Wrong parameter file!")
        sys.exit(1)

    return a


def main():
    load = setting()
    if load.list_R2 is None:
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
                if opt[0] != "PARALLEL_JOB":
                    PARAMETERS.OPTIONS[opt[0]] = opt[1]
                else:
                    THREAD.GROUP = int(opt[1])

    # Parallel running
    pool = mul.Pool(THREAD.GROUP)
    print("Parallel running", THREAD.GROUP, "jobs ...")
    # chunksize = 1 may lead to lower efficiency, but ordered results
    bam = pool.map(helper_quant,
                   it.zip_longest(f_list, r_list,
                                  it.repeat(load.strand, len(f_list)),
                                  it.repeat(load.mode, len(f_list)),
                                  it.repeat(str(load.threads), len(f_list))),
                   chunksize=1)

    # featureCounts
    if load.mode == 2 and bam is not None:
        PARAMETERS.BAM_LIST = bam
        d = os.path.dirname(f_list[0])
        if d == '':
            d = os.getcwd()
        do_featureCounts(PARAMETERS.BAM_LIST, str(load.threads), load.strand,
                         '/'.join([d, "featureCount/counts.txt"]))


if __name__ == "__main__":
    main()
