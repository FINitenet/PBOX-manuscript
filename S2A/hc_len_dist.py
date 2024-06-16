#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @FileName  :hc_len_dist.py
# @Time      :2023/09/27 14:42:07
# @Author    :Yuchen@rlab

import os   
import pysam
import pandas as pd
import glob
from Bio.Seq import Seq
from multiprocessing import Process


def seq_len_5bias_counter(input_file,filename):
    seq_len_5bias = {}
    with pysam.AlignmentFile(input_file, "rb") as input_f:
        for line in input_f:
            stats = line.get_tag('XS')
            if stats == 'Assigned':
                if line.is_reverse:
                    read_seq = line.query_sequence
                    read_seq = Seq(read_seq).reverse_complement()[0]
                else:
                    read_seq = line.query_sequence[0]

                read_len = line.query_alignment_length

                key = str(read_len) + "_" + str(read_seq)
                seq_len_5bias[key] = seq_len_5bias.get(key, 0) + 1
            else:
                continue
    df = pd.DataFrame(seq_len_5bias.items(), columns=['Key', 'Value'])
    df[['length', 'bias']] = df['Key'].str.split('_', expand=True)
    df.rename(columns={'Value': 'count'}, inplace=True)
    df = df[['length', 'bias', 'count']]
    df = df.sort_values(by='length')
    df['library'] = filename
    df.to_csv("/home/chenyc/Bioinformatics/chenyc/Project/TRM-sRNA-seq/QC/2_srna_analysis/" + str(filename)+ "_" + str(rnatype)+"_seq_len_5bias.csv", index=False)

if __name__ == "__main__":

    rnatype = "lncRNA"
    bamfiles = glob.glob("/home/chenyc/Bioinformatics/chenyc/Project/TRM-sRNA-seq/QC/2_srna_analysis/Annotation-level/*" + rnatype + ".bam")
    print(bamfiles)
    progress = []
    for bamfile in bamfiles:
        filename = os.path.basename(bamfile).split('.')[0]
        print(filename)
        p = Process(target=seq_len_5bias_counter, args=(bamfile,filename))
        p.start()
        progress.append(p)
    for p in progress:
        p.join()