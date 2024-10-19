#!/usr/bin/env python3
# script: filterGenomes.py
# author: Daniel Desiro'


import sys
import re
import os
from numpy import mean,std

infasta = sys.argv[1]
outd = sys.argv[2]




################################################################################
# main
################################################################################

def main():
    ############################################################################
    # read fasta
    seq_dict = readFasta(infasta)
    ############################################################################
    # extract proteins
    if not os.path.isdir(outd):
        os.mkdir(outd)
    extractProtein(seq_dict, outd)



################################################################################
# functions
################################################################################

def extractProtein(seq_dict, outdir):
    # extract segment
    # NCBI format: Accession, GenBank Title, Species, Length, Segment, Protein
    #seg_dict, unique_dict, dup_set = dict(), dict(), set()
    prot_dict = dict()
    for header,seq in seq_dict.items():
        acc,gbk,spc,lng,seg,pro = header
        if seg in ["10","11"] and ("NS5" in pro or "NSP5" in pro or "nonstructuralprotein5" in pro.lower().replace(" ","").replace("-","")):
            prot_list = prot_dict.get(pro,[])
            prot_list.append(header)
            prot_dict[pro] = prot_list
    with open(os.path.join(outdir, f"rotavirus_A_proteins.txt"), "w") as outprot:
        for pro in prot_dict.keys():
            outprot.write(f"{pro}\n")
    p_all, u_all, d_all, len_list = dict(), dict(), set(), list()
    for pro,prot_list in prot_dict.items():
        print_dict, unique_dict, dup_set = dict(), dict(), set()
        for header in prot_list:
            seq = seq_dict[header]
            h_name = "|".join(list(header))
            if seq not in dup_set:
                unique_dict[h_name] = seq
                dup_set.add(seq)
            print_dict[h_name] = seq
            if seq not in d_all:
                u_all[h_name] = seq
                d_all.add(seq)
                len_list.append(len(seq))
            p_all[h_name] = seq
        writeFasta(os.path.join(outdir, f"rotavirus_A_prot_{pro}.fa"), print_dict)
        writeFasta(os.path.join(outdir, f"rotavirus_A_prot_{pro}_unique.fa"), unique_dict)
    msig2 = mean(len_list)-2*std(len_list)
    psig2 = mean(len_list)+2*std(len_list)
    ua_all = {k:v for k,v in u_all.items()}
    ul_all = {k:v for k,v in u_all.items() if len(v) >= msig2}
    us_all = {k:v for k,v in u_all.items() if len(v) <= psig2}
    um_all = {k:v for k,v in u_all.items() if len(v) >= msig2 and len(v) <= psig2}
    writeFasta(os.path.join(outdir, f"rotavirus_A_NSP5_unique.fa"), ua_all)

def readFasta(infasta):
    # read fasta
    seq_dict = dict()
    with open(infasta, "r") as infa:
        RNA = ""
        for i,line in enumerate(infa):
            line = line.strip()
            if re.match(">",line):
                if RNA:
                    if ("partial" not in header and "truncated" not in header) and ("B" not in RNA and "Z" not in RNA and "X" not in RNA and "J" not in RNA):
                        seq_dict[tuple([h.strip() for h in header.replace("/","_").split("|")])] = RNA
                header, RNA = line, ""
            elif line: RNA += line
        if ("partial" not in header and "truncated" not in header) and ("B" not in RNA and "Z" not in RNA and "X" not in RNA and "J" not in RNA):
            seq_dict[tuple([h.strip() for h in header.replace("/","_").split("|")])] = RNA
    return seq_dict

def writeFasta(outfasta, seq_dict):
    # write fasta
    with open(outfasta, "w") as outfa:
        for header,seq in seq_dict.items():
            outfa.write(f"{header}\n")
            outfa.write(f"{seq}\n")


################################################################################
# call main
################################################################################

main()