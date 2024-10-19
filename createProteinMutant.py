#!/usr/bin/env python3
# script: createProteinMutant.py
# author: Daniel Desiro'
# dependencies: numpy, RNAfold, VARNAv3-93.jar, inkscape
script_usage="""
usage
    createProteinMutant.py -fsa <in fasta> -pfx <out_prefix> [options]

version
    createProteinMutant.py 0.0.1 (alpha)

dependencies
    x

description
    Validates significance of flexible regions in a sequence.

--prefix,-pfx
    output directory and prefix for result files

--fasta,-fsa
    input fasta file, first entry reference RNA sequence, second entry reference 
    protein sequence, third entry target protein sequence; reference and target 
    protein sequences should be aligned

--reverseComplement,-rvc
    creates reverse complement of each strain if set (default: False)

--targetFrame,-frm
    set the target reading frame (default: 0) (choices: 0,1,2)

--overwrite,-ovr
    overwrite data with folder named in prefix (default: False)

--customModel,-cml
    input custom nucleotide substitution matrix (default:
    AA:0.00,AC:0.34,AG:0.35,AT:0.31,
    CA:0.34,CC:0.00,CG:0.31,CT:0.35,
    GA:0.35,GC:0.31,GG:0.00,GT:0.34,
    TA:0.31,TC:0.35,TG:0.34,TT:0.00,
    -A:0.23,-C:0.26,-G:0.27,-T:0.24)

reference
    Reference.
"""

import argparse as ap
import sys
import os
import re
import time
import pickle
from operator import attrgetter, itemgetter
from itertools import product
from math import ceil, floor, log, isnan
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt


try:
    import multiprocessingPickle4
    import multiprocessing as mp
    ctxx = mp.get_context()
    ctxx.reducer = multiprocessingPickle4.MultiprocessingPickle4()
except:
    import multiprocessing as mp



################################################################################
## main
################################################################################

def main(opt):
    time_s = getTime()
    ############################################################################
    ## create output folder
    print(f"Status: Create output directory ...")
    opt = makeDir(**opt)
    time_s = getTime(time_s, f"Create output directory")
    ########################################################################
    ## read fasta file
    print(f"Status: Read fasta file ...")
    RNA_header, RNA, refAA_header, refAA, tarAA_header, tarAA = readFasta(**opt)
    time_s = getTime(time_s, f"Read fasta file")
    ########################################################################
    ## create mutant
    print(f"Status: Create mutant ...")
    createMutant(RNA_header, RNA, refAA_header, refAA, tarAA_header, tarAA, **opt)
    time_s = getTime(time_s, f"Create mutant")
    ############################################################################
    return opt["var_pfx"]




################################################################################
## functions
################################################################################

class aaSeq(object):
    ## RNA object
    def __init__(self, seq, frame, aa_start, aa_end):
        self.seq = seq
        self.frame = frame
        self.aa_start = aa_start
        self.aa_end = aa_end
        self.start = aa_start*3+frame
        self.end = aa_end*3+frame
        self.len = len(seq)

def readFasta(**opt):
    ## read fasta file
    sequence, first = "", 0
    with open(opt["var_fsa"], "r") as infa:
        for line in infa:
            line = line.strip()
            if re.match(r">", line):
                if first == 0:
                    RNA_header = line
                    first += 1
                elif first == 1:
                    refAA_header = line
                    RNA = revComp(sequence, **opt)
                    sequence = ""
                    first += 1
                elif first == 2:
                    tarAA_header = line
                    refAA = sequence.upper()
                    sequence = ""
                    first += 1
                else:
                    break
            else:
                sequence += line
        tarAA = sequence.upper()
    return RNA_header, RNA, refAA_header, refAA, tarAA_header, tarAA

def makeDir(**opt):
    ## create directory
    dir_name, dir_base = opt["var_pfx"], opt["var_pfx"]
    if not opt["var_ovr"]:
        i = 1
        while os.path.isdir(dir_name):
            dir_name = f"{dir_base}_{i}"
            i += 1
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)
    opt["var_pfx"] = dir_name
    return opt

def revComp(RNA, **opt):
    ## complement dictionary, or transform DNA to RNA
    RNA = RNA.upper()
    D2Rc = {"A":"T","U":"A","T":"A","C":"G","G":"C","R":"Y","Y":"R","M":"K",\
            "K":"M","S":"W","W":"S","B":"V","V":"B","D":"H","H":"D","N":"N"}
    if opt["var_rvc"]: RNA = "".join(D2Rc[i] for i in RNA[::-1])
    else:              RNA = RNA.replace("U","T")
    return RNA

def makeCodons():
    ## make mutations in IAV strains
    # protein coding (+): 5'-3' AUG 
    # complement     (-): 3'-5' UAC
    # structure      (-): 5'-3' CAU
    co_dict = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L", "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
               "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*", "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
               "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L", "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
               "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q", "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
               "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M", "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
               "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K", "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
               "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V", "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
               "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E", "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    rc_dict = {} # reverse complement codons
    aa_dict = {} # amino acid dictionary
    ra_dict = {} # reverse complement aa dict
    opt = {"var_rvc":True}
    for codon,aa in sorted(co_dict.items()):
        comp_codon = revComp(codon, **opt)
        rc_dict[comp_codon] = aa
        ra_list, aa_list = [ra_dict.get(aa,list()), aa_dict.get(aa,list())]
        ra_list.append(comp_codon)
        aa_list.append(codon)
        ra_dict[aa] = ra_list
        aa_dict[aa] = aa_list
    return co_dict, aa_dict, rc_dict, ra_dict

def aaSequence(RNA, frame, codon_table):
    ## get amino acid sequence
    slen = len(RNA)
    rlen = int((slen-frame)/3)*3
    tseq = RNA[frame:frame+rlen]
    codons = [tseq[i:i+3] for i in range(0, len(tseq), 3)]
    aa_seq = "".join(codon_table[c] for c in codons)
    return aa_seq

def createMutant(RNA_header, RNA, refAA_header, refAA, tarAA_header, tarAA, **opt):
    # create mutant RNA sequence
    co_dict, aa_dict, rc_dict, ra_dict = makeCodons()
    if opt["var_rvc"]:
        c_dict, a_dict = rc_dict, ra_dict
    else:
        c_dict, a_dict = co_dict, aa_dict
    frame = opt["var_frm"]
    refAA_complete = aaSequence(RNA, frame, c_dict)
    refAA_search = refAA.replace("-", "")
    prot = re.compile(fr"{refAA_search}")
    maa = aaSeq("",0,0,0)
    refAA_extract = [aaSeq(match.group(0), frame, match.start(0), match.end(0)) for match in prot.finditer(refAA_complete)][0]
    refAA_ext = "-"*refAA_extract.aa_start + refAA + "-"*(len(refAA_complete)-refAA_extract.aa_end)
    tarAA_ext = "-"*refAA_extract.aa_start + tarAA + "-"*(len(refAA_complete)-refAA_extract.aa_end)
    RNAmut_list = list()
    r, g = "", 0
    for rA,tA in zip(refAA_ext, tarAA_ext):
        if rA == tA:
            r += RNA[(g)*3+frame:(g+1)*3+frame]
            g += 1
        elif rA == "-":
            if r:
                RNAmut_list.append([r])
                r = ""
            muts = a_dict[tA][:]
            muts.insert(0,"---")
            RNAmut_list.append(muts)
        elif tA == "-":
            g += 1
            continue
        else:
            if r:
                RNAmut_list.append([r])
                r = ""
            muts = a_dict[tA][:]
            muts.insert(0,RNA[(g)*3+frame:(g+1)*3+frame].lower())
            RNAmut_list.append(muts)
            g += 1
    RNAmut_list.append([r])
    if frame:
        RNAmut_list.insert(0,[RNA[0:frame]])
    if len(refAA_complete)*3+frame < len(RNA):
        RNAmut_list.append([RNA[len(refAA_complete)*3+frame:]])
    RNAmut = ""
    model = {"AA":-1.00, "AC":0.34, "AG":0.35, "AT":0.31,
             "CA":0.34, "CC":-1.00, "CG":0.31, "CT":0.35,
             "GA":0.35, "GC":0.31, "GG":-1.00, "GT":0.34,
             "TA":0.31, "TC":0.35, "TG":0.34, "TT":-1.00,
             "-A":0.23, "-C":0.26, "-G":0.27, "-T":0.24}
    # extract substitution matrix
    model =  {m.split(':')[0]:float(m.split(':')[1]) for m in opt["var_cml"].replace(" ","").split(",")}
    print(model)
    for i,RNApart in enumerate(RNAmut_list):
        if len(RNApart) == 1:
            RNAmut += RNApart[0]
        else:
            org = RNApart[0]
            mut_list = list()
            for rna in RNApart[1:]:
                c = 0
                for o,r in zip(org,rna):
                    c += model[f"{o.upper()}{r.upper()}"]
                mut_list.append((c,org,rna))
            mut_list = sorted(mut_list, key=itemgetter(0,2))
            for i,j in zip(mut_list[0][1],mut_list[0][2]):
                if i.upper() == j:
                    RNAmut += j
                else:      
                    RNAmut += j.lower()
            print(mut_list)
    out_name = os.path.basename(os.path.abspath(opt["var_pfx"]))
    out_write = os.path.join(opt["var_pfx"], f"{out_name}_mutant.fa")
    with open(out_write, "w") as outfile:
        outfile.write(f"{RNA_header} WT\n")
        outfile.write(f"{RNA}\n")
        outfile.write(f"{RNA_header} mutant\n")
        outfile.write(f"{RNAmut}\n")

def getTime(time_s=0, name=""):
    if time_s and name:
        time_e = time.time()-time_s
        time_e = time.strftime("%H:%M:%S", time.gmtime(time_e))
        time_c = time.strftime('%x %X')
        print(f"Status: {name} finished at {time_c} in {time_e}")
    return time.time()

def extractFrames(data_dict, var_maa, var_cls, var_p):
    ## extract aa sequences for all 3 reading frames
    # --minAA,-maa
    # set minimum aa length to be considered a protein (default: 70)
    co_dict, aa_dict, rc_dict, ra_dict = makeCodons()
    if var_cls == "ssRNA-":
        c_dict, a_dict = rc_dict, ra_dict
    else:
        c_dict, a_dict = co_dict, aa_dict
    # do for all vRNAs
    prot = re.compile(r"M[A-Z]+?\*")
    i = 0
    with open(var_p, "w") as out:
        for name,vRNA in sorted(data_dict.items()):
            aa_lst = []
            maa = aaSeq("",0,0,0)
            # get maximum length aa sequences
            for frame in range(3):
                aa_seq = aaSequence(vRNA, frame, c_dict)
                aa_lst += [aaSeq(match.group(0),frame,match.start(0),match.end(0)) 
                    for match in prot.finditer(aa_seq) if len(match.group(0)) > var_maa]
            # get all sub aa sequences M to *
            new_lst = []
            for aa in aa_lst:
                match = prot.search(aa.seq[1:])
                while match:
                    if len(match.group(0)) > var_maa:
                        new_lst.append(aaSeq(match.group(0),aa.frame,aa.aa_start+1+match.start(0),aa.aa_end))
                    match = prot.search(match.group(0)[1:])
            ################################
            aa_lst += new_lst
            aa_lst = sorted(aa_lst, key=attrgetter("start"))
            # get maximum non overlapping aa sequences
            new_lst = []
            maa = aaSeq("",0,0,0)
            for aa in aa_lst:
                if aa.start < maa.end:
                    if aa.len > maa.len:
                        maa = aa
                    else:
                        continue
                else:
                    if maa.len:
                        new_lst.append(maa)
                    maa = aa
            if not new_lst: new_lst.append(maa)
            elif maa != new_lst[-1]: new_lst.append(maa)
            for aa in new_lst:
                out.write(f">{name} frame:{aa.frame} start:{aa.start} end:{aa.end} length:{aa.len}\n{aa.seq}\n")




################################################################################
## parser
################################################################################

if __name__ == "__main__":

    ############################################################################
    ## get time and save call
    sscript = sys.argv[0]
    start_time = time.time()
    current_time = time.strftime('%x %X')
    scall = " ".join(sys.argv[1:])
    with open(f"{sscript}.log", "a") as calllog:
        calllog.write(f"Start : {current_time}\n")
        calllog.write(f"Script: {sscript}\n")
        calllog.write(f"Call  : {scall}\n")
    print(f"Call: {scall}")
    print(f"Status: Started at {current_time}")
    ############################################################################
    ## transform string into int, float, bool if possible
    def trans(s):
        if isinstance(s, str):
            try: return int(s)
            except ValueError:
                try: return float(s)
                except ValueError:
                    if s in ["True", "False"]: return s == "True"
                    else: return s
        else: return s
    ############################################################################
    ## save documentation
    rx_text = re.compile(r"\n^(.+?)\n((?:.+\n)+)",re.MULTILINE)
    rx_oneline = re.compile(r"\n+")
    rx_options = re.compile(r"\((.+?)\:(.+?)\)")
    help_dict, type_dict, text_dict, mand_list = {}, {}, {}, []
    for match in rx_text.finditer(script_usage):
        argument = match.groups()[0].strip()
        text = " ".join(rx_oneline.sub("",match.groups()[1].strip()).split())
        argopts = {"action":"store", "help":None, "default":None, "choices":None}
        for option in rx_options.finditer(text):
            key = option.group(1).strip()
            var = option.group(2).strip()
            if var == "False": argopts["action"] = "store_true"
            if var == "True": argopts["action"] = "store_false"
            if key == "choices": var = [vs.strip() for vs in var.split(",")]
            if key == "default": var = trans(var)
            argopts[key] = var
        if argopts["default"]: add_default = f" (default: {str(argopts['default'])})"
        else: add_default = ""
        argopts["help"] = rx_options.sub("",text).strip()+add_default
        argnames = argument.split(",")
        if len(argnames) > 1:
            if argopts["default"] == None:
                mand_list.append(f"var_{argnames[1][1:]}")
            type_dict[f"var_{argnames[1][1:]}"] = argopts["default"]
            argopts["argshort"] = argnames[1]
            help_dict[argnames[0]] = argopts
        else:
            text_dict[argnames[0]] = argopts["help"]
    ############################################################################
    ## get arguments
    if text_dict["dependencies"]:
        desc = f"{text_dict['description']} (dependencies: {text_dict['dependencies']})"
    else:
        desc = text_dict['description']
    p = ap.ArgumentParser(prog=sscript, prefix_chars="-", usage=text_dict["usage"],
                          description=desc, epilog=text_dict["reference"])
    p.add_argument("-v", "--version", action="version", version=text_dict["version"])
    for argname,argopts in help_dict.items():
        argshort = argopts["argshort"]
        if argopts["choices"]:
            p.add_argument(argshort, argname,            dest=f"var_{argshort[1:]}",\
                           action=argopts["action"],     help=argopts["help"],\
                           default=argopts["default"],   choices=argopts["choices"])
        else:
            p.add_argument(argopts["argshort"], argname, dest=f"var_{argshort[1:]}",\
                           action=argopts["action"],     help=argopts["help"],\
                           default=argopts["default"])
    p._optionals.title = "arguments"
    opt = vars(p.parse_args())
    ############################################################################
    ## validate arguments
    if None in [opt[mand] for mand in mand_list]:
        print("Error: Mandatory arguments missing!")
        print(f"Usage: {text_dict['usage']} use -h or --help for more information.")
        sys.exit()
    for key,var in opt.items():
        if key not in mand_list:
            arg_req, arg_in = type_dict[key], trans(var)
            if type(arg_req) == type(arg_in):
                opt[key] = arg_in
            else:
                print(f"Error: Argument {key} is not of type {type(arg_req)}!")
                sys.exit()
    ############################################################################
    ## call main function
    try:
        #saved = main(**opt)
        saved = main(opt)
    except KeyboardInterrupt:
        print("Error: Interrupted by user!")
        sys.exit()
    except SystemExit:
        print("Error: System exit!")
        sys.exit()
    except Exception:
        print("Error: Script exception!")
        traceback.print_exc(file=sys.stderr)
        sys.exit()
    ############################################################################
    ## finish
    started_time = current_time
    elapsed_time = time.time()-start_time
    elapsed_time = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    current_time = time.strftime('%x %X')
    if saved:
        with open(f"{sscript}.log", "a") as calllog,\
             open(os.path.join(saved,f"call.log"), "a") as dirlog:
            calllog.write(f"Save  : {os.path.abspath(saved)}\n")
            calllog.write(f"Finish: {current_time} in {elapsed_time}\n")
            ## dirlog
            dirlog.write(f"Start : {started_time}\n")
            dirlog.write(f"Script: {sscript}\n")
            dirlog.write(f"Call  : {scall}\n")
            dirlog.write(f"Save  : {os.path.abspath(saved)}\n")
            dirlog.write(f"Finish: {current_time} in {elapsed_time}\n")
    print(f"Status: Saved at {saved}")
    print(f"Status: Finished at {current_time} in {elapsed_time}")
    sys.exit(0)
