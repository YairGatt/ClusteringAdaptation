#!/usr/bin/env python

"""
This script translates kegg pathway format to format usable in my scripts.
It can be called from another script or from the command line.
The script can recieves the org.list file and converts it
Also needs to recieve the cds_from_genomic file for conversion
"""

import os,sys
from Bio import SeqIO
import gzip as gz
from kegg_to_NCBI import parse_cds
import warnings

def parse_kegg(kegg_file, cds_dict, organism_code):
    """
    Parse kegg to dict
    """
    #initialzie dict
    kegg_dict = {}
    converter = 0
    #parse file
    with open(kegg_file) as fl: lines = fl.readlines()
    #iterate through lines
    for line in lines:
        split_lines = line.strip().split()
        if not split_lines[1].startswith(organism_code): continue
        pathway = split_lines[0].split(":")[1].replace(organism_code,"")
        gene_name = split_lines[1].split(":")[1]#.split(".")[0]
        try:
            if pathway in kegg_dict: kegg_dict[pathway].append(cds_dict[gene_name])
            else: kegg_dict[pathway] = [cds_dict[gene_name]]
        except KeyError:
            try:
                if pathway in kegg_dict: kegg_dict[pathway].append(cds_dict[gene_name.split(".")[0]])
                else: kegg_dict[pathway] = [cds_dict[gene_name.split(".")[0]]]
            except KeyError: continue
    #return
    return kegg_dict

def main(args):
    #get files
    kegg_file = args[0]
    cds_file = args[1]
    rna_file = args[2]
    #define_organism_name
    organism_code = os.path.basename(kegg_file).split(".")[0]
    #get cds_from genomic and parse to dict
    cds_dict = parse_cds(cds_file, rna_file)
    #create dict
    kegg_dict = parse_kegg(kegg_file, cds_dict, organism_code)
    #return
    return kegg_dict, cds_dict

if __name__ == "__main__":
    exit(main(sys.argv[1:]))
