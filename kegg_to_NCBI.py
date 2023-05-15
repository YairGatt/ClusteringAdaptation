#!/usr/bin/env python

"""
This script translates kegg identifiers to ncbi names.
It can be called from another script or from the command line.
The script can Can recieve a single name, a list of name or a file with names to convert.
Also needs to recieve the KEGG and the cds_from_genomic file
"""

import os,sys
from Bio import SeqIO
import gzip as gz
import warnings

def parse_kegg(kegg_file):
    """
    Parse kegg to dict
    """
    #initialzie dict
    kegg_dict = {}
    #parse file
    with open(kegg_file) as fl: lines = fl.readlines()
    #iterate through lines
    for line in lines:
        split_lines = line.strip().split()
        gene_name = split_lines[0].split(":")[1]
        gene_id = split_lines[1].split(":")[1]
        kegg_dict[gene_name] = gene_id
    #return
    return kegg_dict

def parse_cds(cds_file, rna_file):
    """
    Parse kegg to dict
    """
    #initialzie dict
    cds_dict = {}
    #parse file
    with gz.open(cds_file) as fl: genes = list(SeqIO.parse(fl,"fasta"))
    with gz.open(rna_file) as fl: genes += list(SeqIO.parse(fl,"fasta"))
    #iterate through genes
    for gene in genes:
        description = gene.description.split()
        if "locus_tag" not in gene.description:
            #warnings.warn("No locus tag in CDS file %s" % gene.description)
            continue
        if "cds" in description[0]: protein_name = description[0].split("cds_")[1].rsplit("_",1)[0]#.split(".")[0]
        if "rna" in description[0]: protein_name = description[0].split("rna_")[1].rsplit("_",1)[0]#.split(".")[0]
        #gene_id = [i for i in description if "GeneID" in i][0].replace("[","").replace("]","").split(":")[1]
        locus_tag = [i for i in description if "locus_tag=" in i][0].replace("[","").replace("]","").split("=")[1]
        #cds_dict[gene_id] = protein_name
        cds_dict[locus_tag] = protein_name
    #return
    return cds_dict

def main(args):
    #get input
    input = args[0]
    genes = []
    #determine type
    print(input, type(input) == str)
    if os.path.isfile(input): #input is file
        with open(input) as fl:
            genes = [i.strip() for i in fl.readlines()]
    elif type(input) == str:
        genes = [input]
    elif type(input) == list:
        genes = input
    #get files
    #kegg_file = args[1]
    cds_file = args[1]
    rna_file = args[2]
    #create dict
    #kegg_dict = parse_kegg(kegg_file)
    #get cds_from genomic and parse to dict
    cds_dict = parse_cds(cds_file, rna_file)
    #convert
    #print genes
    #normal_names = [cds_dict[kegg_dict[i]] for i in genes]
    normal_names = [cds_dict[i] for i in genes]
    #print normal_names
    return normal_names

if __name__ == "__main__":
    exit(main(sys.argv[1:]))
