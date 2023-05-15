#!/usr/bin/env python

"""
This script translates kegg pathways to their descriptions
"""

import os,sys
from Bio import SeqIO
import gzip as gz

class Converter(object):
    """
    Class converting GO categories to their name using the go.obo file from GO
    """
    def __init__(self,descriptions_file):
        """
        Initialize object
        """
        self.descriptions_dict = parse_descriptions(descriptions_file)

    def convert_pathways(self, pathway_list):
        """
        Get list of pathways or a single pathway and convert to names
        """
        if type(pathway_list) == list:
            return [self.descriptions_dict[i] for i in pathway_list]
        elif type(pathway_list) == str:
            if os.path.isfile(pathway_list): #input is file
                with open(pathways) as fl:
                    lines = fl.readlines()
                return [self.descriptions_dict[i.strip()] for i in lines]
            else:
                return self.descriptions_dict[pathway_list]
        else:
            raise Exception("%s is not of type list or str and cannot be converted" % go_list)

def parse_descriptions(descriptions_file):
    """
    Parse kegg to dict
    """
    #initialzie dict
    descriptions_dict = {}
    #parse file
    with open(descriptions_file) as fl:
        lines = fl.readlines()
    #iterate through lines
    for line in lines:
        if line.startswith("#"): continue
        split_lines = line.strip().split(None,1)
        pathway_num = split_lines[0]
        description = split_lines[1]
        descriptions_dict[pathway_num] = description
    #return
    return descriptions_dict

def main(args):
    #get input
    input = args[0]
    pathways = []
    #determine type
    if type(input) == str:
        if os.path.isfile(input): #input is file
            with open(input) as fl:
                pathways = [i.strip() for i in fl.readlines()]
        else:
            pathways = [input]
    elif type(input) == list:
        pathways = input
    #get files
    descriptions_file = args[1]
    #get cds_from genomic and parse to dict
    descriptions_dict = parse_descriptions(descriptions_file)
    #convert
    #print genes
    #normal_names = [cds_dict[kegg_dict[i]] for i in genes]
    normal_names = [descriptions_dict[i] for i in pathways]
    #print normal_names
    return normal_names

if __name__ == "__main__":
    exit(main(sys.argv[1:]))
