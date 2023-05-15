#!/usr/bin/env python

"""
Investigate metabolic enrichment in organism
"""

import sys, os
import argparse
import errno
import create_pathways
import write_metabolic_pathways
import pathway_dict
from Bio import SeqIO
import gzip as gz
import random
import numpy as np
from statsmodels.stats.multitest import multipletests
import itertools
import metabolic_clustering

def process_command_line(argv):
	"""
	Return an args list
	`argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
	"""
	if argv is None:
		argv = sys.argv[1:]
		
	# initialize the parser object:
	parser = argparse.ArgumentParser(description='Process input.', add_help=False)
	
	#define options here
	parser.add_argument(
		'-w', '--workdir', default="./temp/",
		help='Workdir where results will be written.')
	
	parser.add_argument(
		'-l', '--lost_genes_files', nargs='+',
		help='Lists of genes that were lost from each patient.')

	parser.add_argument(
		'-d', '--reference_directory',
		help='Reference directory with NCBI files.')
	
	parser.add_argument(
		'-k', '--kegg_file',
		help='Kegg list of pathways and genes involved in pathway, format is org.list.')
	
	parser.add_argument(
		'--locus_file', default=None,
		help='File converting gebank to refseq annotations to genbank locus tags, available for Acinetobacter baumannii.')
	
	parser.add_argument(
		'-a', '--alpha', type=float, default= 0.05,
		help='Alpha of statistical test.')

	parser.add_argument(
		'--metabolic', action='store_true', default=False,
		help='Metabolic genes only mode.')

	parser.add_argument(
		'-p', '--pathway_descriptions', default= "/home/hosts/disk20/metabolic_pathways/map_title.tab",
		help='Kegg list of pathways and their descriptions.')

	parser.add_argument(# customized description; put --help last
		'-h', '--help', action='help',
		help='Show this help message and exit.')
	
	settings = parser.parse_args(argv)
	
	return settings

def proteins_list(protein_file):
	#get names of proteins from protein file
	with gz.open(protein_file) as fl:
		genes_list = [i.id for i in list(SeqIO.parse(fl, "fasta"))] # i.id.split(".")[0]
	#return
	return genes_list

def simulate_pathways(gene_list, gene_loss_vector, pathways, gene_pathway_dict, number_of_simulations=100000):
	"""
	Simualte gene loss in similar numbers to patients and see how many times the pathways are affected
	"""
	background_patient_distribution = {}
	background_genes_distribution = {}
	#initialize dicts
	for pathway in pathways:
		background_patient_distribution[pathway] = [0 for i in xrange(number_of_simulations)]
		background_genes_distribution[pathway] = [0 for i in xrange(number_of_simulations)]
	#start simulations
	for i in xrange(number_of_simulations):
		#simulation totals
		simulation = []
		#choose lost genes at random
		for patient_vector in gene_loss_vector:
			patient_lost_genes = []
			for fraction in patient_vector:
				gene = random.choice(gene_list)
				patient_lost_genes.append((gene,fraction))
			#append to simulation totals
			simulation.append(patient_lost_genes)
		#count how many times pathways were affected in this simulation
		for patient_vector in simulation:
			#flag symbolizing gene in pathways affected already
			patient_pathway_flag_dict = {}
			#iterate through genes
			for gene in patient_vector:
				#check if gene in pathway
				if gene[0] not in gene_pathway_dict: continue
				for pathway in gene_pathway_dict[gene[0]]:
					#add gene damage
					background_genes_distribution[pathway][i] += float(gene[1])
					#add to patient if not counted yet
					if pathway not in patient_pathway_flag_dict:
						background_patient_distribution[pathway][i] += 1
						patient_pathway_flag_dict[pathway] = 1
	#return
	return background_patient_distribution, background_genes_distribution

def get_p_vals(pathway_patients_dict, pathway_genes_dict, background_patient_distribution, background_genes_distribution):
	"""
	Calculate a p-value for each pathway based on its distribution
	"""
	patient_p_vals = {}
	genes_p_vals = {}
	#iterate through pathways
	for pathway in pathway_patients_dict:
		#start with patients p-value
		#get experiment
		number = pathway_patients_dict[pathway]
		#get background
		affected_at_least = sum([1 for i in background_patient_distribution[pathway] if i >= float(number)])
		#calculate p_value
		patient_p_vals[pathway] = float(affected_at_least) / len(background_patient_distribution[pathway])
		#continue to genes p-value
		#get experiment
		number = pathway_genes_dict[pathway]
		#get background
		affected_at_least = sum([1 for i in background_genes_distribution[pathway] if i >= float(number)])
		#calculate p_value
		genes_p_vals[pathway] = float(affected_at_least) / len(background_genes_distribution[pathway])
	#return results
	return patient_p_vals, genes_p_vals

def correct_multiple_hypotheses(p_value_dict, alpha):
        #print p_value_dict
        #exit()
	#initialize adjusted p values dict
	padjust_dict = {}
	#run multiple hypothesis correction in python using the Benjamini-Hochberg method
	truths, padjust, correctedp1, correctedp2 = multipletests(p_value_dict.values(), alpha*2, method = 'fdr_bh')
	#define p_adjust dict
	for index, pathway in enumerate(p_value_dict.keys()):
		padjust_dict[pathway] = padjust[index]
		
	#return results
	return padjust_dict

def read_metabolic_pathways(lost_genes_files, workdir, metabolic=False):
	"""
	Read all needed objects from files
	"""
	#get patietns dict
	pathway_patients_dict_file = os.path.join(workdir,"pathways_patients.txt")
	pathway_patients_dict = {}
        #print pathway_patients_dict_file
	with open(pathway_patients_dict_file) as fl: content = fl.readlines()
        #print content
        #exit()
	for line in content:
		split_line = line.strip().split()
		pathway_patients_dict[split_line[0]] = int(split_line[1])
        #print pathway_patients_dict
        #exit()
	#get genes dict
	pathway_genes_dict_file = os.path.join(workdir,"pathways_genes.txt")
	pathway_genes_dict = {}
	with open(pathway_genes_dict_file ) as fl: content = fl.readlines()
	for line in content:
		split_line = line.strip().split()
		pathway_genes_dict[split_line[0]] = float(split_line[1])
	#get gene loss vector
	if metabolic: gene_loss_vector_file = os.path.join(workdir,"metabolic_gene_loss_vector.txt")
	else: gene_loss_vector_file = os.path.join(workdir,"gene_loss_vector.txt")
	gene_loss_vector = []
	with open(gene_loss_vector_file ) as fl: content = fl.readlines()
	for line in content:
		split_line = line.strip().replace("[","").replace("]","").split(", ")
		if not split_line[0]:
			gene_loss_vector.append([])
			continue
		gene_loss_vector.append([float(i) for i in split_line])
	
	#return
	return pathway_patients_dict, pathway_genes_dict, gene_loss_vector

def switch_dict(the_dict):
	"""
	Change dict in format x1:[y1,y2,y3] to format y1:[x1,x2,x3]
	"""
	#initialzie
	new_dict = {}
	#iterate through items
	for x in the_dict:
		#iterate through values
		for y in the_dict[x]:
			#add to new dict
			if y in new_dict: new_dict[y] += [x]
			else: new_dict[y] = [x]
	#return
	return new_dict

def main(argv=None):
	#process command line
	settings = process_command_line(argv)
	#list file in reference directory
	onlyfiles = [os.path.join(settings.reference_directory, f) for f in os.listdir(settings.reference_directory) if os.path.isfile(os.path.join(settings.reference_directory, f))]
	#get ncbi files
	cds_file = [i for i in onlyfiles if i.endswith("_cds_from_genomic.fna.gz")][0]
	rna_file = [i for i in onlyfiles if i.endswith("_rna_from_genomic.fna.gz")][0]
	protein_file = [i for i in onlyfiles if i.endswith("_protein.faa.gz")][0]
	#get gene list
	genes_list = proteins_list(protein_file)
	#get pathways
	if settings.locus_file: pathways, cds_dict = create_pathways.main([settings.kegg_file, cds_file, rna_file, settings.locus_file])
	else: pathways, cds_dict = create_pathways.main([settings.kegg_file, cds_file, rna_file])
	#get gene: pathway dict
	gene_pathway_dict = switch_dict(pathways)
	#get list of metabolic proteins
	pathway_proteins = pathways.values()
	metabolic_proteins = list(set(itertools.chain.from_iterable(pathway_proteins)))
	#count how many times each pathway was affected in two ways: number of patients and number of genes
	#run write_metabolic_pathways with skip
	arguements = argparse.Namespace(workdir=settings.workdir, lost_genes_files=settings.lost_genes_files, reference_directory=settings.reference_directory, kegg_file=settings.kegg_file, locus_file=settings.locus_file, skip=False)
	write_metabolic_pathways.main(arguements)
	#read results
	pathway_patients_dict, pathway_genes_dict, gene_loss_vector = read_metabolic_pathways(settings.lost_genes_files, settings.workdir, settings.metabolic)
	#simulated background distribution
	if not settings.metabolic: background_patient_distribution, background_genes_distribution = simulate_pathways(genes_list, gene_loss_vector, pathways, gene_pathway_dict)
	else: background_patient_distribution, background_genes_distribution = simulate_pathways(metabolic_proteins, gene_loss_vector, pathways, gene_pathway_dict)
	#get p-values for all pathways
        #print pathway_patients_dict
        #exit()
	patient_p_vals, genes_p_vals = get_p_vals(pathway_patients_dict, pathway_genes_dict, background_patient_distribution, background_genes_distribution)
        #print patient_p_vals
        #exit()
	#correct for multiple hypotheses
	corrected_patient_p_vals = correct_multiple_hypotheses(patient_p_vals, settings.alpha)
	corrected_genes_p_vals = correct_multiple_hypotheses(genes_p_vals, settings.alpha)
	#oputput significnat pathways
	print "patients"
	for pathway in corrected_patient_p_vals:
		if corrected_patient_p_vals[pathway] <= 2*settings.alpha:
			print pathway_dict.main([pathway, settings.pathway_descriptions]), corrected_patient_p_vals[pathway]
	print "genes"
	for pathway in corrected_genes_p_vals:
		if corrected_genes_p_vals[pathway] <= 2*settings.alpha:
			print pathway_dict.main([pathway, settings.pathway_descriptions]), corrected_genes_p_vals[pathway]

if __name__ == "__main__":
		exit(main())
