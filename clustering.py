#!/usr/bin/env python

"""
Clustering patients of a single species based on the metabolic pathways that were affected
"""

import sys, os
import argparse
import metabolic_clustering
import write_metabolic_pathways
import pathway_dict

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
		help='Lists of genes that underwent loss of function from each patient.')

	parser.add_argument(
		'-d', '--reference_directory',
		help='Reference directory with NCBI files.')
	
	parser.add_argument(
		'-D', '--dimensions', type=int, default=2,
		help='Dimensions for PCA plot, can be 2 (save file) or 3 (show interactivelly).')
	
	parser.add_argument(
		'-k', '--kegg_file',
		help='Kegg list of pathways and genes involved in pathway, format is org.list.')
	
	parser.add_argument(
		'-c', '--clustering_methods', default=[], nargs='+',
		help='List of clustering methods, can be kmeans, gmm, hierarchical or spectral.')
	
	parser.add_argument(
		'-p', '--pathway_descriptions', default= "/home/hosts/disk20/metabolic_pathways/map_title.tab",
		help='Kegg list of pathways and their descriptions.')
	
	parser.add_argument(# customized description; put --help last
		'-h', '--help', action='help',
		help='Show this help message and exit.')
	
	settings = parser.parse_args(argv)
	
	return settings

def read_metabolic_pathways(lost_genes_files, workdir, descriptions_file):
	"""
	Read all needed objects from files
	"""
	#initialize descriptions
	converter = pathway_dict.Converter(descriptions_file)
	#initialize
	clustering_patients_dict = {}
	clustering_genes_dict = {}
	#get patient workdir
	name = os.path.basename(os.path.dirname(os.path.dirname(workdir)))
	#get patients results
	for patient_lost_file in lost_genes_files:
		patient_dir = os.path.dirname(patient_lost_file)
		patient = os.path.basename(patient_dir)
		metabolics_dir = os.path.join(patient_dir,name)
		#initialize dicts
		clustering_patients_dict[patient] = {}
		clustering_genes_dict[patient] = {}
		#add from files
		with open(os.path.join(metabolics_dir,"pathways.txt")) as fl: content = fl.readlines()
		for line in content:
			split_line = line.strip().split()
			clustering_patients_dict[patient][converter.convert_pathways(split_line[0])] = int(split_line[1])
		with open(os.path.join(metabolics_dir,"pathways_full.txt")) as fl:
			content = fl.readlines()
		for line in content:
			split_line = line.strip().split()
			clustering_genes_dict[patient][converter.convert_pathways(split_line[0])] = float(split_line[1])
	#return
	return clustering_patients_dict, clustering_genes_dict

def main(argv=None):
	#process command line
	settings = process_command_line(argv)
	#count how many times each pathway was affected in two ways: number of patients and number of genes
	#run write_metabolic_pathways with skip
	arguements = argparse.Namespace(workdir=settings.workdir, lost_genes_files=settings.lost_genes_files, reference_directory=settings.reference_directory, kegg_file=settings.kegg_file, skip=True)
	write_metabolic_pathways.main(arguements)
	#read results
	clustering_patients_dict, clustering_genes_dict = read_metabolic_pathways(settings.lost_genes_files, settings.workdir, settings.pathway_descriptions)
	if not sum([1 for i in clustering_patients_dict if clustering_patients_dict[i]]): raise Exception("No metabolic pathways affected for organsim %s" % os.path.basename(os.path.dirname(os.path.dirname(settings.lost_genes_files[0]))))
	#create directory
	metabolic_dir = os.path.join(settings.workdir,"clustering")
	write_metabolic_pathways.mkdir(metabolic_dir)
	#clustering
	if settings.clustering_methods != []:
		if settings.dimensions <= 2.5: dimensions = 2
		else: dimensions = 3
		metabolic_clustering.main(clustering_patients_dict, clustering_genes_dict, metabolic_dir=metabolic_dir, D=dimensions, clustering_methods=settings.clustering_methods)

if __name__ == "__main__":
		exit(main())
