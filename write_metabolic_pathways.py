#!/usr/bin/env python

"""
Create patient metabolic pathways files, called from command line or python script
"""

import sys, os
import argparse
import errno
import create_pathways
from Bio import SeqIO
import gzip as gz
import itertools

def process_call_from_script(argv):
	command_line_settings = process_command_line(["-p", "/home/hosts/disk20/metabolic_pathways/map_title.tab"])
	
	for i in dir(command_line_settings):
		try:
			getattr(argv, i)
		except AttributeError:
			setattr(argv, i, getattr(command_line_settings, i))
	return argv

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
		'--skip', action='store_true', default=False,
		help='Get existing files if already written.')

	parser.add_argument(
		'-p', '--pathway_descriptions', default= "/home/hosts/disk20/metabolic_pathways/map_title.tab",
		help='Kegg list of pathways and their descriptions.')

	parser.add_argument(# customized description; put --help last
		'-h', '--help', action='help',
		help='Show this help message and exit.')

	settings = parser.parse_args(argv)

	return settings

def mkdir(path):
	#create directory and don't crash if it already exists
	try:
		os.mkdir(path)
	except OSError as exc:
		if exc.errno != errno.EEXIST:
			raise exc
			pass

def count_lost(lost_genes_files, pathways, metabolic_proteins, metabolic=False):
	"""
	Count how many time each pathway was affected in two ways:
	1. In how many patients
	2. How many genes
	"""
	#initialize dicts
	pathway_patients_dict = {}
	pathway_genes_dict = {}
	gene_loss_vector = []
	#for clustering
	clustering_patients_dict = {}
	clustering_genes_dict = {}
	#initialize all pathways
	for pathway in pathways:
		pathway_patients_dict[pathway] = 0
		pathway_genes_dict[pathway] = 0
	#iterate through lost files
	for patient_lost_file in lost_genes_files:
		#define patient name
		patient = os.path.basename(os.path.dirname(patient_lost_file))
		#initialize clustering dict patient entry
		if patient in clustering_patients_dict: raise Exception("Patient %s appears twice in files" % patient)
		clustering_patients_dict[patient] = {}
		clustering_genes_dict[patient] = {}
		#parse lost genes from file
		with open(patient_lost_file) as fl:
			patient_genes = []
			for i in fl.readlines():
				gene = i.strip().split()
				if not metabolic or gene[0] in metabolic_proteins: patient_genes.append(gene)
			#add to vector
			gene_loss_vector.append([float(i[2]) for i in patient_genes])
		#iterate through pathways and see what genes they lost
		for pathway in pathways:
			normalization_factor = len(pathways[pathway])
			clustering_genes_dict[patient][pathway] = 0
			clustering_patients_dict[patient][pathway] = 0
			#flag symbolizing gene in pathways affected already
			patient_pathway_flag = 0
			#iterate through genes
			for gene in patient_genes:
				#check if gene in pathway
				if gene[0] in pathways[pathway]: #gene[0].split(".")[0]
					#add gene damage
					pathway_genes_dict[pathway] += float(gene[2])
					clustering_genes_dict[patient][pathway] += float(gene[2]) / normalization_factor
					#add to patient if not counted yet
					if patient_pathway_flag == 0:
						pathway_patients_dict[pathway] += 1
						clustering_patients_dict[patient][pathway] += 1
						patient_pathway_flag = 1
	#return
	return pathway_patients_dict, pathway_genes_dict, gene_loss_vector, clustering_patients_dict, clustering_genes_dict

def write_results(workdir, lost_genes_files, clustering_patients_dict, clustering_genes_dict, pathway_patients_dict, pathway_genes_dict, gene_loss_vector, metabolic_gene_loss_vector):
	"""
	Write results for each patient to dir in patient
	Write results to results dir in project for each pathway
	"""
	#get patient workdir
	name = os.path.basename(os.path.dirname(os.path.dirname(workdir)))
	#patient results
	#make directory
	for patient_lost_file in lost_genes_files:
		patient_dir = os.path.dirname(patient_lost_file)
		patient = os.path.basename(patient_dir)
		metabolics_dir = os.path.join(patient_dir,name)
		mkdir(metabolics_dir)
		#write cluster_patients_dict entry for patient
		with open(os.path.join(metabolics_dir,"pathways.txt"), "w") as outfl:
			for pathway in clustering_patients_dict[patient]:
				line = pathway + "\t" + str(clustering_patients_dict[patient][pathway]) + "\n"
				outfl.write(line)
		with open(os.path.join(metabolics_dir,"pathways_full.txt"), "w") as outfl:
			for pathway in clustering_genes_dict[patient]:
				line = pathway + "\t" + str(clustering_genes_dict[patient][pathway]) + "\n"
				outfl.write(line)
	#pathway results
	#make directory for bacteria if it doesn't exist
	mkdir(workdir)
	with open(os.path.join(workdir,"gene_loss_vector.txt"), "w") as outfl:
		for i in gene_loss_vector:
			outfl.write(str(i) + "\n")
	with open(os.path.join(workdir,"metabolic_gene_loss_vector.txt"), "w") as outfl:
		for i in gene_loss_vector:
			outfl.write(str(i) + "\n")
	with open(os.path.join(workdir,"pathways_patients.txt"), "w") as outfl:
		for pathway in pathway_patients_dict:
			line = pathway + "\t" + str(pathway_patients_dict[pathway]) + "\n"
			outfl.write(line)
	with open(os.path.join(workdir,"pathways_genes.txt"), "w") as outfl:
		for pathway in pathway_genes_dict:
			line = pathway + "\t" + str(pathway_genes_dict[pathway]) + "\n"
			outfl.write(line)


def filter_skipped(workdir, lost_genes_files,skip):
	#remove skipped patients from lost genes_files if skip==True
	if not skip: return lost_genes_files
	#get patient workdir
	name = os.path.basename(os.path.dirname(os.path.dirname(workdir)))
	#initialize non-skipped
	lost_genes_files_nonskipped = []
	#iterate through files
	for patient_lost_file in lost_genes_files:
		patient_dir = os.path.dirname(patient_lost_file)
		metabolics_dir = os.path.join(patient_dir,name)
		#if no directory patient is not to be skipped
		if not os.path.isdir(metabolics_dir) or not os.path.isfile(os.path.join(metabolics_dir,"pathways.txt")):
			lost_genes_files_nonskipped.append(patient_lost_file)
	#return
	return lost_genes_files_nonskipped

def main(argv=None):
	#process command line
	if argv == None: settings = process_command_line(argv)
	else: settings = process_call_from_script(argv)
	#filter skipped
	write_lost_genes_files = filter_skipped(settings.workdir, settings.lost_genes_files, settings.skip)
	#skip if no patients are left
	if len(write_lost_genes_files) == 0 and os.path.isfile(os.path.join(settings.workdir,"metabolic_gene_loss_vector.txt")): return 0
	#otherwise calculate how many times each pathway is affected for all patients and write files for the ones that are not to be skipped
	#list files in reference directory
	onlyfiles = [os.path.join(settings.reference_directory, f) for f in os.listdir(settings.reference_directory) if os.path.isfile(os.path.join(settings.reference_directory, f))]
	#get ncbi files
	cds_file = [i for i in onlyfiles if i.endswith("_cds_from_genomic.fna.gz")][0]
	rna_file = [i for i in onlyfiles if i.endswith("_rna_from_genomic.fna.gz")][0]
	protein_file = [i for i in onlyfiles if i.endswith("_protein.faa.gz")][0]
	#get pathways
	pathways, cds_dict = create_pathways.main([settings.kegg_file, cds_file, rna_file])
	#get list of metabolic proteins
	pathway_proteins = pathways.values()
	metabolic_proteins = list(set(itertools.chain.from_iterable(pathway_proteins)))
	#count how many times each pathway was affected in two ways: number of patients and number of genes
	metabolic_pathway_patients_dict, metabolic_pathway_genes_dict, metabolic_gene_loss_vector, metabolic_clustering_patients_dict, metabolic_clustering_genes_dict = count_lost(settings.lost_genes_files, pathways, metabolic_proteins, metabolic=True)
	pathway_patients_dict, pathway_genes_dict, gene_loss_vector, clustering_patients_dict, clustering_genes_dict = count_lost(settings.lost_genes_files, pathways, metabolic_proteins, metabolic=False)
	if metabolic_pathway_genes_dict != pathway_genes_dict or metabolic_clustering_genes_dict != clustering_genes_dict or metabolic_clustering_patients_dict != clustering_patients_dict or metabolic_pathway_patients_dict != pathway_patients_dict: raise Exception("Something is weird in the behaviour of metaboilc mode. Results differ from baseline.")
	#write results
	write_results(settings.workdir, write_lost_genes_files, clustering_patients_dict, clustering_genes_dict, pathway_patients_dict, pathway_genes_dict, gene_loss_vector, metabolic_gene_loss_vector)

if __name__ == "__main__":
		exit(main())
