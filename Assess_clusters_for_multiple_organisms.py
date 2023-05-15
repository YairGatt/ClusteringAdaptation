#!/usr/bin/env python

"""
Search for patients with same organism
"""

import sys, os
import argparse

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
		'-i', '--input',
		help='Input file with clusters crated by metabolic_clustering.py.')
	
	parser.add_argument(
		'-o', '--outfile', default="./temp.txt",
		help='Outfile where results will be written.')
	
	parser.add_argument(
		'-p', '--patients', nargs="+",
		help='All patients from all disks.')
	
	parser.add_argument(# customized description; put --help last
		'-h', '--help', action='help',
		help='Show this help message and exit.')
	
	settings = parser.parse_args(argv)
	
	return settings

def parse_patients(patients):
	"""
	Pase all patients and get a dict with patient:organism
	"""
	#initialize
	patients_dict = {}
	#iterate through patients
	for patient in patients:
		#get organism
		if patient.endswith("/"):
			organism = os.path.basename(os.path.dirname(os.path.dirname(patient)))
			patient_base = os.path.basename(os.path.dirname(patient))
		else:
			organism = os.path.basename(os.path.dirname(patient))
			patient_base = os.path.basename(patient)
		#add to dict
		patients_dict[patient_base] = organism
	#return
	return patients_dict

def parse_clusters(input_file, patients_dict):
	"""
	Parse clustering results and replace patient names with actual organisms. Also get pathways if available
	"""
	#initialize
	clusters_dict = {}
	pathways_dict = {}
	#open file
	with open(input_file) as fl:
		#iterate
		while True:
			line = fl.readline()
			if not line: break
			#if line represents a cluster
			if line.startswith("Cluster_"):
				#get cluster number
				cluster_name = line.split("\t")[0].split("_")[1]
				#initialize cluster list
				clusters_dict[method][cluster_name] = []
				#iterate through patients and add
				for patient in line.strip().split("\t")[1:]:
					if patient in patients_dict: clusters_dict[method][cluster_name] += [patients_dict[patient]]
			#if line represents pathways
			elif line.startswith("Pathways_"):
				#initialize pathways list
				pathways_dict[method][cluster_name] = []
				#iterate through pathwyas and add
				for pathway in line.strip().split("\t")[1:]:
					pathways_dict[method][cluster_name] += [pathway]
			#if line represents a new method
			else:
				method = line.strip()
				clusters_dict[method] = {}
				pathways_dict[method] = {}
	#return
	return clusters_dict, pathways_dict

def assess_multiplicity(clusters_dict):
	"""
	Search for clusters with multiple organisms
	"""
	#initialize
	multiples = {}
	#iterate
	for method in clusters_dict:
		#iterate
		for cluster in clusters_dict[method]:
			#assess number of organisms
			organisms = list(set(clusters_dict[method][cluster]))
			#see if multiples
			if len(organisms) > 1:
				#add method to dict
				if method not in multiples: multiples[method] = {}
				#add cluster to dict
				multiples[method][cluster] = organisms
	#return
	return multiples

def write_to_file(outfile, multiples, pathways_dict):
	"""
	Write clusters with multiple organisms to file
	"""
	#open
	with open(outfile, "w") as outfl:
		#add header
		header ="%s\t%s\t%s\t%s\n" % ("Clustering_Method", "Cluster", "Organisms", "Pathways")
		outfl.write(header)
		#iterate
		for method in multiples:
			print method
			#iterate
			for cluster in multiples[method]:
				#add organisms
				organisms = multiples[method][cluster][0]
				for organism in multiples[method][cluster][1:]:
					organisms += ",%s" % organism
				#add pathways if available
				if method in pathways_dict and cluster in pathways_dict[method] and len(pathways_dict[method][cluster]) >= 1:
					pathways = pathways_dict[method][cluster][0]
					for pathway in pathways_dict[method][cluster][1:]:
						pathways += ",%s" % pathway
				else:
					pathways = "-"
				#construct lines
				line = "%s\tCluster_%s\t%s\t%s\n" % (method, cluster, organisms, pathways)
				#write
				outfl.write(line)

def main(argv=None):
	#process command line
	settings = process_command_line(argv)
	#parse patients
	patients_dict = parse_patients(settings.patients)
	#parse file and replace patients with organisms, also get pathways if available
	clusters_dict, pathways_dict = parse_clusters(settings.input, patients_dict)
	#exit()
	#search for clusters with multiple organisms
	multiples = assess_multiplicity(clusters_dict)
	#exit()
	#write results if there are any to write
	if multiples: write_to_file(settings.outfile, multiples, pathways_dict)

if __name__ == "__main__":
		exit(main())