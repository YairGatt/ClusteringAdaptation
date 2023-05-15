#!/usr/bin/env python

"""
See which cluster have both enriched treatments/something else and pathways. Output file of clusters with enriched treatments, enriched tissues, enriched conditions and enriched references.
Also stringent files with clusters with treatments and not references
"""

import sys, os
import argparse
import errno
from Assess_clusters_for_multiple_organisms import parse_clusters

class Cluster(object):
	"""
	Class representing cluster with full qualities
	"""
	def __init__(self, method, number, pathways, properties):
		"""
		Initialize object
		"""
		self.method = method
		self.number = number
		self.properties = properties
		self.pathways = pathways

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
		'-p', '--cluster_properties', default="./Clusters_properties.txt",
		help='File with enrichment of clusters in different treatments, conditions, references and tissues. Output of Get_cluster_properties.py or Cluster_properties.sh')
	
	parser.add_argument(
		'-c', '--clusters', default="./Clusters.txt",
		help='Clusters for enrichment, output of metabolic_clustering.py or go_clustering.py.')
	
	parser.add_argument(
		'-w', '--workdir', default="./temp/",
		help='Workdir where results will be written.')
	
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

def parse_properties(properties):
	"""
	Parse file with properties of clusters
	"""
	#initialize
	properties_dict = {}
	#open file
	with open(properties) as fl:
		#iterate through lines
		for line in fl.readlines()[1:]:
			#split line
			split_line = line.strip().split()
			#get properties
			method = split_line[0]
			cluster_num = split_line[1].split("_")[1]
			category = split_line[2]
			specific = split_line[3]
			#add to dict
			if method not in properties_dict: properties_dict[method] = {}
			if cluster_num not in properties_dict[method]: properties_dict[method][cluster_num] = {}
			if category not in properties_dict[method][cluster_num]: properties_dict[method][cluster_num][category] = []
			properties_dict[method][cluster_num][category].append(specific)
	#return
	return properties_dict
		
def get_cluster_objects(pathways, properties):
	"""
	Initialize clusters with both pathways and properties as objects.
	Keep a list of all objects, object which do not have reference properties, objects that have only one property category
	"""
	#initialize
	cluster_objects = []
	stringent_cluster_objects = []
	super_stringent_cluster_objects = []
	#iterate
	for method in pathways:
		#iterate
		for cluster_name in pathways[method]:
			#check if pathways
			if pathways[method][cluster_name] and method in properties and cluster_name in properties[method]:
				cluster = Cluster(method, cluster_name, pathways[method][cluster_name], properties[method][cluster_name])
				cluster_objects.append(cluster)
				#stringent
				if "references" not in properties[method][cluster_name]: stringent_cluster_objects.append(cluster)
				if len(properties[method][cluster_name]) == 1: super_stringent_cluster_objects.append(cluster)
	#return
	return cluster_objects, stringent_cluster_objects, super_stringent_cluster_objects

def write_to_files(workdir, cluster_objects, stringent_cluster_objects, super_stringent_cluster_objects):
	"""
	Write results to separate files
	"""
	#iterate through lists
	for object_list, filebase in zip([cluster_objects, stringent_cluster_objects, super_stringent_cluster_objects], ["all", "no_ref", "single_category"]):
		#define cateogries
		if filebase != "no_ref": categories = ["treatments", "conditions", "tissues", "references"]
		else: categories = ["treatments", "conditions", "tissues"]
		#iterate through categories
		for category in categories:
			#define filename
			filename = os.path.join(workdir, "%s_%s.txt" % (filebase, category))
			#open file
			with open(filename, "w") as outfl:
				#write header
				header = "%s\t%s\t%s\t%s\t%s\n" % ("Method", "Cluster", "Category", "Specific", "Pathways")
				outfl.write(header)
				#iterate through clusters and see if they fit
				for cluster in object_list:
					if category in cluster.properties:
						for specific in cluster.properties[category]:
							line = "%s\tCluster_%s\t%s\t%s\t%s\n" % (cluster.method, cluster.number, category, specific, cluster.pathways)
							outfl.write(line)

def main(argv=None):
	#process command line
	settings = process_command_line(argv)
	#create workdir
	mkdir(settings.workdir)
	#parse clusters
	clusters, pathways = parse_clusters(settings.clusters, {})
	#parse properties
	properties = parse_properties(settings.cluster_properties)
	#get cluster objects
	cluster_objects, stringent_cluster_objects, super_stringent_cluster_objects = get_cluster_objects(pathways, properties)
	#write to files
	write_to_files(settings.workdir, cluster_objects, stringent_cluster_objects, super_stringent_cluster_objects)
	
if __name__ == "__main__":
		exit(main())