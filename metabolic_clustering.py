#!/usr/bin/env python

"""
Clustering patients by metabolic pathways that were affected during their 
"""

import sys, os
import numpy as np
import itertools
#Clustering utilities from sklearn
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_samples, silhouette_score
import scipy.stats as stats
#Plotting and clustering personal scripts
import spectral_clustering
from metabolic_plotting import plot_pca, plot_clusters, clustering_pathways_by_jaccard
#general utilities
from FDR import correct_multiple_hypotheses

def convert_to_lists(clustering_patients_dict):
	#convert patients dict to format of dict[patient] = list_of_pathways
	jaccard_dict = {}
	for patient in clustering_patients_dict:
		jaccard_dict[patient] = [i for i in clustering_patients_dict[patient] if clustering_patients_dict[patient][i] == 1]
		#remove empties
		if jaccard_dict[patient] == []: del jaccard_dict[patient]
	#return
	return jaccard_dict

def convert_to_binary_vectors(clustering_patients_dict):
	#convert patients dict to format of dict[patient] = vector_of_0s_and_1s_corresponding_to_pathways
	binary_dict = {}
	#get sorted pathways
	all_keys = []
	for patient in clustering_patients_dict:
		all_keys += [i for i in clustering_patients_dict[patient].keys() if clustering_patients_dict[patient][i] != 0]
	sorted_keys = sorted(set(all_keys))
	#created vector for each patient
	for patient in clustering_patients_dict:
		#initialize
		binary_dict[patient] = []
		#add all pathways
		for i in sorted_keys:
			if i in clustering_patients_dict[patient] and clustering_patients_dict[patient][i] != 0: binary_dict[patient].append(clustering_patients_dict[patient][i])
			else: binary_dict[patient].append(0)
		#another possiblity is just using the keys that everyone has! Maybe that's more fair
		#binary_dict[patient] = [clustering_patients_dict[patient][i] for i in sorted_keys]
		if 1 not in binary_dict[patient]: del binary_dict[patient]
	#return
	return binary_dict,sorted_keys

def convert_to_vectors(clustering_genes_dict):
	#convert patients dict to format of dict[patient] = vector_of_floats_corresponding_to_pathways
	vector_dict = {}
	#get sorted pathways
	all_keys = []
	for patient in clustering_genes_dict:
		all_keys += [i for i in clustering_genes_dict[patient].keys() if clustering_genes_dict[patient][i] != 0]
	sorted_keys = sorted(set(all_keys))
	#sorted_keys = sorted(clustering_genes_dict[patient].keys())
	#created vector for each patient
	for patient in clustering_genes_dict:
		#initialize
		vector_dict[patient] = []
		#add all pathways
		for i in sorted_keys:
			if i in clustering_genes_dict[patient] and clustering_genes_dict[patient][i] != 0: vector_dict[patient].append(clustering_genes_dict[patient][i])
			else: vector_dict[patient].append(0)
		#vector_dict[patient] = [clustering_genes_dict[patient][i] for i in sorted_keys]
		if len(set(vector_dict[patient])) == 1: del vector_dict[patient]
	#keys
	clear_keys = [i for i in sorted_keys]
	#return
	return vector_dict,clear_keys

def clustering_pathways_by_vectors(vector_array, max_num, method='kmeans', k=0):
	"""
	Cluster numeric or binary vector by either kmeans or gmm clustering method
	"""
	#limit to 100 clusters
	max_num = min(max_num, 50)
	#check different numbers of clusters to choose best
	range_n_clusters = list(range(2,max_num + 1))
	#initialize bic if that's what we're gonna use
	best = np.infty
	#initialize list oof bics or silhouettes of the different models
	measures = []
	#cluster
	for n_clusters in range_n_clusters:
		if method == "kmeans":
			clusterer = KMeans(n_clusters=n_clusters, random_state=10)
			cluster_labels = clusterer.fit_predict(vector_array)
			measure = silhouette_score(vector_array, cluster_labels)*-1
			#print("For n_clusters =", n_clusters, "The average silhouette_score is :", -measure)
		elif method == "gmm":
			clusterer = GaussianMixture(n_components=n_clusters, max_iter=1000, n_init=42, covariance_type='full').fit(vector_array)
			cluster_labels = clusterer.predict(vector_array)
			measure2 = clusterer.bic(vector_array)
			#measure = clusterer.aic(vector_array)
			measure = silhouette_score(vector_array, cluster_labels)*-1
			print(method, n_clusters, measure2, measure)
			#print("For n_clusters =", n_clusters, "The BIC score is :", measure)
		elif method == "spectral":
			#np.save("va.npy", vector_array)
			#exit()
			cluster_labels = SpectralClustering(n_clusters=n_clusters, n_init=42).fit_predict(vector_array)
			#measure = clusterer.bic(vector_array)
			#measure = clusterer.aic(vector_array)
			measure = silhouette_score(vector_array, cluster_labels)*-1
			#print("For n_clusters =", n_clusters, "The BIC score is :", measure)

		#add to list
		measures.append(measure)
		
		if measures[-1] < best:
			best = measures[-1]
			if method != "spectral": best_model = clusterer
			best_labels = cluster_labels

	#get n_cluster that are not much worse, up to 5% measure difference
	#intialzie
	ok_measures = []
	#iterate through measures of all n_clusters
	for n,measure in enumerate(measures):
		#check if they are good enough
		if abs(best - measure) < abs(best / 20):
			ok_measures.append((range_n_clusters[n], measure))
	#output to user
	for n_clusters in ok_measures:
		print("%s is acceptable: %s" % (n_clusters[0], n_clusters[1]))
	
	#cluster
	if k != 0:
		if method == "kmeans":
			best_model = KMeans(n_clusters=k, random_state=0).fit(vector_array)
			best_labels = best_model.predict(vector_array)
		elif method == "gmm":
			best_model = GaussianMixture(n_components=k, max_iter=1000, n_init=42, covariance_type='full').fit(vector_array)
			best_labels = best_model.predict(vector_array)
		elif method == "spectral":
			best_labels = SpectralClustering(n_clusters=k, n_init=42).fit_perdict(vector_array)
		
	
	#Predict probabilities		
	if method == "gmm":
		probs = best_model.predict_proba(vector_array)
		#print probs.round(2)
	else:
		probs = 0
	#return
	return best_labels, probs

def to_array(a_dict):
	#initialize
	array = 0
	#iterate
	for i in a_dict.keys():
		if type(array) == int: array = np.array(list(a_dict[i]))
		else: array = np.vstack((array, a_dict[i]))
	#return
	return array

def external_to_labels(labels, patient_list):
	"""
	Convert format of labels from {label:[patients]} dict to [label, label] numpy array
	"""
	#initialize list
	new_labels = []
	#iterate through patients
	for patient in patient_list:
		#iterate thorugh labels
		for label in labels:
			#check if patient has the label and if so append
			if patient in labels[label]: new_labels.append(label)
	#convert to array
	new_labels = np.array(new_labels)
	#return
	return new_labels
	

def clustering(metabolic_dir, array, a_dict, keys, label, patient_list, clustering_method="kmeans", D=3):
	"""
	Cluster array kmeans and plot PCA and heatmaps
	"""
	#initialize name for external labels
	name = "%s_%s" % (clustering_method, label)
	labels = None
	#mode
	if "binary" in label:
		mode = "binary %s" % clustering_method
		metric = "jaccard"
	else:
		mode = clustering_method
		metric = "euclidean"
	#label
	if clustering_method == "kmeans":
		#get label for filename
		clustering_label = "KMeans"
		#if no external labels provided, run clustering
		if type(labels) != np.ndarray:
			#cluster
			labels, probabilities = clustering_pathways_by_vectors(array, len(a_dict) - 1, clustering_method)
		#probabilities are not calculated by kmeans
		probabilities = 0
	elif clustering_method == "gmm":
		clustering_label = "GaussianMixture"
		#if no external labels provided, run clustering
		if type(labels) != np.ndarray: labels, probabilities = clustering_pathways_by_vectors(array, len(a_dict) - 1, clustering_method)
	elif clustering_method == "spectral":
		clustering_label = "Spectral"
		#if no external labels provided, run clustering
		#if type(labels) != np.ndarray: labels = spectral_clustering.spectral_clustering(array, len(a_dict) - 1, 2, True, metric=metric)
		if type(labels) != np.ndarray: labels, probabilities = clustering_pathways_by_vectors(array, len(a_dict) - 1, clustering_method)
		#probabilities are not calculated by spectral clustering
		probabilities = 0
	#pca
	outfile = os.path.join(metabolic_dir,"%s_%s_PCA.pdf" % (clustering_label, label))
	print("plotting")
	#plot_pca(outfile, a_dict.keys(), array, labels, probabilities, D, mode)
	#plot clusters
	outfile_clustering = os.path.join(metabolic_dir,"%s_%s_clustering.pdf" % (clustering_label, label))
	#plot_clusters(labels, array, keys, patient_list, outfile=outfile_clustering, method=mode, metric=metric)
	#return
	return labels

def get_cluster_pathways(cluster, clustering_genes_dict, maximum, THRESHOLD=0.5):
	"""
	Get pathways that appear above THRESHOLD of patients in cluster
	"""
	#initialize pathway:number of times it appears
	pathway_dict = {}
	#iterate through patients
	for patient in cluster:
		#iterate through patwhays
		for pathway in clustering_genes_dict[patient]:
			#if the effect is too small, continue
			if clustering_genes_dict[patient][pathway] < maximum / 20: continue
			#add to dict
			if pathway in pathway_dict: pathway_dict[pathway] += 1
			else: pathway_dict[pathway] = 1
	#normalize
	for pathway in pathway_dict:
		pathway_dict[pathway] = float(pathway_dict[pathway]) / len(cluster)
	#final pathways
	pathways = [pathway for pathway in pathway_dict if pathway_dict[pathway]> THRESHOLD]
	#return
	return pathways

def get_cluster_pathways_stat(cluster, clustering_genes_dict):
        """
        Get pathways that appear in statistically significant number of patients in cluster by fisher's exact test or hypergeom
        """
        #initialize pathway:number of times it appears
        cluster_pathway_dict = {}
        #get all relevant pathways
        for patient in cluster:
            for pathway in clustering_genes_dict[patient]:
                if pathway not in cluster_pathway_dict: cluster_pathway_dict[pathway] = [clustering_genes_dict[patient][pathway]]
                else: cluster_pathway_dict[pathway] += [clustering_genes_dict[patient][pathway]]
        #iterate through pathways and get final nums
        noncluster_pathway_dict = {}
        for pathway in cluster_pathway_dict:
            noncluster_pathway_dict[pathway] = []
            for patient in clustering_genes_dict:
                if patient not in cluster: noncluster_pathway_dict[pathway] += [clustering_genes_dict[patient][pathway]]
        #perform wilcoxon test
        pdict = {}
        for pathway in cluster_pathway_dict:
            #table is [cluster-pathway, cluster-not-pathway], [not cluster-pathway, not cluster-not pathway]
            #cluster_pathway = cluster_pathway_dict[pathway]
            #cluster_not_pathway = len(cluster) - cluster_pathway
            #not_cluster_pathway = noncluster_pathway_dict[pathway]
            #not_cluster_not_pathway = len(clustering_genes_dict) - len(cluster) - noncluster_pathway_dict[pathway]
            #oddsratio, pvalue = stats.fisher_exact([[cluster_pathway, cluster_not_pathway], [not_cluster_pathway, not_cluster_not_pathway]], alternative='greater')
            #pdict[pathway] = pvalue
            print(pathway, cluster_pathway_dict[pathway], noncluster_pathway_dict[pathway], len(set(cluster_pathway_dict[pathway])), cluster_pathway_dict[pathway][0])
            if len(set(cluster_pathway_dict[pathway])) == 1 and cluster_pathway_dict[pathway][0] == 0.0: continue
            pdict[pathway] = stats.mannwhitneyu(cluster_pathway_dict[pathway], noncluster_pathway_dict[pathway], alternative='greater')[1]
        print(pdict)
        #correct
        pdict = correct_multiple_hypotheses(pdict)
        #final pathways
        pathways = [pathway for pathway in pdict if pdict[pathway] <= 0.1]
        #return
        return pathways

def write_results(metabolic_dir, all_labels, patient_list, maximum, clustering_genes_dict, THRESHOLD=0.6):
	"""
	Write results to file
	"""
	with open(os.path.join(metabolic_dir,"Clusters.txt"),"a") as outfl:
		#iterate through all clustering methods
		for method in all_labels:
			#define clusters
			clusters = []
			cluster_pathways = []
			#which cluster are we looking at
			for cluster_name in set(all_labels[method]):
				print(cluster_name)
				#initializepatients in cluster
				cluster = []
				#iterate through all patients
				for n,label in enumerate(all_labels[method]):
					#add patients in cluster to cluster
					if label == cluster_name: cluster.append(patient_list[n])
				#add cluster to list of clusters
				clusters.append(cluster)
				#get cluster pathways
				#pathways = get_cluster_pathways(cluster, clustering_genes_dict, maximum, THRESHOLD)
				pathways = get_cluster_pathways_stat(cluster, clustering_genes_dict)
				cluster_pathways.append(pathways)
			#now we have the clusters for the label and we can write them to the file
			#write clustering method
			outfl.write(str(method) + "\n")
			#write each cluster
			for n, cluster in enumerate(clusters):
				line = "Cluster_%s" % str(n + 1) #initialize line
				#add all patients
				for patient in cluster:
					line += "\t%s" % patient
				#finalize line
				line += "\n"
				#write line
				outfl.write(line)
				#pathway line
				pathway_line = "Pathways_%s" % str(n + 1) #initialize line
				#add all pathways
				for pathway in cluster_pathways[n]:
					pathway_line += "\t%s" % pathway
				#finalize line
				pathway_line += "\n"
				#write line
				outfl.write(pathway_line)

def main(clustering_patients_dict, clustering_genes_dict, metabolic_dir, D=2, clustering_methods=['kmeans']):
	#remove excluded
	excluded = []
	for i in excluded:
		if i in clustering_patients_dict: del clustering_patients_dict[i]
		if i in clustering_genes_dict: del clustering_genes_dict[i]
	#run clsutering in different methods
	clustering_methods = [i.lower() for i in clustering_methods]
	#prepare data
	#binary 0,1 vector for pathways
	binary_dict, binary_keys = convert_to_binary_vectors(clustering_patients_dict)
	#numeric vector of floats for pathways
	vector_dict, vector_keys = convert_to_vectors(clustering_genes_dict)
	#define arrays
	vector_array = to_array(vector_dict)
	binary_array = to_array(binary_dict)
	patient_list = vector_dict.keys()
	#initialize labels
	all_labels = {}
	#cluster by jaccard index in hierarchical clustering
	if "hierarchical" in clustering_methods:
		outfile_jaccard = os.path.join(metabolic_dir,"Heirarchical_Jaccard.pdf")
		jaccard_dict = convert_to_lists(clustering_patients_dict)
		#clustering_pathways_by_jaccard(jaccard_dict, outfile_jaccard)
	#cluster all others
	for array, a_dict, keys, label in zip([vector_array, binary_array], [vector_dict, binary_dict], [vector_keys, binary_keys], ["numeric", "binary"]):
		#iterate through methods
		for clustering_method in [i for i in clustering_methods if i != "hierarchical"]:
			if label == "binary": continue
			#run
			all_labels["%s_%s" % (clustering_method, label)] = clustering(metabolic_dir, array, a_dict, keys, label, patient_list, clustering_method, D)
	#get maximum
	maximum = max(list(itertools.chain.from_iterable(vector_array.tolist())))
	#output results to file
	write_results(metabolic_dir, all_labels, patient_list, maximum, clustering_genes_dict)
	exit()
	with open(os.path.join(metabolic_dir,"Numeric_array.csv"),"w") as outfl:
		#add header
		line = ""
		for i in vector_keys:
			line += "\t%s" % i
		line += "\n"
		outfl.write(line)
		#add patients
		for patient in vector_dict:
			line = patient
			for i in vector_dict[patient]:
				line += "\t%s" % i
			line += "\n"
			outfl.write(line)
	
	with open(os.path.join(metabolic_dir,"Binary_array.csv"),"w") as outfl:
		#add header
		line = ""
		for i in binary_keys:
			line += "\t%s" % i
		line += "\n"
		outfl.write(line)
		#add patients
		for patient in binary_dict:
			line = patient
			for i in binary_dict[patient]:
				line += "\t%s" % i
			line += "\n"
			outfl.write(line)		

if __name__ == "__main__":
		exit(main())
