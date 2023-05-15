#!/usr/bin/env python

"""
This script recieves a list of clusters and Table S2, and assess if any cluster is enriched in any properties
"""

import sys, os
import argparse
import csv
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import itertools

class Patient(object):
    """
    This class is a row in the database with the relevent properties
    """
    def __init__(self, tissues, treatments, additional_conditions, reference):
        self.tissues = tissues
        self.additional_conditions = additional_conditions
        self.treatments = treatments
        self.reference = reference
        self.parameters = {"tissues" : tissues, "conditions" : additional_conditions, "treatments" : treatments, "references" : reference}

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
        '-o', '--outfile', default="./temp.txt",
        help='Outfile where results will be written.')
    
    parser.add_argument(
        '-db', '--database',
        help='Database table.')
    
    parser.add_argument(
        '-c', '--clusters',
        help='Clusters for enrichment, output of metabolic_clustering.py or go_clustering.py.')
    
    parser.add_argument(# customized description; put --help last
        '-h', '--help', action='help',
        help='Show this help message and exit.')
    
    settings = parser.parse_args(argv)
    #return
    return settings

def parse_clusters(clusters_file):
    """
    Parse file with clusters to dict of
    clustering_method : cluster : patients
    also keep all patients for background
    """
    #initialize
    clusters_dict = {}
    pathways_dict = {}
    patients = []
    #open file
    with open(clusters_file) as fl:
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
                clusters_dict[method][cluster_name] += line.strip().split("\t")[1:]
                patients += line.strip().split("\t")[1:]
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
    #remove patients appearing multiple times
    patients = list(set(patients))
    #return
    return clusters_dict, patients

def parse_database(database_file, patients):
    """
    Parse database to list of patients, each defined with relevant properites
    """
    #initialize
    database = {}
    #open file
    with open(database_file) as fl:
        reader = csv.DictReader(fl, delimiter="\t")
        #iterate through rows
        current = 0
        for row in reader:
            #get trial
            trial = row["Experiment"]
            if trial == current: patient_num += 1
            else: patient_num = 1
            #get patient
            patient = "trial_%s_patient_%s" % (trial, patient_num)
            #define trial as currrent
            current = trial
            #continue if patient not in patients
            if patient not in patients: continue
            #if patient in patients add properties to db
            tissues = [i.lower() for i in row["Tissue"].split("+") if i]
            if len(tissues) == 1: tissues = tissues[0].split(" and ")
            #get antibiotics
            treatments = [i.lower() for i in row["Treatment"].split("+") if i]
            #additioanl conditions
            additional_conditions = [i.lower() for i in row["Additional_Condition"].split("+") if i]
            if len(additional_conditions) == 1: additional_conditions = additional_conditions[0].split(", ")
            if len(additional_conditions) == 1: additional_conditions = additional_conditions[0].split(" and ")
            #reference
            reference = row["Reference"]
            #add to Database
            database[patient] = Patient(tissues, treatments, additional_conditions, reference)
    #return
    return database

def cluster_enrichment(clusters, database, patients):
    """
    Check enrichment of each cluster by all properties in db
    """
    #initialize
    enriched_clusters = {}
    #get totals
    general, tissue_dict, treatment_dict, condition_dict, reference_dict = get_totals(database)
    #iterate through each cluster
    for method in clusters:
        enriched_clusters[method] = {}
        for cluster in clusters[method]:
            enriched_clusters[method][cluster] = {}
            #initialize dicts
            tissues = {}
            conditions = {}
            treatments = {}
            references = {}
            #iterate
            for patient in clusters[method][cluster]:
                #add all categories
                for category, dict in zip(["tissues", "treatments", "conditions"], [tissues, treatments, conditions]):
                    for specific in database[patient].parameters[category]:
                        if specific in dict: dict[specific] += 1
                        else: dict[specific] = 1
                    #add reference
                if database[patient].parameters["references"] in references: references[database[patient].parameters["references"]] += 1
                else: references[database[patient].parameters["references"]] = 1
            #now that we have the dicts check for significance
            for specific_dict, general_dict, category in zip([tissues, treatments, conditions, references], [tissue_dict, treatment_dict, condition_dict, reference_dict], ["tissues", "treatments", "conditions", "references"]):
                #print cluster, method, category
                #initialize
                enriched_clusters[method][cluster][category] = {}
                #get num of patients in cluster with anything in category
                total = len([i for i in clusters[method][cluster] if database[i].parameters[category]])
                for possibility in specific_dict:
                    #how many in the cluster came from the tissue
                    x = specific_dict[possibility]
                    if x <= 1: continue
                    #how many in the cluster did not come from the tissue
                    y = total - x #replace with patients in cluster with cateogry
                    #how many outside the cluster came from the tissue
                    t = general_dict[possibility] - x
                    #how many outside the cluster did not come from the tissue
                    n = general[category] -(x + y + t)
                    if cluster == "15" and method == "gmm_numeric" and category == "conditions" and possibility == "adult icu": print x, y, t, n
                    continue
                    #run fisher's exact test
                    od, p = fisher_exact(np.array(([x,y],[t,n])), alternative="greater")
                    #add p to enriched clusters
                    enriched_clusters[method][cluster][category][possibility] = p
    exit()
    #return
    return enriched_clusters

def get_totals(database):
    """
    Get total number of patients with different stuffs
    """
    #initialize
    general = {}
    alls = {}
    for i in ["tissues", "treatments", "conditions", "references"]:
        general[i] = 0
        alls[i] = []
        #iterate
        for patient in database:
            if database[patient].parameters[i]: general[i] += 1
            if i != "references": alls[i] += database[patient].parameters[i]
            else: alls[i] += [database[patient].parameters[i]]
    #convert to totals
    tissue_dict = {}
    treatment_dict = {}
    condition_dict = {}
    reference_dict = {}
    #iterate
    for tissue in alls["tissues"]: tissue_dict[tissue] = alls["tissues"].count(tissue)
    for treatment in alls["treatments"]: treatment_dict[treatment] = alls["treatments"].count(treatment)
    for condition in alls["conditions"]: condition_dict[condition] = alls["conditions"].count(condition)
    for reference in alls["references"]: reference_dict[reference] = alls["references"].count(reference)
    #return
    return general, tissue_dict, treatment_dict, condition_dict, reference_dict 

def correct_clusters(enriched_clusters, alpha=0.05):
    """
    Correct clusters for multiple hypotheses using FDR
    """
    for method in enriched_clusters:
        #initialize p values
        p_values = []
        #get p values
        for cluster in enriched_clusters[method]:
            for category in enriched_clusters[method][cluster]:
                for specific in enriched_clusters[method][cluster][category]:
                    p_values.append(enriched_clusters[method][cluster][category][specific])
        #correct
        if p_values: truths, padjust, correctedp1, correctedp2 = multipletests(p_values, alpha*2, method = 'fdr_bh')
        else: continue
        #add back
        index = 0
        for cluster in enriched_clusters[method]:
            for category in enriched_clusters[method][cluster]:
                for specific in enriched_clusters[method][cluster][category]:
                    enriched_clusters[method][cluster][category][specific] = padjust[index]
                    index += 1
    #return
    return enriched_clusters
    
    
def write_to_file(outfile, enriched_clusters):
    """
    Write results to file
    """
    #open file
    with open(outfile, "w")  as outfl:
        #write header
        header = "Method\tCluster\tCategory\tSpecific\tPval\n"
        outfl.write(header)
        #iterate through lines
        for method in enriched_clusters:
            for cluster in enriched_clusters[method]:
                for category in enriched_clusters[method][cluster]:
                    for specific in enriched_clusters[method][cluster][category]:
                        #if line is worthy add it
                        if enriched_clusters[method][cluster][category][specific] <= 0.1:
                            line = "%s\tCluster_%s\t%s\t%s\t%s\n" % (method, cluster, category, specific, enriched_clusters[method][cluster][category][specific])
                            outfl.write(line)

def main(argv=None):
    #process command line
    settings = process_command_line(argv)
    #parse clusters
    clusters, patients = parse_clusters(settings.clusters)
    #parse database properties of relevant patients
    database = parse_database(settings.database, patients)
    #check enrichment
    enriched_clusters = cluster_enrichment(clusters, database, patients)
    #correct for multiple hypotheses
    enriched_clusters = correct_clusters(enriched_clusters)
    exit()
    #write to file
    write_to_file(settings.outfile, enriched_clusters)

if __name__ == "__main__":
        exit(main())
