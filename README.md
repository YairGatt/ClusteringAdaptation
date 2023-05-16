# ClusteringAdaptation
Main scripts and utility scripts used for the paper "Clinical characteristics and pre-existing genetic variation determine some adaptive strategies employed by Staphylococcus aureus strains during prolonged infection and carriage" by Gatt and Margalit (unpublished)

Main scripts are written in tcsh and can be run from the shell:
-Clustering.sh: This script recieves the results of the Within Host Adaptation pipeline (https://github.com/YairGatt/WithinHostAdaptation) and clusters the different strains within each chosen organism based on KEGG pathways including genes undergoing mutation during host adaptation. Different clustering methods are applied and can be compared.
-Enrichment.sh: This script assess the enrichment of the different clusters outputted by Clustering.sh by different KEGG pathways, in order to clearly define which pathways undergo mutation in each cluster.
-Cluster_properties.sh: This script assess the enrichment of the strains included in the different clusters by clinical properties including antibiotic treatment, tissue from which samples were isolated, and more.
