#! /bin/tcsh

#Assess enrichment of clusters by metabolic pathways

set BINDIR = $1 #home directory of the github
set organism = $2 #organism dir with results of Within Host Adaptation pipeline
set the_organism = `basename $organism`
set workdir = $3 #directory for writing results
set reference_db = $4 #directory with genome in fasta format
set pathways = $5 #KEGG pathways .list file for the organism
#results dir
mkdir -p $workdir/high/${the_organism}
mkdir -p $workdir/changed/${the_organism}
#initialize results vector
set changed_vector = ""
set high_vector = ""
#iterate through relevant trial
foreach patient ($organism/trial_*_patient_*/)
  #define lost genes file
  set changed = `ls $patient | grep "changed_genes_no_repeats.txt" | grep -v "rereading"`
  set high = `ls $patient | grep "high_lost_genes_no_repeats.txt" | grep -v "rereading"`
  #check if loss assessed
  if ("$high" != "") then
    set high_vector = `echo $high_vector ${patient}/${high}`
  endif
  if ("$changed" != "") then
    set changed_vector = `echo $changed_vector ${patient}/${changed}`
  endif
end
#initialize files
echo $the_organism > $workdir/changed/${the_organism}/significant_genes.txt
echo $the_organism > $workdir/high/${the_organism}/significant_genes.txt
echo $the_organism > $workdir/changed/${the_organism}/significant_genes_metabolics_only.txt
echo $the_organism > $workdir/high/${the_organism}/significant_genes_metabolics_only.txt

#run regular mode
python $BINDIR/bin/Clustering/enrichment.py -d $reference_db -k $pathways -l $high_vector -w $workdir/high/${the_organism}/ >> $workdir/high/${the_organism}/significant_genes.txt
python $BINDIR/bin/Clustering/enrichment.py -d $reference_db -k $pathways -l $changed_vector -w $workdir/changed/${the_organism}/ >> $workdir/changed/${the_organism}/significant_genes.txt
#run metabolic only mode
python $BINDIR/bin/Clustering/enrichment.py -d $reference_db -k $pathways -l $high_vector -w $workdir/high/${the_organism}/ --metabolic >> $workdir/high/${the_organism}/significant_genes_metabolic_only.txt
python $BINDIR/bin/Clustering/enrichment.py -d $reference_db -k $pathways -l $changed_vector -w $workdir/changed/${the_organism}/ --metabolic >> $workdir/changed/${the_organism}/significant_genes_metabolics_only.txt
