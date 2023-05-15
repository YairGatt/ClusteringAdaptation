#! /bin/tcsh

#Clustering pipeline for the reuslts of the Within Host Adaptation pipeline
#Input the results directory of Within Host Adaptation pipeline for a single organism, and the different strains will be clustered using different methods.

set BINDIR = $1 #home directory of the github
set organism = $2 #organism dir with results of Within Host Adaptation pipeline
set the_organism = `basename $organism`
set workdir = $3 #directory for writing results
set reference_db = $4 #directory with genome in fasta format
set pathways = $5 #KEGG pathways .list file for the organism
#results dir
mkdir -p $workdir/high/${the_organism}
mkdir -p $workdir/changed/${the_organism}
mkdir -p $workdir/high/${the_organism}/clustering/
mkdir -p $workdir/changed/${the_organism}/clustering/
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
    set high_vector = `echo $high_vector ${patient}/high_lost_genes_no_repeats.txt`
  endif
  if ("$changed" != "") then
    set changed_vector = `echo $changed_vector ${patient}/changed_genes_no_repeats.txt`
  endif
end
#high_lost_genes
python $BINDIR/bin/Clustering/clustering.py -d $reference_db -k $pathways -l $high_vector -c kmeans -w $workdir/high/${the_organism}
python $BINDIR/bin/Clustering/clustering.py -d $reference_db -k $pathways -l $high_vector -c gmm -w $workdir/high/${the_organism}
python $BINDIR/bin/Clustering/clustering.py -d $reference_db -k $pathways -l $high_vector -c hierarchial -w $workdir/high/${the_organism}
python $BINDIR/bin/Clustering/clustering.py -d $reference_db -k $pathways -l $high_vector -c spectral -w $workdir/high/${the_organism}
#changed_genes
python $BINDIR/bin/Clustering/clustering.py -d $reference_db -k $pathways -l $changed_vector -c kmeans -w $workdir/changed/${the_organism}
python $BINDIR/bin/Clustering/clustering.py -d $reference_db -k $pathways -l $changed_vector -c gmm -w $workdir/changed/${the_organism}
python $BINDIR/bin/Clustering/clustering.py -d $reference_db -k $pathways -l $changed_vector -c hierarchical -w $workdir/changed/${the_organism}
python $BINDIR/bin/Clustering/clustering.py -d $reference_db -k $pathways -l $changed_vector -c spectral -w $workdir/changed/${the_organism}
