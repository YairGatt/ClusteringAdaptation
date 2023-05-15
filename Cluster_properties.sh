#!/bin/tcsh

#set workdir from Clustering.sh as workdir to get the properties of the different clusters

set BINDIR = $1 #Directory where all scripts are found
set workdir = $2 #Results directory from Clustering.sh
set Table_S1 = $3 #Table S1 from Gatt and Margalit 2020, or from results of Within Host Adaptation pipeline
foreach clusters ($workdir/high/*/clustering/Clusters.txt $workdir/changed/*/clustering/Clusters.txt)
  set outfile = `echo $clusters | sed 's/\.txt/_properties.txt/'`
  set preoutdir = `dirname $clusters`
  set outdir = $preoutdir/intersection/
  python $BINDIR/bin/Clustering/Get_cluster_properties.py -db $Table_S1 -c $clusters -o $outfile
  python $BINDIR/bin/Clustering/Properties_and_pathways.py -c $clusters -p $outfile -w $outdir
end
