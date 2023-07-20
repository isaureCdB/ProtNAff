cd database/trilib
for m1 in A C; do 
  for m2 in A C; do 
    for m3 in A C; do 
      m=$m1$m2$m3; 
      python3 ../../build-library.py ../fragments.json $m 0.5 1.0 \
      $m-aa-fit-clust0.2 $m-aa-fit-clust0.2.npy  \
      $m-recluster-0.5-all.list $m-recluster-0.5-closest-cluster.rmsd \
      $m-lib \
      | tee build-library-$m.log
      echo $m; 
    done
  done
done