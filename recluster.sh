for m1 in A C; do
  for m2 in A C; do
    for m3 in A C; do
      m=$m1$m2$m3
      python3  recluster.py database/trilib/$m-aa-fit-clust0.2.npy 1.0 database/trilib/$m-recluster-1.0-bigcluster.list database/trilib/$m-recluster-prmsd.npy; 
    done
  done
done