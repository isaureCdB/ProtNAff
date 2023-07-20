for m1 in A C; do
  for m2 in A C; do
    for m3 in A C; do
      for t in 1.0 0.5; do
        m=$m1$m2$m3
        python3  rerecluster.py  database/trilib/$m-aa-fit-clust0.2.npy database/trilib/$m-recluster-1.0-bigcluster.list database/trilib/$m-recluster-prmsd.npy 1.0 $t database/trilib/$m-recluster-$t-all.list database/trilib/$m-recluster-$t-closest.list database/trilib/$m-recluster-$t-closest-cluster.rmsd
      done
    done
  done
done
