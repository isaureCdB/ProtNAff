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
      mutm=$(echo $m | sed 's/A/G/g' | sed 's/C/U/g')
      for domut1 in 0 1; do
        for domut2 in 0 1; do
          for domut3 in 0 1; do
            domut=$domut1$domut2$domut3
            if [ $domut -eq 000 ]; then
              continue
            fi
            mut=$(awk -v m=$m -v domut=$domut -v mutm=$mutm \
              'BEGIN{
                r="" 
                for (n=0;n<length(m); n++) {
                  if (int(substr(domut,n+1,1))) v=mutm; else v = m
                  r = r substr(v,n+1,1)
                } 
                print(r)
              }'
            )
            echo $m $mut
            for f in '.npy' '-conformer.npy' '-intracluster.npy'; do
              python3 ../../create_frag_library/mutate_in_heavy_atom.py $m $mut $m-lib$f $mut-lib$f
              echo 
            done
            for f in '.list' '-conformer.list' '-intracluster.list' '-secondary.list' '-tertiary.list' '-quaternary.list' '-quasi-unique.list'; do
              rm -f $mut-lib$f
              ln -s $m-lib$f $mut-lib$f
            done
          done
        done
      done      

    done
  done
done