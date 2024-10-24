Querying library from the camembert project ( https://github.com/sjdv1982/camembert)
Author: Sjoerd de Vries.

How to run:

# re-run assign_clusters.py based on the reclustering
python3 ../create_frag_library/assign_clusters.py ../database/fragments.json \
    --clustdir ../database/trilib \
    --clustfiles "aa-fit-clust0.2" "recluster-0.5-closest.list" \
    --clustnames clust0.2A clust0.5A \
    --fast \
    ../database/fragments_clust.json rna
python3 ../create_frag_library/assign_clusters.py ../database/fragments_clust.json \
    --clustdir ../database/trilib \
    --clustfiles "aa-fit-clust0.2" "recluster-1.0-closest.list" \
    --clustnames clust0.2A clust1A \
    --fast \
    ../database/fragments_clust.json rna

# convert fragments_clust to npy format
python3 json2npy.py ../database/fragments_clust.json ../database/fragments_clust.npy

query with chaindata=structures.json and frag=fragments_clust.npy (loaded)