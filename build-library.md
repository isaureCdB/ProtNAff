# Converting a fragment clustering to a fragment library for docking

## Definitions

**Cis-modeling**: docking PDB codes that were used to build the library
**Trans-modeling**: docking PDB codes that were not used to build the library

### Colored precision

**Blue** is considered within clustering range (0.5 A).

**Green** is outside clustering range, but is considered as "related". When a continuous region of conformations is clustered, the neighbouring clusters will have an RMSD that is green (and close to blue if there are many trinucleotides in that region). In addition, for docking, it is considered that clusters within "green" RMSD can replace each other reasonably well. Proposed RMSD threshold for "green" : 1.0 A

**Red** is considered as "(possibly) unrelated". For docking, clusters with a "red" RMSD between them don't replace each other very well. Clusters whose closest neighbour cluster is at a "red" RMSD distance are potentially an "island" rather than a continuum. "Island" clusters that come from a single PDB are potentially unique to that PDB.

**Singletons**: clusters that come from a single PDB code.

## Classification

Take all 0.5 A clusters. First, divide the clusters into singletons and non-singletons. 
Then, classify a cluster as "red" or "green"  based on the RMSD towards the closest neighbour cluster.

***Special case***: singletons whose closest neighbour cluster is a singleton from the same PDB code. In this case, search the next-closest neighbour cluster (that is not a singleton from the same PDB code). In any case, replace the closest neighbour with the next-closest and process it using the main procedure. However, if this changes the singleton from green to red, put it aside in a special "quasi-unique" category, which is explained below.

### Main procedure

The goal is to produce several cluster lists for docking. Trans-modeling uses only the primary list. For cis-modeling, elements from the primary list are replaced by elements from the secondary list, tertiary list, quaternary list, or the quasi-unique list, in order to make the docking blind, i.e. not using the conformation of the PDB that is being modeled. In cis-modeling, for the PDB that is being modeled, the secondary list is always being used. One may choose to use the other lists, or one may choose to simply remove the cluster if the cluster heart comes from the PDB that is being modeled.
The lists do not contain any coordinates, only indices from fragments.json. There is a master.list that contains all fragment indices from all lists (eliminating duplicates) except intra-cluster replacements (secondary list and part of the tertiary) and a master.npy containing their coordinates. There is a intracluster pair list and a intracluster.npy with the coordinates of the intracluster replacement.

The clusters are classified as:

- Red non-singletons.
- Green non-singletons.
- Type I green singletons. These are within green RMSD of at least one non-singleton cluster.
- Type II green singletons. Only other singletons are within green RMSD.
- Red singletons.

## Red non-singletons

In this case, it is imperative to have a replacement available for cis-modeling.
Add the fragment index (from fragments.json) to the primary list.
Identify the PDB code of the cluster heart. Identify the cluster member closest to the heart with a different PDB code. Add the fragment index of the cluster member to the secondary list.

## Green non-singletons

The same as for red singletons, but add to the tertiary list instead. The tertiary list also contains the RMSD, the cluster size (at the 0.2 A level), and the number of PDB codes.
Using the tertiary list is optional (you can just discard this cluster with no replacement).

## Type I green singletons

These singletons are close to a larger cluster. Clearly, they belong to a region of conformational space that is reasonably well-sampled. The larger cluster represents this region already. Therefore, the singleton can simply be discarded, implicitly represented by its larger neighbour, both in cis- and trans-modeling.
In practical terms, there is no scenario where keeping this cluster would make cis-modeling any more accurate. For trans-modeling, keeping the cluster would make the result more accurate only if there is a significant chance of unknown conformations that are rather close to the singleton cluster but rather far from its larger neighbour. If this is anticipated, one may choose to treat clusters with a neighbour cluster at the high end of green RMSD (or with a neighbour cluster that is rather small) as type II green singletons instead.
The eliminated singletons are stored in eliminated.list, with the same format as the tertiary list.

## Type II green singletons

These are the most difficult to interpret. Since they are in a sparsely populated region of conformational space, it is possible that such a cluster is unique to the PDB where it came from (like red singletons, see below). Alternatively, it is possible that the cluster and its neighbour are both one-point samples of a continuous region of conformational space. It is most prudent to add the cluster to the primary list and to add its nearest-neighbour to the tertiary list. Alternatively, especially if the RMSD is near the lower end, one may choose to treat such cluster as a type I green singleton, i.e. eliminate it, on the condition that the closest neighbour is not eliminated.

## Red singletons

For cis-modeling, it has no practical value to add the cluster itself to any list: only the nearest-neighbour cluster will ever be relevant. However, as an alternative, one might consider that red singletons on the lower end of RMSD are more appropriately treated as type II green singletons, i.e. as a larger cluster that has been undersampled. If this alternative is rejected, then the singleton is considered to be unique to the PDB where it came from. In that case, it will .not be used in trans-modeling. For cis-modeling, proceed as follows:

Add the index of the nearest-neighbour cluster and the RMSD to the quaternary list. Indices are being stored so that duplicates are being avoided (the added cluster may or may not be in the primary list as well). In cis-modeling, one may use to choose the quaternary list if the RMSD is considered low enough for docking to succeed. Else, one must give up altogether on modeling the fragment.

## Quasi-unique singletons

Note that these have already been treated as red singletons. However, for cis-modeling, where the secondary replacement is too imprecise, one may choose to use a neighbour cluster that comes from the same PDB (pseudo-blind docking). Therefore, create a quasi-unique list with (cluster, closest neighbour) pairs for each PDB code (as fragment indices). When docking, only add the closest neighbour, not the cluster itself (although the reciprocal pair may also exist).
