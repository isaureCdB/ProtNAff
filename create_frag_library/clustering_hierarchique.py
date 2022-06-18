#!/usr/bin/env python3

"""
This script is doing a hierarchical clustering based on the fact that we want
specific clusters, there are 2 rules:
    - first, all values need to be at most at a threshold from on center
    - all center must be further from that threshold

INPUT:
    - a numpy matrix of coordinates

OUTPUT:
    - a dictionnary of clusters containing the list of values in those clusters.
"""

import numpy as np
import argparse
# Those libraries are ours, check rmsdlib.py and smallest_ball_quadra.py
#from smallest_ball_quadra import return_center_radius
from smallest_ball import run_calculation
import rmsdlib
from multiprocessing import Pool

def calculate_rmsd_matrix(values, processor):
    '''
    This function is calculating the matrix of scalar product given:
    values: a matrix n.dim, the matrix of values of interest
    n: number of points
    return: K: a n.n matrix, represents the scalar product of the points
    '''

    runargs = []
    n = values.shape[0]
    for ii in range(n):
        runargs.append([ii, values])

    pool = Pool(processor)

    result = pool.map(run_parallel, runargs)
        
    return np.concatenate([r[1] for r in result]).reshape(n, n)
    

def run_parallel(runarg):
    ref = runarg[0]
    coordinates = runarg[1]
    n = coordinates.shape[0]
    values_ref = coordinates[ref]
    results = rmsdlib.multifit(coordinates, values_ref)[2]
    return ref, results

def mean_prototype(prototype1, prototype2, size1, size2):
    """
    This function is calculating the average structures taking into account
    2 centers of cluster and the size of those clusters.
    INPUT: prototype1 and prototype2 are the prototypes of the 2 clusters
    size1 is the number of structures in the cluster 1
    size2 is the number of structures in the cluster 2
    OUTPUT: An average structure
    """
    mean = np.copy(prototype1)
    sum_size = float(size1 + size2)
    for ii in range(prototype1.shape[0]):
        for jj in range(prototype1.shape[1]):
            mean[ii, jj] = float(size1/sum_size)*prototype1[ii, jj] + \
                            float(size2/sum_size)*prototype2[ii, jj]
    return mean

def align_values_clusters(prototype, values_cluster):
    """
    This function is doing a "all on one" alignement of structures.
    INPUT: prototype is the structure you will align on (mean because it is
    the prototype of a cluster)  (shape = (m,3))
    values_cluster is a set of structures that you will align on mean
    (shape = (n',m,3))
    OUTPUT: a set of structures that was moved (translation + rotation) to be
    alugned on mean (shape = (n',m,3))
    """
    #Here we do an alignement of type "one on all"
    rotation, translation, RMSD = rmsdlib.multifit(values_cluster, prototype)
    #Then the rotations are transposed
    rot = np.transpose(rotation, axes=(0, 2, 1))
    #Calculation of the Center Of Mass (COM)
    COM = values_cluster.sum(axis=1)/values_cluster.shape[1]
    #Definition of new centers of structures
    centered = values_cluster - COM[:, None, :]
    rotated = np.einsum('...ij,...jk->...ik', centered, rot)
    fitted = rotated + COM[:, None, :]
    translated = fitted - translation[:, None, :]
    return translated

def check_fusion(values, prototype, threshold):
    """
    This function is testing if the fusion is kept or not, it means that all
    values in the cluster are at less than the threshold of the prototype.
    INPUT: values a set of structures contained into the cluster
    (shape = (n',m,3))
    prototype the center of the cluster (shape = (m,3))
    threshold the "width" of the cluster
    OUTPUT: a boolean, True if the fusion can be done, False else
    """
    result = True
    for ii in range(values.shape[0]):
        tmp = rmsdlib.fit(values[ii], prototype)
        if tmp[-1] > threshold:
            result = False
            break
    return result

def calculate_new_rmsd(matrix_rmsd, min_clust, max_clust, values_prototype):
    """
    This function is modifing the matrix_rmsd to take into account the fusion
    of clusters, it is putting the line and the column max_clust at the value
    None, and calculatng the line and the column of min_clust by aligning
    the prototype on all the others prototypes.
    INPUT: matrix_rmsd, the actual rmsd matrix (shape = (n,n))
    min_clust, max_clust the two clusters that are merged
    values_prototype the set of prototype structures (shape = (n',m,3))
    OUTPUT: The new rmsd matrix
    """

    n = values_prototype.shape[0]
    for ii in range(n):
        if ii != min_clust and ii != max_clust:
            #If the value is NaN, then it means that this cluster is not a
            #cluster, so do nothing
            if not np.isnan(matrix_rmsd[ii, min_clust]):
                matrix_rmsd[ii, min_clust] = rmsdlib.fit(values_prototype[ii],\
                                                values_prototype[min_clust])[-1]
                matrix_rmsd[min_clust, ii] = matrix_rmsd[ii, min_clust]
        else:
            matrix_rmsd[ii, ii] = None
        matrix_rmsd[ii, max_clust] = None
        matrix_rmsd[max_clust, ii] = None
    return matrix_rmsd

def calculate_new_clusters(clusters_frag):
    """
    This funtion is removing from the dictionnary the empty clusters.
    INPUT: clusters_frag a dictionnary of clusters containing lists of values
    OUTPUT: the same dictionnary but with empty clusters removed
    """
    to_pop = []
    for key in clusters_frag.keys():
        if not clusters_frag[key]:
            to_pop.append(key)
    for item in to_pop:
        clusters_frag.pop(item)
    return clusters_frag

def try_fusion_parallel(runarg):
    x = runarg[0]
    y = runarg[1]
    values_prototype = runarg[2]
    values = runarg[3]
    clusters_frag = runarg[4]
    d = runarg[5]
    t = runarg[6]
    mean_prot = mean_prototype(values_prototype[x], values_prototype[y], 1, 1)
    values_x = align_values_clusters(mean_prot, values[clusters_frag[x]])
    values_y = align_values_clusters(mean_prot, values[clusters_frag[y]])
    values_cluster = np.concatenate((values_x, values_y))

    # Regarder difference entre centre de la boule et prototype moyen
    center, radius = run_calculation(values_cluster, 5)
    values_cluster.shape = (values_cluster.shape[0], values_cluster.shape[1]//d, d)
    center.shape = (int(center.shape[0]/d), d)
    smallest_fusion = check_fusion(values_cluster, center, t)
    return x, y, smallest_fusion, center

def main():
    """
    This is the main script, it is made to be run as it is.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("values", help="The numpy matrix of values, should have \
    the shape (n,m,3)")
    parser.add_argument("threshold1", help="The width of the first final clusters.")
    parser.add_argument("threshold2", help="The width of the second final clusters.")
    parser.add_argument("-c", "--cpu", help="The number of CPUS to use", dest="cpu")
    parser.add_argument("-o", "--output", help="The base name for the motif", dest="output")
    args = parser.parse_args()
    i = 0
    values = np.load(args.values)
    values_prototype = np.copy(values)
    m = values.shape[1]
    d = values.shape[2]
    count = values.shape[0]
    t = float(args.threshold1)
    t2 = float(args.threshold2)
    matrix_rmsd = calculate_rmsd_matrix(values, int(args.cpu))
    # matrix_rmsd = np.load("matrix_rmsd_poster.npy")
    # np.save("matrix_rmsd_poster.npy", matrix_rmsd)
    # np.save("values_poster0.npy", values)
    #Initialisation of clusters as each value is a prototype.
    clusters_frag = {}
    for ii in range(values_prototype.shape[0]):
        clusters_frag[ii] = [ii]

    possible_fusion = True
    copy_rmsd = np.copy(matrix_rmsd)

    while possible_fusion:
        tag = True

        while not np.all((np.isnan(copy_rmsd))) and tag == True:
            i += 1
            try:
                minimum_indices = np.nanargmin(copy_rmsd)
            except ValueError:
                possible_fusion = False
                print("End of clustering.")
                break

            x = minimum_indices % matrix_rmsd.shape[1]
            y = minimum_indices // matrix_rmsd.shape[1]
            if copy_rmsd[x][y] > 2*t:
                possible_fusion = False
                print("End of clustering.")
                break
            mean_prot = mean_prototype(values_prototype[x], values_prototype[y], 1, 1)
            values_x = align_values_clusters(mean_prot, values[clusters_frag[x]])
            values_y = align_values_clusters(mean_prot, values[clusters_frag[y]])
            values_cluster = np.concatenate((values_x, values_y))

            center, radius = run_calculation(values_cluster, int(args.cpu))
            values_cluster.shape = (values_cluster.shape[0], values_cluster.shape[1]//d, d)
            center.shape = (int(center.shape[0]//d), d)
            smallest_fusion = check_fusion(values_cluster, center, t)

            if smallest_fusion:
                count -= 1
                # print("Nb clusters = ", count)
                min_clust = min(x, y)
                max_clust = max(x, y)
                clusters_frag[min_clust] = clusters_frag[x] + clusters_frag[y]
                clusters_frag[max_clust] = []
                values_prototype[min_clust] = center
                matrix_rmsd = calculate_new_rmsd(matrix_rmsd, min_clust, max_clust, \
                                values_prototype)
                copy_rmsd[max_clust] = None
                copy_rmsd[:,max_clust] = None
                for ii in range(values_prototype.shape[0]):
                    if not np.isnan(matrix_rmsd[ii, min_clust]):
                        copy_rmsd[ii,min_clust] = matrix_rmsd[ii, min_clust]
                        copy_rmsd[min_clust,ii] = matrix_rmsd[ii, min_clust]
                tag = False
            else:
                # Tester toutes les fusions possibles.
                # Voir le gain possible en utilisant un seuil pour ne pas faire la fusion
                copy_rmsd[x,y] = None
                copy_rmsd[y,x] = None
                if np.all((np.isnan(copy_rmsd))):
                    possible_fusion = False

        if np.all((np.isnan(copy_rmsd))):
            possible_fusion = False


    clusters_frag_1A = calculate_new_clusters(clusters_frag)

    values_prototype_1A = values_prototype[list(clusters_frag.keys())]

    np.save("{}-dr0.2r-clust{:.1f}-aa".format(args.output, t), values_prototype_1A)
    np.save("{}_clusters".format(args.output), clusters_frag_1A)

    res = []
    for ii in range(1, len(list(clusters_frag_1A.keys()))+1):
        string_res = "Cluster {} ->".format(ii)
        for jj in clusters_frag_1A[list(clusters_frag_1A.keys())[ii-1]]:
            string_res += " {}".format(jj+1)
        res.append(string_res+'\n')
    with open("{}-dr0.2r-clust{}".format(args.output, t), 'w') as ff:
        for line in res:
            ff.write(line)

    possible_fusion = True

    while possible_fusion:
        tag = True

        while not np.all((np.isnan(copy_rmsd))) and tag == True:
            i += 1
            try:
                minimum_indices = np.nanargmin(copy_rmsd)
            except ValueError:
                possible_fusion = False
                print("End of clustering.")
                break

            x = minimum_indices % matrix_rmsd.shape[1]
            y = minimum_indices // matrix_rmsd.shape[1]
            if copy_rmsd[x][y] > 2*t2:
                possible_fusion = False
                print("End of clustering.")
                break
            mean_prot = mean_prototype(values_prototype[x], values_prototype[y], 1, 1)
            values_x = align_values_clusters(mean_prot, values[clusters_frag[x]])
            values_y = align_values_clusters(mean_prot, values[clusters_frag[y]])
            values_cluster = np.concatenate((values_x, values_y))

            center, radius = run_calculation(values_cluster, int(args.cpu))
            values_cluster.shape = (values_cluster.shape[0], values_cluster.shape[1]//d, d)
            center.shape = (int(center.shape[0]//d), d)
            smallest_fusion = check_fusion(values_cluster, center, t2)

            if smallest_fusion:
                count -= 1
                # print("Nb clusters = ", count)
                min_clust = min(x, y)
                max_clust = max(x, y)
                clusters_frag[min_clust] = clusters_frag[x] + clusters_frag[y]
                clusters_frag[max_clust] = []
                values_prototype[min_clust] = center
                matrix_rmsd = calculate_new_rmsd(matrix_rmsd, min_clust, max_clust, \
                                values_prototype)
                copy_rmsd[max_clust] = None
                copy_rmsd[:,max_clust] = None
                for ii in range(values_prototype.shape[0]):
                    if not np.isnan(matrix_rmsd[ii, min_clust]):
                        copy_rmsd[ii,min_clust] = matrix_rmsd[ii, min_clust]
                        copy_rmsd[min_clust,ii] = matrix_rmsd[ii, min_clust]
                tag = False
            else:
                # Tester toutes les fusions possibles.
                # Voir le gain possible en utilisant un seuil pour ne pas faire la fusion
                copy_rmsd[x,y] = None
                copy_rmsd[y,x] = None
                if np.all((np.isnan(copy_rmsd))):
                    possible_fusion = False

        if np.all((np.isnan(copy_rmsd))):
            possible_fusion = False

    clusters_frag_3A = calculate_new_clusters(clusters_frag)

    values_prototype_3A = values_prototype[list(clusters_frag.keys())]

    np.save("{}-dr0.2r-clust{:.1f}-clust{:.1f}".format(args.output, t, t2), values_prototype_3A)
    np.save("{}_clusters_3A".format(args.output), clusters_frag_3A)

    res = []
    for ii in range(1, len(list(clusters_frag_3A.keys()))+1):
        string_res = "Cluster {} ->".format(ii)
        for jj in clusters_frag_3A[list(clusters_frag_3A.keys())[ii-1]]:
            string_res += " {}".format(jj+1)
        res.append(string_res+'\n')
    with open("{}-dr0.2r-clust{}-clust{}".format(args.output, t, t2), 'w') as ff:
        for line in res:
            ff.write(line)

if __name__ == "__main__":
    main()
