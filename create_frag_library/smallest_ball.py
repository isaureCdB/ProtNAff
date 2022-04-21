'''
This script is calculating the smallest ball covering all the given points.
n is the number of points and dim is the dimension of points.
'''
from math import sqrt
import argparse
# from scipy.optimize import linprog
import numpy as np
from multiprocessing import Pool
from multiprocessing import get_context
import rmsdlib
import cProfile, pstats, io

def find_r(values, center, n, dim):
    '''
    This function calculate the r of the ball given:
    values: a matrix n.dim, the matrix of values of interest
    center: a vector of dimension dim, the center of the ball
    return: a value for r
    '''
    radius = float('-inf')
    for ii in range(n):
        r_tmp = 0
        for jj in range(dim):
            r_tmp += (values[ii][jj]-center[jj])**2
        if r_tmp > radius:
            radius = r_tmp
            # print(values[ii])
    return sqrt(radius)

def calculate_O(values, alpha, n, dim):
    '''
    This function is calculating the center of the ball given:
    values: a matrix n.dim, the matrix of values of interest
    alpha: a vector of dimension n, represents alpha
    n: number of points
    dim: dimension of points
    return: center: a vector of dimension dim, the center of the ball
    '''
    center = [0]*dim
    for ii in range(n):
        for jj in range(dim):
            center[jj] += alpha[ii]*values[ii][jj]
    return center

def calculate_rmsd(values):
    n = values.shape[0]
    matrix = np.zeros((n, n))
    for ii in range(n):
        matrix[ii, ii] = 0
        for jj in range(ii +1, n):
            tmp = rmsdlib.fit(values[ii], values[jj])[-1]
            matrix[ii, jj] = tmp
            matrix[jj, ii] = tmp
    return matrix

def calculate_K(values, n, processor):
    '''
    This function is calculating the matrix of scalar product given:
    values: a matrix n.dim, the matrix of values of interest
    n: number of points
    return: K: a n.n matrix, represents the scalar product of the points
    '''
    runargs = []
    for ii in range(n):
        runargs.append([ii, values])

    pool = get_context("spawn").Pool(processor)

    result = pool.map(run_parallel, runargs)
    pool.close()
    pool.join()
    # Concatenate the results into a matrix
    matrix_tmp = np.zeros((n, n))
    for jj in range(len(result)):
        line = result[jj][0]
        matrix_tmp[line, :] = result[jj][1]

    # As we calculated only the upper triangle, we need to copy it to the lower one.
    matrix_out = np.triu(matrix_tmp) + np.tril(matrix_tmp.T, 1) - np.diag(np.diag(matrix_tmp))
    return matrix_out

def run_parallel(runarg):
    ref = runarg[0]
    coordinates = runarg[1]
    n = coordinates.shape[0]
    atoms_ref = coordinates[ref]
    results = np.zeros(n)
    # Calculating only the upper triangle
    for ii in range(ref,n):
        atoms_ii = coordinates[ii]
        # Fit atoms_ii on atoms_ref
        results[ii] = scalar_product(atoms_ref, atoms_ii)
    return ref, results

def scalar_product(atoms1, atoms2):
    '''
    This function is calculating the calar product between two structures.
    INPUT: atoms1 and atoms2, the structures of shape (n,m)
    OUTPUT: the scalar product named sp
    '''
    sp = 0
    for ii in range(atoms1.shape[0]):
        sp += atoms1[ii]*atoms2[ii]
    return sp

def calculate_dual(alpha, K, kappa):
    '''
    This function is calculating the dual problem given:
    alpha: a vector of dimension n, represents alpha
    K: a n.n matrix, represents the scalar product of the points
    kappa: a vector of dimension n, the diagonale of K
    n: number of points
    return: a value for the dual
    '''
    dual = -(np.dot(alpha, np.dot(K, alpha)) - np.dot(kappa, alpha))
    return dual

def run_calculation(values, processor):
    values.shape = (values.shape[0], values.shape[1]*values.shape[2])
    n = values.shape[0]
    dim = values.shape[1]
    processor = min([processor, n])
    # print("Number of points: ", n)
    # print("Dimension of points: ", dim)
    # K = np.load("K_test.npy")
    K = calculate_K(values, n, processor)
    # np.save("K_test.npy", K)
    kappa = K.diagonal()
    #np.save("K{}.npy".format(i), K)

    #Initialisation of alpha
    alpha = np.array([1/n]*n, dtype=np.float32)
    #Initialisation ratio and A amtrix used for linear optimisation
    ratio = 0
    A = np.array([[1]*n])
    i = 0

    while ratio < 0.95:
        gradient = np.dot(alpha, K) * 2 - kappa
        # Pour les vi > 0 dj fait par la fonction linprog
        # res = linprog(gradient, A_eq=A, b_eq=1)
        # v = res.x
        v = np.zeros(n)
        v[np.argmin(gradient)] = 1

        tmp = np.dot(np.transpose(alpha-v), K) * 2
        teta = (np.dot(gradient, (alpha - v))/np.dot(tmp, (alpha - v)))
        teta = np.abs(teta.item())
        teta = min(1, teta)
        alpha_t1 = alpha * (1 - float(teta)) + v * float(teta)
        if i == 100:
            O = calculate_O(values, alpha_t1, n, dim)
            R = find_r(values, O, n, dim)
            dual = calculate_dual(alpha_t1, K, kappa)
            ratio = dual/(R**2)
            # print(R)
            # print(ratio)
            i = 0
            # print("Current ratio: ", ratio)
            # print("")
        #for j in range(alpha.shape[0]):
        #    if alpha_t1[j] < 0.0001:
        #        alpha_t1[j] = 0

        i += 1
        alpha = alpha_t1
    return np.array(O), R

def main():
    '''
    The main function to run
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument("values", help="The numpy matrix of values, should have \
    the shape (n,p,3)")
    parser.add_argument("--np", help="The number of processors to use.", type=int, default=1)
    args = parser.parse_args()
    #number of values
    values = np.load(args.values)
    processor = args.np

    pr = cProfile.Profile()
    pr.enable()
    O, R = run_calculation(values, processor)
    print("Radius = ", R)
    print("center: ", O)
    pr.disable()
    s= io.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print(s.getvalue())

if __name__ == "__main__":
    main()
