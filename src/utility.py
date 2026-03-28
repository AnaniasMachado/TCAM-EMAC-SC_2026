import numpy as np
import os
import scipy.io
import scipy.sparse
import re

epsilon = 10 ** -5

def vec_1_norm(x):
    return np.linalg.norm(x, ord=1)

def matrix_vec_1_norm(H):
    return np.linalg.norm(H.flatten(), ord=1)

def matrix_vec_0_norm(H):
    total = 0
    for x in H.flatten():
        if x > epsilon:
            total += 1
    return total

def matrix_vec_inf_norm(H):
    return max([np.abs(h) for h in H.flatten()])

def matrix_frobenius_norm(H):
    return np.linalg.norm(H, ord="fro")

def matrix_2_1_norm(H):
    total = 0
    for i in range(H.shape[0]):
        total += np.linalg.norm(H[i, :], ord=2)
    return total

def matrix_2_0_norm(H):
    total = 0
    for i in range(H.shape[0]):
        if np.linalg.norm(H[i, :], ord=2) > epsilon:
            total += 1
    return total

def matrix_rank(H):
    U, S, VT = np.linalg.svd(H)
    rank = 0
    for i in range(S.shape[0]):
        if S[i] > epsilon:
            rank += 1
    return rank

def get_instances(matrices_folder):
    matrices_filepath = []
    for filename in os.listdir(matrices_folder):
        filepath = os.path.join(matrices_folder, filename)
        if os.path.isfile(filepath):
            matrix_filepath = f"{matrices_folder}/{filename}"
            matrices_filepath.append(matrix_filepath)
    return matrices_filepath

def read_instance(instance, name):
    # Loads file .mat
    mat = scipy.io.loadmat(instance)

    # Reads matrix
    matrix = mat[name]

    if scipy.sparse.issparse(matrix):
        # Converts to dense matrix
        dense_matrix = matrix.toarray()
        return dense_matrix
    return matrix

def read_maragal_instance(instance):
    # Loads file .mat
    mat = scipy.io.loadmat(instance)

    # Reads matrix
    data = mat["Problem"]
    matrix = data["A"][0, 0]

    if scipy.sparse.issparse(matrix):
        # Converts to dense matrix
        dense_matrix = matrix.toarray()
        return dense_matrix
    return matrix

def get_data_from_instance(instance):
    # Using regex to find values of m, n, r and d
    match = re.search(r'm(\d+)_n(\d+)_r(\d+)_d(\d+)_idx(\d+)', instance)

    if match:
        m = int(match.group(1))
        n = int(match.group(2))
        r = int(match.group(3))
        d = int(match.group(4))
        idx = int(match.group(5))
        
        return [m, n, r, d / 100, idx]
    else:
        raise Exception("ReadError: Could not read all values.")