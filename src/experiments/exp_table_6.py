import time
import os
import numpy as np
import pandas as pd
from scipy.io import savemat

from ..utility import *
from ..methods.local_search import *

problems = ["P123", "P1Sym"]
problem = problems[1]

column_names = ["m", "n", "r", "d", "||A||_0", "||A^+||_0", "||A^+||_1", "||H||_0" , "||H||_1", "time"]

exp = "table_6"
time_limit = 86400

matrices_folder = f"./instances/{exp}"
results_folder = f"./results/{exp}"
solutions_folder = f"./solutions/{exp}"

methods = {
    "P123": "LSLAFI_Det_P3",
    "P1" : "LSLAFI_Det",
    "P1Sym" : "LSLAFI_Det_Symmetric"
}
method = methods[problem]

df = pd.DataFrame(columns=column_names)

instances = get_instances(matrices_folder)

min_unsolvable_m = np.inf

for i in range(1, 6):
    A_norm_0 = 0
    AMP_norm_0 = 0
    AMP_norm_1 = 0
    H_norm_0 = 0
    H_norm_1 = 0
    runtime = 0

    instance = f"{matrices_folder}/Maragal_{i}"
    print(f"Solving for matrix: {instance}")
    A = read_maragal_instance(instance)
    m, n = A.shape
    r = matrix_rank(A)
    d = matrix_vec_0_norm(A) / (m * n)
    d = int(d * 100)
    hatA_flag = 1
    if m < min_unsolvable_m:

        if hatA_flag:
            A = A.T @ A
            r = matrix_rank(A)

        A_MP = np.linalg.pinv(A)

        start_time = time.time()
        H = ls_maragal(exp, solutions_folder, problem, instance, method, hatA_flag, time_limit)
        end_time = time.time()

        if not isinstance(H, int):
            A_norm_0 = matrix_vec_0_norm(A)
            AMP_norm_0 = matrix_vec_0_norm(A_MP)
            AMP_norm_1 = matrix_vec_1_norm(A_MP)
            H_norm_0 = matrix_vec_0_norm(H)
            H_norm_1 = matrix_vec_1_norm(H)
            runtime = end_time - start_time
        else:
            print("TimeLimit: Could not solve instance within the time limit.")
            min_unsolvable_m = m
    instance_results = {
        "m": m,
        "n": n,
        "r": r,
        "d": d,
        "||A||_0": A_norm_0,
        "||A^+||_0": AMP_norm_0,
        "||A^+||_1": AMP_norm_1,
        "||H||_0": H_norm_0,
        "||H||_1": H_norm_1,
        "time": runtime
    }

    df.loc[len(df)] = instance_results

results_filename = f"results_{problem}_LS.csv"
results_filepath = f"{results_folder}/{results_filename}"
df.to_csv(results_filepath, index=False)