import time
import os
import pandas as pd
from scipy.io import savemat

from ..utility import *
from ..methods.local_search import *

problems = ["P123", "P1Sym"]
problem = problems[0]

column_names = ["m", "n", "r", "d", "||A||_0", "||A^+||_0", "||A^+||_1", "||H||_0" , "||H||_1", "time"]

exp = "table_3"

matrices_folder = f"./instances/{exp}"
results_folder = f"./results/{exp}"
solutions_folder = f"./solutions/{exp}"

methods = {
    "P123": "LSFI_Det_P3",
    "P1" : "LSFI_Det",
    "P1Sym" : "LSFI_Det_Symmetric"
}
method = methods[problem]

df = pd.DataFrame(columns=column_names)

instances = get_instances(matrices_folder)

A_norm_0_list = []
AMP_norm_0_list = []
AMP_norm_1_list = []
H_norm_0_list = []
H_norm_1_list = []
time_list = []

count = 0
for instance in instances:
    print(f"Solving for matrix: {instance}")
    count += 1
    A = read_instance(instance, "A")
    m, n, r, d, idx = get_data_from_instance(instance)
    d = int(d * 100)
    hatA_flag = 0

    if hatA_flag:
        A = A.T @ A

    A_MP = np.linalg.pinv(A)

    start_time = time.time()
    H = local_search_procedure(exp, solutions_folder, problem, instance, method, hatA_flag)
    end_time = time.time()

    A_norm_0 = matrix_vec_0_norm(A)
    AMP_norm_0 = matrix_vec_0_norm(A_MP)
    AMP_norm_1 = matrix_vec_1_norm(A_MP)
    H_norm_0 = matrix_vec_0_norm(H)
    H_norm_1 = matrix_vec_1_norm(H)
    runtime = end_time - start_time

    A_norm_0_list.append(A_norm_0)
    AMP_norm_0_list.append(AMP_norm_0)
    AMP_norm_1_list.append(AMP_norm_1)
    H_norm_0_list.append(H_norm_0)
    H_norm_1_list.append(H_norm_1)
    time_list.append(runtime)

    if (count % 5 == 0):
        A_norm_0 = sum(A_norm_0_list) / len(A_norm_0_list)
        AMP_norm_0 = sum(AMP_norm_0_list) / len(AMP_norm_0_list)
        AMP_norm_1 = sum(AMP_norm_1_list) / len(AMP_norm_1_list)
        H_norm_0 = sum(H_norm_0_list) / len(H_norm_0_list)
        H_norm_1 = sum(H_norm_1_list) / len(H_norm_1_list)
        runtime = sum(time_list) / len(time_list)

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

        A_norm_0_list = []
        AMP_norm_0_list = []
        AMP_norm_1_list = []
        H_norm_0_list = []
        H_norm_1_list = []
        time_list = []

results_filename = f"results_{problem}_LS.csv"
results_filepath = f"{results_folder}/{results_filename}"
df.to_csv(results_filepath, index=False)