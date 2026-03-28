import matlab.engine
import os

from ..utility import read_instance, read_maragal_instance, get_data_from_instance, matrix_rank, matrix_vec_0_norm

def local_search_procedure(exp, solutions_folder, problem, instance, method, hatA_flag):
    if method not in ["LSFI_Det", "LSFI_Det_Symmetric", "LSFI_Det_P3"]:
        raise Exception("Error: Function does not exist.")

    eng = matlab.engine.start_matlab()

    eng.addpath(eng.genpath("."), nargout=0)

    eng.cd(".", nargout=0)

    A = read_instance(instance, "A")

    m, n, r, d, idx = get_data_from_instance(instance)
    d = int(d * 100)

    save_path = f"{solutions_folder}/LS/problem_{problem}_m_{m}_n_{n}_d_{d}_idx_{idx}"

    eng.call_local_search_procedure(instance, r, m, n, method, save_path, hatA_flag, nargout=0)

    H = read_instance(save_path, "H")

    eng.quit()

    return H

def local_search_procedure_maragal(exp, solutions_folder, problem, instance, method, hatA_flag, time_limit):
    if method not in ["LSFI_Det", "LSFI_Det_Symmetric", "LSFI_Det_P3"]:
        raise Exception("Error: Function does not exist.")

    eng = matlab.engine.start_matlab()

    eng.addpath(eng.genpath("."), nargout=0)

    eng.cd(".", nargout=0)

    A = read_maragal_instance(instance)

    m, n = A.shape
    r = matrix_rank(A)
    d = matrix_vec_0_norm(A) / (m * n)
    d = int(d * 100)

    save_path = f"{solutions_folder}/LS/problem_{problem}_m_{m}_n_{n}_d_{d}"

    success = eng.call_local_search_procedure_maragal(instance, r, m, n, method, save_path, hatA_flag, float(time_limit), nargout=1)

    eng.quit()

    if success:
        H = read_instance(save_path, "H")
        return H
    else:
        return 0