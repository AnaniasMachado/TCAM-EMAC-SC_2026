using CSV
using DataFrames
using MAT
using Base.GC
using Statistics

include("../types.jl")
include("../utility.jl")

exp = "table_6"
matrices_folder = "./instances/$(exp)"
results_folder = "results/$(exp)"
hatA_flag = 1

df = DataFrame()

for idx in 1:5
    mat_file = "Maragal_$(idx).mat"

    mat_path = joinpath(matrices_folder, mat_file)
    mat_data = matread(mat_path)

    println("Solving for matrix: $mat_path")
    
    problem_data = mat_data["Problem"]
    A = Matrix(problem_data["A"])

    if hatA_flag == 1
        A = A' * A
    end

    AMP = pinv(A)

    m, n = size(A)
    r = calculate_rank(A)
    d = matrix_norm_0(A) / (m * n)

    A_norm_0 = matrix_norm_0(A)
    AMP_norm_0 = matrix_norm_0(AMP)
    AMP_norm_1 = norm(AMP, 1)

    data = DataInst(A, m, n, r, AMP=AMP)

    H_norm_0 = -1.0
    H_norm_1 = -1.0
    time = -1.0

    GC.gc()

    result = DataFrame(
        m = [m],
        n = [n],
        r = [r],
        d = [d],
        A_norm_0 = [A_norm_0],
        AMP_norm_0 = [AMP_norm_0],
        AMP_norm_1 = [AMP_norm_1]
    )

    append!(df, result)

    GC.gc()
end

results_filename = "instance_stats_$(exp).csv"
results_filepath = joinpath(results_folder, results_filename)
CSV.write(results_filepath, df)