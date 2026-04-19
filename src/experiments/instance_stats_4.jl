using CSV
using DataFrames
using MAT
using Base.GC
using Statistics

include("../types.jl")
include("../utility.jl")

exp = "table_4"
matrices_folder = "./instances/$(exp)"
m = 1000
n_values = [100 * i for i in 1:5]
d_values = [0.1, 0.25]

results_folder = "results/$(exp)"

df = DataFrame()

A_norm_0_list = []
AMP_norm_0_list = []
AMP_norm_1_list = []

max_idx = 5

for n in n_values
    for d_val in d_values
        d = Int(100 * d_val)
        for idx in 1:max_idx
            r = Int(3 * n / 4)
            mat_file = "A_m$(m)_n$(n)_r$(r)_d$(d)_idx$(idx).mat"

            mat_path = joinpath(matrices_folder, mat_file)
            mat_data = matread(mat_path)

            println("Solving for matrix: $mat_path")
            
            A = mat_data["A"]
            A = Matrix(A)
            A = A' * A
            AMP = pinv(A)

            A_norm_0 = matrix_norm_0(A)
            AMP_norm_0 = matrix_norm_0(AMP)
            AMP_norm_1 = norm(AMP, 1)

            push!(A_norm_0_list, A_norm_0)
            push!(AMP_norm_0_list, AMP_norm_0)
            push!(AMP_norm_1_list, AMP_norm_1)

            GC.gc()

            if idx == max_idx
                A_norm_0_mean = -1.0
                AMP_norm_0_mean = -1.0
                AMP_norm_1_mean = -1.0

                if !(-1.0 in A_norm_0_list)
                    A_norm_0_mean = mean(A_norm_0_list)
                    AMP_norm_0_mean = mean(AMP_norm_0_list)
                    AMP_norm_1_mean = mean(AMP_norm_1_list)
                end

                result = DataFrame(
                    m = [m],
                    n = [n],
                    r = [r],
                    d = [d],
                    A_norm_0 = [A_norm_0_mean],
                    AMP_norm_0 = [AMP_norm_0_mean],
                    AMP_norm_1 = [AMP_norm_1_mean]
                )

                append!(df, result)

                empty!(A_norm_0_list)
                empty!(AMP_norm_0_list)
                empty!(AMP_norm_1_list)

                GC.gc()
            end
        end
    end
end

results_filename = "instance_stats_$(exp).csv"
results_filepath = joinpath(results_folder, results_filename)
CSV.write(results_filepath, df)