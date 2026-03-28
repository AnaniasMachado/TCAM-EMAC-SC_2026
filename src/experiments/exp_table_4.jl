using CSV
using DataFrames
using MAT
using Base.GC
using Statistics

include("../types.jl")
include("../utility.jl")
include("../methods/solvers.jl")
include("../methods/solvers_cal.jl")
include("../methods/drs.jl")

methods = ["Gurobi", "Gurobi_Cal", "DRS"]
method = methods[3]

# Mixed parameters
problems = ["P1Sym"]
problem = problems[1]
epsilon = 10^(-5)
eps_abs = epsilon
eps_rel = 10^(-3)
fixed_tol = false
eps_opt = epsilon
time_limit = 1200

# Gurobi parameters
constraints_set = [["P1", "Sym"], ["P1Sym"]]
constraints = constraints_set[2]

# DRS parameters
lambda = 10^(-2)

stop_crits = ["Opt", "Fixed_Point"]
stop_crit = stop_crits[2]

exp = "table_4"
matrices_folder = "./instances/$(exp)"
m = 1000
n_values = [100 * i for i in 1:5]
d_values = [0.1, 0.25]

results_folder = "results/$(exp)"

solutions_folder = "./solutions/$(exp)"

df = DataFrame()

A_norm_0_list = []
AMP_norm_0_list = []
AMP_norm_1_list = []
H_norm_0_list = []
H_norm_1_list = []
time_list = []

max_idx = 5
min_unsolvable_m = Inf

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
            AMP_norm_0 = matrix_norm_0(A)
            AMP_norm_1 = norm(AMP, 1)

            data = DataInst(A, m, n, r, AMP=AMP)

            H_norm_0 = -1.0
            H_norm_1 = -1.0
            time = -1.0
            if (m < min_unsolvable_m)
                if method == "Gurobi"
                    try
                        time = @elapsed begin
                            H = gurobi_solver(data, constraints, eps_opt, time_limit)
                        end
                        H_norm_0 = matrix_norm_0(H)
                        H_norm_1 = norm(H, 1)
                        H_rank = calculate_rank(H)

                        problem_label = join(constraints, "_")

                        solution_filename = "Gurobi/problem_$(problem_label)_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                        solution_filepath = joinpath(solutions_folder, solution_filename)
                        matwrite(solution_filepath, Dict("H" => H, "time" => time))
                    catch e
                        if isa(e, ErrorException)
                            global min_unsolvable_m = min(m, min_unsolvable_m)
                        else
                            throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
                        end
                    end
                elseif method == "Gurobi_Cal"
                    try
                        time = @elapsed begin
                            H = gurobi_solver_cal(data, problem, eps_opt, time_limit)
                        end
                        H_norm_0 = matrix_norm_0(H)
                        H_norm_1 = norm(H, 1)
                        H_rank = calculate_rank(H)

                        solution_filename = "Gurobi_Cal/problem_$(problem)_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                        solution_filepath = joinpath(solutions_folder, solution_filename)
                        matwrite(solution_filepath, Dict("H" => H, "time" => time))
                    catch e
                        if isa(e, ErrorException)
                            global min_unsolvable_m = min(m, min_unsolvable_m)
                        else
                            throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
                        end
                    end
                elseif method == "DRS"
                    time = @elapsed begin
                        H, k = drs(A, lambda, eps_abs, eps_rel, problem, fixed_tol, eps_opt, stop_crit, time_limit)
                    end
                    if H == "-"
                        global min_unsolvable_m = min(m, min_unsolvable_m)
                    else
                        H_norm_0 = matrix_norm_0(H)
                        H_norm_1 = norm(H, 1)
                        H_rank = calculate_rank(H)

                        if fixed_tol && stop_crit == "Opt"
                            solution_filename = "DRS_Opt_Eps/problem_$(problem)_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                            solution_filepath = joinpath(solutions_folder, solution_filename)
                            matwrite(solution_filepath, Dict("H" => H, "time" => time, "k" => k))
                        elseif !fixed_tol && stop_crit == "Opt"
                            solution_filename = "DRS_Opt_r0/problem_$(problem)_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                            solution_filepath = joinpath(solutions_folder, solution_filename)
                            matwrite(solution_filepath, Dict("H" => H, "time" => time, "k" => k))
                        elseif fixed_tol && stop_crit == "Fixed_Point"
                            solution_filename = "DRS_FP_Eps/problem_$(problem)_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                            solution_filepath = joinpath(solutions_folder, solution_filename)
                            matwrite(solution_filepath, Dict("H" => H, "time" => time, "k" => k))
                        elseif !fixed_tol && stop_crit == "Fixed_Point"
                            solution_filename = "DRS_FP_r0/problem_$(problem)_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                            solution_filepath = joinpath(solutions_folder, solution_filename)
                            matwrite(solution_filepath, Dict("H" => H, "time" => time, "k" => k))
                        end
                    end
                else
                    throw(ErrorException("Invalid method chose."))
                end
            end

            push!(A_norm_0_list, A_norm_0)
            push!(AMP_norm_0_list, AMP_norm_0)
            push!(AMP_norm_1_list, AMP_norm_1)
            push!(H_norm_0_list, H_norm_0)
            push!(H_norm_1_list, H_norm_1)
            push!(time_list, time)

            GC.gc()

            if idx == max_idx
                A_norm_0_mean = -1.0
                AMP_norm_0_mean = -1.0
                AMP_norm_1_mean = -1.0
                H_norm_0_mean = -1.0
                H_norm_1_mean = -1.0
                time_mean = -1.0

                if !(-1.0 in A_norm_0_list)
                    A_norm_0_mean = mean(A_norm_0_list)
                    AMP_norm_0_mean = mean(AMP_norm_0_list)
                    AMP_norm_1_mean = mean(AMP_norm_1_list)
                    H_norm_0_mean = mean(H_norm_0_list)
                    H_norm_1_mean = mean(H_norm_1_list)
                    time_mean = mean(time_list)
                end

                result = DataFrame(
                    m = [m],
                    n = [n],
                    r = [r],
                    d = [d],
                    A_norm_0 = [A_norm_0_mean],
                    AMP_norm_0 = [AMP_norm_0_mean],
                    AMP_norm_1 = [AMP_norm_1_mean],
                    H_norm_0_mean = [H_norm_0_mean],
                    H_norm_1_mean = [H_norm_1_mean],
                    time_mean = [time_mean]
                )

                append!(df, result)

                empty!(A_norm_0_list)
                empty!(AMP_norm_0_list)
                empty!(AMP_norm_1_list)
                empty!(H_norm_0_list)
                empty!(H_norm_1_list)
                empty!(time_list)

                GC.gc()
            end
        end
    end
end

if method == "Gurobi"
    problem_label = join(constraints, "_")
    results_filename = "results_$(problem_label)_Gurobi.csv"
    results_filepath = joinpath(results_folder, results_filename)
    CSV.write(results_filepath, df)
elseif method == "Gurobi_Cal"
    results_filename = "results_$(problem)_Gurobi_Cal.csv"
    results_filepath = joinpath(results_folder, results_filename)
    CSV.write(results_filepath, df)
elseif method == "DRS"
    if fixed_tol && stop_crit == "Opt"
        results_filename = "results_$(problem)_DRS_Opt_Eps.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    elseif !fixed_tol && stop_crit == "Opt"
        results_filename = "results_$(problem)_DRS_Opt_r0.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    elseif fixed_tol && stop_crit == "Fixed_Point"
        results_filename = "results_$(problem)_DRS_FP_Eps.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    else
        results_filename = "results_$(problem)_DRS_FP_r0.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    end
else
    throw(ErrorException("Invalid method chose."))
end