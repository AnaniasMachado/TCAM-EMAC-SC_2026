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
problems = ["P123"]
problem = problems[1]
epsilon = 10^(-5)
eps_abs = epsilon
eps_rel = 10^(-3)
fixed_tol = false
eps_opt = epsilon
time_limit = 86400

# Gurobi parameters
constraints_set = [["P1", "Sym"], ["P1Sym"]]
constraints = constraints_set[2]

# DRS parameters
lambda = 10^(-2)

stop_crits = ["Opt", "Fixed_Point"]
stop_crit = stop_crits[2]

exp = "table_5"
matrices_folder = "./instances/$(exp)"

results_folder = "results/$(exp)"

solutions_folder = "./solutions/$(exp)"

df = DataFrame()

min_unsolvable_m = Inf

for idx in 1:5
    mat_file = "Maragal_$(idx).mat"

    mat_path = joinpath(matrices_folder, mat_file)
    mat_data = matread(mat_path)

    println("Solving for matrix: $mat_path")

    # println(keys(mat_data))
    
    problem_data = mat_data["Problem"]
    A = Matrix(problem_data["A"])
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

    GC.gc()

    result = DataFrame(
        m = [m],
        n = [n],
        r = [r],
        d = [d],
        A_norm_0 = [A_norm_0],
        AMP_norm_0 = [AMP_norm_0],
        AMP_norm_1 = [AMP_norm_1],
        H_norm_0 = [H_norm_0],
        H_norm_1 = [H_norm_1],
        time = [time]
    )

    append!(df, result)

    GC.gc()
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