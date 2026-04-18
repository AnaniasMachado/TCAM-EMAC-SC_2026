% exp_table_6.m

clc;
clear;

current_folder = fileparts(mfilename('fullpath'));
methods_path = fullfile(current_folder, '..', 'methods');
addpath(methods_path);

problems = ["P123", "P1Sym"];
problem = problems(2);

column_names = ["m", "n", "r", "d", "||A||_0", "||A^+||_0", "||A^+||_1", "||H||_0", "||H||_1", "time"];

exp = "table_6";
time_limit = 86400;
tol = 1e-5;

matrices_folder = "../../instances/" + exp;
results_folder = "../../results/" + exp;
solutions_folder = "../../solutions/" + exp + "/LS/";

methods = containers.Map;
methods("P123") = "LSLAFI_Det_P3";
methods("P1") = "LSLAFI_Det";
methods("P1Sym") = "LSLAFI_Det_Symmetric";

method = methods(problem);

% Initialize results table
df = table('Size',[0 length(column_names)], ...
           'VariableTypes', repmat("double",1,length(column_names)), ...
           'VariableNames', column_names);

min_unsolvable_m = inf;

for i = 1:5

    A_norm_0 = 0;
    AMP_norm_0 = 0;
    AMP_norm_1 = 0;
    H_norm_0 = 0;
    H_norm_1 = 0;
    runtime = 0;

    instance = matrices_folder + "/Maragal_" + i;
    fprintf("Solving for matrix: %s\n", instance);

    data = load(instance);
    A = data.Problem.A;

    if issparse(A)
        A = full(A);
    end

    [m, n] = size(A);
    r = rank(A, tol);

    d = sum(abs(A(:)) > tol) / (m * n);

    hatA_flag = 1;

    if m < min_unsolvable_m

        if hatA_flag == 1
            A = A' * A;
            r = rank(A, tol);
        end

        A_MP = pinv(A);

        savePath = solutions_folder + "/H_Maragal_" + i + ".mat";

        % ===== TIME MEASUREMENT =====
        t_start = tic;
        success = ls_maragal(instance, r, m, n, method, savePath, hatA_flag, time_limit);
        runtime = toc(t_start);
        % ============================

        if success ~= false

            sol = load(savePath);
            H = sol.H;

            A_norm_0 = sum(abs(A(:)) > tol);
            AMP_norm_0 = sum(abs(A_MP(:)) > tol);
            AMP_norm_1 = sum(sum(abs(A_MP)));
            H_norm_0 = sum(abs(H(:)) > tol);
            H_norm_1 = sum(sum(abs(H)));

        else
            fprintf("TimeLimit: Could not solve instance within the time limit.\n");
            min_unsolvable_m = m;
        end
    end

    % Store results
    new_row = {m, n, r, d, A_norm_0, AMP_norm_0, AMP_norm_1, H_norm_0, H_norm_1, runtime};
    df = [df; new_row];

end

% Save CSV
results_filename = "results_" + problem + "_LS.csv";
results_filepath = results_folder + "/" + results_filename;

writetable(df, results_filepath);

fprintf("Results saved to %s\n", results_filepath);