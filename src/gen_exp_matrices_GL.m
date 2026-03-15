d_values = [0.1, 0.25];
n_mtx = 5;
output_dir = "../instances/table_3";

for i = 1:5
    m = 1000;
    n = 100 * i;
    r_values = n / 4;
    gen_multiple_matrices_GL(m, n, r_values, d_values, n_mtx, output_dir);
end