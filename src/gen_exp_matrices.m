m_val = 20;
d_values = 1.0;
n_mtx = 5;
output_dir = "../instances/table_1";

for i = 1:4
    m = m_val * i;
    n = m;
    r_values = m / 4;
    gen_multiple_matrices(m, n, r_values, d_values, n_mtx, output_dir);
end