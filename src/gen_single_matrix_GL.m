function A = gen_single_matrix_GL(m, n, r, density)
    B = sprand(m, r, density);
    k = n - r;
    if k > 0
        W = randn(r, k);
        C = B * W;
        A = [B, C];
    else
        A = B;
    end
    A = full(A);

end