function call_local_search_procedure_maragal(A_path, r, m, n, funcName, savePath, hatA_flag)
    % Load file
    data = load(A_path);

    % Extract matrix from nested struct
    if isfield(data, 'Problem')
        Problem = data.Problem;

        if isfield(Problem, 'A')
            A = Problem.A;
        else
            error('Field A not found inside Problem');
        end
    else
        error('Field Problem not found in .mat file');
    end

    % Convert sparse → dense if needed
    if issparse(A)
        A = full(A);
    end

    if (hatA_flag == 1)
        A = A' * A;
        [m, n] = size(A);
    end

    addpath('src/methods/local_search/nsub');
    [R,C,time] = nsub(A,r);
    
    % Calls an appropriate function with given data parameters
    switch funcName
        case 'LSFI_Det'
            addpath('local_search/P1');
            [normLSFI_Det, timeLSFI_Det, swaps_LSFI_Det, R_LSFI_Det, C_LSFI_Det] = LSFI_Det(A, r, m, n, R, C);
            A_inv = inv(A(R_LSFI_Det, C_LSFI_Det));
            H = zeros(n, m);
            H(C_LSFI_Det, R_LSFI_Det) = A_inv;
        case 'LSFI_Det_Symmetric'
            addpath('local_search/P1sym');
            [norm, time, swaps, R] = LSFI_Det_Symmetric(A, r, m, R);
            A_hat = A(R, R);
            H = zeros(n, m);
            H(R, R) = inv(A_hat);
        case 'LSFI_Det_P3'
            addpath('local_search/P13');
            [norm, time, swaps, C_out] = LSFI_Det_P3(A, r, n, R, C);
            A_hat = A(:, C_out);
            H_hat = inv(A_hat' * A_hat) * A_hat';
            H = zeros(n, m);
            H(C_out, :) = H_hat;
        otherwise
            error('Function not recognized. Use "LSFI_Det", "LSFI_Det_Symmetric" or "LSFI_Det_P3".');
    end
    
    % Saves resulting matrix in the specified path
    save(savePath, 'H');
end