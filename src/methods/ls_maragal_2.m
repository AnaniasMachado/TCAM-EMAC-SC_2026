function success = ls_maragal_2(A_path, r, m, n, funcName, savePath, hatA_flag, time_limit)

    success = true;

    %% Load file
    data = load(A_path);

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

    if issparse(A)
        A = full(A);
    end

    if (hatA_flag == 1)
        A = A' * A;
        [m, n] = size(A);
    end

    t0 = tic;

    % addpath('src/methods/local_search/nsub');
    [R,C,time] = nsub(A,r);

    while true

        switch funcName

            %% =====================================
            %% LSLAFI_Det (inline)
            %% =====================================
            case 'LSLAFI_Det'

                [normLSLAFI_Det,timeLSLAFI_Det,swaps_LSLAFI_Det,R_LSLAFI_Det,C_LSLAFI_Det] = LSLAFI_Det(A,r,m,n,R,C)

                % TIME CHECK
                if toc(t0) > time_limit
                    success = false;
                    return;
                end

                A_inv = inv(A(R,C));
                H = zeros(n,m);
                H(C,R) = A_inv;

            %% =====================================
            %% LSLAFI_Det_Symmetric (inline)
            %% =====================================
            case 'LSLAFI_Det_Symmetric'

                [norm,time,swaps,R,H] = LSLAFI_Det_Symmetric (A,r,m,R)

                % TIME CHECK
                if toc(t0) > time_limit
                    success = false;
                    return;
                end

                A_hat = A(R,R);
                H = zeros(n,m);
                H(R,R) = inv(A_hat);

            %% =====================================
            %% LSLAFI_Det_P3 (inline)
            %% =====================================
            case 'LSLAFI_Det_P3'

                [norm,time,swaps,C_out] = LSLAFI_Det_P3 (A,r,n,R,C)

                % TIME CHECK
                if toc(t0) > time_limit
                    success = false;
                    return;
                end

                A_hat = A(:, C);
                H_hat = (A_hat' * A_hat) \ A_hat';
                H = zeros(n,m);
                H(C,:) = H_hat;

            otherwise
                error('Function not recognized.');
        end

        % Saves resulting matrix in the specified path
        save(savePath, 'H');

        break;
    end
end