function success = ls_general(A_path, r, m, n, funcName, savePath, hatA_flag, time_limit)

    success = true;

    %% Load file
    A = load(A_path);

    if isfield(A, 'A')
        A = A.A;
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
    % [R,C,time] = nsub(A,r);

    [~, R_qr, Erow] = qr(A', 'vector');
    r = sum(abs(diag(R_qr)) > 1e-5);
    R = sort(Erow(1:r));
    [~, ~, Ecol] = qr(A(R,:), 'vector');
    C = sort(Ecol(1:r));
    [max(svd(A(R,C))),min(svd(A(R,C))),length(R),length(C),r];
    assert(rank(A(R,C), 1e-6) == r, 'A(R,C) is rank-deficient');

    switch funcName

        %% =====================================
        %% LSLAFI_Det (inline)
        %% =====================================
        case 'LSLAFI_Det'

            swaps = 0;
            Rb = []; Cb = [];

            for i = 1:m
                if isempty(find(R==i,1))
                    Rb = [Rb;i];
                end
            end
            for j = 1:n
                if isempty(find(C==j,1))
                    Cb = [Cb;j];
                end
            end

            Ar = A(R,C);
            flag = 1;
            eps = 1e-6;

            while flag > 0 && swaps <= 1000

                % TIME CHECK
                if toc(t0) > time_limit
                    success = false;
                    return;
                end

                flag = 0;

                % --- ROWS ---
                [L1,U1,P1] = lu(Ar');
                for i = 1:m-r
                    % TIME CHECK
                    if toc(t0) > time_limit
                        success = false;
                        return;
                    end

                    b = P1 * A(Rb(i), C)';
                    y = L1\b;
                    alfa = U1\y;

                    [biggest_alfa, local_alfa] = max(abs(alfa));
                    if abs(biggest_alfa) > (1+eps)
                        swaps = swaps + 1;
                        el_save = R(local_alfa);
                        R(local_alfa) = Rb(i);
                        Rb(i) = el_save;
                        flag = flag + 1;
                        Ar = A(R,C);
                        [L1,U1,P1] = lu(Ar');
                    end
                end

                % --- COLUMNS ---
                [L2,U2,P2] = lu(Ar);
                for i = 1:n-r
                    % TIME CHECK
                    if toc(t0) > time_limit
                        success = false;
                        return;
                    end

                    b = P2 * A(R, Cb(i));
                    y = L2\b;
                    alfa = U2\y;

                    [biggest_alfa, local_alfa] = max(abs(alfa));
                    if abs(biggest_alfa) > (1+eps)
                        swaps = swaps + 1;
                        el_save = C(local_alfa);
                        C(local_alfa) = Cb(i);
                        Cb(i) = el_save;
                        flag = flag + 1;
                        Ar = A(R,C);
                        [L2,U2,P2] = lu(Ar);
                    end
                end
            end

            A_inv = inv(A(R,C));
            H = zeros(n,m);
            H(C,R) = A_inv;

        %% =====================================
        %% LSLAFI_Det_Symmetric (inline)
        %% =====================================
        case 'LSLAFI_Det_Symmetric'

            swaps = 0;
            Rb = [];

            for i = 1:m
                if isempty(find(R==i,1))
                    Rb = [Rb;i];
                end
            end

            flag = 1;
            det_ref = abs(det(A(R,R)));
            eps = 1e-6;

            while flag > 0

                % TIME CHECK
                if toc(t0) > time_limit
                    success = false;
                    return;
                end

                flag = 0;
                flag_change = 0;

                for i = 1:r
                    % TIME CHECK
                    if toc(t0) > time_limit
                        success = false;
                        return;
                    end

                    R_aux = R;

                    for j = 1:m-r
                        % TIME CHECK
                        if toc(t0) > time_limit
                            success = false;
                            return;
                        end

                        R_aux(i) = Rb(j);
                        value_det = abs(det(A(R_aux,R_aux)));
                        if value_det > (det_ref+eps)
                            det_ref = value_det;
                            j_aux = j;
                            flag_change = 1;
                        end
                        R_aux(i) = R(i);
                    end

                    if flag_change == 1
                        swaps = swaps + 1;
                        el_save = R(i);
                        R(i) = Rb(j_aux);
                        Rb(j_aux) = el_save;
                        flag = flag + 1;
                        flag_change = 0;
                    end
                end
            end

            A_hat = A(R,R);
            H = zeros(n,m);
            H(R,R) = inv(A_hat);

        %% =====================================
        %% LSLAFI_Det_P3 (inline)
        %% =====================================
        case 'LSLAFI_Det_P3'

            swaps = 0;
            Cb = [];

            for j = 1:n
                if isempty(find(C==j,1))
                    Cb = [Cb;j];
                end
            end

            Ar = A(R,C);
            flag = 1;
            eps = 1e-6;

            while flag > 0
                % TIME CHECK
                if toc(t0) > time_limit
                    success = false;
                    return;
                end

                flag = 0;
                Ar = A(R,C);

                [L2,U2,P2] = lu(Ar);
                for i = 1:n-r
                    % TIME CHECK
                    if toc(t0) > time_limit
                        success = false;
                        return;
                    end

                    b = P2 * A(R,Cb(i));
                    y = L2\b;
                    alfa = U2\y;

                    [biggest_alfa, local_alfa] = max(abs(alfa));
                    if abs(biggest_alfa) > (1+eps)
                        swaps = swaps + 1;
                        el_save = C(local_alfa);
                        C(local_alfa) = Cb(i);
                        Cb(i) = el_save;
                        Ar = A(R,C);
                        [L2,U2,P2] = lu(Ar);
                        flag = flag + 1;
                    end
                end
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
end