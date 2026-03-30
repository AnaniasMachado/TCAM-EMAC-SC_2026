function success = ls_maragal(A_path, r, m, n, funcName, savePath, hatA_flag, time_limit)

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

    addpath('src/methods/local_search/nsub');
    [R,C,time] = nsub(A,r);

    while true

        switch funcName

            %% =====================================
            %% LSFI_Det (inline)
            %% =====================================
            case 'LSFI_Det'

                swaps = 0;

                Rb = [];
                Cb = [];

                for i = (1:m)
                    m1 = find(R==i);
                    if isempty(m1)==1
                        Rb = [Rb;i];
                    end
                end

                Ar = A(R,C);

                for j = (1:n)
                    n1 = find(C==j);
                    if isempty(n1)==1
                        Cb = [Cb;j];
                    end
                end

                flag = 1;

                while flag > 0  && swaps <= 1000

                    % TIME CHECK
                    if toc(t0) > time_limit
                        success = false;
                        return;
                    end

                    flag = 0;

                    %% FOR ROWS
                    [L1,U1,P1] = lu(Ar');

                    for i = (1:m-r)

                        b = P1 * A(Rb(i),C)';
                        y = L1\b;
                        alfa = U1\y;

                        local_alfa = find(abs(alfa)>1.000001);
                        T = isempty(local_alfa);

                        if T==0

                            local_alfa = local_alfa(1);

                            swaps = swaps + 1;

                            el_save = R(local_alfa);
                            R(local_alfa) = Rb(i);
                            Rb(i) = el_save;

                            flag = flag + 1;

                            Ar = A(R,C);
                            [L1,U1,P1] = lu(Ar');
                        end
                    end

                    %% FOR COLUMNS
                    [L2,U2,P2] = lu(Ar);

                    for i = (1:n-r)

                        b = P2 * A(R,Cb(i));
                        y = L2\b;
                        alfa = U2\y;

                        local_alfa = find(abs(alfa)>1.000001);
                        T = isempty(local_alfa);

                        if T==0

                            local_alfa = local_alfa(1);

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
                H = zeros(n, m);
                H(C, R) = A_inv;


            %% =====================================
            %% LSFI_Det_Symmetric (inline)
            %% =====================================
            case 'LSFI_Det_Symmetric'

                swaps = 0;

                Rb = [];

                for i = (1:m)
                    m1 = find(R==i);
                    if isempty(m1)==1
                        Rb = [Rb;i];
                    end
                end

                flag = 1;
                det_ref = abs(det(A(R,R)));

                while flag > 0

                    % TIME CHECK
                    if toc(t0) > time_limit
                        success = false;
                        return;
                    end

                    flag = 0;

                    for i = (1:r)

                        R_aux = R;

                        for j = (1:m-r)

                            R_aux(i) = Rb(j);

                            value_det = abs(det(A(R_aux,R_aux)));

                            if value_det > det_ref

                                det_ref = value_det;

                                swaps = swaps + 1;

                                el_save = R(i);
                                R(i) = Rb(j);
                                Rb(j) = el_save;

                                flag = flag + 1;

                                break

                            else
                                R_aux(i) = R(i);
                            end
                        end
                    end
                end

                A_hat = A(R,R);
                H = zeros(n, m);
                H(R, R) = inv(A_hat);


            %% =====================================
            %% LSFI_Det_P3 (inline)
            %% =====================================
            case 'LSFI_Det_P3'

                swaps = 0;

                Cb = [];

                Ar = A(R,C);

                for j = (1:n)
                    n1 = find(C==j);
                    if isempty(n1)==1
                        Cb = [Cb;j];
                    end
                end

                flag = 1;

                while flag > 0

                    % TIME CHECK
                    if toc(t0) > time_limit
                        success = false;
                        return;
                    end

                    flag = 0;

                    Ar = A(R,C);

                    [L2,U2,P2] = lu(Ar);

                    for i = (1:n-r)

                        b = P2 * A(R,Cb(i));
                        y = L2\b;
                        alfa = U2\y;

                        local_alfa = find(abs(alfa)>1);
                        T = isempty(local_alfa);

                        if T==0

                            local_alfa = local_alfa(1);

                            swaps = swaps + 1;

                            el_save = C(local_alfa);
                            C(local_alfa) = Cb(i);
                            Cb(i) = el_save;

                            Ar = A(R,C);
                            [L2,U2,P2] = lu(Ar);

                            flag = flag+1;
                        end
                    end
                end

                A_hat =  A(:,C);
                H_hat = inv(A_hat' * A_hat) * A_hat';

                H = zeros(n, m);
                H(C, :) = H_hat;

            otherwise
                error('Function not recognized.');
        end

        break;
    end
end