function success = call_local_search_procedure_maragal( ...
    A_path, r, m, n, funcName, savePath, hatA_flag, timelimit)

    success = true;

    % Ensure output directory exists
    [folder, ~, ~] = fileparts(savePath);
    if ~isempty(folder) && ~exist(folder, 'dir')
        mkdir(folder);
    end

    % Start pool if needed
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool;
    end

    % Launch async job
    f = parfeval(pool, ...
        @call_local_search_procedure_maragal_inner, ...
        0, ...
        A_path, r, m, n, funcName, savePath, hatA_flag);

    % -------------------------
    % TIME LIMIT CONTROL LOOP
    % -------------------------
    t0 = tic;

    while true
        pause(30);  % avoid busy waiting

        % Case 1: finished
        if strcmp(f.State, 'finished')
            try
                fetchOutputs(f);  % propagate errors
                if ~isempty(f.Error)
                    success = false;
                    warning("Worker error: %s", f.Error.message);
                end
            catch ME
                success = false;
                warning("Execution failed: %s", ME.message);
            end
            return;
        end

        % Case 2: timeout
        if toc(t0) > timelimit
            cancel(f);
            success = false;
            warning('Time limit exceeded. Task cancelled.');
            return;
        end
    end
end