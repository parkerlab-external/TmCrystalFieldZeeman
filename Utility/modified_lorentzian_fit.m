function [params, res, paramErr] = modified_lorentzian_fit(wn_data, scan_data, num_curves, options)
    arguments
        wn_data
        scan_data
        num_curves = 1
        options.ForceShape = false
        options.FitBackground = true
        options.NumRandSamples = 100
        options.AmpRatio = []
        options.FixedPositions = []
        options.scan_err;
    end

    if ~isequal(size(scan_data), size(wn_data))
        error("I am saving you a HUGE headache by throwing you this error. PLEASE make sure the array sizes are equal (transpose might be necessary)")
    end

    opts = optimoptions('lsqnonlin', "Display","None",'SpecifyObjectiveGradient',false);

    
    params = [];
    res = Inf;
    
    lb = [min(wn_data) mean(diff(sort(wn_data))) 0];
    ub = [max(wn_data) (max(wn_data)-min(wn_data))/3 Inf];

    lb = repelem(lb, num_curves);
    ub = repelem(ub, num_curves);

    if options.FitBackground
        lb = [lb 0];
        ub = [ub Inf];
    end

    

    %Allow for fixed positions
    if ~isempty(options.FixedPositions)
        num_fixed = numel(options.FixedPositions);
        lb(1:num_fixed) = options.FixedPositions - mean(abs(diff(options.FixedPositions)))/3;
        ub(1:num_fixed) = options.FixedPositions + mean(abs(diff(options.FixedPositions)))/3;
    end

    %Constratint Matrix:
    Aeq = zeros(numel(lb), numel(lb));

    %Force position differences to be maintained
    Nc = num_fixed-1;
    Aeq(1:Nc,1:Nc) = eye(Nc);

    %Force widths to match too
    Aeq((num_curves+1):(num_curves+Nc),(num_curves+1):(num_curves+Nc)) = eye(Nc);
    Aeq = Aeq - circshift(Aeq, 1, 2);

    %Since the lower bounds contain this info, find Beq:
    Beq = Aeq*lb';
    
    obj_fun = @(params) multi_lorentz_fun( params , wn_data, num_curves) - scan_data;
    
    for i = 1:options.NumRandSamples
        paramsGuess = lb + 0.5*( rand(size(lb)) + rand(size(lb)) ).*(ub - lb);
        paramsGuess(1:num_fixed) = options.FixedPositions;
        paramsGuess(2*num_curves+1:end) = min(scan_data) + rand(1)*(max(scan_data) - min(scan_data));

        [testparams, testres, ~, ~, ~, ~, J] = lsqnonlin(obj_fun, paramsGuess,lb,ub, [], [], Aeq, Beq, [],opts);
        if testres < res
            res = testres;
            params = testparams;
            if isfield(options, "scan_err")
                J = J(:,1:num_curves);
                T = inv(J.'*J);
                paramErr = abs(T*J.'*diag(options.scan_err(:)));
            end
        end
    end

    % if isfield(options, "scan_err")
    %     [~,J] = multi_lorentz_fun(params, wn_data, num_curves);
    %     J = J(:,1:num_curves);
    %     T = inv(J.'*J);
    %     paramErr = abs(T*J.'*diag(options.scan_err(:)));
    % end
    % function [eval, J] = obj_fun(params)
    %     [eval, J] = multi_lorentz_fun( params , wn_data, num_curves);
    %     eval = eval - scan_data;
    % end
end
