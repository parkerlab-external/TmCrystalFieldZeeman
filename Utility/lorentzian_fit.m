function [params, res, paramErr] = lorentzian_fit(wn_data, scan_data, num_curves, options)
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

    opts = optimoptions('lsqnonlin', "Display","None");

    obj_fun = @(params) multi_lorentz_fun(params, wn_data, num_curves) - scan_data;
    if ~isempty(options.AmpRatio)
        obj_fun = @(params) scaled_multi_lorentz_fun(params, wn_data, num_curves, options.AmpRatio) - scan_data;
    end
    
    params = [];
    res = Inf;
    
    lb = [min(wn_data) mean(diff(sort(wn_data))) 0];
    ub = [max(wn_data) (max(wn_data)-min(wn_data)) Inf];
    %If shape is forced keep amp and width constrained
    if options.ForceShape
        lb = [repmat(lb(1), 1, num_curves) lb(2:end)];
        ub = [repmat(ub(1), 1, num_curves) ub(2:end)];
    else
        lb = repelem(lb, num_curves);
        ub = repelem(ub, num_curves);
    end

    %Allow for fixed positions
    if ~isempty(options.FixedPositions)
        num_fixed = numel(options.FixedPositions);
        lb(1:num_fixed) = options.FixedPositions;
        ub(1:num_fixed) = options.FixedPositions;
    end


    if options.FitBackground
        lb = [lb 0];
        ub = [ub Inf];
    end
        
    for i = 1:options.NumRandSamples
        if options.ForceShape
            paramsGuess = [min(wn_data) + rand(1,num_curves)*(max(wn_data) - min(wn_data))...
                (max(wn_data)-min(wn_data))/5 ...
                mean(scan_data)];
        else
            paramsGuess = [min(wn_data) + rand(1,num_curves)*(max(wn_data) - min(wn_data))...
                repmat((max(wn_data)-min(wn_data))/5, 1, num_curves) ...
                repmat(mean(scan_data), 1, num_curves)];
        end
        if options.FitBackground
            paramsGuess = [paramsGuess min(scan_data)];
        end
        [testparams, testres, ~, ~, ~, ~, J] = lsqnonlin(obj_fun, paramsGuess, lb, ub, opts);
        if testres < res
            res = testres;
            params = testparams;
            if isfield(options, "scan_err")
                J = J(:,1:num_curves); %Position subspace only
                T = inv(J.'*J);
                paramErr = T*(J.')*spdiags(options.scan_err(:).^2, 0 , numel(scan_data), numel(scan_data))*J*T;
                paramErr = full(sqrt(diag(paramErr)));
            end
        end
    end
end
