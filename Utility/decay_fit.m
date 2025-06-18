function [amps,tcs,bg,paramVar] = decay_fit(decay_axis,decay_data,num_decays,options)
%DECAY_FIT Fit a single or multiple exponential decays to 1 or 2D input data
%   Detailed explanation goes here
arguments
    decay_axis
    decay_data
    num_decays = 1;
    options.FixedBackground
    options.tc_lb
    options.tc_ub
    options.decay_errors
end
%Try to match up the dimensions of the decay_axis and the decay_data
if size(decay_axis,2) == 1
    decay_axis = decay_axis';
end
if size(decay_axis,2) ~= size(decay_data,2)
    decay_data = decay_data';
end
if size(decay_axis,2) == size(decay_data,1)
    warning("Decay axis and num_lines same size, check that axis is correct");
end


num_lines = size(decay_data,1);
tcGuess = decay_axis(round(end/3));
paramsGuess = [repmat(decay_data(:,1),[num_decays,1])/num_decays;
    repmat(tcGuess, [num_decays, 1]);
    decay_data(1,end);
    ];

%Upper Bounds
ub = [repmat(decay_data(:,1),[num_decays,1])*2;
    repmat(5*decay_axis(end), [num_decays, 1]);
    Inf;
    ];

%Lower Bounds
lb = [zeros(num_lines*num_decays,1);
    mean(abs(diff(decay_axis)))*ones(num_decays, 1);
    0;
    ];

if isfield(options, "FixedBackground")
    lb(end) = options.FixedBackground;
    ub(end) = options.FixedBackground;
end


opts = optimset('Display','off');
[fitParams,resnorm,residual,~,~,~,J] = lsqcurvefit(@exp_fit_fun, paramsGuess, decay_axis, decay_data, lb, ub, opts);

amps = fitParams(1:num_decays*num_lines);
tcs = fitParams(num_decays*num_lines + (1:num_decays));
bg = fitParams(end);

%Testing covariance estimates
if isfield(options, "decay_errors")
    %Here, covParams is actually the standard devation and not the variance
    covParams = zeros(numel(fitParams), 4);
    T = inv(J.'*J);

%     %Linear estimate using standard deviations
%     covParams(:,1) = T*(J.')*options.decay_errors(:);
%     %Using variance
    temp = T*(J.')*spdiags(options.decay_errors(:).^2, 0 , numel(residual(:)), numel(residual(:)))*J*T;
    covParams(:,2) = sqrt(diag(temp));
    
%     %Using residuals estimate
%     temp = T*resnorm/(numel(residual) - numel(fitParams));
%     covParams(:,3) = sqrt(diag(temp));
%     %Linear version
%     covParams(:,4) = T*(J.')*abs(residual(:));
    paramVar = covParams(:,2);
else
    covParams = inv(J'*J)*resnorm/(numel(residual) - numel(fitParams));
    paramVar = full(diag(covParams));
end

    function output = exp_fit_fun(params, decay_axis)
        param_amps = reshape(params(1:num_lines*num_decays),num_lines,num_decays);
        param_tcs = params(num_lines*num_decays+(1:num_decays));
        param_bg = params(end);
        output = zeros(num_lines, numel(decay_axis));
        for i = 1:num_decays
            output = output + param_amps(:,i)*exp(-decay_axis/param_tcs(i));
        end
        output = output + param_bg;
    end
end

