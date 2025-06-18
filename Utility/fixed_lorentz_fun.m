function total = fixed_lorentz_fun(scan_wns, fixed_positions, params)
    %params is [amps x num_fixed offset width bg]
    total = 0;    
    for i = 1:numel(fixed_positions)
        new_params = [fixed_positions(i)+params(end-2), params(end-1), params(i)];
        total = total + multi_lorentz_fun(new_params, scan_wns, 1);
    end
    total = total + params(end);
end