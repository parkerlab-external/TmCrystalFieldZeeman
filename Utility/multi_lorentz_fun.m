function [total, J] = multi_lorentz_fun(params, x, num_curves)
    lorentz_fun = @(params, x) params(3)./( 1 + (( x-params(1) )/params(2) ).^2 );
    num_params = numel(params);
    addBackground = false;
    forceShape = false;
    if num_params == 3 + num_curves || num_params == 3*num_curves + 1
        addBackground = true;
    end
    if num_params == 3 + num_curves || num_params == 2 + num_curves
        forceShape = true;
    end

    
    total = 0;
    for i = 1:num_curves
        if forceShape
            total = total + lorentz_fun(params([i 1+num_curves 2+num_curves]),x);
        else
            total = total + lorentz_fun(params([i i+num_curves i+2*num_curves]),x);
        end
    end
    if addBackground
        total = total + params(end);
    end

    %Compute Jacobian, not for use with ForceShape
    if nargout > 1
        dpos = @(x, p) 2*p(3)*(x-p(1))/(p(2)^2)./( 1+((x-p(1))/p(2)).^2 ).^2;
        dw = @(x, p) 2*p(3)*((x-p(1)).^2)/(p(2)^3)./( 1+((x-p(1))/p(2)).^2 ).^2;

        J = zeros(numel(x), num_params);
        for i = 1:num_curves
            pos = params(i);
            w = params(i+num_curves);
            a = params(i+2*num_curves);

            p = [pos w a];
            J(:,i) = dpos(x, p);
            J(:,i+num_curves) = dw(x,p);
            J(:,i+2*num_curves) = lorentz_fun(p, x)/a;
            if addBackground
                J(:,end) = 1;
            end
        end
    end
end