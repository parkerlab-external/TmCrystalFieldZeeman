function cg = clebsch_gordan(j1, m1,j2, m2, j, m)
    if j1 < 0 || j2 < 0 || j < 0
        cg = 0;
        return;
    end
    if m1+m2 ~= m
        cg = 0;
        return;
    end
    if (j > j1 + j2) || (j < abs(j1 - j2))
        cg = 0;
        return;
    end
    
    %Rename factorial
    f = @(x) factorial(x);

    %Explicit expression for Clebsch-Gordan Coeffs used:
    first_term = sqrt((2*j+1)*f(j+j1-j2)*f(j-j1+j2)*f(j1+j2-j)/f(j1+j2+j+1));
    second_term = sqrt(f(j+m)*f(j-m)*f(j1-m1)*f(j1+m1)*f(j2-m2)*f(j2+m2));

    t = [j1+j2-j, j1-m1, j2+m2, j-j2+m1, j-j1-m2]; %These will always be positive or 0
    max_k = min(t(1:3));
    min_k = max([-t(4:5) 0]);

    third_term = 0;
    for k = min_k:max_k
        third_term = third_term + ((-1)^k)/(f(k)*f(t(1)-k)*f(t(2)-k)*f(t(3)-k)*f(t(4)+k)*f(t(5)+k));
    end

    cg = first_term*second_term*third_term;
end
