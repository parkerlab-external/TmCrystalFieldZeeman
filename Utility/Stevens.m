function operator = Stevens(ang_mom, J2,Jz,Jplus)
    %function for transforming angular indices (j,m) into array indices
    %(j-m (m>0), j+m (m<0))
    ang_ind = @(ang_mom) sub2ind([8 8], ang_mom(1)+1-(ang_mom(2) > 0)*ang_mom(2), ang_mom(1)+1+(ang_mom(2) < 0)*ang_mom(2));
    
    Jminus = Jplus';
    Id = eye(size(J2));
    
    %Set of operator functions
    persistent O
    if isempty(O)
        O{ang_ind([1 1])} = @(Id, J2,Jz,Jplus,Jminus) (Jplus + Jminus)*0.5;
        O{ang_ind([1 0])} = @(Id, J2,Jz,Jplus,Jminus) Jz;
        O{ang_ind([1 -1])} = @(Id, J2,Jz,Jplus,Jminus) -1i*(Jplus - Jminus)*0.5;
        O{ang_ind([2 0])} = @(Id, J2,Jz,Jplus,Jminus)  3*Jz^2 - J2;
        O{ang_ind([2 1])} = @(Id, J2,Jz,Jplus,Jminus)  0.25*(Jz*(Jplus + Jminus) + (Jplus + Jminus)*Jz);
        O{ang_ind([2 -1])} = @(Id, J2,Jz,Jplus,Jminus)  -1i*0.25*(Jz*(Jplus - Jminus) + (Jplus - Jminus)*Jz);
        O{ang_ind([2 -2])} = @(Id, J2,Jz,Jplus,Jminus)  -1i*(Jplus^2 - Jminus^2)*0.5;
        O{ang_ind([2 2])} = @(Id, J2,Jz,Jplus,Jminus)  (Jplus^2 + Jminus^2)*0.5;
        O{ang_ind([4 0])} = @(Id, J2,Jz,Jplus,Jminus)  35*Jz^4 - (30*J2 - 25*Id)*Jz^2 + 3*J2^2 - 6*J2;
        O{ang_ind([4 4])} = @(Id, J2,Jz,Jplus,Jminus)  (Jplus^4 + Jminus^4)*0.5;
        O{ang_ind([4 -4])} = @(Id, J2,Jz,Jplus,Jminus)  -1i*(Jplus^4 - Jminus^4)*0.5;
        O{ang_ind([4 3])} = @(Id, J2,Jz,Jplus,Jminus)  0.25*((Jplus^3 + Jminus^3)*Jz + Jz*(Jplus^3 + Jminus^3));
        O{ang_ind([4 -3])} = @(Id, J2,Jz,Jplus,Jminus)  -1i*0.25*((Jplus^3 - Jminus^3)*Jz + Jz*(Jplus^3 - Jminus^3));
        O{ang_ind([6 0])} = @(Id, J2,Jz,Jplus,Jminus)  231*Jz^6 - (315*J2 - 735*Id)*Jz^4 + (105*J2^2 - 525*J2 + 294*Id)*Jz^2 - 5*J2^3 + 40*J2^2 - 60*J2;
        O{ang_ind([6 6])} = @(Id, J2,Jz,Jplus,Jminus)  (Jplus^6 + Jminus^6)*0.5;
        O{ang_ind([6 -6])} = @(Id, J2,Jz,Jplus,Jminus)  -1i*(Jplus^6 - Jminus^6)*0.5;
    end


    operatorFunc = O{ang_ind(ang_mom)};
    operator = operatorFunc(Id, J2, Jz, Jplus, Jminus);

    
end

