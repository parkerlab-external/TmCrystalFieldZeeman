load(fullfile("Saves","coarse_levels.mat"));
gE = repelem(gE, 2);
eE = repelem(eE - mean(eE), 2);
fit_to_operators = [
    2 0;
%     2 1;
%     2 2;
%     4 0;
%     4 3;
%     6 0;
%     6 6;
    ];
num_operators = size(fit_to_operators,1);

J2 = @(J) J*(J+1)*eye(2*J+1);
Jz = @(J) diag(J:-1:-J);
Jplus = @(J) diag( sqrt(J*(J+1) - ((J-1):-1:-J).*(J:-1:(-J+1))) , 1 );


J2g = 7/2*9/2*eye(8);
Jzg = diag(7/2:-1:-7/2);
Jplusg = diag(sqrt(7/2*9/2 - (5/2:-1:-7/2).*(7/2:-1:-5/2)),1);
testOp = Stevens([2 0], J2g, Jzg, Jplusg);

J2e = 7/2*5/2*eye(6);
Jze = diag(5/2:-1:-5/2);
Jpluse = diag(sqrt(7/2*5/2 - (3/2:-1:-5/2).*(5/2:-1:-3/2)),1);

gCFparts = zeros(num_operators, size(testOp,1), size(testOp,2));
eCFparts = zeros(num_operators,6,6);
for i = 1:num_operators
    gCFparts(i,:,:) = Stevens(fit_to_operators(i,:), J2g, Jzg, Jplusg);
    eCFparts(i,:,:) = Stevens(fit_to_operators(i,:), J2e, Jze, Jpluse);
end

%Fitting to both the excited and ground states
obj_fun = @(params) norm(sort(gE)' - sort(eig(squeeze(sum(params.*gCFparts, 1))))) + norm(sort(eE)' - sort((9/5)*eig(squeeze(sum(params.*eCFparts, 1)))));

%Fitting to just the ground states
% obj_fun = @(params) norm(sort(gE)' - sort(eig(squeeze(sum(params.*gCFparts, 1)))));

paramsGuess = ones(num_operators,1);
[params] = fminsearch(obj_fun,paramsGuess);

res = obj_fun(params);

gCF = 0;
eCF = 0;
for i = 1:num_operators
    gCF = gCF + params(i)*Stevens(fit_to_operators(i,:),J2g, Jzg, Jplusg);
    eCF = eCF + (9/5)*params(i)*Stevens(fit_to_operators(i,:),J2e, Jze, Jpluse);
end



figure("Name","Crystal Field Fitting")
subplot(2,2,1)
xline(eE)
xline(eig(eCF),"blue")
title("Excited From J")
subplot(2,2,2)
xline(gE)
xline(eig(gCF),"blue")
title("Ground from J")


LCF = 0;
for i = 1:num_operators
    LCF = LCF + params(i)*Stevens(fit_to_operators(i,:), J2(3), Jz(3), Jplus(3));
end

LCF = kron(LCF, eye(2));
Sz = kron(eye(7), Jz(1/2));
Lz = kron(Jz(3), eye(2));
Lplus = kron(Jplus(3), eye(2));
Splus = kron(eye(7), Jplus(1/2));
Lminus = Lplus';
Sminus = Splus';
LdotS = Lz*Sz + (Lplus*Sminus + Lminus*Splus)/2;
J2 = kron(J2(3),eye(2)) + 2*LdotS + kron(eye(7), J2(1/2));

[U, ~] = eig(J2);
J2 = U'*J2*U;
LCF = 1.4*U'*LCF*U;

tol = 1e-12;
subplot(2,2,3);
xline(eE)
m = round(diag(J2),12) == (5/2)*(7/2);
xline(eig(LCF(m,m)),"blue")
title("Excited From L")
subplot(2,2,4);
xline(gE)
m = round(diag(J2),12) == (7/2)*(9/2);
xline(eig(LCF(m,m)),"blue")
title("Ground From L")

CFparams = 1.4*params;
CFlabels = fit_to_operators;
save(fullfile("Saves", "CFJparameters.mat"), "CFparams", "CFlabels")

