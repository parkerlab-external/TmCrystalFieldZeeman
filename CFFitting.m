%Load the level fit levels
load(fullfile("Saves","coarse_levels.mat"));
gE = repelem(gE,2);
eE = repelem(eE - mean(eE),2);



%Stevens Operators
fitting_stevens = [
    2 0;
%     2 2;
%     2 1;
%     2 -2;
%     2 -1;
%     4 0;
%     4 4;
%     4 -4;
%     6 6;
    ];
num_fitting = size(fitting_stevens, 1);
params = zeros(1,num_fitting);


%Constants
tol = 1e-12;
tol_dig = 12;
L = 3;
S = 1/2;
I = 1/2;

%Angular Operators
AI = @(J) eye(2*J+1);
Az = @(J) diag(J:-1:-J);
A2 = @(J) J*(J+1)*eye(2*J+1);
Aplus = @(J) diag( sqrt(J*(J+1) - ((J-1):-1:-J).*(J:-1:(-J+1))) , 1 );
Aminus = @(J) Aplus(J)';

%LS Coupling
Lz = kron(Az(L), AI(S));
L2 = kron(A2(L), AI(S));
Lplus = kron(Aplus(L), AI(S));
Lminus = Lplus';

Sz = kron(AI(L),Az(S));
S2 = kron(AI(L),A2(S));
Splus = kron(AI(L), Aplus(S));
Sminus = Splus';

LdotS = Lz*Sz + (Lplus*Sminus + Lminus*Splus)/2;
J2 = L2 + 2*LdotS + S2;



%% First Pass
%Make an ndgrid for the parameter space (brute force)
if num_fitting == 1
    param_steps = 100;
elseif num_fitting == 2
    param_steps = 40;
elseif num_fitting == 3
    param_steps = 20;
else
    param_steps = 10;
end
param_ss = 2;
percent = 1/param_ss^num_fitting; %Look for bottom precentage of residuals

param_max = zeros(1,num_fitting);
paramRanges = cell(1,num_fitting);
paramGrid = cell(1, num_fitting);
max_diff = 1.01*max([gE eE] - [gE eE]',[],"all");
for i = 1:num_fitting
    param_max(i) = max_diff/norm(Stevens(fitting_stevens(i,:), L2, Lz, Lplus));
    paramRanges{i} = linspace(-param_max(i),param_max(i),param_steps);
end
[paramGrid{:}] = ndgrid(paramRanges{:});

%Diagonalize LdotS
[U, ~] = eig(LdotS);
J2 = round(diag(U'*J2*U),tol_dig);


num_passes = 6;
figure("Name", "Brute Force CF Fitting");
tiledlayout('flow')
for passes = 1:num_passes

residuals = zeros(size(paramGrid{1}));

for fitnum = 1:numel(residuals)
    CF = 0;
    for i = 1:num_fitting
        CF = CF + paramGrid{i}(fitnum)*Stevens(fitting_stevens(i,:),L2,Lz,Lplus);
    end
    CF = U'*CF*U;

    % Isolate Each J manifold
    for subJ2 = unique(J2)'
        level_inds = J2 == subJ2;
        CFlevels = sort(round(eig(CF(level_inds, level_inds)),tol_dig))';
        if ~isreal(CFlevels)
            error("Eigenvalues of Crystal Field are not Real")
        end

        if subJ2 == 5/2*7/2
            residuals(fitnum) = residuals(fitnum) + norm(sort(eE) - CFlevels);
        elseif subJ2 == 7/2*9/2
            residuals(fitnum) = residuals(fitnum) + norm(sort(gE) - CFlevels);
        end
    end
end

%Plotting
if ismember(num_fitting, [1 2 3])
    nexttile;
end
if num_fitting == 1
    [~,sortind] = sort(paramGrid{1});
    plot(sort(paramGrid{1}), residuals(sortind));
end
if num_fitting == 2
    scatter3(paramGrid{1}(:), paramGrid{2}(:), residuals(:), 'filled');
end
if num_fitting == 3
    minmask = residuals < prctile(residuals(:), percent*100);
    plotresiduals = (residuals(minmask) - min(residuals(minmask)));
    plotresiduals = plotresiduals./max(plotresiduals);
    scatter3(paramGrid{1}(minmask), paramGrid{2}(minmask), paramGrid{3}(minmask), plotresiduals*(-100) + 120, 'filled')
    hold on
    minmask = residuals == min(residuals, [], "all");
    scatter3(paramGrid{1}(minmask), paramGrid{2}(minmask), paramGrid{3}(minmask), 125, 'filled', 'red')
    hold off
    xlabel(sprintf("O%d,%d", fitting_stevens(1,1), fitting_stevens(1,2)));
    ylabel(sprintf("O%d,%d", fitting_stevens(2,1), fitting_stevens(2,2)));
    zlabel(sprintf("O%d,%d", fitting_stevens(3,1), fitting_stevens(3,2)));

    xlim([min(paramGrid{1},[],"all") max(paramGrid{1},[],"all")])
    ylim([min(paramGrid{2},[],"all") max(paramGrid{2},[],"all")])
    zlim([min(paramGrid{3},[],"all") max(paramGrid{3},[],"all")])
end

if ismember(num_fitting, [1 2 3])
    title(sprintf("Step %d", passes));
end
%End plotting

%Find minimum
[~,minmask] = min(residuals, [], "all");
for i = 1:num_fitting
    params(i) = paramGrid{i}(minmask);
end



%Alter the search ranges for the second pass
%Create an ND grid around each point
% minmask = residuals < (1-percent)*min(residuals, [],"all") + percent*max(residuals, [], 'all');
minmask = residuals < prctile(residuals(:), percent*100);
fprintf("Choosing %d out of %d points\n", sum(minmask,"all"), numel(residuals));

for i = 1:num_fitting
    if passes == 1
        param_max(i) = param_max(i)/(param_steps-1)*(1 - 1/param_ss);
    else
        param_max(i) = param_max(i)/(param_ss-1)*(1 - 1/param_ss);
    end
    paramRanges{i} = linspace(-param_max(i), param_max(i), param_ss);
end
[paramRanges{:}] = ndgrid(paramRanges{:});
%Add the small grid to the previous grid centers
for i = 1:num_fitting
    paramGrid{i} = paramGrid{i}(minmask) + paramRanges{i}(:)';
    paramGrid{i} = paramGrid{i}(:);
end
end



CFparams = params;
CFlabels = fitting_stevens;
save(fullfile("Saves", "CFparameters.mat"), "CFparams", "CFlabels")