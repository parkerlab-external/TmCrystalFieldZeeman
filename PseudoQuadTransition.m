tol_dig = 10;

%Steal The operators from HFCoupling
load(fullfile("Saves","full_hamiltonian.mat")); %Includes SO CFL and HFL, and Jz

%We want everything in the J Basis (NOT necessarily eignebasis for crystal field):
[V, ~] = eig(J.sq);

%Manually create the basis?
L=repelem(3:-1:-3, 2);
S=repmat(1/2:-1:-1/2, 1,7);
Jtemp = [repelem(7/2, 8) repelem(5/2, 6)];
mJ = [7/2:-1:-7/2 5/2:-1:-5/2];
V = zeros(14);
for j = 1:14
    for i = 1:14
        V(i,j) = clebsch_gordan(3,L(i),1/2,S(i),Jtemp(j),mJ(j));
    end
end



Jz = round(diag(V'*J.z*V), tol_dig);
J2 = round(diag(V'*J.sq*V), tol_dig);
[~, sortind] = sortrows([J2  abs(Jz) Jz], [1 2 3], 'descend');

V = V(:,sortind);
Jz = Jz(sortind);
J2 = J2(sortind);
J = J2;
J(J == 7/2*9/2) = 7/2;
J(J == 5/2*7/2) = 5/2;

%Coupling operators
%PQuad
op1 = 2;
dir1 = [-1 2 1];

%M1
op2 = 1;
dir2 = [1 0 -1];

%% Prepare transition matrices
transition1 = cell(1,numel(dir1));
transition1(:) = {zeros(size(V))};

transition2 = cell(1,numel(dir2));
transition2(:) = {zeros(size(V))};

for g = find(J == 7/2)'
    for e = find(J == 5/2)'
        %PseudoQuadrupole as transition1
        for x = 1:numel(dir1)
            if dir1(x) == 2
                transition1{x}(g, e) = transition1{x}(g,e)...
                                + 1/sqrt(2)*clebsch_gordan(J(e), Jz(e), op1, dir1(x), J(g), Jz(g))...
                                + 1/sqrt(2)*clebsch_gordan(J(e), Jz(e), op1, -dir1(x), J(g), Jz(g));
            else
                transition1{x}(g, e) = transition1{x}(g,e)...
                                + clebsch_gordan(J(e), Jz(e), op1, dir1(x), J(g), Jz(g));
            end
        end

        %M1 as transition2
        for x = 1:numel(dir2)
            transition2{x}(g, e) = transition2{x}(g,e)...
                                + clebsch_gordan(J(e), Jz(e), op2, dir2(x), J(g), Jz(g));
        end
    end
end


%% Convert to Eigenbasis
%Convert transitions from current basis (J), back into LSI basis through V,
%then into eigenbasis
fullH = SO + CFL + HFL;
[T, fullE] = eig(round(fullH,tol_dig));
[fullE, sortind] = sort(diag(fullE));
T = T(:,sortind);

for i = 1:numel(transition1)
    transition1{i} = T'*kron(V*transition1{i}*V', eye(2))*T;
end
for i = 1:numel(transition2)
    transition2{i} = T'*kron(V*transition2{i}*V', eye(2))*T;
end

%% Extract Magnetic and Electric Transitions
%transition2 is the magnetic transition
mag_transitions = cell(size(transition2));
full_mag_transitions = mag_transitions;
keep_inds = [];
for i = 1:size(transition2,1)
    temp_cell = cell(1, size(transition2, 2));
    keep = false;
    for x = 1:size(transition2,2)
        cross_terms = transition2{i,x}(1:16, 17:28)';
        mag_transitions{i,x} = cross_terms;

        %Remove non-transition terms
        transition2{i,x}(:, 1:16) = 0;
        transition2{i,x}(17:28, 17:28) = 0;

        full_mag_transitions{i,x} = T*transition2{i,x}*T';
        if any(cross_terms, "all")
            keep = true;
        end
    end
    if keep
        keep_inds = [keep_inds i];
    end
end
mag_transitions = mag_transitions(keep_inds,:);
full_mag_transitions = full_mag_transitions(keep_inds, :);


%transition1 is the electric transition
elec_transitions = cell(size(transition1));
full_elec_transitions = elec_transitions;
keep_inds = [];
for i = 1:size(transition1, 1)
    keep = false;

    for x = 1:size(transition1, 2)
        cross_terms = transition1{i,x}(1:16, 17:28)';
        elec_transitions{i,x} = cross_terms;

        %Remove non-transition terms
        transition1{i,x}(:, 1:16) = 0;
        transition1{i,x}(17:28, 17:28) = 0;

        full_elec_transitions{i,x} = T*transition1{i,x}*T';
        if any(cross_terms, "all")
            keep = true;
        end
    end
    if keep
        keep_inds = [keep_inds i];
    end
end
elec_transitions = elec_transitions(keep_inds, :);
full_elec_transitions = full_elec_transitions(keep_inds,:);



%% Fitting to known levels
load(fullfile("Saves", "HFampfitting.mat")); %Contains both the positions and amplitudes
line_positions = amp_fitted_positions;
load(fullfile("Saves", "coarse_levels.mat"));


pretempscale = -repmat(repelem(gE,4),[12 1]);
% kb = 0.69503;
% tempscaling = line_positions - mean(line_positions,"all");
% tempscaling = tempscaling*0+1; %temporary

[~, uinds, ind_mask] = uniquetol(line_positions, 1e-10);
ind_mask = reshape(ind_mask, size(line_positions));

%Customize some indices which should be averaged
%Close levels at 8765.97
ind_mask(7:8,15:16) = ind_mask(5,13);
%Entire line at 8778
ind_mask(9:12,1:4) = ind_mask(9,1);
%Entire line at 8774.28
ind_mask(5:8,1:4) = ind_mask(5,1);
%Close lines at 8769.9
%These lines seem huge, try removing them from the fit completely
ind_mask(7:8,11:12) = ind_mask(5,9);
ind_mask(5:8,9:12) = -1; %-1 removes these lines from the fitting

%% Save Transitions
save(fullfile("Saves","PQTransition.mat"), "mag_transitions", "full_mag_transitions", "elec_transitions", "full_elec_transitions", "pretempscale","ind_mask");