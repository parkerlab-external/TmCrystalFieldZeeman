%% Setup
%Constants
tol = 1e-12;
tol_dig = 10;
L = 3;
S = 1/2;
I = 1/2;


%SpOrbC = -8771.45/(7/2);
load(fullfile("Saves", "coarse_levels.mat"));
SpOrbC = -mean(eE)/(7/2);
clear eE gE
c = 299792458*1e2; %speed of light in cm/s

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
Lx = (Lplus + Lminus)/2;
Ly = (Lplus - Lminus)/(2i);

Sz = kron(AI(L),Az(S));
S2 = kron(AI(L),A2(S));
Splus = kron(AI(L), Aplus(S));
Sminus = Splus';

LdotS = Lz*Sz + (Lplus*Sminus + Lminus*Splus)/2;
J2 = L2 + 2*LdotS + S2;
Jz = Lz + Sz;
Jplus = Lplus + Splus;
Jminus = Lminus + Sminus;

%Hyperfine (to be used later)
Iplus = Aplus(I);
Iminus = Aminus(I);


%% Crystal Field
load(fullfile("Saves","CFparameters.mat"));
CF = 0;
for i = 1:numel(CFparams)
    CF = CF + CFparams(i)*Stevens(CFlabels(i,:),L2,Lz,Lplus);
end

%% Full Basis Setup
fullBasis = [];
HFL = 0; %Hyperfine interaction reconstructed in L
CFL = 0; %Crystal Field reconstructed in L
IZL = 0;

%% Spin Orbit
SO = SpOrbC*LdotS;

%Diagonalize the spin-orbit:
[U, SOE] = eig(SO);
[SOE, sortind] = sort(diag(round(SOE,tol_dig)));
U = U(:,sortind);

%Look for degenerate s-o levels:
for subSOE = unique(SOE)'
    SOlevel_inds = find(subSOE == SOE);
    USO = U(:,SOlevel_inds);
    
    %Degeneracies should correspond to a single J2
    subJ2 = unique(round(diag(USO'*J2*USO), tol_dig));
    if numel(subJ2) ~= 1
        error("Spin Orbit is not degenerate in J2, problem")
    end

    %Constants for J:
    if subJ2 == 7/2*9/2
        AHF = -(1e6)/c*(1497/4);
    elseif subJ2 == 5/2*7/2
        AHF = -(1e6)/c*(2115/3);
    end

    %Back out the crystal field:
    CFL = CFL + USO*(USO'*CF*USO)*USO';

    %Diagonalize the crystal field
    [V, CFE] = eig(USO'*CF*USO);
    [CFE, sortind] = sort(diag(round(CFE,tol_dig)));
    V = V(:,sortind);
    for subCFE = unique(CFE)'
        CFlevel_inds = find(CFE == subCFE);
        VCF = V(:,CFlevel_inds);

        % IdotJ = Iz*Jz + (Iplus*Jminus + Jminus*Iplus)/2
        CFI = eye(numel(CFlevel_inds));
        HFbasis = @(A) VCF'*USO'*A*USO*VCF;
        IdotJ = kron(CFI, Az(I))*kron(HFbasis(Jz), AI(I)) + (kron(CFI,Iplus)*kron(HFbasis(Jminus),AI(I)) + kron(HFbasis(Jplus),AI(I))*kron(CFI, Iminus))/2;
        HF = AHF*IdotJ;

        %Hyperfine Energies
        [T, HFE] = eig(HF);
        [~, sortind] = sort(diag(round(HFE,tol_dig)));
        T = T(:,sortind);

        %Back out the basis and hyperfine operator
        Tout =  kron(USO,AI(I))*kron(VCF,AI(I))*T;
        HFL = HFL + Tout*HFE*Tout';
        fullBasis = [fullBasis Tout];
    end
end

%% Full Energies
%We have two operators which differ slightly: One says that Crystal field
%is a perturbation on the spin-orbit, the other says that it is not
SO = kron(SO, AI(I));
CFL = kron(CFL, AI(I));
CF = kron(CF, AI(I));

%CF as perturbation
fullHL = SO + CFL + HFL;

%CF from L purely
fullH = SO + CF + HFL;

%Make Hermitian
fullH = (fullH + fullH')/2;
fullHL = (fullHL + fullHL')/2;

%% Perturbative Crystal Field Output
fullE = sort(eig(fullHL));
gE = fullE(1:16);
eE = fullE(17:end);
fullLines = eE - gE';

line_positions = sort(fullLines(:));
% save("Lines\hf_lines_CFJ.mat", "line_positions");

%Get fit levels
[U, E] = eig(fullHL);
[~,sortind] = sort(diag(E));
U = U(:,sortind);
CFLevels = round(diag(U'*(SO + CFL)*U), tol_dig);
HFShifts = round(diag(U'*HFL*U), tol_dig);
mF = kron(Jz, eye(2)) + kron(eye(size(Jz)), I*[1 0; 0 -1]);
mF = diag(U'*mF*U);

save(fullfile("Saves","CFHFLevels"), "CFLevels", "HFShifts", "mF");

eE = [];
gE = [];
eHFShifts = {};
gHFShifts = {};
for CFLevel = unique(CFLevels)'
    level_inds = CFLevel == CFLevels;

    if CFLevel > 0
        eE = [eE CFLevel];
        eHFShifts = [eHFShifts {HFShifts(level_inds)}];
    else
        gE = [gE CFLevel];
        gHFShifts = [gHFShifts {HFShifts(level_inds)}];
    end
end
shift = mean(gE);
gE = gE - shift;
eE = eE - shift;


% save("Levels\fit_levels.mat", "eE", "gE", "eHFShifts", "gHFShifts");

%We want to save the CF positions and HF shifts in a 12x16 matrix
%Just going to hard code the numbers because I'm tired
CF_positions = CFLevels(17:end) - CFLevels(1:16)';
HF_shifts = HFShifts(17:end) - HFShifts(1:16)';

% save("Lines\ordered_lines.mat", "CF_positions", "HF_shifts");

%Save Operators too
J = [];
J.z = Jz;
J.sq = J2;
L = [];
L.z = Lz;
L.x = Lx;
L.y = Ly;

S = [];
S.z = Sz;
S.x = (Splus + Sminus)/2;
S.y = (Splus - Sminus)/(2i);

save(fullfile("Saves","full_hamiltonian.mat"), "SO", "CFL", "HFL", "L", "J", "S", "I", "fullBasis");


%% Non-perturbative Crystal Field Output
%Warning: This is not the right way to do it, deleting this part by
%commenting

% fullE = sort(real(eig(fullH)));
% gE = fullE(1:16);
% eE = fullE(17:end);
% fullLines = eE - gE';
% 
% 
% line_positions = sort(fullLines(:));
% % save("Lines\hf_lines_CFL.mat", "line_positions");

%% Plotting
figure("Name", "Crystal Field Fitting")
load(fullfile("Saves", "rough_position_scan.mat"));
plot(scan_wvs, scan);
hold on
%Crystal Field Fit positions:
xline(1e7./uniquetol(line_positions, tol));

%Level Fit lines with hyperfine:
load(fullfile("Saves", "coarse_levels.mat"));
coarse_positions = repelem(sort(eE),4)' - repelem(sort(gE),4);
% coarse_positions = repelem(line_positions, 4, 4);
line_positions = coarse_positions + HF_shifts;
xline(1e7./uniquetol(line_positions(:), tol), "blue")





