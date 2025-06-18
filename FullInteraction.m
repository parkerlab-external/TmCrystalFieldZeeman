function [xaxis, predicted] = FullInteraction(Bfieldmag, linewidth, type, sourcefile, Mweight, Eweight, options)
arguments
    Bfieldmag = 0; %tesla
    linewidth = 0;
    type = "unpolarized1";
    sourcefile ="fittedPQ.mat";
    Mweight = 1;
    Eweight = 1;
    options.xaxis;
    options.N = 500;
end

func_dir = fileparts(mfilename('fullpath'));
load(fullfile(func_dir,"Saves","full_hamiltonian.mat"));
%Massaging
Iz = kron(eye(size(L.x)), [1 0; 0 -1]/2);
Ix = kron(eye(size(L.x)), [0 1; 1 0]/2);
Iy = (Iz*Ix - Ix*Iz)/(1i);
Lx = kron(L.x,eye(2));
Lz = kron(L.z,eye(2));
Sz = kron(S.z,eye(2));
Sx = kron(S.x,eye(2));
Ly = (Lz*Lx - Lx*Lz)/(1i);
Sy = (Sz*Sx - Sx*Sz)/(1i);


load(fullfile(func_dir, "Saves", sourcefile));
temp = reshape(temp, [12, 16]);
load(fullfile(func_dir, "Saves", "rough_position_scan"));
load(fullfile(func_dir, "Saves", "HFampfitting"));

if ~isfield(options,"xaxis")
    %Set up the plotting axis
    scan_wns = 1e7./scan_wvs;
    wn_range = max(scan_wns) - min(scan_wns);
    xaxis = linspace(min(scan_wns) - 0.05*wn_range , max(scan_wns) + 0.05*wn_range, 100000);
    xaxis_coarse = linspace(min(scan_wns) - 0.05*wn_range , max(scan_wns) + 0.05*wn_range, 5000);
    xaxis(min(abs(xaxis - amp_fitted_positions(:))) > 0.1) = [];
    xaxis = sort([xaxis xaxis_coarse]);
else
    xaxis = options.xaxis;
end

%We need the correction from crystal field levels to the amplitude fitted
%levels
load(fullfile(func_dir, "Saves", "CFHFLevels.mat"));
CFLevels = CFLevels + HFShifts;
CFLevels = CFLevels(17:28) - CFLevels(1:16)';
CF_correction = amp_fitted_positions - CFLevels;

amp_fitted_widths = amp_fitted_widths';


predicted = zeros(size(xaxis));
predicted_mag = predicted;
predicted_elec = predicted;

% Constants
N = options.N;
Bmag = Bfieldmag;
Btheta = acos(2*rand(1,N)-1);
Bphi = 2*pi*rand(1,N);
Kphi = 2*pi*rand(1,N);
stdwidth = linewidth;


tol = 1e-11;
tol_dig = 10;
muB = 9.27e-24/(6.36e-34)/(3e10); %cm-1 per T

%turn off E1 or M1
elec_weight = Eweight*elec_weight;
mag_weight = Mweight*mag_weight;


E = cell(1,3);
E(:) = {zeros(size(Lz))};
M = cell(1,3);
M(:) = {zeros(size(Lz))};
for i = 1:size(full_elec_transitions,1)
    for x = 1:size(full_elec_transitions, 2)
        E{x} = E{x} + elec_weight(i)*(full_elec_transitions{i,x});
    end
end
for i = 1:size(full_mag_transitions,1)
    for x = 1:size(full_mag_transitions,2)
        M{x} = M{x} + mag_weight(i)*(full_mag_transitions{i,x});
    end
end

%Selecting polarization for E and M
Ep = @(pvec) pvec(1)*E{1} + pvec(2)*E{2} + pvec(3)*E{3};
Mp = @(pvec) pvec(1)*M{1} + pvec(2)*M{2} + pvec(3)*M{3};

M1 = Mp([1i 0 1i]/sqrt(2));
E2 = Ep([0 0 0]);

%Polarization
% Bpol = [0;0;1];
% Bpol = [0;0;1];
% Epol = [1;0;0];
% Epol = [0;1;0];

switch type
    case "unpolarized1"
        Ktheta = acos(2*rand(1,N)-1);
        Bpol = [0;0;1];
        Epol = [1;0;0];
    case "unpolarized2"
        Ktheta = acos(2*rand(1,N)-1);
        Bpol = [1;0;0];
        Epol = [0;0;1];
    case "parallel1" %Electric field along B field
        Ktheta = 0*acos(2*rand(1,N)-1);
        Bpol = [1;0;0];
        Epol = [0;0;1];
    case "transverse1" %Electric field perpendicular to B field
        Ktheta = 0*acos(2*rand(1,N)-1);
        Bpol = [0;0;1];
        Epol = [1;0;0];
    case "parallel2"
        Ktheta = pi/2 + 0*acos(2*rand(1,N)-1);
        Epol = [0;1;0];
        Bpol = [0;0;1];
    case "transverse2"
        Ktheta = pi/2 + 0*acos(2*rand(1,N)-1);
        Epol = [0;0;1];
        Bpol = [0;1;0];
    otherwise %Unpolarized
        Ktheta = acos(2*rand(1,N)-1);
        Bpol = [0;0;1];
        Epol = [1;0;0];
end

%Magnetic field and full hamiltonian
B = Bmag*muB*(Lz + 2*Sz);
H = SO + CFL + HFL;


%Now we have all operators point "up". To proceed, we can keep the
%crystal field pointed up and rotate the magnetic field and excitation

%Rotating the magnetic field and excitation together
Rn = @(theta, n) expm(-1i*theta*(n(1)*(Lx+Sx+Ix) + n(2)*(Ly+Sy+Iy) + n(3)*(Lz+Sz+Iz)));
Rx = @(theta) Rn(theta, [1 0 0]);
Ry = @(theta) Rn(theta, [0 1 0]);
Rz = @(theta) Rn(theta, [0 0 1]);

%Basis transform from vectors to operators
U = [-1 1i 0;
    0 0 sqrt(2);
    1 1i 0]/sqrt(2);

%Rotations of vectors, I want axis angle repn:
genx = [0 0 0;
    0 0 1;
    0 -1 0];
geny = [0 0 -1;
    0 0 0;
    1 0 0];
genz = [0 1 0;
    -1 0 0;
    0 0 0];
Rvec = @(vec, theta) expm(theta*(vec(1)*genx + vec(2)*geny + vec(3)*genz));
Rvecx = @(theta) Rvec([1 0 0], theta);
Rvecy = @(theta) Rvec([0 1 0], theta);
Rvecz = @(theta) Rvec([0 0 1], theta);

%Simulation part
for i = 1:N
    %Rotate B
    Brotated = Rz(Bphi(i))'*Rx(Btheta(i))'*B*Rx(Btheta(i))*Rz(Bphi(i));

   
    [Trotated, rotatedE] = eig(H + Brotated);
    [rotatedE, sortinds] = sort(diag(round(rotatedE, tol_dig)));
    Trotated = Trotated(:,sortinds);
    full_positions = rotatedE' - rotatedE;
    line_positions = full_positions(1:16,17:28) + CF_correction';
    
    %Rotate the transitions - Rotate lagging behind B, then with B, then
    %around B. Then convert to eigenbasis to get amplitudes
    RX1 = Rx(Btheta(i) - Ktheta(i));
    Bdirection = [sin(Bphi(i))*sin(Btheta(i)) cos(Bphi(i))*sin(Btheta(i))  cos(Btheta(i))];
%     Bdirection2 = Rvecz(Bphi(i))*Rvecx(Btheta(i))*[0;0;1]; %Check
    RB = expm(-1i*Kphi(i)*(Bdirection(1)*(Lx+Sx+Ix) + Bdirection(2)*(Ly+Sy+Iy) + Bdirection(3)*(Lz+Sz+Iz)));
    
    %Rotate Magnetic Operators
%     Mrotated = RB'*Rz(Bphi(i))'*RX1'*M1*RX1*Rz(Bphi(i))*RB;
    Mrotated = Mp(U*Rvec(Bdirection, Kphi(i))*Rvecz(Bphi(i))*Rvecx(Btheta(i)-Ktheta(i))*Bpol);
        
    %Sum of squares
    Mamps = abs(Trotated'*Mrotated*Trotated).^2;
    Mamps = Mamps(1:16, 17:28);

    %Rotate Electric Operators
%     Erotated = RB'*Rz(Bphi(i))'*RX1'*E2*RX1*Rz(Bphi(i))*RB;
    Erotated = Ep(U*Rvec(Bdirection, Kphi(i))*Rvecz(Bphi(i))*Rvecx(Btheta(i)-Ktheta(i))*Epol);

    Eamps = abs(Trotated'*Erotated*Trotated).^2;
    Eamps = Eamps(1:16, 17:28);

    Eamps = Eamps.*temp';
    Mamps = Mamps.*temp';

    line_amplitudes = Mamps + Eamps;

    for j = 1:numel(line_amplitudes)
        if linewidth == 0
            stdwidth = amp_fitted_widths(j);
        end
        predicted_mag = predicted_mag + multi_lorentz_fun([line_positions(j) stdwidth Mamps(j)/stdwidth], xaxis, 1);
        predicted_elec = predicted_elec + multi_lorentz_fun([line_positions(j) stdwidth Eamps(j)/stdwidth], xaxis, 1);
    end
end
predicted = (predicted_mag + predicted_elec)/N;
end







