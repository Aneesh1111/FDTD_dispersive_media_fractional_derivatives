close all; clear variables; clc;

%% constants
eps0  = 8.8541878128e-12;
mu0 = 1.256637062e-6;
c0 = 1/sqrt(eps0*mu0);
imp0 = sqrt(mu0/eps0);

tic 

%% FDTD parameters
dt = 1e-12;%1.768e-12;
dx = 1e-3;%1.1e-3 ;%/ 1.2;
dy = dx;
imax = 500;  % max discretization in space
nmax = 2001; % max discretization in time
C = c0*dt*sqrt(2)/dx;  % Courant factor must be less than 1

%% dielectric permittivity parameters

tau = 318e-12;
eps_inf = 2;

alpha = 0.2;
beta = [0.3, 0.5, 0.9];
a = tau.^alpha;
b = tau.^beta .* [9, 2, 10];
eps_rs = 60;
% conductivity = 0.035;

% alpha = 0;
% beta = [0.9];
% a = 0;%tau.^alpha;
% b = tau.^beta .* [1];
% eps_rs = 60;

conductivity = 0.0;
N = length(alpha);  % E
M = length(beta);   % P

%% GPOF method (calculates coefficients for FDTD scheme)
NN = 2000;
gpof_max = 8;
L = gpof_max;

wk = zeros(1,NN);
if beta(1) == 0
    wk(1) = 0;
else
    wk(1) = 1;
end
for k = 1:NN-1
    wk(k+1) = wk(k)*(1 - (1+beta(1))/k);
end
k = 0:1:NN-1;
[b_gpof_P1, a_gpof_P1] = gpof(wk,k,gpof_max,L+1);


if length(beta) >= 2
    wk = zeros(1,NN);
    if beta(2) == 0
        wk(1) = 0;
    else
        wk(1) = 1;
    end
    for k = 1:NN-1
        wk(k+1) = wk(k)*(1 - (1+beta(2))/k);
    end
    k = 0:1:NN-1;
    [b_gpof_P2, a_gpof_P2] = gpof(wk,k,gpof_max,L+1);
end

if length(beta) >= 3
    wk = zeros(1,NN);
    if beta(3) == 0
        wk(1) = 0;
    else
        wk(1) = 1;
    end
    for k = 1:NN-1
        wk(k+1) = wk(k)*(1 - (1+beta(3))/k);
    end
    k = 0:1:NN-1;
    [b_gpof_P3, a_gpof_P3] = gpof(wk,k,gpof_max,L+1);
end

if length(beta) >= 4
    wk = zeros(1,NN);
    if beta(4) == 0
        wk(1) = 0;
    else
        wk(1) = 1;
    end
    for k = 1:NN-1
        wk(k+1) = wk(k)*(1 - (1+beta(4))/k);
    end
    k = 0:1:NN-1;
    [b_gpof_P4, a_gpof_P4] = gpof(wk,k,gpof_max,L+1);
end

if length(beta) >= 5
    wk = zeros(1,NN);
    if beta(5) == 0
        wk(1) = 0;
    else
        wk(1) = 1;
    end
    for k = 1:NN-1
        wk(k+1) = wk(k)*(1 - (1+beta(5))/k);
    end
    k = 0:1:NN-1;
    [b_gpof_P5, a_gpof_P5] = gpof(wk,k,gpof_max,L+1);
end

if length(beta) == 5
    a_gpof_P = [a_gpof_P1, a_gpof_P2, a_gpof_P3, a_gpof_P4, a_gpof_P5];
    b_gpof_P = [b_gpof_P1, b_gpof_P2, b_gpof_P3, b_gpof_P4, b_gpof_P5];
end

if length(beta) == 3
    a_gpof_P = [a_gpof_P1, a_gpof_P2, a_gpof_P3];
    b_gpof_P = [b_gpof_P1, b_gpof_P2, b_gpof_P3];
end

if length(beta) == 2
    a_gpof_P = [a_gpof_P1, a_gpof_P2];
    b_gpof_P = [b_gpof_P1, b_gpof_P2];
end

if length(beta) == 1
    a_gpof_P = a_gpof_P1;
    b_gpof_P = b_gpof_P1;
end

wk = zeros(1,NN);
if alpha(1) == 0
    wk(1) = 0;
else
    wk(1) = 1;
end
for k = 1:NN-1
    wk(k+1) = wk(k)*(1 - (1+alpha)/k);
end
k = 0:1:NN-1;
[b_gpof_E, a_gpof_E] = gpof(wk,k,gpof_max,L+1);

%% plotting constants
x = 1:imax;
y = 1:imax;

%% update parametes
Ex = zeros(imax, imax);
Ey = zeros(imax, imax);
Px = zeros(imax, imax);
Py = zeros(imax, imax);
Px_prev = zeros(imax, imax);
Py_prev = zeros(imax, imax);
Ex_prev = zeros(imax, imax);
Ey_prev = zeros(imax, imax);

Hz = zeros(imax, imax);
Hzx = zeros(imax, imax);
Hzy = zeros(imax, imax);

vj_Ex = zeros(imax, imax, gpof_max*N);
vj_Ey = zeros(imax, imax, gpof_max*N);
vj_Px = zeros(imax, imax, gpof_max*M);
vj_Py = zeros(imax, imax, gpof_max*M);

vj_Ex_prev = zeros(imax, imax, gpof_max*N);
vj_Ey_prev = zeros(imax, imax, gpof_max*N);
vj_Px_prev = zeros(imax, imax, gpof_max*M);
vj_Py_prev = zeros(imax, imax, gpof_max*M);

vj_Ex_sum = zeros(imax,imax,N);
vj_Ey_sum = zeros(imax,imax,N);
vj_Px_sum = zeros(imax,imax,M);
vj_Py_sum = zeros(imax,imax,M);

%% PML
PML_thickness = 50;  % PML thickness
sigma_x = zeros(imax, imax);
sigma_y = zeros(imax, imax);

m = 1;
sigma_x_max = 0.8*(m + 1)/(imp0*dx);
sigma_y_max = 0.8*(m + 1)/(imp0*dy);
for i = 1:imax
    if i < PML_thickness
        sigma_x(i,:) = ((PML_thickness-i)/PML_thickness)^3 * sigma_x_max;
    elseif i > imax-PML_thickness
        sigma_x(i,:) = ((i-(imax-PML_thickness))/PML_thickness)^3 * sigma_x_max;
    end
end
for i = 1:imax
    if i < PML_thickness
        sigma_y(:,i) = ((PML_thickness-i)/PML_thickness)^3 * sigma_y_max;
    elseif i > imax-PML_thickness
        sigma_y(:,i) = ((i-(imax-PML_thickness))/PML_thickness)^3 * sigma_y_max;
    end
end
sigma_x_star = sigma_x * imp0^2;
sigma_y_star = sigma_y * imp0^2;

factor = 1;
eps2 = eps0*factor;
mu2 = mu0*factor;

%% TFSF
tfsf_pos = 25;
pml_offset = 15;

tfsf_end = 75;

tfsf_left = PML_thickness + pml_offset;
tfsf_right = imax - PML_thickness - pml_offset;

tfsf_pos = 25;
tfsf_end = 75;
tfsf_left = 20;
tfsf_right = 80;

tfsf_pos = 100;
tfsf_end = 400;
tfsf_left = 100;
tfsf_right = 400;

tfsf_x = tfsf_right - tfsf_left;
tfsf_y = tfsf_end - tfsf_pos;

Hxz_tfsf_in_pml1 = 1:PML_thickness+pml_offset-1;
Hxz_tfsf_in_tfsf = PML_thickness+pml_offset:imax-1-PML_thickness-pml_offset;
Hxz_tfsf_in_pml2 = imax-PML_thickness-pml_offset:imax-1;

Ey_tfsf_in_pml1 = 2:PML_thickness+pml_offset;
Ey_tfsf_in_tfsf = PML_thickness+pml_offset+1:imax-PML_thickness-pml_offset;
Ey_tfsf_in_pml2 = imax-PML_thickness-pml_offset+1:imax;

%% circular dielectric
eps = eps0*ones(imax,imax);
simple_eps = ones(imax,imax);

dielectric = zeros(imax,imax);
radius = 30;  % in units of [dx]
center = imax/2;

for i = 1:imax
    for j = 1:imax
        if ((i-center)^2 + (j-center)^2) <= radius^2
            dielectric(j,i) = 1;
            simple_eps(j,i) = eps_inf;
        end
    end
end

x1 = center - radius;
x2 = center + radius;
y1 = center - radius;
y2 = center + radius;
media_start_x = x1;
media_end_x = x2;
media_start_y = y1;
media_end_y = y2;

% epsr = 2;
% eps(media_start_y:media_end_y,media_start_x:media_end_x) = epsr*eps(media_start_y:media_end_y,media_start_x:media_end_x);

%% PML constants
CaEx = (1 - sigma_y*dt/(2*eps2)) ./ (1 + sigma_y*dt/(2*eps2));
CbEx = (dt./(eps*dy)) ./ (1 + sigma_y*dt/(2*eps2));
CaEy = (1 - sigma_x*dt/(2*eps2)) ./ (1 + sigma_x*dt/(2*eps2));
CbEy = (dt./(eps*dx)) ./ (1 + sigma_x*dt/(2*eps2));
DaHzx = (1 - sigma_x_star*dt/(2*mu2)) ./ (1 + sigma_x_star*dt/(2*mu2));
DbHzx = (dt/(mu0*dx)) ./ (1 + sigma_x_star*dt/(2*mu2));
DaHzy = (1 - sigma_y_star*dt/(2*mu2)) ./ (1 + sigma_y_star*dt/(2*mu2));
DbHzy = (dt/(mu0*dy)) ./ (1 + sigma_y_star*dt/(2*mu2));

%% Electric field capture
Ex_trans = zeros(1, nmax);
Ex_pos2 = zeros(1, nmax);
Ex_1D = zeros(1, nmax);
H_tfsf = zeros(tfsf_y+1, tfsf_x+1, nmax);
Ex_tfsf = zeros(tfsf_y+1, tfsf_x+1, nmax);
Ey_tfsf = zeros(tfsf_y+1, tfsf_x+1, nmax);
Hz_ss = zeros(imax,imax);
Hz_ss_2pi = zeros(imax,imax);
Hz_ss_quarterpi = zeros(imax,imax);
Hz_ss_pi = zeros(imax,imax);
Hz_ss_pi_3quarter = zeros(imax,imax);

posx = 50;
posy = media_end_y;
posx2 = posx;
posy2 = media_start_y;

%% 1D source
source_pos = 5;

Ey_1D = zeros(1, imax+source_pos);
Hz_1D = zeros(1, imax+source_pos);
ty_mu = dt/(mu0*dy);

epsr_arr = eps0*ones(1, imax+source_pos);
% epsr_arr(media_start_y+source_pos:media_end_y+source_pos) = eps0;

mur_boundary_const_first = (c0*dt - dy)/(c0*dt + dy);

%% FDTD constants
C5 = conductivity/2 + eps_inf* eps0/dt;
C4 = conductivity/2 - eps_inf* eps0/dt;
C2 = sum(a./(dt.^alpha));
C1 = sum(b./(dt.^beta));
C3 = (eps0*eps_rs*(1 + C2)/C5);

%% GIF animation
numFrames = 100;   % Number of frames in the animation
frames(numFrames) = struct('cdata',[],'colormap',[]);
frame = 1;

fig = figure;
axis tight manual;
axis([0 imax -0.8, 0.8]);

%% FDTD scheme
for n = 2:nmax-1
    
    %% 1D source update equations
        
    % boudnary condition save parameters
    Ey_1D_prev_imax = Ey_1D(end);
    Ey_1D_prev_imax1 = Ey_1D(end-1);
    
    Hz_1D_prev_2 = Hz_1D(2);
    Hz_1D_prev_1 = Hz_1D(1);

    % update electric field
    Ey_1D(1:end-1) = Ey_1D(1:end-1) + dt./(dy*eps0) .* diff(Hz_1D);  % simple dielectric 
    % update electric boundaries
    Ey_1D(end) = Ey_1D_prev_imax1 + mur_boundary_const_first*(Ey_1D(end-1) - Ey_1D_prev_imax);
    % electric field source
    Ey_1D(source_pos) = Ey_1D(source_pos) + source_function(n, dt);
    
    % magnetic field update
    Hz_1D(2:end) = Hz_1D(2:end) + ty_mu*diff(Ey_1D);
    % update magnetic boundaries
    Hz_1D(1) = Hz_1D_prev_2 + mur_boundary_const_first*(Hz_1D(2) - Hz_1D_prev_1);
    
    
    %% update vj (auxilliary vectors)
    for i = 1:M
        % vj Px update
        vec = (a_gpof_P(:,i).*exp(b_gpof_P(:,i)));
        expanded_vec = reshape(vec, [1, 1, numel(vec)]);
        Px_dielectric = Px(y1:y2,x1:x2) .* dielectric(y1:y2,x1:x2);
        Px_gpof_3d = expanded_vec.*Px_dielectric;
        
        vec2 = exp(b_gpof_P(:,i));
        expanded_vec2 = reshape(vec2, [1, 1, numel(vec2)]);
        vj_Px_gpof_3d = expanded_vec2.*vj_Px_prev(y1:y2,x1:x2, 1+gpof_max*(i-1):gpof_max*i);
        
        vj_Px(y1:y2,x1:x2, 1+gpof_max*(i-1):gpof_max*i) = Px_gpof_3d + vj_Px_gpof_3d;
        vj_Px_sum(:,:,i) = sum(vj_Px(:,:,1+gpof_max*(i-1):gpof_max*i), 3);
        
        % vj Py update
        Py_dielectric = Py(y1:y2,x1:x2) .* dielectric(y1:y2,x1:x2);
        Py_gpof_3d = expanded_vec.*Py_dielectric;
        vj_Py_gpof_3d = expanded_vec2.*vj_Py_prev(y1:y2,x1:x2, 1+gpof_max*(i-1):gpof_max*i);

        vj_Py(y1:y2,x1:x2, 1+gpof_max*(i-1):gpof_max*i) = Py_gpof_3d + vj_Py_gpof_3d;
        vj_Py_sum(:,:,i) = sum(vj_Py(:,:,1+gpof_max*(i-1):gpof_max*i), 3);
    end
    vj_Px_prev = vj_Px;
    vj_Py_prev = vj_Py;
    
    % update vj for Ex and Ey
    for i = 1:N
        % vj Ex
        vec = (a_gpof_E(:,i).*exp(b_gpof_E(:,i)));
        expanded_vec = reshape(vec, [1, 1, numel(vec)]);
        Ex_dielectric = Ex(y1:y2,x1:x2) .* dielectric(y1:y2,x1:x2);
        Ex_gpof_3d = expanded_vec.*Ex_dielectric;
        
        vec2 = exp(b_gpof_E(:,i));
        expanded_vec2 = reshape(vec2, [1, 1, numel(vec2)]);
        vj_Ex_gpof_3d = expanded_vec2.*vj_Ex_prev(y1:y2,x1:x2, 1+gpof_max*(i-1):gpof_max*i);
        
        vj_Ex(y1:y2,x1:x2, 1+gpof_max*(i-1):gpof_max*i) = Ex_gpof_3d + vj_Ex_gpof_3d;
        vj_Ex_sum(:,:,i) = sum(vj_Ex(:,:,1+gpof_max*(i-1):gpof_max*i), 3);
        
        % vj Ey
        Ey_dielectric = Ey(y1:y2,x1:x2) .* dielectric(y1:y2,x1:x2);
        Ey_gpof_3d = expanded_vec.*Ey_dielectric;
        vj_Ey_gpof_3d = expanded_vec2.*vj_Ey_prev(y1:y2,x1:x2, 1+gpof_max*(i-1):gpof_max*i);

        vj_Ey(y1:y2,x1:x2, 1+gpof_max*(i-1):gpof_max*i) = Ey_gpof_3d + vj_Ey_gpof_3d;
        vj_Ey_sum(:,:,i) = sum(vj_Ey(:,:,1+gpof_max*(i-1):gpof_max*i), 3);
    end
    vj_Ex_prev = vj_Ex;
    vj_Ey_prev = vj_Ey;
    
    beta_vec = b./(dt.^beta);
    expanded_beta_vec = reshape(beta_vec, [1, 1, numel(beta_vec)]);
    alpha_vec = a./(dt.^alpha);
    expanded_alpha_vec = reshape(alpha_vec, [1, 1, numel(alpha_vec)]);
    
    psi_lPx = expanded_beta_vec .* vj_Px_sum;
    psi_lPy = expanded_beta_vec .* vj_Py_sum;
    psi_kEx = expanded_alpha_vec .* vj_Ex_sum;
    psi_kEy = expanded_alpha_vec .* vj_Ey_sum;
    
    psi_lPx_sum = sum(psi_lPx,3);
    psi_lPy_sum = sum(psi_lPy,3);
    psi_kEx_sum = sum(psi_kEx,3);
    psi_kEy_sum = sum(psi_kEy,3);
    
    %% update P (polarisation)
    Px(:,1:imax-1) = 1/(1 + C1 + C3/dt) * ( C3*( (Hz(1:imax,2:imax) - Hz(1:imax,1:imax-1))/dx - C4*Ex(:,1:imax-1) + Px(:,1:imax-1)/dt ) - psi_lPx_sum(:,1:imax-1) + eps0*eps_rs*psi_kEx_sum(:,1:imax-1) );
    Py(1:imax-1,:) = 1/(1 + C1 + C3/dt) * ( C3*( -(Hz(2:imax,1:imax) - Hz(1:imax-1,1:imax))/dy - C4*Ey(1:imax-1,:) + Py(1:imax-1,:)/dt ) - psi_lPy_sum(1:imax-1,:) + eps0*eps_rs*psi_kEy_sum(1:imax-1,:) );
    
    %% magnetic field update
    Hzy(1:imax,1:imax-1) = DaHzy(1:imax,1:imax-1).*Hzy(1:imax,1:imax-1) + DbHzy(1:imax,1:imax-1).*(Ex(1:imax,2:imax) - Ex(1:imax,1:imax-1));%diff(Ex(1:imax-1,:),1,2);
    Hzx(1:imax-1,1:imax) = DaHzx(1:imax-1,1:imax).*Hzx(1:imax-1,1:imax) - DbHzx(1:imax-1,1:imax).*(Ey(2:imax,1:imax) - Ey(1:imax-1,1:imax));%diff(Ey(:,1:imax-1),1,1);
    
    Hzx(tfsf_pos,tfsf_left:tfsf_right) = Hzx(tfsf_pos,tfsf_left:tfsf_right) + (dt/(mu0*dx)).*Ey_1D(source_pos);
    Hzx(tfsf_end,tfsf_left:tfsf_right) = Hzx(tfsf_end,tfsf_left:tfsf_right) - (dt/(mu0*dx)).*Ey_1D(tfsf_end-tfsf_pos+source_pos);
    
    Hz = Hzx + Hzy;
    
    %% electric field update
    Ex(1:imax,2:imax) = CaEx(1:imax,2:imax).*Ex(1:imax,2:imax) + CbEx(1:imax,2:imax).*(Hz(1:imax,2:imax) - Hz(1:imax,1:imax-1));%diff(Hz(2:imax,:),1,2);
    Ey(2:imax,1:imax) = CaEy(2:imax,1:imax).*Ey(2:imax,1:imax) - CbEy(2:imax,1:imax).*(Hz(2:imax,1:imax) - Hz(1:imax-1,1:imax));%diff(Hz(:,2:imax),1,1);
    
    %% simple dielectric 
%     Ex(1:imax,2:imax) = CaEx(1:imax,2:imax).*Ex(1:imax,2:imax) + CbEx(1:imax,2:imax).*(Hz(1:imax,2:imax) - Hz(1:imax,1:imax-1))./simple_eps(1:imax,2:imax);
%     Ey(2:imax,1:imax) = CaEy(2:imax,1:imax).*Ey(2:imax,1:imax) - CbEy(2:imax,1:imax).*(Hz(2:imax,1:imax) - Hz(1:imax-1,1:imax))./simple_eps(2:imax,1:imax);

    %% complex media electric update
    ex_y1 = y1;
    ex_y2 = y2;
    ex_x1 = x1+1;
    ex_x2 = x2;
    
    ey_y1 = y1+1;
    ey_y2 = y2;
    ey_x1 = x1;
    ey_x2 = x2;
    
    for i = ex_x1+1:ex_x2
        for j = ex_y1:ex_y2
            if ((i-center)^2 + (j-center)^2) <= radius^2
                Ex(j,i) = 1/C5 .*( (Hz(j,i) - Hz(j,i-1))/dx - C4*Ex_prev(j,i) - (Px(j,i) - Px_prev(j,i))/dt );
            else
                Ex(j,i) = CaEx(j,i).*Ex_prev(j,i) + CbEx(j,i).*(Hz(j,i) - Hz(j,i-1));
            end
        end
    end
    for i = ey_x1:ey_x2
        for j = ey_y1+1:ey_y2
            if ((i-center)^2 + (j-center)^2) <= radius^2
                Ey(j,i) = 1/C5 .*( -(Hz(j,i) - Hz(j-1,i))/dy - C4*Ey_prev(j,i) - (Py(j,i) - Py_prev(j,i))/dt );
            else
                Ey(j,i) = CaEy(j,i).*Ey_prev(j,i) - CbEy(j,i).*(Hz(j,i) - Hz(j-1,i));
            end
        end
    end

%     Ex(ex_y1:ex_y2,ex_x1+1:ex_x2) = 1/C5 .*( (Hz(ex_y1:ex_y2,ex_x1+1:ex_x2) - Hz(ex_y1:ex_y2,ex_x1:ex_x2-1))/dx - C4*Ex_prev(ex_y1:ex_y2,ex_x1+1:ex_x2) - (Px(ex_y1:ex_y2,ex_x1+1:ex_x2) - Px_prev(ex_y1:ex_y2,ex_x1+1:ex_x2))/dt );
%     Ey(ey_y1+1:ey_y2,ey_x1:ey_x2) = 1/C5 .*( -(Hz(ey_y1+1:ey_y2,ey_x1:ey_x2) - Hz(ey_y1:ey_y2-1,ey_x1:ey_x2))/dy - C4*Ey_prev(ey_y1+1:ey_y2,ey_x1:ey_x2) - (Py(ey_y1+1:ey_y2,ey_x1:ey_x2) - Py_prev(ey_y1+1:ey_y2,ey_x1:ey_x2))/dt );
%     Ex(ex_y1:ex_y2,ex_x1:ex_x2) = CaEx(ex_y1:ex_y2,ex_x1:ex_x2).*Ex_prev(ex_y1:ex_y2,ex_x1:ex_x2) + 1/eps_inf * CbEx(ex_y1:ex_y2,ex_x1:ex_x2).*(Hz(ex_y1:ex_y2,ex_x1:ex_x2) - Hz(ex_y1:ex_y2,ex_x1-1:ex_x2-1));
%     Ey(ey_y1:ey_y2,ey_x1:ey_x2) = CaEy(ey_y1:ey_y2,ey_x1:ey_x2).*Ey_prev(ey_y1:ey_y2,ey_x1:ey_x2) - 1/eps_inf * CbEy(ey_y1:ey_y2,ey_x1:ey_x2).*(Hz(ey_y1:ey_y2,ey_x1:ey_x2) - Hz(ey_y1-1:ey_y2-1,ey_x1:ey_x2));
    
    
    %% electric field TFSF
    tfsf_pos_Ex = tfsf_pos;
    tfsf_end_Ex = tfsf_end-1;
    
    % TFSF offset for complex dielectric
    tfsf_complex_offset = ones(tfsf_end_Ex-tfsf_pos_Ex+1,1);
    tfsf_complex_offset(y1-tfsf_pos_Ex+1:y2-tfsf_pos_Ex+1) = tfsf_complex_offset(y1-tfsf_pos_Ex+1:y2-tfsf_pos_Ex+1) ;%.* 1/eps_inf;
    
    Hz_correct_dimensions = reshape(Hz_1D(source_pos+1:tfsf_end_Ex+source_pos-tfsf_pos_Ex+1), length(tfsf_pos_Ex:tfsf_end_Ex), 1);
    Ex(tfsf_pos_Ex:tfsf_end_Ex,tfsf_left) = Ex(tfsf_pos_Ex:tfsf_end_Ex,tfsf_left) + tfsf_complex_offset.*CbEx(tfsf_pos_Ex:tfsf_end_Ex,tfsf_left).*Hz_correct_dimensions;
    Ex(tfsf_pos_Ex:tfsf_end_Ex,tfsf_right+1) = Ex(tfsf_pos_Ex:tfsf_end_Ex,tfsf_right+1) - tfsf_complex_offset.*CbEx(tfsf_pos_Ex:tfsf_end_Ex,tfsf_right+1).*Hz_correct_dimensions;
    
    Ey(tfsf_pos,tfsf_left:tfsf_right) = Ey(tfsf_pos,tfsf_left:tfsf_right) + CbEy(tfsf_pos,tfsf_left:tfsf_right).* Hz_1D(source_pos);
    Ey(tfsf_end,tfsf_left:tfsf_right) = Ey(tfsf_end,tfsf_left:tfsf_right) + CbEy(tfsf_end,tfsf_left:tfsf_right).* Hz_1D(tfsf_end-tfsf_pos+source_pos+1);

%     %% magnetic field update
%     Hzy(1:imax,1:imax-1) = DaHzy(1:imax,1:imax-1).*Hzy(1:imax,1:imax-1) + DbHzy(1:imax,1:imax-1).*(Ex(1:imax,2:imax) - Ex(1:imax,1:imax-1));%diff(Ex(1:imax-1,:),1,2);
%     Hzx(1:imax-1,1:imax) = DaHzx(1:imax-1,1:imax).*Hzx(1:imax-1,1:imax) - DbHzx(1:imax-1,1:imax).*(Ey(2:imax,1:imax) - Ey(1:imax-1,1:imax));%diff(Ey(:,1:imax-1),1,1);
%     
%     tfsf_pos_Hzx = tfsf_pos;
%     Hzx(tfsf_pos_Hzx,tfsf_left:tfsf_right) = Hzx(tfsf_pos_Hzx,tfsf_left:tfsf_right) + (dt/(mu0*dx)).*Ey_1D(source_pos);
%     Hzx(tfsf_end,tfsf_left:tfsf_right) = Hzx(tfsf_end,tfsf_left:tfsf_right) - (dt/(mu0*dx)).*Ey_1D(tfsf_end-tfsf_pos_Hzx+source_pos);
%     
%     Hz = Hzx + Hzy;
    
    %% update P_prev and E_prev
    Px_prev = Px;
    Py_prev = Py;
    
    Ex_prev = Ex;
    Ey_prev = Ey;
    
    %% store fields
%     Ex_trans(n) = Ey_1D(source_pos);%Ey(posy,posx);
%     Ex_pos2(n) = Ey(posy2,posx2);
%     Ex_1D(n) = Ey_1D();
%     H_tfsf(:,:,n) = Hz(tfsf_pos:tfsf_end,tfsf_left:tfsf_right);
%     Ex_tfsf(:,:,n) = Ex(tfsf_pos:tfsf_end,tfsf_left:tfsf_right);
%     Ey_tfsf(:,:,n) = Ey(tfsf_pos:tfsf_end,tfsf_left:tfsf_right);
    if n == 1000
        Hz_ss = Hz;
    end
    if n == 5000
        Hz_ss_2pi = Hz;
    end
    if n == 2000
        Hz_ss_quarterpi = Hz;
    end
    if n == 3000
        Hz_ss_pi = Hz;
    end
    if n == 4000
        Hz_ss_pi_3quarter = Hz;
    end
    
    %% plotting
    
    if mod(n,20) == 0
        figure(1)
%         surf(x,y,real(Ey));
        clims = [-2,2];
        imagesc(x,y,real(Hz),clims);
        colormap('jet')
        colorbar
        
        xline(PML_thickness)
        xline(imax-PML_thickness)
        yline(PML_thickness)
        yline(imax-PML_thickness)
        viscircles([center,center],radius,'Color','k');
        xlabel('$\Delta x$', 'Interpreter', 'latex')
        ylabel('$\Delta y$', 'Interpreter', 'latex')
        set(gca,'YDir','normal')
        title('$\Re e\{H_z\}$', 'Interpreter', 'latex')
        frames(frame) = getframe(fig);
        pause(0.1);
        frame = frame + 1;
    end
% 
%     if n == 1000
%         figure(2)
% %         surf(x,y,real(Ey));
%         clims = [-0.1,0.1];
%         imagesc(x,y,real(Ey));
%         colormap('jet')
%         colorbar
%         
%         xline(PML_thickness)
%         xline(imax-PML_thickness)
%         yline(PML_thickness)
%         yline(imax-PML_thickness)
%         rectangle('Position', [media_start_x,media_start_y,media_end_x-media_start_x,media_end_y-media_start_y])
%         xlabel('x')
%         ylabel('y')
%         getframe;
%     end
    
    disp(n)
end
toc
%%
% Create an animated GIF from the captured frames
filename = '2D_FDTD_animation_slower_2.gif';
for frame = 1:numFrames
    im = frame2im(frames(frame));
    [imind, cm] = rgb2ind(im, 256);
    if frame == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%
% % Post-processing
% %%%%%%%%%%%%%%%%%%%%%%%%%%


% figure(1)
% clims = [-0.1,0.1];
% imagesc(x(1:300),y(1:300),real(Hz_ss_pi_3quarter(tfsf_pos:tfsf_end-1,tfsf_left:tfsf_right-1)));
% colormap('jet')
% colorbar
% 
% % xline(PML_thickness)
% % xline(imax-PML_thickness)
% % yline(PML_thickness)
% % yline(imax-PML_thickness)
% viscircles([150,150],radius,'Color','k');
% xlabel('x')
% ylabel('y')
% 
% 
% %%
% 
% % figure(1)
% % plot(Ex_trans)
% 
% close all;
% 
% figure(1)
% imagesc(x,y,real(Hz_ss));
% % colormap('jet')
% colorbar
% xline(PML_thickness)
% xline(imax-PML_thickness)
% yline(PML_thickness)
% yline(imax-PML_thickness)
% % viscircles([center,center],radius,'Color','k');
% xlabel('x')
% ylabel('y')
% 
% figure(2)
% imagesc(x,y,real(Hz_ss_2pi));
% % colormap('jet')
% colorbar
% xline(PML_thickness)
% xline(imax-PML_thickness)
% yline(PML_thickness)
% yline(imax-PML_thickness)
% % viscircles([center,center],radius,'Color','k');
% xlabel('x')
% ylabel('y')
% 
% figure(3)
% imagesc(x,y,real(Hz_ss_quarterpi));
% % colormap('jet')
% colorbar
% xline(PML_thickness)
% xline(imax-PML_thickness)
% yline(PML_thickness)
% yline(imax-PML_thickness)
% % viscircles([center,center],radius,'Color','k');
% xlabel('x')
% ylabel('y')
% 
% error = Hz_ss_2pi - Hz_ss;
% figure(5)
% imagesc(x,y,abs(error));
% colorbar
% 
% abs_Hz = abs(sqrt(Hz_ss_quarterpi.^2 + Hz_ss.^2));
% figure(6)
% imagesc(x,y,abs_Hz);
% colorbar
% xline(PML_thickness)
% xline(imax-PML_thickness)
% yline(PML_thickness)
% yline(imax-PML_thickness)
% % viscircles([center,center],radius,'Color','k');
% xlabel('x')
% ylabel('y')
% 
% abs_Hz = abs(sqrt(Hz_ss_quarterpi.^2 + Hz_ss_pi.^2));
% figure(7)
% imagesc(x,y,abs_Hz);
% colorbar
% xline(PML_thickness)
% xline(imax-PML_thickness)
% yline(PML_thickness)
% yline(imax-PML_thickness)
% % viscircles([center,center],radius,'Color','k');
% xlabel('x')
% ylabel('y')
% 
% abs_Hz = abs(sqrt(Hz_ss_pi_3quarter.^2 + Hz_ss_pi.^2));
% figure(8)
% imagesc(x,y,abs_Hz);
% colorbar
% xline(PML_thickness)
% xline(imax-PML_thickness)
% yline(PML_thickness)
% yline(imax-PML_thickness)
% % viscircles([center,center],radius,'Color','k');
% xlabel('x')
% ylabel('y')
% 
% abs_Hz = abs(sqrt(Hz_ss_pi_3quarter.^2 + Hz_ss_2pi.^2));
% abs_Hz_correct = abs_Hz(tfsf_pos:tfsf_end-1,tfsf_left:tfsf_right-1)./1.66666;%(5/3);
% figure(9)
% imagesc(x(1:300),y(1:300),abs_Hz_correct);
% colorbar
% xlabel('x')
% ylabel('y')
% title("FDTD abs Hz")
% 
% analytical_abs_Hz = readmatrix("Hz.txt");
% figure(11)
% imagesc(x(1:300),y(1:300),analytical_abs_Hz);
% colorbar
% xlabel('x')
% ylabel('y')
% title("analytical abs Hz")
% 
% error_relative = (abs_Hz_correct - analytical_abs_Hz)./analytical_abs_Hz;
% figure(12)
% imagesc(x(1:300),y(1:300),error_relative);
% colorbar
% xlabel('$\Delta x$', 'Interpreter', 'latex')
% ylabel('$\Delta y$', 'Interpreter', 'latex')
% % title("relaltive error")
% 
% 
% figure(4)
% imagesc(x,y,real(Hz_ss_pi));
% % colormap('jet')
% colorbar
% xline(PML_thickness)
% xline(imax-PML_thickness)
% yline(PML_thickness)
% yline(imax-PML_thickness)
% % viscircles([center,center],radius,'Color','k');
% xlabel('x')
% ylabel('y')
% 
% figure(10)
% imagesc(x,y,real(Hz_ss_pi_3quarter));
% % colormap('jet')
% colorbar
% xline(PML_thickness)
% xline(imax-PML_thickness)
% yline(PML_thickness)
% yline(imax-PML_thickness)
% % viscircles([center,center],radius,'Color','k');
% xlabel('x')
% ylabel('y')
% 
% %% 
% average_error = mean(error_relative, 'all');
% 
% %%
% writematrix(error_relative, "relative_error_1111MHz.txt");
% writematrix(analytical_abs_Hz, "analytical_abs_Hz_1111MHz.txt");
% writematrix(abs_Hz_correct, "abs_Hz_correct_1111MHz.txt");
% writematrix(average_error, "average_error_1111MHz.txt");
% 
% %%
% f = (10/9)*1e9;
% omega = 2*pi*f;
% eps_r = eps_inf + eps_rs*(1 + a(1)*(1i*omega).^alpha(1))./(1 + b(1)*(1i*omega).^beta(1) + b(2)*(1i*omega).^beta(2) + b(3)*(1i*omega).^beta(3));
% 
% % %%
% % close all;
% % 
% % abs_Hz_correct = readmatrix("abs_Hz_correct_2000MHz.txt");
% % figure(9)
% % imagesc(x(1:300),(y(1:300)),abs_Hz_correct);
% % colorbar
% % colormap('jet')
% % xlabel('$\Delta x$', 'Interpreter', 'latex')
% % ylabel('$\Delta y$', 'Interpreter', 'latex')
% % set(gca,'YDir','normal')
% % % title("FDTD abs Hz")
% % 
% % analytical_abs_Hz = readmatrix("analytical_abs_Hz_2000MHz.txt");
% % figure(11)
% % imagesc(x(1:300),(y(1:300)),analytical_abs_Hz);
% % colorbar
% % colormap('jet')
% % xlabel('$\Delta x$', 'Interpreter', 'latex')
% % ylabel('$\Delta y$', 'Interpreter', 'latex')
% % set(gca,'YDir','normal')
% % % title("analytical abs Hz")
% % 
% % error_relative = readmatrix("relative_error_2000MHz.txt");
% % figure(12)
% % imagesc(x(1:300),y(1:300),abs(error_relative));
% % colorbar
% % colormap('jet')
% % xlabel('$\Delta x$', 'Interpreter', 'latex')
% % ylabel('$\Delta y$', 'Interpreter', 'latex')
% % set(gca,'YDir','normal')
% % % title("relaltive error")
% 
% 
% %%
% close all;
% 
% h1 = readmatrix("abs_Hz_correct_500MHz.txt");
% h2 = readmatrix("abs_Hz_correct_667MHz.txt");
% h3 = readmatrix("abs_Hz_correct_833MHz.txt");
% h4 = readmatrix("abs_Hz_correct_1000MHz.txt");
% h5 = readmatrix("abs_Hz_correct_1111MHz.txt");
% h6 = readmatrix("abs_Hz_correct_1250MHz.txt");
% h7 = readmatrix("abs_Hz_correct_1429MHz.txt");
% h8 = readmatrix("abs_Hz_correct_1667MHz.txt");
% h9 = readmatrix("abs_Hz_correct_2000MHz.txt");
% 
% c1 = readmatrix("analytical_abs_Hz_500MHz.txt");
% c2 = readmatrix("analytical_abs_Hz_667MHz.txt");
% c3 = readmatrix("analytical_abs_Hz_833MHz.txt");
% c4 = readmatrix("analytical_abs_Hz_1000MHz.txt");
% c5 = readmatrix("analytical_abs_Hz_1111MHz.txt");
% c6 = readmatrix("analytical_abs_Hz_1250MHz.txt");
% c7 = readmatrix("analytical_abs_Hz_1429MHz.txt");
% c8 = readmatrix("analytical_abs_Hz_1667MHz.txt");
% c9 = readmatrix("analytical_abs_Hz_2000MHz.txt");
% 
% a1 = sum(abs(h1 - c1), 'all')/sum(c1, 'all');
% a2 = sum(abs(h2 - c2), 'all')/sum(c2, 'all');
% a3 = sum(abs(h3 - c3), 'all')/sum(c3, 'all');
% a4 = sum(abs(h4 - c4), 'all')/sum(c4, 'all');
% a5 = sum(abs(h5 - c5), 'all')/sum(c5, 'all');
% a6 = sum(abs(h6 - c6), 'all')/sum(c6, 'all');
% a7 = sum(abs(h7 - c7), 'all')/sum(c7, 'all');
% a8 = sum(abs(h8 - c8), 'all')/sum(c8, 'all');
% a9 = sum(abs(h9 - c9), 'all')/sum(c9, 'all');
% 
% err1 = readmatrix("relative_error_500MHz.txt");
% err2 = readmatrix("relative_error_667MHz.txt");
% err3 = readmatrix("relative_error_833MHz.txt");
% err4 = readmatrix("relative_error_1000MHz.txt");
% err5 = readmatrix("relative_error_1111MHz.txt");
% err6 = readmatrix("relative_error_1250MHz.txt");
% err7 = readmatrix("relative_error_1429MHz.txt");
% err8 = readmatrix("relative_error_1667MHz.txt");
% err9 = readmatrix("relative_error_2000MHz.txt");
% 
% max_err1 = max(abs(err1), [], 'all');
% max_err2 = max(abs(err2), [], 'all');
% max_err3 = max(abs(err3), [], 'all');
% max_err4 = max(abs(err4), [], 'all');
% max_err5 = max(abs(err5), [], 'all');
% max_err6 = max(abs(err6), [], 'all');
% max_err7 = max(abs(err7), [], 'all');
% max_err8 = max(abs(err8), [], 'all');
% max_err9 = max(abs(err9), [], 'all');
% max_err = [max_err1, max_err2, max_err3, max_err4, max_err5, max_err6, max_err7, max_err8, max_err9];
% 
% min_err1 = min(abs(err1), [], 'all');
% min_err2 = min(abs(err2), [], 'all');
% min_err3 = min(abs(err3), [], 'all');
% min_err4 = min(abs(err4), [], 'all');
% min_err5 = min(abs(err5), [], 'all');
% 
% % a1 = mean(abs(err1), 'all');
% % a2 = mean(abs(err2), 'all');
% % a3 = mean(abs(err3), 'all');
% % a4 = mean(abs(err4), 'all');
% % a5 = mean(abs(err5), 'all');
% % a6 = mean(abs(err6), 'all');
% % a7 = mean(abs(err7), 'all');
% % a8 = mean(abs(err8), 'all');
% % a9 = mean(abs(err9), 'all');
% 
% freq = [0.5, 0.667, 0.833, 1, 1.11, 1.25, 1.429, 1.667, 2];  % GHz
% average_err = [a1, a2, a3, a4, a5, a6, a7, a8, a9];
% 
% figure(1)
% loglog(freq, abs(average_err), '.-')
% hold on
% loglog(freq, max_err, 'x-')
% grid on
% ylabel("Relative Error", 'Interpreter', 'latex')
% xlabel("Frequnecy [GHz]", 'Interpreter', 'latex')
% legend("Averaged $RE$", "Maximum $RE$", 'Location', 'northwest', 'Interpreter', 'latex')
% ylim([0.0001 1])


% figure(2)
% imagesc(x(1:300),(y(1:300)),analytical_abs_Hz);
% colorbar
% colormap('jet')


%%
% f = (0.01:0.01:100)*1e9;
% omega = 2*pi*f;
% 
% dz = dx;
% 
% eps_r = eps_inf.*ones(1,length(omega));
% 
% if length(beta) == 3
%     % eps_r = eps_inf + eps_rs*(1)./(1 + b(1)*(1i*omega).^beta(1) + b(2)*(1i*omega).^beta(2) + b(3)*(1i*omega).^beta(3));
%     eps_r = eps_inf + eps_rs*(1 + a(1)*(1i*omega).^alpha(1))./(1 + b(1)*(1i*omega).^beta(1) + b(2)*(1i*omega).^beta(2) + b(3)*(1i*omega).^beta(3));
% end
% if length(beta) == 1
%     eps_r = eps_inf + eps_rs*(1 + a(1)*(1i*omega).^alpha(1))./(1 + b(1)*(1i*omega).^beta(1));
% end
% if length(beta) == 2
%     eps_r = eps_inf + eps_rs*(1 + a(1)*(1i*omega).^alpha(1))./(1 + b(1)*(1i*omega).^beta(1) + b(2)*(1i*omega).^beta(2));
% end
% 
% d = (y2 - y1)*dy;
% 
% T_analytical = exp(-1i*omega.*d.*sqrt(eps_r)/c0) ;
% 
% % fabry-perot comparison
% n1 = 1; n2 = sqrt(eps_r);
% k = omega.*sqrt(eps_r)./c0;
% k0 = omega./c0;
% 
% i_media_start = media_start_y;
% i_media_stop = media_end_y;
% isource = tfsf_pos;
% 
% d1 = 0;
% d2 = 0;
% dm = i_media_stop - i_media_start;
% ds = i_media_start - isource;
% 
% r_air = (n1 - n2)./(n1 + n2);
% r_media = -r_air;
% t_squared = (4*n1*n2./(n1 + n2).^2);
% power_sum = 1./(1 - r_media.^2 .* exp(-1i*2*k*dm*dz));
% numer = exp(-1i*k0*(ds + d2)*dz) .* exp(-1i*k*dm*dz) .* t_squared .* power_sum;
% denom = exp(-1i*k0*(ds - d1)*dz) + r_air .* exp(-1i*k0*(ds + d1)*dz) .* (1 - t_squared.*exp(-1i*k*2*dm*dz).*power_sum);
% T_fp = numer./denom;
% 
% % Ex_source = zeros(1, nmax);
% % for n = 2:nmax-1
% %     Ex_source(n) = source_function(n,dt);
% % end
% % Ex_source_freq = fft(Ex_source);
% 
% Ex_pos2_freq = fft(Ex_pos2);
% 
% Ex_pos1 = Ex_trans;
% Ex_pos1_freq = fft(Ex_pos1);
% 
% T = Ex_pos1_freq./Ex_pos2_freq;
% 
% % pos = 1000;
% % E1 = [Ex_pos1(1:pos), zeros(1,nmax-pos)];
% % E2 = [Ex_pos2(1:pos), zeros(1,nmax-pos)];
% % E1_fft = fft(E1);
% % E2_fft = fft(E2);
% % T2 = E1_fft./E2_fft;
% 
% fff = (1/dt)*(0:nmax-1)/nmax;
% figure(6)
% semilogx(fff*1e-9, real(T), 'bx')
% hold on
% semilogx(fff*1e-9, imag(T), 'r+')
% semilogx(f*1e-9, real(T_fp), 'b')
% semilogx(f*1e-9, imag(T_fp), 'r')
% 
% % semilogx(fff*1e-9, real(T2), 'b+')
% % semilogx(fff*1e-9, imag(T2), 'rx')
% 
% hold off
% xlim([0.1 10])
% % ylim([-1 1])
% legend("Re\{T\} FDTD", "Im\{T\} FDTD", "Re\{T\} Analytical", "Im\{T\} Analytical", 'Location', 'northeast', 'Interpreter', 'latex')
% xlabel("Frequnecy [GHz]")
% title("$\epsilon_r(\omega) = 3$", 'Interpreter', 'latex')
% % title("$\epsilon_r(\omega) = 2 + 60 \frac{1 + (j\omega \tau)^{0.2}}{1 + 9(j\omega \tau)^{0.3} + 2(j\omega \tau)^{0.5} + 10(j\omega \tau)^{0.9}}$", 'Interpreter', 'latex')
% grid minor





%% function definitions

function wave_source = source_function(t, dt)
    fc = 6e9;
    fbw = 5/6;
    [~,yq] = gauspuls((t)*dt, fc, fbw, -3);
    
    f = (5)*1e9;
    eps0  = 8.8541878128e-12;
    mu0 = 1.256637062e-6;
    E0 = sqrt(mu0/eps0);

    
    wave_source = E0*sin(2*pi*f*(t-1)*dt);% yq;
end

function [Bt, At]=gpof(f,t,M,L)
    N=length(f);
    
    if(t(1)~=0)
        disp('t-vector must start at 0!')
        return
    end
    if exist('L')==0
        L=round(N/2);
        fprintf('Choosing default pencil parameter N/2=%d...\n',L);
    end
    if (M>=L)
        fprintf('Requested order %d is larger than pencil parameter %d\n',M,L)
        M=L-1;
        fprintf('Choosing new order %d...\n',M)
    end
    if (L>N-M)
        fprintf('Number of support points too small for order %d\n',N)
        return
    end
    
    Y1=zeros(N-L,L);
    for rw=0:N-L-1
        for cl=0:L-1
            Y1(rw+1,cl+1)=f(rw+cl+1);
        end
    end
    
    Y2=zeros(N-L,L);
    for rw=1:N-L
        for cl=0:L-1
            Y2(rw,cl+1)=f(rw+cl+1);
        end
    end
    
    [U, D, V]=svd(Y1);
    invD=pinv(D);
    Z=invD*U'*Y2*V;
    Z=Z(1:M,1:M);
    zitemp=eig(Z);
    zi=zitemp(1:M);
    
    rankdef=sum(zitemp==0);
    if(rankdef~=0)
        fprintf('Order of approximation %d too high by %d (rank %d)\n', M,rankdef,M-rankdef);
    end
    
    dt=(t(2)-t(1));
    Bt=log(zi)/dt;
    
    Zmat=zeros(M,N);
    for rw=1:M
        for cl=0:N-1
            Zmat(rw,cl+1)=zi(rw)^cl;
        end
    end
    
    At=(f*pinv(Zmat)).';
    return
end
