%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Aneesh Deogan
% Date:   09/06/2023
% Script: 1D FDTD scheme for frequency dispersive media
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clean workspace
close all; clear variables; clc;

%% constants
f = (0.01:0.01:100)*1e9;
eps0  = 8.8541878128e-12;
mu0 = 1.256637062e-6;
c0 = 1/sqrt(eps0*mu0);
imp0 = sqrt(mu0/eps0);
omega = 2*pi*f;
tau = 318e-12;
eps_s = 2;%50;
eps_inf = 2;%2
d_eps = eps_s - eps_inf;

%% FDTD parameters
off = 0.5;
dt = 1.768e-12/off;
dz = 1.1e-3/off;
imax = 500*off;  % max discretization in space
nmax = 10000*off; % max discretization in time
isource = 200*off;
i_media_start = 300*off;
i_media_stop = 330*off;

%% dielectric permittivity parameters

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



% a_gpof_E = [1.917075539048859, -0.560287204184830, -0.198639691967215, -0.085752446603181, -0.045331304427383, -0.018930049994618, -0.006964496770346, -0.001170344468412]';
% b_gpof_E = [-5.347160450874707 + 3.141592653589793i, -2.579206230697158 + 0.00i, -1.356262901949623 + 0.00i, -0.725900394409128 + 0.00i, -0.359939622406301 + 0.00i, -0.151364238040643 + 0.00i, -0.045973859900755 + 0.00i, -0.006388607475665 + 0.00i]';
% a_gpof_P1 = [2.669664779268722, -1.130353379179308, -0.319124863825985, -0.127427146917053, -0.060073567274821, -0.023621097952209, -0.007853612584541, -0.001211110983464]';
% b_gpof_P1 = [-5.016172362040061 + 3.141592653589793i, -2.708543106502772 + 0.0000i, -1.413642612584595 + 0.0000i, -0.760037597871560 + 0.0000i, -0.380999711174307 + 0.0000i, -0.163439794019462 + 0.0000i, -0.051602128858320 + 0.0000i, -0.007889929899061 + 0.00i]';
% a_gpof_P2 = [5.331832122564569, -3.463142726299982, -0.575971123489243, -0.187842028805195, -0.072239715677446, -0.024849262526770, -0.006883150744650, -9.041148289020017e-04]';
% b_gpof_P2 = [-4.682192687845149 + 3.141592653589793i, -3.015762892824395 + 0.00i, -1.535975941329309 + 0.00i, -0.831468861641271 + 0.00i, -0.425035368945261 + 0.00i, -0.188985508593928 + 0.00i, -0.063884248848722 + 0.00i, -0.011465455385803 + 0.00i]';
% a_gpof_P3 = [34.711296294933700, -33.087263703948906, -0.496500152254662, -0.094861637251343, -0.024878882072776, -0.006381535486200, -0.001285233632307, -1.251502861917283e-04]';
% b_gpof_P3 = [-4.830293886144696 + 3.141592653589793i, -4.220728932009782 + 0.0000i, -1.819457273975490 + 0.0000i, -0.988060793726207 + 0.0000i, -0.520664832308133 + 0.0000i, -0.245289842399021 + 0.0000i, -0.092246963047317 + 0.0000i, -0.020870775550631 + 0.0000i]';
% a_gpof_P = [a_gpof_P1, a_gpof_P2, a_gpof_P3];
% b_gpof_P = [b_gpof_P1, b_gpof_P2, b_gpof_P3];


%% FDTD update parameters initialisation
vj_E = zeros(gpof_max*N, imax);
vj_P = zeros(gpof_max*M, imax);
Ex = zeros(nmax, imax);
P = zeros(nmax, imax);
Hy = zeros(nmax, imax-1);
vj_E_prev = zeros(gpof_max*N, imax);
vj_P_prev = zeros(gpof_max*M, imax);
Ex_prev = zeros(nmax, imax);
P_prev = zeros(1, imax);
Hy_prev = zeros(1, imax);
vj_P_sum = zeros(M,imax);
vj_E_sum = zeros(N,imax);

%% FDTD boundary conditions
mur_boundary_const_first = (c0*dt - dz)/(c0*dt + dz);

% PML
PML_thickness = 50*off;  % PML thickness
sigma = zeros(1, imax);
for i = 1:imax
    if i < PML_thickness
        sigma(i) = (2*eps0/(3*dt)) * ((PML_thickness-i)/PML_thickness)^3;
    elseif i > imax-PML_thickness
        sigma(i) = (2*eps0/(3*dt)) * ((i-(imax-PML_thickness))/PML_thickness)^3;
    end
end
PML_const_E = sigma*dt/(2*eps0*1);
PML_const_H = sigma/(eps0*1);

%% FDTD constants
C5 = conductivity/2 + eps_inf* eps0/dt;
C4 = conductivity/2 - eps_inf* eps0/dt;
C2 = sum(a./(dt.^alpha));
C1 = sum(b./(dt.^beta));
C3 = (eps0*eps_rs*(1 + C2)/C5);
tz_mu = dt/(mu0*dz);

%% FDTD air-dielectric-air 
tic

for n = 2:nmax
    
    % update vj (auxilliary vectors)
    for i = 1:M
        vj_P(1+gpof_max*(i-1):gpof_max*i,i_media_start:i_media_stop) = (a_gpof_P(:,i).*exp(b_gpof_P(:,i)))*P(n-1,i_media_start:i_media_stop) + exp(b_gpof_P(:,i)).*vj_P_prev(1+gpof_max*(i-1):gpof_max*i,i_media_start:i_media_stop);
        vj_P_sum(i,:) = sum(vj_P(1+gpof_max*(i-1):gpof_max*i,:));
    end
    vj_P_prev = vj_P;
    
    for i = 1:N
        vj_E(1+gpof_max*(i-1):gpof_max*i,i_media_start:i_media_stop) = (a_gpof_E(:,i).*exp(b_gpof_E(:,i)))*Ex(n-1,i_media_start:i_media_stop) + exp(b_gpof_E(:,i)).*vj_E_prev(1+gpof_max*(i-1):gpof_max*i,i_media_start:i_media_stop);
        vj_E_sum(i,:) = sum(vj_E(1+gpof_max*(i-1):gpof_max*i,:));
    end
    vj_E_prev = vj_E;
    
    psi_lP = (b./(dt.^beta))' .* vj_P_sum;
    psi_kE = (a./(dt.^alpha))' .* vj_E_sum;
    
    psi_lP_sum = sum(psi_lP,1);
    psi_kE_sum = sum(psi_kE,1);
    
    % update P (polarisation)
    P(n,2:end-1) = 1/(1 + C1 + C3/dt) * ( C3*( diff(Hy(n-1,:))/dz - C4*Ex(n-1,2:end-1) + P(n-1,2:end-1)/dt ) - psi_lP_sum(2:end-1) + eps0*eps_rs*psi_kE_sum(2:end-1) );
    
    % update E (electric field)    
    % air
    Ex(n,2:i_media_start-1) = Ex(n-1,2:i_media_start-1).*(1 - PML_const_E(2:i_media_start-1))./(1 + PML_const_E(2:i_media_start-1)) + dt/(dz*eps0) * diff(Hy(n-1,1:i_media_start-1))./(1 + PML_const_E(2:i_media_start-1));
    % media
    Ex(n,i_media_start:i_media_stop) = 1/C5 *( (diff(Hy(n-1,i_media_start-1:i_media_stop)))/dz - C4*Ex(n-1,i_media_start:i_media_stop) - (P(n,i_media_start:i_media_stop) - P(n-1,i_media_start:i_media_stop))/dt );
    % air
    Ex(n,i_media_stop+1:end-1) = Ex(n-1,i_media_stop+1:end-1).*(1 - PML_const_E(i_media_stop+1:end-1))./(1 + PML_const_E(i_media_stop+1:end-1)) + dt/(dz*eps0) * diff(Hy(n-1,i_media_stop:end))./(1 + PML_const_E(i_media_stop+1:end-1));

    % update electric boundaries
    Ex(n,1) = Ex(n-1,2) + mur_boundary_const_first*(Ex(n,2) - Ex(n-1,1));
    Ex(n,imax) = Ex(n-1,imax-1) + mur_boundary_const_first*(Ex(n,imax-1) - Ex(n-1,imax));
    % electric field source
    Ex(n,isource) = Ex(n,isource) + source_function(n, dt);
    
    % update H (magnetic field)
    Hy(n,:) = Hy(n-1,:).*(1 - PML_const_H(1:imax-1)*dt) + tz_mu*diff(Ex(n,:));
        
    if n == 500
        Ex_500 = Ex(n,:);
    end
    if n == 1000
        Ex_1000 = Ex(n,:);
    end
    if n == 1500
        Ex_1500 = Ex(n,:);
    end

%     % plot field
%     if mod(n,50) == 0
%         figure(1)
%         plot(real(Ex(n,:)));
%         hold on
%         xline(i_media_start);
%         xline(i_media_stop);
%         hold off
% %         ylim([-1, 1])
%         xlim([0 imax])
%         getframe;
%         
%     end
    disp(n)
end
toc

%%

% plot(real(Ex(2,:)));
% axis tight
% set(gca,'nextplot','replacechildren','visible','off')
% f = getframe(gcf);
% [im,map] = rgb2ind(f.cdata,256,'nodither');
% im(1,1,1,20) = 0;
% for k = 5:20:500
%   plot(real(Ex(k,:)));
%   f = getframe(gcf);
%   im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
% end
% imwrite(im,map,'DancingPeaks.gif','DelayTime',0,'LoopCount',inf) %g443800

%% post-processing - plot transmittance of electric field

d = (i_media_stop-i_media_start)*dz;
% d = 30*dz;

if length(beta) == 5
    % eps_r = eps_inf + eps_rs*(1)./(1 + b(1)*(1i*omega).^beta(1) + b(2)*(1i*omega).^beta(2) + b(3)*(1i*omega).^beta(3));
    eps_r = eps_inf + eps_rs*(1 + a(1)*(1i*omega).^alpha(1))./(1 + b(1)*(1i*omega).^beta(1) + b(2)*(1i*omega).^beta(2) + b(3)*(1i*omega).^beta(3) + b(4)*(1i*omega).^beta(4) + b(5)*(1i*omega).^beta(5));
end
if length(beta) == 3
    % eps_r = eps_inf + eps_rs*(1)./(1 + b(1)*(1i*omega).^beta(1) + b(2)*(1i*omega).^beta(2) + b(3)*(1i*omega).^beta(3));
    eps_r = eps_inf + eps_rs*(1 + a(1)*(1i*omega).^alpha(1))./(1 + b(1)*(1i*omega).^beta(1) + b(2)*(1i*omega).^beta(2) + b(3)*(1i*omega).^beta(3));
end
if length(beta) == 1
    eps_r = eps_inf + eps_rs*(1 + a(1)*(1i*omega).^alpha(1))./(1 + b(1)*(1i*omega).^beta(1));
end
if length(beta) == 2
    eps_r = eps_inf + eps_rs*(1 + a(1)*(1i*omega).^alpha(1))./(1 + b(1)*(1i*omega).^beta(1) + b(2)*(1i*omega).^beta(2));
end


T_analytical = exp(-1i*omega.*d.*sqrt(eps_r)./c0) ;

% (fabry-perot)
n1 = 1; n2 = sqrt(eps_r);
k = omega.*sqrt(eps_r)./c0;
k0 = omega./c0;

d1 = 0;
d2 = 0;
dm = i_media_stop - i_media_start;
ds = i_media_start - isource;

r_air = (n1 - n2)./(n1 + n2);
r_media = -r_air;
t_squared = (4*n1*n2./(n1 + n2).^2);
power_sum = 1./(1 - r_media.^2 .* exp(-1i*2*k*dm*dz));
numer = exp(-1i*k0*(ds + d2)*dz) .* exp(-1i*k*dm*dz) .* t_squared .* power_sum;
denom = exp(-1i*k0*(ds - d1)*dz) + r_air .* exp(-1i*k0*(ds + d1)*dz) .* (1 - t_squared.*exp(-1i*k*2*dm*dz).*power_sum);
T_fp = numer./denom;

pos1 = i_media_stop+d2+1;
Ex_pos1 = Ex(:, pos1);
Ex_pos1_freq = fft(Ex_pos1);

pos2 = i_media_start-d1+1;
Ex_pos2 = Ex(:, pos2);
Ex_pos2_freq = fft(Ex_pos2);


T = Ex_pos1_freq./Ex_pos2_freq;
fff = (1/dt)*(0:nmax-1)/nmax;


figure(6)
semilogx(fff*1e-9, real(T), 'bx')
hold on
semilogx(fff*1e-9, imag(T), 'rx')

% semilogx(ff4*1e-9, real(T_c4), 'bo')
% semilogx(ff4*1e-9, imag(T_c4), 'ro')
% semilogx(ff1*1e-9, real(T_c1), 'b+')
% semilogx(ff1*1e-9, imag(T_c1), 'r+')

% semilogx(fff*1e-9, abs(T).^2, 'c.-')
% semilogx(f*1e-9, abs(T_fp).^2, 'g')
semilogx(f*1e-9, real(T_fp), 'b')
semilogx(f*1e-9, imag(T_fp), 'r')

xlim([0.1 10])
% ylim([-1 1])

% legend("Re\{T\} FDTD: C = 2", "Im\{T\} FDTD: C = 2", "Re\{T\} FDTD: C = 4", "Im\{T\} FDTD: C = 4", "Re\{T\} FDTD: C = 1.04", "Im\{T\} FDTD: C = 1.04", "Re\{T\} Analytical", "Im\{T\} Analytical", 'Location', 'northeast', 'Interpreter', 'latex')
legend("Re\{$T$\}", "Im\{$T$\}", "Re\{$T_{analytical}$\}", "Im\{$T_{analytical}$\}", 'Location', 'northeast', 'Interpreter', 'latex')
xlabel("Frequnecy [GHz]", 'Interpreter', 'latex')
ylabel("$E_{m2}$/$E_{m1}$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2 + 60 \frac{1 + (j\omega \tau)^{0.2}}{1 + 9(j\omega \tau)^{0.3} + 2(j\omega \tau)^{0.5} + 10(j\omega \tau)^{0.9}}$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2 + 60 \frac{1}{1 + 10(j\omega \tau)^{0.9}}$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2 + 60 \frac{1}{1 + 9(j\omega \tau)^{0.3} + 2(j\omega \tau)^{0.5}}$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2 + 60 \frac{1 + (j\omega \tau)^{0.2}}{1 + 9(j\omega \tau)^{0.3} + 2(j\omega \tau)^{0.5}}$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2 + 60 \frac{1}{1 + (j\omega \tau)^{0.9}}$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2 + 10 \frac{1}{1 + (j\omega \tau)^{0.9}}$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2 + 10 \frac{1}{1 + 10(j\omega \tau)^{0.9}}$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 62$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2 + 60 \frac{1}{1 + 9(j\omega \tau)^{0.3} + 2(j\omega \tau)^{0.5} + 10(j\omega \tau)^{0.9}}$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2 + 60 \frac{1}{1 + 5(j\omega \tau)^{0.9}}$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2 + 60 \frac{1}{1 + (j\omega \tau)^{0.9} + (j\omega \tau)^{0.9} + (j\omega \tau)^{0.9} + (j\omega \tau)^{0.9} + (j\omega \tau)^{0.9}}$", 'Interpreter', 'latex')

grid minor


%% 

omega = 2*pi*fff;

% if length(beta) == 5
%     % eps_r = eps_inf + eps_rs*(1)./(1 + b(1)*(1i*omega).^beta(1) + b(2)*(1i*omega).^beta(2) + b(3)*(1i*omega).^beta(3));
%     eps_r = eps_inf + eps_rs*(1 + a(1)*(1i*omega).^alpha(1))./(1 + b(1)*(1i*omega).^beta(1) + b(2)*(1i*omega).^beta(2) + b(3)*(1i*omega).^beta(3) + b(4)*(1i*omega).^beta(4) + b(5)*(1i*omega).^beta(5));
% end
if length(beta) == 3
    % eps_r = eps_inf + eps_rs*(1)./(1 + b(1)*(1i*omega).^beta(1) + b(2)*(1i*omega).^beta(2) + b(3)*(1i*omega).^beta(3));
    eps_r = eps_inf + eps_rs*(1 + a(1)*(1i*omega).^alpha(1))./(1 + b(1)*(1i*omega).^beta(1) + b(2)*(1i*omega).^beta(2) + b(3)*(1i*omega).^beta(3));
end
% if length(beta) == 1
%     eps_r = eps_inf + eps_rs*(1 + a(1)*(1i*omega).^alpha(1))./(1 + b(1)*(1i*omega).^beta(1));
% end
% if length(beta) == 2
%     eps_r = eps_inf + eps_rs*(1 + a(1)*(1i*omega).^alpha(1))./(1 + b(1)*(1i*omega).^beta(1) + b(2)*(1i*omega).^beta(2));
% end

% fabry-perot comparison
n1 = 1; n2 = sqrt(eps_r);
k = omega.*sqrt(eps_r)./c0;
k0 = omega./c0;

d1 = 0;
d2 = 0;
dm = i_media_stop - i_media_start;
ds = i_media_start - isource;

r_air = (n1 - n2)./(n1 + n2);
r_media = -r_air;
t_squared = (4*n1*n2./(n1 + n2).^2);
power_sum = 1./(1 - r_media.^2 .* exp(-1i*2*k*dm*dz));
numer = exp(-1i*k0*(ds + d2)*dz) .* exp(-1i*k*dm*dz) .* t_squared .* power_sum;
denom = exp(-1i*k0*(ds - d1)*dz) + r_air .* exp(-1i*k0*(ds + d1)*dz) .* (1 - t_squared.*exp(-1i*k*2*dm*dz).*power_sum);
T_fp = numer./denom;

pos1 = i_media_stop+d2+1;
Ex_pos1 = Ex(:, pos1);
Ex_pos1_freq = fft(Ex_pos1);

pos2 = i_media_start-d1+1;
Ex_pos2 = Ex(:, pos2);
Ex_pos2_freq = fft(Ex_pos2);

T = Ex_pos1_freq./Ex_pos2_freq;
fff = (1/dt)*(0:nmax-1)/nmax;

relative_error_real = abs((real(T)' - real(T_fp)));
relative_error_imag = abs((imag(T)' - imag(T_fp)));
relative_error = abs(reshape(T,size(T_fp)) - T_fp)./abs(T_fp);

figure(7)
% semilogx(fff*1e-9, relative_error_real, 'b')
% hold on
% semilogx(fff*1e-9, relative_error_imag, 'r')
loglog(fff*1e-9, relative_error, 'k.-')
% hold on
% loglog(fff_normal*1e-9, relative_error_normal, 'k-.')
xlim([0.1 10])
ylim([0.0001 1])
ylabel("Relative Error", 'Interpreter', 'latex')
xlabel("Frequnecy [GHz]", 'Interpreter', 'latex')
% legend("PML = 150$\Delta z$", "PML = 50$\Delta z$", 'Location', 'northeast', 'Interpreter', 'latex')

% legend("Abs\{Error\}", "Re\{Error\}", "Im\{Error\}", 'Location', 'northeast', 'Interpreter', 'latex')
% legend("Re\{Error\}: C = 2", "Im\{Error\}: C = 2", "Re\{Error\}: C = 4", "Im\{Error\}: C = 4", "Re\{Error\}: C = 1.04", "Im\{Error\}: C = 1.04", 'Location', 'northeast', 'Interpreter', 'latex')

% title("$\epsilon_r(\omega) = 2 + 60 \frac{1 + (j\omega \tau)^{0.2}}{1 + 9(j\omega \tau)^{0.3} + 2(j\omega \tau)^{0.5} + 10(j\omega \tau)^{0.9}}$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2 + 60 \frac{1}{1 + 10(j\omega \tau)^{0.9}}$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2 + 60 \frac{1}{1 + 9(j\omega \tau)^{0.3} + 2(j\omega \tau)^{0.5}}$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2 + 60 \frac{1 + (j\omega \tau)^{0.2}}{1 + 9(j\omega \tau)^{0.3} + 2(j\omega \tau)^{0.5}}$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2 + 60 \frac{1}{1 + (j\omega \tau)^{0.9}}$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2 + 10 \frac{1}{1 + (j\omega \tau)^{0.9}}$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2 + 10 \frac{1}{1 + 10(j\omega \tau)^{0.9}}$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 62$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2 + 60 \frac{1}{1 + 9(j\omega \tau)^{0.3} + 2(j\omega \tau)^{0.5} + 10(j\omega \tau)^{0.9}}$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2 + 60 \frac{1}{1 + 5(j\omega \tau)^{0.9}}$", 'Interpreter', 'latex')
% title("$\epsilon_r(\omega) = 2 + 60 \frac{1}{1 + (j\omega \tau)^{0.9} + (j\omega \tau)^{0.9} + (j\omega \tau)^{0.9} + (j\omega \tau)^{0.9} + (j\omega \tau)^{0.9}}$", 'Interpreter', 'latex')

grid minor

% figure(3)
% semilogx(fff*1e-9, real(T), 'kx')
% hold on
% semilogx(fff*1e-9, imag(T), 'k+')
% semilogx(f*1e-9, real(T_fp), 'b')
% semilogx(f*1e-9, imag(T_fp), 'r--')
% xlim([0.1 10])
% legend("Re\{T\} FDTD", "Im\{T\} FDTD", "Re\{T\} Analytical", "Im\{T\} Analytical", 'Location', 'northeast', 'Interpreter', 'latex')
% xlabel("Frequnecy [GHz]")
% grid minor


%%
function wave_source = source_function(t, dt)
    fc = 6e9;
    fbw = 5/6;
    tau_width = fbw;
    t0 = tau_width;
    w0 = 2*pi*fc;

    [yi,yq] = gauspuls((t-150)*dt, fc, fbw, -3);
    wave_source = yq;

%     wave_source = exp(-(t-t0)^2 / tau_width^2) * sin(w0*t*dt);
%     wave_source = exp(-(dt*t-t0*dt)^2 / tau_width^2) * sin(w0*t*dt);
end


function [Amp, alfa, freq, theta] = polynomial_method(x,p,Ts)
    N = length(x);
    T = toeplitz(x(p:N-1), x(p:-1:1));
    a = tls(T, -x(p+1:N));
    
    indeterminate_form = sum(isnan(a) | isinf(a));
    if indeterminate_form
        return;
    end
    
    c = transpose([1; a]);
    r = roots(c);
    alfa = log(abs(r))/Ts;
    freq = atan2(imag(r), real(r))/(2*pi*Ts);
    
    alfa(isinf(alfa)) = realmax*sign(alfa(isinf(alfa)));
    len_vandermonde = N;
    Z = zeros(len_vandermonde, p);
    
    for i=1:length(r)
        Z(:,i) = transpose(r(i).^(0:len_vandermonde-1));
    end
    
    rZ = real(Z);
    iZ = imag(Z);
    rZ(isinf(rZ)) = realmax*sign(rZ(isinf(rZ)));
    iZ(isinf(iZ)) = realmax*sign(iZ(isinf(iZ)));
    
    Z = rZ + 1i*iZ;
    
    indeterminate_formZ = sum(isnan(Z) | isinf(Z));
    if indeterminate_formZ
        return;
    else
        h = tls(Z, x(1:len_vandermonde));
    end
    Amp = abs(h);
    theta = atan2(imag(h),real(h));
end

function x = tls(A,b)
    [~,n] = size(A);
    C = [A b'];
    [~,~,V] = svd(C);
    VXY = V(1:n,1+n:end);
    VYY = V(1+n:end, 1+n:end);
    if VYY == 0
        x = zeros(n,1);
        disp('Error: not tls solution');
        return;
    end
    x = -VXY/VYY;
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