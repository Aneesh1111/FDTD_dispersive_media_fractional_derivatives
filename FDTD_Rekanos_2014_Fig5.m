close all; clear variables; clc;

% constants
f = (0.1:0.1:10)*1e9;
eps0  = 8.8541878128e-12;
mu0 = 1.256637062e-6;
c0 = 1/sqrt(eps0*mu0);
imp0 = sqrt(mu0/eps0);
omega = 2*pi*f;
tau = 153e-12;
eps_s = 2;%50;
eps_inf = 2;
alpha = 0.9;
d_eps = eps_s - eps_inf;

% analytical solution
eps_r = eps_inf + d_eps./(1 + (1i*omega*tau).^alpha);
R_analytical = abs((1 - sqrt(eps_r))./(1 + sqrt(eps_r)));

% plot permittivity
% figure(4)
% semilogx(f, real(eps_r))
% figure(5)
% semilogx(f, imag(eps_r))

% plotting
% figure(3)
% semilogx(f*1e-9, R_analytical)
% ylim([0.5 0.8])

dz = 1.1e-3;
d = 30*dz;



% Prony's method
N = 2000;
wk = zeros(1,N);
wk(1) = 1;
for k = 1:N-1
    wk(k+1) = wk(k)*(1 - (1+alpha)/k);
end
k = 0:1:N-1;
M = 8;
Ts = 1;
[Amp, alfa, freq, theta] = polynomial_method(wk,M,Ts);

% GPOF method
L = M;
[b_gpof, a_gpof] = gpof(wk,k,M,L+1);

% FDTD parameters
dt = 1.768e-12;
dz = 1.1e-3;
% lm = 350e-9; dz = lm/20; dt = dz/c0;
imax = 500;  % max discretization in space
nmax = 10000; % max discretization in time
isource = 100;
i_media_start = 200;
i_media_stop = 230;

% constant
q = (tau/dt)^alpha;
A = 1/(eps0*eps_inf + eps0*d_eps/(1 + q));
qq = q/(1 + q);
qq_eps = eps0*d_eps/(1 + q);
tz_mu = dt/(mu0*dz);

eps_dielectric = eps0*ones(1, i_media_stop-i_media_start+2)*10;

% update parametes
vj = zeros(M, imax);
Ex = zeros(nmax, imax);
P = zeros(nmax, imax);
Hy = zeros(nmax, imax-1);
vj_prev = zeros(M, imax);
Ex_prev = zeros(nmax, imax);
P_prev = zeros(1, imax);
Hy_prev = zeros(1, imax);

mur_boundary_const_first = (c0*dt - dz)/(c0*dt + dz);
mur_boundary_const_first_dielectric = (c0*dt/sqrt(eps_inf) - dz)/(c0*dt/sqrt(eps_inf) + dz);

% FDTD air-media-air
% for n = 2:nmax-1    
%     % update vj
%     for i = i_media_start:i_media_stop
% %         for j = 1:M
% %             vj(j, i) = -abs(a_gpof(j))*exp(-abs(b_gpof(j)))*P(n-1,i) + exp(-abs(b_gpof(j)))*vj_prev(j, i);
% % %             vj(j, i) = a_gpof(j)*exp(b_gpof(j))*P(n-1,i) + exp(b_gpof(j))*vj_prev(j, i); 
% %             vj_prev(j, i) = vj(j, i);
% %         end
%         vj(:, i) = a_gpof.*exp(b_gpof)*P(n-1,i) + exp(b_gpof).*vj_prev(:, i);
%         vj_prev(:, i) = vj(:, i);
%     end
%     vj_sum = sum(vj);
%     
%     % update electric field in media
%     for i = i_media_start:i_media_stop+1
%         Ex(n,i) = A*(eps0*eps_inf*Ex(n-1,i) + P(n-1,i) + qq*vj_sum(i) + dt*(Hy(n-1,i-1) - Hy(n-1,i))/dz);  % Rekanos
% %         Ex(n,i) = Ex(n-1,i) + dt/(dz*eps0) * (Hy(n-1,i) - Hy(n-1,i-1));  % free space
% %         Ex(n,i) = Ex(n-1,i) + dt/(dz*eps_dielectric(i-i_media_start+1)) * (Hy(n-1,i) - Hy(n-1,i-1));  % simple dielectric
%     end
%     % update electric field in air
%     for i = 2:i_media_start-1
%         Ex(n,i) = Ex(n-1,i) + dt/(dz*eps0) * (Hy(n-1,i) - Hy(n-1,i-1));
%     end
%     % update electric field in air
%     for i = i_media_stop+2:imax-1
%         Ex(n,i) = Ex(n-1,i) + dt/(dz*eps0) * (Hy(n-1,i) - Hy(n-1,i-1));
%     end
%     % update electric boundaries
%     Ex(n,1) = Ex(n-1,2) + mur_boundary_const_first*(Ex(n,2) - Ex(n-1,1));
%     Ex(n,imax) = Ex(n-1,imax-1) + mur_boundary_const_first*(Ex(n,imax-1) - Ex(n-1,imax));
%     % electric field source
%     Ex(n,isource) = Ex(n,isource) + source_function(n, dt);
%         
%     % update polatizarion vector
%     for i = i_media_start:i_media_stop
%         P(n,i) = qq_eps*Ex(n,i) - qq*vj_sum(i);
%     end
%     
%     % update magnetic field in media
%     for i = i_media_start:i_media_stop
%         Hy(n,i) = Hy(n-1,i) - tz_mu*(Ex(n,i+1) - Ex(n,i));  % Rekanos
% %         Hy(n,i) = Hy(n-1,i) + dt/(dz*mu0)*(Ex(n,i+1) - Ex(n,i));  % free space or dielectric
%     end
%     % update magnetic field in air
%     for i = 1:i_media_start-1
%         Hy(n,i) = Hy(n-1,i) + dt/(dz*mu0)*(Ex(n,i+1) - Ex(n,i));
%     end
%     % update magnetic field in air
%     for i = i_media_stop+1:imax-1
%         Hy(n,i) = Hy(n-1,i) + dt/(dz*mu0)*(Ex(n,i+1) - Ex(n,i));
%     end
%     % magnetic field source
%     Hy(n,isource-1) = Hy(n,isource-1) - source_function(n, dt)/imp0;
%     
% %     if mod(n,20) == 0
% %         figure(1)
% %         plot(Ex(n,:));
% %         hold on
% %         xline(i_media_start);
% %         xline(i_media_stop);
% %         hold off
% %         ylim([-2, 2])
% %         xlim([0 500])
% %         getframe;
% %     end
%     disp(n)
% end

% FDTD scheme: media only
for n = 2:nmax-1    
    % update vj
    for i = 1:imax
        for j = 1:M
%             vj(j, i) = -Amp(j)*exp(alfa(j))*P(n-1,i) + exp(alfa(j))*vj_prev(j, i);
            vj(j, i) = -abs(a_gpof(j))*exp(-abs(b_gpof(j)))*P(n-1,i) + exp(-abs(b_gpof(j)))*vj_prev(j, i); 
            vj_prev(j, i) = vj(j, i);
        end
    end
    vj_sum = sum(vj);
    
    % update electric field in media
    for i = 2:imax-1
%         Ex(n,i) = A*(eps0*eps_inf*Ex(n-1,i) + P(n-1,i) + qq*vj_sum(i) + dt*(Hy(n-1,i-1) - Hy(n-1,i))/dz);  % Rekanos
%         Ex(n,i) = Ex(n-1,i) + dt/(dz*eps0) * (Hy(n-1,i) - Hy(n-1,i-1));  % free space
        Ex(n,i) = Ex(n-1,i) + dt/(dz*eps0*eps_inf) * (Hy(n-1,i) - Hy(n-1,i-1));  % constant simple dielectric
%         Ex(n,i) = Ex(n-1,i) + dt/(dz*eps_dielectric(i-i_media_start+1)) * (Hy(n-1,i) - Hy(n-1,i-1));  % simple air-dielectric-air
    end
    % update electric boundaries air
%     Ex(n,1) = Ex(n-1,2) + mur_boundary_const_first*(Ex(n,2) - Ex(n-1,1));
%     Ex(n,imax) = Ex(n-1,imax-1) + mur_boundary_const_first*(Ex(n,imax-1) - Ex(n-1,imax));
    % update electric boundaries simple constant dielectric
    Ex(n,1) = Ex(n-1,2) + mur_boundary_const_first_dielectric*(Ex(n,2) - Ex(n-1,1));
    Ex(n,imax) = Ex(n-1,imax-1) + mur_boundary_const_first_dielectric*(Ex(n,imax-1) - Ex(n-1,imax));
    % electric field source
    Ex(n,isource) = Ex(n,isource) + source_function(n, dt);
        
    % update polatizarion vector
    for i = 1:imax
        P(n,i) = qq_eps*Ex(n,i) - qq*vj_sum(i);
    end
    
    % update magnetic field in media
    for i = 1:imax-1
        Hy(n,i) = Hy(n-1,i) - tz_mu*(Ex(n,i+1) - Ex(n,i));  % Rekanos
%         Hy(n,i) = Hy(n-1,i) + dt/(dz*mu0)*(Ex(n,i+1) - Ex(n,i));  % free space or dielectric
    end
    % magnetic field source
%     Hy(n,isource-1) = Hy(n,isource-1) - source_function(n, dt)/imp0;
    
%     if mod(n,1) == 0
%         figure(1)
%         plot(Ex(n,:));
%         hold off
%         ylim([-2, 2])
%         xlim([0 500])
%         getframe;
%     end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dielectric post-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reflec_coeff = abs((1 - sqrt(eps_dielectric(1)))/(1 + sqrt(eps_dielectric(1))));
% 
% trans_coeff = exp(-1i*omega*d*sqrt(eps_dielectric(1))/c0);
% 
% source_t = 1:1:nmax;
% source_input = gauspuls((source_t-100)*dt, 6e9, 5/6)';
% % source_input = sin(2*pi*6e9*source_t*dt)';
% source_input_fft = fft(source_input);
% 
% pos1 = 290;
% Ex_pos1 = Ex(:, pos1);
% Ex_pos1_freq = fft(Ex_pos1);
% 
% T = Ex_pos1_freq./source_input_fft;
% fff = (1/dt)*(0:nmax-1)/nmax;
% figure(6)
% % semilogx(fff*1e-9, real(Ex_pos1_freq))
% semilogx(fff*1e-9, real(T))
% hold on
% semilogx(fff*1e-9, imag(T), '--')
% % semilogx(f*1e-9, real(T_analytical))
% % semilogx(f*1e-9, imag(T_analytical))
% xlim([0.1 10])
% ylim([-1 1])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rekanos post-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%
T_analytical = exp(-1i*omega.*d.*sqrt(eps_r)/c0) ;

pos1 = i_media_stop;
Ex_pos1 = Ex(:, pos1);
Ex_pos1_freq = fft(Ex_pos1);

pos2 = i_media_start;
Ex_pos2 = Ex(:, pos2);
Ex_pos2_freq = fft(Ex_pos2);

T = Ex_pos1_freq./Ex_pos2_freq;
fff = (1/dt)*(0:nmax-1)/nmax;
figure(6)
% semilogx(fff*1e-9, real(Ex_pos1_freq))
semilogx(fff*1e-9, real(T), 'b-.')
hold on
semilogx(fff*1e-9, imag(T), 'r-.')
semilogx(f*1e-9, real(T_analytical), 'b')
semilogx(f*1e-9, imag(T_analytical), 'r')
xlim([0.1 10])
% ylim([-1 1])
legend("Re\{T\} FDTD", "Im\{T\} FDTD", "Re\{T\} Analytical", "Im\{T\} Analytical", 'Location', 'northeast', 'Interpreter', 'latex')
xlabel("Frequnecy [GHz]")
grid minor


function wave_source = source_function(t, dt)
    fc = 6e9;
    fbw = 5/6;
%     tau_width = 30;
%     t0 = tau_width*3;
    w0 = 2*pi*fc;
%     wave_source = sin(w0*t*dt);
    wave_source = gauspuls((t-100)*dt, fc, fbw);
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