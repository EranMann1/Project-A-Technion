function [Hy, Hz] = find_generated_Hyz_field(y,z,waveguide,k,Lambda,eta,I_0,k_y,I_norm,N,y_0,z_0)
    % calc the magnetic field created by the currents runnig throw the
    % particls. each current corespond to diffrent k_y for a spesific point
    % in space (y,z)
    
    % k_y and dk_y are vectors of the first dimention
    % I_norm is a 2D length(k_y) by number of particles
    % waveguide contains Hs,Ds are in the second dimention as well
    % k_y cannot contain +-k (otherwise E_out will be NaN)
    % the amaunt of memory in use is length(k)*lentgh(I_norm)*N
    
    Hs = waveguide.Hs;
    Ds = waveguide.Ds;
    Ns = zeros(1,1,2*N+1); % the third dimention
    Ns(:) = (-N:N);
    Alphas = 2*pi*Ns/Lambda + k_y;
    Betas = conj(sqrt(k^2-Alphas.^2));
    
    k_z = conj(sqrt(4*pi^2-k_y.^2));
    %sumTensor = exp(-1j*(k_z*abs(z-z_0)+k_y*y_0+Alphas.*(y - Ds)+Betas.*abs(z-Hs)))./k_z.*I_norm./Betas;
    sumTensor = exp(-1j*(k_z*abs(z-z_0)+k_y*(y-y_0)+Alphas.*(y - Ds)+Betas.*abs(z-Hs)))./k_z.*I_norm./Betas;
    %E_out = k^2*eta^2*I_0/4/Lambda*sum(sumTensor,'all');
    Hy = k/(4*Lambda)*trapz(k_y,sum(sumTensor.*(k_z.*sign(z-z_0)+Betas.*sign(z-Hs)),[2,3]),1);
    Hz = -k/(4*Lambda)*trapz(k_y,sum(sumTensor.*(k_y+Alphas),[2,3]),1);
end