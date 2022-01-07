function [WiresElectricField] = PoissonSum(Lambda,k,eta,theta_in,theta_out,Ds,Hs,Is,reff,M,N)
    % calculates the electric field at the m-th wire created by the other
    % wires
    % 
    n =(-N:N);
    Alphas = k*sin(theta_in)+n*k*(abs(sin(theta_out)-sin(theta_in)));
    Betas = sqrt(k^2-Alphas.^2)';
    WiresElectricField = 0;
    for index =1:length(Is)
        if index == M
            vecY = exp(-1j*(reff)*Alphas);
            vecZ = 1./Betas;
        else
            vecY = exp(-1j*(Ds(M)-Ds(index))*Alphas);
            vecZ = exp(-1j*Betas*abs(Hs(M)-Hs(index)))./Betas;
        end
        WiresElectricField = WiresElectricField-Is(index)*k*eta/(2*Lambda)*vecY*vecZ;
    end
end