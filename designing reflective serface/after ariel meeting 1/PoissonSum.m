function [WiresElectricField] = PoissonSum(Lambda,k,eta,theta_in,theta_out,Ds,Hs,Is,reff,M,N)
    % calculates the electric field at the m-th wire created by the other
    % wires
    % 
    n =(-N:N);
    nM=[(-N:-1),(1:N)];
    Alphas = k*sin(theta_in)+n*k*(abs(sin(theta_out)-sin(theta_in)));
    Betas = sqrt(k^2-Alphas.^2)';
    AlphasM = k*sin(theta_in)+nM*k*(abs(sin(theta_out)-sin(theta_in)));
    BetasM = sqrt(k^2-AlphasM.^2)';
 
    WiresElectricField = 0;
    for index =1:length(Is)
        if index == M
            vecS = 1./(Lambda.*BetasM)-1j./(2*pi*abs(nM'));
            WiresElectricField = WiresElectricField - Is(M)*k*eta/2*sum(vecS);
            WiresElectricField = WiresElectricField + 1j*Is(M)*k*eta/2/pi*log(2*pi*reff/Lambda);
            WiresElectricField = WiresElectricField - Is(M)*eta/2/Lambda/cos(theta_in);
        else
            vecY = exp(-1j*(Ds(M)-Ds(index))*Alphas);
            vecZ = exp(-1j*Betas*abs(Hs(M)-Hs(index)))./Betas;
            WiresElectricField = WiresElectricField-Is(index)*k*eta/(2*Lambda)*vecY*vecZ;
        end
        
    end
end