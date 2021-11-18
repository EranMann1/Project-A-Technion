function [WiresElectricField] = HenklesSum(Lambda,k,eta,theta_in,Ds,Hs,Is,M,N)
    % calculates the electric field at the m-th wire created by the other
    % wires
    % 
    
    phase = -k*sin(theta_in)*Lambda;
    WiresElectricField=0;
    nEqPhaseVec = exp(1j*(-N:N)*phase);
    EqPhaseVec = [ exp(1j*(-N:-1)*phase),  exp(1j*(1:N)*phase)];
    for m=1:length(Is)
       if m==M
           distanceVec = abs((Lambda*[(-N:-1),(1:N)]')); 
           HenVec = besselh(0,2,k*distanceVec);
           sum = EqPhaseVec*HenVec;
       else 
           distanceVec = sqrt( (Ds(M)-(Ds(m)+Lambda*(-N:N)')).^2 + (Hs(M)-Hs(m))^2 ); 
           HenVec = besselh(0,2,k*distanceVec);
           sum = nEqPhaseVec*HenVec;
       end
       WiresElectricField = WiresElectricField - k*eta/4*Is(m)*sum;
    end
    
end