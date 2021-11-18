function [ElecricFieldOnWires] = FindElectricFiield(Lambda,k,eta,theta_in,theta_out,E_in,Ds,Hs,Is,reff,N,type)
    % the fufunction calculates the electric fild at each wire
    % vars: pirode wave vector, wave impedance, angle of incident, electric
    % field amplitude, the postions of the wires, th currntes of the wires,
    % the effctive radius of the wires and the order of calculation
    ElecricFieldOnWires = zeros(length(Is),1);
    ElectricIn=@(d,h) E_in*exp(-1j*k*(sin(theta_in)*d +cos(theta_in)*h));
    if strcmp(type,'Henkels')
        SelfFiels=@(selfI)-k*eta/4*selfI*besselh(0,2,k*reff);
        for m=1:length(Is)
            ElecricFieldOnWires(m) = ElectricIn(Ds(m),Hs(m)) + SelfFiels(Is(m)) + HenklesSum(Lambda,k,eta,theta_in,Ds,Hs,Is,m,N);
        end
    end
    if strcmp(type,'Poisson')
         for m=1:length(Is)
            ElecricFieldOnWires(m) = ElectricIn(Ds(m),Hs(m)) + PoissonSum(Lambda,k,eta,theta_in,theta_out,Ds,Hs,Is,reff,m,N);
        end
    end
end