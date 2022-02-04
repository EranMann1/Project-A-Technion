function [ElecricFieldOnWires] = FindElectricFiield(Lambda,k,eta,theta_in,E_in,Ds,Hs,Is,reff,N)
    % the fufunction calculates the electric fild at each wire
    % vars: pirode wave vector, wave impedance, angle of incident, electric
    % field amplitude, the postions of the wires, th currntes of the wires,
    % the effctive radius of the wires and the order of calculation
    ElectricIn=@(d,h) E_in*exp(-1j*k*(sin(theta_in)*d +cos(theta_in)*h));
    SelfFiels=@(selfI)-k*eta/4*selfI*besselh(0,2,k*reff);
    ElecricFieldOnWires = zeros(length(Is),1);
    for m=1:length(Is)
        ElecricFieldOnWires(m) = ElectricIn(Ds(m),Hs(m)) + SelfFiels(Is(m)) + HenklesSum(Lambda,k,eta,theta_in,Ds,Hs,Is,m,N);
    end
    
end