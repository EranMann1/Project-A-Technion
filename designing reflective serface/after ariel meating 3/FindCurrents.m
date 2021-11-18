function [currents] = FindCurrents(Ds,Hs,Lambda,Ein,Eout,Zin,Zout,theta_in,theta_out,k,phi)
    n = (-1:1);
    Alphas = k*sin(theta_in)+n*k*(abs(sin(theta_out)-sin(theta_in)));
    Betas = sqrt(k^2-Alphas.^2);
    EqautionMatrix = [  exp(1j*Alphas(convert(0))*Ds).*exp(1j*Betas(convert(0))*Hs) ;
                        exp(1j*Alphas(convert(1))*Ds).*exp(1j*Betas(convert(1))*Hs) ;
                        exp(1j*Alphas(convert(-1))*Ds).*exp(1j*Betas(convert(-1))*Hs) ;
                        exp(1j*Alphas(convert(0))*Ds).*exp(-1j*Betas(convert(0))*Hs) ;
                        exp(1j*Alphas(convert(1))*Ds).*exp(-1j*Betas(convert(1))*Hs) ;
                        exp(1j*Alphas(convert(-1))*Ds).*exp(-1j*Betas(convert(-1))*Hs) ];
    EqationVector = [Ein*2*Lambda/Zin ; 0; 0; 0; -Eout*exp(1j*phi)*2*Lambda/Zout ; 0];    
    currents = inv(EqautionMatrix)*EqationVector;
end