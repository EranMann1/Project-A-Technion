function [Hy_reg,Hz_reg] = MagneticFieldOfCurrent(I_0,y0,z0,y,z,k)
    rho = sqrt((y-y0).^2+(z-z0).^2);
    H_phi = I_0*1j/4*besselh(1,2,k*rho);
    Hy_reg = -H_phi.*(z-z0)./rho;
    Hz_reg = H_phi.*(y-y0)./rho;
end