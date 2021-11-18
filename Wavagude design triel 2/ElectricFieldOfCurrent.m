function Field = ElectricFieldOfCurrent(I_0,y_0,z_0,ys,zs,k,eta)
% ys,zs are first and second dimention vectors respectivly
Field = -k*eta/4*I_0*besselh(0,2,k*sqrt((ys-y_0).^2+(zs-z_0).^2));
end