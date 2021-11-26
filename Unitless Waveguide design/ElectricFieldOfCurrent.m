function Field = ElectricFieldOfCurrent(I_0,y_0,z_0,ys,zs,k,eta)
% ys,zs are first and second dimention vectors respectivly
Field = -1/2/pi*besselh(0,2,2*pi*sqrt((ys-y_0).^2+(zs-z_0).^2));
end