function Min=MinimumDetA(reflective_surface)
D=0;
lambda = reflective_surface.lambda;
Lambda = reflective_surface.Lambda/lambda;
k = 2*pi;
eta = reflective_surface.eta;
zs = reflective_surface.Zs*lambda/eta;
hs = reflective_surface.Hs/lambda;
ds = reflective_surface.Ds/lambda;
I_0 = 1;
z_0 = 0;
y_0 = 0;
theta_in = reflective_surface.theta_in;
theta_out = reflective_surface.theta_out;
r_eff = reflective_surface.r_eff/lambda;
clear reflective_surface;
numH = 400;
numKy = 400;
N_sumation=100;
Ky = k*linspace(0,1,numKy);
Hs = linspace(0.01,10,numH);
Dets = zeros(numH,numKy);
for index=1:numH
    waveguide = createWaveguide(D,Hs(index),ds,hs,zs);
    Dets(index,:)=detsOfKy(waveguide,k,Ky,N_sumation,eta,r_eff,Lambda);
end
Min=min(abs(Dets(:)));




end