D=0;
innerPath=pwd;
cd('..');
path=pwd;
cd('..');
cd([pwd '\designing reflective serface\after ariel meating 3']);
reflective_surface = load('FinalReflector.mat');
cd(innerPath);
Lambda = reflective_surface.Lambda;
k = reflective_surface.k;
eta = reflective_surface.eta;
zs = reflective_surface.Zs;
hs = reflective_surface.Hs;
ds = reflective_surface.Ds;
I_0 = 1;
z_0 = 0;
y_0 = 0;
lambda = reflective_surface.lambda;
theta_in = reflective_surface.theta_in;
theta_out = reflective_surface.theta_out;
r_eff = reflective_surface.reff;
clear reflective_surface;
numH = 200;
numKy = 100;
N_sumation=200;
Ky = k*linspace(-0.9,0.9,numKy);
Hs = lambda*linspace(0.01,10,numH);
Dets = zeros(numH,numKy);
for index=1:numH
    index
    waveguide = createWaveguide(D,Hs(index),ds,hs,zs);
    Dets(index,:)=detsOfKy(waveguide,k,Ky,N_sumation,eta,r_eff,Lambda);
end

surf(Hs/lambda,Ky/k,abs(Dets)');