D=0;
innerPath=pwd;
cd('..');
path=pwd;
cd('..');
cd([pwd '\designing reflective serface\after ariel meating 3']);
reflective_surface = load('FinalReflector.mat');
cd(innerPath);
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
r_eff = reflective_surface.reff/lambda;
clear reflective_surface;
numH = 200;
numKy = 100;
N_sumation=200;
Ky = k*linspace(-1,1,numKy);
Hs = linspace(0.01,10,numH);
Dets = zeros(numH,numKy);
for index=1:numH
    index
    waveguide = createWaveguide(D,Hs(index),ds,hs,zs);
    Dets(index,:)=detsOfKy(waveguide,k,Ky,N_sumation,eta,r_eff,Lambda);
end
figure;
imagesc(Hs,Ky/k,log(abs(Dets)')/log(10));
colorbar;
xlabel('Hs/lambda');
ylabel('ky*lambda/2pi');
