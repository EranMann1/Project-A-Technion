clear;
close all;
%% large H pull
D=0;
innerPath=pwd;
cd('..');
path=pwd;
cd('..');
cd([pwd '\designing reflective serface\after ariel meating 3']);
% reflective_surface = load('FinalReflector.mat');
reflective_surface = load('realy small thickness surface.mat');
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
numH = 400;
numKy = 1000;
N_sumation=100;
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


%% small H
H_const = cos(theta_in)*cos(theta_out)/2/(cos(theta_in)+cos(theta_out));

numH = 400;
numKy = 5000;
N_sumation=100;
Ky = k*linspace(-1,1,numKy);
Hs = linspace(0.01,4*H_const,numH);
Dets = zeros(numH,numKy);
for index=1:numH
    index
    waveguide = createWaveguide(D,Hs(index),ds,hs,zs);
    Dets2(index,:)=detsOfKy(waveguide,k,Ky,N_sumation,eta,r_eff,Lambda);
end

figure;
imagesc(Hs/H_const,Ky/k,log(abs(Dets2)')/log(10));
colorbar;
xlabel('Hs/H_const');
ylabel('ky*lambda/2pi');
zlabel('log(det(A))');

%% D changes H=lambda
H_const = cos(theta_in)*cos(theta_out)/2/(cos(theta_in)+cos(theta_out));

H=Lambda;

numD = 100;
numKy = 400;
N_sumation=100;
Ky = k*linspace(-1,1,numKy);
Ds = linspace(0,Lambda,numD);
Dets = zeros(numD,numKy);
for index=1:numD
    index
    waveguide = createWaveguide(Ds(index),H,ds,hs,zs);
    Dets5(index,:)=detsOfKy(waveguide,k,Ky,N_sumation,eta,r_eff,Lambda);
end

figure;
imagesc(Ds/Lambda,Ky/k,log(abs(Dets5)')/log(10));
colorbar;
xlabel('D/Lambda');
ylabel('ky*lambda/2pi');
zlabel('log(det(A))');