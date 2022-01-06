clear;
close all;
%% D=0 H=H_dist
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
H = 1/2*cos(theta_in)*cos(theta_out)/2/(cos(theta_in)+cos(theta_out));% H_dist
waveguide = createWaveguide(D,H,ds,hs,zs);
N_sumation=100;

reKy=k*linspace(-2,2,101);
imagKy=k*linspace(-1.5,1.5,101);
for index1=1:length(reKy)
    disp(index1)
    for index2=1:length(imagKy)
        Dets(index1,index2)=detsOfKy(waveguide,k,reKy(index1)+1j*imagKy(index2),N_sumation,eta,r_eff,Lambda);
    end
end

figure;
imagesc(reKy/k,imagKy/k,log(abs(Dets'))/log(10));
colorbar;
xlabel('Re{K_y}/k');
ylabel('IM{K_y}/k');
zlabel('log(det(A))');