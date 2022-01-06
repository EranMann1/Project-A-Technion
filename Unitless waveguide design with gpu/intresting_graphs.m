clear all;
%% louding relevent informatoin
cd('C:\Users\user\OneDrive - Technion\Documents\MATLAB\technion\semester 5\project a\designing reflective serface\after ariel meating 3');
reflective_surface = load('FinalReflector.mat');
cd('C:\Users\user\OneDrive - Technion\Documents\MATLAB\technion\semester 5\project a\Wavagude design triel 2');
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

%% calculating the relevent waveguides
H_const = lambda/(2*(1/cos(theta_in)+1/cos(theta_out)));
D=0;
constructive_waveguide = createWaveguide(D,H_const,ds,hs,zs);
y_0=0;
z_0=0;

%% parameters of grid:
y_runs = linspace(-50*Lambda,50*Lambda,1000);
z_middle_waveguide = 0;
z_out_of_waveguide = 2*max(constructive_waveguide.Hs);
z_runs = linspace(-2*max(constructive_waveguide.Hs),2*max(constructive_waveguide.Hs),100);
y_middle = 0;
y_far = 50*Lambda;
k_y=zeros(402,1);
k_y(:) = (-2.005*k:0.01*k:2.005*k);
dk_y = k_y;
dk_y(:) = 0.01*k*ones(length(k_y),1);
N=100;

%% calculations
E_middle_waveguide = zeros(length(y_runs),1);
E_out_of_waveguide = zeros(length(y_runs),1);
E_cross_waveguide = zeros(length(z_runs),1);
E_side_waveguide = zeros(length(y_runs),1);
I_const = currentsOfPlaneWave(constructive_waveguide,k,k_y,N,eta,r_eff,Lambda);

for index=1:length(y_runs)
    E_middle_waveguide(index) = find_generated_field(y_runs(index),z_middle_waveguide,constructive_waveguide,k,Lambda,eta,I_0,k_y,dk_y,I_const,N,y_0,z_0);
    E_out_of_waveguide(index) = find_generated_field(y_runs(index),z_out_of_waveguide,constructive_waveguide,k,Lambda,eta,I_0,k_y,dk_y,I_const,N,y_0,z_0); 
end

for index=1:length(z_runs)
    E_cross_waveguide(index) = find_generated_field(y_middle,z_runs(index),constructive_waveguide,k,Lambda,eta,I_0,k_y,dk_y,I_const,N,y_0,z_0);
    E_side_waveguide(index) = find_generated_field(y_far,z_runs(index),constructive_waveguide,k,Lambda,eta,I_0,k_y,dk_y,I_const,N,y_0,z_0); 
end

 E_middle_waveguide =  E_middle_waveguide + ElectricFieldOfCurrent(I_0,y_0,z_0,y_runs,z_middle_waveguide,k,eta);
 E_out_of_waveguide =  E_out_of_waveguide + ElectricFieldOfCurrent(I_0,y_0,z_0,y_runs,z_out_of_waveguide,k,eta);
 E_cross_waveguide =  E_cross_waveguide + ElectricFieldOfCurrent(I_0,y_0,z_0,y_middle,z_runs,k,eta);
 E_side_waveguide =  E_side_waveguide + ElectricFieldOfCurrent(I_0,y_0,z_0,y_far,z_runs,k,eta);
 
 %% plloting
 y_divided = y_runs/Lambda;
 z_devided = z_runs/min(abs(constructive_waveguide.Hs));
 
 figure(1);
 plot(y_divided,abs(E_middle_waveguide));
 hold on;
 plot(y_divided,abs(E_out_of_waveguide));
 title('electric field in intresting crossections as a function of y');
 xlabel('y [Lambda]');
 ylabel('Electric Field');
 legend('electric field in middle of the waveguide', 'electric field out of the waveguide');
 hold off;
 
 figure(2);
 plot(z_devided,abs(E_cross_waveguide));
 hold on;
 plot(z_devided,abs(E_side_waveguide));
 title('electric field in intresting crossections as a function of z');
 xlabel('z [H]');
 ylabel('Electric Field');
 legend('electric field cross waveguide on sourse', 'electric field cross waveguide far from sourse');
 hold off;
 
 