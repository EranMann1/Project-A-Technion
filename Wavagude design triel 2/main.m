
% takes the reflective surface and makes it a waveguide, calculates the
% electric field in a 2D grid in space for 2 diffrent waveguides - first 
% constructive interfirance and first distructive interfirance. 
% plots the resulting fields

clear all;
Title = '2D constructive and distructive interfirance very large space with wide gap big n';
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
H_const = 20*lambda/(2*(1/cos(theta_in)+1/cos(theta_out)));
H_dist = 20*lambda/(4*(1/cos(theta_in)+1/cos(theta_out)));
D=0;
constructive_waveguide = createWaveguide(D,H_const,ds,hs,zs);
distructive_waveguide = createWaveguide(D,H_dist,ds,hs,zs);

%% parameters of grid:
z_const = linspace(-2*max(constructive_waveguide.Hs),2*max(constructive_waveguide.Hs),100);
z_dist = linspace(-2*max(constructive_waveguide.Hs),2*max(constructive_waveguide.Hs),100);
y = linspace(-5*Lambda,5*Lambda,60); %on source
% y = linspace(20*Lambda,30*Lambda,100); %far from sorce
k_y=zeros(402,1);
k_y(:) = (-2.005*k:0.01*k:2.005*k);
N=200;

%% calculations
E_const = zeros(length(y),length(z_const));
E_dist = zeros(length(y),length(z_dist));
I_const = currentsOfPlaneWave(constructive_waveguide,k,k_y,N,eta,r_eff,Lambda);
I_dist =  currentsOfPlaneWave(distructive_waveguide,k,k_y,N,eta,r_eff,Lambda);
for indexy=1:length(y)
    for indexz=1:length(z_const)
        E_const(indexy,indexz) = find_generated_field(y(indexy),z_const(indexz),constructive_waveguide,k,Lambda,eta,I_0,k_y,I_const,N,y_0,z_0);
        E_dist(indexy,indexz) = find_generated_field(y(indexy),z_dist(indexz),distructive_waveguide,k,Lambda,eta,I_0,k_y,I_dist,N,y_0,z_0);
    end
    indexy % to see where we are in the simulation
end
E_reg_const = ElectricFieldOfCurrent(I_0,y_0,z_0,y',z_const,k,eta);
E_reg_dist = ElectricFieldOfCurrent(I_0,y_0,z_0,y',z_dist,k,eta);
E_const_tot = E_const + E_reg_const;
E_dist_tot = E_dist + E_reg_dist;

%% ploting and saving
z_const_norm = z_const/H_const;
z_dist_norm = z_dist/H_dist;
y_norm = y/Lambda;

a1=figure(1);
surf(z_const_norm,y_norm,real(E_const));
set(get(gca,'Children'),'EdgeColor','None');
ylabel('y [Lambda]');
xlabel('z [H]');
zlabel('generated electric field');
title('generated electric field in constructive interfirance waveguide');
view([-37.5,30]);

a2=figure(2);
surf(z_const_norm,y_norm,real(E_const_tot));
set(get(gca,'Children'),'EdgeColor','None'); 
ylabel('y [Lambda]');
xlabel('z [H]');
zlabel('total electric field');
title('toral electric field in constructive interfirance waveguide');
view([-37.5,30]);

a3=figure(3);
surf(z_dist_norm,y_norm,real(E_dist));
set(get(gca,'Children'),'EdgeColor','None'); 
ylabel('y [Lambda]');
xlabel('z [H]');
zlabel('generated electric field');
title('generated electric field in distructive interfirance waveguide');
view([-37.5,30]);

a4=figure(4);
surf(z_dist_norm,y_norm,real(E_dist_tot));
set(get(gca,'Children'),'EdgeColor','None'); 
ylabel('y [Lambda]');
xlabel('z [H]');
zlabel('total electric field');
title('total electric field in distructive interfirance waveguide');
view([-37.5,30]);

NewSave