%% louding relevent informatoin
path=pwd;
cd('..');
cd([pwd '\designing reflective serface\after ariel meating 3']);
reflective_surface = load('FinalReflector.mat');
cd(path);
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

%% calculating the relevent waveguides
H_const = 1/(2*(1/cos(theta_in)+1/cos(theta_out)));
H_dist = 1/(4*(1/cos(theta_in)+1/cos(theta_out)));
D=0;
constructive_waveguide = createWaveguide(D,H_const,ds,hs,zs);
distructive_waveguide = createWaveguide(D,H_dist,ds,hs,zs);

%% parameters of grid:
z_const = linspace(-max(constructive_waveguide.Hs),max(constructive_waveguide.Hs),60);
y = linspace(-5*Lambda,5*Lambda,60); %on source
% y = linspace(20*Lambda,30*Lambda,100); %far from sorce
k_y=zeros(402,1);
k_y(:) = (-2.005*k:0.01*k:2.005*k);
N=200;

%% calculations
E_const = zeros(length(y),length(z_const));
I_const = currentsOfPlaneWave(constructive_waveguide,k,k_y,N,eta,r_eff,Lambda);
H_const_y = E_const;
H_const_z = E_const;
for indexy=1:length(y)
    for indexz=1:length(z_const)
        E_const(indexy,indexz) = find_generated_field(y(indexy),z_const(indexz),constructive_waveguide,k,Lambda,eta,I_0,k_y,I_const,N,y_0,z_0);
        [H_const_y(indexy,indexz), H_const_z(indexy,indexz)] = find_generated_Hyz_field(y(indexy),z_const(indexz),constructive_waveguide,k,Lambda,eta,I_0,k_y,I_const,N,y_0,z_0);% to write
    end
    indexy % to see where we are in the simulation
end
E_reg_const = ElectricFieldOfCurrent(I_0,y_0,z_0,y',z_const,k,eta);
[Hy_reg,Hz_reg] = MagneticFieldOfCurrent(I_0,y_0,z_0,y',z_const,k);
E_tot = E_const + E_reg_const;
Hy_tot=H_const_y+Hy_reg;
Hz_tot=H_const_z+Hz_reg;
S_y=-1/2*real(E_tot.*conj(Hz_tot));
S_z=1/2*real(E_tot.*conj(Hy_tot));

%% plots and results
P_z = sum(S_z,1);
figure;
quiver(y,z_const,S_y,S_z);
