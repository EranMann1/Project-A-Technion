clear
% simultion of 4 metagriding
%% constants
eta = 377; %ohm
f= 10*10^9; %Hz
c = 3*10^8; %m/s
lambda = c/f;
k = 2*pi/lambda;
omega = 2*pi*f;
w = 0.0254e-3; % meters, 1 mil
r_eff = w/4;

%% paramters
theta_in = 15*pi/180; %rads
theta_out = 60*pi/180;% rads
phi = 25*pi/180; % rads
Ein = 1;

sinin = sin(theta_in);
cosin = cos(theta_in);
sinout = sin(theta_out);
cosout = cos(theta_out);

Zin = eta/cos(theta_in);
Zout = eta/cos(theta_out);
Lambda = lambda/abs(sinout-sinin);

%% vars
I = sym('I', [4 1]);
d = sym('d', [1,4]);
h = sym('h', [1,4]);
syms sin_in sin_out cos_in cos_out Z_in Z_out phase Lam K;
% phase = exp(1j*phi)


%% equations  - current as a function of everything
eq1 = exp(1j*K*(sin_in*d+cos_in*h))*I == 2*Lam/Z_in;
eq2 = exp(1j*K*(sin_out*d+cos_out*h))*I == 0;
eq3 = exp(1j*K*(sin_in*d-cos_in*h))*I == 0;
eq4 = exp(1j*K*(sin_out*d-cos_out*h))*I == 2*Lam*phase/Z_out;
eqns = [eq1,eq2,eq3,eq4];
Currents = solve(eqns,I);

%% iterations
N = 10;
d2 = linspace(0,Lambda,N+1);
d3 = d2;
d4 = d2;
h2 = linspace(0,lambda,N+1);
h3 = h2;
h4 = h2;

matrix = [d2 ; d3 ; d4 ; h2 ; ];

for [D2,D3,D4,H2,H3,H4]=[d2,d3,d4,h2,h3,h4]
   if H3>H2 && H4>H3
       
   end
end

for [i1,i2]=[1:10,1:10]
    i1
    i2
end