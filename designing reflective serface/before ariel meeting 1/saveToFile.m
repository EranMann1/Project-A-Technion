%%
filename='N=505.mat';

%%
optimInput=solution;
D2 = optimInput(1);
D3 = optimInput(2);
D4 = optimInput(3);
H2 = optimInput(4);
H3 = optimInput(5);
H4 = optimInput(6);

Is2 = double(subs(Is, [d h],[[0 D2 D3 D4] [0 H2 H3 H4]]));
Fields = FindElectricFiield(Lambda,k,eta,theta_in,Ein,[0 D2 D3 D4],[0 H2 H3 H4],Is2,reff,N);
Impadences= Fields'./Is2;
Ds = [0, solution(1:3)];
Hs = [0, solution(4:6)];
save(filename, 'Ds','Hs' ,'Impadences');
