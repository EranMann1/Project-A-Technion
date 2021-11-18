function time = timeit(type,N)
load('N=500 good results.mat');
if type=='H'
    tic;
    FindElectricField(Lambda,k,eta,theta_in,theta_out,Ein,Ds,Hs,Is,reff,N,'Henkels');
    time = toc;
end
if type=='P'
    tic;
    FindElectricField(Lambda,k,eta,theta_in,theta_out,Ein,Ds,Hs,Is,reff,N,'Poisson');
    time = toc;
end
