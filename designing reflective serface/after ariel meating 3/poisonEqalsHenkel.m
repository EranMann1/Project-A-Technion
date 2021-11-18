load('N=500 good results.mat');
n=(10:2000);
H = zeros(6,length(n));
P = zeros(6,length(n));
for index=1:length(n)
    H(:,index) = FindElectricField(Lambda,k,eta,theta_in,theta_out,Ein,Ds,Hs,Is,reff,n(index),'Henkels');
    P(:,index) = FindElectricField(Lambda,k,eta,theta_in,theta_out,Ein,Ds,Hs,Is,reff,n(index),'Poisson');
end
figure(1);
hold on;
grid on;
plot(n,real(H(1,:)));
plot(n,real(P(1,:)));
plot(n,imag(H(1,:)));
plot(n,imag(P(1,:)));
legend('r(H)','r(p)','i(H)','i(p)')
hold off;
