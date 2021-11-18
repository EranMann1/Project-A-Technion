
 Is2 = double(subs(Is, [d h],[[0 solution(1:3)],[0 solution(4:6)]]));
 N = (1:1000);
 Fieldsrts = zeros(4,1000);
for n = N
    Fieldsrts(:,n) = FindElectricFiield(Lambda,k,eta,theta_in,Ein,[0 solution(1:3)],[0 solution(4:6)],Is2,reff,n);
end
figure(1);
plot(N,abs(Fieldsrts(1:4,:)))
grid on;