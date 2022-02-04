function [TotalElectricFieldMatrix, GeneratedElectricFieldMatrix] = ElectricFieldYZ2(Zs,Ds,Hs,k,Lambda,y,z,theta_in,Ein,N,eta)
    % the function calculte and sidply the elctric field Ex as a function
    % of Y and Z
    % Y is a colum vector Z is a row vector !!!
    sinin = sin(theta_in);
    cosin = cos(theta_in);
    EinYZ = Ein*exp(-1j*k*(sinin*y+cosin*z));
    %EIsYz = zeros(length(y), length(z), length(Is));
    EIsYz2 = zeros(length(y), length(z));
    
    Is = Ein*exp(-1j*k*(sinin*Ds'+cosin*Hs'))./Zs;
    
    n = (-N:N); % 
    phase = k*sinin*Lambda;
    for index = 1:length(Is)
        betaN = sqrt(k^2-(k*sinin+n*2*pi/Lambda).^2)'; 
        matrixY = exp(-1j*(y-Ds(index))*(2*pi*n-phase)/Lambda);
        matrixZ = exp(-1j*betaN*(abs(z-Hs(index))))./betaN;
        EIsYz2 = EIsYz2 - eta*k/(2*Lambda)*Is(index)*matrixY*matrixZ;
    end
    GeneratedElectricFieldMatrix = EIsYz2;
    TotalElectricFieldMatrix = EinYZ + GeneratedElectricFieldMatrix;
end