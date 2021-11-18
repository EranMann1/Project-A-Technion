function [TotalElectricFieldMatrix, GeneratedElectricFieldMatrix] = ElectricFieldYZ(Is,Ds,Hs,k,Lambda,y,z,theta_in,Ein,N,eta,theta_out)
    % the function calculte and sidply the elctric field Ex as a function
    % of Y and Z
    % Y is a colum vector Z is a row vector !!!
    sinin = sin(theta_in);
    cosin = cos(theta_in);
    sinout = sin(theta_out);
    cosout = cos(theta_out);
    EinYZ = Ein*exp(-1j*k*(sinin*y+cosin*z));
    
    n =(-N:N);
    Alphas = k*sinin+n*k*(sinout-sinin); 
    Betas = sqrt(k^2-Alphas.^2)';
    
    GeneratedElectricFieldMatrix = zeros(length(y),length(z));
    
    for index =1:length(Is)
        matrixY = exp(-1j*(y-Ds(index))*Alphas);
        matrixZ = exp(-1j*Betas*abs(z-Hs(index)))./Betas;
        GeneratedElectricFieldMatrix = GeneratedElectricFieldMatrix-Is(index)*k*eta/(2*Lambda)*matrixY*matrixZ;
    end
    TotalElectricFieldMatrix = EinYZ + GeneratedElectricFieldMatrix;
    
    %EIsYz = zeros(length(y), length(z), length(Is));
%     EIsYz2 = zeros(length(y), length(z));
%     n = (-N:N); % 
%     phase = -k*sinin*Lambda;
%     for index = 1:length(Is)
%         betaN = sqrt(k^2-(k*sinin+(n*2*pi-phase)/Lambda).^2)'; 
%         matrixY = exp(-1j*(y-Ds(index))*(2*pi*n-phase)/Lambda);
%         matrixZ = exp(-1j*betaN*(abs(z-Hs(index))))./betaN;
%         EIsYz2 = EIsYz2 - eta*k/(2*Lambda)*Is(index)*matrixY*matrixZ;
%     end
%     GeneratedElectricFieldMatrix = EIsYz2;
%     TotalElectricFieldMatrix = EinYZ + GeneratedElectricFieldMatrix;
end