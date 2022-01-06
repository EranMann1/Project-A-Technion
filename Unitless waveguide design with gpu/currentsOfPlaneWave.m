function [Currents] = currentsOfPlaneWave(waveguide,k,k_y,N,eta,r_eff,Lambda) 
    % calc the currnets in each particle for every ky in k_y (vector 1D);
    % Currents is a 2D length(k_y) by number of particles
    % waveguide contains Zs, Hs, Ds are in the second dimention
    % k_y is in the first dimention
    % N is the order of summation
    
    Zs = waveguide.Zs;
    Ds = waveguide.Ds;
    Hs = waveguide.Hs;
    Ns = zeros(1,1,1,2*N+1);
    Ns_no_zero = zeros(1,1,1,2*N);
    Ky=zeros(1,1,length(k_y));
    Ky(:)=k_y; % 3th dim
    MatrixA = zeros(length(Ds),length(Hs),length(k_y));
    Currents = zeros(length(k_y),length(Ds));
   
    Ns(1,1,1,:) =(-N:N);
    Ns_no_zero(1,1,1,:) = [(-N:-1) (1:N)];
    
    Alphas = 2*pi*Ns/Lambda+Ky; % 4D (1,1,Ky,Ns);
    Alphas_no_zeros = 2*pi*Ns_no_zero/Lambda+Ky; % 4D (1,1,Ky,Ns_no_zero);
    
    Betas = conj(sqrt(k^2-Alphas.^2)); % 4D (1,1,Ky,Ns);
    Betas_no_zeros = conj(sqrt(k^2-Alphas_no_zeros.^2)); % 4D (1,1,Ky,Ns);
    
    Sum4D = sum(1./(Lambda*Betas_no_zeros)-1j./(2*pi*abs(Ns_no_zero)),4);
    
    Vec_diag = Zs-1j*pi*(1/pi*log(2*pi*r_eff/Lambda)+(1./(Lambda*conj(sqrt(4*pi^2-Ky.^2)))+Sum4D));
    for index = 1:12
        MatrixA_diag(index,index,:) = Vec_diag(1,index,:);
    end
    
    MatrixDs = transpose(Ds)-Ds;
    MatrixHs = abs(transpose(Hs)-Hs);
    MatrixA_nondiag = pi/Lambda*sum(exp(-1j*Alphas.*MatrixDs).*exp(-1j*Betas.*MatrixHs)./Betas,4);
    %MatrixA_nondiag(isnan(MatrixA_nondiag)||isinf(MatrixA_nondiag))= Vec_diag;
    MatrixA_nondiag(isinf(MatrixA_nondiag)) = 0;
    MatrixA_nondiag(isnan(MatrixA_nondiag)) = 0;
    
    MatrixA = MatrixA_diag+MatrixA_nondiag;
    PW_vec = exp(-1j*(conj(sqrt(4*pi^2-Ky.^2)).*transpose(Hs)+Ky.*transpose(Ds)));
    
    
    for index = 1:length(k_y)
        temp_mat(:,:) = MatrixA(:,:,index);
        temp_PW_vec(:) = PW_vec(:,1,index);
        Currents(index,:) = temp_PW_vec/temp_mat;
    end

%     Ns_no_zero = zeros(1,1,2*N);
%     Ns(:) = (-N:1:N);
%     Ns_no_zero(:) = [(-N:1:-1) (1:1:N)];
%     Alphas = 2*pi*Ns/Lambda + k_y;
%     Alphas_no_zero = 2*pi*Ns_no_zero/Lambda + k_y;
%     Betas = conj(sqrt(k^2-Alphas.^2));
%     Betas_no_zero = conj(sqrt(k^2-Alphas_no_zero.^2));
%     k_z=sqrt(k^2-k_y.^2);
%     Matrix = zeros(length(Zs),length(Zs));
%     Vector = zeros(1,length(Zs));
%     for index0=1:length(k_y)
%        Vector(:) = exp(-1j.*(k_z(index0).*Hs+k_y(index0).*Ds));
%        for index1=1:length(Zs)
%            Matrix(index1,index1) = Zs(index1) + k*eta/2*(-1j/pi*log(2*pi*r_eff/Lambda)+1/Lambda/Betas(N+1)+sum(1/Lambda/Betas_no_zero(index0,1,:)-1j/2/pi/Ns_no_zero,'all'));
%            for index2=1:length(Zs)
%                Matrix(index1,index2) = k*eta/2/Lambda*sum(exp(-1j*(Alphas(index0,1,:)*(Ds(index1)-Ds(index2))+Betas(index0,1,:)*abs(Hs(index1)-Hs(index2))))./Betas(index0,1,:),'all');
%            end
%        end
%        Currents(index0,:) = Vector*inv(Matrix); 
%     end
end