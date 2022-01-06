function [dets] = detsOfKy(waveguide,k,k_y,N,eta,r_eff,Lambda) 
    % calc the currnets in each particle for every ky in k_y (vector 1D);
    % dets is a 2D length(k_y) by number of particles
    % waveguide contains Zs, Hs, Ds are in the second dimention
    % k_y is in the first dimention
    % N is the order of summation
    
    Zs = gpuArray(waveguide.Zs);
    Ds = gpuArray(waveguide.Ds);
    Hs = gpuArray(waveguide.Hs);
    Ns = gpuArray(zeros(1,1,1,2*N+1));
    Ns_no_zero = gpuArray(zeros(1,1,1,2*N));
    Ky=gpuArray(zeros(1,1,length(k_y)));
    Ky(:)=k_y; % 3th dim
    MatrixA = gpuArray(zeros(length(Ds),length(Hs),length(k_y)));
    Currents = gpuArray(zeros(length(k_y),length(Ds)));
   
    Ns(1,1,1,:) =gpuArray((-N:N));
    Ns_no_zero(1,1,1,:) = gpuArray([(-N:-1) (1:N)]);
    
    Alphas = 2*pi*Ns/Lambda+Ky; % 4D (1,1,Ky,Ns);
    Alphas_no_zeros = 2*pi*Ns_no_zero/Lambda+Ky; % 4D (1,1,Ky,Ns_no_zero);
    
    Betas = conj(sqrt(complex(k^2-Alphas.^2))); % 4D (1,1,Ky,Ns);
    Betas_no_zeros = conj(sqrt(complex(k^2-Alphas_no_zeros.^2))); % 4D (1,1,Ky,Ns);
    
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
    
    for index = 1:length(k_y)
        temp_mat(:,:) = MatrixA(:,:,index);
        dets(index)=det(temp_mat);
    end

end