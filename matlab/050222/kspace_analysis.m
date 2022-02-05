function [detMatrix] = detCurrentMatrixKyH(k,ky,reff,Lambda,eta,Norder,Hs,Ds,Zs,H)
%% DETCURRENTMATRIXKYH k-space analysis of PEC waveguide
    NumLayers = length(Hs);
    Matrix = zeros(NumLayers,NumLayers);
    phi = -ky*Lambda;
    for Mindex = 1:NumLayers
        for mindex = 1:NumLayers
            if (Mindex == mindex)
                Matrix(Mindex,mindex) = Zs(Mindex)-k*eta/(2*Lambda)*sumElemnts123(phi,Norder,reff,Lambda,Hs,Ds,H);
            else
                Matrix(Mindex,mindex) = -k*eta/(2*Lambda)*sumElemnts4(phi,Norder,Hs,Ds,H);
            end 
        end
    end
    detMatrix = det(Matrix);
end
%%
function [Alphas,Betas] = calcAlphasBetas(k,Lambda,phi,IndexVector)
    Alphas = (2*pi*IndexVector-phi)/Lambda;
    Betas  = conj(sqrt(k^2-Alphas.^2));
end