function WaveGuide = createWaveguide(D,H,ds,hs,zs)
    % create a sturct that represent the waveguide
    % gets the parameters of the reflective surface' and of the spacing.
    % ds, hs, zs are in the first dimention
    % Ds, Hs, Zs are in the second dimention
    
    s=length(ds);
    Ds = zeros(1,2*s);
    Hs = zeros(1,2*s);
    Zs = zeros(1,2*s);
    for index = 1:s
        Ds(index) = ds(s+1-index);
        Hs(index) = H + hs(s+1-index);
        Zs(index) = zs(s+1-index);
        Ds(s+index) = D - ds(index);
        Hs(s+index) = -H - hs(index);
        Zs(s+index) = zs(index);
    end
    WaveGuide.Ds=Ds;
    WaveGuide.Hs=Hs;
    WaveGuide.Zs=Zs;
end