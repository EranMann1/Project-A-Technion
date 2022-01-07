
function best=reflective(thicknes)
    %% constants
    eta = 377; %ohm
    f= 10*10^9; %Hz
    c = 3*10^8; %m/s
    lambda = c/f;
    k = 2*pi/lambda;
    omega = 2*pi*f;
    w = 3*0.0254e-3; % meters, 3 mil
    reff = w/4;
    epsilon=8.85e-12;
    t=thicknes;

    %% paramters

    theta_in = 10*pi/180; %rads
    theta_out = 70*pi/180;% rads
    phi = 25*pi/180; % rads
    Ein = 1;

    sinin = sin(theta_in);
    cosin = cos(theta_in);
    sinout = sin(theta_out);
    cosout = cos(theta_out);

    Zin = eta/cos(theta_in);
    Zout = eta/cos(theta_out);
    Eout = Ein*sqrt(cosin/cosout);
    Lambda = lambda/abs(sinout-sinin);

    deltaRes=18.3e-3*eta/lambda; %resistence of copper wire


    %% Optimization parameters
    numOfIter=5;% how maney random starting conditions to check
    maxNumOfItter=50;
    N=1000; % to what degree to calculate the infinate sums 
    weigt=0.5;%how much to consider power lost and how much 
    save('var.mat','t','Lambda','Ein','Eout','Zin','Zout','theta_in','theta_out','k','phi','eta','reff','N','weigt','deltaRes');
    best.Zs=zeros(6,1);% the best results with only capacitance tipe impidances
    best.PowerLost=inf;
    best.CopperPowerLost=inf;
    best.Ds=zeros(1,6);
    best.Hs=zeros(1,6);
    best.Is=FindCurrents(best.Ds,best.Hs,Lambda,Ein,Eout,Zin,Zout,theta_in,theta_out,k,phi);
    best.lambda=lambda;
    best.eta=eta;
    best.theta_in=theta_in;
    best.theta_out=theta_out;
    best.r_eff=reff;
    best.Lambda=Lambda;

    %% Optimization
    i=1; %index
    total=0;% for percentage of good ones
    while and(i<=numOfIter,or(i==1,total<=maxNumOfItter)) 
        total=total+1;
        % initial position
        x0=[Lambda*rand(1,5),lambda*rand(1,5)];
        Ds= [0 x0(1:5)];
        Hs= [0 x0(6:10)];
        Is =  FindCurrents(Ds,Hs,Lambda,Ein,Eout,Zin,Zout,theta_in,theta_out,k,phi);
        Fields = FindElectricField(Lambda,k,eta,theta_in,theta_out,Ein,Ds,Hs,Is,reff,N,'Poisson');
        Zs = Fields./Is;
        if imag(Zs)<0
            % Set nondefault solver options
            options = optimoptions('fmincon','Display','off');

            % Solve
            [solution,objectiveValue] = fmincon(@objectiveFcn,x0,[],[],[],[],[],[],...
                @constraintFcn,options);

            % Clear variables
            clearvars options
            Ds = [0 solution(1:5)];
            Hs = [0 solution(6:10)];
            Is =  FindCurrents(Ds,Hs,Lambda,Ein,Eout,Zin,Zout,theta_in,theta_out,k,phi);
            Fields = FindElectricField(Lambda,k,eta,theta_in,theta_out,Ein,Ds,Hs,Is,reff,N,'Poisson');
            Zs = Fields./Is;
            objective = objectiveFcn(solution);
            if imag(Zs)<0
                 i=i+1;%go on with the index
            end
            if objective<best.PowerLost+weigt*best.CopperPowerLost % is it the best result so far?

              if imag(Zs)<0
                  best.PowerLost=1/2*(abs(Is).^2).*real(Zs);
                  best.CopperPowerLost=1/2*abs(Is).^2*deltaRes;
                  best.Ds=[0,solution(1:5)];
                  best.Hs=[0,solution(6:10)];
                  best.Is=FindCurrents(best.Ds,best.Hs,Lambda,Ein,Eout,Zin,Zout,theta_in,theta_out,k,phi);
                  best.Zs=FindElectricField(Lambda,k,eta,theta_in,theta_out,Ein,best.Ds,best.Hs,best.Is,reff,N,'Poisson')./best.Is;
              end
            end
        end
    end


end

    %% optimization function
    function f = objectiveFcn(optimInput)

    load('var.mat');
    % Edit the lines below with your calculation

    Ds= [0 optimInput(1:5)];
    Hs= [0 optimInput(6:10)];
    Is = FindCurrents(Ds,Hs,Lambda,Ein,Eout,Zin,Zout,theta_in,theta_out,k,phi);
    Fields = FindElectricField(Lambda,k,eta,theta_in,theta_out,Ein,Ds,Hs,Is,reff,N,'Poisson');
    Zs= Fields./Is;
    f= sum(abs((1/2*abs(Is).^2).*real(Zs)))+weigt*sum(abs((1/2*abs(Is).^2).*deltaRes));% minimum power lost (\gain)
    %f= sum(abs(real(Zs)./imag(Zs)));
    end

    %% costrains

    function [c,ceq] = constraintFcn(optimInput)
    load('var.mat');

    % Note, if no equality constraints, specify ceq = []
    D2 = optimInput(1);
    D3 = optimInput(2);
    D4 = optimInput(3);
    D5 = optimInput(4);
    D6= optimInput(5);
    H2 = optimInput(6);
    H3 = optimInput(7);
    H4 = optimInput(8);
    H5 = optimInput(9);
    H6 = optimInput(10);


    c(1) = H2 - H3;
    c(2) = H3 - H4;
    c(3) = H4 - H5;
    c(4) = H5 - H6;
    c(5) = - H2;
    c(6) = - D2;
    c(7) = - D3;
    c(8) = - D4;
    c(9) = - D5;
    c(10) = - D6;
    c(11) = D2 - Lambda;
    c(12) = D3 - Lambda;
    c(13) = D4 - Lambda;
    c(14) = D5 - Lambda;
    c(15) = D6 - Lambda;
    c(16) = H6 - t;
    ceq = [];
    end



