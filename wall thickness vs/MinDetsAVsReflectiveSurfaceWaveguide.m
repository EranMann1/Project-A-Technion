
clear;
close all;

load('C:\Users\user\OneDrive - Technion\Documents\GitHub\Project-A-Technion\designing reflective serface\after ariel meating 3\realy small thickness surface.mat','lambda');
title='faster.mat';
thickness=(lambda/10:lambda/10:3*lambda);
real_thickness=thickness;
detA=zeros(1,length(thickness));
load('surf.mat');
for i=1:length(thickness)
    disp(['starting iteration ' num2str(i) ' out of ' num2str(length(thickness))]);
    surface = reflective(thickness(i));
    surfaces.num=surfaces.num+1;
    surfaces.array(surfaces.num)=surface;
    save('surf.mat','surfaces');
    disp('Second function of itteration');
    detA(i) = MinimumDetA(surface);
    surfaces.results(surfaces.num)=detA(i);
    real_thickness(i)=max(surface.Hs);
    save(title);
end

figure
plot(real_thickness,detA);
xlabel('max reflective layer thickness')
ylabel('det(A) of best conductive mode');
    
save(title);