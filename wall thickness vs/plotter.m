load('surf');
j=1;
for i=1:length(surfaces.results)
   if surfaces.results(i)>0
      thickness(j)=max(surfaces.array(i).Hs);
      result(j)=surfaces.results(i);
      j=j+1;
   end 
end
figure;
plot(thickness,log(result),'.');