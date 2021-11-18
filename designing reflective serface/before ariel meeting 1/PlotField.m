function PlotField(Field,y,z)

    figure;
    hold on;
    surf(y,z,real(Field)');
    set(get(gca,'Children'),'EdgeColor','None'); colormap(jet(256)); colorbar;
    grid on;
    xlabel('y','Interpreter','latex'); ylabel('z','Interpreter','latex'); zlabel('$|E_x(y,z)|^2$','Interpreter','latex'); 
    title('Electric Field','Interpreter','latex');
    hold off;
end