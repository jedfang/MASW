%  MASWaves Dispersion
%  Version: 07.2017
%% 
%  MASWaves_plot_data(u,N,dx,x1,L,T,Tmax,du,FigWidth,FigHeight,FontSize)
%% 
%  The function MASWaves_plot_data plots recorded multichannel surface wave data 
%  in the offset-time domain
%
%% Input
%  u             u(x,t) offset-time shot gather
%  N             Number of receivers
%  dx            Receiver spacing [m]
%  x1            Source offset [m]
%  L             Length of receiver spread [m]
%  T             Time of individual recordings [s]
%  Tmax          Total recording time [s] 
%  du            Scale factor for offset between traces 
%  FigWidth      Width of figure [cm]
%  FigHeight     Height of figure [cm]
%  FigFontSize   Font size for axis labels [pt]
%
%% Output
%  Plot of recorded surface wave data in the offset-time domain
%
%% Subfunctions
%  (None)
%
%%
function MASWaves_plot_data(u,N,dx,x1,L,T,Tmax,du,FigWidth,FigHeight,FigFontSize)

for j = 1:N
    signal_plot(:,j) = u(:,j) + (j-1)*du*dx;
    plot(signal_plot(:,j),T,'k');
    hold on
end
set(gca,'YDir','reverse')
set(gca, 'FontSize', FigFontSize)

% Axis limits and axis labels
axis([-5*dx*du,L*du+5*dx*du,0,Tmax]) 
set(gca,'XTick',0:6*dx*du:L*du,'XTickLabel',x1:(6*dx):x1+L)
set(gca,'TickDir','out')
box off
xlabel('Distance from source [m]','FontSize',FigFontSize)
ylabel('Time [s]','FontSize',FigFontSize)

% Size of figure
set(gcf,'units','centimeters')
pos = [2, 2, FigWidth, FigHeight]; 
set(gcf,'Position',pos)

hold off
end
