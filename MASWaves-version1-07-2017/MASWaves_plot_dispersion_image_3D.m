% MASWaves Dispersion
% Version: 07.2017
%% 
%  [fplot,cplot,Aplot] = MASWaves_plot_dispersion_image_3D...
%    (f,c,A,fmin,fmax,FigWidth,FigHeight,FigFontSize)
%
%% 
%  The function MASWaves_plot_dispersion_image_3D plots the three-dimensional 
%  dispersion image of the recorded wavefield. The slant-stacked amplitude (A) 
%  is presented in the frequency - phase velocity - normalized summed 
%  amplitude domain. 
%
%  MASWaves_plot_dispersion_image_3D plots the dispersion image between the 
%  limits [f_min,f_max,cT_min,cT_max]
%  
%% Input
%  f             Frequency range used for analysis [Hz]
%  c             Velocity range used for analysis [m/s]
%  A             Summed (slant-stack) amplitude corresponding to 
%                different combinations of f and cT
%  fmin          Lower limit of the frequency axis [Hz]
%  fmax          Upper limit of the frequency axis [Hz]
%  FigWidth      Width of figure [cm]
%  FigHeight     Height of figure [cm]
%  FigFontSize   Font size for axis labels [pt]
% 
%% Output
%  Three dimensional dispersion image with limits [f_min,f_max,cT_min,cT_max] 
%
%  fplot         Frequency range of dispersion image [Hz] 
%  cplot         Velocity range of dispersion image [m/s] 
%  Aplot         Summed (slant-stack) amplitudes corresponding to 
%                fplot and cplot
%
%% Subfunctions 
%  (None)
%
%%
function [fplot,cplot,Aplot] = MASWaves_plot_dispersion_image_3D(f,c,A,fmin,fmax,FigWidth,FigHeight,FigFontSize)
        
% Limits of frequency axis
[~,no_fmin] = (min(abs(f(:,1)-fmin)));
[~,no_fmax] = (min(abs(f(:,1)-fmax)));

% Select data corresponding to frequnecy range [fmin,fmax]
% Compute absolute value (length) of complex numbers
Aplot = A(no_fmin:no_fmax,:); 
fplot = f(no_fmin:no_fmax,:);
cplot = c(no_fmin:no_fmax,:);

% Plot the 3D dispersion image
surf(fplot,cplot,Aplot);
colormap(jet)
shading interp
grid on

% Axis limits and axis labels
set(gca,'FontSize',FigFontSize);
xlabel('Frequency [Hz]','FontSize',FigFontSize,'Fontweight','normal')
ylabel('Phase velocity [m/s]','FontSize',FigFontSize,'Fontweight','normal')
zlabel('Normalized amplitude','FontSize',FigFontSize,'Fontweight','normal')

% Size of figure
set(gcf,'units','centimeters')
pos = [2, 2, FigWidth, FigHeight]; 
set(gcf,'Position',pos)
box off
set(gca,'TickDir','out')

% Colorbar
c = colorbar('location','NorthOutside','color','k');
c.FontSize = FigFontSize;

% Width and location of colorbar
axpos = get(gca,'Position');
cpos = get(c,'Position');
cpos(2) = cpos(2) - 0.5*cpos(4);
cpos(4) = 0.5*cpos(4);
set(c,'Position',cpos)
set(gca,'Position',axpos)

% View
az = 45;
el = 45;
view(az,el)
end