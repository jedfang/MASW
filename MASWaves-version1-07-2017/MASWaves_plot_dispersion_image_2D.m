% MASWaves Dispersion
% Version: 07.2017
%%
%  [fplot, cplot, Aplot] = MASWaves_plot_dispersion_image_2D...
%     (f,c,A,fmin,fmax,resolution,FigWidth,FigHeight,FigFontSize)
%%
%  The function MASWaves_plot_dispersion_image_2D plots the two-dimensional 
%  dispersion image of the recorded wavefield. The slant-stacked amplitude 
%  (A) is presented in the frequency - phase velocity - normalized summed 
%  amplitude domain using a color scale.
%
%  MASWaves_plot_dispersion_image_2D plots the dispersion image between the 
%  limits [f_min,f_max,cT_min,cT_max]
%
%% Input
%  f             Frequency range used for analysis [Hz]
%  c             Velocity range used for analysis [m/s]
%  A             Summed (slant-stack) amplitude corresponding to
%                different combinations of f (omega = 2*pi*f) and cT
%  fmin          Lower limit of the frequency axis [Hz]
%  fmax          Upper limit of the frequency axis [Hz]
%  resoluton     Number of contour lines
%                - resolution = 100 is generally recommended
%  FigWidth      Width of figure [cm]
%  FigHeight     Height of figure [cm]
%  FigFontSize   Font size for axis labels [pt]
%
%% Output
%  Two dimensional dispersion image with limits [f_min,f_max,cT_min,cT_max]
%
%  fplot         Frequency range of the dispersion image [Hz] 
%  cplot         Velocity range of the dispersion image [m/s] 
%  Aplot         Summed (slant-stack) amplitudes corresponding to 
%                fplot and cplot
%
%% Subfunctions
%  (None)
%
%%
function [fplot, cplot, Aplot] = MASWaves_plot_dispersion_image_2D(f,c,A,fmin,fmax,resolution,FigWidth,FigHeight,FigFontSize)

% Limits of frequency axis
[~,no_fmin] = (min(abs(f(:,1)-fmin)));
[~,no_fmax] = (min(abs(f(:,1)-fmax)));

% Select data corresponding to frequnecy range [fmin,fmax]
% Compute absolute value (length) of complex numbers
Aplot = A(no_fmin:no_fmax,:); 
fplot = f(no_fmin:no_fmax,:);
cplot = c(no_fmin:no_fmax,:);

% Plot the 2D dispersion image
[~,ch] = contourf(fplot,cplot,Aplot,resolution);
set(ch,'LineStyle','none');
set(ch,'edgecolor','none');
colormap(jet)
shading flat
grid on

% Axis limits and axis labels
set(gca,'XTick',0:10:fmax+0.01);
set(gca,'FontSize',FigFontSize);
xlim([fmin fmax])
xlabel('Frequency [Hz] ','FontSize',FigFontSize,'Fontweight','normal','color','k')
ylabel('Phase velocity [m/s] ','FontSize',FigFontSize,'Fontweight','normal','color','k')

% Size of figure
set(gcf,'units','centimeters')
pos = [2, 2, FigWidth, FigHeight]; 
set(gcf,'Position',pos)
box off
set(gca,'TickDir','out')

% Colorbar
c = colorbar('location','NorthOutside','color','k');
c.FontSize = FigFontSize;
c.Label.String = 'Normalized amplitude';
c.Label.FontSize = FigFontSize;

% Width and location of colorbar
axpos = get(gca,'Position');
cpos = get(c,'Position');
cpos(2) = cpos(2) - 0.5*cpos(4);
cpos(4) = 0.5*cpos(4);
set(c,'Position',cpos)
set(gca,'Position',axpos)

end