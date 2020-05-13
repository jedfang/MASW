% MASWaves Inversion
% Version: 07.2017
%%
%  MASWaves_plot_theor_exp_dispersion_curves(c_t,lamda_t,...
%    c_curve0,lambda_curve0,c_curve0_up,lambda_curve0_up,c_curve0_low,...
%    lambda_curve0_low,up_low_boundaries,FigWidth,FigHeight,FigFontSize)
%
%%
%  The function MASWaves_plot_theor_exp_dispersion_curves is used to plot
%  the theoretical and experimental fundamental mode dispersion curves,
%  with or without the upper/lower experimental boundaries.
%
%  The dispersion curve is presented as Rayleigh wave phase velocity
%  vs. wavelength.
%
%% Input
%
%  Theoretical fundamental mode dispersion curve
%  c_t               Phase velocity [m/s]
%  lambda_t          Wavelength [m]
%
%  Experimental fundamental mode dispersion curve
%  c_curve0          Phase velocity [m/s]
%  lambda_curve0     Wavelength [m]
%  c_curve0_up       Phase velocity, upper bound curve [m/s]
%                    (Can be assigned as 'nan' or [] if upper/lower boundaries
%                    are not wanted.)
%  lambda_curve0_up  Wavelength, upper bound curve [m]
%                    (Can be assigned as 'nan' or [] if upper/lower boundaries
%                    are not wanted.)
%  c_curve0_low      Phase velocity, upper bound curve [m/s]
%                    (Can be assigned as 'nan' or [] if upper/lower boundaries
%                    are not wanted.)
%  lambda_curve0_low Wavelength, upper bound curve [m]
%                    (Can be assigned as 'nan' or [] if upper/lower boundaries
%                    are not wanted.)
%
%  up_low_boundaries - 'yes'      Upper/lower boundaries for the experimental 
%                                 fundamental mode dispersion curve are wanted.
%                    - 'no'       Upper/lower boundaries for the experimental 
%                                 fundamental mode dispersion curve are not wanted.
%  FigWidth          Width of figure [cm]
%  FigHeight         Height of figure [cm]
%  FigFontSize       Font size for axis labels [pt]
%
%% Output
%  (None)
%
%% Subfunctions
%  (None)
%
%%
function  MASWaves_plot_theor_exp_dispersion_curves(c_t,lambda_t,...
    c_curve0,lambda_curve0,c_curve0_up,lambda_curve0_up,c_curve0_low,...
    lambda_curve0_low,up_low_boundaries,FigWidth,FigHeight,FigFontSize)


% With upper/lower boundaries
if strcmp(up_low_boundaries,'yes')
    hold on
    plot(c_curve0,lambda_curve0,'ko-','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k')
    plot(c_curve0_up,lambda_curve0_up,'k+--','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k')
    plot(c_t,lambda_t,'r-','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1.5)
    plot(c_curve0_low,lambda_curve0_low,'k+--','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k')
    legend({'Exp.','Exp. up/low','Theor.'},'location','southwest','FontSize',FigFontSize)
    hold on
end

% Without upper/lower boundaries
plot(c_curve0,lambda_curve0,'ko-','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
plot(c_t,lambda_t,'r-','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1.5)

if strcmp(up_low_boundaries,'no')
    legend({'Exp.','Theor.'},'location','southwest','FontSize',FigFontSize)
end

% Axis labels and axis limits
set(gca, 'FontSize', FigFontSize)
axis ij
grid on
xlabel('Rayleigh wave velocity [m/s]','FontSize',FigFontSize,'Fontweight','normal')
ylabel('Wavelength [m]','FontSize',FigFontSize,'Fontweight','normal')

% Size of figure
set(gcf,'units','centimeters')
pos = [2, 2, FigWidth, FigHeight];
set(gcf,'Position',pos)
box off
set(gca,'TickDir','out')
hold off
end
