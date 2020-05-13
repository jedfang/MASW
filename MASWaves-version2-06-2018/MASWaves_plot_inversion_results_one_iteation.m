% MASWaves Inversion
% Version: 07.2017
%%
%  MASWaves_plot_inversion_results_one_iteation(c_t,lambda_t,...
%    c_curve0,lambda_curve0,c_curve0_up,lambda_curve0_up,c_curve0_low,...
%    lambda_curve0_low,n,beta,h,e,up_low_boundaries,FigWidth,FigHeight,FigFontSize)
%
%%
%  The function MASWaves_plot_inversion_results_one_iteation is used to plot
%  (i) the theoretical and experimental fundamental mode dispersion curves,
%  with or without the upper/lower experimental boundaries and (ii) the
%  inverted shear wave velocity profile
%
%  The dispersion curve is presented as Rayleigh wave phase velocity
%  vs. wavelength.
%
%  The shear wave velocity profile is presented as shear wave velocity
%  vs. depth
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
%  n                 Number of finite thickness layers
%  beta              Shear wave velocity [m/s] (vector of length n+1)
%  h                 Layer thicknesses [m] (vector of length n)
%  e                 Misfit [%]
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
function MASWaves_plot_inversion_results_one_iteation(c_t,lambda_t,...
    c_curve0,lambda_curve0,c_curve0_up,lambda_curve0_up,c_curve0_low,...
    lambda_curve0_low,n,beta,h,e,up_low_boundaries,FigWidth,FigHeight,FigFontSize)


% Plot experimental and theoretical dispersion curves
subplot(1,2,1), hold on

% With upper/lower boundaries
if strcmp(up_low_boundaries,'yes')
    obs = plot(c_curve0,lambda_curve0,'ko-','MarkerSize',3,...
        'MarkerFaceColor','k','MarkerEdgeColor','k');
    obs_up = plot(c_curve0_up,lambda_curve0_up,'k+--','MarkerSize',3,...
        'MarkerFaceColor','k','MarkerEdgeColor','k');
    calc = plot(c_t,lambda_t,'r-','MarkerSize',3,...
        'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1.5);
    obs_low = plot(c_curve0_low,lambda_curve0_low,'k+--','MarkerSize',3,...
        'MarkerFaceColor','k','MarkerEdgeColor','k');
    legend([obs,obs_up,calc],{'Exp.','Exp. up/low','Theor.'},'location','southwest','FontSize',FigFontSize)
end

% Without upper/lower boundaries
if strcmp(up_low_boundaries,'no')
    obs = plot(c_curve0,lambda_curve0,'ko-','MarkerSize',3,...
        'MarkerFaceColor','k','MarkerEdgeColor','k');
    calc = plot(c_t,lambda_t,'r-','MarkerSize',3,...
        'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1.5);
    legend([obs,calc],'Exp','Theor.','location','southwest','FontSize',FigFontSize)
end

% Figure title
set(gca,'fontsize',FigFontSize,'fontweight','normal')
title(['Misfit: \epsilon = ', num2str(e), '%'])

% Axis labels and axis limits
xlabel('Rayleigh wave velocity [m/s]','fontsize',FigFontSize,'fontweight','normal','color','k')
ylabel('Wavelength [m]','fontsize',FigFontSize,'fontweight','normal','color','k')
axis ij
grid on,  box off
set(gca, 'XColor', 'k'); set(gca, 'YColor', 'k'); set(gca,'XMinorTick','on','YMinorTick','on'); set(gca,'TickDir','out')
xlim([100*(floor(min([c_curve0; c_t])/100)) 100*(ceil(max([c_curve0; c_t])/100))]), ylim([5*(floor(min(lambda_curve0)/5)) 5*(ceil(max(lambda_curve0)/5))])


% Compute the depth vector z
z = zeros(1,(n+1));
for i = 1:n
    z(i+1) = sum(h(1:i));
end

% Plot shear wave velocity profile
subplot(1,2,2), hold on
for i = 1:n
    plot([beta(i),beta(i)], [z(i) z(i+1)],'k-','LineWidth',1)
    plot([beta(i),beta(i+1)], [z(i+1),z(i+1)],'k-','LineWidth',1)
end
plot([beta(n+1),beta(n+1)], [z(n+1) z(n+1)+5],'k--','LineWidth',1)

% Axis labels and axis limits
set(gca,'fontsize',FigFontSize,'fontweight','normal')
axis ij
grid on, box off
set(gca, 'XColor', 'k'); set(gca, 'YColor', 'k'); set(gca,'XMinorTick','on','YMinorTick','on'); set(gca,'TickDir','out')
xlabel('Shear wave velocity, [m/s]','FontSize',FigFontSize,'Fontweight','normal')
ylabel('Depth, [m]','FontSize',FigFontSize,'Fontweight','normal')
xlim([100*(floor(min(beta)/100)) 100*(ceil(max(beta)/100))]), ylim([5*(floor(min(z)/5)) 5*(ceil(max(z)/5))])

% Size of figure
set(gcf,'units','centimeters');
pos=[5, 5, FigWidth, FigHeight];
set(gcf,'Position',pos);
end