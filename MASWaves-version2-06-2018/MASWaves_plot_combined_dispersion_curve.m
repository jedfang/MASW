% MASWaves Combination
% Version: 06.2018
%%
%  [c_m,c_ps,c_ms,lambda_m] = MASWaves_plot_combined_dispersion_curve...
%    (no_measurements,c_all,lambda_all,c_mean,c_plus_std,c_min_std,lambda_mean,...
%    PlotAll,FigWidth,FigHeight,FigFontSize)
%
%%
%  The function MASWaves_plot_combined_dispersion_curve plots the combined
%  mean dispersion curve along with upper and lower bound curves.
%
%  In addition, the elementary dispersion curve data points can be shown
%  along with the composite curves.
%
%% Input
%  Elementary dispersion data
%  no_measurements  Number of measurements
%  lambda_all       A cell array of length no_measurements, where
%                   - cell no. 1 contains the wavelength values of
%                     elementary dispersion curve no. 1
%                   - cell no. 2 contains the wavelength values of
%                     elementary dispersion curve no. 2
%                   - etc.
%  c_all            A cell array of length no_measurements, where
%                   - cell no. 1 contains the Rayleigh wave phase velocity
%                     values of elementary dispersion curve no. 1
%                   - cell no. 2 contains the Rayleigh wave phase velocity
%                     values of elementary dispersion curve no. 2
%                   - etc.
%
%  Combined dispersion curve
%  c_mean           Rayleigh wave velocity [m/s]
%  lambda_mean      Wavelength [m]
%  c_plus_std       Upper bound Rayleigh wave velocity
%                   [m/s] (Mean value plus one standard deviation)
%  c_min_std        Lower bound Rayleigh wave velocity
%                   [m/s] (Mean value minus one standard deviation)
%
%  PlotAll          - '0' Show mean/upper bound/lower bound curves.
%                   - '1' Show mean/upper bound/lower bound curves
%                         and the elementary dispersion data

%% Output
%  (none)
%
%% Subfunctions
%  (none)
%%
function MASWaves_plot_combined_dispersion_curve(no_measurements,c_all,lambda_all,...
    c_mean,c_plus_std,c_min_std,lambda_mean,PlotAll,FigWidth,FigHeight,FigFontSize)

figure, hold on

% Show elementary dispersion curves
if PlotAll == 1
    data = plot(c_all{1},lambda_all{1},'o','Color',[0.75 0.75 0.75],'MarkerFaceColor',[0.75 0.75 0.75],'MarkerSize',3);
    for i = 2:no_measurements
        plot(c_all{i},lambda_all{i},'o','Color',[0.75 0.75 0.75],'MarkerFaceColor',[0.75 0.75 0.75],'MarkerSize',3)
    end
end

mean_curve = plot(c_mean,lambda_mean,'ko-','MarkerFaceColor','k','MarkerSize', 2,'LineWidth',1);
boundary_curve = plot(c_plus_std,lambda_mean,'k--', 'LineWidth',1);
plot(c_min_std,lambda_mean,'k--', 'LineWidth',1)

% Axis labels and axis limits
set(gca, 'FontSize', FigFontSize)
axis ij
grid on
xlabel('Rayleigh wave velocity [m/s]','FontSize',FigFontSize,'Fontweight','normal')
ylabel('Wavelength [m]','FontSize',FigFontSize,'Fontweight','normal')

% Legend
if PlotAll == 0
    hleg1=legend([mean_curve,boundary_curve],'Mean','Mean \pm std','Location','SouthWest');
elseif PlotAll == 1
    hleg1=legend([mean_curve,boundary_curve data],'Mean','Mean \pm std', 'Data','Location','SouthWest');
end

set(hleg1,'Fontsize',FigFontSize)

% Size of figure
set(gcf,'units','centimeters')
pos = [2, 2, FigWidth, FigHeight];
set(gcf,'Position',pos)
box off
set(gca,'TickDir','out')

end