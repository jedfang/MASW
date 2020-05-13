% MASWaves Dispersion
% Version: 07.2017
%%
%  MASWaves_plot_dispersion_curve(f_curve0,c_curve0,lambda_curve0,...
%     f_curve0_up,c_curve0_up,lambda_curve0_up,f_curve0_low,c_curve0_low,lambda_curve0_low,...
%     type,up_low_boundaries,FigWidth,FigHeight,FigFontSize)
%
%%
%  The function MASWaves_plot_dispersion_curve is used to plot the
%  fundamental mode dispersion curve, with or without upper/lower
%  boundaries. The dispersion curve is either presented as frequency vs.
%  Rayleigh wave velocity or as Rayleigh wave velocity vs. wavelength.
%
%% Input
%  Fundamental mode dispersion curve
%  f_curve0          Frequency [Hz]
%  c_curve0          Rayleigh wave velocity [m/s]
%  lambda_curve0     Wavelength [m]
%
%  f_curve0_up       Frequency, upper bound curve [Hz]
%                    (Can be assigned as 'nan' or [] if upper/lower boundaries
%                    are not wanted.)
%  c_curve0_up       Rayleigh wave velocity, upper bound curve [m/s]
%                    (Can be assigned as 'nan' or [] if upper/lower boundaries
%                    are not wanted.)
%  lambda_curve0_up  Wavelength, upper bound curve [m]
%                    (Can be assigned as 'nan' or [] if upper/lower boundaries
%                    are not wanted.)
%  f_curve0_low      Frequency, upper bound curve [Hz]
%                    (Can be assigned as 'nan' or [] if upper/lower boundaries
%                    are not wanted.)
%  c_curve0_low      Rayleigh wave velocity, upper bound curve [m/s]
%                    (Can be assigned as 'nan' or [] if upper/lower boundaries
%                    are not wanted.)
%  lambda_curve0_low Wavelength, upper bound curve [m]
%                    (Can be assigned as 'nan' or [] if upper/lower boundaries
%                    are not wanted.)
%  type              Controls how the dispersion curve is presented
%                    - 'f_c'      Frequency vs. Rayleigh wave velocity
%                    - 'c_lambda' Rayleigh wave velocity vs. wavelength
%  up_low_boundaries - 'yes'      Upper/lower boundaries for the fundamental
%                                 mode dispersion curve are wanted.
%                    - 'no'       Upper/lower boundaries for the fundamental
%                                 mode dispersion curve are not wanted.
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
function  MASWaves_plot_dispersion_curve(f_curve0,c_curve0,lambda_curve0,...
    f_curve0_up,c_curve0_up,lambda_curve0_up,f_curve0_low,c_curve0_low,lambda_curve0_low,...
    type,up_low_boundaries,FigWidth,FigHeight,FigFontSize)

%Frequency vs. Rayleigh wave phase velocity
if strcmp(type,'f_c')
    
    % With upper/lower boundaries
    if strcmp(up_low_boundaries,'yes')
        hold on
        plot(f_curve0,c_curve0,'ko-','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k')
        plot(f_curve0_up,c_curve0_up,'r+--','MarkerSize',3,'MarkerFaceColor','r','MarkerEdgeColor','r')
        plot(f_curve0_low,c_curve0_low,'r+--','MarkerSize',3,'MarkerFaceColor','r','MarkerEdgeColor','r')
        legend({'Exp.','Exp. up/low'},'location','northeast','FontSize',FigFontSize)
        hold on
    end
    
    % Without upper/lower boundaries
    plot(f_curve0,c_curve0,'ko-','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k')
    if strcmp(up_low_boundaries,'no')
        legend('Exp.','location','northeast','FontSize',FigFontSize)
    end
    
    % Axis labels and axis limits
    set(gca, 'FontSize', FigFontSize)
    grid on
    xlabel('Frequency [Hz]','FontSize',FigFontSize,'Fontweight','normal')
    ylabel('Rayleigh wave velocity [m/s]','FontSize',FigFontSize,'Fontweight','normal')
    
    % Size of figure
    set(gcf,'units','centimeters')
    pos = [2, 2, FigWidth, FigHeight];
    set(gcf,'Position',pos)
    box off
    set(gca,'TickDir','out')
end

% Rayleigh wave phase velocity vs. wavelength                     
if strcmp(type,'c_lambda')
    
    % With upper/lower boundaries
    if strcmp(up_low_boundaries,'yes')
        hold on
        plot(c_curve0,lambda_curve0,'ko-','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k')
        plot(c_curve0_up,lambda_curve0_up,'r+--','MarkerSize',3,'MarkerFaceColor','r','MarkerEdgeColor','r')
        plot(c_curve0_low,lambda_curve0_low,'r+--','MarkerSize',3,'MarkerFaceColor','r','MarkerEdgeColor','r')
        legend({'Exp.','Exp. up/low'},'location','southwest','FontSize',FigFontSize)
        hold on
    end
    
    % Without upper/lower boundaries
    plot(c_curve0,lambda_curve0,'ko-','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k')
        if strcmp(up_low_boundaries,'no')
        legend('Exp.','location','southwest','FontSize',FigFontSize)
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
end
hold off
end
