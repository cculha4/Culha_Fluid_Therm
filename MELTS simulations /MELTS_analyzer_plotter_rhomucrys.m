%%%%%%%%%%%%%%%%%%%%%%%%
%%MELTS plotter
%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
load('data.mat')
%% prepare formating options
HA = {'HorizontalAlignment','left','center','right'};
VA = {'VerticalAlignment','bottom','middle','top'};
UN = {'Units','Normalized','Inches'};
TX = {'Interpreter','Latex'};
TL = {'TickLabelInterpreter','Latex'};
LW = {'LineWidth',1,1.25,1.5,2};
FS = {'FontSize',10,15,18,21,24};
MS = {'MarkerSize',6,8,12};
LS = {'LineStyle','-','--','-.',':'};
printing = 1; %print 1 not print 0

%density of the melt phase
figure
hold on
plot(T_liq(index_h:index_c),v_rho_liq,'--k','linewidth',2)
hold on
plot(T_liq(index_h:index_c),rho_liq(index_h:index_c),'.k','linewidth',10)
hold on
% title('\rho')
text(0.98,0.9,'Melt density, $\rho_m$ [kg/cm$^3$]',UN{[1,2]},TX{:},FS{[1,4]},VA{[1,2]},HA{[1,4]},'Color',[0 0 0]);%FigDemoColorMap(30,:));
xlabel('Temperature in $^{\circ}$C','interpreter','latex','FontSize',18)
 ylabel(' ','interpreter','latex','FontSize',18)
set(gca,'FontSize', 18,TL{:},FS{[1,3]});
if printing
print('rho_melt','-dpng')
end

figure
plot(T_liq(index_h:index_c),mu_liq(index_h:index_c),'.k','LineWidth',10)
hold on
plot(T_liq(index_h:index_c),v_mu_liq,'--k','LineWidth',2)
hold on
% title('\mu')
text(0.98,0.9,'Melt viscosity, $\mu$ [log Pa $\cdot$ s]',UN{[1,2]},TX{:},FS{[1,4]},VA{[1,2]},HA{[1,4]},'Color',[0 0 0]);%FigDemoColorMap(30,:));
xlabel('Temperature in $^{\circ}$C','interpreter','latex','FontSize',18)
 ylabel(' ','interpreter','latex','FontSize',18)
set(gca,'FontSize', 18,TL{:},FS{[1,3]});
if printing
print('mu_melt','-dpng')
end


figure
 plot(T_liq(index_h:index_c),vol_sol_n(index_h:index_c),'.k','LineWidth',10)
 hold on
plot(T_liq(index_h:index_c),v_vol_sol,'--k','LineWidth',2)
% title('Volume','interpreter','latex','FontSize',18)
text(0.98,0.9,'Crystallinity, $\phi$ ',UN{[1,2]},TX{:},FS{[1,4]},VA{[1,2]},HA{[1,4]},'Color',[0 0 0]);%FigDemoColorMap(30,:));
xlabel('Temperature in $^{\circ}$C','interpreter','latex','FontSize',18)
ylabel(' ','interpreter','latex','FontSize',18)
set(gca,'FontSize', 18,TL{:},FS{[1,3]});
if printing
print('vol_crystallinity','-dpng')
end
