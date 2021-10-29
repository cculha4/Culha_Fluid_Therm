%%%%%%%%%
%%MELTS Analyzer Plotter
%%%%%%%%%
close all
clear all
load('data.mat')
%% plotting
figure;
%subplot(2,3,1)
hold on;
plot(T_liq(index_h:index_c),rho_thermal(index_h:index_c),'m','linewidth',5)
plot(T_liq(index_h:index_c),rho_chemical(index_h:index_c),'r','linewidth',5)
plot(T_liq(index_h:index_c),rho_bub,'b','linewidth',5)
plot(T_liq(index_h:index_c),rho_cry,'g','linewidth',5)
plot(T_liq(index_h:index_c),RHO,'k','linewidth',5)
title('bulk density due to different contributions', 'interpreter','latex','FontSize',18)
xlabel('Temperature in $^{\circ}$ C','interpreter','latex','FontSize',18)
ylabel('$\overline{\rho}$ g/cm$^3$','interpreter','latex','FontSize',18)
th = sprintf('thermal'); %d\rho/dT = %.1E', drho_thermal);
ch = sprintf('chemical');% d\rho/dT = %.1E', drho_chemical);
bb = sprintf('bubble');% d\rho/dT = %.1E', drho_bub);
cry = sprintf('crystal');% d\rho/dT = %.1E', drho_cry);
r = sprintf('bulk density');% d\rho/dT = %.1E', dRHO);
legend(th,ch,bb,cry,r)
legend('Location','southeast')
set(gca,'FontSize', 18);
%tightfig;
%set(gca,'FontSize','14')

print('bulk_rho','-dpng')



%plotting
figure
% subplot(2,3,2)
hold on
plot(T_liq(index_h:index_c),v_rho_liq,'-r','linewidth',2)
plot(T_liq(index_h:index_c),v_rho_gas,'-b','linewidth',2)
plot(T_liq(index_h:index_c),v_rho_sol,'-k','linewidth',2)

hold on
plot(T_liq(index_h:index_c),rho_liq(index_h:index_c),'.r','linewidth',10)
plot(T_liq(index_g_h:index_c),rho_gas(index_g_h:index_c),'.b','linewidth',10)
plot(T_liq(index_h:index_c),rho_sol_tbl(index_h:index_c),'.k','linewidth',10)
hold on
title('\rho')
xlabel('Temperature in $^{\circ}$C','interpreter','latex','FontSize',18)
ylabel('$\rho$ g/cm$^3$','interpreter','latex','FontSize',18)
set(gca,'FontSize', 18);
legend('liquid','gas','solid')
legend('Location','southeast')
%tightfig;
%set(gca,'FontSize','14')
print('rho','-dpng')

figure
% subplot(2,3,[3 6])
plot(T_liq(index_h:index_c),mu_liq(index_h:index_c),'.r')
hold on
plot(T_liq(index_h:index_c),v_mu_liq)
hold on
title('\mu')
xlabel('Temperature in $^{\circ}$C','interpreter','latex','FontSize',18)
ylabel('log_{10} Pa $\cdot$s','interpreter','latex','FontSize',18)
set(gca,'FontSize', 18);
%tightfig;
% set(gca,'FontSize','14')
print('mu','-dpng')

figure
% subplot(2,3,4)
plot(T_liq(index_h:index_c),mass_liq(index_h:index_c),'.r')
hold on
plot(T_liq(index_h:index_c),v_mass_liq)
hold on
plot(T_liq(index_g_h:index_c),mass_gas(index_g_h:index_c),'.b')
hold on
plot(T_liq(index_h:index_c),v_mass_gas)
plot(T_liq(index_h:index_c),mass_sol_tbl(index_h:index_c),'.k')
plot(T_liq(index_h:index_c),v_mass_sol)
title('Mass')
legend('liquid','','gas','','solid','')
xlabel('Temperature in $^{\circ}$ C','interpreter','latex','FontSize',18)
ylabel('gm','interpreter','latex','FontSize',18)
set(gca,'FontSize', 18);
%tightfig;
print('mass','-dpng')
% set(gca,'FontSize','14')

figure
% subplot(2,3,5)
% plot(T_liq(index_h:index_c),vol_liq_n(index_h:index_c),'.r','LineWidth',5)
hold on
plot(T_liq(index_h:index_c),v_vol_liq,'-r','LineWidth',3)
hold on
% plot(T_liq(index_g_h:index_c),vol_gas_n(index_g_h:index_c),'.b','LineWidth',5)
hold on
plot(T_liq(index_h:index_c),v_vol_gas,'-b','LineWidth',3)
% plot(T_liq(index_h:index_c),vol_sol_n(index_h:index_c),'.k','LineWidth',5)
plot(T_liq(index_h:index_c),v_vol_sol,'-k','LineWidth',3)
title('Volume','interpreter','latex','FontSize',18)
legend('liquid','gas','solid')
legend('Location','northwest')
xlabel('Temperature in $^{\circ}$C','interpreter','latex','FontSize',18)
ylabel('Volume Fraction','interpreter','latex','FontSize',18)
%set(gca,'FontSize','14')
%set(gcf, 'Units', 'Inches', 'Position', [0, 0, 30, 20], 'PaperUnits', 'Inches', 'PaperSize', [30, 20])
set(gca,'FontSize', 18);
% tightfig;
print('vol','-dpng')

figure
plot(T_liq(index_h:index_c),H2O_liq(index_h:index_c),'.r')
hold on
plot(T_liq(index_h:index_c),v_H2O_liq)
hold on
plot(T_liq(index_g_h:index_c),H2O_gas(index_g_h:index_c),'.b')
hold on
plot(T_liq(index_h:index_c),v_H2O_gas)
title('H2O')
xlabel('Temperature in C')
ylabel('%')
print('h2o','-dpng')

figure
plot(T_liq(index_h:index_c),CO2_liq(index_h:index_c),'.r')
hold on
plot(T_liq(index_h:index_c),v_CO2_liq)
hold on
plot(T_liq(index_g_h:index_c),CO2_gas(index_g_h:index_c),'.b')
hold on
plot(T_liq(index_h:index_c),v_CO2_gas)
title('CO2')
xlabel('Temperature in C')
ylabel('%')
print('co2','-dpng')

figure
plot(T_liq(index_h:index_c),SiO2_liq(index_h:index_c),'.r')
hold on
plot(T_liq(index_h:index_c),v_SiO2_liq)
title('SiO2')
xlabel('Temperature in C')
ylabel('%')
print('sio2','-dpng')
% save('data.mat')

Colorss = [0.0317 0.6911 0.633;
    0.642 0.0344 0.9307;
    0.56687,0.488,0.97776;
    0.376409,0.97139,0.09359;
    0.212548,0.11245,0.66173;
    0.7921568,0.74321,0.6027728;
    0.14544,0.6385413,0.473817;
    0.4891424,0.59413,0.35625;
    0.83822,0.0762,0.0890;
    0.75, 0.75, 0;
    0.75, 0.75, 0.75;
    0, 0.75, 0.75;
    1, 0, 0; 
    1, 0, 1];
figure
hold on
plot(T_liq(index_h:index_c),mass_ol(index_h:index_c),'linewidth',10,'Color',Colorss(1,:))
plot(T_liq(index_h:index_c),mass_q(index_h:index_c),'linewidth',10,'Color',Colorss(2,:))
plot(T_liq(index_h:index_c),mass_f(index_h:index_c),'linewidth',10,'Color',Colorss(3,:))
plot(T_liq(index_h:index_c),mass_sp(index_h:index_c),'linewidth',10,'Color',Colorss(4,:))
plot(T_liq(index_h:index_c),mass_am(index_h:index_c),'linewidth',10,'Color',Colorss(5,:))
plot(T_liq(index_h:index_c),mass_ap(index_h:index_c),'linewidth',10,'Color',Colorss(6,:))
plot(T_liq(index_h:index_c),mass_bi(index_h:index_c),'linewidth',10,'Color',Colorss(7,:))
plot(T_liq(index_h:index_c),mass_cli(index_h:index_c),'linewidth',10,'Color',Colorss(8,:))
plot(T_liq(index_h:index_c),mass_or(index_h:index_c),'linewidth',10,'Color',Colorss(9,:))
plot(T_liq(index_h:index_c),mass_cm(index_h:index_c),'linewidth',10,'Color',Colorss(10,:))
plot(T_liq(index_h:index_c),mass_fa(index_h:index_c),'linewidth',10,'Color',Colorss(11,:))
plot(T_liq(index_h:index_c),mass_ga(index_h:index_c),'linewidth',10,'Color',Colorss(12,:))
plot(T_liq(index_h:index_c),mass_gr(index_h:index_c),'linewidth',10,'Color',Colorss(13,:))
plot(T_liq(index_h:index_c),mass_wh(index_h:index_c),'linewidth',10,'Color',Colorss(14,:))

legend('olivine','quartz','feldspar','spinel','amphibole','apatite','biotite','clinopyroxene','orthopyroxen',...
    'cummingtonite','fayalite','garnet','graphite','whitlockite')
title('Crystals Mass')
xlabel('Temperature in C')
ylabel('mass in cc')
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 15, 10], 'PaperUnits', 'Inches', 'PaperSize', [15, 10])
set(gca,'FontSize', 18);
print('crystal_mass','-dpng')

figure
hold on
plot(T_liq(index_h:index_c),rho_ol(index_h:index_c),'linewidth',10,'Color',Colorss(1,:))
plot(T_liq(index_h:index_c),rho_q(index_h:index_c),'linewidth',10,'Color',Colorss(2,:))
plot(T_liq(index_h:index_c),rho_f(index_h:index_c),'linewidth',10,'Color',Colorss(3,:))
plot(T_liq(index_h:index_c),rho_sp(index_h:index_c),'linewidth',10,'Color',Colorss(4,:))
plot(T_liq(index_h:index_c),rho_am(index_h:index_c),'linewidth',10,'Color',Colorss(5,:))
plot(T_liq(index_h:index_c),rho_ap(index_h:index_c),'linewidth',10,'Color',Colorss(6,:))
plot(T_liq(index_h:index_c),rho_bi(index_h:index_c),'linewidth',10,'Color',Colorss(7,:))
plot(T_liq(index_h:index_c),rho_cli(index_h:index_c),'linewidth',10,'Color',Colorss(8,:))
plot(T_liq(index_h:index_c),rho_or(index_h:index_c),'linewidth',10,'Color',Colorss(9,:))
plot(T_liq(index_h:index_c),rho_cm(index_h:index_c),'linewidth',10,'Color',Colorss(10,:))
plot(T_liq(index_h:index_c),rho_fa(index_h:index_c),'linewidth',10,'Color',Colorss(11,:))
plot(T_liq(index_h:index_c),rho_ga(index_h:index_c),'linewidth',10,'Color',Colorss(12,:))
plot(T_liq(index_h:index_c),rho_gr(index_h:index_c),'linewidth',10,'Color',Colorss(13,:))
plot(T_liq(index_h:index_c),rho_wh(index_h:index_c),'linewidth',10,'Color',Colorss(14,:))

legend('olivine','quartz','feldspar','spinel','amphibole','apatite','biotite','clinopyroxene','orthopyroxen',...
    'cummingtonite','fayalite','garnet','graphite','whitlockite')
title('Crystals density')
xlabel('Temperature in C')
ylabel('density in gm/cc')
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 15, 10], 'PaperUnits', 'Inches', 'PaperSize', [15, 10])
set(gca,'FontSize', 18);
print('crystal_rho','-dpng')

