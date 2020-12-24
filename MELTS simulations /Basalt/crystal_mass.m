%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Data input and analysis%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Details of each output found http://melts.ofm-research.org/Manual/UnixManHtml/Output.html
close all
clear all
clc
%constants
type = 1; %1=mafic 0=felsic
Tcold = 800;
Thot = 1040;
Pressure = 200; %MPa
gas = 1; %1 if there is gas present
%Data input general analysis
str = {'melts-liquid.tbl'};
deliniator = {','};
filename=char(strcat(str(1)));
deliniator=char(strcat(deliniator(1)));
input = importdata(filename,deliniator);
T_liq = input.data(:,2);
mu_liq = (input.data(:,50)-1); %in log10 poise to log10 Pa s
mass_liq = input.data(:,5);
mass_sol = input.data(:,51);
rho_liq = input.data(:,6);
rho_sol = input.data(:,52);
H2O_liq = input.data(:,21);
CO2_liq = input.data(:,22);
SiO2_liq = input.data(:,7);
visc_liq = .1*10.^input.data(:,50);
vol_sys = input.data(:,61);
vol_liq = input.data(:,29);
vol_sol = input.data(:,56);
dVdT_liq = input.data(:,67);
dVdT_sys = input.data(:,63);

nm = 'fluid.tbl';
if exist(nm,'file')==2
    str = {nm};
    deliniator = {','};
    filename=char(strcat(str(1)));
    deliniator=char(strcat(deliniator(1)));
    input = importdata(filename,deliniator);
end
T_gas = input.data(:,2).*gas;
rho_gas = [zeros(numel(vol_liq)-numel(input.data(:,6).*gas),1); input.data(:,6).*gas];
mass_gas = [zeros(numel(vol_liq)-numel(input.data(:,5).*gas),1);input.data(:,5).*gas];
H2O_gas = [zeros(numel(vol_liq)-numel(input.data(:,21).*gas),1); input.data(:,21).*gas];
CO2_gas = [zeros(numel(vol_liq)-numel(input.data(:,22).*gas),1); input.data(:,22).*gas];
vol_gas = [zeros(numel(vol_liq)-numel(input.data(:,29).*gas),1); input.data(:,29).*gas];

%% dealing with crystal density

T_liq_sol = [max(T_liq):(T_liq(2)-T_liq(1)):min(T_liq)]';
vol_sol_tbl = zeros(size(T_liq_sol));
mass_sol_tbl = zeros(size(T_liq_sol));
nm='amphibole.tbl';
if exist(nm,'file')==2
    str = {nm};
    deliniator = {','};
    filename=char(strcat(str(1)));
    deliniator=char(strcat(deliniator(1)));
    input = importdata(filename,deliniator);
    if isempty(input)==0
        T_am = input.data(:,2);
        index = zeros(size(T_am));
        for i = 1:numel(T_am)
            index(i) = find(T_liq_sol==T_am(i));
        end
        rho_am = zeros(size(T_liq_sol));
        mass_am = zeros(size(T_liq_sol));
        vol_am = zeros(size(T_liq_sol));
        rho_am(index) = input.data(:,6);
        mass_am(index) = input.data(:,5);
        vol_am(index) = input.data(:,29);
        %calculate the volume and mass and hence the density of the solids
        mass_sol_tbl = mass_sol_tbl+mass_am;
        vol_sol_tbl = vol_sol_tbl+vol_am;
    end
end
nm='apatite.tbl';
if exist(nm,'file')==2
    str = {nm};
    deliniator = {','};
    filename=char(strcat(str(1)));
    deliniator=char(strcat(deliniator(1)));
    input = importdata(filename,deliniator);
    if isempty(input)==0
        T_ap = input.data(:,2);
        index = zeros(size(T_ap));
        for i = 1:numel(T_ap)
            index(i) = find(T_liq_sol==T_ap(i));
        end
        rho_ap = zeros(size(T_liq_sol));
        mass_ap = zeros(size(T_liq_sol));
        vol_ap = zeros(size(T_liq_sol));
        rho_ap(index) = input.data(:,6);
        mass_ap(index) = input.data(:,5);
        vol_ap(index) = input.data(:,29);
        %calculate the volume and mass and hence the density of the solids
        mass_sol_tbl = mass_sol_tbl+mass_ap;
        vol_sol_tbl = vol_sol_tbl+vol_ap;
    end
end

nm='biotite.tbl';
if exist(nm,'file')==2
    str = {nm};
    deliniator = {','};
    filename=char(strcat(str(1)));
    deliniator=char(strcat(deliniator(1)));
    input = importdata(filename,deliniator);
    if isempty(input)==0
        T_bi = input.data(:,2);
        index = zeros(size(T_bi));
        for i = 1:numel(T_bi)
            index(i) = find(T_liq_sol==T_bi(i));
        end
        rho_bi = zeros(size(T_liq_sol));
        mass_bi = zeros(size(T_liq_sol));
        vol_bi = zeros(size(T_liq_sol));
        rho_bi(index) = input.data(:,6);
        mass_bi(index) = input.data(:,5);
        vol_bi(index) = input.data(:,29);
        %calculate the volume and mass and hence the density of the solids
        mass_sol_tbl = mass_sol_tbl+mass_bi;
        vol_sol_tbl = vol_sol_tbl+vol_bi;
    end
end

nm='clinopyroxene.tbl';
if exist(nm,'file')==2
    str = {nm};
    deliniator = {','};
    filename=char(strcat(str(1)));
    deliniator=char(strcat(deliniator(1)));
    input = importdata(filename,deliniator);
    if isempty(input)==0
        T_cli = input.data(:,2);
        index = zeros(size(T_cli));
        for i = 1:numel(T_cli)
            if isempty(find(T_liq_sol==T_cli(i)))==0
            index(i) = find(T_liq_sol==T_cli(i));
            end
        end
        %index=index(index~=0);
        rho_cli = zeros(size(T_liq_sol));
        mass_cli = zeros(size(T_liq_sol));
        vol_cli = zeros(size(T_liq_sol));
        rho_cli(index) = input.data(:,6);
        mass_cli(index) = input.data(:,5);
        vol_cli(index) = input.data(:,29);
        %calculate the volume and mass and hence the density of the solids
        mass_sol_tbl = mass_sol_tbl+mass_cli;
        vol_sol_tbl = vol_sol_tbl+vol_cli;
    end
end
nm = 'feldspar.tbl';
if exist(nm,'file')==2
    str = {nm};
    deliniator = {','};
    filename=char(strcat(str(1)));
    deliniator=char(strcat(deliniator(1)));
    input = importdata(filename,deliniator);
    if isempty(input)==0
        T_f = input.data(:,2);
        index = zeros(size(T_f));
        for i = 1:numel(T_f)
            index(i) = find(T_liq_sol==T_f(i));
        end
        rho_f = zeros(size(T_liq_sol));
        mass_f = zeros(size(T_liq_sol));
        vol_f = zeros(size(T_liq_sol));
        rho_f(index) = input.data(:,6);
        mass_f(index) = input.data(:,5);
        vol_f(index) = input.data(:,29);
        %calculate the volume and mass and hence the density of the solids
        mass_sol_tbl = mass_sol_tbl+mass_f;
        vol_sol_tbl = vol_sol_tbl+vol_f;
    end
end
nm = 'olivine.tbl';
if exist(nm,'file')==2
    str = {nm};
    deliniator = {','};
    filename=char(strcat(str(1)));
    deliniator=char(strcat(deliniator(1)));
    input = importdata(filename,deliniator);
    if isempty(input)==0
        T_ol = input.data(:,2);
        index = zeros(size(T_ol));
        for i = 1:numel(T_ol)
            index(i) = find(T_liq_sol==T_ol(i));
        end
        rho_ol = zeros(size(T_liq_sol));
        mass_ol = zeros(size(T_liq_sol));
        vol_ol = zeros(size(T_liq_sol));
        rho_ol(index) = input.data(:,6);
        mass_ol(index) = input.data(:,5);
        vol_ol(index) = input.data(:,29);
        %calculate the volume and mass and hence the density of the solids
        mass_sol_tbl = mass_sol_tbl+mass_ol;
        vol_sol_tbl = vol_sol_tbl+vol_ol;
    end
end
nm='orthopyroxene.tbl';
if exist(nm,'file')==2
    str = {nm};
    deliniator = {','};
    filename=char(strcat(str(1)));
    deliniator=char(strcat(deliniator(1)));
    input = importdata(filename,deliniator);
    if isempty(input)==0
        T_or = input.data(:,2);
        index = zeros(size(T_or));
        for i = 1:numel(T_or)
            index(i) = find(T_liq_sol==T_or(i));
        end
        rho_or = zeros(size(T_liq_sol));
        mass_or = zeros(size(T_liq_sol));
        vol_or = zeros(size(T_liq_sol));
        rho_or(index) = input.data(:,6);
        mass_or(index) = input.data(:,5);
        vol_or(index) = input.data(:,29);
        %calculate the volume and mass and hence the density of the solids
        mass_sol_tbl = mass_sol_tbl+mass_or;
        vol_sol_tbl = vol_sol_tbl+vol_or;
    end
end
nm='quartz.tbl';
if exist(nm,'file')==2
    str = {nm};
    deliniator = {','};
    filename=char(strcat(str(1)));
    deliniator=char(strcat(deliniator(1)));
    input = importdata(filename,deliniator);
    if isempty(input)==0
        T_q = input.data(:,2);
        index = zeros(size(T_q));
        for i = 1:numel(T_q)
            index(i) = find(T_liq_sol==T_q(i));
        end
        rho_q = zeros(size(T_liq_sol));
        mass_q = zeros(size(T_liq_sol));
        vol_q = zeros(size(T_liq_sol));
        rho_q(index) = input.data(:,6);
        mass_q(index) = input.data(:,5);
        vol_q(index) = input.data(:,29);
        %calculate the volume and mass and hence the density of the solids
        mass_sol_tbl = mass_sol_tbl+mass_q;
        vol_sol_tbl = vol_sol_tbl+vol_q;
    end
end
nm='spinel.tbl';
if exist(nm,'file')==2
    str = {nm};
    deliniator = {','};
    filename=char(strcat(str(1)));
    deliniator=char(strcat(deliniator(1)));
    input = importdata(filename,deliniator);
    if isempty(input)==0
        T_sp = input.data(:,2);
        index = zeros(size(T_sp));
        for i = 1:numel(T_sp)
            index(i) = find(T_liq_sol==T_sp(i));
        end
        rho_sp = zeros(size(T_liq_sol));
        mass_sp = zeros(size(T_liq_sol));
        vol_sp = zeros(size(T_liq_sol));
        rho_sp(index) = input.data(:,6);
        mass_sp(index) = input.data(:,5);
        vol_sp(index) = input.data(:,29);
        %calculate the volume and mass and hence the density of the solids
        mass_sol_tbl = mass_sol_tbl+mass_sp;
        vol_sol_tbl = vol_sol_tbl+vol_sp;
    end
end
nm = 'rhm-oxide.tbl';
if exist(nm,'file')==2
    str = {nm};
    deliniator = {','};
    filename=char(strcat(str(1)));
    deliniator=char(strcat(deliniator(1)));
    input = importdata(filename,deliniator);
end

rho_sol_tbl = mass_sol_tbl./vol_sol_tbl;
rho_sol_tbl (isnan(rho_sol_tbl))=0;
%%
%%%Index
[c_c, index_c] = min(abs(T_liq-Tcold));
[c_h, index_h] = min(abs(T_liq-Thot));
[c_c_g, index_c_g] = min(abs(T_liq_sol-Tcold));
[c_h_g, index_h_g] = min(abs(T_liq_sol-Thot));
[c_c_c, index_c_c] = min(abs(T_liq_sol-Tcold));
[c_h_c, index_h_c] = min(abs(T_liq_sol-Thot));
rho_C = rho_liq(index_c);
rho_H = rho_liq(index_h);
mu_C = mu_liq(index_c);
mu_H = mu_liq(index_h);

%fitting functions
p_rho_liq = polyfit(T_liq(index_h:index_c),rho_liq(index_h:index_c),10);
p_rho_gas = polyfit(T_liq(index_h:index_c),rho_gas(index_h:index_c),10);
p_rho_sol = polyfit(T_liq_sol(index_h_c:index_c_c),rho_sol_tbl(index_h_c:index_c_c),10);
p_mu_liq = polyfit(T_liq(index_h:index_c),mu_liq(index_h:index_c),10);
p_mass_liq = polyfit(T_liq(index_h:index_c),mass_liq(index_h:index_c),10);
p_mass_gas = polyfit(T_liq(index_h:index_c),mass_gas(index_h:index_c),10);
p_mass_sol = polyfit(T_liq_sol(index_h_c:index_c_c),mass_sol_tbl(index_h_c:index_c_c),10);
p_H2O_liq = polyfit(T_liq(index_h:index_c),H2O_liq(index_h:index_c),10);
p_H2O_gas = polyfit(T_liq(index_h:index_c),H2O_gas(index_h:index_c),10);
p_CO2_liq = polyfit(T_liq(index_h:index_c),CO2_liq(index_h:index_c),10);
p_CO2_gas = polyfit(T_liq(index_h:index_c),CO2_gas(index_h:index_c),10);
p_SiO2_liq = polyfit(T_liq(index_h:index_c),SiO2_liq(index_h:index_c),10);


%extrapolating points
v_rho_liq=polyval(p_rho_liq,T_liq(index_h:index_c));
v_rho_gas=polyval(p_rho_gas,T_liq(index_h:index_c));
v_rho_sol=polyval(p_rho_sol,T_liq(index_h:index_c));
v_mu_liq=polyval(p_mu_liq,T_liq(index_h:index_c));
v_mass_liq = polyval(p_mass_liq,T_liq(index_h:index_c));
v_mass_gas = polyval(p_mass_gas,T_liq(index_h:index_c));
v_mass_sol = polyval(p_mass_sol,T_liq(index_h:index_c));
v_H2O_liq = polyval(p_H2O_liq,T_liq(index_h:index_c));
v_H2O_gas = polyval(p_H2O_gas,T_liq(index_h:index_c));
v_CO2_liq = polyval(p_CO2_liq,T_liq(index_h:index_c));
v_CO2_gas = polyval(p_CO2_gas,T_liq(index_h:index_c));
v_SiO2_liq = polyval(p_SiO2_liq,T_liq(index_h:index_c));

%density change coefficients
alpha = (polyval(p_rho_liq,Tcold) - polyval(p_rho_liq,Thot))/polyval(p_rho_liq,Thot);
beta = (polyval(p_rho_gas,(Thot+Tcold)/2) - polyval(p_rho_liq,Thot))/polyval(p_rho_liq,Thot);


%volume percent


vol_gass = [zeros(numel(vol_liq)-numel(vol_gas),1); vol_gas];
vol_sys = vol_liq + vol_sol_tbl+vol_gass;
vol_liq_n = vol_liq./vol_sys;
vol_sol_n = vol_sol_tbl./vol_sys;
vol_gas_n = vol_gass./vol_sys;
% for i = 1:numel(vol_gas);
%     [c, index] = min(abs(T_liq-T_liq(i)));
%     vol_gas_n(i) = vol_gas(i)./vol_sys(index);
% end
%%%%Removing Nan values
vol_liq_n1 = vol_liq_n(index_h:index_c);
trash = find(isnan(vol_liq_n1));
vol_liq_n1(trash) = mean([vol_liq_n1(trash+1) vol_liq_n1(trash-1)]);

vol_gas_n1 = vol_gas_n(index_h:index_c);
trash = find(isnan(vol_gas_n1));
vol_gas_n1(trash) = mean([vol_gas_n1(trash+1) vol_gas_n1(trash-1)]);

vol_sol_tbl1 = vol_sol_n(index_h:index_c);
trash = find(isnan(vol_sol_tbl1));
vol_sol_tbl1(trash) = mean([vol_sol_tbl1(trash+1) vol_sol_tbl1(trash-1)]);


p_vol_liq = polyfit(T_liq(index_h:index_c),vol_liq_n1,10);
p_vol_gas = polyfit(T_liq(index_h:index_c),vol_gas_n1,10);
p_vol_sol = polyfit(T_liq(index_h:index_c),vol_sol_tbl1,3);
v_vol_liq = polyval(p_vol_liq,T_liq(index_h:index_c));
v_vol_gas = polyval(p_vol_gas,T_liq(index_h:index_c));
v_vol_sol = polyval(p_vol_sol,T_liq(index_h:index_c));


%% Calculating the bulk density change due to
%thermal expansivity
dpdT = dVdT_liq.*(rho_liq(2)-rho_liq(1))./(vol_liq(2)-vol_liq(1));
rho_thermal = dpdT.*(T_liq-1100)+rho_liq(1);
drho_thermal = (rho_thermal(index_h)-rho_thermal(index_c))/(Thot-Tcold);

%chemical alteration
rho_chemical = rho_liq - dpdT.*T_liq;
drho_chemical = (rho_chemical(index_h)-rho_chemical(index_c))/(Thot-Tcold);

%bubble presence
rho_bub = vol_gas_n.*rho_gas+ (1-vol_gas_n)*rho_liq(1);
drho_bub = (rho_bub(index_h)-rho_bub(index_c))/(Thot-Tcold);

%crystal presence
rho_cry = vol_sol_n.*rho_sol_tbl + (1-vol_sol_n).*rho_liq(1);
drho_cry = (rho_cry(index_h)-rho_cry(index_c))/(Thot-Tcold);

%bulk density
RHO = rho_liq.*(1-vol_sol_n-vol_gas_n)+rho_gas.*vol_gas_n+vol_sol_n.*rho_sol_tbl;%+dpdT.*T_liq;
dRHO = (RHO(index_h)-RHO(index_c))/(Thot-Tcold);

%% plotting
figure;
subplot(2,3,1)
hold on;
plot(T_liq(index_h:index_c),rho_thermal(index_h:index_c),'m','linewidth',10)
plot(T_liq(index_h:index_c),rho_chemical(index_h:index_c),'r','linewidth',10)
plot(T_liq(index_h:index_c),rho_bub(index_h:index_c),'b','linewidth',10)
plot(T_liq(index_h:index_c),rho_cry(index_h:index_c),'g','linewidth',10)
plot(T_liq(index_h:index_c),RHO(index_h:index_c),'.k','linewidth',10)
title('Bulk \rho due to different components')
xlabel('Temperature in C')
ylabel('in gm/cc or x1000 kg/m^3')
th = sprintf('thermal d\rho/dT = %.1E', drho_thermal);
ch = sprintf('chemical d\rho/dT = %.1E', drho_chemical);
bb = sprintf('bubble d\rho/dT = %.1E', drho_bub);
cry = sprintf('crystal d\rho/dT = %.1E', drho_cry);
r = sprintf('\rho d\rho/dT = %.1E', dRHO);
legend(th,ch,bb,cry,r)
legend('Location','southeast')
set(gca,'FontSize', 18);
%tightfig;
%set(gca,'FontSize','14')

%print('bulk_rho','-dpdf')



%plotting
%figure
subplot(2,3,2)
plot(T_liq(index_h:index_c),rho_liq(index_h:index_c),'.r','linewidth',10)
hold on
plot(T_liq(index_h:index_c),v_rho_liq)
hold on
plot(T_liq(index_h:index_c),rho_gas(index_h:index_c),'.b','linewidth',10)
plot(T_liq(index_h:index_c),rho_sol_tbl(index_h:index_c),'.k','linewidth',10)
hold on
plot(T_liq(index_h:index_c),v_rho_gas)
plot(T_liq(index_h:index_c),v_rho_sol)
title('\rho')
xlabel('Temperature in C')
ylabel('in gm/cc or x1000 kg/m^3')
set(gca,'FontSize', 18);
%tightfig;
%set(gca,'FontSize','14')
%print('rho','-dpdf')

%figure
subplot(2,3,[3 6])
plot(T_liq(index_h:index_c),mu_liq(index_h:index_c),'.r')
hold on
plot(T_liq(index_h:index_c),v_mu_liq)
hold on
title('\mu')
xlabel('Temperature in C')
ylabel('10^ Pa s')
set(gca,'FontSize', 18);
%tightfig;
%set(gca,'FontSize','14')
%print('mu','-dpdf')

%figure
subplot(2,3,4)
plot(T_liq(index_h:index_c),mass_liq(index_h:index_c),'.r')
hold on
plot(T_liq(index_h:index_c),v_mass_liq)
hold on
plot(T_liq(index_h:index_c),mass_gas(index_h:index_c),'.b')
hold on
plot(T_liq(index_h:index_c),v_mass_gas)
plot(T_liq(index_h:index_c),mass_sol_tbl(index_h:index_c),'.k')
plot(T_liq(index_h:index_c),v_mass_sol)
title('Mass')
legend('liquid','','gas','','solid','')
xlabel('Temperature in C')
ylabel('gm')
set(gca,'FontSize', 18);
%tightfig;
%print('mass','-dpdf')
%set(gca,'FontSize','14')

%figure
subplot(2,3,5)
plot(T_liq(index_h:index_c),vol_liq_n(index_h:index_c),'.r','LineWidth',5)
hold on
plot(T_liq(index_h:index_c),v_vol_liq)
hold on
plot(T_liq(index_h:index_c),vol_gas_n(index_h:index_c),'.b','LineWidth',5)
hold on
plot(T_liq(index_h:index_c),v_vol_gas)
plot(T_liq(index_h:index_c),vol_sol_n(index_h:index_c),'.k','LineWidth',5)
plot(T_liq(index_h:index_c),v_vol_sol)
title('Volume')
%legend('liquid','','gas','','solid','')
xlabel('Temperature in C')
ylabel('Volume Fraction')
%set(gca,'FontSize','14')
%set(gcf, 'Units', 'Inches', 'Position', [0, 0, 30, 20], 'PaperUnits', 'Inches', 'PaperSize', [30, 20])
set(gca,'FontSize', 18);
tightfig;
print('vol','-dpng')

figure
plot(T_liq(index_h:index_c),H2O_liq(index_h:index_c),'.r')
hold on
plot(T_liq(index_h:index_c),v_H2O_liq)
hold on
plot(T_liq(index_h:index_c),H2O_gas(index_h:index_c),'.b')
hold on
plot(T_liq(index_h:index_c),v_H2O_gas)
title('H2O')
xlabel('Temperature in C')
ylabel('%')
print('h2o','-dpdf')

figure
plot(T_liq(index_h:index_c),CO2_liq(index_h:index_c),'.r')
hold on
plot(T_liq(index_h:index_c),v_CO2_liq)
hold on
plot(T_liq(index_h:index_c),CO2_gas(index_h:index_c),'.b')
hold on
plot(T_liq(index_h:index_c),v_CO2_gas)
title('CO2')
xlabel('Temperature in C')
ylabel('%')
print('co2','-dpdf')

figure
plot(T_liq(index_h:index_c),SiO2_liq(index_h:index_c),'.r')
hold on
plot(T_liq(index_h:index_c),v_SiO2_liq)
title('SiO2')
xlabel('Temperature in C')
ylabel('%')
print('sio2','-dpdf')
save('data.mat')

figure
plot(T_liq(index_h:index_c),mass_ol(index_h:index_c),'linewidth',10)
hold on
plot(T_liq(index_h:index_c),mass_f(index_h:index_c),'linewidth',10)
plot(T_liq(index_h:index_c),mass_am(index_h:index_c),'linewidth',10)
plot(T_liq(index_h:index_c),mass_ap(index_h:index_c),'linewidth',10)
plot(T_liq(index_h:index_c),mass_bi(index_h:index_c),'linewidth',10)
plot(T_liq(index_h:index_c),mass_cli(index_h:index_c),'linewidth',10)
plot(T_liq(index_h:index_c),mass_or(index_h:index_c),'linewidth',10)
legend('olivine','feldspar','amphibole','apatite','biotite','clinopyroxene','orthopyroxen')
title('Crystals Mass')
xlabel('Temperature in C')
ylabel('mass in cc')
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 15, 10], 'PaperUnits', 'Inches', 'PaperSize', [15, 10])
set(gca,'FontSize', 35);
print('crystal_mass','-dpng')

figure
plot(T_liq(index_h:index_c),rho_ol(index_h:index_c),'linewidth',10)
hold on
plot(T_liq(index_h:index_c),rho_f(index_h:index_c),'linewidth',10)
plot(T_liq(index_h:index_c),rho_am(index_h:index_c),'linewidth',10)
plot(T_liq(index_h:index_c),rho_ap(index_h:index_c),'linewidth',10)
plot(T_liq(index_h:index_c),rho_bi(index_h:index_c),'linewidth',10)
plot(T_liq(index_h:index_c),rho_cli(index_h:index_c),'linewidth',10)
plot(T_liq(index_h:index_c),rho_or(index_h:index_c),'linewidth',10)
%legend('olivine','feldspar','amphibole','apatite','biotite','clinopyroxene','orthopyroxen')
title('Crystals density')
xlabel('Temperature in C')
ylabel('density in gm/cc')
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 15, 10], 'PaperUnits', 'Inches', 'PaperSize', [15, 10])
set(gca,'FontSize', 35);
print('crystal_rho','-dpng')



