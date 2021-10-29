%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Data input and analysis%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Details of each output found http://melts.ofm-research.org/Manual/UnixManHtml/Output.html
close all
clear all
clc
%constants
type = 0; %1=mafic 0=felsic
Tcold = 700;
Thot = 940;
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
minTsol = min(T_liq_sol);
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
        T_am = input.data(:,2); minTsol = min(T_am(T_am>652));iminTsol = find(round(T_am,2)==round(minTsol,2));
        index = zeros(size(T_am(1:iminTsol)));
        for i = 1:numel(index)
            index(i) = find(round(T_liq_sol,2)==round(T_am(i),2));
        end
        rho_am = zeros(size(T_liq_sol));
        mass_am = zeros(size(T_liq_sol));
        vol_am = zeros(size(T_liq_sol));
        rho_am(index) = input.data(1:iminTsol,6);
        mass_am(index) = input.data(1:iminTsol,5);
        vol_am(index) = input.data(1:iminTsol,29);
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
        T_ap = input.data(:,2);  minTsol = min(T_ap(T_ap>652));iminTsol = find(round(T_ap,2)==round(minTsol,2));
        index = zeros(size(T_ap(1:iminTsol)));
        for i = 1:numel(index)
            index(i) = find(round(T_liq_sol,2)==round(T_ap(i),2));
        end
        rho_ap = zeros(size(T_liq_sol));
        mass_ap = zeros(size(T_liq_sol));
        vol_ap = zeros(size(T_liq_sol));
        rho_ap(index) = input.data(1:iminTsol,6);
        mass_ap(index) = input.data(1:iminTsol,5);
        vol_ap(index) = input.data(1:iminTsol,29);
        %calculate the volume and mass and hence the density of the solids
        mass_sol_tbl = mass_sol_tbl+mass_ap;
        vol_sol_tbl = vol_sol_tbl+vol_ap;
    end
end
nm='cummingtonite.tbl';
if exist(nm,'file')==2
    str = {nm};
    deliniator = {','};
    filename=char(strcat(str(1)));
    deliniator=char(strcat(deliniator(1)));
    input = importdata(filename,deliniator);
    if isempty(input)==0
        T_cm = input.data(:,2); minTsol = min(T_cm(T_cm>652)); iminTsol = find(round(T_cm,2)==round(minTsol,2));
        index = zeros(size(T_cm(1:iminTsol)));
        for i = 1:numel(index)
            index(i) = find(round(T_liq_sol,2)==round(T_cm(i),2));
        end
        rho_cm = zeros(size(T_liq_sol));
        mass_cm = zeros(size(T_liq_sol));
        vol_cm = zeros(size(T_liq_sol));
        rho_cm(index) = input.data(1:iminTsol,6);
        mass_cm(index) = input.data(1:iminTsol,5);
        vol_cm(index) = input.data(1:iminTsol,29);
        %calculate the volume and mass and hence the density of the solids
        mass_sol_tbl = mass_sol_tbl+mass_cm;
        vol_sol_tbl = vol_sol_tbl+vol_cm;
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
        T_ap = input.data(:,2); minTsol = min(T_ap(T_ap>652)); iminTsol = find(round(T_ap,2)==round(minTsol,2));
        index = zeros(size(T_ap(1:iminTsol)));
        for i = 1:numel(index)
            index(i) = find(round(T_liq_sol,2)==round(T_ap(i),2));
        end
        rho_ap = zeros(size(T_liq_sol));
        mass_ap = zeros(size(T_liq_sol));
        vol_ap = zeros(size(T_liq_sol));
        rho_ap(index) = input.data(1:iminTsol,6);
        mass_ap(index) = input.data(1:iminTsol,5);
        vol_ap(index) = input.data(1:iminTsol,29);
        %calculate the volume and mass and hence the density of the solids
        mass_sol_tbl = mass_sol_tbl+mass_ap;
        vol_sol_tbl = vol_sol_tbl+vol_ap;
    end
end


nm='fayalite.tbl';
if exist(nm,'file')==2
    str = {nm};
    deliniator = {','};
    filename=char(strcat(str(1)));
    deliniator=char(strcat(deliniator(1)));
    input = importdata(filename,deliniator);
    if isempty(input)==0
       
        T_fa = input.data(:,2);minTsol = min(T_fa(T_fa>652));iminTsol = find(round(T_fa,2)==round(minTsol,2));
        index = zeros(size(T_fa(1:iminTsol)));
        for i = 1:numel(index)
            index(i) = find(round(T_liq_sol,2)==round(T_fa(i),2));
        end
        rho_fa = zeros(size(T_liq_sol));
        mass_fa = zeros(size(T_liq_sol));
        vol_fa = zeros(size(T_liq_sol));
        rho_fa(index) = input.data(1:iminTsol,6);
        mass_fa(index) = input.data(1:iminTsol,5);
        vol_fa(index) = input.data(1:iminTsol,29);
        %calculate the volume and mass and hence the density of the solids
        mass_sol_tbl = mass_sol_tbl+mass_fa;
        vol_sol_tbl = vol_sol_tbl+vol_fa;
    end
end
nm='garnet.tbl';
if exist(nm,'file')==2
    str = {nm};
    deliniator = {','};
    filename=char(strcat(str(1)));
    deliniator=char(strcat(deliniator(1)));
    input = importdata(filename,deliniator);
    if isempty(input)==0
        T_ga = input.data(:,2); minTsol = min(T_ga(T_ga>652)); iminTsol = find(round(T_ga,2)==round(minTsol,2));
        index = zeros(size(T_ga(1:iminTsol)));
        for i = 1:numel(index)
            index(i) = find(round(T_liq_sol,2)==round(T_ga(i),2));
        end
        rho_ga = zeros(size(T_liq_sol));
        mass_ga = zeros(size(T_liq_sol));
        vol_ga = zeros(size(T_liq_sol));
        rho_ga(index) = input.data(1:iminTsol,6);
        mass_ga(index) = input.data(1:iminTsol,5);
        vol_ga(index) = input.data(1:iminTsol,29);
        %calculate the volume and mass and hence the density of the solids
        mass_sol_tbl = mass_sol_tbl+mass_ga;
        vol_sol_tbl = vol_sol_tbl+vol_ga;
    end
end

nm='graphite.tbl';
if exist(nm,'file')==2
    str = {nm};
    deliniator = {','};
    filename=char(strcat(str(1)));
    deliniator=char(strcat(deliniator(1)));
    input = importdata(filename,deliniator);
    if isempty(input)==0
        T_gr = input.data(:,2); minTsol = min(T_gr(T_gr>652)); iminTsol = find(round(T_gr,2)==round(minTsol,2));
        index = zeros(size(T_gr(1:iminTsol)));
        for i = 1:numel(index)
            index(i) = find(round(T_liq_sol,2)==round(T_gr(i),2));
        end
        rho_gr = zeros(size(T_liq_sol));
        mass_gr = zeros(size(T_liq_sol));
        vol_gr = zeros(size(T_liq_sol));
        rho_gr(index) = input.data(1:iminTsol,6);
        mass_gr(index) = input.data(1:iminTsol,5);
        vol_gr(index) = input.data(1:iminTsol,29);
        %calculate the volume and mass and hence the density of the solids
        mass_sol_tbl = mass_sol_tbl+mass_gr;
        vol_sol_tbl = vol_sol_tbl+vol_gr;
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
        T_cli = input.data(:,2); minTsol = min(T_cli(T_cli>652)); iminTsol = find(round(T_cli,2)==round(minTsol,2));
        index = zeros(size(T_cli(1:iminTsol)));
        for i = 1:numel(index)
            index(i) = find(round(T_liq_sol,2)==round(T_cli(i),2));
        end
        %index=index(index~=0);
        rho_cli = zeros(size(T_liq_sol));
        mass_cli = zeros(size(T_liq_sol));
        vol_cli = zeros(size(T_liq_sol));
        rho_cli(index) = input.data(1:iminTsol,6);
        mass_cli(index) = input.data(1:iminTsol,5);
        vol_cli(index) = input.data(1:iminTsol,29);
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
        T_f = input.data(:,2); minTsol = min(T_f(T_f>652)); iminTsol = find(round(T_f,2)==round(minTsol,2));
        index = zeros(size(T_f(1:iminTsol)));
        for i = 1:numel(index)
            index(i) = find(round(T_liq_sol,2)==round(T_f(i),2));
        end
        rho_f = zeros(size(T_liq_sol));
        mass_f = zeros(size(T_liq_sol));
        vol_f = zeros(size(T_liq_sol));
        rho_f(index) = input.data(1:iminTsol,6);
        mass_f(index) = input.data(1:iminTsol,5);
        vol_f(index) = input.data(1:iminTsol,29);
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
        T_ol = input.data(:,2); minTsol = min(T_ol(T_ol>652)); iminTsol = find(round(T_ol,2)==round(minTsol,2));
        index = zeros(size(T_ol(1:iminTsol)));
        for i = 1:numel(index)
            index(i) = find(round(T_liq_sol,2)==round(T_ol(i),2));
        end
        rho_ol = zeros(size(T_liq_sol));
        mass_ol = zeros(size(T_liq_sol));
        vol_ol = zeros(size(T_liq_sol));
        rho_ol(index) = input.data(1:iminTsol,6);
        mass_ol(index) = input.data(1:iminTsol,5);
        vol_ol(index) = input.data(1:iminTsol,29);
        %calculate the volume and mass and hence the density of the solids
        mass_sol_tbl = mass_sol_tbl+mass_ol;
        vol_sol_tbl = vol_sol_tbl+vol_ol;
    end
end
nm='whitlockite.tbl';
if exist(nm,'file')==2
    str = {nm};
    deliniator = {','};
    filename=char(strcat(str(1)));
    deliniator=char(strcat(deliniator(1)));
    input = importdata(filename,deliniator);
    if isempty(input)==0
        T_wh = input.data(:,2); minTsol = min(T_wh(T_wh>652)); iminTsol = find(round(T_wh,2)==round(minTsol,2));
        index = zeros(size(T_wh(1:iminTsol)));
        for i = 1:numel(index)
            index(i) = find(round(T_liq_sol,2)==round(T_wh(i),2));
        end
        rho_wh = zeros(size(T_liq_sol));
        mass_wh = zeros(size(T_liq_sol));
        vol_wh = zeros(size(T_liq_sol));
        rho_wh(index) = input.data(1:iminTsol,6);
        mass_wh(index) = input.data(1:iminTsol,5);
        vol_wh(index) = input.data(1:iminTsol,29);
        %calculate the volume and mass and hence the density of the solids
        mass_sol_tbl = mass_sol_tbl+mass_wh;
        vol_sol_tbl = vol_sol_tbl+vol_wh;
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
        T_or = input.data(:,2); minTsol = min(T_or(T_or>652)); iminTsol = find(round(T_or,2)==round(minTsol,2));
        index = zeros(size(T_or(1:iminTsol)));
        for i = 1:numel(index)
            index(i) = find(round(T_liq_sol,2)==round(T_or(i),2));
        end
        rho_or = zeros(size(T_liq_sol));
        mass_or = zeros(size(T_liq_sol));
        vol_or = zeros(size(T_liq_sol));
        rho_or(index) = input.data(1:iminTsol,6);
        mass_or(index) = input.data(1:iminTsol,5);
        vol_or(index) = input.data(1:iminTsol,29);
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
        T_q = input.data(:,2); minTsol = min(T_q(T_q>652)); iminTsol = find(round(T_q,2)==round(minTsol,2));
        
        index = zeros(size(T_q(1:iminTsol)));
        for i = 1:numel(index)
            index(i) = find(round(T_liq_sol,2)==round(T_q(i),2));
        end
        rho_q = zeros(size(T_liq_sol));
        mass_q = zeros(size(T_liq_sol));
        vol_q = zeros(size(T_liq_sol));
        rho_q(index) = input.data(1:iminTsol,6);
        mass_q(index) = input.data(1:iminTsol,5);
        vol_q(index) = input.data(1:iminTsol,29);
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
        T_sp = input.data(:,2); minTsol = min(T_sp(T_sp>652)); iminTsol = find(round(T_sp,2)==round(minTsol,2));
        
        index = zeros(size(T_sp(1:iminTsol)));
        for i = 1:numel(index)
            index(i) = find(round(T_liq_sol,2)==round(T_sp(i),2));
        end
        rho_sp = zeros(size(T_liq_sol));
        mass_sp = zeros(size(T_liq_sol));
        vol_sp = zeros(size(T_liq_sol));
        rho_sp(index) = input.data(1:iminTsol,6);
        mass_sp(index) = input.data(1:iminTsol,5);
        vol_sp(index) = input.data(1:iminTsol,29);
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
if ~exist('mass_ol','var')
mass_ol = zeros(size(T_liq));
rho_ol = NaN(size(T_liq));
end
hold on
if ~exist('mass_f','var')
    mass_f = zeros(size(T_liq));
    rho_f = NaN(size(T_liq));
end
if ~exist('mass_sp','var')
    mass_sp = zeros(size(T_liq));
    rho_sp = NaN(size(T_liq));
end
if ~exist('mass_q','var')
    mass_q = zeros(size(T_liq));
    rho_q = NaN(size(T_liq));
end
if ~exist('mass_am','var')
    mass_am = zeros(size(T_liq));
    rho_am = NaN(size(T_liq));
end
if ~exist('mass_cm','var')
    mass_cm = zeros(size(T_liq));
    rho_cm = NaN(size(T_liq));
end
if ~exist('mass_fa','var')
    mass_fa = zeros(size(T_liq));
    rho_fa = NaN(size(T_liq));
end
if ~exist('mass_ga','var')
    mass_ga = zeros(size(T_liq));
    rho_ga = NaN(size(T_liq));
end
if ~exist('mass_gr','var')
    mass_gr = zeros(size(T_liq));
    rho_gr = NaN(size(T_liq));
end
if ~exist('mass_wh','var')
    mass_wh = zeros(size(T_liq));
    rho_wh = NaN(size(T_liq));
end
if ~exist('mass_ap','var')
    mass_ap = zeros(size(T_liq));
    rho_ap = NaN(size(T_liq));
end
if ~exist('mass_bi','var')
    mass_bi = zeros(size(T_liq));
    rho_bi = NaN(size(T_liq));
end
if ~exist('mass_cli','var')
    mass_cli = zeros(size(T_liq));
    rho_cli = NaN(size(T_liq));
end
if ~exist('mass_or','var')
    mass_or = zeros(size(T_liq));
    rho_or = NaN(size(T_liq));
end
rho_sol_tbl = mass_sol_tbl./vol_sol_tbl;
% rho_sol_tbl (isnan(rho_sol_tbl))=0;

%%
%%%Index
[c_c, index_c] = min(abs(T_liq-Tcold));
[c_h, index_h] = min(abs(T_liq-Thot));
[c_c_g, index_c_g] = min(abs(T_liq_sol-Tcold));
[c_h_g, index_h_g] = min(abs(T_liq_sol-Thot));
[c_c_c, index_c_c] = min(abs(T_liq_sol-Tcold));
[c_h_c, index_h_c] = min(abs(T_liq_sol-Thot));
begingas = find(rho_gas>0.03); c_g_h = max(T_liq(begingas));index_g_h = find(T_liq==c_g_h);
% index_g_h = find(T_liq==c_g_h);
rho_C = rho_liq(index_c);
rho_H = rho_liq(index_h);
mu_C = mu_liq(index_c);
mu_H = mu_liq(index_h);

%fitting functions
p_rho_liq = polyfit(T_liq(index_h:index_c),rho_liq(index_h:index_c),2);
p_rho_gas = polyfit(T_liq(index_g_h:index_c),rho_gas(index_g_h:index_c),10);
p_rho_sol = polyfit(T_liq_sol(index_h_c:index_c_c),rho_sol_tbl(index_h_c:index_c_c),10);
p_mu_liq = polyfit(T_liq(index_h:index_c),mu_liq(index_h:index_c),1);
p_mass_liq = polyfit(T_liq(index_h:index_c),mass_liq(index_h:index_c),10);
p_mass_gas = polyfit(T_liq(index_g_h:index_c),mass_gas(index_g_h:index_c),10);
p_mass_sol = polyfit(T_liq_sol(index_h_c:index_c_c),mass_sol_tbl(index_h_c:index_c_c),10);
p_H2O_liq = polyfit(T_liq(index_h:index_c),H2O_liq(index_h:index_c),10);
p_H2O_gas = polyfit(T_liq(index_g_h:index_c),H2O_gas(index_g_h:index_c),10);
p_CO2_liq = polyfit(T_liq(index_h:index_c),CO2_liq(index_h:index_c),10);
p_CO2_gas = polyfit(T_liq(index_g_h:index_c),CO2_gas(index_g_h:index_c),10);
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
vol_sys = vol_liq + vol_sol+vol_gass; %tbl
vol_liq_n = vol_liq./vol_sys;
vol_sol_n = vol_sol./vol_sys; %tbl
vol_gas_n = vol_gass./vol_sys;
% for i = 1:numel(vol_gas);
%     [c, index] = min(abs(T_liq-T_liq(i)));
%     vol_gas_n(i) = vol_gas(i)./vol_sys(index);
% end
%%%%Removing Nan values
rho_sol_n1 = rho_sol_tbl(index_h:index_c); 
trash = find(isnan(rho_sol_n1)|isinf(rho_sol_n1));L = isnan(rho_sol_n1) + isinf(rho_sol_n1);
nonzero = rho_sol_n1(~L);
if size(trash)>0
    for i = 1:numel(trash)
        if trash(i) ==1
            rho_sol_n1(trash(i)) = nonzero(1);
        else
        rho_sol_n1(trash(i)) = rho_sol_n1(trash(i)-1);
        end
    end
end
rho_sol_tbl(index_h:index_c) = rho_sol_n1;
vol_liq_n1 = vol_liq_n(index_h:index_c);
trash = find(isnan(vol_liq_n1));
if size(trash)>0
    if trash == numel(vol_liq_n1)
        vol_liq_n1(trash) = vol_liq_n1(trash-1);
    else
        vol_liq_n1(trash) = mean([vol_liq_n1(trash+1) vol_liq_n1(trash-1)]);
    end
end
vol_gas_n1 = vol_gas_n(index_g_h:index_c);
trash = find(isnan(vol_gas_n1));
if size(trash)>0
    if trash == numel(vol_gas_n1)
        vol_gas_n1(trash) = vol_gas_n1(trash-1);
    else
        vol_gas_n1(trash) = mean([vol_gas_n1(trash+1) vol_gas_n1(trash-1)]);
    end
end
vol_sol_tbl1 = vol_sol_n(index_h:index_c);
trash = find(isnan(vol_sol_tbl1));
if size(trash)>0
    if trash == numel(vol_sol_tbl1)
        vol_sol_tbl1(trash) = vol_sol_tbl1(trash-1);
    else
        vol_sol_tbl1(trash) = mean([vol_sol_tbl1(trash+1) vol_sol_tbl1(trash-1)]);
    end
end
vol_sol_tbl = [zeros(index_h-1,1); vol_sol_tbl1; vol_sol_tbl(end).*ones(numel(T_liq)-index_c,1)];

p_vol_liq = polyfit(T_liq(index_h:index_c),vol_liq_n1,3);
p_vol_gas = polyfit(T_liq(index_g_h:index_c),vol_gas_n1,3);
p_vol_sol = polyfit(T_liq(index_h:index_c),vol_sol_tbl(index_h:index_c),3);
p_rho_sol = polyfit(T_liq_sol(index_h:index_c),rho_sol_tbl(index_h:index_c),10);
v_vol_liq = polyval(p_vol_liq,T_liq(index_h:index_c));
v_vol_gas = polyval(p_vol_gas,T_liq(index_h:index_c));
% v_vol_gas = [zeros(index_g_h,1); v_vol_gas];
v_vol_gas(v_vol_gas<0)=0;
v_vol_sol = polyval(p_vol_sol,T_liq(index_h:index_c));
v_rho_sol=polyval(p_rho_sol,T_liq(index_h:index_c));


%% Calculating the bulk density change due to
%thermal expansivity
dpdT = dVdT_liq.*(rho_liq(2)-rho_liq(1))./(vol_liq(2)-vol_liq(1));
rho_thermal = dpdT.*(T_liq-1100)+rho_liq(1);
drho_thermal = (rho_thermal(index_h)-rho_thermal(index_c))/(Thot-Tcold);

%chemical alteration
rho_chemical = rho_liq - dpdT.*T_liq;
drho_chemical = (rho_chemical(index_h)-rho_chemical(index_c))/(Thot-Tcold);

%bubble presence 
v_rho_gas = [zeros(numel(v_vol_gas)-numel(v_rho_gas),1);v_rho_gas];
% v_vol_gas(v_vol_gas<0)=0;

rho_bub = v_vol_gas.*v_rho_gas+ (1-v_vol_gas)*rho_liq(1);
drho_bub = (rho_bub(1)-rho_bub(end))/(Thot-Tcold);

%crystal presence
rho_cry = v_vol_sol.*v_rho_sol + (1-v_vol_sol).*rho_liq(1);
drho_cry = (rho_cry(1)-rho_cry(end))/(Thot-Tcold);

%bulk density
RHO = rho_liq(index_h:index_c).*(1-v_vol_sol-v_vol_gas)+v_rho_gas.*v_vol_gas+v_vol_sol.*v_rho_sol;%+dpdT.*T_liq;
dRHO = (RHO(1)-RHO(end))/(Thot-Tcold);
 save('datatest.mat')


