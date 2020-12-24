%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Felspar equilibrium%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Details of each output found http://melts.ofm-research.org/Manual/UnixManHtml/Output.html
close all
clear all

type = 0; %1=mafic 0=felsic
Tcold = 702;
Thot = 940;
Pressure = 200; %MPa



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
        vol_f = zeros(size(T_f));
        rho_f = input.data(:,6);
        mass_f = input.data(:,5);
        vol_f = input.data(:,29);
        anorthite = input.data(:,32);
        albite = input.data(:,31);
    end
end

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
% LC = {'Color',color};
printing = 1; %print 1 not print 0
%% plotting
figure;
%subplot(2,3,1)
hold on;
plot(anorthite,T_f,'-k','linewidth',5)
%plot(1-albite,T_f,'--k','linewidth',5)
title('Phase Diagram', 'interpreter','latex','FontSize',18)
xlabel('Anorthite','interpreter','latex','FontSize',18)
ylabel('Temperature [$^\circ$C]','interpreter','latex','FontSize',18)

if printing
print('anorthite','-dpng')
end
