%self created Phase diagram
close all;
clear all;

An = [0:.10:1];
Temp = [ 750, 755, 760, 770, 785, 800, 825, 850, 950, 1065, 1200];

p_an = polyfit(Temp,An,4);
v_an = polyval(p_an,Temp);

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
plot(An,Temp,'-k','linewidth',5)
plot(v_an,Temp,'--r','linewidth',5)
%plot(1-albite,T_f,'--k','linewidth',5)
title('Phase Diagram', 'interpreter','latex','FontSize',18)
xlabel(' Anorthite [mol $\%$] ','interpreter','latex','FontSize',18)
ylabel(' Temperature [$^{\circ}$C]','interpreter','latex','FontSize',18)
set(gca,'Ylim',[750, 1200],TL{:},FS{[1,5]});
set(gca,'Xlim',[0, 1],TL{:},FS{[1,5]});

if printing
print('anorthite_syn_model','-dpng')
end
save('anorthite_syn')