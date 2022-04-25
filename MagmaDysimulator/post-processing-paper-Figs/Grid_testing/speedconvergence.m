% 119 <111>, 126 <211>, 146 <311>, 149 <411>

clear all
close all
stamps = [119, 127, 146, 149];
folder{1} = '/Volumes/My Mac Passport/Stanford_2015/Research_main/Multiphase-2dCode/MultiphaseNS/Magma/Double Diffusion/Crystals and Temperature/Cooling/Unzen/domainsize/111';
folder{2} = '/Volumes/My Mac Passport/Stanford_2015/Research_main/Multiphase-2dCode/MultiphaseNS/Magma/Double Diffusion/Crystals and Temperature/Cooling/Unzen/800Ctop/Run 8';
folder{3} = '/Volumes/My Mac Passport/Stanford_2015/Research_main/Multiphase-2dCode/MultiphaseNS/Magma/Double Diffusion/Crystals and Temperature/Cooling/Unzen/domainsize/311';
folder{4} = '/Volumes/My Mac Passport/Stanford_2015/Research_main/Multiphase-2dCode/MultiphaseNS/Magma/Double Diffusion/Crystals and Temperature/Cooling/Unzen/domainsize/411';
inty = [111, 211, 311, 411];
numb = [811, 800, 831, 841];
intx = inty;
xtl =0;
perc = 101; 
tmax = 1000; maximg = 200;
%nondimensional parameters
a  = 0.001;
L  = 0.1; 
n = 1.8;
vv_d = (2*(2600-2170)*a^2*9.8)/(9*1.5E+04*n);
vv_b = (2*(3000-2223)*a^2*9.8)/(9*9.5E+01*n);
v_crystal       = vv_b; %8.0692e-6;
v_temperature  = 1e-7/L;
ttc = L/v_crystal; 
ttt = L/v_temperature;
v  = v_crystal; 
tt = ttc;
Tcd = 940;%1000;
Tsd = 700;%737; %700;
Tcb = 1040;%1100;
Tsb = 800;%936; %700; 
Tc = Tcb;
Ts = Tsb;
dT = Tc-Ts;

meanvel = zeros(1,4);
maxvel = meanvel;
for j = 1:4
    cd(folder{j})
    str = {'data_magma_chamber_trc_ev_' mat2str(perc) '_' mat2str(numb(j)) '_' mat2str(stamps(j)) '.dat'};
    filename=char(strcat(str(1),str(2),str(3),str(4),str(5),str(6),str(7)));
 
    importdata(filename,' ');
    input=ans;
    close all
    if isempty(input)
        continue
    end
    Uin=input.data(:,4)./v;
    Vin=input.data(:,5)./v;
    PS1in = input.data(:,9);
    P1in = input.data(:,10);
    
    locDownwell = find(Uin<0);
    meanvel(j) = mean(Uin(locDownwell));
    minvel(j)  = min(Uin(locDownwell));
    stdvel(j) = std(Uin(locDownwell));
    
end




%% set printing options
format    = '-dpng';
resl      = '-r200';
rend      = '-opengl';
printfig  = true;
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




f = figure;
axh = 5;
axw = axh;
ahs = 0.15;
avs = 0.15;
axb = 0.7;
axt = 0.5;
axl = 1.2;
axr = 0.5;
cbh = axh; cbw = 0.2; ht = 0.05;
fh = axb + 1*axh + 1*avs + 1*ht + axt ;
fw = axl + axw + ahs +  axr;
set(f,'Units','Inches','Position',[axl 11 fw fh]);
set(f,'PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f,'Color','w','InvertHardcopy','off', 'MenuBar','none');
set(f,'Resize','off','Toolbar','none');
ax = axes('Units','Inches','position',[axl axb axw axh]);



errorbar(inty,abs(meanvel),stdvel,'LineWidth',2,'Color','b');
hold on
errorbar(inty,abs(minvel),stdvel,'LineWidth',2,'Color','r');
legend('Mean Downwelling Speed', 'Max Downwelling Speed',FS{[1,3]},TX{[1,2]})
xticks([111 211 311 411])
 xticklabels({'111x111','211x211','311x311','411x411'})


xlabel('Grid Resolution',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Speed [Non-dimensional]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);
set(gca,TL{:},FS{[1,4]});
if printfig
    print(f,format,resl,rend,'convergence_test','-loose');

end