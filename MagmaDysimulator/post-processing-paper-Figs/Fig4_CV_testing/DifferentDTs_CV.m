%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Looking at how dT changes with CV size%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
%%%% at CV size of 10
folder{1} = '/Volumes/My Mac Passport/Stanford_2015/Research_main/Multiphase-2dCode/MultiphaseNS/Magma/Double Diffusion/Crystals and Temperature/Cooling/Unzen/CV test/CV10/Run 2';
 mainfolder = folder{1};      %where the figure is saved
 cd(mainfolder)
input = load('Crystaltracker.mat','BigDT','BigDT_end');
BigDT = input.BigDT./(1040-800); %nondimensionalize
BigDT_end = input.BigDT_end./(1040-800); %nondimensionalize
mBigDT(1) = mean(BigDT(~isnan(BigDT)));
stdBigDT(1) = std(BigDT(~isnan(BigDT)));
mBigDT_end(1) = mean(BigDT_end(~isnan(BigDT_end)));
stdBigDT_end(1) = std(BigDT_end(~isnan(BigDT_end)));

folder{2} = '/Volumes/My Mac Passport/Stanford_2015/Research_main/Multiphase-2dCode/MultiphaseNS/Magma/Double Diffusion/Crystals and Temperature/Cooling/Unzen/800Ctop/Run 8';
mainfolder = folder{2};      %where the figure is saved
cd(mainfolder)
input = load('Crystaltracker.mat','BigDT','BigDT_end');
BigDT = input.BigDT./(1040-800);
BigDT_end = input.BigDT_end./(1040-800);
mBigDT(2) = mean(BigDT(~isnan(BigDT)));
stdBigDT(2) = std(BigDT(~isnan(BigDT)));
mBigDT_end(2) = mean(BigDT_end(~isnan(BigDT_end)));
stdBigDT_end(2) = std(BigDT_end(~isnan(BigDT_end)));

folder{3} = '/Volumes/My Mac Passport/Stanford_2015/Research_main/Multiphase-2dCode/MultiphaseNS/Magma/Double Diffusion/Crystals and Temperature/Cooling/Unzen/CV test/CV30/Run 2';
 mainfolder = folder{3};      %where the figure is saved
 cd(mainfolder)
input = load('Crystaltracker.mat','BigDT','BigDT_end');
BigDT = input.BigDT./(1040-800);
BigDT_end = input.BigDT_end./(1040-800);
mBigDT(3) = mean(BigDT(~isnan(BigDT)));
stdBigDT(3) = std(BigDT(~isnan(BigDT)));
mBigDT_end(3) = mean(BigDT_end(~isnan(BigDT_end)));
stdBigDT_end(3) = std(BigDT_end(~isnan(BigDT_end)));



%% set printing options
% figname   = 'Case3_loss';
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
% LC = {'Color',color};




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
errorbar([10 20 30], mBigDT,stdBigDT,'LineWidth',2,'Color',[230,97,1]./255);
hold on
errorbar([10 20 30], mBigDT_end,stdBigDT,'LineWidth',2,'Color',[94,60,153]./255);
legend('In $\&$ Out of Equilibrium', 'In Equilibrium',FS{[1,3]},TX{[1,2]})
xlim([5 35])  

 xticks([10 20 30])
 xticklabels({'10x10','20x20','30x30'})
 ylim([-0.15 0.075])

xlabel('Control Volume Size',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Average Change in Temperature',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);
set(gca,TL{:},FS{[1,4]});
if printfig
    print(f,format,resl,rend,'CV_test','-loose');
    % %
end