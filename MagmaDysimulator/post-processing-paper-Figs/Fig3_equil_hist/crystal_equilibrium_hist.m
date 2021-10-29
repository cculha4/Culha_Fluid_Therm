clear all
close all

load('Crystaltracker.mat')
NOCrystals = Tracker(end);
%non-dims
Tcd = 940;%1000;
Tsd = 700;%737; %700;
Tcb = 1100;
Tsb = 936; %700; 
Tc = Tcb;
Ts = Tsb;
dT = Tc-Ts;


printfig = 0;
f3 = figure;
axh = 5;
axw = axh;
ahs = 0.15;
avs = 0.15;
axb = 0.7;
axt = 0.5;
axl = 1.3;
axr = 0.5;
cbh = axh; cbw = 0.2; ht = 0.05;
fh = axb + 1*axh + 1*avs + 1*ht + axt ;
fw = axl + simnum_C*axw + (simnum_C-1)*ahs + (simnum_C-1)*cbw + axr;
set(f3,'Units','Inches','Position',[0.7 12 fw fh]);
set(f3,'PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f3,'Color','w','InvertHardcopy','off', 'MenuBar','none');
set(f3,'Resize','off','Toolbar','none');
ax = axes('Units','Inches','position',[axl axb axw axh]);
[freq_count, bin_value] = hist(BigDT./(max(Tempin)-min(Tempin)),30);
h1 = bar(bin_value,freq_count./NOCrystals,'FaceColor',[230,97,1]./255,'FaceAlpha',.9,'BarWidth', 1)
hold on
[freq_count, bin_value] = hist(BigDT_end./(max(Tempin)-min(Tempin)),30);
h2 = bar(bin_value,freq_count./NOCrystals,'FaceColor',[94,60,153]./255,'FaceAlpha',.9,'BarWidth', 1)

xlabel('Relative Max Temperature Gradient [non-dim]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Number of crystals',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);
set(gca,'XLim',[-.2, .2],TL{:},FS{[1,4]});
if printfig
    print(f3,format,resl,rend,'LargeTHist_nondim','-loose');
end
