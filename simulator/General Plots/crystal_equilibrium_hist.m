clear all
close all

load('Crystaltracker.mat')

f3 = figure;
axh = 5;
axw = axh;
ahs = 0.15;
avs = 0.15;
axb = 0.7;
axt = 0.5;
axl = 0.8;
axr = 0.5;
cbh = axh; cbw = 0.2; ht = 0.05;
fh = axb + 1*axh + 1*avs + 1*ht + axt ;
fw = axl + simnum_C*axw + (simnum_C-1)*ahs + (simnum_C-1)*cbw + axr;
set(f3,'Units','Inches','Position',[0.7 12 fw fh]);
set(f3,'PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f3,'Color','w','InvertHardcopy','off', 'MenuBar','none');
set(f3,'Resize','off','Toolbar','none');
ax = axes('Units','Inches','position',[axl axb axw axh]);
h1 = histogram(BigDT./(max(Tempin)-min(Tempin)),30,'FaceColor',[230,97,1]./255,'FaceAlpha',.7)
% h1.Normalization = 'probability';
h1.BinWidth = .05;
hold on
h2 = histogram(BigDT_end./(max(Tempin)-min(Tempin)),30,'FaceColor',[94,60,153]./255,'FaceAlpha',.7)
% h2.Normalization = 'probability';
h2.BinWidth = .05;

xlabel('Relative Max Temperature Gradient [non-dim]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Number of crystals',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);
set(gca,TL{:},FS{[1,4]});
if printfig
%     cd(mainfolder);
    print(f3,format,resl,rend,'LargeTHist_nondim','-loose');
    % %
end
