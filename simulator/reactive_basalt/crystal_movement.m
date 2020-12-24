%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Crystal Movement%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load data
clear all
close all
load('Crystaltracker.mat')
%clear variables that are repeated
small = []; big =[]; depth = [];

%define the final region (advecting closer)
small = xx(end,:)<0.08;
big = xx(end,:)>0.06;
L   = yy(1,:)<0.05;
R   = yy(1,:)>0.05;
depthL = find(small + big +L ==3);
depthR = find(small + big  +R==3);

%pick some crystals that end in the same area
wanted_crystalsL = [];
for k=depthL%1:numel(TempHist(1,:))
  if sum(~isnan(TempHist(:,k)))>22
       wanted_crystalsL = [wanted_crystalsL k];
   end
end
wanted_crystalsR = [];
for k=depthR%1:numel(TempHist(1,:))
  if sum(~isnan(TempHist(:,k)))>22
       wanted_crystalsR = [wanted_crystalsR k];
   end
end

%clear variables that are repeated
small = []; big =[]; depth = [];

%define the final region (advecting)
small = xx(end,:)<0.04;
big = xx(end,:)>0.02;
depth = find(small + big ==2);

%pick some crystals that end in the same area
wanted_crystalsT = [];
for k=depth%1:numel(TempHist(1,:))
  if sum(~isnan(TempHist(:,k)))>22
       wanted_crystalsT = [wanted_crystalsT k];
   end
end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%FIGURES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set printing options

% figname   = 'Case3_loss';
format    = '-dpng';
resl      = '-r200';
rend      = '-opengl';
printfig  = true;


% prepare formating options
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
inty = 211; intx = 211;

% prepare axes/borders dimensions
axh = 3*2;
axw = axh*inty/intx;%4.5;
ahs = 0.15;
avs = 0.15;
axb = 0.7;
axt = 0.2;
axl = 0.8;
axr = 0.5;
cbh = axh; cbw = 0.2;
fh = axb + 1*axh + 1*avs +           axt;
fw = axl + 1*axw + 1*ahs + 1.5*cbw + axr;

% initialize figure and axes
for l = timestep(1,:)
str       = {'Step_c_' mat2str(l)};
figname   = char(strcat(str(1),str(2)));

f = figure;
set(f,'Units','Inches','Position',[0.7 12 fw fh]);
set(f,'PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f,'Color','w','InvertHardcopy','off', 'MenuBar','none');
set(f,'Resize','off','Toolbar','none');
ax(1) = axes('Units','Inches','position',[axl         axb         axw axh]);

pcolcrystals = [55,126,184; 77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153]./255;
colcrystals = [pcolcrystals; pcolcrystals; pcolcrystals];

plot(yy(l-100,:),xx(l-100,:),'.k','MarkerSize',32)
hold on
for k = 1:numel(wanted_crystalsL)
    plot(yy(l-100,wanted_crystalsL(k)),xx(l-100,wanted_crystalsL(k)),'.','MarkerSize',32,'Color',colcrystals(1,:))
end
for k = 1:numel(wanted_crystalsR)
    plot(yy(l-100,wanted_crystalsR(k)),xx(l-100,wanted_crystalsR(k)),'.','MarkerSize',32,'Color',colcrystals(2,:))
end
for k = 1:numel(wanted_crystalsT)
    plot(yy(l-100,wanted_crystalsT(k)),xx(l-100,wanted_crystalsT(k)),'.','MarkerSize',32,'Color',colcrystals(3,:))
end
xlabel('Width [m]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Height [m]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);
str       = {'t = ' mat2str(floor((l-100)*tmax/maximg)) 's'};
time   = char(strcat(str(1),str(2),str(3)));
text(1,1,time,UN{[1,2]},TX{:},FS{[1,4]},VA{[1,2]},HA{[1,4]},'Color',[0 0 0]);%FigDemoColorMap(30,:));
set(gca,'XLim',[0,.1],'Ylim',[0,.1],TL{:},FS{[1,4]});
hold off; drawnow;




if printfig
    print(f,format,resl,rend,figname,'-loose');
 end
    
end
