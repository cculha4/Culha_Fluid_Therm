%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Crystal Movement%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load data
clear all
close all
load('Crystaltracker.mat')
v_T = sqrt(1e-7/750);
um = 4.2E-05;
a   = 0.002;
%clear variables that are repeated
small = []; big =[]; depth = [];

%define the final region (advecting closer)
small = xx(end,:)<0.3;
big = xx(end,:)>0.15;
L   = yy(1,:)<0.15;
R   = yy(1,:)>0.15;
depthL = find(small + big +L ==3);
depthR = find(small + big  +R==3);

%pick some crystals that end in the same area
wanted_crystalsL = [];
for k=depthL%1:numel(TempHist(1,:))
  if sum(~isnan(TempHist(:,k)))>50
       wanted_crystalsL = [wanted_crystalsL k];
   end
end
wanted_crystalsR = [];
for k=depthR%1:numel(TempHist(1,:))
  if sum(~isnan(TempHist(:,k)))>50
       wanted_crystalsR = [wanted_crystalsR k];
   end
end

%clear variables that are repeated
small = []; big =[]; depth = [];

%define the final region (advecting)
small = xx(end,:)<0.15;
big = xx(end,:)>0.00;
depth = find(small + big ==2);

%pick some crystals that end in the same area
wanted_crystalsT = [];
for k=depth%1:numel(TempHist(1,:))
  if sum(~isnan(TempHist(:,k)))>50
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
printfig  = false;


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
ahs = 0.1;
avs = 0.1;
axb = 0.7;
axt = 0.2;
axl = 0.8;
axr = 0.5;
cbh = axh; cbw = 0.2;
fh = axb + 1*axh + 1*ahs +           axt;
fw = axl + 1*axw + 1*avs + 1.5*cbw + axr;

% initialize figure and axes
for l = timestep(1,:)
str       = {'Step_bw_' mat2str(l)};
figname   = char(strcat(str(1),str(2)));

f = figure;
set(f,'Units','Inches','Position',[0.7 12 fw fh]);
set(f,'PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f,'Color','w','InvertHardcopy','off', 'MenuBar','none');
set(f,'Resize','off','Toolbar','none');
ax(1) = axes('Units','Inches','position',[axl         axb         axw axh]);

pcolcrystals = [55,126,184; 77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153]./255;
colcrystals = [pcolcrystals; pcolcrystals; pcolcrystals];

plot(yy(l-100,:),xx(l-100,:),'.','Color','k','MarkerSize',10)%[.65 .65 .65], 'MarkerSize',10)
hold on
xlabel('Width [m]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Height [m]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2-avs,axh/2]);
str       = {'t = ' mat2str(floor((l-100)*tmax/maximg*10)/10) 's'};
timetx   = char(strcat(str(1),str(2),str(3)));
text(1,1,timetx,UN{[1,2]},TX{:},FS{[1,4]},VA{[1,2]},HA{[1,4]},'Color',[0 0 0]);%FigDemoColorMap(30,:));
set(gca,'XLim',[0,.3],'Ylim',[0,.3],TL{:},FS{[1,4]});
hold off; drawnow;




if printfig
    print(f,format,resl,rend,figname,'-loose');
 end
    
end



