%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Crystal Movement Profile%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load data
clear all
close all
load('Crystaltracker_old.mat')
v_T = sqrt(1e-7/750);
%nondimensional parameters
a  = 0.001;
L  = 0.3; 
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
%clear variables that are repeated
small = []; big =[]; depth = [];

%define the final region (advecting closer)
small = xx(end,:)<0.3;
big = xx(end,:)>0.15;
L   = yy(1,:)<0.15;
R   = yy(1,:)>0.15;
depthL = find(small + big +L ==3);
depthR = find(small + big  +R==3);
pcolcrystals = [55,126,184; 77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153]./255;
colcrystals = [pcolcrystals; pcolcrystals; pcolcrystals];

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
inty = 211; intx = 211;


%Plotting profiles
axh = 8;
axw = axh;
ahs = 0.15;
avs = 0.15;
axb = 0.7;
axt = 0.5;
axl = 1;
axr = 0.5;
cbh = axh; cbw = 0.2; ht = 0.05;
fh = axb + 1*axh + 1*avs + 1*ht + axt ;
fw = axl + simnum_C*axw + (simnum_C-1)*ahs + (simnum_C-1)*cbw + axr;
p = 0;

str_C       = {'CrystalProfile_samedepthLRT_ND_all' };
figname   = char(strcat(str_C(1)));
f2 = figure;
set(f2,'Units','Inches','Position',[0.7 12 fw fh]);
set(f2,'PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f2,'Color','w','InvertHardcopy','off', 'MenuBar','none');
set(f2,'Resize','off','Toolbar','none');
ax = axes('Units','Inches','position',[axl axb axw axh]);




for k=depthL%[8 15 50])
    if sum(~isnan(TempHist(:,k)))>50
        %smoothing
        for l = 1:size(TempHist(:,k))
            if l == 1
                TempHistsmooth(l) = (TempHist(1,k)*2 + TempHist(2,k))/3;
            elseif l == length(TempHist(:,k))
                TempHistsmooth(l) = (TempHist(end-1,k) + TempHist(end,k)*2)/3;
            else
                TempHistsmooth(l) = (TempHist(l-1,k)+TempHist(l,k)+TempHist(l+1,k))/3;
            end
        end
       plot(time./tt,(TempHist(:,k)-Ts)./dT,'Color',[colcrystals(1,:) .2],'LineWidth',2)
       hold on
    end
end

for k=depthR%[5 20 6])
    if sum(~isnan(TempHist(:,k)))>50
%smoothing
        for l = 1:size(TempHist(:,k))
            if l == 1
                TempHistsmooth(l) = (TempHist(1,k)*2 + TempHist(2,k))/3;
            elseif l == length(TempHist(:,k))
                TempHistsmooth(l) = (TempHist(end-1,k) + TempHist(end,k)*2)/3;
            else
                TempHistsmooth(l) = (TempHist(l-1,k)+TempHist(l,k)+TempHist(l+1,k))/3;
            end
        end
        plot(time./tt,(TempHist(:,k)-Ts)./dT,'Color',[colcrystals(2,:) .2],'LineWidth',2)
       hold on
    end
end

for k=depth   % feel free to specify these crystals [ 3 4        5        6  ])
    if sum(~isnan(TempHist(:,k)))>50
%smoothing
       for l = 1:size(TempHist(:,k))
            if l == 1
                TempHistsmooth(l) = (TempHist(1,k)*2 + TempHist(2,k))/3;
            elseif l == length(TempHist(:,k))
                TempHistsmooth(l) = (TempHist(end-1,k) + TempHist(end,k)*2)/3;
            else
                TempHistsmooth(l) = (TempHist(l-1,k)+TempHist(l,k)+TempHist(l+1,k))/3;
            end
        end
        plot(time./tt,(TempHist(:,k)-Ts)./dT,'Color',[colcrystals(3,:) .2],'LineWidth',2)
       hold on
    end
end



xlabel('time',TX{:},FS{[1,5]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Temperature',TX{:},FS{[1,6]},UN{[1,3]},'Position',[-axl/2,axh/2]);
set(gca,'Ylim',[0, 1],TL{:},FS{[1,6]});
        
if printfig
     print(f2,format,resl,rend,'All_BGP_crystals','-loose');
end


f3 = figure;
set(f3,'Units','Inches','Position',[0.7 12 fw fh]);
set(f3,'PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f3,'Color','w','InvertHardcopy','off', 'MenuBar','none');
set(f3,'Resize','off','Toolbar','none');
ax = axes('Units','Inches','position',[axl axb axw axh]);




for k=depthL([15 50])
    if sum(~isnan(TempHist(:,k)))>50
        %smoothing
        for l = 1:size(TempHist(:,k))
            if l == 1
                TempHistsmooth(l) = (TempHist(1,k)*2 + TempHist(2,k))/3;
            elseif l == length(TempHist(:,k))
                TempHistsmooth(l) = (TempHist(end-1,k) + TempHist(end,k)*2)/3;
            else
                TempHistsmooth(l) = (TempHist(l-1,k)+TempHist(l,k)+TempHist(l+1,k))/3;
            end
        end
       plot(time./tt,(TempHist(:,k)-Ts)./dT,'Color',colcrystals(1,:),'LineWidth',2)
       hold on
    end
end


for k=depthR([6 20])
    if sum(~isnan(TempHist(:,k)))>50
%smoothing
        for l = 1:size(TempHist(:,k))
            if l == 1
                TempHistsmooth(l) = (TempHist(1,k)*2 + TempHist(2,k))/3;
            elseif l == length(TempHist(:,k))
                TempHistsmooth(l) = (TempHist(end-1,k) + TempHist(end,k)*2)/3;
            else
                TempHistsmooth(l) = (TempHist(l-1,k)+TempHist(l,k)+TempHist(l+1,k))/3;
            end
        end
        plot(time./tt,(TempHist(:,k)-Ts)./dT,'Color',colcrystals(2,:),'LineWidth',2)
       hold on
    end
end
for k=depth([3 5])
    if sum(~isnan(TempHist(:,k)))>50
%smoothing
       for l = 1:size(TempHist(:,k))
            if l == 1
                TempHistsmooth(l) = (TempHist(1,k)*2 + TempHist(2,k))/3;
            elseif l == length(TempHist(:,k))
                TempHistsmooth(l) = (TempHist(end-1,k) + TempHist(end,k)*2)/3;
            else
                TempHistsmooth(l) = (TempHist(l-1,k)+TempHist(l,k)+TempHist(l+1,k))/3;
            end
        end
        plot(time./tt,(TempHist(:,k)-Ts)./dT,'Color',colcrystals(3,:),'LineWidth',2)
       hold on
    end
end


xlabel('time',TX{:},FS{[1,5]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Temperature',TX{:},FS{[1,6]},UN{[1,3]},'Position',[-axl/2,axh/2]);
set(gca,'Ylim',[0, 1],TL{:},FS{[1,6]});
        
if printfig
     print(f3,format,resl,rend,figname,'-loose');
end