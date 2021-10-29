 close all
clear all
load('crystaltracker.mat')
printfig = 1;
%% Plot crystal Temperature profile
%% Non-dimensional form
a  = 0.001;
n = 1.8;
L = .1;
vv_d = (2*(2600-2170)*a^2*9.8)/(9*1.5E+04*n);
vv_b = (2*(3000-2223)*a^2*9.8)/(9*9.5E+01*n);
v_crystal       = vv_d; %8.0692e-6;
v_temperature  = 1e-7;
ttc = L/v_crystal; 
ttt = L^2/v_temperature;
v  = v_crystal; 
tt = ttt;
Tcd = 940;%1000;
Tsd = 700;%737; %700;
Tcb = 1040;%1100;
Tsb = 800;%936; %700; 
Tc = Tcd;
Ts = Tsd;
dT = Tc-Ts;
%% prepare axes/borders dimensions
axh = 5;
axw = axh*1.25;
ahs = 0.15;
avs = 0.15;
axb = 0.7;
axt = 0.5;
axl = 0.8;
axr = 0.5;
cbh = axh; cbw = 0.2; ht = 0.05;
fh = axb + 1*axh + 1*avs + 1*ht + axt ;
fw = axl + simnum_C*axw + (simnum_C-1)*ahs + (simnum_C-1)*cbw + axr;
p = 0;
%% initialize figure and axes
str_C       = {'CrystalProfile' };
figname   = char(strcat(str_C(1)));
f = figure;
set(f,'Units','Inches','Position',[0.7 12 fw fh]);
set(f,'PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f,'Color','w','InvertHardcopy','off', 'MenuBar','none');
set(f,'Resize','off','Toolbar','none');
ax = axes('Units','Inches','position',[axl axb axw axh]);


%cooling
for k=1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)<-10
            plot(time./tt,(TempHist(:,k)-Ts)/dT,'-b','LineWidth',2)
            hold on
        end
    end
end
%neutral tempearture
for k = 1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)>-10 && BigDT(k)<10
            plot(time./tt,(TempHist(:,k)-Ts)/dT,'-','LineWidth',2,'Color',[.5 .5 .5])
            hold on
        end
        
    end

end
%heating
for k = 1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)>10
            plot(time./tt,(TempHist(:,k)-Ts)/dT,'-','LineWidth',2,'Color',[1 0 0])
            hold on
        end
    end
end

xlabel('time',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Temperature',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);
set(gca,'XLim',[0,max(time)./tt],'Ylim',([Ts,max(Tempin)]-Ts)./dT,TL{:},FS{[1,4]});
        
if printfig
%cd(mainfolder);
    print(f,format,resl,rend,figname,'-loose');
    % %
end
%% initialize figure and axes
%% prepare axes/borders dimensions
axh = 5;
axw = axh*1.25;
ahs = 0.15;
avs = 0.15;
axb = 0.7;
axt = 0.5;
axl = 0.8;
axr = 0.5;
cbh = axh; cbw = 0.2; ht = 0.05;
fh = axb + 1*axh + 1*avs + 1*ht + axt ;
fw = axl + simnum_C*axw + (simnum_C-1)*ahs + (simnum_C-1)*cbw + axr;
p = 0;

str_C       = {'CrystalProfile_samedepth_6_8' };
figname   = char(strcat(str_C(1)));
f2 = figure;
set(f2,'Units','Inches','Position',[0.7 12 fw fh]);
set(f2,'PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f2,'Color','w','InvertHardcopy','off', 'MenuBar','none');
set(f2,'Resize','off','Toolbar','none');
ax = axes('Units','Inches','position',[axl axb axw axh]);

%identify the crystals that are in the same dpeth range
small = xx(end,:)<0.08;
big = xx(end,:)>0.06;
depth = find(small + big ==2);

%cooling
for k=depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)<-10
            plot(time./tt,(TempHist(:,k)-Ts)./dT,'-b','LineWidth',2)
            hold on
        end
    end
end
%neutral tempearture
for k = depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)>-10 && BigDT(k)<10
            plot(time./tt,(TempHist(:,k)-Ts)./dT,'-','LineWidth',2,'Color',[.5 .5 .5])
            hold on
        end
        
    end

end
%heating
for k = depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)>10
            plot(time./tt,(TempHist(:,k)-Ts)./dT,'-','LineWidth',2,'Color',[1 0 0])
            hold on
        end
    end
end

xlabel('time',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Temperature',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);
set(gca,'XLim',[0,max(time)./tt],'Ylim',([Ts,max(Tempin)]-Ts)./dT,TL{:},FS{[1,4]});
        
if printfig
  %  cd(mainfolder);
    print(f2,format,resl,rend,figname,'-loose');
    % %
end






axh = 5;
axw = axh*1.25;
ahs = 0.15;
avs = 0.15;
axb = 0.7;
axt = 0.5;
axl = 0.8;
axr = 0.5;
cbh = axh; cbw = 0.2; ht = 0.05;
fh = axb + 1*axh + 1*avs + 1*ht + axt ;
fw = axl + simnum_C*axw + (simnum_C-1)*ahs + (simnum_C-1)*cbw + axr;
p = 0;

str_C       = {'CrystalProfile_samedepth_4_6' };
figname   = char(strcat(str_C(1)));
f2 = figure;
set(f2,'Units','Inches','Position',[0.7 12 fw fh]);
set(f2,'PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f2,'Color','w','InvertHardcopy','off', 'MenuBar','none');
set(f2,'Resize','off','Toolbar','none');
ax = axes('Units','Inches','position',[axl axb axw axh]);

%identify the crystals that are in the same dpeth range
small = xx(end,:)<0.06;
big = xx(end,:)>0.04;
depth = find(small + big ==2);

%cooling
for k=depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)<-10
            plot(time./tt,(TempHist(:,k)-Ts)./dT,'-b','LineWidth',2)
            hold on
        end
    end
end
%neutral tempearture
for k = depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)>-10 && BigDT(k)<10
            plot(time./tt,(TempHist(:,k)-Ts)./dT,'-','LineWidth',2,'Color',[.5 .5 .5])
            hold on
        end
        
    end

end
%heating
for k = depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)>10
            plot(time./tt,(TempHist(:,k)-Ts)./dT,'-','LineWidth',2,'Color',[1 0 0])
            hold on
        end
    end
end

xlabel('time',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Temperature',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);
set(gca,'XLim',[0,max(time)./tt],'Ylim',([Ts,max(Tempin)]-Ts)./dT,TL{:},FS{[1,4]});
        
if printfig
  %  cd(mainfolder);
    print(f2,format,resl,rend,figname,'-loose');
    % %
end




axh = 5;
axw = axh*1.25;
ahs = 0.15;
avs = 0.15;
axb = 0.7;
axt = 0.5;
axl = 0.8;
axr = 0.5;
cbh = axh; cbw = 0.2; ht = 0.05;
fh = axb + 1*axh + 1*avs + 1*ht + axt ;
fw = axl + simnum_C*axw + (simnum_C-1)*ahs + (simnum_C-1)*cbw + axr;
p = 0;

str_C       = {'CrystalProfile_samedepth_2_4' };
figname   = char(strcat(str_C(1)));
f2 = figure;
set(f2,'Units','Inches','Position',[0.7 12 fw fh]);
set(f2,'PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f2,'Color','w','InvertHardcopy','off', 'MenuBar','none');
set(f2,'Resize','off','Toolbar','none');
ax = axes('Units','Inches','position',[axl axb axw axh]);

%identify the crystals that are in the same dpeth range
small = xx(end,:)<0.04;
big = xx(end,:)>0.02;
depth = find(small + big ==2);

%cooling
for k=depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)<-10
            plot(time./tt,(TempHist(:,k)-Ts)./dT,'-b','LineWidth',2)
            hold on
        end
    end
end
%neutral tempearture
for k = depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)>-10 && BigDT(k)<10
            plot(time./tt,(TempHist(:,k)-Ts)./dT,'-','LineWidth',2,'Color',[.5 .5 .5])
            hold on
        end
        
    end

end
%heating
for k = depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>15
        if BigDT(k)>10
            plot(time./tt,(TempHist(:,k)-Ts)./dT,'-','LineWidth',2,'Color',[1 0 0])
            hold on
        end
    end
end

xlabel('time',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Temperature',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);
set(gca,'XLim',[0,max(time)./tt],'Ylim',([Ts,max(Tempin)]-Ts)./dT,TL{:},FS{[1,4]});
        
if printfig
  %  cd(mainfolder);
    print(f2,format,resl,rend,figname,'-loose');
    % %
end








axh = 5;
axw = axh*1.25;
ahs = 0.15;
avs = 0.15;
axb = 0.7;
axt = 0.5;
axl = 0.8;
axr = 0.5;
cbh = axh; cbw = 0.2; ht = 0.05;
fh = axb + 1*axh + 1*avs + 1*ht + axt ;
fw = axl + simnum_C*axw + (simnum_C-1)*ahs + (simnum_C-1)*cbw + axr;
p = 0;

str_C       = {'CrystalProfile_samedepth_8_10' };
figname   = char(strcat(str_C(1)));
f2 = figure;
set(f2,'Units','Inches','Position',[0.7 12 fw fh]);
set(f2,'PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f2,'Color','w','InvertHardcopy','off', 'MenuBar','none');
set(f2,'Resize','off','Toolbar','none');
ax = axes('Units','Inches','position',[axl axb axw axh]);

%identify the crystals that are in the same dpeth range
small = xx(end,:)<0.1;
big = xx(end,:)>0.08;
depth = find(small + big ==2);

%cooling
for k=depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)<-10
            plot(time./tt,(TempHist(:,k)-Ts)./dT,'-b','LineWidth',2)
            hold on
        end
    end
end
%neutral tempearture
for k = depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)>-10 && BigDT(k)<10
            plot(time./tt,(TempHist(:,k)-Ts)./dT,'-','LineWidth',2,'Color',[.5 .5 .5])
            hold on
        end
        
    end

end
%heating
for k = depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)>10
            plot(time./tt,(TempHist(:,k)-Ts)./dT,'-','LineWidth',2,'Color',[1 0 0])
            hold on
        end
    end
end

xlabel('time ',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Temperature',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);
set(gca,'XLim',[0,max(time)./tt],'Ylim',([Ts,max(Tempin)]-Ts)./dT,TL{:},FS{[1,4]});
        
if printfig
  %  cd(mainfolder);
    print(f2,format,resl,rend,figname,'-loose');
    % %
end





axh = 5;
axw = axh*1.25;
ahs = 0.15;
avs = 0.15;
axb = 0.7;
axt = 0.5;
axl = 0.8;
axr = 0.5;
cbh = axh; cbw = 0.2; ht = 0.05;
fh = axb + 1*axh + 1*avs + 1*ht + axt ;
fw = axl + simnum_C*axw + (simnum_C-1)*ahs + (simnum_C-1)*cbw + axr;
p = 0;

str_C       = {'CrystalProfile_samedepth_0_2' };
figname   = char(strcat(str_C(1)));
f2 = figure;
set(f2,'Units','Inches','Position',[0.7 12 fw fh]);
set(f2,'PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f2,'Color','w','InvertHardcopy','off', 'MenuBar','none');
set(f2,'Resize','off','Toolbar','none');
ax = axes('Units','Inches','position',[axl axb axw axh]);

%identify the crystals that are in the same dpeth range
small = xx(end,:)<0.02;
big = xx(end,:)>0.00;
depth = find(small + big ==2);

%cooling
for k=depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)<-10
            plot(time./tt,(TempHist(:,k)-Ts)./dT,'-b','LineWidth',2)
            hold on
        end
    end
end
%neutral tempearture
for k = depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)>-10 && BigDT(k)<10
            plot(time./tt,(TempHist(:,k)-Ts)./dT,'-','LineWidth',2,'Color',[.5 .5 .5])
            hold on
        end
        
    end

end
%heating
for k = depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)>10
            plot(time./tt,(TempHist(:,k)-Ts)./dT,'-','LineWidth',2,'Color',[1 0 0])
            hold on
        end
    end
end

xlabel('time',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Temperature',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);
set(gca,'XLim',[0,max(time)./tt],'Ylim',([Ts,max(Tempin)]-Ts)./dT,TL{:},FS{[1,4]});
        
if printfig
  %  cd(mainfolder);
    print(f2,format,resl,rend,figname,'-loose');
    % %
end




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
h1 = histogram((BigDT)./dT,30,'FaceColor',[230,97,1]./255,'FaceAlpha',.7)

hold on
h2 = histogram((BigDT_end)./dT,30,'FaceColor',[94,60,153]./255,'FaceAlpha',.7)


xlabel('Max Temperature Gradient ',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Number of crystals',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);
set(gca,TL{:},FS{[1,4]});

