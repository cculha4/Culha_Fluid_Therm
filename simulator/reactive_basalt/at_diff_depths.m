close all
clear all

load('Crystaltracker')
totalxtls = numel(Xo);
f = figure;

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
totalxtls_at_depth = [0];

%cooling
for k=depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>22
        if BigDT(k)<-10
            plot(time,TempHist(:,k),'-b','LineWidth',2)
            hold on
        end
    end
    totalxtls_at_depth = totalxtls_at_depth + 1;
end
%neutral tempearture
for k = depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)>-10 && BigDT(k)<10
            plot(time,TempHist(:,k),'-','LineWidth',2,'Color',[.5 .5 .5])
            hold on
        end
        
    end

end
%heating
for k = depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)>10
            plot(time,TempHist(:,k),'-','LineWidth',2,'Color',[1 0 0])
            hold on
        end
    end
end

xlabel('time [s]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Temperature [$\mathrm{^{\circ}}~$C]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);
set(gca,'XLim',[0,max(time)],'Ylim',[800,max(Tempin)],TL{:},FS{[1,4]});
str       = {'Number of crystals = ' mat2str(totalxtls_at_depth)};
NOxtls   = char(strcat(str(1),str(2)));
text(1,1,NOxtls,UN{[1,2]},TX{:},FS{[1,4]},VA{[1,2]},HA{[1,4]},'Color',[0 0 0]);
        
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
totalxtls_at_depth = [0];

%cooling
for k=depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)<-10
            plot(time,TempHist(:,k),'-b','LineWidth',2)
            hold on
        end
    end
totalxtls_at_depth = totalxtls_at_depth + 1;
end
%neutral tempearture
for k = depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)>-10 && BigDT(k)<10
            plot(time,TempHist(:,k),'-','LineWidth',2,'Color',[.5 .5 .5])
            hold on
        end
        
    end

end
%heating
for k = depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)>10
            plot(time,TempHist(:,k),'-','LineWidth',2,'Color',[1 0 0])
            hold on
        end
    end
end

xlabel('time [s]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Temperature [$\mathrm{^{\circ}}~$C]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);
set(gca,'XLim',[0,max(time)],'Ylim',[800,max(Tempin)],TL{:},FS{[1,4]});
        str       = {'Number of crystals = ' mat2str(totalxtls_at_depth)};
        NOxtls   = char(strcat(str(1),str(2)));

        text(1,1,NOxtls,UN{[1,2]},TX{:},FS{[1,4]},VA{[1,2]},HA{[1,4]},'Color',[0 0 0]);
                
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
totalxtls_at_depth = [0];

%cooling
for k=depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)<-10
            plot(time,TempHist(:,k),'-b','LineWidth',2)
            hold on
        end
    end
totalxtls_at_depth = totalxtls_at_depth + 1;
end
%neutral tempearture
for k = depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)>-10 && BigDT(k)<10
            plot(time,TempHist(:,k),'-','LineWidth',2,'Color',[.5 .5 .5])
            hold on
        end
        
    end

end
%heating
for k = depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>15
        if BigDT(k)>10
            plot(time,TempHist(:,k),'-','LineWidth',2,'Color',[1 0 0])
            hold on
        end
    end
end

xlabel('time [s]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Temperature [$\mathrm{^{\circ}}~$C]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);
set(gca,'XLim',[0,max(time)],'Ylim',[800,max(Tempin)],TL{:},FS{[1,4]});
        str       = {'Number of crystals = ' mat2str(totalxtls_at_depth)};
        NOxtls   = char(strcat(str(1),str(2)));
        text(1,1,NOxtls,UN{[1,2]},TX{:},FS{[1,4]},VA{[1,2]},HA{[1,4]},'Color',[0 0 0]);
                
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
totalxtls_at_depth = [0];

%cooling
for k=depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)<-10
            plot(time,TempHist(:,k),'-b','LineWidth',2)
            hold on
        end
    end
totalxtls_at_depth = totalxtls_at_depth + 1;
end
%neutral tempearture
for k = depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)>-10 && BigDT(k)<10
            plot(time,TempHist(:,k),'-','LineWidth',2,'Color',[.5 .5 .5])
            hold on
        end
        
    end

end
%heating
for k = depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)>10
            plot(time,TempHist(:,k),'-','LineWidth',2,'Color',[1 0 0])
            hold on
        end
    end
end

xlabel('time [s]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Temperature [$\mathrm{^{\circ}}~$C]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);
set(gca,'XLim',[0,max(time)],'Ylim',[800,max(Tempin)],TL{:},FS{[1,4]});
        str       = {'Number of crystals = ' mat2str(totalxtls_at_depth)};
        NOxtls   = char(strcat(str(1),str(2)));
        text(1,1,NOxtls,UN{[1,2]},TX{:},FS{[1,4]},VA{[1,2]},HA{[1,4]},'Color',[0 0 0]);
                
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
totalxtls_at_depth = [0];

%cooling
for k=depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)<-10
            plot(time,TempHist(:,k),'-b','LineWidth',2)
            hold on
        end
    end
totalxtls_at_depth = totalxtls_at_depth + 1;
end
%neutral tempearture
for k = depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)>-10 && BigDT(k)<10
            plot(time,TempHist(:,k),'-','LineWidth',2,'Color',[.5 .5 .5])
            hold on
        end
        
    end

end
%heating
for k = depth%1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if BigDT(k)>10
            plot(time,TempHist(:,k),'-','LineWidth',2,'Color',[1 0 0])
            hold on
        end
    end
end

xlabel('time [s]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Temperature [$\mathrm{^{\circ}}~$C]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);
set(gca,'XLim',[0,max(time)],'Ylim',[800,max(Tempin)],TL{:},FS{[1,4]});
        str       = {'Number of crystals = ' mat2str(totalxtls_at_depth)};
        NOxtls   = char(strcat(str(1),str(2)));
        text(1,1,NOxtls,UN{[1,2]},TX{:},FS{[1,4]},VA{[1,2]},HA{[1,4]},'Color',[0 0 0]);
                
if printfig
  %  cd(mainfolder);
    print(f2,format,resl,rend,figname,'-loose');
    % %
end
