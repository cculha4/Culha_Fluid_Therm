%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Crystal Tracker%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

%% Folder Names
%folder{1} = '/Volumes/My Mac Passport/Stanford_2015/Research_main/Multiphase-2dCode/MultiphaseNS/Magma/Double Diffusion/Crystals and Temperature/Constant Crystal Volume/r001phi10mu95/Run 8';
%mainfolder = folder{1};      %where the figure is saved
startfilename_C{1} = 'data_magma_chamber_xtl_trc_ev_102_800_';
startfilename_T{1} = 'data_magma_chamber_trc_ev_102_800_';
simnum_C = numel(startfilename_C);
simnum_T = numel(startfilename_T);
alphabet = {'A.','B.','C.','D.','E.','F.','G.','H.','I.','J.','K.','L.','M.','N.','O.','P.','Q.','R.','S.','T.','U.','V.','W.','X.','Y.','Z.'};
%% domain properties
tmax = 1000; maximg = 200;
timestep = 101:167;
maxposxtls = 100000;
Tempx    = NaN(max(timestep-100),maxposxtls);
xx       = NaN(max(timestep-100),maxposxtls);
yy       = NaN(max(timestep-100),maxposxtls);
TempHist = NaN(max(timestep-100),maxposxtls);
TempHist = NaN(max(timestep-100),maxposxtls);
%% set printing options
% figname   = 'Case3_loss';
format    = '-dpng';
resl      = '-r200';
rend      = '-opengl';
printfig  = false;
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
%% prepare axes/borders dimensions
axh = 5;
axw = axh*2;
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


%% entering into each type of simulation
for i = 1:numel(startfilename_C)
    %foldername = folder{i};%char(strcat(fn(1),fn(2)));
    %cd(foldername)
    %% enter into each time
    for j = timestep(i,:)
        p = p + 1;
        str_C = {startfilename_C{i}  mat2str(j) '.dat'};
        str_T = {startfilename_T{i}  mat2str(j) '.dat'};
        filename_C=char(strcat(str_C(1),str_C(2),str_C(3)));
        filename_T=char(strcat(str_T(1),str_T(2),str_T(3)));
        importdata(filename_C,' ');
        input_C=ans;
        importdata(filename_T,' ');
        input_T=ans;
        if isempty(input_C)
            continue
        end
        if isempty(input_T)
            continue
        end
        Xo = 0; Yo = 0;
        Xo=input_C.data(:,1);
        Yo=input_C.data(:,2);
        Tracker=input_C.data(:,5)+1;
        Tempin=input_T.data(:,8);
        Xin = input_T.data(:,1);
        Yin=input_T.data(:,2);

        for k = 1:numel(Xo)
            distanceX = 0;
            distanceX = sqrt((Xin-Xo(k)).^2+(Yin-Yo(k)).^2);
            [orderx, nn] = sort(distanceX);
            closest_ns = nn(1:10); 
            closest_n = find(distanceX == min(distanceX));
            xx(j-100,Tracker(k)) = Xin(closest_n); yy(j-100,Tracker(k)) = Yin(closest_n); TempX(j-100,Tracker(k)) = mean(Tempin(closest_ns));
            TempHist(j-100,Tracker(k)) = TempX(j-100,Tracker(k));
        end

%             TempHistAv(j-100,:) = (TempHist(j-100,1:end-3)+TempHist(j-100,2:end-2)+TempHist(j-100,3:end-1)+TempHist(j-100,4:end))./4;
    end
end

time = (1:max(timestep-100))*tmax/maximg;
DTempHist = TempHist(2:end,:)-TempHist(1:end-1,:);
DTempHistDt = DTempHist/ (time(2)-time(1));
meanDT = 0;

for k = 1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        meanDT(k) = mean(DTempHist(~isnan(DTempHist(:,k)),k));
        firstvala = TempHist(~isnan(TempHist(:,k)),k);
        firstval(k) = firstvala(1);
        lastvala = TempHist(~isnan(TempHist(:,k)),k);
        lastval(k) = lastvala(end);
        BigDT(k) = lastval(k)- firstval(k);

    else
        BigDT(k) = NaN;
    end
end
%identify the crystals that made it to the end
for k = 1:numel(TempHist(1,:))
    if sum(~isnan(TempHist(:,k)))>2
        if ~isnan(TempHist(end,k))
        meanDT_end(k) = mean(DTempHist(~isnan(DTempHist(:,k)),k));
        firstvala = TempHist(~isnan(TempHist(:,k)),k);
        firstval_end(k) = firstvala(1);
        lastvala = TempHist(~isnan(TempHist(:,k)),k);
        lastval_end(k) = lastvala(end);
        BigDT_end(k) = lastval_end(k)- firstval_end(k);
        
        else
            BigDT_end(k) = NaN;
        end
    else
        BigDT_end(k) = NaN;
    end
end

save('Crystaltracker_old')
