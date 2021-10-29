%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Crystal Anorthite Profile%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load data
clear all
close all
load('Crystaltracker_old.mat')
load('anorthite_syn.mat')


printfig  = false;
v_T = sqrt(1e-7/750);
um  = 8.0692e-6;
a   = 0.001;
tau = a/um;
growthrate = 10^(-3); %mu m/s
composition_conversion = 10^2/(10^4); %%only use this if your viscosity is O(2) different
%clear variables that are repeated
small = []; big =[]; depth = [];

%define the final region (advecting closer)
small = xx(end,:)<0.15;
big = xx(17,:)>0.15;
depth = find(small  == 1);

pcolcrystals = [166,206,227;
31,120,180;
178,223,138;
51,160,44;
251,154,153;
227,26,28;
253,191,111;
255,127,0;
202,178,214;
106,61,154;
255,255,153;
177,89,40]./255;
colcrystals = [pcolcrystals; pcolcrystals; pcolcrystals];

%pick some crystals that end in the same area
wanted_crystalsL = [];
for k=depth%1:numel(TempHist(1,:))
  if sum(~isnan(TempHist(:,k)))>50
       wanted_crystalsL = [wanted_crystalsL k];
   end
end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%FIGURES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set printing options

format    = '-dpng';
resl      = '-r200';
rend      = '-opengl';



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

plotmax = 6;
%Plotting profiles
axh = 1.25;
axw = 6;
ahs = 0.15;
avs = 0.25;
axb = 0.7;
axt = 0.5;
axl = 1;
axr = 0.5;
cbh = axh; cbw = 0.2; ht = 0.05;
fh = axb + plotmax*axh + plotmax*avs + axt ;
fw = axl + 1*axw + 1*ahs + 1*axr;
p = 0;

str_C       = {'crystal_temp_ND' };
figname   = char(strcat(str_C(1)));
f2 = figure;
set(f2,'Units','Inches','Position',[0.7 12 fw fh]);
set(f2,'PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f2,'Color','w','InvertHardcopy','off', 'MenuBar','none');
set(f2,'Resize','off','Toolbar','none');

for cnting = 1:plotmax
    ax(cnting) = axes('Units','Inches','position',[axl axb+(cnting-1)*(axh+avs) axw axh]);
end



cnt = 0;
v_an = zeros(100,numel(time));
for k=depth
    if sum(~isnan(TempHist(:,k)))>30
        cnt = cnt +1;
        if cnt >6
            continue
        end
        %smoothing
        for l = 1:size(~isnan(TempHist(:,k)))
            if l == 1
                TempHistsmooth(cnt,l) = (TempHist(1,k)*2 + TempHist(2,k))/3;
            elseif l == length(TempHist(:,k))
                TempHistsmooth(cnt,l) = (TempHist(end-1,k) + TempHist(end,k)*2)/3;
            else
                TempHistsmooth(cnt,l) = (TempHist(l-1,k)+TempHist(l,k)+TempHist(l+1,k))/3;
            end
        end

        axes(ax(cnt)); hold on;
        %shift TempHistsmooth to cooler temperatures
         Tempnew = (TempHistsmooth(cnt,:)-800)*(800-750)/(1040-800)+750;
        v_an(cnt,:) = polyval(p_an,Tempnew);
        v_an_xtl  = v_an(cnt,~isnan(v_an(cnt,:)));
        thickness = time./composition_conversion.*growthrate;
        lengthy = numel(v_an(~isnan(v_an(cnt,:))));
        thicknessy = thickness(1:lengthy);
        plot(time'./tau,(TempHistsmooth(cnt,:)-800)./(1040-800),'Color','k','LineWidth',2)
        set(gca,'Xlim',[0, max(time./tau)],TL{:},FS{[1,5]});
        if cnt == 1
            xlabel('time [non-dim]',TX{:},FS{[1,5]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
        else
            set(gca,'XTickLabel',[]);
        end
        set(gca,'Ylim',[0.5, 0.9],TL{:},FS{[1,5]});
        if cnt == plotmax/2
            ylabel('temperature [non-dim]',TX{:},FS{[1,5]},UN{[1,3]},'Position',[-axl/2,axh/2]);
        end
    end
end
        
if printfig
  %  cd(mainfolder);
     print(f2,format,resl,rend,figname,'-loose');
    % %
end



