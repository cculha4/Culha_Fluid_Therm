clear all
close all
inty= 211; intx= 211; %x and y values
 xtl =0;
  perc = 101; numb = 850;
%workdir = '/Volumes/My Mac Passport/Stanford_2015/Research_main/Multiphase-2dCode/MultiphaseNS/Magma/Double Diffusion/Crystals and Temperature/Constant Crystal Number/InitialMeltsDist';
% cd(workdir)
tmax = 5000; maximg = 200;
for l=260
 
           str = {'data_magma_chamber_trc_ev_' mat2str(perc) '_' mat2str(numb) '_' mat2str(l) '.dat'};
          filename=char(strcat(str(1),str(2),str(3),str(4),str(5),str(6),str(7)));
              str = {'data_magma_chamber_xtl_trc_ev_' mat2str(perc) '_' mat2str(numb) '_' mat2str(l) '.dat'};
          filenamextl=char(strcat(str(1),str(2),str(3),str(4),str(5),str(6),str(7)));

    importdata(filename,' ');
    input=ans;
    importdata(filenamextl,' ');
    inputxtls = ans;
    if ~isstruct(inputxtls)
        xo = []; yo = []; tkr = [];
    else
        xo = inputxtls.data(:,1);
        yo = inputxtls.data(:,2);
        tkr = inputxtls.data(:,5);
    end

    importdata(filename,' ');
    input=ans;
    close all
    if isempty(input)
        continue
    end
    Xin=input.data(:,1);
    Yin=input.data(:,2);
    Tempin=input.data(:,8);
    rhoin1=input.data(:,6);
    muin1 = input.data(:,7);
    %rhoin2=input.data(:,5);
    %muin2 = input.data(:,6);
    Uin=input.data(:,4);
    Vin=input.data(:,5);
    PS1in = input.data(:,9);
    P1in = input.data(:,10);
    
    [m,n]=size(input(1,:));
    
    
    
    for j=1:inty
        for i=1:intx
            U(i,j)=Uin((j-1)*intx+i);
            V(i,j)=Vin((j-1)*intx+i);
            X(i,j)=Xin((j-1)*intx+i);
            Y(i,j)=Yin((j-1)*intx+i);
            temp(i,j)=Tempin((j-1)*intx+i);
            rho1(i,j)=rhoin1((j-1)*intx+i);
            mu1(i,j) = muin1((j-1)*intx+i);
            %rho2(i,j)=rhoin2((j-1)*intx+i);
            %mu2(i,j) = muin2((j-1)*intx+i);
            phiS(i,j) = PS1in((j-1)*intx+i);
            phi(i,j) = P1in((j-1)*intx+i);
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%FIGURES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set printing options
str       = {'Step_BW_' mat2str(l)};
figname   = char(strcat(str(1),str(2)));
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
f = figure;
set(f,'Units','Inches','Position',[0.7 12 fw fh]);
set(f,'PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f,'Color','w','InvertHardcopy','off', 'MenuBar','none');
set(f,'Resize','off','Toolbar','none');
ax(1) = axes('Units','Inches','position',[axl         axb         axw axh]);

% newmap = colormein(temp,[[1 1 0]; [.5 .5 .5]; .95 .75 0; [.5 .5 .5];.8 .45 0;],...
%                         [[.95 .75 0];[.5 .5 .5]; .8 .45 0; [.5 .5 .5];.6 0 0],...
%                         [800 875 885 965 975 1040]);%, [800   990  1040]) ; %0 0 .5; 0 0 .5; 1 0 0;
% [254,240,217]./255 [253,204,138]./255 [252,141,89]./255 [215,48,31]./255
  newmap = colormein(temp,[1 1 0; .8 0.4 0],[.8 .4 0; .6 0 0],[800 990 1040]);%, [800   990  1040]) ; %0 0 .5; 0 0 .5; 1 0 0;
 minT = 700; maxT = 1040;

axes(ax(1)); hold on;
%contourf(Y,X,temp,100,'edgecolor','none')%100,'edgecolor','none')
 %contour(Y,X,temp,[800 850 900 950 1000 1040],'LineWidth',4,'edgecolor',[.7 .7 .7])%100,'edgecolor','none')

 axis equal tight
%  contour(Y,X,phi+990,[990,990],'LineWidth',1)
%  contour(Y,X,phiS+990,[990,990],'LineWidth',1)
for z = 1:numel(yo)
    if numel(yo)>z
    if yo(z) >max(max(Y)) || yo(z) <min(min(Y)) || xo(z)>max(max(X)) || xo(z) <min(min(X))
        yo(z) = []; xo(z) =[];
    end
    end
end

crystals = [45 79 94 95];
colcrystals = [55,126,184; 77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153]./255;
 plot(yo,xo,'.k','MarkerSize',35)

%colormap(newmap); box on;
xlabel('Width [m]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Height [m]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);
%text(0.98,0.02,'Temperature [$\mathrm{^{\circ}}~$C]',UN{[1,2]},TX{:},FS{[1,4]},VA{[1,2]},HA{[1,4]},'Color',[1 1 1]);%FigDemoColorMap(30,:));
str       = {'t = ' mat2str(floor((l-100)*tmax/maximg)) 's'};
time   = char(strcat(str(1),str(2),str(3)));
text(1,1,time,UN{[1,2]},TX{:},FS{[1,4]},VA{[1,2]},HA{[1,4]},'Color',[0 0 0]);%FigDemoColorMap(30,:));
%caxis([minT maxT])
%colorbar('Limits',[minT,maxT],'Ticks',[ceil(minT),800, 900, 1000,maxT],TL{:},FS{[1,3]},UN{[1,3]},'Position',[axl + 1*axw + 1*ahs ,axb,cbw,cbh]);
%     axl+axw+ahs+cbw,axb+axh+avs,cbw,cbh]);
set(gca,'XLim',[min(Yin), max(Yin)], 'YLim', [min(Xin) max(Xin)],FS{[1,3]})
hold off; drawnow;




if printfig
%      cd(workdir);
    print(f,format,resl,rend,figname,'-loose');   
 end
    
end
