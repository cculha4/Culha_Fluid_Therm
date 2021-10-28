%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%Bulk density of the magma%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear all
close all
inty= 211; intx= 211; %x and y values
 xtl =0;
 rdx = 0.001;
 rhox = 3000;
 nav = 20;
  perc = 101; numb = 800;
workdir = '/Volumes/My Mac Passport/Stanford_2015/Research_main/Multiphase-2dCode/MultiphaseNS/Magma/Double Diffusion/Crystals and Temperature/Cooling/Unzen/800Ctop/Run 8';
 cd(workdir)
  meandens = 2.4127e+03; %reactive basalt
% meandens = 2.4024e+03; %nonreactive basalt
tmax = 1000; maximg = 200;
for l=[101:130]%101:100+maximg
 
           str = {'data_magma_chamber_trc_ev_' mat2str(perc) '_' mat2str(numb) '_' mat2str(l) '.dat'};
          filename=char(strcat(str(1),str(2),str(3),str(4),str(5),str(6),str(7)));
              str = {'data_magma_chamber_xtl_trc_ev_' mat2str(perc) '_' mat2str(numb) '_' mat2str(l) '.dat'};
          filenamextl=char(strcat(str(1),str(2),str(3),str(4),str(5),str(6),str(7)));

    importdata(filename,' ');
    input=ans;
    importdata(filenamextl,' ');
    inputxtls = ans;
    xo = inputxtls.data(:,1);
    yo = inputxtls.data(:,2);
    importdata(filename,' ');
    input=ans;
    close all
    if isempty(input)
        continue
    end
    Xin=input.data(:,1);
    Yin=input.data(:,2);
    dx = (max(Xin)-min(Xin))/intx; dy = (max(Yin)-min(Yin))/inty;
    Tempin=input.data(:,8);
    rhoin1=input.data(:,6);
    muin1 = input.data(:,7);
    %rhoin2=input.data(:,5);
    %muin2 = input.data(:,6);
    Uin=input.data(:,4);
    Vin=input.data(:,5);
    PS1in = input.data(:,9);
    P1in = input.data(:,10);
    
    [m,n]=size(input.data(1,:));
    
    
    
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
    
    x = X(:,1);
    y = Y(1,:);
    Solid_ratio = zeros(intx-1,inty-1);
    s_part = zeros(intx-1,inty-1);
    
    
    for k = 1:numel(xo)
        s_part_in = -rdx + sqrt((Xin-xo(k)).^2 + (Yin - yo(k)).^2);
        for j = 1:inty
            for i = 1:intx
                s_part(i,j) = s_part_in((j-1)*intx+i);
            end
        end
        Solid_ratio_update = volume_frac(s_part,1,intx,1,inty,x,y,rdx,xo(k),yo(k),dx,dy);
        Solid_ratio = Solid_ratio + Solid_ratio_update;
    end
    Solid_ratioF = zeros(intx,inty);
    Solid_ratioF(2:end,2:end) = Solid_ratio;
    density = rho1.*(1-Solid_ratioF) + rhox.*Solid_ratioF;
    Ker = 1/(nav*nav).*ones(nav,nav);
    covdens = conv2(density,Ker,'valid');
    avdens = covdens(1:nav:end,1:nav:end);
    [gradavdens_x, gradavdens_y] = gradient(avdens);
    perc_deviation_neutral = (avdens-meandens)./meandens.*100;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%FIGURES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set printing options
str       = {'Bulkdensity_pdev_' mat2str(l)};
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

% rmin = 2200;
% rmax=3000;
rmin = -10;%min(min(perc_deviation_neutral));
rmax = 10;%max(max(perc_deviation_neutral));
rmean = mean([rmin rmax]);
newmap = colormein(perc_deviation_neutral,[230,97,1; 247,247,247]./255,[247,247,247; 94,60,153]./255,[rmin rmean  rmax]);%, [800   990  1040]) ; %0 0 .5; 0 0 .5; 1 0 0;
 
 axes(ax(1)); hold on;
 perc_deviation_neutral(floor(intx/nav)+1,:) = rmin;
 perc_deviation_neutral(floor(intx/nav)+1,9) = rmax;
imagesc(y(1:end-nav)+nav*dy/2,x(1:end-nav)+nav*dx/2,perc_deviation_neutral)

colormap(newmap); box on;
 xlabel('Width [m]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
 ylabel('Height [m]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);
%  text(0.98,0.075,'$\%$ deviation from average density',UN{[1,2]},TX{:},FS{[1,4]},VA{[1,2]},HA{[1,4]},'Color',[0 0 0]);%FigDemoColorMap(30,:));
 strrho       = {'$\Delta_y \rho_{CV}$  = ' mat2str(floor(max(max(gradavdens_y(2:end-1,2:end-1))))) ' [kg/m$^3$]'};
 rhoprint   = char(strcat(strrho(1),strrho(2),strrho(3)));
%  text(0.98,0.02,rhoprint,UN{[1,2]},TX{:},FS{[1,4]},VA{[1,2]},HA{[1,4]},'Color',[0 0 0]);%FigDemoColorMap(30,:));

 str       = {'t = ' mat2str(floor((l-100)*tmax/maximg)) 's'};
 time   = char(strcat(str(1),str(2),str(3)));
 text(1,1,time,UN{[1,2]},TX{:},FS{[1,4]},VA{[1,2]},HA{[1,4]},'Color',[0 0 0]);%FigDemoColorMap(30,:));
%  caxis([2200  3000])
 colorbar('Limits',[rmin,rmax],'Ticks',[rmin rmean rmax],TL{:},FS{[1,3]},UN{[1,3]},'Position',[axl + 1*axw + 1*ahs ,axb,cbw,cbh]);
xlim([y(1) y(end-nav)])
ylim([x(1) x(end-nav)])




if printfig
    print(f,format,resl,rend,figname,'-loose');   
 end
    
end