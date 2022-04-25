clear all
close all
inty= 210+1; intx= 210+1; 
 xtl =0;
  perc = 110; numb = 110;
tmax = 1000; maximg = 200;
%nondimensional parameters
a  = 0.001;
L  = 0.1; 
n = 1.8;
vv_d = (2*(2600-2170)*a^2*9.8)/(9*1.5E+04*n);
vv_b = (2*(3000-2223)*a^2*9.8)/(9*9.5E+01*n);
v_crystal       = vv_b; %8.0692e-6;
v_temperature  = 1e-7/L;
ttc = a/v_crystal; 
ttt = L/v_temperature;
v  = v_temperature; 
tt = ttt;
Tcd = 940;%1000;
Tsd = 700;%737; %700;
Tcb = 1040;%1100;
Tsb = 800;%936; %700; 
Tc = Tcb;
Ts = Tsb;
dT = Tc-Ts;

for l=[ 141]%[173]%101:100+maximg
 
           str = {'data_magma_chamber_trc_ev_' mat2str(perc) '_' mat2str(numb) '_' mat2str(l) '.dat'};
          filename=char(strcat(str(1),str(2),str(3),str(4),str(5),str(6),str(7)));
              str = {'data_magma_chamber_xtl_trc_ev_' mat2str(perc) '_' mat2str(numb) '_' mat2str(l) '.dat'};
          filenamextl=char(strcat(str(1),str(2),str(3),str(4),str(5),str(6),str(7)));

    importdata(filename,' ');
    input=ans;
    importdata(filenamextl,' ');
    inputxtls = ans;
    xo = inputxtls.data(:,1)./L;
    yo = inputxtls.data(:,2)./L;
    importdata(filename,' ');
    input=ans;
    close all
    if isempty(input)
        continue
    end
    Xin=input.data(:,1)./L;
    Yin=input.data(:,2)./L;
    Tempin=-(Ts-input.data(:,8))./dT;
    rhoin1=input.data(:,6);
    muin1 = input.data(:,7);
    %rhoin2=input.data(:,5);
    %muin2 = input.data(:,6);
    Uin=input.data(:,4)./v;
    Vin=input.data(:,5)./v;
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
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%FIGURES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set printing options
str       = {'Step_Thotyellow_nondimTemperature_Tdif_prof_' mat2str(l)};
figname   = char(strcat(str(1),str(2)));
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
  newmap = colormein(temp,[.6 0 0 ; .8 0.4 0],[.8 .4 0;1 1 0],([Ts mean([800 1040]) Tc]-Ts)./dT);%, [800   990  1040]) ; %0 0 .5; 0 0 .5; 1 0 0;
 minT = (Ts-Ts)./dT; maxT = -(Ts-Tc)./dT;

axes(ax(1)); hold on;
 contourf(Y,X,temp,100,'edgecolor','none')%100,'edgecolor','none')
 contour(Y,X,temp,-(Ts-[800 850 900 950 1000 1040])./dT,'LineWidth',4,'edgecolor',[.7 .7 .7])

 axis equal tight
%       contour(Y,X,phiS+990,[990,990],'LineWidth',1)
for z = 1:numel(yo)
    if numel(yo)>z
    if yo(z) >max(max(Y)) || yo(z) <min(min(Y)) || xo(z)>max(max(X)) || xo(z) <min(min(X))
        yo(z) = []; xo(z) =[];
    end
    end
end
        
 plot(yo,xo,'.k','MarkerSize',32)
colormap(newmap); box on;
xlabel('Width [NonDim]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Height [NonDim]',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);
str       = {'t = ' mat2str(floor((l-100)*tmax/maximg./tt*100)/100) '[NonDim]'};
time   = char(strcat(str(1),str(2),str(3)));
text(1,1,time,UN{[1,2]},TX{:},FS{[1,4]},VA{[1,2]},HA{[1,4]},'Color',[0 0 0]);%FigDemoColorMap(30,:));
caxis([minT maxT])
colorbar('Limits',[minT,maxT],'Ticks',[0 .2 .4 .6 .8 1],TL{:},FS{[1,3]},UN{[1,3]},'Position',[axl + 1*axw + 1*ahs ,axb,cbw,cbh]);
set(gca,TL{:},FS{[1,4]});
hold off; drawnow;




if printfig
    print(f,format,resl,rend,figname,'-loose');   
 end
    
end
