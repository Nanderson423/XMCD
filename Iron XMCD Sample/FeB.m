% For Plotting the Hysteresis

clearvars
close all

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% CONSTANTS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
nh = 10-6.61;%2; % Number of holes 
L = 0; %Angular Momentum
S = 5/2; %Spin
J = L+S; %The shell is more than half filled so J=L+S
g = 1+(J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1)); %Lande factor given by g=1+(J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1))
kB = 8.6173324e-5; %Boltzmann constant eV/K
mu_B = 5.7883818066e-5; %Bohr magneton eV/T
%%%%%%%%%%%%%%%%%%%%%%%

%Scan numbers

%12 degrees
Fe6_12 = [213 230 238 243 248 253 258 266 271 276 281];
Fe1_12 = [214 231 239 244 249 254 259 267 272 277 282];

%90 degrees
Fe6_90 = [];
Fe1_90 = [];

%Initializes fileName
fileName = {};

for scans = Fe1_12
    fileName{end+1} = strcat('\APS 2018 Feb 4IDC\NAA_Feb18.',num2str(sprintf('%04d',scans)),'.dat');
end

%Gets fileLength
fileLength = length(fileName);

skipPlot = 0;
scanType = 0;
%StepE = [700 714;714 728];
StepE = [703 710; 714 728];
backscale=.5;

%Preallocate space for figure handle holder figs
figs = cell(1,fileLength+3);
%Preallocates space for storing data from the files
filedata = cell(1,fileLength);
%Preallocates space for Bfield
Bfield = zeros(1,fileLength);
%Preallocates space for Temp
Temp = zeros(1,fileLength);
%Preallocates space for zPos
zPos = zeros(1,fileLength);
%Preallocates space for sampleType
sampleType = cell(1,fileLength);
%Gets color order
color = get(groot,'DefaultAxesColorOrder');
color(3,:) = [0 .5 0];
%Gets monitor positions
set(0,'units','pixels');
monPos = get(0,'MonitorPosition');
printplot = 0;

if scanType == 0
    scanTypeStr = 'TEY';
elseif scanType == 1
    scanTypeStr = 'TFY';
elseif scanType == 2
    scanTypeStr = 'REF';
end


for n=1:fileLength
    
    data=LoadData(fileName{n});
    
    %Sets the long format to prevent/reduce rounding errors
    format long e
    
    %Pulls data from specific column into indentifiers 5/9 for TEY, 6/10
    %for TFY and 7/11 for REF
    data.Energy = flip(data.raw{2});
    data.monR = flip(data.raw{4});
    data.teyR = flip(data.raw{scanType + 5});
    data.monL = flip(data.raw{8});
    data.teyL = flip(data.raw{scanType + 9});
        
    %Manipulate/normalize data before plotting
    data.leftXAS = data.teyL./data.monL; %Normalize left polarized
    data.rightXAS = data.teyR./data.monR; %Normalize right polarized
    data.leftXAS = data.leftXAS-data.leftXAS(1);
    data.rightXAS = data.rightXAS-data.rightXAS(1);
            
    %Addition
    data.aveXAS_wBG = (data.rightXAS + data.leftXAS);
    data.shift = data.aveXAS_wBG(1);
    data.aveXAS_wBG = data.aveXAS_wBG - data.shift;
    data.ScaleFactor = max(data.aveXAS_wBG);
    data.aveXAS_wBG = data.aveXAS_wBG/data.ScaleFactor;
    
    %Subtraction
    data.XMCD = data.rightXAS - data.leftXAS;
    data.XMCD = data.XMCD/data.ScaleFactor;
    
    data.rightXAS = (data.rightXAS-data.shift/2)/data.ScaleFactor;
    data.leftXAS = (data.leftXAS-data.shift/2)/data.ScaleFactor;
    
    %%%%%%%%%%%%%%%%%
    
    %Removing background (steps)
    [data.background.y,~,~] = background(data.Energy, data.aveXAS_wBG, StepE(1,:), StepE(2,:), backscale);
    %[data.background.y,~,~] = background(data.Energy, data.aveXAS_wBG, [703 710] );
    data.aveXAS = data.aveXAS_wBG - data.background.y;
    
    %%%%%%%%%%%%%%%%%%%%
                
    %Finds the index for the lower energy bound of M4 - 1315 to 1350
    [~,qLow] = min(abs(data.Energy-StepE(2,1)));
    [~,qHigh] = min(abs(data.Energy-StepE(2,2)));
    
    %Finds the index for the lower energy bound of p (M5) - 1280 - 1315
    [~,pLow] = min(abs(data.Energy-StepE(1,1)));
    [~,pHigh] = min(abs(data.Energy-StepE(1,2)));
    
    %Integrates over p (M5)
    data.pInt = trapz(data.Energy(pLow:pHigh),data.XMCD(pLow:pHigh));
    data.pSInt = data.XMCD(95);
    
    %Integrates over q (M4)
    data.qInt = trapz(data.Energy(qLow:qHigh),data.XMCD(qLow:qHigh));
    
    %Integrates over full range aveXAS
    data.rInt = trapz(data.Energy(pLow:qHigh),data.aveXAS(pLow:qHigh));
    %data.rInt = trapz(data.Energy(pLow:pHigh),data.aveXAS(pLow:pHigh)) + trapz(data.Energy(qLow:qHigh),data.aveXAS(qLow:qHigh));
    data.rSInt = data.aveXAS(95);
    
    %Integrates over just q (M4) aveXAS
    data.rqInt = trapz(data.Energy(qLow:qHigh),data.aveXAS(qLow:qHigh));
    
    %Calculates m_L
    data.m_L = -(4*(data.pInt+data.qInt)/(3*data.rInt))*nh; % (4(p+q)/3r)*nh
    
    %Calculates m_S effective
    data.m_S = -(2*(data.pInt-2*data.qInt)/(data.rInt))*nh; %nh*(p-2q)/r
    
    %Calculates Tz
    data.Tz = -(5/6)*(data.qInt/data.rInt)*nh;
    
    
    %Plotting%
    %%%%%%%%%%
    if skipPlot==0 && (data.B==5)% || data.B==-5)
        xlims=[700 745];
        figure
        %set(figs{n},'name',['Temperature:',' ',num2str(Temp(n)),' K @ ',num2str(Bfield(n)),' T'],'NumberTitle','Off');
        %Plots left and right on subplot 1
        %subplot(3,1,1)
        plot(data.Energy,data.rightXAS,'Color',color(1,:),'LineWidth',3,'DisplayName','\mu_+')
        hold on
        plot(data.Energy,data.leftXAS,'Color',color(2,:),'LineWidth',3,'DisplayName','\mu_-')
        title(['T = ',num2str(data.T),' K and B = ',num2str(data.B),' T'],'FontSize',40);
        leg = legend('show','location','NorthEast');
        set(leg,'FontSize',30);
        legend('boxoff');
        axis([xlims(1) xlims(2) -0.03 .73]);
        hax = gca;
        set(hax,'FontSize',30,'LineWidth',2,'XTick',[xlims(1):10:xlims(2)],'YTick',[0:.2:1.4],'Linewidth',4);
        hax.XAxis.MinorTick='on';
        hax.YAxis.MinorTick='on';
        hax.XAxis.MinorTickValues=[xlims(1):5:xlims(2)];
        hax.YAxis.MinorTickValues=[-2:.1:2];
        ylabel('XAS (arb. units)','FontSize',50);
        xlabel('Energy (eV)','FontSize',50);
        axis square
        
        %Plots the different of L-R in subplot 2
        %subplot(3,1,2)
        figure
        plot(data.Energy,data.XMCD,'Color',color(1,:),'LineWidth',3,'DisplayName','\mu_+ - \mu_-')
        %title([sampleType{n} ' on Graphene at T = ',num2str(Temp(n)),' K and B = ',num2str(Bfield(n)),' T'],'FontSize',40);
        leg = legend('show','location','NorthEast');
        set(leg,'FontSize',30);
        legend('boxoff');
        axis([xlims(1) xlims(2) -.05 .05]);
        hax = gca;
        set(hax,'FontSize',30,'LineWidth',2,'XTick',[xlims(1):10:xlims(2)],'YTick',[-1:.1:1],'Linewidth',4);
        hax.XAxis.MinorTick='on';
        hax.YAxis.MinorTick='on';
        hax.XAxis.MinorTickValues=[xlims(1):5:xlims(2)];
        hax.YAxis.MinorTickValues=[-2:.01:2];
        ylabel('XMCD (arb. units)','FontSize',50);
        xlabel('Energy (eV)','FontSize',50);
        axis square
        
        %Plots the mean (L+R)/2 in subplot 3
        %subplot(3,1,3)
        figure
        plot(data.Energy,data.aveXAS_wBG,'Color',color(1,:),'LineWidth',3,'DisplayName','\mu_+ + \mu_-')
        hold on
        plot(data.Energy,data.background.y,'Color',color(2,:),'LineWidth',3,'DisplayName','Background Line');
        %title([sampleType{n} ' on Graphene at T = ',num2str(Temp(n)),' K and B = ',num2str(Bfield(n)),' T'],'FontSize',40);
        leg = legend('show','location','NorthEast');
        set(leg,'FontSize',30);
        legend('boxoff');
        axis([xlims(1) xlims(2) -0.03 1.04]);
        hax = gca;
        set(hax,'FontSize',30,'LineWidth',2,'XTick',[xlims(1):10:xlims(2)],'YTick',[0:.2:1.2],'Linewidth',4);
        hax.XAxis.MinorTick='on';
        hax.YAxis.MinorTick='on';
        hax.XAxis.MinorTickValues=[xlims(1):5:xlims(2)];
        hax.YAxis.MinorTickValues=[0:.1:1.2];
        ylabel('Total XAS (arb. units)','FontSize',50);
        xlabel('Energy (eV)','FontSize',50);
        axis square
        
        if printplot == 1
            
            figure
            plot(data.Energy,data.rightXAS,'Color',color(1,:),'LineWidth',1.2,'DisplayName','\mu_+')
            hold on
            plot(data.Energy,data.leftXAS,'Color',color(2,:),'LineWidth',1.2,'DisplayName','\mu_-')
            hax = gca;
            set(hax,'FontSize',8,'LineWidth',1,'XLim',[xlims(1) xlims(2)],'XTick',[xlims(1):10:xlims(2)],'YLim',[-0.03, .63],'YTick',[0:.2:1.4]);
            hax.XAxis.MinorTick='on';
            hax.YAxis.MinorTick='on';
            hax.XAxis.MinorTickValues=[xlims(1):5:xlims(2)];
            hax.YAxis.MinorTickValues=[-2:.05:2];
            leg = legend('show','location','NorthEast');
            set(leg,'FontSize',8);
            legend('boxoff');
            ylabel('XAS (arb. units)','FontSize',13);
            xlabel('Energy (eV)','FontSize',13);
            axis square
            fig = gcf;
            set(gcf, 'color', 'none','inverthardcopy', 'off');
            set(gca, 'xcolor', [0 0 0],'ycolor', [0 0 0],'color', 'none');
            fig.Units = 'centimeters';
            fig.Position = [0 0 8.58924 8.58924];
            pbaspect([4 4 1]);
            %printeps(fig.Number,'FeXAS');


            figure
            plot(data.Energy,data.XMCD,'Color',color(1,:),'LineWidth',1.2,'DisplayName','\mu_+ - \mu_-')
            hax = gca;
            set(hax,'FontSize',8,'LineWidth',1,'XLim',[xlims(1) xlims(2)],'XTick',[xlims(1):10:xlims(2)],'YLim',[-.063 .023],'YTick',[-1:.02:1]);
            hax.XAxis.MinorTick='on';
            hax.YAxis.MinorTick='on';
            hax.XAxis.MinorTickValues=[xlims(1):5:xlims(2)];
            hax.YAxis.MinorTickValues=[-2:.01:2];
            leg = legend('show','location','NorthEast');
            set(leg,'FontSize',8);
            legend('boxoff');
            ylabel('XAS (arb. units)','FontSize',13);
            xlabel('Energy (eV)','FontSize',13);
            axis square
            fig = gcf;
            set(gcf, 'color', 'none','inverthardcopy', 'off');
            set(gca, 'xcolor', [0 0 0],'ycolor', [0 0 0],'color', 'none');
            fig.Units = 'centimeters';
            fig.Position = [0 0 8.93 8.93];
            pbaspect([4 4 1]);
            %printeps(fig.Number,'FeXMCD');


            figure
            plot(data.Energy,data.aveXAS_wBG,'Color',color(1,:),'LineWidth',1.2,'DisplayName','\mu_+ + \mu_-')
            hold on
            plot(data.Energy,data.background.y,'Color',color(2,:),'LineWidth',1.2,'DisplayName','Background Line');
            hax = gca;
            set(hax,'FontSize',8,'LineWidth',1,'XLim',[xlims(1) xlims(2)],'XTick',[xlims(1):10:xlims(2)],'YLim',[-0.03, 1.04],'YTick',[0:.2:1.2]);
            hax.XAxis.MinorTick='on';
            hax.YAxis.MinorTick='on';
            hax.XAxis.MinorTickValues=[xlims(1):5:xlims(2)];
            hax.YAxis.MinorTickValues=[0:.1:1.2];
            leg = legend('show','location','NorthEast');
            set(leg,'FontSize',8);
            legend('boxoff');
            ylabel('XAS (arb. units)','FontSize',13);
            xlabel('Energy (eV)','FontSize',13);
            axis square
            fig = gcf;
            set(gcf, 'color', 'none','inverthardcopy', 'off');
            set(gca, 'xcolor', [0 0 0],'ycolor', [0 0 0],'color', 'none');
            fig.Units = 'centimeters';
            fig.Position = [0 0 8.58924 8.58924];
            pbaspect([4 4 1]);
            %printeps(fig.Number,'FeTotalXAS');
            
        end
        
    end
    
    %Saves data into filedata
    filedata{n}=data;
    

end %Ending loop for multiple files

%Collects the data from all the files
%Preallocates matrix space
Bfield=zeros(1,fileLength);
m_L=zeros(1,fileLength);
m_S=zeros(1,fileLength);
Tz=zeros(1,fileLength);
pSInt = Tz;
rSInt = Tz;
for n=1:fileLength
    Bfield(n) = filedata{n}.B;
    Temp(n) = filedata{n}.T;
    m_L(n) = filedata{n}.m_L;
    m_S(n) = filedata{n}.m_S;
    pSInt(n) = filedata{n}.pSInt;
    rSInt(n) = filedata{n}.rSInt;
end


%Creates figure and plots m_S vs B
figure%figs{length(fileName)+1} = figure(length(fileName)+1);
set(figs{length(fileName)+1},'name',strcat('Moment vs B @ ',num2str(Temp(1)),' K'),'NumberTitle','Off')
plot(Bfield,m_S,'Linewidth',4,'DisplayName','\langleS_Z^{eff}\rangle Experimental','Marker','o','MarkerSize',12,'Color',color(1,:),'MarkerFaceColor',color(1,:));
hold on
plot(Bfield,m_L,'Linewidth',4,'DisplayName','\langleL_Z\rangle Experimental','Marker','s','MarkerSize',12,'Color',color(3,:),'MarkerFaceColor',color(3,:));
%plot(Bfield,m_S_calc,'Linewidth',2,'DisplayName','m_S Calculated','Marker','o','MarkerSize',10,'Color',color(2,:),'MarkerFaceColor',color(2,:));
%plot(Bfield,m_tot,'Linewidth',2,'DisplayName','m_{total} Calculated','Marker','o','MarkerSize',10,'Color',color(4,:),'MarkerFaceColor',color(4,:));
hax = gca;
set(hax,'XLim',[-5.5 5.5],'XTick',[-6:2:6],'YLim',[-.6 .6],'YTick',[-2:.1:2],'FontSize',30,'LineWidth',4);
hax.XAxis.MinorTick='on';
hax.YAxis.MinorTick='on';
hax.XAxis.MinorTickValues=[-6:1:6];
hax.YAxis.MinorTickValues=[-6:1:6];
x = [-7:.01:-.01 .01:.01:7];
y = Brillouin(x,Temp(1),J,g);
plot(x,y,'LineWidth',2,'DisplayName','Brillouin','Marker','none','Color','k');
legID = legend('show','location','NorthWest');
set(legID,'FontSize',30);
%title([sampleType{1} ' at T = ',num2str(Temp(1)),' K'],'FontSize',40);
hline=line([-6 6],[0 0],'LineWidth',2,'Color','k','LineStyle','--');
uistack(hline,'bottom');
xlabel('Magnetic Field (T)','FontSize',50);
ylabel('Magnetic Moments (\mu_B)','FontSize',50);


% figure;
% for n=1:fileLength
%     plot(filedata{n}.Energy,filedata{n}.aveXAS_wBG+n*1.5-1,'DisplayName',strcat(sampleType{n},' @T=',num2str(Temp(n)),' &B=',num2str(Bfield(n))),'linewidth',2);
%     hold on
%     
% end
% %+leg = legend('show','location','NorthEast');
% set(gca,'fontsize',18)
% 
% figure;
% for n=1:fileLength
%     plot(filedata{n}.Energy,filedata{n}.XMCD+(n-1)*.6,'DisplayName',strcat(sampleType{n},' @T=',num2str(Temp(n)),' &B=',num2str(Bfield(n))),'linewidth',2)
%     hold on
%     
% end
% %leg = legend('show','location','NorthEast');
% set(gca,'fontsize',18)