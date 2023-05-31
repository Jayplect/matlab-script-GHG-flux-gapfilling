addpath 'C:\Users\grad_student\OneDrive\Documents\Matlab files\Elkle'
%non-growing season Resp
%import NEE and Tair as table from combineNGS located in finalfolder 
%combineNGS=xlsread("C:\Users\grad_student\OneDrive\Documents\finalfolder\combineNGS.csv");

R_all = combineNGS(:,[1]); 
R_all = table2array (R_all);
R_all(1,:) = [];
T_all = combineNGS(:,[2]);
T_all = table2array (T_all);
T_all(1,:) = [];
T_all = fillmissing(T_all,'linear');
find(isnan(T_all));

%combineNGS=array2table(combineNGS);
dframe = rmmissing(combineNGS,'DataVariables',{'NEE'});
Rmeas = dframe(:,[1]);
Rmeas = table2array (Rmeas);
Tinput = dframe(:,[2]);
Tinput = table2array (Tinput);
Tinput = fillmissing(Tinput,'linear');
find(isnan(Tinput));

[R_nonGS] = gapfill_resp(Rmeas,Tinput,R_all,T_all);
ngsT = T_all;
ngsR = R_all;
%fit the curve and save in fitting
respiration = fittype('r1./(1+exp(r2.*(r3-Tinput)))','independent',{'Tinput'},'dependent',{'Rmeas'});
startpar = [6.3 0.19 19]
fitting = fit(Tinput,Rmeas,respiration,'Lower',[0,0,0],'Startpoint',startpar)

% growing season
%import NEEnightflux and Tair as table from nighttimeGS located in finalfolder 
%nighttimeGS=xlsread('C:\Users\grad_student\OneDrive\Documents\finalfolder\nighttimeGS.csv');
T_all = nighttimeGS(:,[1]);
T_all = table2array (T_all);
T_all(1,:) = [];
R_all = nighttimeGS(:,[2]);
R_all = table2array (R_all);
R_all(1,:) = [];
T_all = fillmissing(T_all,'linear');
find(isnan(T_all))

dframe = rmmissing(nighttimeGS,'DataVariables',{'NEEnightflux'});
Tinput = dframe(:,[1]);
Tinput = table2array (Tinput);
Rmeas = dframe(:,[2]);
Rmeas = table2array (Rmeas);
Tinput = fillmissing(Tinput,'linear');
find(isnan(Tinput))

[R_GS1] = gapfill_resp(Rmeas,Tinput,R_all,T_all);
gsT = T_all;
gsR = R_all;
%fit the curve and save in fitting
respiration = fittype('r1./(1+exp(r2.*(r3-Tinput)))','independent',{'Tinput'},'dependent',{'Rmeas'});
startpar = [6.3 0.19 19];
fitting1 = fit(Tinput,Rmeas,respiration,'Lower',[0,0,0],'Startpoint',startpar);

% MAIN GROWING SEASON 
%load NEEgrowingseason from dataframe, nighttimeGS, and Rg as COLUMN for entire period under consideration
NEEgrowingseason(1,:) = [];
GPP_GS = NEEgrowingseason - R_GS1 ;
Rg(1,:) = [];
%plot Rg to see if missing values and fill using chunk in next line
Rg = fillmissing(Rg,'linear', 'EndValues','none');
find(isnan(Rg))
tablex = table(Rg,GPP_GS);
writetable(tablex,'finalfolder/dataGS.csv','Delimiter',',')
%next step is to move to R file- "matlab continuation script" and run line
%97 to 103


%plot graph non growing season
hold off
p=plot(ngsT,ngsR, '.r',"color","red",'MarkerSize',40)
hold on
%plot graph growing season
p1=plot(gsT,gsR, '.r',"color","blue",'MarkerSize',40)
h1=get(gca,'Children')

%rename legend according to year
legendstr1={'2019/2020 Growing period fluxes','2019/2020 Non-growing period fluxes'};
legend(h1([2 1]),legendstr1{[2 1]})
plot(fitting)
plot(fitting1)
ylabel('Respiration CO_2 flux (\mumol m^-^2 s^-^1)');
xlabel('Temperature (^oC)');
legend('Location','northwest')
legend('Box','off')
%inspect graph before applying limits
set(gca,'FontSize',18)

xlim([-25 30])
xticks(-30:10:30)

ylim([-10 20])
yticks(-10:10:20)


%CC growing season Resp
%import NEEnightCC and Tair from nighttimeCC located in finalfolder
T_all = nighttimeCC(:,[1]);
T_all = table2array (T_all);
T_all(1,:) = [];
R_all = nighttimeCC(:,[2]);
R_all = table2array (R_all);
R_all(1,:) = []
dframe = rmmissing(nighttimeCC,'DataVariables',{'NEEnightCC'})
Tinput = dframe(:,[1]);
Tinput = table2array (Tinput);
Rmeas = dframe(:,[2]);
Rmeas = table2array (Rmeas);
[R_CC] = gapfillrespCC(Rmeas,Tinput,R_all,T_all);
covercropT = T_all;
covercropR = R_all;
%fit the curve and save in fitting
respiration = fittype('r1./(1+exp(r2.*(r3-Tinput)))','independent',{'Tinput'},'dependent',{'Rmeas'});
startpar = [6.3 0.19 19];
fitting2 = fit(Tinput,Rmeas,respiration,'Lower',[0,0,0],'Startpoint',startpar);

% COVERCROP GROWING PERIODS in fall/spring 
%load NEEcovercrop from dataframe (i.e., nighttimeCC) and Rg for entire period under consideration
NEEcovercrop(1,:) = []
GPP_CC = NEEcovercrop - R_CC 
Rg1(1,:) = []
%plot Rg to see if missing values and fill using chunk in next line
%Rg = fillmissing(Rg,'linear', 'EndValues','none')
tablex = table(Rg1,GPP_CC)
writetable(tablex,'finalfolder/dataCC.csv','Delimiter',',')
%next step is to move to R file- "matlab continuation script" and run line 54
%to 59

%plot graph for growing period, non-growing and post growing period
hold off
p=plot(ngsT,ngsR, '.r',"color","red",'MarkerSize',40)
hold on
p1=plot(gsT,gsR, '.r',"color","blue",'MarkerSize',40)
p2=plot(covercropT,covercropR, '.r',"color","green",'MarkerSize',40)
h=get(gca,'Children')
legendstr={'2020/2021 Growing period fluxes','2020/2021 Non-growing period fluxes', '2020/2021 Post-growing period fluxes'};
legend(h([3 2 1]),legendstr{[2 1 3]})
plot(fitting)
plot(fitting1)
plot(fitting2)
ylabel('Respiration CO_2 flux (\mumol m^-^2 s^-^1)');
xlabel('Temperature (^oC)');
legend('Location','northwest')
legend('Box','off')
%inspect graph before applying limits
set(gca,'FontSize',18)

xlim([-25 30])
xticks(-30:10:30)

ylim([-10 20])
yticks(-10:10:20)

