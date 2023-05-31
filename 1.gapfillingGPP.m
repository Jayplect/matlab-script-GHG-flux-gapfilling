addpath 'C:\Users\grad_student\OneDrive\Documents\Matlab files\Elkle'

% MAIN GROWING SEASON
% load dataframe-'dataGS' from 'finalfolder' 
% Make sure each column is formated to number
Kdown_all = dataGS(:,[2]);
Kdown_all = table2array (Kdown_all);
GPP_all = dataGS(:,[4]);
GPP_all = table2array (GPP_all);
dframe = rmmissing(dataGS,'DataVariables',{'GPPmod'});
Kdown = dframe(:,[2]);
Kdown = table2array (Kdown);
GPPmeas = dframe(:,[4]);
GPPmeas = table2array (GPPmeas);

%if you want to save output, change the plus sign in the fit of gapfill_GPP
%to negative before invoking next line. Make sure the next line are
%+vegppmeas and gpp_all
[GPP_GS1] = gapfill_GPP(GPPmeas,Kdown,GPP_all,Kdown_all);
gsKdown_all = Kdown_all;
gsGPP_all = GPP_all;

%fit the curve and save in photofitting for plotting
photosynthesis = fittype('(a*Px*Kdown)/(a*Kdown-Px)','independent',{'Kdown'},'dependent',{'GPP'});
% define start points for the parameters a and Px 
startpar = [0.05 101];
photofitting = fit(Kdown,+GPPmeas,photosynthesis,'Startpoint',startpar);

%plot curve: only use this for P12
hold off
p3=plot(gsKdown_all,-gsGPP_all, '.r',"color","green",'MarkerSize',40);
hold on

plot(photofitting)
set(p3,'linewidth',2);
legend('Location','northwest')
legend('Box','off')
legend('2018/2019 Growing season fluxes','Fitted curve')
ylabel('GEP (\mumol m^-^2 s^-^1)');
xlabel('Incoming shortwave radiation (W m^-^2)');
set(gca,'FontSize',18);

%set x and y limits %inspect graph before applying limits
xlim([0 1020]);
xticks(0:200:1020);

ylim([-12 70]);
yticks(-15:15:70)


% load dataframe-'dataCC' from 'finalfolder'
% CC growing season GPP
Kdown_all = dataCC(:,[2]);
Kdown_all = table2array (Kdown_all);
GPP_all = dataCC(:,[4]);
GPP_all = table2array (GPP_all);
dframe = rmmissing(dataCC,'DataVariables',{'GPPmod'});
Kdown = dframe(:,[2]);
Kdown = table2array (Kdown);
GPPmeas = dframe(:,[4]);
GPPmeas = table2array (GPPmeas);

%if you want to save output, change the plus sign in the fit of gapfill_GPP
%to negative before invoking next line
[GPP_CC] = gapfillGppCC(GPPmeas,Kdown,GPP_all,Kdown_all);
ccKdown_all = Kdown_all;
ccGPP_all = GPP_all;

%fit the curve and save in fitting for plotting
photosynthesis1 = fittype('(a*Px*Kdown)/(a*Kdown+Px)','independent',{'Kdown'},'dependent',{'GPP'});
% define start points for the parameters a and Px 
startpar = [0.05 101];
photofitting1 = fit(Kdown,-GPPmeas,photosynthesis1,'Startpoint',startpar);


%plot curve Use this chunk below for P3
hold off
p3=plot(gsKdown_all,-gsGPP_all, '.r',"color","green",'MarkerSize',40)
hold on
p8 = plot(ccKdown_all,-ccGPP_all, '.r',"color","red",'MarkerSize',40)
h8 = get(gca,'Children')
legendstr8={'2019/2020 Growing period fluxes', '2019/2020 Post-growing period fluxes'};
legend(h8([2 1]),legendstr8{[1 2]})
plot(photofitting)
plot(photofitting1)
legend('Location','northwest')
legend('Box','off')
set(gca,'FontSize',20);
ylabel('GEP (\mumol m^-^2 s^-^1)');
xlabel('Incoming shortwave radiation (W m^-^2)');

%set x and y limits %inspect graph before applying limits
xlim([0 1020]);
xticks(0:200:1020);

ylim([-12 70]);
yticks(-15:15:70)
