close all; clear all;
%% DEFINIZIONE COSTANTI e caricamento nomi delle cartelle
FOLDER_NAME = 'A2MB';
data =load_data(FOLDER_NAME); %1 patient per struct and each struct has the data
TR=2.6; %s
THRES=0.05; %sogliatura per i p-value
numPazienti = size(data,2);
numROI  = size(data(1).ROI,2);
numSamplRoi = size(data(1).ROI(1).tac,1);

%% Load alternativo
load HW9_data.mat
TR=2.6; %s
THRES=0.05; %sogliatura per i p-value
numPazienti = size(data,2);
numROI  = size(data(1).ROI,2);
numSamplRoi = size(data(1).ROI(1).tac,1);

%% task 2 - deriva temporale -- da sistemare

% filtraggio segnali
ROIfiltrate = zeros(numPazienti, size(data(1).ROI,2),size(data(1).ROI(1).tac,1));
varianza_thres=5;
for paziente =1:1:numPazienti %per ogni paziente, filtra le ROI
    ROIfiltrate(paziente,:,:) = dataFilter(data(paziente).ROI,'butter'); %check errors
    [media,varianza] = analisiDeriva(squeeze(ROIfiltrate(paziente,:,:)));
    if varianza > varianza_thres
        disp('non vi è deriva temporale lineare')
    end
end




%% task 3
%derivata motion e padding
motion_diff = [diff(data(1).motion,1,1) ; zeros(1,size(data(1).motion,2))];
idxCSF = find(data(1).explVarCSF>70,1,'first');
idxWM = find(data(1).explVarWM>50,1,'first');

%spiegazione sogg1 della varianza: ~~ 4 per white; ~~ per CSF
%indici per la spiegazione della varianza
explVarWM_idx = zeros(numPazienti,1);
explVarCSF_idx = zeros(numPazienti,1);

for paziente =1:1:numPazienti
    explVarWM_idx(paziente) = find(data(paziente).explVarWM,1,'first');
    explVarCSF_idx(paziente) = find(data(paziente).explVarCSF,1,'first');
end
%matrice per regressione 225*24
%        1                                  6           6   
% X= [ones(size(diff(data(1).motion,1),1)),data(1).motion,motion_diff, data(1).CSF(:,[1,idxCSF(1)]),data(1).WM(:,[1,idxWM(1)])];
betas = cell(numPazienti,1);
risultati = cell(numPazienti,numROI);
for paziente = 1:1:numPazienti
    motion_diff = [diff(data(paziente).motion,1,1) ; zeros(1,size(data(paziente).motion,2))];
   X=[data(paziente).motion, motion_diff, data(paziente).CSF(:,[1,1:explVarCSF_idx(paziente)]), ...
       data(paziente).WM(:,[1,1:explVarWM_idx(paziente)])]; 
   for acq=1:1:numROI
       Y=data(paziente).ROI(acq).tac_filtered;
       betas{paziente} =  (X'*X)\X'*Y;
       data(paziente).regr(acq,:)= data(paziente).ROI(acq).tac - X*betas{paziente};
       data(paziente).regrFilt(acq,:)= data(paziente).ROI(acq).tac_filtered - X*betas{paziente};
   end
   
end

%% Task 4

%A2MB20 è il paziente 15
tmpPaziente = 15;
ROI_OP=60;
subplot(1,2,1)
 plot(data(tmpPaziente).regr(ROI_OP,:))
title('Segnale regredito filtrato')

subplot(1,2,2)
plot(data(tmpPaziente).regrFilt(ROI_OP,:))
title('Segnale regredito non filtrato')


%% task 5
%Pearson correlation
pearsCorr(numPazienti) = struct('FC',0,'FC_parz',0,...
            'signif',0,'signifParz',0);
parfor paziente =1:1:numPazienti
    curROI = zeros(numROI,numSamplRoi);
    for roi =1:1:numROI
        curROI(roi,:)=data(paziente).ROI(roi).tac_filtered;
    end
    [FC_paz, signifIDpaz] = corr(curROI);
    [FC_paz_part,signifIDpaz_part] = partialcorr(curROI,curROI);
    pearsCorr(paziente) = struct('FC',FC_paz,'FC_parz',FC_paz_part,...
        'signif',signifIDpaz,'signifParz',signifIDpaz_part);

    
end
clear signifIDpaz signifIDpaz_part FC_paz FC_paz_part

figure(5)
for paziente =1:1:numPazienti
   subplot(1,2,1)
   imagesc(pearsCorr(paziente).FC);colormap jet;colorbar
   subplot(1,2,2)
   imagesc(pearsCorr(paziente).FC_parz);colormap jet;colorbar
   pause(0.2)
end
%% task6 - copiato
alpha0=0.05;
for paziente=1:numPazienti
    triangl=triu(pearsCorr(paziente).signif);
    vett_pval=triangl(:);
    pos_pval=vett_pval>0;
    vett_pval=vett_pval(pos_pval==1);
    vett_pval=sort(unique(vett_pval));
    j=1;
    temp=(j*alpha0)/length(vett_pval);
    while vett_pval(j)<temp
        j=j+1;
        temp=(j*alpha0)/length(vett_pval);
    end
    alpha(paziente)=vett_pval(j);
end

%% task 6 - sogliatura hard 
for paziente=1:1:numPazienti
    pearsCorr(paziente).sogliatura = pearsCorr(paziente).signif.*alpha;
 
end

%% task 7
figure(7)
subplot(2,2,1)
imagesc(pearsCorr(tmpPaziente).FC)
colormap jet; colorbar
subplot(2,2,2)
imagesc(pearsCorr(tmpPaziente).signif)
colormap jet; colorbar

subplot(2,2,3)
imagesc(pearsCorr(tmpPaziente).FC_parz)
colormap jet; colorbar

subplot(2,2,4)
plot(pearsCorr(tmpPaziente).signifParz)
colormap jet; colorbar


%% task 8
FC_gruppo = zeros(size(pearsCorr(1).FC));
for paziente = 1:1:numPazienti
    FC_gruppo = FC_gruppo + pearsCorr(paziente).FC;
end
FC_gruppo = FC_gruppo/numPazienti;
figure(8)
imagesc(FC_gruppo); colormap jet; colorbar
%% task 9
for paz=1:1:numPazienti
    [pearsCorr(paz).dist]=abs(distance_wei(inv(pearsCorr(tmpPaziente).FC)));
    [pearsCorr(paz).dist_parz]=abs(distance_wei(inv(pearsCorr(tmpPaziente).FC_parz)));
    pearsCorr(paz).ge = efficiency_wei(dist.*pearsCorr(tmpPaziente).FC);
    pearsCorr(paz).ge_parz = efficiency_wei(dist_parz.*pearsCorr(tmpPaziente).FC_parz);
    [pearsCorr(paz).lambda,pearsCorr(paz).efficiency] = charpath(pearsCorr(paz).dist);
    [pearsCorr(paz).lambda_parz,pearsCorr(paz).efficiency_parz] = charpath(pearsCorr(paz).dist_parz);
end

figure(9)
imagesc(pearsCorr(tmpPaziente).dist); colormap jet ; colorbar
%% task 10 



%% 2 - plot deriva nel tempo
figure(1)
for i=1:1:numPazienti
    for j=1:1:numROI
        curROI =data(i).ROI(j).tac; 
        plot(curROI)
        pause(0.1)
        [mean, var] = analisiDeriva(curROI);
    end
end



