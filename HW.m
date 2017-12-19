%% ANDREA PITTARO
% homework finale di Neuroingegneria - gruppo 9

% Lo script è suddiviso in sezioni, una per ogni richiesta di conti e, nel
% caso vi fossero parti di codice importanti che venivano ripetute, è stata
% creata una funzione. Si è adottato il parfor invece che il for ove
% possibile per far eseguire i calcoli a tutti i cor di cui il processore è
% dotato in parallelo, riducendo il tempo di calcolo di circa un 30%

close all; clear all;
%% DEFINIZIONE COSTANTI e caricamento nomi delle cartelle
FOLDER_NAME = 'A2MB';
data =load_data(FOLDER_NAME); %1 patient per struct and each struct has the data
TR=2.6; %s
THRES=0.05; %sogliatura per i p-value
numPazienti = size(data,2);
numROI  = size(data(1).ROI,2);
numSamplRoi = size(data(1).ROI(1).tac,1);
ThresWM = 50; % 50% della varianza spiegata
ThresCSF = 70;% 70% della varianza spiegata
%% Load alternativo
load HW9_data.mat
TR=2.6; %s
THRES=0.05; %sogliatura per i p-value
numPazienti = size(data,2);
numROI  = size(data(1).ROI,2);
numSamplRoi = size(data(1).ROI(1).tac,1);
ThresWM = 50; % 50% della varianza spiegata
ThresCSF = 70;% 70% della varianza spiegata
%struct con i risultati
ris(numPazienti) = struct('beta',0,'derivaM',0,'derivaVar',0,'ROIfilt',0);

%% task 2 - deriva temporale -- da sistemare

% filtraggio segnali
varianza_thres=5;
for paziente =1:1:numPazienti %per ogni paziente, filtra le ROI
    ris(paziente).ROIfilt = dataFilter(data(paziente).ROI,'butter'); %check errors
%     [ris(paziente).derivaM,ris(paziente).derivaVar] = analisiDeriva(ris(paziente).ROIfilt);
%     if ris(paziente).derivaVar > varianza_thres
%         disp('non vi è deriva temporale lineare')
%     end
end

%% filtro manu Ws=0.004/(Fs/2);
Fs=1/TR;
Ws=0.004/(Fs/2);
Wp=0.007/(Fs/2);
Rs=0.05;
Rp=0.95;
Rp_db=-20*log10(Rp);
Rs_db=-20*log10(Rs);
[n,Wn]=buttord(Wp,Ws,Rp_db,Rs_db);
[b,a]=butter(n,Wn,'high');
N=512;
[H,F] = freqz(b,a,N,Fs);
plot(F,abs(H))
title('Risposta in frequenza del filtro')
xlabel('Frequenza [Hz]')
ylabel('Modulo');

%FILTRAGGIO DATI
filtromanu = zeros(numROI, numPazienti,225);
for j=1:numPazienti
    for i=1:numROI
        filtromanu(i,j,:)=filtfilt(b,a,data(j).ROI(i).tac);
    end
end

%% test filtri
figure(2)
hold on 
plot(data(15).ROI(2).tac )

plot(data(15).ROI(2).tac_filtered )

plot(squeeze(filtromanu(2,15,:) ))
plot(ris(paziente).ROIfilt(2,:))
legend('orig data','filtered given','manu','mio')

%% task 3
%indici per la spiegazione della varianza
explVar_idxCSF = zeros(numPazienti,1);
explVar_idxWM = zeros(numPazienti,1);
numRowCSF=size(data(1).CSF,1);
parfor paziente = 1:1:numPazienti %sceglie il minimo indice 
    explVar_idxWM(paziente) = find(data(paziente).explVarWM>ThresWM,1,'first');
    explVar_idxCSF(paziente) = find(data(paziente).explVarCSF>ThresCSF,1,'first');
end

explVar_idx = round([mean(explVar_idxCSF),mean(explVar_idxWM)]);

for paziente = 1:1:numPazienti
    motion_diff = [diff(data(paziente).motion,1,1) ; zeros(1,size(data(paziente).motion,2))];
     X=[data(paziente).motion, motion_diff, data(paziente).CSF(:,1:explVar_idx(1)),...
         data(paziente).WM(:,1:explVar_idx(2))]; 
     ris(paziente).beta = zeros(35,numROI);  
   for acq=1:1:numROI
       Y=ris(paziente).ROIfilt(acq,:)';
       ris(paziente).beta(:,acq)= (X'*X)\X'*Y;
       model = X*ris(paziente).beta(:,acq);
       data(paziente).regrFilt(acq,:)= ris(paziente).ROIfilt(acq,:) - model';
   end  
end



%% Task 4

%A2MB20 è il paziente 15
tmpPaziente = 15;
ROI_OP=60;
figure(4)
subplot(1,2,1)
plot(data(tmpPaziente).ROI(ROI_OP).tac)
title('Segnale regredito non filtrato')

subplot(1,2,2)
plot(data(tmpPaziente).regrFilt(ROI_OP,:))
title('Segnale regredito filtrato')
hold on 
plot(data(tmpPaziente).ROI(ROI_OP).tac)

figure(3)
hold on
plot(data(tmpPaziente).regrFilt(ROI_OP,:))
plot(data(tmpPaziente).ROI(ROI_OP).tac)
plot(ris(tmpPaziente).ROIfilt(ROI_OP,:))

legend('filtr&regr','orig','filtr NO regr')
hold off

%% task 5
%Pearson correlation
pearsCorr(numPazienti) = struct('FC',0,'FC_parz',0,...
            'signif',0,'signifParz',0,'sogliatura',0);
parfor paziente =1:1:numPazienti
    [FC_paz, signifIDpaz] = corr(data(paziente).regrFilt');
    [FC_paz_part,signifIDpaz_part] = partialcorr(data(paziente).regrFilt');
    pearsCorr(paziente) = struct('FC',FC_paz,'FC_parz',FC_paz_part,...
        'signif',signifIDpaz,'signifParz',signifIDpaz_part,'sogliatura',0);

    
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
alpha = zeros(numPazienti,1);
parfor paziente=1:1:numPazienti
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
alpha0=0.05;
alpha = zeros(numPazienti,1);
pValueSize = size(pearsCorr(1).signif,1);
Q_s = cell(numPazienti);
FDR = cell(numPazienti);
%FDR perché meno restrittivo del permutation test
for paziente=1:1:numPazienti
%     %matrice simmetrica, estraggo una matrice triangolare superiore
    upTri = triu(pearsCorr(paziente).signif) - eye(pValueSize);
    curAnalisys=upTri(upTri(:)>0);
    if isempty(curAnalisys)
        break;
    end
    [FDR{paziente},Q_s{paziente}] = mafdr(curAnalisys);

end

%% task 7
figure(7)
suptitle('Matrici FC per il soggetto 20')
tmpPaziente=12;

subplot(2,2,1)
imagesc(pearsCorr(tmpPaziente).FC); colormap jet; colorbar
title('matrice FC corr')

subplot(2,2,2)
imagesc(pearsCorr(tmpPaziente).signif); colormap jet; colorbar
title('matrice p-values corr')

subplot(2,2,3)
imagesc(pearsCorr(tmpPaziente).FC_parz); colormap jet; colorbar
title('matrice FC corr parziale')

subplot(2,2,4)
imagesc(pearsCorr(tmpPaziente).signifParz); colormap jet; colorbar
title('matrice p-values  corr parziale')

%% task 8
FC_gruppo = zeros(size(pearsCorr(1).FC));
for paziente = 1:1:numPazienti
    FC_gruppo = FC_gruppo + pearsCorr(paziente).FC;
end
FC_gruppo = FC_gruppo/numPazienti;
figure(8)
imagesc(FC_gruppo); colormap jet; colorbar
%% task 9
parfor paz=1:1:numPazienti
    posCorr =  pearsCorr(paz).signif(pearsCorr(paz).signif>=0);
    [pearsCorr(paz).dist]=distance_wei(posCorr);
%     [pearsCorr(paz).dist_parz]=abs(distance_wei(posCorr));
    pearsCorr(paz).ge = efficiency_wei(pearsCorr(paz).dist);
%     pearsCorr(paz).ge_parz = efficiency_wei(pearsCorr(paz).FC_parz);
    [pearsCorr(paz).lambda,pearsCorr(paz).efficiency] = charpath(pearsCorr(paz).dist);
%     [pearsCorr(paz).lambda_parz,pearsCorr(paz).efficiency_parz] = charpath(pearsCorr(paz).dist_parz);
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



