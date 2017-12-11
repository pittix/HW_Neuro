close all; clear all;
%% DEFINIZIONE COSTANTI e caricamento nomi delle cartelle
FOLDER_NAME = 'A2MB';
data =load_data(FOLDER_NAME); %1 patient per struct and each struct has the data
TR=2.6; %s
numPazienti = size(data,2);
numROI  = size(data(1).ROI,1);
numSamplRoi = size(data(1).ROI(1).tac,1);
%% Load alternativo
load HW9_data.mat


%% task 2 - deriva temporale

% filtraggio segnali
ROIfiltrate = zeros(numPazienti, size(data(1).ROI,2),size(data(1).ROI(1).tac,1));
for paziente =1:1:numPazienti %per ogni paziente, filtra le ROI
    ROIfiltrate(paziente,:,:) = dataFilter(data(paziente).ROI,'butter'); %check errors
    [media,varianza] = analisiDeriva(squeeze(ROIfiltrate(paziente,:,:)));
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
    explVarWM_idx = find(data(paziente).explVarWM,1,'first');
    explVarCSF_idx = find(data(paziente).explVarCSF,1,'first');
end
%matrice per regressione 225*24
%        1                                  6           6   
% X= [ones(size(diff(data(1).motion,1),1)),data(1).motion,motion_diff, data(1).CSF(:,[1,idxCSF(1)]),data(1).WM(:,[1,idxWM(1)])];
betas = cell(numPazienti,1);
for paziente = 1:1:numPazienti
    motion_diff = [diff(data(paziente).motion,1,1) ; zeros(1,size(data(paziente).motion,2))];
   X=[data(paziente).motion, motion_diff, data(paziente).CSF(:,[1,1:explVarCSF_idx(paziente)]), ...
       ata(paziente).WM(:,[1,1:explVarWM_idx(paziente)])]; 
   for acq=1:1:size(data(paziente).ROI,1)
       Y=data(paziente).ROI(acq).tac_filtered;
       betas{paziente} =  (X'*X)\X'*Y;
   end
    % manca stima di parametri dati i beta
end

%% Task 4

%A2MB20 è il paziente 15
tmpPaziente = 15;
ROI_OP=60;
subplot(1,2,1)
 plot(data(tmpPaziente).ROI(ROI_OP).tac)
title('Segnale regredito filtrato')

subplot(1,2,2)
% plot(plot(data(tmpPaziente).ROI(ROI_OP).tac_filtered))
title('Segnale regredito non filtrato')

processed_fMRI = zeros(numROI,numSamplROI,numPazienti);
%devo applicare i beta.....come? :((

%% task 5
%valid volumes?


%% task 6


%% task 7


%% task 8

%% task 9

%% task 10 


%% 2 - plot deriva nel tempo
figure(1)
for i=1:1:numFolders
    curROI =malf_ROI{i}; 
    for j=1:1:size(malf_ROI{1},1)
        plot(curROI(j).tac)
        pause(0.5)
        [mean, var] = analisiDeriva(curROI);
    end
end



