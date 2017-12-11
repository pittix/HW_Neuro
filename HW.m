close all; clear all;
%% DEFINIZIONE COSTANTI e caricamento nomi delle cartelle
FOLDER_NAME = 'A2MB';
data =load_data(FOLDER_NAME); %1 patient per struct and each struct has the data
TR=2.6; %s
numPazienti = size(data,2);
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

%matrice per regressione 225*24
%        1                                  6           6   
% X= [ones(size(diff(data(1).motion,1),1)),data(1).motion,motion_diff, data(1).CSF(:,[1,idxCSF(1)]),data(1).WM(:,[1,idxWM(1)])];
betas = cell(numPazienti,size(data(1).ROI,1));
for paziente = 1:1:numPazienti
    motion_diff = [diff(data(paziente).motion,1,1) ; zeros(1,size(data(paziente).motion,2))];
   X=[data(paziente).motion, data(paziente).motion.^2, motion_diff, motion_diff.^2]; 
   for acq=1:1:size(data(paziente).ROI,1)
       Y=data(paziente).ROI(acq).tac_filtered;
       betas{paziente,acq} =  (X'*X)\X'*Y;
   end
   
end

beta;

%spiegazione sogg1 della varianza: ~~ 4 per white; ~~ per CSF


%% task 5
%valid volumes?


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



