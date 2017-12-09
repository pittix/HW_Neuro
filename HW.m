close all; clear all;
%% DEFINIZIONE COSTANTI e caricamento nomi delle cartelle
FOLDER_NAME = 'A2MB';
data=load_data(FOLDER_NAME); %1 patient per struct and each struct has the data
TR=2.6; %s
A=data(1).ROI;

filtered = dataFilter(data(1).ROI); %check errors


%% task 3
%derivata motion e padding
motion_diff = [diff(data(1).motion,1,1) ; zeros(1,size(data(1).motion,2))];
idxCSF = find(data(1).explVarCSF>70);
idxWM = find(data(1).explVarWM>50);

%matrice per regressione 225*24
%        1                                  6           6   
X= [ones(size(diff(data(1).motion,1),1)),data(1).motion,motion_diff, data(1).CSF(:,[1,idxCSF(1)]),data(1).WM(:,[1,idxWM(1)])];

Y= 0 ;

beta;

%spiegazione sogg1 della varianza: ~~ 4 per white; ~~ per CSF


%% PUNTO 5
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



%%
Fs = 0.38461538462;  % Sampling Frequency
Fstop = 0.005;       % Stopband Frequency
Fpass = 0.006;       % Passband Frequency
Astop = 80;          % Stopband Attenuation (dB)
Apass = 1;           % Passband Ripple (dB)
match = 'passband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, Fs);
Hd = design(h, 'butter', 'MatchExactly', match);

filt = filter(Hd,data(1).ROI(1).tac);


%% OK - stesso prof
Fs = 0.38461538462;  % Sampling Frequency

Fstop = 0.005;   % Stopband Frequency
Fpass = 0.006;   % Passband Frequency
Astop = 7;      % Stopband Attenuation (dB)
Apass = 1;       % Passband Ripple (dB)
match = 'both';  % Band to match exactly

% Construct an FDESIGN object and call its ELLIP method.
h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, Fs);
Hd = design(h, 'ellip', 'MatchExactly', match);

filt = filter(Hd,data(3).ROI(2).tac);

figure(1)
plot(filt)
hold on
plot(data(3).ROI(2).tac_filtered)
hold off
%% 

% %% plot all series of all data translation and rotation
% for i=1:1:numFolders
%     subplot(2,3,1)
%     plot(motionRegressors{i}(:,1))
%     subplot(2,3,2)
%     plot(motionRegressors{i}(:,2))
%     subplot(2,3,3)
%     plot(motionRegressors{i}(:,3))
%     subplot(2,3,4)
%     plot(motionRegressors{i}(:,4))
%     subplot(2,3,5)
%     plot(motionRegressors{i}(:,5))
%     subplot(2,3,6)
%     plot(motionRegressors{i}(:,6))
%     pause(1)
% 
% %     for j=1:1:size(motionRegressors{1},2)
% %         
% %     end
% end