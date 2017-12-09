function [data] = load_data(folder)
%Carica I dati e li ritorna nella struct data
%Carica I dati e li ritorna nelle celle di variabili

TIMESERIES_NAME='atlas_sampled_timeseries';
MOT_EST_NAME='estimated_motion';
folder_list = dir(folder);
%Il -2 è per escludere la cartella corrente (A2MB) e la cartella precedente
%(sotto Linux)
numFolders = length(folder_list)-2; %numero di soggetti. 

% caricamento variabili
% CSF_princ = cell(numFolders,1);
% motionRegressors = cell(numFolders,1);
% WM_princ = cell(numFolders,1);
% explVarCSF = cell(numFolders,1);
% explVarWM = cell(numFolders,1);
% malf_ROI = cell(numFolders,1);

% data = zeros(numFolders,1);
for i=1:1:numFolders
    folder_list_path = strcat(folder_list(i+2).folder, '/', folder_list(i+2).name, '/');   
    load(strcat(folder_list_path,TIMESERIES_NAME));
    load(strcat(folder_list_path,MOT_EST_NAME));
%inutile se usiamo la struct
%     malf_ROI{i}         = malfROI;
%     CSF_princ{i}        = CSFprinc;
%     WM_princ{i}         = WMprinc;
%     explVarCSF{i}       = explained_variance_CSF;
%     explVarWM{i}        = explained_variance_WM;
%     motionRegressors{i} = Motion_regressors;
    
    dataTmp = struct('ROI',malfROI, 'CSF',CSFprinc, 'WM', WMprinc,...
        'explVarCSF',explained_variance_CSF, 'explVarWM',explained_variance_WM,...
        'motion',Motion_regressors);
    data(i) = dataTmp;
    
    clear malfROI CSFprinc WMprinc explained_variance_CSF explained_variance_WM Motion_regressors
end

% data = struct('ROI','malf_ROI', 'CSF','CSF_princ', 'WM', 'WM_princ');
% data.add

end

