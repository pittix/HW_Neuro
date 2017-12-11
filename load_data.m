function [data] = load_data(folder)
%Carica I dati e li ritorna nella struct data
%Carica I dati e li ritorna nelle celle di variabili

TIMESERIES_NAME='atlas_sampled_timeseries';
MOT_EST_NAME='estimated_motion';
folder_list = dir(folder);
%Il -2 è per escludere la cartella corrente (A2MB) e la cartella precedente
%(sotto Linux)
numFolders = length(folder_list)-2; %numero di soggetti. 

%inizializzazione struct
data(numFolders) = struct('ROI',0,'CSF',0,'WM',0,'explVarCSF',0,'explVarWM',0,'motion',0);

 for i=1:1:numFolders
    folder_list_path = strcat(folder_list(i+2).folder, '/', folder_list(i+2).name, '/');   
    load(strcat(folder_list_path,TIMESERIES_NAME));
    load(strcat(folder_list_path,MOT_EST_NAME));
   
    dataTmp = struct('ROI',malfROI, 'CSF',CSFprinc, 'WM', WMprinc,...
        'explVarCSF',explained_variance_CSF, 'explVarWM',explained_variance_WM,...
        'motion',Motion_regressors);
    data(i) = dataTmp;
    
    clear malfROI CSFprinc WMprinc explained_variance_CSF explained_variance_WM Motion_regressors
end

end

