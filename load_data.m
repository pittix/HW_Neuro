function [data] = load_data(folder)
%Carica I dati e li ritorna nella struct data
%Carica I dati e li ritorna nelle celle di variabili

TIMESERIES_NAME='atlas_sampled_timeseries.mat';
MOT_EST_NAME='estimated_motion.mat';
folder_list = dir(folder);
matlVer_str = version('-release');
matlVer = str2num(matlVer_str(1:end-1));
%Il -2 è per escludere la cartella corrente (A2MB) e la cartella precedente
%(sotto Linux)
numFolders = length(folder_list)-2; %numero di soggetti. 

%inizializzazione struct
data(numFolders) = struct('ROI',0,'CSF',0,'WM',0,'explVarCSF',0,'explVarWM',0,'motion',0);

 for i=1:1:numFolders
    %percorso completo per la cartella corrente
    folder_list_path = strcat(pwd,'/',folder,'/',folder_list(i+2).name,'/');  
    %caricamento delle matrici di quella cartella
    load(strcat(folder_list_path,TIMESERIES_NAME));
    load(strcat(folder_list_path,MOT_EST_NAME));
    %creazione della struct per il soggetto corrente
    dataTmp = struct('ROI',malfROI, 'CSF',CSFprinc, 'WM', WMprinc,...
        'explVarCSF',explained_variance_CSF, 'explVarWM',explained_variance_WM,...
        'motion',Motion_regressors);
    data(i) = dataTmp;
    %eliminazione delle matrici appartenenti al soggetto corrente
    clear malfROI CSFprinc WMprinc explained_variance_CSF explained_variance_WM Motion_regressors
end

end

