%% ANDREA PITTARO
% homework finale di Neuroingegneria - gruppo 9

% Lo script è suddiviso in sezioni, una per ogni richiesta di conti e, nel
% caso vi fossero parti di codice importanti che venivano ripetute, è stata
% creata una funzione. Si è adottato il parfor invece che il for ove
% possibile per far eseguire i calcoli a tutti i core di cui il processore è
% dotato in parallelo, riducendo il tempo di calcolo di circa un 30%

close all; clear all;
%% DEFINIZIONE COSTANTI e caricamento nomi delle cartelle
FOLDER_NAME = 'A2MB';
data =load_data(FOLDER_NAME); %1 patient per struct and each struct has the data
TR=2.6; %s
THRES=0.05; %sogliatura per i p-value
numSoggetti = size(data,2);
numROI  = size(data(1).ROI,2);
numSamplRoi = size(data(1).ROI(1).tac,1);
ThresWM = 50; % 50% della varianza spiegata
ThresCSF = 70;% 70% della varianza spiegata
time = TR * (1:numSamplRoi);
%struct dei risultati. A fine file verranno generate le matrici richieste
ris(numSoggetti) = struct('beta',0,'ROIfilt',0,'regrFilt',0,'FC',0, ...
        'FC_parz',0,'signif',0,'signif_parz',0,'sogliatura',0,'CPL_corr',0,...
        'CPL_parz',0,'GE_corr',0,'GE_parz',0,'eff_corr',0,'eff_parz',0,...
        'FCThres',0,'FC_parzThres',0,'signifThres',0,'signifThres_parz',0);
    
sogg_rif = 15; %soggetto A2MB20
ROI_rif  = 60; 

%% task 2 - deriva temporale -- da sistemare

% filtraggio segnali
varianza_thres=5;
drift = zeros(numSoggetti,numROI);
for sogg =1:1:numSoggetti %per ogni paziente, filtra le ROI
    ris(sogg).ROIfilt = dataFilter(data(sogg).ROI,'butter'); %check errors
    [drift(sogg,:)] = analisiDeriva(ris(sogg).ROIfilt);
    drift_medio=mean(abs(drift(:))); %valore assoluto per confrontarli
    drift_varianza=var(abs(drift(:))); 

end

%% test filtri
f2=figure(2);
set(f2,'Name','Confronto dati originali e filtrati','NumberTitle','off',...
    'units','normalized','outerposition',[0 0 1 1])
hold on 
plot(time,data(sogg_rif).ROI(ROI_rif).tac )

plot(time,data(sogg_rif).ROI(ROI_rif).tac_filtered )

plot(time,ris(sogg_rif).ROIfilt(ROI_rif,:))
legend('dati originali','dati filtrati forniti','dati filtrati')
hold off
title('Confronto segnali pre e post filtraggio')
xlabel('tempo')
%% task 3
%indici per la spiegazione della varianza
explVar_idxCSF = zeros(numSoggetti,1);
explVar_idxWM = zeros(numSoggetti,1);
numRowCSF=size(data(1).CSF,1);
parfor sogg = 1:1:numSoggetti %sceglie il minimo indice 
    explVar_idxWM(sogg) = find(data(sogg).explVarWM>ThresWM,1,'first');
    explVar_idxCSF(sogg) = find(data(sogg).explVarCSF>ThresCSF,1,'first');
end

explVar_idx = round([mean(explVar_idxCSF),mean(explVar_idxWM)]);

parfor sogg = 1:1:numSoggetti
    motion_diff = [diff(data(sogg).motion,1,1) ; zeros(1,size(data(sogg).motion,2))];
    X=[data(sogg).motion, motion_diff, data(sogg).CSF(:,1:explVar_idx(1)),...
         data(sogg).WM(:,1:explVar_idx(2))]; 
     ris(sogg).beta = zeros(35,numROI);  
     ris(sogg).regrFilt = zeros(numROI,225);  
   for roi=1:1:numROI
       Y=ris(sogg).ROIfilt(roi,:)';
       ris(sogg).beta(:,roi)= (X'*X)\X'*Y;
       model = X*ris(sogg).beta(:,roi);
       ris(sogg).regrFilt(roi,:) = ris(sogg).ROIfilt(roi,:) - model';
   end  
end



%% Task 4

%A2MB20 è il paziente 15
f3=figure(3);
set(f3,'Name','Segnali filtrati e regrediti','NumberTitle','off',...
    'units','normalized','outerposition',[0 0 1 1])
hold on
plot(ris(sogg_rif).regrFilt(ROI_rif,:))
plot(data(sogg_rif).ROI(ROI_rif).tac)
plot(ris(sogg_rif).ROIfilt(ROI_rif,:))

legend('filtrato & regredito','originale','filtrato NON regredito')
hold off

%% task 5 calcolo correlazione e correlazione parziale
parfor sogg =1:1:numSoggetti
    [ris(sogg).FC, ris(sogg).signif] = corr(ris(sogg).regrFilt');
    [ris(sogg).FC_parz,ris(sogg).signif_parz] = partialcorr(ris(sogg).regrFilt');    
end

f5=figure(5);
set(f5,'Name','Confronto corr e corr parziale','NumberTitle','off',...
    'units','normalized','outerposition',[0 0 1 1])
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
t=text( 0.5, 0,'', 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
for sogg =1:1:numSoggetti
   set(t,'String',['Connettività funzionale soggetto ',num2str(sogg)]);
   subplot(1,2,1)
   imagesc(ris(sogg).FC);colormap jet;colorbar
   title('correlazione')
   subplot(1,2,2)
   imagesc(ris(sogg).FC_parz);colormap jet;colorbar
   title('correlazione parziale')
   pause(0.25)
end
%% task6 
alpha0=0.05;
alpha = zeros(numSoggetti,1);
%correlazione totale
for sogg=1:1:numSoggetti
    tri=triu(ris(sogg).signif);
    vett_pval=tri(:);
    pos_pval=vett_pval>0;
    vett_pval=vett_pval(pos_pval==1);
    vett_pval=sort(unique(vett_pval));
    j=1;
    temp=(j*alpha0)/length(vett_pval);
    while vett_pval(j)<temp 
        j=j+1;
        temp=(j*alpha0)/length(vett_pval);
    end
    alpha(sogg)=vett_pval(j);
    mask=ris(sogg).signif<alpha(sogg);
    ris(sogg).signifThres=mask .* ris(sogg).signif;
    ris(sogg).FCThres=mask .* ris(sogg).FC;
end
alpha_parz = zeros(numSoggetti,1);
%correlazione totale
for sogg=1:1:numSoggetti
    tri=triu(ris(sogg).signif_parz);
    vett_pval=tri(:);
    pos_pval=vett_pval>0;
    vett_pval=vett_pval(pos_pval==1);
    vett_pval=sort(unique(vett_pval));
    j=1;
    temp=(j*alpha0)/length(vett_pval);
    while vett_pval(j)<temp 
        j=j+1;
        temp=(j*alpha0)/length(vett_pval);
    end
    alpha_parz(sogg)=vett_pval(j);
    mask=ris(sogg).signif_parz<alpha(sogg);
    ris(sogg).signifThres_parz=mask .* ris(sogg).signif_parz;
    ris(sogg).FCThres_parz=mask .* ris(sogg).FC_parz;
end

%% task 6 - sogliatura hard - prova 
alpha0=0.05;
alpha = zeros(numSoggetti,1);
pValueSize = size(ris(1).signif,1);
Q_s = cell(numSoggetti);
FDR = cell(numSoggetti);
%FDR perché meno restrittivo del permutation test
for sogg=1:1:numSoggetti
    %matrice simmetrica, estraggo una matrice triangolare superiore
    upTri = triu(ris(sogg).signif) - eye(pValueSize);
    curAnalisys=upTri(upTri(:)>0);
    upTri_parz = triu(ris(sogg).signif_parz) - eye(pValueSize);
    curAnalisys_parz=upTri(upTri_parz(:)>0);
    
    if isempty(curAnalisys)
        break;
    end
    [FDR{sogg},Q_s{sogg}] = mafdr(curAnalisys);
    [FDR_parz{sogg},Q_s_parz{sogg}] = mafdr(curAnalisys);

end

%% task 7
f7=figure(7);
set(f7,'Name','confronto correlazioni e p-values ','NumberTitle','off',...
    'units','normalized','outerposition',[0 0 1 1])
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, 'Matrici FC per il soggetto 20', 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;

subplot(2,2,1)
imagesc(ris(sogg_rif).FC); colormap jet; colorbar
title('matrice FC corr')

subplot(2,2,2)
imagesc(ris(sogg_rif).signif); colormap jet; colorbar
title('matrice p-values corr')

subplot(2,2,3)
imagesc(ris(sogg_rif).FC_parz); colormap jet; colorbar
title('matrice FC corr parziale')

subplot(2,2,4)
imagesc(ris(sogg_rif).signif_parz); colormap jet; colorbar
title('matrice p-values  corr parziale')

%% task 8 media dei dati per ottenere la matrice di connettività di gruppo
FC_gruppo = zeros(size(ris(1).FC));
FC_gruppo_parz = zeros(size(ris(1).FC_parz));
for sogg = 1:1:numSoggetti
    FC_gruppo = FC_gruppo + ris(sogg).FC;
    FC_gruppo_parz = FC_gruppo_parz + ris(sogg).FC_parz;
end
FC_gruppo = FC_gruppo/numSoggetti;
FC_gruppo_parz = FC_gruppo_parz/numSoggetti;

%plot dei risultati
f9 = figure(9);
set(f9,'Name','NumberTitle','off','Connettività funzionale di gruppo',...
    'units','normalized','outerposition',[0 0 1 1])
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, 'Connettività funzionale di gruppo', 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;

subplot(1,2,1)
imagesc(FC_gruppo); colormap jet; colorbar
title('Correlazione')
subplot(1,2,2)
imagesc(FC_gruppo_parz); colormap jet; colorbar
title('Correlazione parziale')
%% task 9
parfor sogg=1:1:numSoggetti
    mask=ris(sogg).signifThres>=0;
    inver_corr=1./(mask.*ris(sogg).FCThres);
    [D_corr, B_corr]=distance_wei(inver_corr);
    [ris(sogg).CPL_corr,ris(sogg).eff_corr,~]=charpath(D_corr,0,0);
    ris(sogg).GE_corr=efficiency_wei(mask.*ris(sogg).FCThres);

    mask_parcorr=ris(sogg).FCThres_parz>=0;
    inver_parcorr=1./(mask_parcorr.*ris(sogg).FCThres_parz);
    [D_parcorr, B_parcorr]=distance_wei(inver_parcorr);
    [ris(sogg).CPL_parz,ris(sogg).eff_parz,~]=charpath(D_parcorr,0,0);
    ris(sogg).GE_parz=efficiency_wei(mask_parcorr.*ris(sogg).FCThres_parz);
end
 
%% Task 10
 ascissa= 1:1:numSoggetti;
 f10=figure(10);
 set(f10,'NumberTitle','off','Name','G.E. e C.P.L','units',...
     'normalized','outerposition',[0 0 1 1])
 subplot(4,1,1)
 hold on
 plot(ascissa,[ris(1:numSoggetti).eff_corr],'-bo')
 plot(ascissa,[ris(1:numSoggetti).GE_corr],'-ro')
 title('Global efficiency per correlazione')
 legend('charpath','efficiencywei','Location','NorthEast')
 hold off
 xlim([0,23])
 subplot(4,1,2)
 plot(ascissa,[ris(1:numSoggetti).CPL_corr],'--om')
 title('Characteristic path length per correlazione')
 xlim([0,23])
 subplot(4,1,3)
 hold on
 plot(ascissa,[ris(1:numSoggetti).eff_parz],'ob-')
 plot(ascissa,[ris(1:numSoggetti).GE_parz],'or-')
 title('Global efficiency per correlazione parziale')
 legend('charpath','efficiencywei','Location','NorthEast')
 hold off
 xlim([0,23])
 subplot(4,1,4)
 plot(ascissa,[ris(1:numSoggetti).CPL_parz],'-om')
 title('Characteristic path length per correlazione parziale')


 %% Creazione matrici finali per la consegna
 
 
 %salvataggio
 nome_file =  'results_ANDREA_PITTARO.mat'
 save(nome_file,'')