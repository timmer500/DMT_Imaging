%% Lempel-Ziv analysis
% Chris Timmermann 2022

% Do this to unpack the code
mex LZ76.cpp
addpath  '/Users/christophertimmermann/Documents/MATLAB/fieldtrip-20200130'
% subject_list = {'S01TWD1', 'S02WTD2', 'S03CT', 'S06JBD1', 'S07MND2','S08RSD1','S09SRD1','S10RMD1','S11AED2','S12APD1','S13HKD2','S14LLD1','S14llD2','S15LPD1','S16JMD1','S17CSD2','S18SCD2','S19SGD1','S22EKD1','S23LPJD2','S25MMD1','S01TWD2', 'S02WTD1', 'S04CT', 'S06JBD2', 'S07MND1','S08RSD2','S09SRD2','S10RMD2','S11AED1','S12APD2','S13HKD1','S15LPD2','S16JMD2','S17CSD1','S18SCD1','S19SGD2','S22EKD2','S23LPJD1','S25MMD2'};
subject_list = {'S03CT', 'S23LPJD1','S18SCD1','S25MMD2'};
numsubjects = length(subject_list);
parentfolder = '/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed';

% Extract data and apply LZs to each trial

for s=1:numsubjects
    subject = subject_list{s};
    
    Filefold = [parentfolder '/' subject_list{s}];
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', s, subject);
    cd(Filefold);
    load('data_clean.mat');
    
    cfg = [];
    cfg.implicitref   = 'FCz'; %ERr I assume this is the ref for a brain prodcuts cap?
    cfg.reref = 'yes';
    cfg.dftfilter       = 'yes';
    cfg.channel = 'EEG';
    cfg.dftfreq       = [17.5 35];
    cfg.dftreplace    = 'zero';
    cfg.refchannel = 'EEG';
    [dataref] = ft_preprocessing(cfg, data_clean);
    
    % Select only EEG chans
                 
    datamat = fieldtrip2mat(dataref); %Turn into matrix format
    
                for jj = 1:size(datamat,2) %Loop through chans
                fprintf('LZs computation Channel %d .....', jj) ;
                sig = squeeze(datamat(:,jj,:));
                    for ii=1:size(sig,1) % loop through each trial
                        data = sig(ii,:);
                        LZs(jj,ii) = LZ76(data > mean(data)); %obtain the LZs value for that string
                        EntropyRt(jj,ii) = LZ76(data> mean(data))*(log2(length(data))/length(data)) ;% obtain entropy rate
                    end
                     fprintf('..completed \n');
                end

                save LZs  LZs
                save EntropyRt EntropyRt
                clear LZs EntropyRt
end
        
%% Save in averaged structures

data5m=0 ;% 1 if you want to do analysis over 5 mins, 0 if over 8 mins

subject_list = {'S01TWD1', 'S02WTD2', 'S03CT', 'S06JBD1', 'S07MND2','S08RSD1','S09SRD1','S10RMD1','S11AED2','S12APD1','S13HKD2','S14LLD1','S14LLD2','S15LPD1','S16JMD1','S17CSD2','S18SCD2','S19SGD1','S22EKD1','S23LPJD2','S25MMD1','S01TWD2', 'S02WTD1', 'S04CT', 'S06JBD2', 'S07MND1','S08RSD2','S09SRD2','S10RMD2','S11AED1','S12APD2','S13HKD1','S15LPD2','S16JMD2','S17CSD1','S18SCD1','S19SGD2','S22EKD2','S23LPJD1','S25MMD2'};
% subject_list = {'S03CT', 'S23LPJD1','S18SCD1','S25MMD2'};
numsubjects = length(subject_list);
parentfolder = '/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed';


for s=1:numsubjects
    subject = subject_list{s};
    
    Filefold = [parentfolder '/' subject_list{s}];
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', s, subject);
    cd(Filefold);
    load('LZs.mat'); % switch to specAll3b if you want the 3 second data to be done
    load('EntropyRt.mat'); % switch to specAll3b if you want the 3 second data to be done
    load('/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/Analysis/5minData/DMTfinal.mat')
    load('data_clean.mat')
    
    cfg = [];
    cfg.channel     = 'EEG';
    [datatrials] = ft_selectdata(cfg, data_clean);
    
    for ll=1:3 % generate 3 Fieldtrip structures: before, after and end based on an unrelated datafile
        LZ_ER{ll} = DMTfinal{1};
    end
    
     % this little stub below finds the nearest values lower than the timepoints
    % that we want to divide our data for analysis
    
    % first put all times of clean trials in an array
    for i=1:length(datatrials.time)    
        dattimes(i,:) = datatrials.time{i};
    end
    trials = dattimes(:,1);
    
    % define here times of interest
    m5 = (60*5);
    m8 = (60*8);
    m8s30 = (60*8)+30;% Add 30 secs corresponding to injection confound
    m13s30 = (13*60)+30; % add 30 secs to end of acute as in pilot
    m16 = (16*60);
    m23 = (23*60);
    m28 = (60*28);

    % extract the index of trials relevant for analysis

    times = [m5 m8 m8s30 m13s30 m16 m23 m28];
    for t=1:length(times)
        tt = times(t);

        k = bsxfun(@minus,trials,tt);
        k(k >= 0) = -inf;
        [~,ii] = max(k);
        ii(all(k == -inf)) = 0;
        indx(t) = ii;
    end
    
    if data5m == 1
    time(1,:) = [1 indx(1)];
    time(2,:) = [indx(3)+1 indx(4)];
    time(3,:) = [indx(6)+1 indx(7)];
    else
    time(1,:) = [1 indx(2)];
    time(2,:) = [indx(3)+1 indx(5)];
    time(3,:) = [indx(6)+1 indx(7)];
    end
    
    % assign relevant trials to times of interest
    
    
    for nn=1:length(time)   % loop through segments, before post and after    
            LZsav =  mean(LZs(:,[time(nn,1):time(nn,2)]),2); % take the mean out of each segment
            ERav =  mean(EntropyRt(:,[time(nn,1):time(nn,2)]),2); % take the mean out of each segment
        LZ_ER{nn}.powspctrm = [LZsav ERav];
        LZ_ER{nn}.freq = [1 2];
 
        clear dat newdat
    end
       
    save LZ_ER8min LZ_ER
    clear LZ_ER8min  freq* dat trials freqOsc freqall indx ii time dattimes datatrials

end

    %% Put it all on seprate DMT and plcebo conditions and grandaverage them

    data5m=0
    
PCBlist = {'S01TWD1', 'S02WTD2', 'S04CT', 'S06JBD1', 'S07MND2','S08RSD1','S09SRD1','S10RMD1','S11AED2','S12APD1','S13HKD2','S14LLD2','S15LPD1','S16JMD1','S17CSD2','S18SCD2','S19SGD1','S22EKD1','S23LPJD2','S25MMD1'};
DMTlist = {'S01TWD2', 'S02WTD1', 'S03CT', 'S06JBD2', 'S07MND1','S08RSD2','S09SRD2','S10RMD2','S11AED1','S12APD2','S13HKD1','S14LLD1','S15LPD2','S16JMD2','S17CSD1','S18SCD1','S19SGD2','S22EKD2','S23LPJD1','S25MMD2'};
% DMTlistalt = {'S01TWD2', 'S02WTD1', 'S04CT', 'S06JBD2', 'S07MND1','S08RSD2','S09SRD2','S10RMD2','S11AED1','S12APD2','S13HKD1','S14LLD1','S15LPD2','S16JMD2','S17CSD1','S18SCD1','S19SGD2','S21FHD1','S22EKD2','S23LPJD1','S25MMD2'};
% alt version has S21FHD1 which is a DMT sessions with no placebo control.
% which is ood for corrs
numsubjects = length(PCBlist);
parentfolder = '/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/';

for s=1:numsubjects
    subject = PCBlist{s};
    
    PCBfile = [parentfolder PCBlist{s}];
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', s, subject);
    cd(PCBfile);
    if data5m==1
    load('LZ_ER.mat');
    else
        load('LZ_ER8min.mat');
    end
    PCBpre{s} = LZ_ER{1};
    PCBpost{s} = LZ_ER{2};
    PCBfinal{s} = LZ_ER{3};
       
    DMTfile = [parentfolder DMTlist{s}];
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', s, subject);
    cd(DMTfile);
    if data5m==1
    load('LZ_ER.mat');
    else
        load('LZ_ER8min.mat');
    end
    load('LZ_ER.mat');
    DMTpre{s} = LZ_ER{1};
    DMTpost{s} = LZ_ER{2};
    DMTfinal{s} = LZ_ER{3};
    
end

if data5m==1
    cd '/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/Analysis/LZs/5mindata'
else
cd '/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/Analysis/LZs/8mindata'
end

save PCBpre PCBpre
save PCBpost PCBpost
save PCBfinal PCBfinal
save DMTpre DMTpre
save DMTpost DMTpost
save DMTfinal DMTfinal

%% GAs and substract Powspctrm

cfg = [];
cfg.keepindividual = 'yes';
GAPCBpre= ft_freqgrandaverage(cfg, PCBpre{:});
GAPCBpost= ft_freqgrandaverage(cfg, PCBpost{:});
GAPCBfinal= ft_freqgrandaverage(cfg, PCBfinal{:});
GADMTpre= ft_freqgrandaverage(cfg, DMTpre{:});
GADMTpost= ft_freqgrandaverage(cfg, DMTpost{:});
GADMTfinal= ft_freqgrandaverage(cfg, DMTfinal{:});

GAPCBdiff = GAPCBpre;
GAPCBdiff.powspctrm = GAPCBpost.powspctrm - GAPCBpre.powspctrm;
GADMTdiff = GADMTpre;
GADMTdiff.powspctrm = GADMTpost.powspctrm - GADMTpre.powspctrm;

GAPCBdiffFinal = GAPCBpre;
GAPCBdiffFinal.powspctrm = GAPCBfinal.powspctrm - GAPCBpre.powspctrm;
GADMTdiffFinal = GADMTpre;
GADMTdiffFinal.powspctrm = GADMTfinal.powspctrm - GADMTpre.powspctrm;

%% Do this to extract values


DMTdiff = nanmean(nanmean(GADMTdiff.powspctrm(:,:,1),3),2);
PCBdiff = nanmean(nanmean(GAPCBdiff.powspctrm(:,:,1),3),2);

DMTpost = nanmean(nanmean(GADMTpost.powspctrm(:,:,1),3),2);
PCBpost = nanmean(nanmean(GAPCBpost.powspctrm(:,:,1),3),2);


%% Now lets run the stats

dirtysubs = [9 14 17];
removedirty = 1; %Remove dirty subs? S13HK, S

LZ = 1 ;
ER = 2;
freqbands = {LZ, ER};

numfreq = length(freqbands);

nperms = 7500;
xdata = GADMTdiff;
ydata = GAPCBdiff;


if removedirty ==1
    xdata.powspctrm(dirtysubs,:,:) = [];
    ydata.powspctrm(dirtysubs,:,:) = [];
end

% mean(GADMTpost.powspctrm(:,:,[[25:31] [37:59]]),3);

for ff=1:numfreq
    fr = freqbands{ff};

% fr = 1;
    cfg = [];
    cfg.frequency        =fr;
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_depsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;
    cfg.channel = {'EEG'};
    cfg.avgoverfreq = 'yes';
    cfg.clusterstatistic = 'maxsum';
%     cfg.minnbchan        = 2;
    cfg.tail             = 0; % 2 sided
    cfg.clustertail      = 0;
    cfg.alpha            = 0.025; % this is a 2 tailed
    cfg.numrandomization = nperms;
    % specifies with which sensors other sensors can form clusters
    cfg_neighb.method    = 'template';
    cfg_neighb.layout    =   'easycapMR32final.mat';
    % cfg_neighb.layout    =   'EEG1020.lay';
    cfg_neighb.template = 'easycapMR32rad_neighb.mat';
    % cfg_neighb.template = 'elec1020_neighb.mat';
    cfg_neighb.feedback = 'no'; % Yes if you want to check the neighberhood structure
    cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, xdata);

    subj = size(xdata.powspctrm,1);
    design = zeros(2,2*subj);
    for i = 1:subj
      design(1,i) = i;
    end
    for i = 1:subj
      design(1,subj+i) = i;
    end
    design(2,1:subj)        = 1;
    design(2,subj+1:2*subj) = 2;

    cfg.design   = design;
    cfg.ivar = 2;
    cfg.uvar = 1;

    stat{ff} = ft_freqstatistics(cfg, xdata, ydata);
end



%% Figs for cluster stats results

rat = stat;
cfg = [];
cfg.alpha  = 0.025;
cfg.parameter = 'stat';
cfg.zlim   = [-8 8];
% cfg.contournum = 10;
% cfg.colormap = coolwarm;
cfg.layout = 'easycapMR32final.mat';
 cfg.highlightsymbolseries     = ['.','o'];
cfg.subplotsize    = [1 1];
cfg.marker = 'off'; 
cfg.comment = 'no';
cfg.highlightsizeseries  = [21 7];
cfg.highlightcolorpos  = [1 0 0];
cfg.highlightcolorneg  = [1 0 0];
cfg.gridscale = 1000;
if length(rat)==1
    A= isfieldRecursive(rat, 'negclusters');
    if A==0
        
       figure; ft_topoplotER(cfg, rat);
    else
        ft_clusterplot(cfg, rat);
    end  
title('', 'Position', [0 0.64 1], 'FontSize', 30,'FontName', 'Arial')
 dum = get(gcf);
  set(gcf,'Position',[dum.Position(1) dum.Position(2) dum.Position(3)/1.8 dum.Position(4)/1.8] , 'color', [1 1 1])  
  else
            for nn=1:length(rat)
                dat = rat{nn};
                
                if nn==3
                cfg.zlim   = [-4 4];
                else
                    cfg.zlim   = [-5 5];
                end

                
     
            B= isfieldRecursive(dat, 'negclusters');
                if B==0;
                    figure;ft_topoplotER(cfg, dat);
                else 
                     try
                   ft_clusterplot(cfg, dat);
                         catch ME
                         if (strcmp(ME.identifier,'no clusters present with a p-value lower than the specified alpha, nothing to plot'))
                             
                         end
                           figure;ft_topoplotER(cfg, dat);
                     end
                end
               
            title('', 'Position', [0 0.64 1], 'FontSize', 30,'FontName', 'Arial')
   dum = get(gcf);
  set(gcf,'Position',[dum.Position(1) dum.Position(2) dum.Position(3)/2.8 dum.Position(4)/2.8] , 'color', [1 1 1])  
  h = gcf
%   print(h,sprintf('DMTvsPLA%d',nn),'-dpng')
if nn==3
  export_fig(sprintf('DMTdiffvsPCBdiffT4%d.png',nn),'-m2')
else
      export_fig(sprintf('DMTdiffvsPCBdiffT4%d.png',nn),'-m2')
end
            end
     
end
         %% Run this to plot separate chans
        
        stat{2}.stat([1:14]) = 0;
        stat{2}.stat([17:end]) = 0;
        
        %% Get max absolute to report results
        
rat = stat;
 
if length(rat)==1
         vec = rat.stat;
       [Y,I] = max(abs(vec));
       maxtval=vec(I);
       
        A= isfieldRecursive(rat, 'negclusters','prob');
        B= isfieldRecursive(rat, 'posclusters','prob');
        if A==1 && B==0
            wa = rat.negclusters.prob;
            pr= min(wa);
        elseif A==0 && B==1
            pr= min(rat.posclusters.prob);   
             else
            pr = nan;
        end
else
    
   for d=1:length(rat)
       vec = rat{d}.stat;
       [Y,I(d)] = max(abs(vec));
       maxtval(d)=vec(I(d));
            dummy = find(rat{d}.prob<0.025);
       if isempty(dummy)
           mask{d}=nan;
       else
           mask{d}=dummy;
       end
       
        A= isfieldRecursive(rat{d}, 'negclusters','prob');
        B= isfieldRecursive(rat{d}, 'posclusters','prob');
        if A==1 && B==0
            pr(d)= min([rat{d}.negclusters.prob]);
        elseif A==0 && B==1
            pr(d)= min([rat{d}.posclusters.prob]);
        elseif A==0 && B==0
            pr(d) = nan;
        end
   end
end
   
   maxtval
   pr*2
%    pr*2   
save maskdiff mask

%% Plot the stats

removedirty=1;
import iosr.statistics.boxPlot

dirtysubs = [9 14 17];

PCBdiffmean = mean(GAPCBdiff.powspctrm(:,:,2),2);

PCBpostmean = mean(GAPCBpost.powspctrm(:,:,2),2);
PCBpremean = mean(GAPCBpre.powspctrm(:,:,2),2);

DMTdiffmean = mean(GADMTdiff.powspctrm(:,:,2),2);
DMTpostmean = mean(GADMTpost.powspctrm(:,:,2),2);
DMTpremean = mean(GADMTpre.powspctrm(:,:,2),2);
xdata = DMTdiffmean;
ydata = PCBdiffmean;
    
if removedirty ==1
    xdata(dirtysubs,:) = [];
    ydata(dirtysubs,:) = [];
end

    data = [ ydata xdata];

    figure;
    name = {'PCB','DMT'};
    nameelse = {'',''};
      h=   boxPlot([1 2],data,'boxColor','k', 'boxAlpha', 0.2,'lineWidth',2 ,'scatterMarker','.','meanSize',40);
      h.scatterSize = 800;
      h.scatterColor = [0.7 0 0];
      h.boxWidth = 0.4;
      h.scatterLayer = 'bottom';
        h.meanMarker = '.';
          h.meancolor = 'y';
      h.showScatter = true;
      h.showMean = true;
    %   h.showViolin = true;
    % h.notch = true; 
%     ylim([0 20])
       set(gcf, 'color', [1 1 1])
    set(gca,'fontsize', 25)
    
    set(gca,'XTickLabel' , name)
         ylabel('LZs')


    d = computeCohen_d(data(:,1),data(:,2), 'paired');
    d=abs(round(d,2));
    [H,P,CI,STATS] =  ttest(data(:,1),data(:,2));
    abs(round(STATS.tstat,2))
    H
    pr= round(P,3)
    n=0.1
    tstat = abs(round(STATS.tstat,2))
    df = STATS.df
%%
    
    y = [ydata;xdata ];

x = [ones(size(xdata,1),1) ; ones(size(ydata,1),1)+1];

figure
h = boxplot([ydata, xdata],'Labels',{'PCB','DMT'},'colors' ,'r')
% title('LZc')
set(h,'LineWidth',4)

hold on
scatter(x,y,40, 'k', 'filled','MarkerFaceAlpha',0.5)
line([x(1:size(ydata,1)) x(size(ydata,1)+1:end)]',[y(1:size(ydata,1)) y(size(ydata,1)+1:end)]','LineWidth',2,'Color', [0 0 0 0.3])

set(gca, 'linewidth', 2,'Fontsize',24,'Box', 'on')
set(gcf, 'color', [1 1 1],'Position',[744 784 338 266])
yl = ylabel('', 'FontWeight', 'bold','FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',20);

%   text(1.25 ,0.6,['p = ' num2str(round(pTotal,3))],'FontSize',20,'HorizontalAlignment','left','Color', 'k');
   
    
%     plot([1 2], [max([xdata; ydata])+n max([xdata; ydata])+n],'k','LineWidth',2)
    text(1.5,(max([xdata; ydata]+0.5)),['d = ' num2str(d)],'FontSize',25,'HorizontalAlignment','center','Color', 'r');
    if pr<0.001
        text(1.5,(max([xdata; ydata])),['p < 0.001'],'FontSize',25,'HorizontalAlignment','center','Color', 'k');
        else
        text(1.5,(max([xdata; ydata])),['p = ' num2str(pr)],'FontSize',25,'HorizontalAlignment','center','Color', 'k');
    end
    
     export_fig(sprintf('DMTvsPCBLZc.png'),'-m2')