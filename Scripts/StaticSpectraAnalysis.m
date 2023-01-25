%% EEG averaged analysis (Fig 5a)

% Put it all on a structure 

data5m = 0 ; % determine if 5 mins analysis oor 8
PCBlist = {'S01D1', 'S02D2', 'S04D1', 'S06D1', 'S07D2','S08D1','S09D1','S10D1','S11D2','S12D1','S13D2','S14D2','S15D1','S16D1','S17D2','S18D2','S19D1','S22D1','S23D2','S25D1'};
DMTlist = {'S01D2', 'S02D1', 'S04D2', 'S06D2', 'S07D1','S08D2','S09D2','S10D2','S11D1','S12D2','S13D1','S14D1','S15D2','S16D2','S17D1','S18D1','S19D2','S22D2','S23D1','S25D2'};
%DMT list alt has one more subject
numsubjects = length(DMTlist);
parentfolder = '/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/';

for s=1:numsubjects
    subject = PCBlist{s};
    
    PCBfile = [parentfolder PCBlist{s}];
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', s, subject);
    cd(PCBfile);
    load('data_clean.mat');
    
    cfg = [];
    cfg.implicitref   = 'FCz'; %ERr I assume this is the ref for a brain prodcuts cap?
    cfg.reref = 'yes';
    cfg.dftfilter       = 'yes';
    cfg.dftfreq       = [17.5 35];
    cfg.dftreplace    = 'zero';
    cfg.refchannel = 'EEG';
    [dataref] = ft_preprocessing(cfg, data_clean);
%     
    st = 0;
    m5 = (60*5);
    m8 = (60*8);
    m8s30 = (60*8)+30;% Add 30 secs corresponding to injection confound
    m13s30 = (13*60)+30; % add 30 secs to end of acute as in pilot
    m16s30 = (16*60)+30;
    m20 = (60*20);
    m23 = (23*60);
    m28 = (60*28);

    times = [m5 m8 m8s30 m13s30 m16s30 m20 m23 m28];
    % for nn=1:length(times)
    %     
    % [val,idx(nn)]=min(abs(trs(:,1)-times(nn)));
    % end
    if data5m == 1
    time(1,:) = [st times(1)-0.004];
    time(2,:) = [times(3) times(4)-0.004];
    time(3,:) = [times(7) times(8)-0.004];
    else
    time(1,:) = [st times(2)-0.004];
    time(2,:) = [times(3) times(5)-0.004];
    time(3,:) = [times(6) times(8)-0.004];
    end

    for nn=1:length(time)

        cfg             = [];           % segment the continuous data into two-second pieces
        cfg.toilim    = time(nn,:);
        data_trialscut     = ft_redefinetrial(cfg, dataref);

        % FFT for low frequencies:
        cfg             = [];
        cfg.method      = 'mtmfft';     % analyses entire spectrum for the entire data length
        cfg.output      = 'pow';
        cfg.channel     = 'EEG';
        cfg.taper       = 'hanning';       % default = dpss
        cfg.foi         = 1: 0.5: 30;    % frequency band of interest (min:stepsize:max)
        cfg.keeptrials = 'no';
        flowr = ft_freqanalysis(cfg, data_trialscut);

        cfg             = [];
        cfg.method      = 'mtmfft';     % analyses entire spectrum for the entire data length
        cfg.output      = 'pow';
        cfg.channel     = 'EEG';
        cfg.taper       = 'dpss';       % default = dpss
        cfg.foilim      = [30 45];    % frequency band of interest
        cfg.tapsmofrq   = 3;          % frequency smoothing
        cfg.keeptrials = 'no';
        fhighr = ft_freqanalysis(cfg, data_trialscut);

        % concatenate high and low frequencies:
        fallr{nn} = flowr;
        fallr{nn}.freq = cat(2, flowr.freq, fhighr.freq);
        fallr{nn}.powspctrm = cat(2, flowr.powspctrm, fhighr.powspctrm);


    end
    
    PCBpre{s} = fallr{1};
    PCBpost{s} = fallr{2};
    PCBfinal{s} = fallr{3};
    
    DMTfile = [parentfolder DMTlist{s}];
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', s, subject);
    cd(DMTfile);
    load('data_clean.mat');
    
    cfg = [];
    cfg.implicitref   = 'FCz'; %ERr I assume this is the ref for a brain prodcuts cap?
    cfg.reref = 'yes';
    cfg.dftfilter       = 'yes';
    cfg.dftfreq       = [17.5 35];
    cfg.dftreplace    = 'zero';
    cfg.refchannel = 'EEG';
    [dataref] = ft_preprocessing(cfg, data_clean);
    
    times = [m5 m8 m8s30 m13s30 m16s30 m20 m23 m28];
    % for nn=1:length(times)
    %     
    % [val,idx(nn)]=min(abs(trs(:,1)-times(nn)));
    % end
    if data5m == 1
    time(1,:) = [st times(1)-0.004];
    time(2,:) = [times(3) times(4)-0.004];
    time(3,:) = [times(7) times(8)-0.004];
    else
    time(1,:) = [st times(2)-0.004];
    time(2,:) = [times(3) times(5)-0.004];
    time(3,:) = [times(6) times(8)-0.004];
    end

    for nn=1:length(time)

        cfg             = [];           % segment the continuous data into two-second pieces
        cfg.toilim    = time(nn,:);
        data_trialscut     = ft_redefinetrial(cfg, dataref);

        % FFT for low frequencies:
        cfg             = [];
        cfg.method      = 'mtmfft';     % analyses entire spectrum for the entire data length
        cfg.output      = 'pow';
        cfg.channel     = 'EEG';
        cfg.taper       = 'hanning';       % default = dpss
        cfg.foi         = 1: 0.5: 30;    % frequency band of interest (min:stepsize:max)
        cfg.keeptrials = 'no';
        flowr = ft_freqanalysis(cfg, data_trialscut);

        cfg             = [];
        cfg.method      = 'mtmfft';     % analyses entire spectrum for the entire data length
        cfg.output      = 'pow';
        cfg.channel     = 'EEG';
        cfg.taper       = 'dpss';       % default = dpss
        cfg.foilim      = [30 45];    % frequency band of interest
        cfg.tapsmofrq   = 3;          % frequency smoothing
        cfg.keeptrials = 'no';
        fhighr = ft_freqanalysis(cfg, data_trialscut);

        % concatenate high and low frequencies:
        fallr{nn} = flowr;
        fallr{nn}.freq = cat(2, flowr.freq, fhighr.freq);
        fallr{nn}.powspctrm = cat(2, flowr.powspctrm, fhighr.powspctrm);


    end   
    DMTpre{s} = fallr{1};
    DMTpost{s} = fallr{2};
    DMTfinal{s} = fallr{3};  
    
end

%% Save this shiat

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

%% Extract values for correlations

x =1;[d, dbix] = min(abs(x-GADMTfinal.freq));
x=4;[d, deix] = min(abs(x-GADMTfinal.freq)) ;
x=8;[d, teix] = min(abs(x-GADMTfinal.freq)) ;
% x=10;[d, aleix] = min(abs(x-DMTsub{1}.freq)) ;
x =13;[d, aheix] = min(abs(x-GADMTfinal.freq)) ;
x=30;[d, beix] = min(abs(x-GADMTfinal.freq)) ;
x=45;[d, geix] = min(abs(x-GADMTfinal.freq));


DeltaDMTdiff = nanmean(nanmean(GADMTdiff.powspctrm(:,:,[dbix deix]),3),2);
ThetaDMTdiff = nanmean(nanmean(GADMTdiff.powspctrm(:,mask{2},[deix teix]),3),2);
AlphaDMTdiff = nanmean(nanmean(GADMTdiff.powspctrm(:,mask{3},[teix aheix]),3),2);
BetaDMTdiff = nanmean(nanmean(GADMTdiff.powspctrm(:,mask{4},[aheix beix]),3),2);
GammaDMTdiff = nanmean(nanmean(GADMTdiff.powspctrm(:,mask{5},[beix geix]),3),2);

DeltaPCBdiff = nanmean(nanmean(GAPCBdiff.powspctrm(:,:,[dbix deix]),3),2);
ThetaPCBdiff = nanmean(nanmean(GAPCBdiff.powspctrm(:,mask{2},[deix teix]),3),2);
AlphaPCBdiff = nanmean(nanmean(GAPCBdiff.powspctrm(:,mask{3},[teix aheix]),3),2);
BetaPCBdiff = nanmean(nanmean(GAPCBdiff.powspctrm(:,mask{4},[aheix beix]),3),2);
GammaPCBdiff = nanmean(nanmean(GAPCBdiff.powspctrm(:,mask{5},[beix geix]),3),2);

DeltaDMT = nanmean(nanmean(GADMTpost.powspctrm(:,:,[dbix deix]),3),2);
ThetaDMT = nanmean(nanmean(GADMTpost.powspctrm(:,mask{2},[deix teix]),3),2);
AlphaDMT = nanmean(nanmean(GADMTpost.powspctrm(:,mask{3},[teix aheix]),3),2);
BetaDMT = nanmean(nanmean(GADMTpost.powspctrm(:,mask{4},[aheix beix]),3),2);
GammaDMT = nanmean(nanmean(GADMTpost.powspctrm(:,mask{5},[beix geix]),3),2);

DeltaPCB = nanmean(nanmean(GAPCBpost.powspctrm(:,:,[dbix deix]),3),2);
ThetaPCB = nanmean(nanmean(GAPCBpost.powspctrm(:,mask{2},[deix teix]),3),2);
AlphaPCB = nanmean(nanmean(GAPCBpost.powspctrm(:,mask{3},[teix aheix]),3),2);
BetaPCB = nanmean(nanmean(GAPCBpost.powspctrm(:,mask{4},[aheix beix]),3),2);
GammaPCB = nanmean(nanmean(GAPCBpost.powspctrm(:,mask{5},[beix geix]),3),2);

DeltaDMTmean = nanmean(nanmean(GADMTpost.powspctrm(:,:,[dbix deix]),3),2);
ThetaDMTmean = nanmean(nanmean(GADMTpost.powspctrm(:,:,[deix teix]),3),2);
AlphaDMTmean = nanmean(nanmean(GADMTpost.powspctrm(:,:,[teix aheix]),3),2);
BetaDMTmean = nanmean(nanmean(GADMTpost.powspctrm(:,:,[aheix beix]),3),2);
GammaDMTmean = nanmean(nanmean(GADMTpost.powspctrm(:,:,[beix geix]),3),2);

DeltaPCBmean = nanmean(nanmean(GAPCBpost.powspctrm(:,:,[dbix deix]),3),2);
ThetaPCBmean = nanmean(nanmean(GAPCBpost.powspctrm(:,:,[deix teix]),3),2);
AlphaPCBmean = nanmean(nanmean(GAPCBpost.powspctrm(:,:,[teix aheix]),3),2);
BetaPCBmean = nanmean(nanmean(GAPCBpost.powspctrm(:,:,[aheix beix]),3),2);
GammaPCBmean = nanmean(nanmean(GAPCBpost.powspctrm(:,:,[beix geix]),3),2);


DeltaDMTdiffmean = nanmean(nanmean(GADMTdiff.powspctrm(:,:,[dbix deix]),3),2);
ThetaDMTdiffmean = nanmean(nanmean(GADMTdiff.powspctrm(:,:,[deix teix]),3),2);
AlphaDMTdiffmean = nanmean(nanmean(GADMTdiff.powspctrm(:,:,[teix aheix]),3),2);
BetaDMTdiffmean = nanmean(nanmean(GADMTdiff.powspctrm(:,:,[aheix beix]),3),2);
GammaDMTdiffmean = nanmean(nanmean(GADMTdiff.powspctrm(:,:,[beix geix]),3),2);

DeltaPCBdiffmean = nanmean(nanmean(GAPCBdiff.powspctrm(:,:,[dbix deix]),3),2);
ThetaPCBdiffmean = nanmean(nanmean(GAPCBdiff.powspctrm(:,:,[deix teix]),3),2);
AlphaPCBdiffmean = nanmean(nanmean(GAPCBdiff.powspctrm(:,:,[teix aheix]),3),2);
BetaPCBdiffmean = nanmean(nanmean(GAPCBdiff.powspctrm(:,:,[aheix beix]),3),2);
GammaPCBdiffmean = nanmean(nanmean(GAPCBdiff.powspctrm(:,:,[beix geix]),3),2);


%% Run stats

removedirty = 1; %Remove dirty subs?  S11AE, S7CS,  S16JM (due to having alpha peak at 7 Hz) S25MM, 
noslice = 0; % average beta and gamma bands avoiding the slice artifact?
delta = [1 4];
theta = [4 8];
alpha = [8 13];


if noslice==1
    beta = 46;
    gamma=47;
else
    beta = [13 30];
    gamma = [30 45];
end
betalow = [13 16.5];
betahigh = [18.5 30];
gammalow = [30 34];
gammahigh = [36 45];
freqbands = {delta, theta, alpha, beta, gamma};

numfreq = length(freqbands);

nperms = 7500;
xdata = GADMTdiff;
ydata = GAPCBdiff;

if noslice==1
xdata.powspctrm(:,:,91) =  mean(xdata.powspctrm(:,:,[[25:31] [37:59]]),3);
xdata.powspctrm(:,:,92) =  mean(xdata.powspctrm(:,:,[[59:67] [73:90]]),3);
ydata.powspctrm(:,:,91) =  mean(ydata.powspctrm(:,:,[[25:31] [37:59]]),3);
ydata.powspctrm(:,:,92) =  mean(ydata.powspctrm(:,:,[[59:67] [73:90]]),3);
xdata.freq(end+1:end+2) = [46 47];
ydata.freq(end+1:end+2) = [46 47];
end

if removedirty ==1
    
    xdata.powspctrm(dirtysubs,:,:) = [];
    ydata.powspctrm(dirtysubs,:,:) = [];
end



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
    cfg.minnbchan        = 1;
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
cfg.zlim   = [-4 4];
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
                    cfg.zlim   = [-4 4];
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
        
        stat{2}.stat([1:6]) = 0;
        stat{2}.stat([9:end]) = 0;
        
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
        end
   end
end
   
   maxtval
   pr*2
   
   save maskdiffvsdiff mask
%    pr*2   

%% Plot Spectra
dirtysubs = [9 14 17]; % just remove S12 LL and S16 JM as they are excluded from fMRI
logspectra=0;
removedirty=1;
long=0; % plot all conditions or just 2?

if logspectra==1
  GADMTpost.powspctrm =   log(GADMTpost.powspctrm);
    GAPCBpost.powspctrm =   log(GAPCBpost.powspctrm);
  GADMTpre.powspctrm =   log(GADMTpre.powspctrm);
  GAPCBpre.powspctrm =   log(GAPCBpre.powspctrm);
end

xdata = GADMTpost;
ydata = GAPCBpost;
zdata = GADMTpre;
idata = GAPCBpre;

if removedirty ==1
    
    xdata.powspctrm(dirtysubs,:,:) = [];
    ydata.powspctrm(dirtysubs,:,:) = [];
    zdata.powspctrm(dirtysubs,:,:) = [];
    idata.powspctrm(dirtysubs,:,:) = [];
end
numsubs = size(xdata.powspctrm,1);

lim = [1 45]; % select freqs to plot
cols = linspecer(3);
%     m = log(squeeze(nanmean(DMT.powspctrm,1))); % mean across trials
    m = squeeze(nanmean(xdata.powspctrm,1)); % mean across trials
    x = mean(m);
    y = xdata.freq;
    z = std(m)/sqrt(numsubs);
    
    m2 = squeeze(nanmean(ydata.powspctrm,1)); % mean across trials
    x2 = mean(m2);
    z2 = std(m2)/sqrt(numsubs);   
    
    m3 = squeeze(nanmean(zdata.powspctrm,1)); % mean across trials
    x3 = mean(m3);
    z3 = std(m3)/sqrt(numsubs);   
    
    m4 = squeeze(nanmean(idata.powspctrm,1)); % mean across trials
    x4 = mean(m4);
    z4 = std(m4)/sqrt(numsubs);   
    
    figure; 
    if long==1
    [l,p] = boundedline(y, x, z, y, x2, z2, y, x3, z3, y, x4, z4, 'alpha', 'transparency', 0.28);
    else
    [l,p] = boundedline(y, x, z, y, x2, z2, 'alpha', 'transparency', 0.28);
    end
    
    set(l(1), 'linewidth', 2, 'color', cols(2,:));
    set(p(1), 'facecolor', cols(2,:));
    set(l(2), 'linewidth', 2, 'color', cols(1,:));
    set(p(2), 'facecolor', cols(1,:));
    if long==1
    set(l(3), 'linewidth', 1.5, 'color', [0 0 0]);
    set(p(3), 'facecolor', [0 0 0]);
    set(l(4), 'linewidth', 1.5, 'color', [0.3 0.3 0.3]);
    set(p(4), 'facecolor', [0.3 0.3 0.3]);
    end

    allfigs = allchild(gcf);
    set(allfigs(10), 'linewidth', 1,'Fontsize',25, 'Box', 'on')
    ylim = [-0.2 3.5];

    set(gcf, 'color', [1 1 1],'Position',[539 555 373 307]);
    set(gca, 'xlim', lim,'ylim', ylim,'XTick',[4 8 13 30]); % here change the frequencies to explore [10:10:40]
 % But put labales on the first one
 if long==1
        legend('DMT post', 'DMT pre', 'PCB post', 'PCB pre');
 else
     legend('DMT','PCB');
 end
        allfigs = allchild(gcf);
        set(allfigs(10), 'fontsize',24, 'FontWeight', 'bold', 'Box', 'off', 'FontName', 'Arial', 'location', 'northeast', 'LineWidth', 10);
        xl = xlabel('Frequency (Hz)', 'FontWeight', 'bold','FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',24);
        yl = ylabel('Power (\muv^2)', 'FontWeight', 'bold', 'FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',24);
%          yl = ylabel('log Power', 'FontWeight', 'bold', 'FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',24);


export_fig(sprintf('SpectraDMTvsPCB.png'),'-m2')