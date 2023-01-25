%% Dynamic Spectra analysis (Fig 5c)
% Chris Timmermann 2022

% Determine a regressor to correlate EEG data
truerat=0;
ranking=0;
interp=1; % intepolate across 840 points?

numDMT = xlsread('/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/CorrMatrix/IntensityRatings.xlsx','DMT');

if truerat==1
numDMT(numDMT==11)=10;
end

numPCB = xlsread('/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/CorrMatrix/IntensityRatings.xlsx','PCB');

S = numDMT - numPCB;

S = [zeros(20,1) S];

if ranking==1
    for ss=1:size(S,1)
         data = S(ss,:);
        X = data;
        [temp,S(ss,:)]  = ismember(X,unique(X));
    end
end

S = S';    

if interp==1
    % define query points for interpolation
    qv = 0:28/840:(28-0/840);
    qv(1) = [];

    for i=1:size(S,2)
        int = S(:,i);
        Sinterp(i,:) = interp1(0:28, int, qv);
    end
    clear S
    S = Sinterp;
end

badsubs = [1 6 7 12 14 18];
S(badsubs,:) = [];

%% Cluster corrected temporal correlations
% Replacing LME as a more specific measure that takes into account cluster
% activity. This makes more sense to use.
%     DMTlist = {'S02WTD1', 'S03CT', 'S06JBD2', 'S07MND1','S10RMD2','S11AED1','S12APD2','S13HKD1','S15LPD2','S17CSD1','S18SCD1','S19SGD2','S23LPJD1','S25MMD2', ...  
% Set parameters
IRASA = 1;

DMTlist = {'S02WTD1', 'S03CT', 'S06JBD2', 'S07MND1','S10RMD2','S11AED1','S12APD2','S13HKD1','S15LPD2','S17CSD1','S18SCD1','S19SGD2','S23LPJD1','S25MMD2'};  
PCBlist = {'S02WTD2', 'S04CT', 'S06JBD1', 'S07MND2','S10RMD1','S11AED2','S12APD1','S13HKD2','S15LPD1','S17CSD2','S18SCD2','S19SGD1','S23LPJD2','S25MMD1'}; % use S21FHD1 if wnat to extend DMT listparentfolder = '/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/';
parentfolder = '/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/';

if interp==1
    time4corr = [1:840];
else
    time4corr = [1:28];
end
%load some dummy variables
load('/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/Analysis/IRASA/Corrs/dummfreq')

if IRASA==1
    
load('/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/S01TWD1/spec2s.mat')
    freqvar = spec2s{1}.freq;
else
    load('/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/Analysis/8minData/DMTpost.mat')
    freqvar = DMTpost{1}.freq;
end
% Int = Int(2:17,:);
dumm.freq = (1:6);
dumm.powspctrm = [];
% We get indeces from frequency bands of interest
x =1;[d, dbix] = min(abs(x-freqvar));
x=4;[d, deix] = min(abs(x-freqvar)) ;
x=8;[d, teix] = min(abs(x-freqvar)) ;
x =13;[d, aheix] = min(abs(x-freqvar)) ;
x=30;[d, beix] = min(abs(x-freqvar)) ;
x=45;[d, geix] = min(abs(x-freqvar));


for s=1:length(DMTlist)
    subject = DMTlist{s};

    DMTfile = [parentfolder DMTlist{s}];
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', s, subject);
    cd(DMTfile);
    if IRASA==1
    load('IRall.mat');
    var = IRall;
    else
        load('allSpectra.mat');
        var = allSpectra;
    end
    load('cleantrs.mat')
    load('LZs.mat')
    delta = squeeze(mean(var(dbix:deix,:,:)));
    theta  = squeeze(mean(var(deix:teix,:,:)));
    alpha  = squeeze(mean(var(teix:aheix,:,:)));
    beta  = squeeze(mean(var(aheix:beix,:,:)));
    gamma  = squeeze(mean(var(beix:geix,:,:)));
    lz = LZs';

    int = S(s,:);

    deltaALL = nan(840,32);
    thetaALL = nan(840,32);
    alphaALL = nan(840,32);
    betaALL = nan(840,32);
    gammaALL = nan(840,32);
    lzALL = nan(840,32);
    deltaALL(cleantrs,:) = delta;
    thetaALL(cleantrs,:) = theta;
    alphaALL(cleantrs,:) = alpha;
    betaALL(cleantrs,:) = beta;
    gammaALL(cleantrs,:) = gamma;
    lzALL(cleantrs,:) = lz;
    
       PCBfile = [parentfolder PCBlist{s}];
    cd(PCBfile);
    
    if IRASA==1
    load('IRall.mat');
    var = IRall;
    else
        load('allSpectra.mat');
        var = allSpectra;
    end
    load('cleantrs.mat')
    load('LZs.mat')
    
    delta = squeeze(mean(var(dbix:deix,:,:)));
    theta  = squeeze(mean(var(deix:teix,:,:)));
    alpha  = squeeze(mean(var(teix:aheix,:,:)));
    beta  = squeeze(mean(var(aheix:beix,:,:)));
    gamma  = squeeze(mean(var(beix:geix,:,:)));
    lz = LZs';

    deltaALLp = nan(840,32);
    thetaALLp = nan(840,32);
    alphaALLp = nan(840,32);
    betaALLp = nan(840,32);
    gammaALLp = nan(840,32);
    lzALLp = nan(840,32);
    deltaALLp(cleantrs,:) = delta;
    thetaALLp(cleantrs,:) = theta;
    alphaALLp(cleantrs,:) = alpha;
    betaALLp(cleantrs,:) = beta;
    gammaALLp(cleantrs,:) = gamma;
    lzALLp(cleantrs,:) = lz;
    

    tps = [1:30:840]; % determine timepoints, average minutes
    for e=1:32 % electrodes
           if interp==1
                deltam(1,:) = deltaALL(:,e);
                thetam(1,:) = thetaALL(:,e);
                alpham(1,:) = alphaALL(:,e);        
                betam(1,:) = betaALL(:,e);
                gammam(1,:) = gammaALL(:,e);
                lzm(1,:) = lzALL(:,e);
                
                
                deltamp(1,:) = deltaALLp(:,e);
                thetamp(1,:) = thetaALLp(:,e);
                alphamp(1,:) = alphaALLp(:,e);        
                betamp(1,:) = betaALLp(:,e);
                gammamp(1,:) = gammaALLp(:,e);
                lzmp(1,:) = lzALLp(:,e);
                
           else
                 for ii=1:length(tps)
                     tim = tps(ii);
                    %do correlations per channel. First build a matrix between channel averaged
                    %freq band and intensity for each subject         
                        % Build matrix per subject and per electrode
                        deltam(ii) = nanmean(deltaALL(tim:(tim+29),e),1);
                        thetam(ii) = nanmean(thetaALL(tim:(tim+29),e),1);
                        alpham(ii) = nanmean(alphaALL(tim:(tim+29),e),1);        
                        betam(ii) = nanmean(betaALL(tim:(tim+29),e),1);
                        gammam(ii) = nanmean(gammaALL(tim:(tim+29),e),1);
                        lzm(ii) = nanmean(lzALL(tim:(tim+29),e),1);


                        deltamp(ii) = nanmean(deltaALLp(tim:(tim+29),e),1);
                        thetamp(ii) = nanmean(thetaALLp(tim:(tim+29),e),1);
                        alphamp(ii) = nanmean(alphaALLp(tim:(tim+29),e),1);        
                        betamp(ii) = nanmean(betaALLp(tim:(tim+29),e),1);
                        gammamp(ii) = nanmean(gammaALLp(tim:(tim+29),e),1);
                        lzmp(ii) = nanmean(lzALLp(tim:(tim+29),e),1);
                 end
           end
     
         % now subtract DMT form placebo before doing correlations
         deltaf = deltam - deltamp;
         thetaf = thetam - thetamp;
         alphaf = alpham - alphamp;
         betaf = betam - betamp;
         gammaf = gammam - gammamp;
         lzf = lzm - lzmp;

    % do our correlation and normalize using atanh
    deltacorr(e) =atanh(corr(deltaf(time4corr)',int(time4corr)','type','Pearson','rows','complete'));
    thetacorr(e) =atanh(corr(thetaf(time4corr)',int(time4corr)','type','Pearson','rows','complete'));
%         alphalcorr(e,sub) = corr(alphal',Int(:,sub),'type','Pearson');
    alphacorr(e) =atanh(corr(alphaf(time4corr)',int(time4corr)','type','Pearson','rows','complete'));
    betacorr(e) =atanh(corr(betaf(time4corr)',int(time4corr)','type','Pearson','rows','complete'));
    gammacorr(e) =atanh(corr(gammaf(time4corr)',int(time4corr)','type','Pearson','rows','complete'));
    lzcorr(e) =atanh(corr(lzf(time4corr)',int(time4corr)','type','Pearson','rows','complete'));

    end
    % put correlation info in a FT structur for cluster stats
    DMTcorrs{s} = dummfreq;
    DMTcorrs{s}.powspctrm(:,1) = deltacorr';
    DMTcorrs{s}.powspctrm(:,2) = thetacorr';
    %     DMT1stlevel{sub}.powspctrm(:,3) = alphalcorr(:,sub);
    DMTcorrs{s}.powspctrm(:,3) = alphacorr';
    DMTcorrs{s}.powspctrm(:,4) = betacorr';
    DMTcorrs{s}.powspctrm(:,5) = gammacorr';
    DMTcorrs{s}.powspctrm(:,6) = lzcorr';
    DMTcorrs{s}.freq = 1:6;

    % we generate our control variable putting zeros
    DMTcorrsControl{s} = DMTcorrs{1};
    DMTcorrsControl{s}.powspctrm = zeros(32,6);

    deltas(:,:,s) = deltaALL -deltaALLp ;
    thetas(:,:,s) =thetaALL -thetaALLp;
    alphas(:,:,s) =alphaALL -alphaALLp;
    betas(:,:,s) = betaALL -betaALLp;
    gammas(:,:,s) = gammaALL -gammaALLp;
    lzs(:,:,s) = lzALL -lzALLp;
end


%% now run stats
% GA for stats
addpath  '/Users/christophertimmermann/Documents/MATLAB/fieldtrip-20200130'
cfg.keepindividual = 'yes';
GADMTcorrs= ft_freqgrandaverage(cfg, DMTcorrs{:});
GADMTcorrsControl= ft_freqgrandaverage(cfg, DMTcorrsControl{:});
%% Do stats

dirtysubs = [6 12];
removedirty = 1; %Remove dirty subs?  S11AE, S7CS, S25MM,  S16JM (due to having alpha peak at 7 Hz)

nperms = 7500;
xdata = GADMTcorrs;
ydata = GADMTcorrsControl;

if removedirty ==1
    
    xdata.powspctrm(dirtysubs,:,:) = [];
    ydata.powspctrm(dirtysubs,:,:) = [];
end
numfreq = length(xdata.freq);

for ff=1:numfreq
    fr = ff;
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
    masksigNormal{ff} = find(stat{ff}.mask);
end

%% Figs for cluster stats results
%    colors = cbrewer('div', 'RdBu', 64);
% colors = flipud(colors); % puts red on top, blue at the bottom
% colormap(colors);


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
                cfg.zlim   = [-8 8];
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
  export_fig(sprintf('DMTdiffvsPCBdiffT8%d.png',nn),'-m2')
else
      export_fig(sprintf('DMTdiffvsPCBdiffT4%d.png',nn),'-m2')
end
            end
     
end

 
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
  
  
 %% Professional
 mkdir FigsDMTvsDumm
cd FigsDMTvsDumm

rat = stat;
cfg = [];
cfg.alpha  = 0.025;
cfg.parameter = 'stat';
cfg.zlim   = [-6 6];
% cfg.contournum = 10;
% cfg.colormap = colors;
cfg.layout = 'easycapMR32final.mat';
cfg.subplotsize    = [1 1];
cfg.marker = 'off'; 
cfg.comment = 'no';
cfg.highlightsizeseries  = [15 15 15 15];
cfg.highlightcolorpos  = [1 1 1];
cfg.highlightcolorneg  = [0 0 0];
% cfg.gridscale = 2000;

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
            for nn=1:length(stat)
                dat = rat{nn};
                
     
            dat = rat{nn};
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
  set(gcf,'Position',[dum.Position(1) dum.Position(2) dum.Position(3)/1.8 dum.Position(4)/1.8] , 'color', [1 1 1])  
  h = gcf
%   print(h,sprintf('DMTvsPLA%d',nn),'-dpng')
%   export_fig(sprintf('DMTvsPLAT5%d.png',nn))
            end
     
end
         export_fig('AftervsAfterT63.png') 

         
 %% Plot the timecourses of stuff
 
%  save('IRASAplot.mat', 'betas', '-append')
truerat=0;
ranking=0;
interp=1; % intepolate across 840 points?

numDMT = xlsread('/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/CorrMatrix/IntensityRatings.xlsx','DMT');

if truerat==1
numDMT(numDMT==11)=10;
end

numPCB = xlsread('/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/CorrMatrix/IntensityRatings.xlsx','PCB');

S = numDMT - numPCB;

S = [zeros(20,1) S];

if ranking==1
    for ss=1:size(S,1)
         data = S(ss,:);
        X = data;
        [temp,S(ss,:)]  = ismember(X,unique(X));
    end
end

S = S';    

if interp==1
    % define query points for interpolation
    qv = 0:28/840:(28-0/840);
    qv(1) = [];

    for i=1:size(S,2)
        int = S(:,i);
        Sinterp(i,:) = interp1(0:28, int, qv);
    end
    clear S
    S = Sinterp;
end

badsubs = [1 6 7 12 14 18];
S(badsubs,:) = [];

load('/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/Analysis/EEGvsIntensity/IRASAplot.mat','deltas','masksig')
deltaIrasa = deltas;
% maskdeltaIrasa = masksig{1};

% clear deltas masksig;

% load('/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/Analysis/EEGvsIntensity/IRASAplotNormal.mat')


 deltasig = squeeze(mean(deltaIrasa(:,masksig{1},:),2));
%   alphasig = squeeze(mean(alphas(:,masksigNormal{3},:),2));
%   betasig = squeeze(mean(betas(:,masksigNormal{4},:),2));
%     lzsig = squeeze(mean(lzs(:,masksigNormal{6},:),2));
  alphasig = squeeze(mean(alphas(:,masksig{3},:),2));
  betasig = squeeze(mean(betas(:,masksig{4},:),2));
    lzsig = squeeze(mean(lzs(:,masksig{6},:),2));

deltainterp = fillmissing(deltasig,'linear',1);
alphasinterp = fillmissing(alphasig,'linear',1);
betasinterp = fillmissing(betasig,'linear',1);
lzsinterp = fillmissing(lzsig,'linear',1);

% convolve and obtain mean with sliding window as in MRI for display
% purposes only
alldata = cat(3,deltainterp,alphasinterp,betasinterp,lzsinterp);

dirtysubs = [6 12];
alldata(:,dirtysubs,:) = [];

% slidwin=22;
% sigmas=3;
% padded = slidwin/2;
% 
% for ff=1:size(alldata,3) % go through each measure 
%     for ss=1:size(deltainterp,2)
% %         alpharegpad = zeros(size(alldata,2),size(alldata,1)+ slidwin) ;% initialize padded variable
%         forconv = alldata(:,:,ff);
%         final = [zeros(size(alldata,2), padded)  forconv' zeros(size(alldata,2), padded)] ;
%         [c , alldataconv(:,:,ff)] = tapered_sliding_window(final',slidwin, sigmas);
%     end
%     
% end

for ff=1:size(alldata,3)
    alldataff = alldata(:,:,ff);
        alldataff = [zeros(11,size(alldataff,2)) ; alldataff ;zeros(11,size(alldataff,2)) ];
    alldatalowp= ft_preproc_lowpassfilter(alldataff',0.5,0.01); % take mean of chans of interest, detrend and lowpass as in MRI
                alldatalowp(:,[1:11 end-10:end]) = [];
                    alldatalowpass(:,:,ff) = alldatalowp;

end




for ff=1:size(alldatalowpass,3)
        alldatalowpassnorm(:,:,ff) =  rescale(alldatalowpass(:,:,ff)');
end
% alldataconvnorm =  normalize(alldataconv,2,'range');

alldatalowpassnormbased = alldatalowpassnorm - mean(alldatalowpassnorm(1:240,:,:),1);

% for ff=1:size(alldatalowpassnormbased,3)
%     plot(mean(alldatalowpassnormbased(:,:,ff),2))
%     hold on
% end
% 
Slowpass = ft_preproc_lowpassfilter(S,0.5,0.01); % take mean of chans of interest, detrend and lowpass as in MRI
% yyaxis right
% 
% plot(mean(Slowpass,1))

% Do ProPlot

Slowpass([6 12],:) = [];
ddata = alldatalowpassnormbased(:,:,1);
adata = alldatalowpassnormbased(:,:,2);
bdata = alldatalowpassnormbased(:,:,3);
ldata = alldatalowpassnormbased(:,:,4);

zdata = Slowpass';

numsubs = size(ddata,2);

lim = [-29 870]; % select time limit to plot in axis

% determine time vector
timemin = zeros(1,840);
timemin(1:30) = 1;

for tt=2:28
    timemin((tt-1)*30:(tt-1)*30+30) = tt;
end

y = 1:840;
x = mean(ddata,2);
z = std(ddata,[],2)/sqrt(numsubs);

x1 = mean(adata,2);
z1 = std(adata,[],2)/sqrt(numsubs);

x2 = mean(bdata,2);
z2 = std(bdata,[],2)/sqrt(numsubs);

x3 = mean(ldata,2);
z3 = std(ldata,[],2)/sqrt(numsubs);

x4 = mean(zdata,2);
z4 = std(zdata,[],2)/sqrt(numsubs);
fig = figure;

clrs = get(gca,'colororder');

% set(fig,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);%  yyaxis left
yyaxis left

[l,p] = boundedline(y, x, z, y, x1, z1, y, x2, z2, y, x3, z3,'alpha', 'transparency', 0.20);

for jj=1:4
set(l(jj), 'linewidth', 1.5,'color',clrs(jj,:));
set(p(jj), 'facecolor', clrs(jj,:));
end

xl = xlabel('Time (minutes)', 'FontWeight', 'bold','FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',12);
yl = ylabel({'\DeltaEEG'; '(normalized)'}, 'FontWeight', 'bold', 'FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',12);

yyaxis right

[l,p] = boundedline(y, x4, z4,'alpha', 'transparency', 0.20);

set(l(1), 'linewidth', 1.5,'color',clrs(5,:));
set(p(1), 'facecolor', clrs(5,:));
ax = gca;
ax.YAxis(1).Color = 'k';

ax.YAxis(2).Color = 'k';
ylim([0 10]);
yl = ylabel('\DeltaIntensity', 'FontWeight', 'bold', 'FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',12,'color',clrs(5,:));
set(gca, 'linewidth', 2,'Fontsize',24, 'Box', 'on')
set(gcf, 'color', [1 1 1],'Position', [1000 918 858 356]);
%     set(gca, 'xlim', lim,'ylim',ylim); % here change the frequencies to explore
set(gca, 'xlim', lim); % here change the frequencies to explore
set(gca,'XTick',[0:120:840] );
set(gca,'XTickLabel',[-8:4:20]);
set(gcf,'Position',[706 688 533 215])
% But put labales on the first one
legend('Delta', 'Alpha','Beta','LZs','fontsize',14, 'Box', 'on', 'FontName', 'Arial', 'location', 'northwest', 'LineWidth', 1,'AutoUpdate','off');
left_color = [0, 0, 1];
right_color =  	clrs(5,:);
% yyaxis left
oldylim = ylim;
hold on
yli = ylim;
    plot([240 240],[yli(1) yli(2)],'k.-.','linewidth', 3)
    ylim(oldylim);     
          
    hold off
    

export_fig(sprintf('TimeCourseSignificantEEGIRASA.png'),'-m2.5')