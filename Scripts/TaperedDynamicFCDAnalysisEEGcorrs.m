%% Do Dynamic FCD agianst EEG metrics using a Tapaered Window

slidwin = 22;
sigmas = 3;

cd('/Users/christophertimmermann/Documents/Imaging_fMRI/Leor_Preproc/RS_lowpass/CompleteSchaeff100Ext')

allD = dir( '*DMT.mat');
allP = dir( '*PCB.mat');
allfiles = [allD;allP];

for ii=1:size(allD,1)
    file = load(allD(ii).name);
    final = [zeros(size(file.BOLD_AAL,1),slidwin/2) file.BOLD_AAL zeros(size(file.BOLD_AAL,1),slidwin/2)];
    DMTCorrTime(:,:,:,ii) = tapered_sliding_window(final',slidwin, sigmas);
    clear file final
    file = load(allP(ii).name);   
    final = [zeros(size(file.BOLD_AAL,1),slidwin/2) file.BOLD_AAL zeros(size(file.BOLD_AAL,1),slidwin/2)];
    PCBCorrTime(:,:,:,ii) = tapered_sliding_window(final',slidwin, sigmas);
    clear file final
end

% FCD for correlations

DMTCorrTime(:,:,:,badsubs) = [];
PCBCorrTime(:,:,:,badsubs) = [];

  
for ss=1:size(DMTCorrTime,4) %loop through subs
    for tt=1:size(DMTCorrTime,3) %loop through time       
        for aa=1:size(DMTCorrTime,1) %loop through area
        PCBCorrTime(aa,aa,tt,ss) = 0;
        DMTCorrTime(aa,aa,tt,ss) = 0;
        end
    end
end
        
% Normalize
PCBCorrTime = atanh(PCBCorrTime);
DMTCorrTime = atanh(DMTCorrTime);

% Obtain FCD
for ss=1:size(DMTCorrTime,4) %loop through subs
    for tt=1:size(DMTCorrTime,3) %loop through time
        for aal=1:size(DMTCorrTime,1) %loop through area
            % Determine mean of correlations. ie get FCD
            FCDDMTtime(aal,tt,ss) = nanmean(DMTCorrTime(aal,:,tt,ss),2 );
            FCDPCBtime(aal,tt,ss) = nanmean(PCBCorrTime(aal,:,tt,ss),2 );
        end
    end
end

FCDdiffTime = FCDDMTtime - FCDPCBtime;

%% Dynamic GFC per network

xdata = FCDDMTtime;
ydata = FCDPCBtime;

 vis = [1:9      51:58];
 sm = [10 :15    59:66]; 
 dan = [16:23    67:73];
 van = [24:30    74:78];
 limbic = [31:33 79:80];
 fp = [34 :37 81: 89];
 dmn = [38: 50 90: 100];
 subc = [101:112];
 
%  nets= {vis, sm , dan ,van , limbic , fp, dmn};
  nets= {vis, sm ,limbic , dan ,van ,  fp, dmn};

 for tt=1:length(nets)
     net = nets{tt};

    GlobalnetDMT(:,tt,:) = nanmean(xdata(net,:,:),1);
    GlobalnetPCB(:,tt,:) = nanmean(ydata(net,:,:),1);
    WithinDelta(:,:,tt) =  squeeze(nanmean(nanmean(xdata(net,net,:,:),1),2)) - squeeze(nanmean(nanmean(ydata(net,net,:,:),1),2));
 end
  
 Globaldiff = GlobalnetDMT -  GlobalnetPCB;
 
 Globaldiff = permute(Globaldiff,[2,1,3]);
 GlobalnetDMT = permute(GlobalnetDMT,[2,1,3]);
 GlobalnetPCB = permute(GlobalnetPCB,[2,1,3]);

%% Generate regressor
% load the mask here to see the channels from the corresponding regressors
load('/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/Analysis/IRASA/5minData/maskdiffmaxpaired.mat')
clear regressor Reg*

bands = {'Delta','Theta','Alpha','Beta','Gamma'};
areas = {'Occipital','Parietal','Temporal','Central','Frontal'};

% for bb=1:length(bands)
    for bb=1

    for aa=1:length(areas)

%determine options and parameters
printfig = 0; % Print figure?
fmri = 1; % fMRI variable for analysis 1 = FCD per area, 2 = FCD per network
generalregs = 0;
ovrlap = 1;
delaylag = 0;
noconv = 0;
clearbad = 1; % clear bad EEG trials?
regre =bb; %determine regressor of interest. 1 alpha, 2 beta, 3 lz. 3 = gamma and 4 = lz in regressorsused =5. If IRASA is selected 1 delta, 2 theta, 3 alpha, 4 beta, 5 gamma; If travelling waves 1= BW travel, 2 = FW travel
diffz = 0; %use difference between DMT and placebo?
irasa = 0;
del = 0; % how many TRs will the lag have?
bonf = 0;
convolved = 2; % use here if using already convolved regressors
areaa = aa; %applies for general regressors  0 =  old regs, 1=occipital, 2=parietal, 3=temporal, 4=central, 5=frontal
regressorsused = 8; % 8 = IRASA per area channels, % 9 = LZc across channs

slidwind=22;
padded = slidwind/2;

load('/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/Analysis/badtrs4regs.mat')

    if regressorsused==8 && areaa ==1
%     load('/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/Analysis/regressorsHilbAlphBetGammaLzsInterpscrubbedIRASAOccipital.mat')
    load('/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/Analysis/RegressorsInterpscrubbedIRASA_Occipital.mat')
    elseif regressorsused==8 && areaa ==2
%     load('/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/Analysis/regressorsHilbAlphBetGammaLzsInterpscrubbedIRASAParietal.mat')
    load('/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/Analysis/RegressorsInterpscrubbedIRASA_Parietal.mat')
    elseif regressorsused==8 && areaa ==3
%     load('/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/Analysis/regressorsHilbAlphBetGammaLzsInterpscrubbedIRASATemporal.mat')
    load('/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/Analysis/RegressorsInterpscrubbedIRASA_Temporal.mat')
    elseif regressorsused==8 && areaa ==4
%     load('/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/Analysis/regressorsHilbAlphBetGammaLzsInterpscrubbedIRASACentral.mat')
    load('/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/Analysis/RegressorsInterpscrubbedIRASA_Central.mat')
    elseif regressorsused==8 && areaa ==5
%     load('/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/Analysis/regressorsHilbAlphBetGammaLzsInterpscrubbedIRASAFrontal.mat')
    load('/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/Analysis/RegressorsInterpscrubbedIRASA_Frontal.mat')
    elseif regressorsused==9
    load('/Users/christophertimmermann/Documents/MATLAB/DMT_Imaging/EEG_analysis/ChrisClean/Processed/Analysis/regressorsTravellingWavesInterp.mat')
  end
        

if diffz==1
    regressor = regdiff;
else
    regressor = regDMT;
end

addpath '/Users/christophertimmermann/Documents/Imaging_programs/spm12'

for ii=1:size(regressor,1)
    
     if regressorsused==10
         out = regressor(ii,:)';
     else
     out = regressor(ii,:,regre)';
     end
  
    if noconv==1
        alphareg(:,ii) = out;
    else
        xBF = spm_get_bf(struct('dt',2,'name','hrf')); % convolve with HRF
        alphareg(:,ii) = spm_Volterra(struct('u',out,'name',{{'task'}}),xBF.bf);
    end
    
    clear out
end


 if ovrlap==0
    for ii=1:size(alphareg,2)
        reg = alphareg(:,ii);
        tps = [1:30:840];
        for pp=1:length(tps)
            tim = tps(pp);
                AvgRegWing   =  nanmean(reg([tim:(tim+29)])); % Extract time window of interest
                RegWind(pp,ii) = AvgRegWing; % do correlation matrix, with time on 3rd dimension, subjects in the 4th
        end
    end
 else
            
  slidwin=22;
  sigmas=3;
padded = slidwin/2;


alpharegpad = nan(size(alphareg,1)+slidwin,size(alphareg,2)) ;% initialize padded variable
alpharegpad(padded+1:(size(alpharegpad,1)-padded),:) = alphareg;

% if already convolved, don't convolve again
if convolved==1
    RegWind = alphareg;
elseif convolved==2 % not previously convolved but using the sliding window. Requires interpolated data
        alpharegpad = nan(size(alphareg',1),size(alphareg',2)+ slidwin) ;% initialize padded variable
        alpharegpad(:,padded+1:(size(alpharegpad',1)-padded)) = alphareg';
        [c , RegWind] = tapered_sliding_window(alpharegpad',slidwin, sigmas);
        RegWind = RegWind';

 end
 end


if clearbad==1
    if regressorsused==9
        badregDMT(:,badsubs2) = []; 
    end
            cleanRegWind = RegWind; % clear bad segments

    RegWind(find(badregDMT)) = nan; % clear bad segments
    else
       
end

if delaylag==1
    RegWind2 = nan(del,size(RegWind,2));
    RegWind2 = [RegWind2 ;RegWind];
    RegWind = RegWind2;
    if size(RegWind,1)>840
    RegWind(841:end,:) = [];
    end
end

if regressorsused~=9
    RegWind(:,badsubs2)=[]; %remove guys with bad EEG data
end

inf = (237-slidwind/2:240+slidwind/2); % we eliminate the first 30 seconds
toreplotRegWind = RegWind(inf,:) ;
RegWind(inf,:) = nan;


% figure;

% plot(nanmean(RegWind,2))

clearvars -except fmri cleanRegWind printfig bands bb areas aa RegWind FCD* *CorrTime badsubs* generalregs regre diffz irasa del bonf slidwin sigmas inf toreplot* irasa areaa regressorsused Global*

% MLE for each correlation point in time



if fmri==1
    if diffz==1
    FCDinterest = FCDdiffTime;
    else 
        FCDinterest = FCDDMTtime;
    end
elseif fmri==2
    if diffz==1
    FCDinterest = Globaldiff;
    else
    FCDinterest = GlobalnetDMT;
    end
end




FCDinterest(:,:,badsubs2) = []; %remove guys with bad EEG data

% LME analysis.
%Stack all subjects' Intensity in single vector
start = 0; % Minute to start correlation
ends = 28; % Minute to end correlation
mins = ends-start; % determine how many mins for correlations
time = (30*mins);
clear FCstring FCd FCDfinal All

if start >0  
    if ovrlap==0
        FCDfinal = FCDinterest(:,(start):(ends),:); % adjust to whatever time you are interested
alpharegfinal = RegWind((start):(ends),:);
    else
FCDfinal = FCDinterest(:,(start*30):(ends*30)-1,:); % adjust to whatever time you are interested
alpharegfinal = RegWind((start*30):(ends*30)-1,:);
    end
else
   FCDfinal = FCDinterest;
   alpharegfinal = RegWind;
end

Points = size(FCDfinal,2); % DONT FORGET TO CHANGE THIS BACK
All = [];
count = 1;
   for i = 1 : size(FCDfinal,3)
	All(count:count+Points-1,1) = repmat(i,Points,1);
	count = count+Points;
   end



% rerrange matrix so that each area is a column to loop through LMEs   
for aa=1:size(FCDfinal,1) % loop through area 
    FCd = squeeze(FCDfinal(aa,:,:));% 1 area at a time
    FCstring(:,aa) = reshape(FCd,1,[]);
end


Intensity = reshape(alpharegfinal,1,[]);% reshape p vaues into an array for FDR
All(:,2) = Intensity';
All(:,3) = 1:size(All,1);

%select different lme's to be run

% mod{1} = 'FCD ~ Intensity + (1 | Subject) + (-1 + Intensity | Subject)';
mod{1} = 'Intensity ~ FCD + (1 | Subject) + (-1 + FCD | Subject)';


 for mm=1:length(mod)
     model = mod{mm};
     
    fprintf('*** Running on model %d \n', mm)

 % Run LME looping through each AAl area FCD value  
             
     for ii=1:size(FCstring,2) % loop through AAL area 
        fprintf('*** Running on area %d \n', ii)
        
         All(:,4) = FCstring(:,ii); % select AAL area
         tbl = table(All(:,1),All(:,2), All(:,3), All(:,4), 'VariableNames',{'Subject','Intensity','Time','FCD'});
       %Model with slope and intercept free Beta fixed subject random
       lme = fitlme(tbl, model);

       tstat(ii,1) =  lme.Coefficients(2,4);
       pval(ii,1) = lme.Coefficients(2,6);
%          pval(ii,1) = lme.Coefficients.pValue(2);
       rsqrd(ii,1) = lme.Rsquared.Adjusted;
       betasInt(ii,1) = lme.Coefficients(2,2);           
       All(:,4) = [];
         
     end
            tvl(:,mm) = double(tstat); pvl(:,mm) = double(pval);rvl(:,mm) = double(rsqrd); bvl(:,mm) = double(betasInt);              
end
   
% Correct across models

mm=1;

if bonf==1
    pcor=pvl*112;
else
[thr2 pcor adj2] = fdr(pvl(:,mm));
end
sigaal = find(pcor<0.05)

clearvars -except fmri cleanRegWind printfig bands bb areas aa RegWind FCD* *CorrTime badsubs* generalregs regre diffz irasa del bonf slidwin sigmas inf toreplot* *vl sigaal areaa regressorsused lme Global* pcor

% Plot networks

if fmri==2
    colorst = linspecer(5);
    y= [tvl];

    figure;
    name = {'VIS','SM','LIM','DAN','SAL','FP','DMN'};

    x = 1:7;

    positiveIndexes = find(y > 0);
    negativeIndexes = find(y <  0);

    b = bar(y,'EdgeColor','k','LineWidth',2, 'facecolor', 'flat');
    % b.FaceColor = [.5 0 .5];
    % b(1).FaceColor = [.5 .5 .5];
     hold on
    clr = nan(length(x),3);
    for ii= 1:length(positiveIndexes)
        idx = positiveIndexes(ii);
        clr(idx,:) = colorst(2,:);
    end

    for ii= 1:length(negativeIndexes)
        idx = negativeIndexes(ii);
        clr(idx,:) = colorst(1,:);
    end

    b.CData = clr;


    set(gca, 'XTick',x, 'XTickLabel', name,'Fontsize',24,'XTickLabelRotation',45); % set x-axis labels
    ylabel('T','Fontsize',24)

    % ylabel('\beta','Fontsize',20,'FontWeight', 'bold')


    % h = legend( 'Remainers','Converted','AutoUpdate','off','Fontsize',20)

    % legend( 'Remainers','Converted','AutoUpdate','off')
    set(gcf, 'color', [1 1 1],'Position',[1000 1026 471 312]);    
    sigs = find(pcor<0.05);


    % add significance stuff with the amazing function sigstar

    for tt = 1:length(sigs)
        rr=sigs(tt);
        sigline = [rr rr];
          H = sigstar({sigline'},pcor(rr))
            set(H(1),'LineWidth',0.0000000001)
            set(H(2),'Fontsize',27)

    end
        set(gca, 'fontsize',25, 'Box', 'on');


    %     if min(bvl)>0
    %             ylim([min(bvl)-1  max(bvl)+2  ])
    %     else
    % %         set(gca, 'yDir', 'reverse')

    % if isempty(sigs)==1
    %     if mean(tvl)>0
    %         ylim([0 4 ])
    %     else
    %         ylim([-4 0])
    %     end
    % end

     ylim([-5 5])

    if printfig==1
    % determine name to plot
    if regressorsused==8
        yname = [bands{regre} areas{areaa} ];
    elseif regressorsused==10
        yname = 'LZc';
    elseif regressorsused==9
        if bb==1
        yname = 'BWTravWaves';
        elseif bb==2
            yname = 'FWTravWaves';
        end
    end
    export_fig(sprintf([yname 'Bar.png']),'-m2.5')

    end
     hold off


    close all
    
else


    % Replot in AAL
    aal90=0;
    schaef100=1;

    if aal90==1
    path2mask = '/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/CorrMatrix/roi_mnifsl(2).hdr';
    elseif schaef100==1
            path2mask = '/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/CorrMatrix/Schaeff1007NExt.nii.gz'; % path to Shcaeffer 100 areas extended
    end
    Mask = load_nii(path2mask);

    Maskfcd = Mask;

    Maskfcd.img = zeros(size(Maskfcd.img));

    for tt =1 :length(sigaal) % For each significant aal
        sal = sigaal(tt);
             idx = find(Mask.img == sal); % obtain 3d index for each significant aal.       
              Maskfcd.img(idx) = tvl(sal);
              clear idx
    end

    % determine name to plot
    if regressorsused==8
        yname = ['FCD' bands{regre} areas{areaa} '_Irasa' ];
    elseif regressorsused==10
        yname = 'FCDLZc';
    elseif regressorsused==9
        if areaa==1
        yname = 'BWTravWaves';
        elseif areaa==2
            yname = 'FWTravWaves';
        end
    end


    rmpath '/Users/christophertimmermann/Documents/Imaging_programs/spm12'
    rmpath  '/Users/christophertimmermann/Documents/MATLAB/fieldtrip-20200130'
    % cd  '/Users/christophertimmermann/Documents/Imaging_fMRI/Leor_Preproc/RS_NoHighpass_motionD_Detrend1/DMT_LongSchaeff100Ext/EEG_fMRI_nii_results'
    save_nii(Maskfcd, [yname '.nii.gz']);

    BrainNet_MapCfg('BrainMesh_ICBM152.nv',[yname '.nii.gz'])

end
            end
    end
    


%% Plot timecourse

% Do whole plot vs Intensity

toplot=5;

DMTintplot = FCDinterest;

RegWind(inf,:) = toreplotRegWind;

xdata = squeeze(mean(DMTintplot(toplot,:,:),1));
% ydata = squeeze(mean(PCBintplot(toplot,:,:),1));
% zdata = zscore(fillmissing(cleanRegWind));

zdata = [ nan(11,14) ; zscore(cleanRegWind(12:830,:)) ;nan(10,14) ];

numsubs = size(xdata,2);

lim = [-29 870]; % select time limit to plot in axis

% determine time vector
timemin = zeros(1,840);
timemin(1:30) = 1;
for tt=2:28
    timemin((tt-1)*30:(tt-1)*30+30) = tt;
end

x = mean(xdata,2);
z = std(xdata,[],2)/sqrt(numsubs);

y = 1:840;


x2 = mean(zdata,2);
z2 = std(zdata,[],2)/sqrt(numsubs);

clrs = linspecer(3);

% get(gca,'colororder');
fig = figure;
set(fig,'defaultAxesColorOrder',[clrs(1,:); clrs(2,:)]);%  yyaxis left

yyaxis left

[l,p] = boundedline(y, x, z ,'alpha', 'transparency', 0.20);

for jj=1
set(l(jj), 'linewidth', 2,'color',clrs(jj,:));
set(p(jj), 'facecolor', clrs(jj,:));
end

xl = xlabel('Time (minutes)', 'FontWeight', 'bold','FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',30);
yl = ylabel({['GFC ' name{toplot}]; '(Fisher Z)'}, 'FontWeight', 'bold', 'FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',30);

yyaxis right

[l,p] = boundedline(y, x2, z2,'alpha', 'transparency', 0.20);

set(l(1), 'linewidth', 2,'color',clrs(2,:));
set(p(1), 'facecolor', clrs(2,:));

yl = ylabel({[bands{regre} ' ' areas{areaa}]; '(z)'}, 'FontWeight', 'bold', 'FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',30);
set(gca, 'linewidth', 2,'Fontsize',28, 'Box', 'on')
set(gcf, 'color', [1 1 1],'Position', [1000 918 858 356]);
%     set(gca, 'xlim', lim,'ylim',ylim); % here change the frequencies to explore
set(gca, 'xlim', lim); % here change the frequencies to explore
set(gca,'XTick',[0:60:840] );
set(gca,'XTickLabel',[-8:2:20]);
% But put labales on the first one
% legend('DMT', 'Placebo','fontsize',24, 'Box', 'off', 'FontName', 'Arial', 'location', 'northeast', 'LineWidth', 10,'AutoUpdate','off');
left_color = [0, 0, 1];
right_color =  	clrs(2,:);

oldylim = ylim;
hold on
% yyaxis left
plot([240 240],[-2 2],'k.-.','linewidth', 3)
ylim(oldylim);

export_fig(sprintf([yname ' ' name{toplot} 'Timecourse.png']),'-m2.5')


