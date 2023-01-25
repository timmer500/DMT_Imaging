%% Dynamic Resting state Functional Connectivity (dRSFC)
% Chris Timmermann 2022

% Do Dynamic GFC agianst Intensity using a Tapaered Window

% Then  use tapered sliding window to do sliding window correlation matrix

slidwin = 22;
sigmas = 3;

cd('/Users/christophertimmermann/Documents/Imaging_fMRI/Leor_Preproc/RS_NoHighpass_motionD_Detrend1/DMT_LongSchaeff100Ext')

% select here parcellated data 
allD = dir( '*DMT.mat');
allP = dir( '*PCB.mat');
allfiles = [allD;allP];

for ii=1:size(allD,1)
    file = load(allD(ii).name);
    final = [zeros(size(file.BOLD_AAL,1),slidwin/2) file.BOLD_AAL zeros(size(file.BOLD_AAL,1),slidwin/2)];
    DMTCorrTime(:,:,:,ii) = tapered_sliding_window(final',slidwin, sigmas); % file provided
    clear file final
    file = load(allP(ii).name);   
    final = [zeros(size(file.BOLD_AAL,1),slidwin/2) file.BOLD_AAL zeros(size(file.BOLD_AAL,1),slidwin/2)];
    PCBCorrTime(:,:,:,ii) = tapered_sliding_window(final',slidwin, sigmas);
    clear file final
end
    

     
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

% Obtain GFC
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


%% Dynamic GFC per network. 
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
 
  nets= {vis, sm ,limbic , dan ,van , fp, dmn};

 for tt=1:length(nets)
     net = nets{tt};

    GlobalnetDMT(:,tt,:) = nanmean(xdata(net,:,:),1);
    GlobalnetPCB(:,tt,:) = nanmean(ydata(net,:,:),1);
 end
  
 Globaldiff = GlobalnetDMT -  GlobalnetPCB;
 
 Globaldiff = permute(Globaldiff,[2,1,3]);
 GlobalnetDMT = permute(GlobalnetDMT,[2,1,3]);
 GlobalnetPCB = permute(GlobalnetPCB,[2,1,3]);


%% Correlate these metrix (integrity or GFC over time with intensity ratings over time or whatever regressor

clear alphareg
printfig = 0;
truerat=0;
ranking=0;
ovrlap=1;
regressor = 1; % 1= Intensity, 2 = Body, 3 = Visual, 4 = Emo, 5 = Plasmalevels
diffz=1; % use difference?
if regressor>1
    % load microphenomenology
        var = Intensity
     if regressor ==2
        var = Bodyalt;
     elseif regressor==3
        var = Visual;
     elseif regressor==4        
         var = Emo;
     elseif regressor==5
  
         var = fillmissing(plasma,'linear')';
     end

    S = mean(var,2);
     if regressor~=5
    S = [zeros(8,1) ; S];
     end
    addpath '/Users/christophertimmermann/Documents/MATLAB/spm12'
    % define query points for interpolation
    qv = 0:28/840:(28-0/840);
    qv(1) = [];

        for i=1:size(S,2)
            if regressor==5
              intintp  = S(:,i)';
            else
            int = S(:,i);
            intintp = interp1(0:28, int, qv);
            end

            outbp = ft_preproc_lowpassfilter(intintp',0.5,  0.08);
            xBF = spm_get_bf(struct('dt',2,'name','hrf')); % convolve with HRF
            alphareg(:,i) = spm_Volterra(struct('u',outbp,'name',{{'task'}}),xBF.bf);
        end
        fxx = alphareg;
        for ii=1:size(FCDdiffTime,3)
            alphareg(:,ii) = fxx;
        end
else
    
    numDMT = xlsread('/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/CorrMatrix/IntensityRatings.xlsx','DMT');

    if truerat==1
    numDMT(numDMT==11)=10;
    end

    numPCB = xlsread('/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/CorrMatrix/IntensityRatings.xlsx','PCB');

    S = numDMT - numPCB;
%         S = numDMT ;


    S = [zeros(20,1) S];

    if ranking==1
        for ss=1:size(S,1)
             data = S(ss,:);
            X = data;
            [temp,S(ss,:)]  = ismember(X,unique(X));
        end
    end

    S(badsubs,:) = [];

    S = S';

    if ovrlap==0
        S(1,:) = [];
        RegWind = S;
    else
        
    addpath '/Users/christophertimmermann/Documents/MATLAB/spm12'

    % define query points for interpolation
    qv = 0:28/840:(28-0/840);
    qv(1) = [];

        for i=1:size(S,2)
            int = S(:,i);
            intintp = interp1(0:28, int, qv);

            outbp = ft_preproc_lowpassfilter(intintp',0.5,  0.08);
            xBF = spm_get_bf(struct('dt',2,'name','hrf')); % convolve with HRF
            alphareg(:,i) = spm_Volterra(struct('u',outbp,'name',{{'task'}}),xBF.bf);
        end
    end

end
padded = slidwin/2;
addpath '/Users/christophertimmermann/Documents/MATLAB/spm12'

alpharegpad = nan(size(alphareg,1)+slidwin,size(alphareg,2)) ;% initialize padded variable
alpharegpad(padded+1:(size(alpharegpad,1)-padded),:) = alphareg;

[c  AvgRegWing] = tapered_sliding_window(alpharegpad,slidwin, sigmas);
RegWind= AvgRegWing'; % do correlation matrix, with time on 3rd dimension, subjects in the 4th

inf = (237-slidwin/2:240+slidwin/2); % we eliminate the first 30 seconds
toreplotRegWind = RegWind(inf,:) ;
RegWind(inf,:) = nan;

%% MLE analysis for each correlation point in time

% LME analysis.
%Stack all subjects' Intensity in single vector
start = 0; % Minute to start correlation
ends = 28; % Minute to end correlation
mins = ends-start; % determine how many mins for correlations


fmri = 1; % fMRI variable for analysis 1 = GFC per area, 2 = GFC per network, 3 = Integrity of networks

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
elseif fmri==3
    if diffz==1
    FCDinterest = Withindiff;
    else
    FCDinterest = WithinDMT;
    end
end

All = [];
count = 1;

if start >0  
FCDfinal = FCDinterest(:,(start*30):(ends*30)-1,:); % adjust to whatever time you are interested
alpharegfinal = RegWind((start*30):(ends*30)-1,:);
else
   FCDfinal = FCDinterest;
   alpharegfinal = RegWind;
end

time = size(FCDfinal,2);
Points = time; % DONT FORGET TO CHANGE THIS BACK

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

mod{1} = 'FCD ~ Intensity + (1 | Subject) + (-1 + Intensity | Subject)';

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
       rsqrd(ii,1) = lme.Rsquared.Adjusted;
       betasInt(ii,1) = lme.Coefficients(2,2);           
       All(:,4) = [];
         
     end
            tvl(:,mm) = double(tstat); pvl(:,mm) = double(pval);rvl(:,mm) = double(rsqrd); bvl(:,mm) = double(betasInt);              
end
   
% Correct across models

mm=1;
[thr2 pcor adj2] = fdr(pvl(:,mm));

sigaal = find(pcor<0.05)

% Plot network-level results (could be either integrity (Fig 3a/b lower
% right) GFC (Fig 3a/b lower left)

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
    
    ylim([-7 7])
    if printfig==1
    if regressor==1
        yname = 'Intensity';
    elseif regressor==5
        yname = 'Plasma DMT';
   
    end
      export_fig(sprintf([yname 'vsGFCnet.png']),'-m2.5')

    end
close all
%% Replot in ROIS significant correlations

aal90=0; 
aal401=0; 
schaef100=1;


if schaef100==1
        path2mask = '/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/CorrMatrix/Schaeff1007NExt.nii.gz'; % path to Shcaeffer 100 areas extended
end
Mask = load_nii(path2mask);

Maskfcd = Mask;

Maskfcd.img = zeros(size(Maskfcd.img));


for tt = 1:112 % For each significant aal
         idx = find(Mask.img == tt); % obtain 3d index for each significant aal.   
         fds{tt} = idx;
          Maskfcd.img(idx) = bvl(tt);
          clear idx
end

[v,loc] = min(Maskfcd.img(:));

Maskfcd.img(find(Maskfcd.img)) = Maskfcd.img(find(Maskfcd.img)) + abs(v) + 0.0000000001;


save_nii(Maskfcd, 'IntensityFCDnonthresh.nii.gz');

rmpath '/Users/christophertimmermann/Documents/MATLAB/spm12'

BrainNet_MapCfg('BrainMesh_ICBM152.nv','IntensityFCDnonthresh.nii.gz')

[r c] = find(bvl<0)
intersect(r,sigaal)

 export_fig(sprintf('BVLDynamicDMTvsPCBNonthreshHSV.png'),'-m2.5')

%% Here for Neurosynth mapping


aal90=0; 
aal401=0; 
schaef100=1;
schaef1000=0;

if aal90==1
path2mask = '/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/CorrMatrix/roi_mnifsl(2).img';
elseif aal401==1
path2mask = '/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/CorrMatrix/ALL_REGS.nii';
elseif schaef100==1
        path2mask = '/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/CorrMatrix/Schaeff1007NExt.nii.gz'; % path to Shcaeffer 100 areas extended
        elseif schaef1000==1
    path2mask = '/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/CorrMatrix/Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm.nii.gz'; % path to Enzo's 401 areas
end
Mask = load_nii(path2mask);

Maskfcd = Mask;

Maskfcd.img = zeros(size(Maskfcd.img));

for tt = 1:length(tvl) % For each significant aal
         idx = find(Mask.img == tt); % obtain 3d index for each significant aal.       
          Maskfcd.img(idx) = tvl(tt);
          clear idx
end

save_nii(Maskfcd, 'FCDDMTvsPCBDynamicNonthreshTapered.nii.gz');


BrainNet_MapCfg('BrainMesh_ICBM152.nv','FCDDMTvsPCBDynamicNonthreshTapered.nii.gz')


%% Do FCD vs intensity plot (Fig 3C)

colors = cbrewer('div', 'RdBu', 64);
colors = flipud(colors); % puts red on top, blue at the bottom

meanFCD = nanmean(FCDdiffTime,3);

% re-arrange so parcellations for ecah network are contiguous
 vis = [1:9      51:58];
 sm = [10 :15    59:66]; 
 dan = [16:23    67:73];
 van = [24:30    74:78];
 limbic = [31:33 79:80];
 fp = [34 :37 81: 89];
 dmn = [38: 50 90: 100];
 subc = [101:112];
 
toplot = zeros(size(meanFCD));

toplot(1:length(dmn),:) = meanFCD(dmn,:);
toplot(length(dmn)+1:length(dmn)+length(fp),:) = meanFCD(fp,:);
toplot(length(dmn)+length(fp)+1:length(dmn)+length(fp)+length(van),:) = meanFCD(van,:);
toplot(length(dmn)+length(fp)+length(van)+1:length(dmn)+length(fp)+length(van)+length(dan),:) = meanFCD(dan,:);
toplot(length(dmn)+length(fp)+length(van)+length(dan)+1:length(dmn)+length(fp)+length(van)+length(dan)+length(limbic),:) = meanFCD(limbic,:);
toplot(length(dmn)+length(fp)+length(van)+length(dan)+length(limbic)+1:length(dmn)+length(fp)+length(van)+length(dan)+length(limbic)+length(sm),:) = meanFCD(sm,:);
toplot(length(dmn)+length(fp)+length(van)+length(dan)+length(limbic)+length(sm)+1:length(dmn)+length(fp)+length(van)+length(dan)+length(limbic)+length(sm)+length(vis),:) = meanFCD(vis,:);
toplot(length(dmn)+length(fp)+length(van)+length(dan)+length(limbic)+length(sm)+length(vis)+1:length(dmn)+length(fp)+length(van)+length(dan)+length(limbic)+length(sm)+length(vis)+length(subc),:) = meanFCD(subc,:);
    

clrs = linspecer(10);

fig = figure;
imagesc(toplot,[-0.35 0.35]);
set(fig,'defaultAxesColorOrder',[[0 0 0]; clrs(10,:)]);%  yyaxis left

yyaxis left
colormap(parula);
hold on

 % determine the points where the line goes
%  ll=[length(vis) (length(vis)+length(sm)) (length(vis)+length(sm)+ length(dan)) (length(vis)+length(sm) + length(dan)+length(van))  (length(vis)+length(sm)+ length(dan)+length(van)+length(limbic)) (length(vis)+length(sm)+ length(dan)+length(van)+length(limbic)+length(fp)) (length(vis)+length(sm)+ length(dan)+length(van)+length(limbic)+length(fp)+length(dmn))];
 ll=[length(dmn) (length(dmn)+length(fp)) (length(dmn)+length(fp)+ length(van)) (length(dmn)+length(fp) + length(van)+length(dan))  (length(dmn)+length(fp)+ length(van)+length(dan)+length(limbic)) (length(dmn)+length(fp)+ length(van)+length(dan)+length(limbic)+length(sm)) (length(dmn)+length(fp)+ length(van)+length(dan)+length(limbic)+length(sm)+length(vis))];

 
%  name = {'VIS','SM','DAN','SAL','LIM','FP','DMN','SC'};
  name = {'DMN','FP','SAL','DAN','LIM','SM','VIS','SC'};

 set(gca, 'fontsize',20, 'Box', 'on','YTickLabel',[]);
 xl = xlabel('Time (minutes)', 'FontWeight', 'bold','FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',24);
  yl = ylabel('Area', 'FontWeight', 'bold','FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',24);

% caxis([-1 1])
locticks = [(ll(1)-0)/2 ((ll(2)-ll(1))/2)+ll(1) ((ll(3)-ll(2))/2)+ll(2) ((ll(4)-ll(3))/2)+ll(3) ((ll(5)-ll(4))/2)+ll(4) ((ll(6)-ll(5))/2)+ll(5) ((ll(7)-ll(6))/2)+ll(6)   ((112-ll(7))/2)+ll(7)]; % location for where names will be

timedisp = -8:2:20;

set(gcf, 'color', [1 1 1]);   
gbh = colorbar
set(gbh,'fontsize',15,  'FontWeight', 'bold','TickDirection' , 'out','YTick',[-0.3 0 0.3]);
ylabel(gbh, 'Fisher Z','Position', [0 0 0])
set(gcf, 'color', [1 1 1],'Position',[451 733 934 324]);
set(gca,'YTick',locticks, 'Linewidth',2,'YTickLabel', name,'XTick',[0:60:840],'XTickLabel',timedisp,'TickLength',[0 0])
xlim([-1 840])

    for s=1:length(ll)

%     line([ll(s) ll(s)],[0 1000],'Color', [0.25 0.25 0.25 0.5], 'linewidth',2)
    line([0 1000],[ll(s) ll(s)],'Color', [0.25 0.25 0.25 0.5], 'linewidth',2)

    end
    
    names = {'DAN','SAL','FP','DMN'};
    
for tp=1:length(toplot)
    if tp==1
%     nets = readtable('/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/FCDtimecourse/DAN_res.csv') ;
        nets = readtable('/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/FCDtimecourse/DMN_res.csv') ;
    ll2 = [1 ll(1)];
    elseif tp==2
%     nets = readtable('/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/FCDtimecourse/Salience_res.csv') ;
        nets = readtable('/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/FCDtimecourse/FP_res.csv') ;

    ll2 = ll([1 2]);
    elseif tp==3
%     nets = readtable('/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/FCDtimecourse/FP_res.csv') ;
 nets = readtable('/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/FCDtimecourse/Salience_res.csv') ;
    ll2 = ll([2 3]);
    elseif tp==4
%     nets = readtable('/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/FCDtimecourse/DMN_res.csv') ;
nets = readtable('/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/FCDtimecourse/DAN_res.csv');
    ll2 = ll([3 4]);
    end
    sigtime = find(nets.pvalue_plus<0.025);
    
    
    line([sigtime(1) sigtime(end)],[ll2(1) ll2(1)],'Color', [0 0 0 0.8], 'linewidth',3)
    line([sigtime(1) sigtime(end)],[ll2(2) ll2(2)],'Color', [0 0 0 0.8], 'linewidth',3)
    line([sigtime(1) sigtime(1)],[ll2(1) ll2(2)],'Color', [0 0 0 0.8], 'linewidth',3)
    line([sigtime(end) sigtime(end)],[ll2(1) ll2(2)],'Color', [0 0 0 0.8], 'linewidth',3)
    


    
end

%          yyaxis left
    plot([240 240],[0 1000],'k.-.','linewidth', 3)
    
    
    RegWind = alpharegpad;
RegWind(1:11,:) = [];
RegWind(end-10:end,:) = [];

 zdata = RegWind;

    numsubs = 14;

    lim = [-29 870]; % select time limit to plot in axis

    % determine time vector
    timemin = zeros(1,840);
    timemin(1:30) = 1;
    for tt=2:28
        timemin((tt-1)*30:(tt-1)*30+30) = tt;
    end



    % y = timemin;
    y = 1:840;

    x2 = mean(zdata,2);
    z2 = std(zdata,[],2)/sqrt(numsubs);
    
      
    yyaxis right
[l2,p2] = boundedline(y, x2, z2,'alpha', 'transparency', 0.20);
    
    set(l2, 'linewidth', 2,'color',clrs(10,:));
    set(p2, 'facecolor', clrs(10,:));
    
    yl = ylabel('\DeltaIntensity', 'FontWeight', 'bold', 'FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',20,'color',clrs(10,:));
 set(gca, 'fontsize',20, 'Box', 'on','color',clrs(10,:));

 set(gca,'YColor',clrs(10,:));

%     title('\DeltaGFC')
     export_fig(sprintf('FCDTimecourse_sig.png'),'-m2.5')


