%% Very similar to dRSFC script we only now do dynamic analyses per edge as oposed to ROI with GFC
% we need already parcelated data from dRSFC script

slidwin = 22;
sigmas = 3;

cd('/Users/christophertimmermann/Documents/Imaging_fMRI/Leor_Preproc/RS_NoHighpass_motionD_Detrend1/DMT_LongSchaeff100Ext')
% cd('/Users/christophertimmermann/Documents/Imaging_fMRI/Leor_Preproc/RS_lowpass/CompleteSchaeff100Ext')
% cd('/Users/christophertimmermann/Documents/Imaging_fMRI/Leor_Preproc/RS_NoHighpass_motionD_Detrend1_Ratings/LongSchaeff100Ext')
allD = dir( '*DMT.mat');
allP = dir( '*PCB.mat');
allfiles = [allD;allP];

 vis = [1:9      51:58];
 sm = [10 :15    59:66]; 
 dan = [16:23    67:73];
 van = [24:30    74:78];
 limbic = [31:33 79:80];
 fp = [34 :37 81: 89];
 dmn = [38: 50 90: 100];
 subc = [101:112];
 
 
for ii=1:size(allD,1)
    file = load(allD(ii).name);
    file.BOLD_AAL = [file.BOLD_AAL(vis,:) ; file.BOLD_AAL(sm,:);file.BOLD_AAL(dan,:);file.BOLD_AAL(van,:);file.BOLD_AAL(limbic,:);file.BOLD_AAL(fp,:);file.BOLD_AAL(dmn,:);file.BOLD_AAL(subc,:)];
    final = [zeros(size(file.BOLD_AAL,1),slidwin/2) file.BOLD_AAL zeros(size(file.BOLD_AAL,1),slidwin/2)];
%         final = [zeros(size(BOLD_AAL,1),slidwin/2) BOLD_AAL zeros(size(BOLD_AAL,1),slidwin/2)];
    DMTCorrTime(:,:,:,ii) = tapered_sliding_window(final',slidwin, sigmas);
    clear file final
    file = load(allP(ii).name);   
    file.BOLD_AAL = [file.BOLD_AAL(vis,:) ; file.BOLD_AAL(sm,:);file.BOLD_AAL(dan,:);file.BOLD_AAL(van,:);file.BOLD_AAL(limbic,:);file.BOLD_AAL(fp,:);file.BOLD_AAL(dmn,:);file.BOLD_AAL(subc,:)];
    final = [zeros(size(file.BOLD_AAL,1),slidwin/2) file.BOLD_AAL zeros(size(file.BOLD_AAL,1),slidwin/2)];
    PCBCorrTime(:,:,:,ii) = tapered_sliding_window(final',slidwin, sigmas);
    clear file final
end
    

%% Obtain FCD
% get identity diagonal to zero

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


%% Intensity ratings for cleaned subs

clear alphareg
truerat=0;
ranking=0;
ovrlap=1;
regressor = 1; % 1= Intensity, 2 = Body, 3 = Visual, 4 = Emo, 5 = Plasmalevels

if regressor>1
    % load microphenomenology

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
        %      outbp = ft_preproc_lowpassfilter(intintp',0.5, 0.3);
            xBF = spm_get_bf(struct('dt',2,'name','hrf')); % convolve with HRF
            alphareg(:,i) = spm_Volterra(struct('u',outbp,'name',{{'task'}}),xBF.bf);
        end
        fxx = alphareg;
        for ii=1:size(DMTCorrTime,4)
            alphareg(:,ii) = fxx;
        end
else
    
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
    %         intintp(1:30) = [];
    %         intintp(end+1:end+30) = 0;
    %         intintp(1:30) = 0;

            outbp = ft_preproc_lowpassfilter(intintp',0.5,  0.08);
        %      outbp = ft_preproc_lowpassfilter(intintp',0.5, 0.3);
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

figure;
plot(mean(RegWind,2))

%% Do pairwise MLE stats
% hold on
% Fill with nans half of the tables so that we then loop through them
FinalCorrMatTime = DMTCorrTime - PCBCorrTime;

fmrivar = FinalCorrMatTime;
%Stack all subjects' Intensity in single vector
start = 0; % Minute to start correlation
ends = 28; % Minute to end correlation
mins = ends-start; % determine how many mins for correlations

All = [];
count = 1;

if start >0  
FCDfinal = fmrivar(:,(start*30):(ends*30)-1,:); % adjust to whatever time you are interested
alpharegfinal = RegWind((start*30):(ends*30)-1,:);
else
   FCDfinal = fmrivar;
   alpharegfinal = RegWind;
end

time = size(FCDfinal,3);
Points = time; % DONT FORGET TO CHANGE THIS BACK

   for i = 1 : size(FCDfinal,4)
	All(count:count+Points-1,1) = repmat(i,Points,1);
	count = count+Points;
   end
   
   

for ss=1:size(FCDfinal,4)
    for pp=1:size(FCDfinal,2)
        FCDfinal(pp,pp:end,:,ss) = nan;
    end
    corrsubs{ss} = FCDfinal(:,:,:,ss);
end

% reduce dimensionality of matrix so that we don't struggle with 4D stuff
for s=1:length(corrsubs) % loop through each sub
    s
    corrmats = corrsubs{s};
    for t=1:size(corrmats,3) % loop through each tp

        tpcorr = corrmats(:,:,t);
        arrcorr = reshape(tpcorr,1,[]);% reshape p vaues into an array for FDR
        gdx = find(~isnan(arrcorr)); %find the non-nanvalues of the simmettrical thing
        B = arrcorr'; % delete nans from array
        B = arrcorr(~isnan(arrcorr))';
        rvaltimesub(:,t,s) = B;
    end
end


finalrvaltimes = rvaltimesub(:,1:840,:); % reduce time to 16 mins
RegWind2 = alpharegfinal(1:840,:);
Alpha = reshape(alpharegfinal,1,[]);% reshape p vaues into an array for FDR
% All(:,2) = Alpha';

%select different lme's to be run

mod{1} = 'Edge ~ Alpha + (1 | Subject) + (-1 + Alpha | Subject)';
% mod{1} = 'Edge ~ 1 + Alpha  + (1 | Subject)';
% mod{3} = 'Alpha ~  Edge +  Subject*Edge + (1 | Subject)';
% mod{4} = 'CorrAlpha*Edge ~ drug + (1 | Edge) + (1 | Subject)';

sizefinalrvaltimes = size(finalrvaltimes,1);

 for mm=1:length(mod)
     model = mod{mm};
     
    fprintf('*** Running on model %d \n', mm)

              
 % Run LME looping through each edge  
             
     parfor ii=1:sizefinalrvaltimes % loop through edges 
        fprintf('*** Running on edge %d \n', ii)
        newmat = squeeze(finalrvaltimes(ii,:,:)); % work with matrix so as to get an array for a looped edge with all times and subs

          edge = reshape(newmat,1,[]);% reshape p vaues into an array for FDR     
%           All(:,3) = edge';

         tbl = table(All,Alpha', edge', 'VariableNames',{'Subject','Alpha','Edge'});

       %Model with slope and intercept free Beta fixed subject random
%            lme = fitlme(tbl, 'Delta ~ Intensity + (1 | Subject) + (-1 + Intensity | Subject)');
       lme = fitlme(tbl, model);

       tstat(ii,1) =  lme.Coefficients(2,4);
       pval(ii,1) = lme.Coefficients.pValue(2);
       rsqrd(ii,1) = lme.Rsquared.Adjusted;
       betasInt(ii,1) = lme.Coefficients(2,2);           
       edge = [];
     end
            tvl(:,mm) = double(tstat); pvl(:,mm) = double(pval);rvl(:,mm) = double(rsqrd); bvl(:,mm) = double(betasInt);   
           
end
   
   %% Correct across models
   
   aal=0;
   schaef100=1;
   load('mycoolwarm.mat')
   
   colors = cbrewer('div', 'RdBu', 64);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);


   mm=1;
%    for ii=1:length(mod)
   [thr2 pcor adj2] = fdr(pvl(:,mm));
%       pcor = pvl(:,mm)*4005 ; % Bonferroni correction
   

    arrcorr(gdx) = pcor; % fill array with correct values in correct position
    FinalPvals = reshape(arrcorr,112,112); % reshape again
    arrcorr(gdx) = tvl(:,mm);
    FinalTvals = reshape(arrcorr,112,112); % reshape again

   
for jj=1:size(FinalTvals,2)
FinalTvals(jj,jj:end) = 0;
end

[r c] = find(FinalPvals<0.05);
sigp = [r c];
simsigp  =  [c r];
% allsig = [sigp;simsigp];


positive = zeros(90,2);
negative = zeros(90,2);

% separate positive from negative correlations
for ii=1:length(simsigp)
    if FinalTvals(simsigp(ii,2),simsigp(ii,1))>0
       positive(ii,:) = simsigp(ii,:);
       negative(ii,:) = nan; 
    elseif FinalTvals(simsigp(ii,2),simsigp(ii,1))<0
        negative(ii,:) = simsigp(ii,:); 
        positive(ii,:) = nan; 
    end
end

% remove nans and zeros remaining
positive(any(isnan(positive),2),:) = []; % remove nans
negative(any(isnan(negative),2),:) = [];
positive( ~any(positive,2), : ) = [];  % remove zeros
negative( ~any(negative,2), : ) = []; 


for ii=1:length(positive)
    FinalTvals(positive(ii,1),positive(ii,2))=abs(max(FinalTvals(:)));
end

for ii=1:length(negative)
    FinalTvals(negative(ii,1),negative(ii,2))=abs(max(FinalTvals(:)))*-1;
end


 

figure;
imagesc(FinalTvals,[(abs(max(FinalTvals(:)))*-1) abs(max(FinalTvals(:)))])
colormap(colors)
% imagesc(FinalTvals2,[-8 8])

% imagesc(FinalTvals,[-7 7])

if aal==1
    
    name = {'Frontal','Orbital','Limbic','Occipital','Parietal','Subcortical','Temporal'};
    locticks = [(20-0)/2 ((28-20)/2)+20 ((42-28)/2)+28 ((56-42)/2)+42 ((70-56)/2)+56 ((78-70)/2)+70 ((90-78)/2)+78]; % location for where names will be
    colormap(mycoolwarm)
    % set(gca,'XTick',[1:12])
    % set(gca,'YTick',[1:12])
    % set(gca, 'XTickLabel', name); % set x-axis labels
    % set(gca, 'YTickLabel', name); % set y-axis labels
    % allfigs = allchild(gcf);
    set(gca, 'fontsize',15, 'Box', 'on','YTickLabel',[],'XTickLabel',[]);
    % caxis([-1 1])
    set(gcf, 'color', [1 1 1]);   
    gbh = colorbar
    set(gbh,'YTick',[-3:3:3],'fontsize',15,  'FontWeight', 'bold')
    set(gcf, 'color', [1 1 1],'position', [552 818 338 239]);
    set(gca,'YTick',locticks, 'XTick',locticks,'XTickLabel', name,'YTickLabel', name,'XTickLabelRotation',45)
    % 0-20 frontal
    % 21-28 orbital
    % 29-42 limbic
    % 43-56 occipital
    % 57 -  70 parietal
    % 71 - 78 subcortical
    % 79- 90 temporal
    hold on
    line([20 20],[0 90],'Color', [0.25 0.25 0.25 0.5], 'linewidth',2)
    line([28 28],[0 90],'Color', [0.25 0.25 0.25 0.5], 'linewidth',2)
    line([42 42],[0 90],'Color', [0.25 0.25 0.25 0.5], 'linewidth',2)
    line([56 56],[0 90],'Color', [0.25 0.25 0.25 0.5], 'linewidth',2)
    line([70 70],[0 90],'Color', [0.25 0.25 0.25 0.5], 'linewidth',2)
    line([78 78],[0 90],'Color', [0.25 0.25 0.25 0.5], 'linewidth',2)

    line([0 90],[20 20],'Color', [0.25 0.25 0.25 0.5], 'linewidth',2)
    line([0 90],[28 28],'Color', [0.25 0.25 0.25 0.5], 'linewidth',2)
    line([0 90],[42 42],'Color', [0.25 0.25 0.25 0.5], 'linewidth',2)
    line([0 90],[56 56],'Color', [0.25 0.25 0.25 0.5], 'linewidth',2)
    line([0 90],[70 70],'Color', [0.25 0.25 0.25 0.5], 'linewidth',2)
    line([0 90],[78 78],'Color', [0.25 0.25 0.25 0.5], 'linewidth',2)
    
    elseif schaef100==1
        
 vis = [1:9      51:58]; 
 sm = [10 :15    59:66]; % 17 start
 dan = [16:23    67:73]; % starts at (6+8+17) = 31
 van = [24:30    74:78]; % starts at (8+7+31) = 46
 limbic = [31:33 79:80]; % (7 + 5 + 46) = 58
 fp = [34 :37 81: 89]; % (3 + 2 + 58) = 63
 dmn = [38: 50 90: 100]; % (4 + 10 + 63) = 77
 subc = [101:112]; % 3 + 11 + 77 = 91
 
 % determine the points where the line goes
 ll=[length(vis) (length(vis)+length(sm)) (length(vis)+length(sm)+ length(dan)) (length(vis)+length(sm) + length(dan)+length(van))  (length(vis)+length(sm)+ length(dan)+length(van)+length(limbic)) (length(vis)+length(sm)+ length(dan)+length(van)+length(limbic)+length(fp)) (length(vis)+length(sm)+ length(dan)+length(van)+length(limbic)+length(fp)+length(dmn))];

 
    name = {'VIS','SM','DAN','SAL','LIM','FP','DMN','SC'};
    set(gca, 'fontsize',15, 'Box', 'on','YTickLabel',[],'XTickLabel',[]);
% caxis([-1 1])
locticks = [(ll(1)-0)/2 ((ll(2)-ll(1))/2)+ll(1) ((ll(3)-ll(2))/2)+ll(2) ((ll(4)-ll(3))/2)+ll(3) ((ll(5)-ll(4))/2)+ll(4) ((ll(6)-ll(5))/2)+ll(5) ((ll(7)-ll(6))/2)+ll(6)   ((112-ll(7))/2)+ll(7)]; % location for where names will be

set(gcf, 'color', [1 1 1]);   
gbh = colorbar
set(gbh,'YTick',[ceil(-1*(abs(max(FinalTvals(:)))))+1:floor(abs(max(FinalTvals(:))))-1:floor(abs(max(FinalTvals(:))))-1],'fontsize',15,  'FontWeight', 'bold','TickDirection' , 'out')
ylabel(gbh, 'T','Position', [0 0 0])
    set(gcf, 'color', [1 1 1],'position', [552 786 325 271]);
set(gca,'YTick',locticks, 'XTick',locticks,'XTickLabel', name,'Linewidth',2,'YTickLabel', name,'XTickLabelRotation',45)


hold on

    for s=1:length(ll)

    line([ll(s) ll(s)],[0 1000],'Color', [0.25 0.25 0.25 0.5], 'linewidth',2)
    line([0 1000],[ll(s) ll(s)],'Color', [0.25 0.25 0.25 0.5], 'linewidth',2)

    end
end


% imagesc automatically flips the y-axis so that the smallest values go on
% top. Set this right if we want the origin to be in the left bottom
% corner.
set(gca, 'ydir', 'normal');
axis square;
 
% add the colorbar, make it prettier
% handles = colorbar;
% handles.TickDirection = 'out';
% handles.Box = 'off';
% handles.Label.String = 'T-statistic';
% drawnow;

% box on
% h=gca
% gca.Linewidth = 5
 export_fig(sprintf('IntensityVsCorrsMat.png'),'-m2.5')
