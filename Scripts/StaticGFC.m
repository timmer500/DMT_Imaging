%% OBTAIN GFC AT EACH PARCEL THEN AT NETWORK LEVEL.
% Chris Timmermann 2022

% determine paths according to each dataset
    parentfolder = '/Users/christophertimmermann/Documents/Imaging_fMRI/Leor_Preproc/RS_lowpass/';
    conditions = {'DMT','PCB'};
        subjects = {'S01','S02', 'S03', 'S06', 'S07', 'S08', 'S09', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S22', 'S23', 'S25'};
    
path2mask = '/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/CorrMatrix/Schaeff1007NExt.nii.gz'; % path to Schaeffer parcellation with extension to subcrotical areas, found now in a filder GFC for you, change the path
  
% FIRST OBRAIN CORRELATION MATRIX
Mask = load_nii(path2mask);

for ss = 1 :length(subjects)
    suj = subjects{ss};
    tic
    for cc=1:length(conditions)
        cond = conditions{cc};
        
               
            path2data = [parentfolder  cond '/post1/MNI152/scrubbed_FD_040/clean_with_intact/' suj '/' cond '_post1_rdsmffms6FWHM_bd_M_V_DV_WMlocal2_modecorr.nii.gz']; % change to your data
            
            Data = load_nii(path2data);
                        

                maxmask = max(Mask.img(:)); % For each ROI specified in the mask
            
                 tic  % Start loop to adapt Volumes to Mask
                 for jj = 1:maxmask % For each ROI specified in the mask
                    for t=1:size(Data.img,4)  % for each timepoint in the data                        
                        MatT = double(Data.img(:,:,:,t));
                           Avg_Signal(jj,t) = mean(MatT(find(Mask.img==jj)));
                        clear MatT
                    end              
                end
                toc
                
                CorrMat(:,:) = corr( Avg_Signal(:,:)'); % do correlation matrix
                
                for k=1:size(CorrMat,1) % make identity diagnoal 0 instead of 1
                    CorrMat(k,k) = 0;
                end
                
                % Save outcome variables according to respective data
                % selected
                
                       if cc==1
                           DMTpostCorr(:,:,ss) = CorrMat;
                     
                    elseif cc==2
                           PCBpostCorr(:,:,ss) = CorrMat;
                       
                        end
                    
                clear  CorrMat MatT  Data Avg_Signal
                toc
    end
end

%% Global Functional Connectivity at each parcel and total across the brain 

% determine experimental and control conditions matrices
exp = DMTpostCorr;
cont = PCBpostCorr;

% normalize data for stats
ExpFishZ =atanh(exp);
ContFishZ = atanh(cont);

% WE obtain GFC but taking the mean of fisher z values of each row (ie ROI)
for tt=1:size(ExpFishZ,3) %loop through subs
    for aal=1:size(ExpFishZ,1) %loop through area
        FCDDMTPost(aal,tt) = nanmean(ExpFishZ(aal,:,tt),2 );
        FCDPCBPost(aal,tt) = nanmean(ContFishZ(aal,:,tt),2 );
    end
end


% Determine mean of correlations. ie get FCD
for ss=1:size(ExpFishZ,3)
    for aal=1:size(ExpFishZ,1)
        fcdExp(aal,ss) = nanmean(ExpFishZ(aal,:,ss),2);
        fcdCont(aal,ss) = nanmean(ContFishZ(aal,:,ss),2);
        
%         fcdExp(aal,ss) = mean(ExpFishZ(aal,:,ss),2);
%         fcdCont(aal,ss) = mean(ContFishZ(aal,:,ss),2);
    end
end

% Now obtain a global connectivity per condition for each subject
 for tt=1:length(nets)
     net = nets{tt};
    GlobalnetDMTPost(tt,:) = squeeze(nanmean(FCDDMTPost(net,:),1));
    GlobalnetPCBPost(tt,:) = squeeze(nanmean(FCDPCBPost(net,:),1));
 end
 
save GlobalNet  GlobalnetDMTPost GlobalnetPCBPost

% stats on GFC at each ROI between conditions

for tt=1:size(xdata,1)
    
    [h,p(tt),ci,tst] = ttest(xdata(tt,:),ydata(tt,:));
    tstat(tt) = tst.tstat;
   
end

[fdrcorthres, fdrindp] = newFDR(p', 0.05) ; % FDR correction

if fdrcorrection==1
    sigaal = find(p<fdrcorthres) %
else
    sigaal = find(p<0.05) % sigaal provides the significant ROIs of GFC
end


% stats on total GFC across the brain between conditions

FCDtotalEXP = nanmean(xdata,1);
FCDtotalCont = nanmean(ydata,1);
[h,pTotal,ci,tstatTotal] = ttest(FCDtotalEXP,FCDtotalCont)

y = [FCDtotalCont'; FCDtotalEXP' ];

x = [ones(size(xdata,2),1) ; ones(size(ydata,2),1)+1];

figure
h = boxplot([FCDtotalCont', FCDtotalEXP'],'Labels',{'PCB','DMT'},'colors' ,'r')
% title('Functional Connectivity Density')
set(h,'LineWidth',4)

hold on
scatter(x,y,40, 'k', 'filled','MarkerFaceAlpha',0.5)
line([x(1:size(ydata,2)) x(size(ydata,2)+1:end)]',[y(1:size(ydata,2)) y(size(ydata,2)+1:end)]','LineWidth',2,'Color', [0 0 0 0.3])

set(gca, 'linewidth', 2,'Fontsize',26,'Box', 'on')
set(gcf, 'color', [1 1 1])
yl = ylabel('Fisher Z', 'FontWeight', 'bold','FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',30);

  text(1.25 ,0.6,['p = ' num2str(round(pTotal,3))],'FontSize',26,'HorizontalAlignment','left','Color', 'k');
  
 export_fig(sprintf('TotalGFC_GSR.png'),'-m2')
 
 %% Obtain GFC per network. 
 
 % This for the 7 networks
nets{1} = [1:9 51:58]; % Visual network
nets{2} = [10:15 59:66]; % SM
nets{3} = [16:23 67:73]; % DAN
nets{4} = [24:30 74:78]; % SAL
nets{5} = [31:33 79:80]; % Limbic
nets{6} = [34:37 81:89]; % FP
nets{7} = [38:50 90:100]; % DMN 
nets{8} = [101:112]; % SubC subcortical 'network' if that makes sense 


 for tt=1:length(nets)
     net = nets{tt};
       post(:,tt) =    (squeeze(nanmean(xdata(net,:),1)));
       pre(:,tt) =     (squeeze(nanmean(ydata(net,:),1)));
 end

 % now ttests per network
for ii=1:size(post,2)
    [h,pnets(ii),ci,tst] = ttest(post(:,ii),pre(:,ii));
    tstatnets(ii) = tst.tstat;
end

% fdr or bonferroni correction. Need fdr script
[adj pnetscor l] = fdr(pnets);
pbonfcor = pnets*length(pnets)

%% Plot GFC per parcell 

path2mask = '/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/CorrMatrix/Schaeff1007NExt.nii.gz'; % change here
      
Mask = load_nii(path2mask);

Maskfcd = Mask;

Maskfcd.img = zeros(size(Maskfcd.img));

for tt = 1:length(sigaal) % For each significant aal
    sal = sigaal(tt);
         idx = find(Mask.img == sal); % obtain 3d index for each significant aal.       
          Maskfcd.img(idx) = tstat(sal);
          tstat(sal)
          clear idx
end

save_nii(Maskfcd, 'GFCsig.nii.gz');

BrainNet_MapCfg('BrainMesh_ICBM152.nv','GFCsig.nii.gz')
export_fig(sprintf('GFC_DMTvsPCB_SanityCheck_1_6.png'),'-m2.5') % print figure
 
 %% Plot conditions separately 
 
plotsep=1;
Mask = load_nii(path2mask);

if plotsep==1
Maskdmt = Mask;
Maskpcb = Mask;

Maskdmt.img = zeros(size(Maskdmt.img));
Maskpcb.img = zeros(size(Maskpcb.img));


fcdmeanDMT = nanmean(xdata,2);%*10;
fcdmeanPCB = nanmean(ydata,2);%*10;


for tt = 1:length(fcdmeanDMT) % For each significant aal
         idx = find(Mask.img == tt); % obtain 3d index for each significant aal.       
          Maskdmt.img(idx) = fcdmeanDMT(tt);
          Maskpcb.img(idx) = fcdmeanPCB(tt);
          clear idx
end

save_nii(Maskdmt, 'FCDdmt1000.nii');
save_nii(Maskpcb, 'FCDpcb1000.nii');

BrainNet_MapCfg('BrainMesh_ICBM152.nv','FCDpcb1000.nii')
end

 export_fig(sprintf('FCDpcb100_0_0.4.png'),'-m2')


