
%% this scripts relies on functions from the BrainSpace toolbox
% https://brainspace.readthedocs.io/en/latest/

%% load surface-projected DMT data

DMT_dir='/Volumes/MANESH/Psychedelic_data/DMT/FSA5_MGH';
outdir='/Users/manesh/Desktop/psychedelic_gradients/DMT_Chris';

L1_DMTpre='*DMTpre*lh*.mgh';
R1_DMTpre='*DMTpre*rh*.mgh';
DMTpre_data_L1 = deblank(icatb_listFiles_inDir(DMT_dir, L1_DMTpre));
DMTpre_data_R1 = deblank(icatb_listFiles_inDir(DMT_dir, R1_DMTpre));

L1_DMTpost='*DMTpost*lh*.mgh';
R1_DMTpost='*DMTpost*rh*.mgh';
DMTpost_data_L1 = deblank(icatb_listFiles_inDir(DMT_dir, L1_DMTpost));
DMTpost_data_R1 = deblank(icatb_listFiles_inDir(DMT_dir, R1_DMTpost));

L1_PCB='*PLCB*lh*.mgh';
R1_PCB='*PLCB*rh*.mgh';
PLCB_data_L1 = deblank(icatb_listFiles_inDir(DMT_dir, L1_PCB));
PLCB_data_R1 = deblank(icatb_listFiles_inDir(DMT_dir, R1_PCB));

%% compute gradients
for k=1:16

    currFile_pre_L1=deblank(DMTpre_data_L1(k,:));
    currFile_pre_R1=deblank(DMTpre_data_R1(k,:));
    
    currFile_pre_L1=deblank(DMTpost_data_L1(k,:));
    currFile_pre_R1=deblank(DMTpost_data_R1(k,:));

    currFile_PLCB_L1=deblank(PLCB_data_L1(k,:));
    currFile_PLCB_R1=deblank(PLCB_data_R1(k,:));

    curr_Z_DMTpre(k,:,:) = getZmat(DMT_dir,currFile_pre_L1, currFile_pre_R1); %combine hemis and get downsampled FC matrix
    curr_Z_DMTpost(k,:,:)= getZmat(DMT_dir,currFile_post_L1, currFile_post_R1);
    curr_Z_PLCB(k,:,:) = getZmat(DMT_dir,currFile_PLCB_L1, currFile_PLCB_R1);

    gm = GradientMaps();
    curr_gm = gm.fit(curr_Z_DMTpre(k,:,:));
    gradients_DMTpre(k,:).gradients = curr_gm.gradients;
    gradients_DMTpre(k,:).lambda = curr_gm.lambda;
    
    gm = GradientMaps();
    curr_gm = gm.fit(curr_Z_DMTpost(k,:,:));
    gradients_DMTpost(k,:).gradients = curr_gm.gradients;
    gradients_DMTpost(k,:).lambda = curr_gm.lambda;
    
    gm = GradientMaps();
    curr_gm = gm.fit(curr_Z_PLCB(k,:,:));
    gradients_PLCB(k,:).gradients = curr_gm.gradients;
    gradients_PLCB(k,:).lambda = curr_gm.lambda;   
end

mean_Z = mean(cat(1,curr_Z_DMTpre,curr_Z_DMTpostcurr_Z_PLCB),1);

gm = GradientMaps();
curr_gm = gm.fit(mean_Z);
mean_template.gradients = curr_gm.gradients;
mean_template.lambda = curr_gm.lambda;

save(strcat(outdir,'/gradients_DMTpre.mat'),'gradients_DMTpre');
save(strcat(outdir,'/gradients_DMTpost.mat'),'gradients_DMTpost');
save(strcat(outdir,'/gradients_DMTplcb.mat'),'gradients_PLCB');
save(strcat(outdir,'/mean_template.mat'),'mean_template');

%% align gradients (procrustes rotation)

all_embeddings=[gradients_DMTpre gradients_DMTpost gradients_PLCB];

[DMTgradients16_PREPOSTPLCB_aligned.gradients, DMTgradients16_PREPOSTPLCB.xfms] = mica_iterativeAlignment(all_embeddings,5,mean_template);

save('DMTgradients16_PREPOST_aligned.mat','DMTgradients16_PREPOST_aligned')

%% display gradients

load('/Volumes/MANESH/Psychedelic_data/fsa5_to_fsa10k/fsa_10k_midsurface_LR.mat')
SurfStatViewData(DMTgradients16_PREPOSTPLCB_aligned.gradients{1,1}(:,1),fsa_10k)
colormap(parula)



