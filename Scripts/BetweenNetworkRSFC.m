%% BETWEEN RSFC FIG 2B
% Chris Timmermann 2022
% here we do regression to determine between network connectivity

list_ICA=1:7; 
name = {'VIS', 'SM', 'LIM', 'DAN', 'SAL', 'FP', 'DMN'};

numICA=numel(list_ICA);

sbj = 1:34; Number of subjects

for sub=sbj
	sub
    
  path2analysis = ['/Users/christophertimmermann/Documents/Imaging_fMRI/Raw/analysis/Yeo7NStriatbinaryAftvsAft/dr_stage1_subject',sprintf('%05d',sub-1),'.txt']; % select here the output of RSFC_Network.sh file. Each one containing the mean PE values
    
    tcs=load(path2analysis);
    tcs = tcs(:, [1 2 5 3 4 6 7]);
    all(:,:,sub) = tcs;

    for i=1:numICA
        i
        i_num=list_ICA(i);
        for j=1:numICA
            j_num=list_ICA(j);
            if j~=i

                Y=tcs(:,i_num);	
                X=[ones(size(tcs(:,j_num))) tcs(:,j_num)];             

               [TempBeta,CI]=regress(Y,X);
                beta_maineffect(i,j,sub)=TempBeta(2);

            end
        end
    end

    
    for i=1:numICA
        for j=1:numICA
            sym_beta_maineffect(i,j,sub)=(beta_maineffect(i,j,sub)+beta_maineffect(j,i,sub))/2;
            sym_beta_maineffect(j,i,sub)=sym_beta_maineffect(i,j,sub);
        end
    end

end


colors = cbrewer('div', 'RdBu', 2000);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
 
totsub = size(sym_beta_maineffect,3);

%difference for each subject

for sub=1:totsub/2
    for i=1:numICA
        for j=1:numICA
            diff_beta(j,i,sub)=sym_beta_maineffect(i,j,(sub+(totsub/2)))-sym_beta_maineffect(i,j,sub);
        end
    end
end
           
range = [-0.7 0.7];
load('mycoolwarm.mat')

figure(1) % experimental condition betas
after_means=mean(sym_beta_maineffect(:,:,1:size(sym_beta_maineffect,3)/2),3)

% get rid of diagonal
for jj=1:size(after_means,2)
    after_means(jj,jj) = nan;
end

% plot experimental conditions
imagesc(after_means,range)
% colormap(mycoolwarm)
colormap(colors);
set(gca,'XTick',[1:12])
set(gca,'YTick',[1:12])
set(gca, 'XTickLabel', name); % set x-axis labels
set(gca, 'YTickLabel', name); % set y-axis labels
allfigs = allchild(gcf);
set(allfigs(10), 'fontsize',20, 'Box', 'on');
caxis([-1 1])
set(gca,'XTickLabelRotation',45)
tbh = colorbar
set(tbh,'YTick',[-1:1:1],'fontsize',15,  'FontWeight', 'bold')
ylabel(tbh, 'r','Position', [0 0 0])
set(gcf, 'color', [1 1 1],'position', [552 786 325 271]);

title('DMT', 'fontsize',20)

set(gca, 'ydir', 'normal');
axis square;
set(gcf, 'Position', [800 783 497 274]);
export_fig(sprintf('DMTbetween.png'),'-m2.5')


figure(2) % control condition betas
before_means=mean(sym_beta_maineffect(:,:,size(sym_beta_maineffect,3)/2+1:size(sym_beta_maineffect,3)),3);
for jj=1:size(before_means,2)
    before_means(jj,jj) = nan;
end

% plot
imagesc(before_means,range)
colormap(colors);
set(gca,'XTick',[1:12])
set(gca,'YTick',[1:12])
set(gca, 'XTickLabel', name); % set x-axis labels
set(gca, 'YTickLabel', name); % set y-axis labels
allfigs = allchild(gcf);
set(allfigs(10), 'fontsize',20, 'Box', 'on');
caxis([-1 1])
set(gca,'XTickLabelRotation',45)
gbh = colorbar
set(gbh,'YTick',[-1:1:1],'fontsize',15,  'FontWeight', 'bold')
ylabel(gbh, 'r','Position', [0 0 0])
set(gcf, 'color', [1 1 1],'position', [552 786 325 271]);


title('PCB', 'fontsize',20)

set(gca, 'ydir', 'normal');
axis square;
set(gcf, 'Position', [800 783 497 274]);
 export_fig(sprintf('PCBbetween.png'),'-m2.5')

 % Plot the contrast
figure(3) % after minus before connectivity 
[H,P,CI,STATS_interaction_diff]=ttest(sym_beta_maineffect(:,:,1:totsub/2),sym_beta_maineffect(:,:,(totsub/2)+1:totsub),[],[],3)
for jj=1:size(STATS_interaction_diff.tstat,2)
    STATS_interaction_diff.tstat(jj,jj) = nan;
end

maxtstat = abs(max(STATS_interaction_diff.tstat(:)));

imagesc(STATS_interaction_diff.tstat,[maxtstat*-1 maxtstat])
colormap(colors);
set(gca,'XTick',[1:7])
set(gca,'YTick',[1:7])
set(gca, 'XTickLabel', name); % set x-axis labels
set(gca, 'YTickLabel', name); % set y-axis labels
allfigs = allchild(gcf);
set(allfigs(10), 'fontsize',20, 'Box', 'on');
% caxis([-1 1])
set(gcf, 'color', [1 1 1]);    
set(gca,'XTickLabelRotation',45)
cbh = colorbar
set(cbh,'YTick',[floor(maxtstat*-1)+1 0 ceil(maxtstat)-1],'fontsize',15,  'FontWeight', 'bold')
ylabel(cbh, 'T','Position', [0 0 0])
set(gcf, 'color', [1 1 1],'position', [552 786 325 271]);

hold on
% Fill with nans half of the table to do fdr correction
P2= P;

for pp=1:size(P,2)
    P(pp,pp:end) = nan;
end

% now do fdr correction

temp = reshape(P.',1,size(P,2)*size(P,2)).';
idxnonann = find(~isnan(temp));
temp(any(isnan(temp),2),:) = [];
[thr cor adj] = fdr(temp');

 var=cor;

% fill the table with corrected values
pcor=nan(size(P)) ;

pcor(idxnonann) = var;

[r c] = find(pcor<0.05);

plot(r,c,'*','color', 'w','MarkerSize',10)
plot(c,r,'*','color', 'w','MarkerSize',10)

title('DMT vs PCB', 'fontsize',20)

set(gca, 'ydir', 'normal');
axis square;
 export_fig(sprintf('DMTvsPCBbetween.png'),'-m2.5')

