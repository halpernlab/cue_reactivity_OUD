function [zmapthresh] = cluster_permutation_main_1d(cond1,cond2)
% This is only for 1D (LFP) cluster permutation. Not 2D
% Data should be fixed timepoint.
% input: cond1 & cond2 : trial x timepoints
% output: p<0.05 timepoint
%% perform perm testing
% target conds: both_conds, baseline(cond2) + cond1:control 
% Cond1 - Cond2
both_conds = cat(1,cond2(:,:), cond1(:,:));

real_condition_mapping = [-ones(1,size(cond1,1)) ones(1,size(cond2,1))]; % # of each chns % [-1 -1 1 1]
nTimepoints            = size(both_conds,2); % # of time
voxel_pval             = 0.05;
mcc_voxel_pval         = 0.05;
mcc_cluster_pval       = 0.05;
n_permutes             = 1000; % 1000

% compute actual t-test of difference (using unequal N and std)
% tnum is the same with diff_map. (nanmean of the different channels)
% tdenom is nanstd of each power spectrum data. to get a real t-value
tnum   = squeeze(nanmean(both_conds(real_condition_mapping==1,:),1) - nanmean(both_conds(real_condition_mapping==-1,:),1)); % result: freq x time. difference between two condition in mean
tdenom = sqrt((nanstd(both_conds(real_condition_mapping==1,:),0,1).^2)./sum(real_condition_mapping==1)...
    + (nanstd(both_conds(real_condition_mapping==-1,:),0,1).^2)./sum(real_condition_mapping==-1));
real_t = tnum./tdenom;

% initialize null hypothesis matrices
permuted_tvals  = zeros(n_permutes,nTimepoints);
max_clust_info  = zeros(n_permutes,1);

% for example in sub1, shuffle 8 chns and reassign new label (cond1 cond2)
% 8 power-spectrum images will be shuffled
% generate pixel-specific null hypothesis parameter distributions
for permi = 1:n_permutes
    % choose how many subj to perm
    fake_condition_mapping = real_condition_mapping;
    shuf_cntr = [];

%     for ii = 1:size(cond1,1)
%         shuf_cntr = [shuf_cntr randi([0 1],1)];
%     end
%     shuffle = logical([shuf_cntr  shuf_cntr]);
    for ii = 1:length(fake_condition_mapping)
        shuf_cntr = [shuf_cntr randi([0 1],1)];
    end
    shuffle = logical(shuf_cntr);
    fake_condition_mapping(shuffle) = -1*real_condition_mapping(shuffle);
    
    % compute t-map of null hypothesis
    tnum   = squeeze(nanmean(both_conds(fake_condition_mapping==-1,:),1)-nanmean(both_conds(fake_condition_mapping==1,:),1));
    tdenom = sqrt( ((nanstd(both_conds(fake_condition_mapping==-1,:),0,1).^2)./sum(fake_condition_mapping==-1)) ...
                  + ((nanstd(both_conds(fake_condition_mapping==1,:),0,1).^2)./sum(fake_condition_mapping==1)) );
    tmap   = tnum./tdenom;

    % save all permuted values
    permuted_tvals(permi,:,:) = tmap;
    
    % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
    % here, clusters obtained by parametrically thresholding the t-maps
    tmap(abs(tmap)<tinv(1-voxel_pval,size(both_conds,1)-2))=0;

    % get number of elements in largest supra-threshold cluster
    clustinfo = bwconncomp(tmap);
    max_clust_info(permi) = max([0 cellfun(@numel,clustinfo.PixelIdxList)]); % zero accounts for empty maps

%     disp(permi)
end

% compute Z-map
zmap = (real_t-squeeze(nanmean(permuted_tvals,1)))./squeeze(nanstd(permuted_tvals));
zmapthresh = zmap;
% going from Z value to p value % norminv(1-voxel_pval) is the standard def unit that corrsep to p<.05
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=false; 
zmapthresh=logical(zmapthresh);

%% correct for mc
% apply cluster-level corrected threshold
zmapthresh = zmap;
% uncorrected pixel-level threshold
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;
clustinfo  = bwconncomp(zmapthresh);
clust_info = cellfun(@numel,clustinfo.PixelIdxList);
% % MC based on size
clust_size_threshold   = prctile(max_clust_info,100-mcc_cluster_pval*100); % among the cluster sizes, find 95% of the size
% identify clusters to remove
whichclusters2remove   = find(clust_info<=clust_size_threshold); % remove clusters of small size than size threshold

% remove clusters
for i=1:length(whichclusters2remove)
    zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
end




end

