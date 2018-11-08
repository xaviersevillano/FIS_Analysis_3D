%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Facial Improvement Score Analysis              %
%                                                               %
% Function that generates random groups (with different         %
% ammount of treated subjects) and perform the Facial           %
% Improvement Score (FIS) on them.                              %
%                                                               %
% As output a FIS histogram is calculated and saved as figure.  %
%                                                               %
% Input:                                                        %
%   res_path: Result path where the histigrams will be stored.  %
%   type_sample: 0 for humans , 1 for mice.                     %
%   age_range: [min max] vector defining the range of ages.     %
%   nameOutput: name for the stored histogram.                  %
%   treat_level: Treatment dosis level; 1 for low, 2 for high.  %
%                                                               %
% Notes:                                                        %
%   All humans have an predefined treatment dosis level of 1.   %
%   All mice have an predefined age of 1.                       %
%                                                               %
% Example:                                                      %
%   FIS_analysis_3D('results\',1,[0 2],'mice_low',1)            %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FIS_analysis_3D(res_path,type_sample,age_range,nameOutput,treat_level) 
%% Initialization
if ~exist(res_path,'dir'), mkdir(res_path); end
clc; close all;
scale       = 1;
idx         = [];
age         = [];
treatment   = [];
list_dir    = {};
sample_name = {};

if (type_sample==0)
    % Humans
    num_point   = 21;
    treat_level = 1;
    path_det        = 'Humans\Annotations3D\';
    path_dataset    = 'Humans\Dataset\';
    l_dir       = dir(fullfile(path_dataset,'*.ini'));
    % sample_name    =[sample_name; {l_dir.name}'];
    for i=1:length(l_dir)
        list_dir = [list_dir; [path_det l_dir(i).name(1:end-3) 'txt']]; 
        fid  = fopen([path_dataset l_dir(i).name]);
        a = textscan(fid,'%s');
        fclose(fid);
        sample_name = [sample_name {l_dir(i).name(1:end-4)}];
        idx         = [idx str2num(a{1}{13})];
        age         = [age str2num(a{1}{7})];
        treatment   = [treatment str2num(a{1}{16})];
    end 
else % Mice
    num_point   = 12;
    age_range = [0 2];
    path_det        = 'Mice\Annotations3D\';
    path_dataset    = 'Mice\Dataset\';
    l_dir       = dir(fullfile(path_dataset,'*.ini'));
    % sample_name    =[sample_name; {l_dir.name}'];
    for i=1:length(l_dir)
        list_dir = [list_dir; [path_det l_dir(i).name(1:end-3) 'txt']]; 
        fid  = fopen([path_dataset l_dir(i).name]);
        a = textscan(fid,'%s');
        fclose(fid);
        sample_name = [sample_name {l_dir(i).name(1:end-4)}];
        idx         = [idx str2num(a{1}{13})];
        age         = [age str2num(a{1}{7})];
        treatment   = [treatment str2num(a{1}{16})];
    end 
end

id_age      = and(age>=min(age_range),age<max(age_range));
sample_name = sample_name(id_age);
list_dir    = list_dir(id_age);
age         = age(id_age);
idx         = idx(id_age);
treatment   = treatment(id_age);

id_treat    = or(treatment==0,and(treatment== treat_level,idx>1));
sample_name = sample_name(id_treat);
list_dir    = list_dir(id_treat);
age         = age(id_treat);
idx         = idx(id_treat);
treatment   = treatment(id_treat);

idx(treatment>=1) = 4;

%% Read detections
x=[];
y=[];
z=[];
puntos = 1:num_point;
for i=1:length(list_dir)
    [X, Y, Z] = read_det(list_dir{i}, num_point, true);
    x = [x; X(puntos)./scale];
    y = [y; Y(puntos)./scale];
    z = [z; Z(puntos)./scale];
end
num_point = length(puntos);

%% Calculate centroid
C = [mean(x,2), mean(y,2), mean(z,2)];
CS = sum(((x-repmat(C(:,1),1,num_point)).^2+(y-repmat(C(:,2),1,num_point)).^2+(z-repmat(C(:,3),1,num_point)).^2),2).^0.5;

%% Translation and scale centering
x = (x-repmat(C(:,1),1,num_point))./repmat(CS,1,num_point);
y = (y-repmat(C(:,2),1,num_point))./repmat(CS,1,num_point);
z = (z-repmat(C(:,3),1,num_point))./repmat(CS,1,num_point);

%% Calculate Form Matrix
[FM, id_ij] = FM_creator_3D(x, y, z, num_point);

idx_EU = idx==1;
idx_SD = or(idx==2,idx==3);
% idx_SDP = idx==2;
idx_SD_ECGC = idx==4;

mean_dist_EU = mean(FM(:,idx_EU),2);
mean_dist_SD = mean(FM(:,idx_SD),2);
mean_dist_SD_ECGC = mean(FM(:,idx_SD_ECGC),2);

s_dist_EU = var(FM(:,idx_EU),0,2);
s_dist_SD = var(FM(:,idx_SD),0,2);
s_dist_SD_ECGC = var(FM(:,idx_SD_ECGC),0,2);

n_EU = sum(idx_EU);
n_SD = sum(idx_SD);
n_SD_ECGC = sum(idx_SD_ECGC);

diff_EUvsSD         = [];
diff_EUvsSD_ECGC    = [];
diff_EUvsSubGroup   = [];

%% Subgroups Analysis
n_Samples   = sum(idx_SD_ECGC);
if (type_sample==0) 
    n_SubGroups = 50*n_Samples;
else
    n_SubGroups = 10*n_Samples;
end
IDX_SubGroups = [];
CI_EU_SD = {};
CI_EU_ECGC = {};
for rep=1:n_SubGroups   
    fprintf('Generating %i of %i Subgroups\n',rep,n_SubGroups);
    % Subgroups Generation
    idx_SubGroup = Subgroups_Generation(floor((rep-1)/(n_SubGroups/n_Samples)),idx_SD,idx_SD_ECGC,n_Samples);
    IDX_SubGroups = [IDX_SubGroups;idx_SubGroup];
    
	% EDMA-II z test with montecarlo parametric bootstrapp
    alpha = 0.10;
    num_boot = 1000;
    CI_idx = floor([num_boot*(alpha/2), num_boot*(1-(alpha/2))]);

    FDM_EUvsSD      = mean_dist_EU - mean_dist_SD;
    FDM_EUvsSD_ECGC = mean_dist_EU - mean_dist_SD_ECGC;
    FDM_SDvsSD_ECGC = mean_dist_SD - mean_dist_SD_ECGC;

    mean_dist_SubGroup = mean(FM(:,idx_SubGroup),2);
    cov_EU=cov(FM(:,idx_EU)');
    cov_SD=cov(FM(:,idx_SD)');
    cov_SD_ECGC=cov(FM(:,idx_SD_ECGC)');
    cov_SubGroup = cov(FM(:,idx_SubGroup)');

    mean_dist_EU_boot=[];
    mean_dist_SD_boot=[];
    mean_dist_SD_ECGC_boot=[];
    mean_dist_SubGroup_boot=[];
    
    for i = 1:num_boot
        mean_dist_EU_boot       = [mean_dist_EU_boot        ; mean(mvnrnd(mean_dist_EU,cov_EU,size(FM(:,idx_EU),2)))];
        mean_dist_SD_boot       = [mean_dist_SD_boot        ; mean(mvnrnd(mean_dist_SD,cov_SD,size(FM(:,idx_SD),2)))];
        mean_dist_SD_ECGC_boot  = [mean_dist_SD_ECGC_boot   ; mean(mvnrnd(mean_dist_SD_ECGC,cov_SD_ECGC,size(FM(:,idx_SD_ECGC),2)))];
        mean_dist_SubGroup_boot  = [mean_dist_SubGroup_boot   ; mean(mvnrnd(mean_dist_SubGroup,cov_SubGroup,size(FM(:,idx_SubGroup),2)))];
    end
    
    FDM_EUvsSD_boot         = sort(mean_dist_EU_boot - mean_dist_SD_boot);
    FDM_EUvsSD_ECGC_boot    = sort(mean_dist_EU_boot - mean_dist_SD_ECGC_boot);
    FDM_EUvsSubGroup_boot   = sort(mean_dist_EU_boot - mean_dist_SubGroup_boot);


    diff_EUvsSD      = [diff_EUvsSD; or(FDM_EUvsSD_boot(CI_idx(1),:)>0,FDM_EUvsSD_boot(CI_idx(2),:)<0)];
    diff_EUvsSD_ECGC = [diff_EUvsSD_ECGC; or(FDM_EUvsSD_ECGC_boot(CI_idx(1),:)>0,FDM_EUvsSD_ECGC_boot(CI_idx(2),:)<0)];
    diff_EUvsSubGroup   = [diff_EUvsSubGroup; or(FDM_EUvsSubGroup_boot(CI_idx(1),:)>0,FDM_EUvsSubGroup_boot(CI_idx(2),:)<0)];
    
    CI_EU_SD{end+1} = [FDM_EUvsSD_boot(CI_idx(1),:)',FDM_EUvsSD_boot(CI_idx(2),:)'];
    CI_EU_ECGC{end+1} = [FDM_EUvsSD_ECGC_boot(CI_idx(1),:)',FDM_EUvsSD_ECGC_boot(CI_idx(2),:)'];

end

p_EUvsDS=sum(floor(sum(diff_EUvsSD)/n_SubGroups))*100/size(id_ij,2)
p_EUvsECGC=sum(floor(sum(diff_EUvsSD_ECGC)/n_SubGroups))*100/size(id_ij,2)
p_EUvsSubGroup=sum(diff_EUvsSubGroup,2)*100/size(id_ij,2);
FIS=(p_EUvsDS-p_EUvsSubGroup)*100/p_EUvsDS;
FIS_base=(p_EUvsDS-p_EUvsECGC)*100/p_EUvsDS
p_value = sum(FIS>FIS_base)/n_SubGroups

save([res_path nameOutput '.mat'],'diff_EUvsSD','diff_EUvsSD_ECGC','diff_EUvsSubGroup', 'IDX_SubGroups', 'FIS', 'FIS_base', 'CI_EU_SD', 'CI_EU_ECGC')

figure(1)
for i=0:n_Samples-1
subplot(n_Samples+1,1,i+1)
histogram(FIS(((n_SubGroups/n_Samples)*i)+1:(n_SubGroups/n_Samples)*(i+1)),floor(min(FIS)-1):floor(max(FIS)+1),'FaceColor',[0.8-(i*0.1) 0.8-(i*0.1) 0.8-(i*0.1)])
end
subplot(n_Samples+1,1,n_Samples+1)
histogram(FIS,floor(min(FIS)-1):floor(max(FIS)+1),'FaceColor',[0 0 0])
saveas(gcf,[res_path nameOutput],'fig')
saveas(gcf,[res_path nameOutput],'png')

end

function [X, Y, Z] = read_det(name, num_points, DIM)
    fdet=fopen(name,'r');
    data=textscan(fdet,'%f');
    fclose(fdet);
    X   = [];
    Y   = [];
    Z   = [];
    if (~DIM), data    = reshape(data{1},2,size(data{1},1)/2);
    else data    = reshape(data{1},3,size(data{1},1)/3); end
    X       = reshape(data(1,:),num_points,size(data,2)/num_points)';
    Y       = reshape(data(2,:),num_points,size(data,2)/num_points)';
    if (DIM), Z=reshape(data(3,:),num_points,size(data,2)/num_points)'; end
end

function [FM, id_ij] = FM_creator_3D(x, y, z, num_point)
    FM = zeros(sum(1:num_point-1),size(x,1));
    id_ij = [];
    for s=1:size(x,1)
        k = 1;
        for i=1:size(x,2)
            S_i = [x(s,i), y(s,i), z(s,i)];
            for j = i+1:size(x,2)
                S_j     = [x(s,j), y(s,j), z(s,j)];
                FM(k,s) = (sum((S_i-S_j).^2)).^0.5;
                if (s==1); id_ij   = [id_ij, [i;j]]; end;
                k       = k+1;
            end
        end
        FM(:,s) = FM(:,s)./geomean(FM(:,s));
    end
end

function idx_sub = Subgroups_Generation(n_treated,idx_SD,idx_SD_ECGC,n)
ECGC = find(idx_SD_ECGC);
SD = find(idx_SD);

new_ECGC = datasample(ECGC,n_treated,'Replace',false);
if (n-n_treated<length(SD))
    new_SD = datasample(SD,n-n_treated,'Replace',false);
else
    new_SD = datasample(SD,n-n_treated,'Replace',true);
end
idx_sub = false(size(idx_SD));
idx_sub([new_ECGC,new_SD]) = true;
end