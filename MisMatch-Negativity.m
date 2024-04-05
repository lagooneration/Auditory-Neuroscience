%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MISMATCH NEGATIVITY using FIELDTRIP 


close all; clear; clc
spm('defaults', 'eeg');

% set the block
block = 'fdMMN'; % 'fdMMN' 'oMMN'

root            = ['D:\work\MMN\wo_polhemus' '_' block filesep];
path_out         = [ root '\Figures\']; 

prefix          = 'lp_subject';
suffix          = [];
file_extn       = 'png';

% set the list of conditions
if strcmp(block, 'fdMMN')
    
    cond_list = {'f1 STD4' 'f1 STD-L' 'f1 DEV' 'f2 STD4' 'f2 STD-L' 'f2 DEV' ...
                 'd STD4'  'd STD-L'  'd DEV' }; 
    
elseif strcmp(block, 'oMMN')
    
    cond_list = {'o STD4'  'o STD-L'  'o DEV'};
    
end

subject_list    = [3 5:7 9:15 17:22] ;

conds_to_compare = [9 8];

cond_list = cond_list(conds_to_compare);

chans       = '[A-D]';
save_fig    = 1;

alpha       = 0.05;
stat_test_tail = -1;

flag_behav  = 0;
layout_dir  = 'D:\tools\fieldtrip-20170111\template\layout\';

% times       = [0 1.5];


%% Process

% Load behavioural performance
if flag_behav == 1
    load('results_behaviour')
end

k = 1;
avgs = [];
bad_channels = [];

for subj = subject_list
    
    disp(['Processing Subject: ' num2str(subj)]);
    
    path = [root prefix num2str(subj) suffix '.mat'];
    
    D = spm_eeg_load(path);
    
    for jj = 1:numel(cond_list)
        % disp(['Processing Subject: ' num2str(subj_list(s)) ' cond: ' cond_list{jj}]);
        
        if flag_behav == 1
            % select all trials of given condition with correct response
            correct_response = find(strcmp(behaviour{subj},[ cond_list{jj} ' C']));
            
            % remove bad trials
            selected_trials = setdiff(correct_response, D.badtrials);
        else
            % retain good trials of given condition
            selected_trials = D.indtrial(cond_list{jj},'GOOD');
        end
        
        % convert data to field trip format
        % data = ftraw(D,subset,'',selected_trials);
        tmp_data = fttimelock(D,'EEG','',selected_trials,'');
        
        % compute mean per condition
        cfg = [];
        timelock = ft_timelockanalysis(cfg, tmp_data);
        
        data{jj} = timelock;
        
    end
    
    cond_A{k} = data{1};
    
    cond_B{k} = data{2};
    
    k = k + 1;
end

%% perform the statistical test using randomization and a clustering approach
% using the NEW freqstatistics function

% returns 'neighbours' struct required for computation of cluster level metric
load('X:\auditory\EEG\Resources\biosemi128_neighb.mat');

% select data from all EEG channels from temporal lobe
chan_list = selectchannels(D,['regexp_^' chans '.*']);

cfg                  = [];
cfg.neighbours       = neighbours;
cfg.minnbchan        = 2;
cfg.neighbourdist    = 4;
% set the range of latencies to compute the output
cfg.latency          = [0 inf];
% selection of channels - Nx1 cell-array or 'all' (takes a lot of time though)
cfg.channel          = tmp_data.label(chan_list); % see CHANNELSELECTION
% set whether to average the response across time
cfg.avgovertime      = 'no'; % 'yes' or 'no'
% set whether to average the response across channels
cfg.avgoverchan      = 'no'; % 'yes' or 'no'
% set the measure to use for quantifying the effect at the single-sample level
cfg.statistic        = 'ft_statfun_depsamplesT';
% 'maxsum' (default), 'maxsize', or 'wcm'
cfg.clusterstatistic = 'maxsum';
% Method for calculation of the significance probability
cfg.method           = 'montecarlo';
% number of random draws in monte carlo method
cfg.numrandomization = 1000;
% Correction for multiple comparisons - 'no', 'max', 'cluster', 'bonferoni', 'holms', or 'fdr'
cfg.correctm         = 'cluster';
% choose between a one-sided (1, -1) and a two-sided (0) statistical test
cfg.tail             = stat_test_tail;
cfg.clustertail      = stat_test_tail;
% false alarm rate of the permutation test
cfg.alpha            = alpha/(2-abs(stat_test_tail));
% cfg.clusteralpha     = 0.05;
% store information about the Unit of Observations
cfg.design           = [1:numel(subject_list) 1:numel(subject_list); % subject number
                        ones(1,numel(subject_list)) 2*ones(1,numel(subject_list))];  % condition number
% indicate the row of the design matrix that contains the unit variable
cfg.uvar             = 1; % "subject" is unit of observation
% indicate the row of the design matrix that contains the independent variable
cfg.ivar             = 2; % "condition" is the independent variable

stat    = ft_timelockstatistics(cfg, cond_A{:}, cond_B{:});

%%
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
GA_cond_A     = ft_timelockgrandaverage(cfg, cond_A{:});  
GA_cond_B     = ft_timelockgrandaverage(cfg, cond_B{:});

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_diff_cond_A_vs_B = ft_math(cfg, GA_cond_A, GA_cond_B);

% get relevant (significant) values
if stat_test_tail >= 0
    pos_cluster_pvals = [stat.posclusters(:).prob];
    pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
    pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
end

if stat_test_tail <= 0
    neg_cluster_pvals = [stat.negclusters(:).prob];
    neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
    neg = ismember(stat.negclusterslabelmat, neg_signif_clust);
end

%%

% plot

% subplot parameters in combo plot
% p = fnumSubplots(numel(start_times));

disp ' ';
mkdir(path_out);

start_times = [0:0.05:0.2]; 
time_duration = [0.05*ones(1,numel(start_times))];

screen_size = get(0, 'ScreenSize');

for k = 1:numel(start_times);
    
    time_pos = start_times(k);
    time_dur = time_duration(k);
    
    cfg = [];
    cfg.xlim = (time_pos +[0 time_dur]);
    % cfg.zlim = [-5e-14 5e-14];
    
    sample_span = single(ceil(1+cfg.xlim(1)*D.fsample):ceil(cfg.xlim(2)*D.fsample));
    
    if stat_test_tail >= 0
        pos_int = all(pos(:, sample_span), 2);
    end
    if stat_test_tail <= 0
        neg_int = all(neg(:, sample_span), 2);
    end
    if stat_test_tail == 0
        highlight_channel = find(or(pos_int,neg_int));
    elseif stat_test_tail == 1
        highlight_channel = find(pos_int);
    elseif stat_test_tail == -1
        highlight_channel = find(neg_int);
    end
    
    cfg.highlight   = 'labels';
    cfg.highlightsymbol  = 'x';
    cfg.highlightchannel = highlight_channel;
    cfg.highlightcolor   = [1 1 1];
    cfg.comment     = 'xlim';
    cfg.commentpos  = 'title';
    cfg.layout      = [layout_dir filesep 'biosemi128.lay'];
    cfg.colorbar   	= 'yes';
    cfg.colormap    = 'jet';

    h_fig = figure;
    % subplot(p(1),p(2),k);
    
    ft_topoplotER(cfg, GA_diff_cond_A_vs_B);
    
    if save_fig
        h_fig.InvertHardcopy = 'off';
%         set(h_fig, 'Position', [0 0 screen_size(3) screen_size(4) ] ); % set to scren size
%         set(h_fig, 'PaperPositionMode','auto');
        saveas(h_fig, [path_out 'EEG-Scalptopo-select_subj-' strrep(cell2mat(cond_list),'.','') '-diff_cond' '-time_' strrep(num2str(time_pos),'.','p') 's' '-dur_' strrep(num2str(time_dur),'.','p') 's'], file_extn);
    end
    close(h_fig);
end  
