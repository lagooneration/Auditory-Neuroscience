%Calculate ERP response by multiple regression
%   Usage:
%   	rerp_result = rerp(data, rerp_profile);
%           Run analysis defined in rerp_profile against data
%
%   Parameters:
%       data:
%           CONTINUOUS data from EEG.data or EEG.icaact
%
%       rerp_profile:
%           RerpProfile object
%
%   See also:
%       RerpProfile, RerpResult, pop_rerp, pop_rerp_study

function rerp_result = rerp(data, rerp_profile, varargin)
import rerp_dependencies.*
p = rerp_profile;
s = p.settings;

assert(length(s.penalty_func) <= 1, 'rerp: can only calculate one penalty function at a time');

rerp_result=RerpResult(p);
tic;

% Want a tall matrix
if size(data, 2) == rerp_profile.pnts
    data = permute(data,[2 1]);
end

assert(size(data,1) == rerp_profile.pnts, 'rerp: data must have time as first dimension, channels or components as second dimension');

dim = size(data,2);
datalength = size(data,1);
nbins = size(data,3);

if nargin < 2
    help rerp;
    return;
end;

if nbins > 1
    dim = dim*nbins;
    data = reshape(data, [datalength, dim]);
    new_time_series = zeros(1, size(data,2));
    
    if s.type_proc
        for m=1:dim
            new_time_series(m) = p.include_chans(ceil(m/nbins));
        end
        p.include_chans=new_time_series;
        p.include_comps=[];
        
    else
        for m=1:dim
            new_time_series(m) = p.include_comps(ceil(m/nbins));
        end
        p.include_comps=new_time_series;
        p.include_chans=[];
    end
    
    ersp_flag=1;
    
else
    ersp_flag=0;
end

% Assert valid artifact index configuration
assert(~(s.artifact_rejection_enable && isempty(p.computed_artifact_indexes) && isempty(p.variable_artifact_indexes)), 'rerp: either disable artifact rejection or provide artifact index variable in rerp_profile.computed_artifact_indexes OR rerp_profile.variable_artifact_indexes (see pop_rerp doc for details');

% Set up artifact indexes based on profile
if s.artifact_rejection_enable
    artifact_indexes=p.computed_artifact_indexes;
    
    if s.artifact_variable_enable
        if ~isempty(p.variable_artifact_indexes)
            artifact_indexes=p.variable_artifact_indexes;
        else
            disp('rerp: artifact variable not found, using computed indexes');
        end
    end
    
else
    artifact_indexes=false(datalength,1);
end

s.first_phase_lambda=s.first_phase_lambda(:);
if numel(s.first_phase_lambda) < 3
    disp('rerp: rerp_profile.settings.first_phase_lambda has less than 3 elements, resetting to default');
    s.first_phase_lambda=[-log10(1e8); s.first_phase_lambda; log10(1e8)];
end

if s.num_xvalidation_folds < 2
    disp('rerp: rerp_profile.settings.num_xvalidation_folds < 2, using 2 instead');
    s.num_xvalidation_folds=2;
end

%Make sure artifact indexes are valid
assert(length(artifact_indexes)==datalength, 'rerp: artifact index vector must be same length as data');
if size(artifact_indexes, 2) == rerp_profile.pnts
    artifact_indexes = permute(artifact_indexes,[2 1]);
end

%We don't add a row of ones to predictor: assumes no dc bias in data.
fprintf('rerp: generating predictor, time=%f\n', toc);
[predictor, data_pad, parameter_idx_layout, tags] = p.predictor;

%Extend data with zeros to match predictor size/alignment
data = [zeros(data_pad(1), size(data,2)); data; zeros(data_pad(2), size(data,2))];
artifact_indexes = [false(data_pad(1),1); artifact_indexes; false(data_pad(2),1)];

%Discard artifact frames
data = double(data(~artifact_indexes, :));
predictor = predictor(~artifact_indexes, :);

%Don't worry about these
arg_parser=inputParser;
addOptional(arg_parser,'P', []);
addOptional(arg_parser,'q',[]);
addOptional(arg_parser,'L', []);
addOptional(arg_parser,'xval_train_idx', {});
addOptional(arg_parser,'xval_test_idx', {});
addOptional(arg_parser,'xval_P', {});
addOptional(arg_parser,'xval_q', {});
addOptional(arg_parser,'xval_L', {});
parse(arg_parser, varargin{:});

%Cache factors
if isempty(arg_parser.Results.P)
    P = predictor'*predictor;
else
    P = arg_parser.Results.P;
end

if isempty(arg_parser.Results.q)
    q = predictor'*data;
else
    q = arg_parser.Results.q;
end

L={};
if isempty(arg_parser.Results.L)
    for m=1:length(s.penalty_func)
        str = s.penalty_func{m};
        if strcmp(str,'Elastic net')||strcmp(str,'L1 norm')
            fprintf('rerp: caching LU factorization for ADMM (L1 norm or Elastic net), time=%f\n', toc);
            RHO=1.5;
            L = sparse(chol(P + RHO*speye(size(P)), 'lower' ));
            break;
        end
    end
else
    L = arg_parser.Results.L;
end

% When doing cross validation, we need predictors based on subsets of the
% data
xval_train_idx=arg_parser.Results.xval_train_idx;
xval_test_idx=arg_parser.Results.xval_test_idx;
xval_P=arg_parser.Results.xval_P;
xval_q=arg_parser.Results.xval_q;
xval_L=arg_parser.Results.xval_L;

if isempty(xval_train_idx)||isempty(xval_test_idx)||isempty(xval_P)||isempty(xval_q)
    fprintf('rerp: setting up cross validation predictors, time=%f\n', toc');
    setup_xval_predictors;
end

if s.regularization_enable
    if ~isempty(s.penalty_func)
        
        for m=1:length(s.penalty_func)
            this_result = RerpResult(p);
            
            %Finding lambda with cross validation
            if s.cross_validate_enable
                fprintf('rerp: beginning %s grid search, time=%f\n',this_result.analysis_name, toc);
                [~, enidx] = intersect('Elastic net', s.penalty_func);
                
                if strcmp(s.penalty_func{m},'L1 norm')
                    if ~(any(enidx) && s.elasticnet_quick_zoom)
                        this_result.gridsearch.lambda_range= repmat({s.first_phase_lambda}, 1, dim);
                        this_result.analysis_name = 'L1 norm';
                    end
                end
                
                if strcmp(s.penalty_func{m},'L2 norm')
                    
                    if ~(any(enidx) && s.elasticnet_quick_zoom)
                        this_result.gridsearch.lambda_range= repmat({s.first_phase_lambda}, 1, dim);
                        this_result.analysis_name = 'L2 norm';
                    end
                end
                
                if strcmp(s.penalty_func{m},'Elastic net')
                    this_result.gridsearch.lambda_range= repmat({[s.first_phase_lambda, s.first_phase_lambda]}, 1, dim);
                    this_result.analysis_name= 'Elastic net';
                end
                
                %Keep either the full result, or a lightweight copy
                rerp_result = rerp_grid_search(this_result, 1);
                if ~s.save_grid_search
                    rerp_result.gridsearch=[];
                end
                
            else
                
                %Specified Lambda
                fprintf('rerp: beginning %s single lambda, time=%f\n' ,this_result.analysis_name, toc);
                
                if strcmp(s.penalty_func{m},'L1 norm')
                    this_result.lambda=repmat({s.lambda(1)}, 1, dim);
                    this_result.analysis_name = 'L1 norm';
                end
                
                if strcmp(s.penalty_func{m},'L2 norm')
                    this_result.analysis_name = 'L2 norm';
                    this_result.lambda=repmat({s.lambda(2)}, 1, dim);
                end
                
                if strcmp(s.penalty_func{m},'Elastic net')
                    this_result.analysis_name = 'Elastic net';
                    this_result.lambda=repmat({s.lambda}, 1, dim);
                end
                
                rerp_result = rerp_regularized(this_result);
            end
        end
    else
        error('rerp: regularization is enabled, but a penalty function was not specified');
    end
    
else
    this_result = RerpResult(p);
    this_result.analysis_name='Least squares';
    this_result.lambda=cell( 1, size(data, 2));
    rerp_result = rerp_regularized(this_result);
end

%Get rsquare performance on a per event type basis
cross_validate_event_types(rerp_result);

rerp_result.compute_time_seconds=toc;
rerp_result.date_completed=datestr(now,'yyyy-mm-dd-HH:MM:SS');
rerp_result.ersp_flag=ersp_flag;

dsname = regexp(rerp_result.rerp_profile.eeglab_dataset_name,'.*[\\\/](.*)\.set', 'tokens');
if ersp_flag
    type=' rERSP';
else
    type='';
end

rerp_result.name=[dsname{1}{1} ' ' rerp_result.analysis_name type ' ' rerp_result.date_completed];
rerp_result.parameter_idx_layout=parameter_idx_layout; 
rerp_result.tags=tags; 
disp('rerp: done');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regression functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Recursively zoom in on lambda
    function [rerp_result, optimal_result] = rerp_grid_search(rerp_result, level)
        rerp_result.gridsearch.level = level;
        
        if strcmp(rerp_result.analysis_name, 'Elastic net')
            % Quickly find a ballpark range of lambda by searching L1 and
            % L2 norm axes seperately.
            if s.elasticnet_quick_zoom;
                if level==1
                    % Find best L1 and L2 lambda seperately first: get a closer range of lambda to start
                    % the elastic net grid search
                    fprintf('rerp: elastic net quick zoom, starting L1 norm search, time=%f\n',toc);
                    L1_rerp_result = RerpResult({});
                    L1_rerp_result.analysis_name={'L1 norm'};
                    L1_rerp_result.gridsearch.lambda_range = repmat({rerp_result.gridsearch.lambda_range{1}(:,1)}, 1, dim);
                    [L1_rerp_result, L1_optimal_params] = rerp_grid_search(L1_rerp_result, level);
                    rerp_result.gridsearch.L1_norm=L1_rerp_result;
                    
                    fprintf('rerp: elastic net quick zoom, starting L2 norm search, time=%f\n',toc);
                    L2_rerp_result = RerpResult({});
                    L2_rerp_result.analysis_name={'L2 norm'};
                    L2_rerp_result.gridsearch.lambda_range = repmat({rerp_result.gridsearch.lambda_range{1}(:,2)}, 1, dim);
                    [L2_rerp_result, L2_optimal_params] = rerp_grid_search(L2_rerp_result, level);
                    rerp_result.gridsearch.L2_norm=L2_rerp_result;
                    
                    lambda1 = L1_optimal_params.lambda;
                    lambda2 = L2_optimal_params.lambda;
                    
                    opti_lambda = cell(size(lambda1));
                    for i=1:length(lambda1)
                        opti_lambda{i} = [lambda1{i}, lambda2{i}];
                    end
                    
                    rerp_result.gridsearch.lambda_range = get_new_lambda(rerp_result.gridsearch.lambda_range, opti_lambda);
                    level=level+1;
                end
            end
            
            %Setup grids
            lambda_grid_L1 = cell(size(size(rerp_result.gridsearch.lambda_range)));
            lambda_grid_L2 = cell(size(size(rerp_result.gridsearch.lambda_range)));
            for j=1:length(rerp_result.gridsearch.lambda_range)
                this_lambda_range = rerp_result.gridsearch.lambda_range{j};
                [lambda_grid_L1{j}, lambda_grid_L2{j}] = meshgrid(this_lambda_range(:,1), this_lambda_range(:, 2));
            end
            
            %Execute analysis for each lambda
            grid_results = cell(1,numel(rerp_result.gridsearch.lambda_range{1}));
            for j=1:numel(lambda_grid_L1{1})
                fprintf('rerp: elastic net grid search level %d, lambda %d/%d, time=%f\n', level, j, numel(lambda_grid_L1{1}), toc);
                
                this_lambda = cell(size(lambda_grid_L1));
                
                for k=1:length(this_lambda)
                    this_lambda{k} = [lambda_grid_L1{k}(j), lambda_grid_L1{k}(j)];
                end
                
                this_rerp_result= RerpResult({});
                this_rerp_result.analysis_name= 'Elastic net';
                this_rerp_result.lambda= this_lambda;
                grid_results{j}= rerp_regularized(this_rerp_result);
            end
        end
        
        if strcmp(rerp_result.analysis_name, 'L1 norm')
            %Execute analysis for each lambda
            grid_results = cell(1,length(rerp_result.gridsearch.lambda_range{1}));
            for j=1:numel(rerp_result.gridsearch.lambda_range{1})
                fprintf('rerp: L1 norm grid search level %d, lambda %d/%d, time=%f\n' , level, j, numel(rerp_result.gridsearch.lambda_range{1}),toc);
                
                this_lambda = cell(size(rerp_result.gridsearch.lambda_range));
                
                for k=1:length(this_lambda)
                    this_lambda{k} = rerp_result.gridsearch.lambda_range{k}(j);
                end
                
                this_rerp_result = RerpResult({});
                this_rerp_result.analysis_name='L1 norm';
                this_rerp_result.lambda = this_lambda;
                grid_results{j} = rerp_regularized(this_rerp_result);
            end
        end
        
        if strcmp(rerp_result.analysis_name, 'L2 norm')
            %Execute analysis for each lambda
            grid_results = cell(1,length(rerp_result.gridsearch.lambda_range{1}));
            for j=1:numel(rerp_result.gridsearch.lambda_range{1})
                fprintf('rerp: L2 norm grid search level %d, lambda %d/%d, time=%f\n' , level, j, numel(rerp_result.gridsearch.lambda_range{1}),toc);
                this_lambda = cell(size(rerp_result.gridsearch.lambda_range));
                
                for k=1:length(this_lambda)
                    this_lambda{k} = rerp_result.gridsearch.lambda_range{k}(j);
                end
                
                this_rerp_result = RerpResult({});
                this_rerp_result.analysis_name='L2 norm';
                this_rerp_result.lambda = this_lambda;
                grid_results{j} = rerp_regularized(this_rerp_result);
            end
        end
        
        % Extract the highest rsquare result for each channel/IC
        fprintf('rerp: grid search level %d, finding optimal lambda, time=%f\n', level, toc);
        rerp_result.gridsearch.grid_results=grid_results;
        rsq = cellfun(@(x) x.average_total_rsquare', grid_results,'UniformOutput',0);
        [~, maxidx] = max(cell2mat(rsq),[], 2);
        optimal_results = grid_results(maxidx);
        
        % Assign the optimal values found at this level to the result
        fprintf('rerp: grid search level %d, consolidating optimal results, time=%f\n' , level, toc);
        for i=1:dim
            rerp_result.rerp_estimate(:,i) = optimal_results{i}.rerp_estimate(:,i);
            rerp_result.admm_residual(:,i) = optimal_results{i}.admm_residual(:,i);
            rerp_result.average_total_rsquare(i) = optimal_results{i}.average_total_rsquare(i);
            rerp_result.lambda{i} = optimal_results{i}.lambda{i};
            
            for j=1:s.num_xvalidation_folds
                rerp_result.total_xval_folds(j).noise_variance(i) = optimal_results{i}.total_xval_folds(j).noise_variance(i);
                rerp_result.total_xval_folds(j).data_variance(i) = optimal_results{i}.total_xval_folds(j).data_variance(i);
                rerp_result.total_xval_folds(j).num_samples(i) = optimal_results{i}.total_xval_folds(j).num_samples(i);
            end
            
        end
        
        % Either zoom in another level or return the optimal result found
        % at this level
        if level < s.num_grid_zoom_levels
            fprintf('rerp: grid search level %d, entering next level, time=%f\n',level, toc);
            fprintf('rerp: grid search level %d, optimal lambda, time=%f\n', level, toc);
            
            next_rerp_result = RerpResult({});
            next_rerp_result.analysis_name = rerp_result.analysis_name;
            next_rerp_result.gridsearch.lambda_range = get_new_lambda(rerp_result.gridsearch.lambda_range, rerp_result.lambda);
            [rerp_result.gridsearch.next_stage, optimal_result] = rerp_grid_search(next_rerp_result, level+1);
        else
            rerp_result.gridsearch.next_stage = {};
            optimal_result = rerp_result;
        end
        
        if mean(rerp_result.average_total_rsquare) > mean(optimal_result.average_total_rsquare)
            optimal_result= rerp_result;
        end
    end


% Do regularized regression for the specified penalty functions (one value
% of lambda per penalty per channel). Called without fold number
% externally.
    function rerp_result = rerp_regularized(rerp_result, fold_num)
        
        assert(iscell(rerp_result.lambda)&&(length(rerp_result.lambda)==size(data,2)),'rerp: invalid lambda, must be cell array same dim as data');
        
        if ~exist('fold_num', 'var')
            % If we didn't use fold number, that means we called the function
            % from outside and should crossvalidate the result
            xval=1;
            this_predictor= predictor;
            this_data= data;
            this_P= P;
            this_q= q;
            this_L=L;
            
        else
            % Indicates internal cross validation call
            xval=0;
            this_predictor = predictor(xval_train_idx{fold_num}, :);
            this_data = data(xval_train_idx{fold_num}, :);
            this_P= xval_P{fold_num};
            this_q= xval_q{fold_num};
            if ~isempty(xval_L)
                this_L= xval_L{fold_num};
            else
                this_L=L;
            end
        end
        
        %Configure hot start parameters, if not provided
        if isempty(rerp_result.rerp_estimate) || isempty(rerp_result.admm_residual)
            rerp_result.rerp_estimate = zeros(size(predictor,2), dim);
            rerp_result.admm_residual = zeros(size(predictor,2), dim);
        end
        
        %Elastic net
        if strcmp(rerp_result.analysis_name, 'Elastic net')
            assert(length(rerp_result.lambda{1})==2,'rerp: lambda must be a vector of length 2 for elastic net regularization');
            [rerp_result.rerp_estimate, rerp_result.admm_residual] = rerp_elastic_net(this_predictor, this_data, this_q, rerp_result.lambda, rerp_result.rerp_estimate, rerp_result.admm_residual, this_L);
        end
        
        %L1 norm
        if strcmp(rerp_result.analysis_name, 'L1 norm')
            assert(length(rerp_result.lambda{1})==1,'rerp: lambda must be a vector of length 1 for L1 norm regularization');
            [rerp_result.rerp_estimate, rerp_result.admm_residual] = rerp_L1_norm(this_predictor, this_data, this_q, rerp_result.lambda, rerp_result.rerp_estimate, rerp_result.admm_residual, this_L);
        end
        
        %L2 norm
        if strcmp(rerp_result.analysis_name, 'L2 norm')
            assert(length(rerp_result.lambda{1})==1,'rerp: lambda must be a vector of length 1 for L2 norm regularization');
            rerp_result.rerp_estimate = rerp_L2_norm(this_P, this_q, rerp_result.lambda);
        end
        
        %Least squares
        if strcmp(rerp_result.analysis_name, 'Least squares')
            rerp_result.rerp_estimate = rerp_least_squares(this_P, this_q);
        end
        
        % Cross validate the resultant estimate
        if xval
            rerp_result = cross_validate(rerp_result);
        end
        
    end

% ADMM function (Boyd) for L1 norm regularization (LASSO)
    function [rerp_estimate, admm_residual] = rerp_L1_norm(predictor, data, q, lambda, x, u, L)
        rerp_estimate = zeros(size(predictor,2), size(data, 2));
        admm_residual = zeros(size(predictor,2), size(data, 2));
        
        for i=1:size(data, 2)
            [rerp_estimate(:,i), admm_residual(:,i)] = rerp_dependencies.admm_lasso_boyd(predictor, data(:,i), q(:,i), lambda{i}, x(:,i), u(:,i), L);
        end
    end

%Closed form solution to L2 norm regularization (ridge)
    function rerp_estimate = rerp_L2_norm(P, q, lambda)
        rerp_estimate = zeros(size(q));
        I=speye(size(P));
        
        for i=1:size(rerp_estimate,2)
            rerp_estimate(:,i) = (P+lambda{i}*I)\q(:,i);
            
            if nnz(rerp_estimate(:,i))==0
                %Could not estimate wih mldivide, try pinv instead
                rerp_estimate(:,i) = pinv(P+lambda{i}*I, eps)*q(:,i);
                
                if nnz(rerp_estimate(:,i))==0
                    %Failed to estimate the ERP with this lambda, ignore
                    %that lambda
                    rerp_estimate(:,i) = NaN(size(rerp_estimate(:,i)));
                end
            end
            
        end
    end

% ADMM function (Boyd) for L1 + L2 norm regularization (elastic net)
    function [rerp_estimate, admm_residual] = rerp_elastic_net(predictor, data, q, lambda, x, u, L)
        rerp_estimate = zeros(size(predictor,2), size(data, 2));
        admm_residual = zeros(size(predictor,2), size(data, 2));
        for i=1:size(rerp_estimate ,2)
            [rerp_estimate(:,i), admm_residual(:,i)] = rerp_dependencies.admm_elastic_net_boyd(predictor, data(:,i), q(:,i), lambda{i}, x(:,i), u(:,i), L);
        end
    end

% Closed form solution to least squares
    function rerp_estimate = rerp_least_squares(P, q)
        rerp_estimate = P\q;
        
        pinv_P=[]; 
        for i=1:size(rerp_estimate,2)
            if nnz(rerp_estimate(:,i))==0
                %Could not estimate wih mldivide, try pinv instead
                if isempty(pinv_P)
                    pinv_P=pinv(P);
                end
                
                rerp_estimate(:,i) = pinv(P)*q(:,i);

                if nnz(rerp_estimate(:,i))==0
                    %Failed to estimate the ERP, ignore it
                    rerp_estimate(:,i) = NaN(size(rerp_estimate(:,i)));
                end
            end
        end
    end

    function rerp_result = cross_validate(rerp_result)
        %Cross validate a single rerp for the data as a whole, as well as
        %individual variables.
        import rerp_dependencies.meanw;
        import rerp_dependencies.RerpXvalFold;
        
        fprintf('rerp: cross validating parameters\n');
        for i=1:s.num_xvalidation_folds
            this_fold=RerpXvalFold;
            rerp_result.total_xval_folds(i)=this_fold;
            
            fprintf('rerp: %s fold %d / %d, time=%f\n', rerp_result.analysis_name, i, s.num_xvalidation_folds,toc);
            
            % Create copy so we keep the whole data rerp estimate
            res_copy = copy(rerp_result);
            res_copy = rerp_regularized(res_copy, i);
            this_test_predictor = predictor(xval_test_idx{i},:);
            this_test_data = data(xval_test_idx{i},:);
            
            %Remove any time index that does not intersect with this
            %section of predictor.
            keep_idx = sum(this_test_predictor,2)~=0;
            
            this_data_model = this_test_predictor(keep_idx,:)*res_copy.rerp_estimate;
            this_noise = this_test_data(keep_idx,:) - this_data_model;
            
            % Get the statistics of the whole data
            this_fold.noise_variance = var(this_noise)';
            this_fold.data_variance = var(this_test_data(keep_idx, :))';
            this_fold.num_samples = repmat(size(this_noise, 1), [size(this_noise, 2) 1]);
        end
        
        % Calculate rsquare for each channel or IC
        noisevar = [rerp_result.total_xval_folds(:).noise_variance];
        datavar= [rerp_result.total_xval_folds(:).data_variance];
        w= [rerp_result.total_xval_folds(:).num_samples];
        
        rerp_result.average_total_rsquare = 1-meanw(noisevar./datavar, w, 2)';
    end


    function rerp_result = cross_validate_event_types(rerp_result)
        % Do final cross validation of each event type seperately. Called after we
        % have found the optimal lambda
        import rerp_dependencies.meanw;
        import rerp_dependencies.RerpXvalFold;

        fprintf('rerp: cross validating event types\n');
        block=zeros(length(parameter_idx_layout), size(data,2));
        for i=1:s.num_xvalidation_folds
            this_fold=RerpXvalFold;
            rerp_result.event_xval_folds(i)=this_fold;
            
            fprintf('rerp: %s fold %d / %d, time=%f\n', rerp_result.analysis_name, i, s.num_xvalidation_folds, toc);
            rerp_result.event_xval_folds(i).noise_variance=block;
            rerp_result.event_xval_folds(i).data_variance=block;
            
            % Create copy so we keep the whole data rerp estimate
            res_copy = copy(rerp_result);
            res_copy = rerp_regularized(res_copy, i);
            
            for j=1:length(parameter_idx_layout)
                this_test_predictor = predictor(xval_test_idx{i}, parameter_idx_layout{j});
                this_test_data = data(xval_test_idx{i},:);
                
                %Remove any time index that does not intersect with this
                %event type.
                keep_idx = sum(this_test_predictor, 2)~=0;
                
                this_data_model = this_test_predictor(keep_idx,:)*res_copy.rerp_estimate(parameter_idx_layout{j}, :);
                this_noise = this_test_data(keep_idx,:) - this_data_model;
                
                % Get the statistics for this event type
                this_fold.noise_variance(j,:) = var(this_noise);
                this_fold.data_variance(j,:) = var(this_test_data(keep_idx, :));
                this_fold.num_samples(j,:) = repmat(size(this_noise, 1), [1 size(this_noise, 2)]);
            end
        end
        
        % Calculate rsquare for each channel or IC
        noisevar = reshape([rerp_result.event_xval_folds(:).noise_variance],[length(parameter_idx_layout), size(data,2), s.num_xvalidation_folds]);
        datavar= reshape([rerp_result.event_xval_folds(:).data_variance],[length(parameter_idx_layout), size(data,2), s.num_xvalidation_folds]);
        w= reshape([rerp_result.event_xval_folds(:).num_samples],[length(parameter_idx_layout),size(data,2), s.num_xvalidation_folds]);

        rerp_result.average_event_rsquare = 1-meanw(noisevar./datavar, w, 3);
    end

%Cache factors for cross-validation
    function setup_xval_predictors
        chunk_size = floor(size(predictor,1)/s.num_xvalidation_folds);
        for i=1:s.num_xvalidation_folds
            start_idx = (i-1)*chunk_size +1;
            end_idx = start_idx+chunk_size;
            test_indexes = start_idx:(end_idx-1);
            train_indexes = setdiff(1:size(predictor,1), test_indexes);
            
            xval_train_idx{i}=train_indexes;
            xval_test_idx{i}=test_indexes;
            xval_P{i}=predictor(train_indexes,:)'*predictor(train_indexes,:);
            xval_q{i}=predictor(train_indexes,:)'*data(train_indexes,:);
            
            % Need a lot of memory for ADMM functions
            if ~isempty(L)
                xval_L{i}=sparse(chol(xval_P{i} + RHO*speye(size(xval_P{i})), 'lower' ));
            end
        end
    end

% When zooming, get the next grid of lambda (if mesh==1, using meshgrid, sqrt the number of pts)
    function new_lambda = get_new_lambda(lambda, optimal_lambda)
        new_lambda = cell(size(lambda));
        
        % Go through lambda for each time series
        for i=1:length(lambda)
            this_lambda = lambda{i};
            this_optimal_lambda = optimal_lambda{i};
            
            % Go through multivariate lambda (L1, L2, etc)
            for j=1:size(this_lambda, 2)
                this_optimal_lambda_col = this_optimal_lambda(j);
                this_lambda_col = this_lambda(:,j);
                
                %Find the index of the lambda_range which is closest to optimal lambda
                [~, idx] = min(abs(this_lambda_col - this_optimal_lambda_col));
                
                if idx-1 < 1
                    lamstart = this_lambda_col(idx) - this_lambda_col(idx+1);
                else
                    lamstart = this_lambda_col(idx-1);
                end
                
                if idx+1 > length(this_lambda_col)
                    lamend = this_lambda_col(idx) + this_lambda_col(idx-1);
                else
                    lamend = this_lambda_col(idx+1);
                end
                
                newlam = linspace(lamstart, lamend, s.num_grid_points)';
                new_lambda{i} = [new_lambda{i} newlam];
            end
        end
    end

end

