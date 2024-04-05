%Contains all results from rerp function, and methods to plot the results.
%   Usage:
%       rerp_result_gui(rerp_results);
%           Opens result plotting GUI (recommended)
%
%   Parameters:
%       rerp_results:
%           Vector of RerpResult objects
%
%   See doc RerpResult for extensive information on methods and properties.
%
%   Also see:
%       pop_rerp, rerp, RerpProfile, RerpResultStudy, rerp_result_gui
classdef RerpResult < matlab.mixin.Copyable
    properties
        rerp_profile % Profile used to derive this result
        analysis_name % A title assigned by the rerp function
        
        date_completed %Exact date and time analysis was completed
        compute_time_seconds=0; %How many compute seconds this result took
        
        lambda=[] % The optimal lambda if grid search was used, or the only lambda otherwise
        ersp_flag=0 % 1 if this result was computed from time-frequency decomposed data
        
        rerp_estimate % The rerp or rersp esitmates for each channel and frequency. must be reshaped for ersp.
        admm_residual % If regularization required alternating direction method of multipliers
        
        tags={};
        parameter_idx_layout={};  %Gives the index into the parameter vector for each variable
        
        average_total_rsquare % Rsquare as computed on the entire data with full signal estimate
        average_event_rsquare % Rsquare as computed only on the part of the data affected by the variable, with it's part of the signal estimate
        
        total_xval_folds % Cross validation structures for each fold
        event_xval_folds % Cross validation structures for each fold
        
        gridsearch % A complete history of the grid search process, if any
        
        rerp_plot_spec %RerpPlotSpec object
        name %Name of file this result was pulled from
    end
    
    methods
        function obj = RerpResult(rerp_profile, rerp_plot_spec)
            import rerp_dependencies.RerpPlotSpec
            import rerp_dependencies.RerpXvalFold
            
            if nargin < 1
                help RerpResult
                return;
            end
            
            obj.total_xval_folds=RerpXvalFold;
            obj.event_xval_folds=RerpXvalFold;
            
            obj.rerp_profile = rerp_profile;
            if nargin > 1
                obj.rerp_plot_spec=rerp_plot_spec;
            else
                obj.rerp_plot_spec=RerpPlotSpec;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Plotting
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function plotRerpEventTypes(obj, h)
            %Plot rerp estimates for time series on same axis (one plot per event type or hed tag).
            %Option for only plotting ERP waveforms which correspond to
            %statistically significant Rsquare at specified p-value threshold.
            %    Usage:
            %        rerp_result.plotRerpEventTypes(rerp_plot_spec.event_idx, rerp_plot_spec.ts_idx);
            %        rerp_result.plotRerpEventTypes(rerp_plot_spec.event_idx, rerp_plot_spec.ts_idx, rerp_plot_spec.exclude_insignificant);
            %        rerp_result.plotRerpEventTypes(rerp_plot_spec.event_idx, rerp_plot_spec.ts_idx, rerp_plot_spec.exclude_insignificant, rerp_plot_spec.significance_level);
            %        rerp_result.plotRerpEventTypes(rerp_plot_spec.event_idx, rerp_plot_spec.ts_idx, rerp_plot_spec.exclude_insignificant, rerp_plot_spec.significance_level, h);
            %    Parameters:
            %        rerp_plot_spec.exclude_insignificant [true false] - whether to skip plotting statistically insignificant ERP estimates
            %        rerp_plot_spec.significance_level - p-value threshold for determining significance
            import rerp_dependencies.*
            
            if ~exist('h','var')
                h=figure;
            end
            
            if obj.rerp_plot_spec.exclude_insignificant
                significance_label = [' ( significant @ p < ' num2str(obj.rerp_plot_spec.significance_level) ' )'];
            else
                significance_label = '';
            end
            
            rsquare_significance = obj.get_event_rsquare_significance;
            
            hold all;
            assert(~obj.ersp_flag, 'RerpResult: plotRerpEventTypes is invalid, use plotRersp instead');
            
            [tags, estimates, xaxis_ms] = obj.get_plotting_params;
            
            x_label='time (ms)';
            y_label='amplitude (RMS microvolt)';
            
            if isempty(obj.rerp_plot_spec.event_idx)
                obj.rerp_plot_spec.event_idx = 1:length(tags);
            end
            
            if isempty(obj.rerp_plot_spec.ts_idx)
                obj.rerp_plot_spec.ts_idx = 1:size(estimates,2);
            end
            
            datasetname = regexp(obj.rerp_profile.eeglab_dataset_name,'.*[\\\/](.*).set','tokens');
            if ~isempty(datasetname)
                datasetname = {{regexprep(datasetname{1}{1},'[\_]','\\\_')}};
            else
                datasetname = {{''}};
            end
            
            if obj.rerp_plot_spec.constant_scale
                max_y = 1.1*max(max([estimates{obj.rerp_plot_spec.event_idx, obj.rerp_plot_spec.ts_idx}]));
                min_y = 1.1*min(min([estimates{obj.rerp_plot_spec.event_idx, obj.rerp_plot_spec.ts_idx}]));
            end
            
            m=1;
            props=get(findobj(h, 'tag', 'legend'));
            for i=obj.rerp_plot_spec.event_idx
                
                if obj.rerp_plot_spec.exclude_insignificant
                    this_ts_idx = obj.rerp_plot_spec.ts_idx(rsquare_significance(i, obj.rerp_plot_spec.ts_idx)==1);
                else
                    this_ts_idx = obj.rerp_plot_spec.ts_idx;
                end
                
                if ~isempty(obj.rerp_plot_spec.ts_idx)
                    scrollsubplot(4,1,m,h);
                    hold all;
                    this_estimate = [estimates{i, this_ts_idx}];
                    if ~isempty(this_estimate)
                        plot(xaxis_ms{i}', this_estimate);
                    end
                    
                    hcmenu = uicontextmenu;
                    uimenu(hcmenu, 'Label', 'Publish graph', 'Callback', @RerpResult.gui_publish);
                    set(gca,'uicontextmenu', hcmenu);
                    
                    if obj.rerp_plot_spec.constant_scale
                        ylim([min_y, max_y]);
                    end
                    
                    if obj.rerp_profile.settings.hed_enable
                        titl = ['Tag:' tags{i}];
                    else
                        titl = ['Event type:' tags{i}];
                    end
                    
                    if obj.rerp_profile.settings.type_proc
                        leg = [datasetname{1}{1} ' - ' obj.analysis_name ', Channel: '];
                        ts_label=obj.rerp_profile.include_chans(this_ts_idx);
                    else
                        leg = [datasetname{1}{1} ' - ' obj.analysis_name ', Component: '];
                        ts_label=obj.rerp_profile.include_comps(this_ts_idx);
                    end
                    
                    title([titl significance_label]);
                    xlim([min(xaxis_ms{i}) max(xaxis_ms{i})]);
                    
                    legend_idx = cellfun(@(x) [leg x], regexp(num2str(ts_label),'\s*','split'), 'UniformOutput' ,false);
                    
                    if isempty(props)
                        a=legend(legend_idx);
                        leg_prop = get(a);
                        leg_prop.UserData.lstrings=legend_idx;
                        set(a,'UserData', leg_prop.UserData);
                    else
                        new_ls = [props(m).UserData.lstrings legend_idx];
                        a=legend(new_ls);
                        leg_prop = get(a);
                        leg_prop.UserData.lstrings={new_ls};
                        set(a,'UserData', leg_prop.UserData);
                    end
                    
                    pr= get(gca,'UserData');
                    pr.legend=a;
                    set(gca, 'UserData',pr);
                    
                    xlabel(x_label);
                    ylabel(y_label);
                end
                m=m+1;
            end
        end
        
        function plotRerpTimeSeries(obj, h)
            %Plot rerp estimates for event types on same axis (one plot per time series).
            %Option for only plotting ERP waveforms which correspond to
            %statistically significant Rsquare at specified p-value threshold.
            %    Usage:
            %        rerp_result.plotRerpTimeSeries(rerp_plot_spec.event_idx, rerp_plot_spec.ts_idx);
            %        rerp_result.plotRerpTimeSeries(rerp_plot_spec.event_idx, rerp_plot_spec.ts_idx, rerp_plot_spec.exclude_insignificant);
            %        rerp_result.plotRerpTimeSeries(rerp_plot_spec.event_idx, rerp_plot_spec.ts_idx, rerp_plot_spec.exclude_insignificant, rerp_plot_spec.significance_level);
            %        rerp_result.plotRerpTimeSeries(rerp_plot_spec.event_idx, rerp_plot_spec.ts_idx, rerp_plot_spec.exclude_insignificant, rerp_plot_spec.significance_level, h);
            %    Parameters:
            %        rerp_plot_spec.exclude_insignificant [true false] - whether to skip plotting statistically insignificant ERP estimates
            %        rerp_plot_spec.significance_level - p-value threshold for determining significance
            import rerp_dependencies.*
            
            if ~exist('h','var')
                h=figure;
            end
            
            if obj.rerp_plot_spec.exclude_insignificant
                significance_label = [' ( significant @ p < ' num2str(obj.rerp_plot_spec.significance_level) ' )'];
            else
                significance_label = '';
            end
            
            rsquare_significance = obj.get_event_rsquare_significance;
            assert(~obj.ersp_flag, 'RerpResult: plotRerpTimeSeries is invalid, , use plotRersp instead');
            
            [tags, estimates, xaxis_ms] = obj.get_plotting_params;
            
            x_label='time (ms)';
            y_label='amplitude (RMS microvolt)';
            
            if isempty(obj.rerp_plot_spec.event_idx)
                obj.rerp_plot_spec.event_idx = 1:length(tags);
            end
            
            if isempty(obj.rerp_plot_spec.ts_idx)
                obj.rerp_plot_spec.ts_idx = 1:size(estimates,2);
            end
            
            datasetname = regexp(obj.rerp_profile.eeglab_dataset_name,'.*[\\\/](.*).set','tokens');
            datasetname = {{regexprep(datasetname{1}{1},'[\_]','\\\_')}};
            
            if obj.rerp_plot_spec.constant_scale
                max_y = 1.1*max(max([estimates{obj.rerp_plot_spec.event_idx, obj.rerp_plot_spec.ts_idx}]));
                min_y = 1.1*min(min([estimates{obj.rerp_plot_spec.event_idx, obj.rerp_plot_spec.ts_idx}]));
            end
            
            m=1;
            props=get(findobj(h, 'tag', 'legend'));
            new_idx=[];
            for i=obj.rerp_plot_spec.ts_idx
                scrollsubplot(4,1,m,h);
                
                n=1;
                for j=obj.rerp_plot_spec.event_idx
                    if ~obj.rerp_plot_spec.exclude_insignificant||rsquare_significance(j, i)
                        plot(xaxis_ms{j}, estimates{j, i});
                        new_idx(n)=j;
                        hold all;
                        n=n+1;
                    end
                end
                
                if obj.rerp_plot_spec.constant_scale
                    ylim([min_y, max_y]);
                end
                
                if ~isempty(obj.rerp_plot_spec.event_idx)
                    hcmenu = uicontextmenu;
                    uimenu(hcmenu, 'Label', 'Publish graph', 'Callback', @RerpResult.gui_publish);
                    set(gca,'uicontextmenu', hcmenu)
                    
                    if obj.rerp_profile.settings.type_proc
                        titl = ['Channel:' num2str(obj.rerp_profile.include_chans(i))];
                    else
                        titl = ['Component:' num2str(obj.rerp_profile.include_comps(i))];
                    end
                    
                    if obj.rerp_profile.settings.hed_enable
                        leg = [datasetname{1}{1} ' - ' obj.analysis_name ', Tag: '];
                    else
                        leg = [datasetname{1}{1} ' - ' obj.analysis_name ', Event type: '];
                    end
                    
                    title([titl significance_label]);
                    xlim([min(cell2mat(xaxis_ms(new_idx))) max(cell2mat(xaxis_ms(new_idx)))]);
                    legend_idx = cellfun(@(x) [leg x], tags(new_idx) , 'UniformOutput' ,false);
                    
                    if isempty(props)
                        a=legend(strtrim(legend_idx));
                        leg_prop = get(a);
                        leg_prop.UserData.lstrings=legend_idx;
                        set(a,'UserData', leg_prop.UserData);
                    else
                        new_ls = [props(m).UserData.lstrings legend_idx];
                        a=legend(new_ls);
                        leg_prop = get(a);
                        leg_prop.UserData.lstrings={new_ls};
                        set(a,'UserData', leg_prop.UserData);
                    end
                    
                    pr= get(gca,'UserData');
                    pr.legend=a;
                    set(gca, 'UserData',pr);
                    
                    xlabel(x_label);
                    ylabel(y_label);
                    m=m+1;
                end
            end
        end
        
        function plotRerpTotalRsquared(obj, h)
            %Plot the average rsquare for the time series as a whole with significance markings.
            %Time series are ranked on average rsquare, with ttest taken across xvalidation folds.
            %   Usage:
            %       rerp_result.plotRerpTotalRsquared(rerp_plot_spec.ts_idx, rerp_plot_spec.significance_level);
            %       rerp_result.plotRerpTotalRsquared(rerp_plot_spec.ts_idx, rerp_plot_spec.significance_level, h);
            import rerp_dependencies.*
            
            if ~exist('h','var')
                h=figure;
            end
            
            vals = obj.average_total_rsquare(obj.rerp_plot_spec.ts_idx);
            
            tmax = max(vals);
            tmin = min(min(vals), 0);
            
            if ~isempty(obj.rerp_plot_spec.significance_level)
                rsquare_significance = obj.get_total_rsquare_significance;
                rsquare_significance = rsquare_significance(obj.rerp_plot_spec.ts_idx);
            end
            
            %Plot average total R2
            p=plot(0:(length(obj.rerp_plot_spec.ts_idx)-1), vals);
            line_props = get(p);
            set(gca,'xtickmode','manual');
            set(gca, 'xtick', 0:(length(obj.rerp_plot_spec.ts_idx)-1));
            
            %Figure out the legend
            datasetname = regexp(obj.rerp_profile.eeglab_dataset_name,'.*[\\\/](.*).set','tokens');
            if ~isempty(datasetname)
                datasetname = {{regexprep(datasetname{1}{1},'[\_]','\\\_')}};
                legend_idx=[datasetname{1}{1} ' - ' obj.analysis_name];
                
            else
                legend_idx= obj.analysis_name;
            end
            
            %We save handles to the plots in the legend so we can
            %rereference them later when doing overplotting.
            props=get(findobj(h, 'tag', 'legend'));
            if isempty(props)
                leg = legend(strtrim(legend_idx));
                props = get(leg);
                props.UserData.plotHandles = p;
                props.UserData.lstrings={legend_idx};
                set(gca, 'xticklabel', obj.rerp_plot_spec.ts_idx);
                
            else
                plotHandles = [props.UserData.plotHandles p];
                leg = legend(plotHandles, {props.UserData.lstrings{:} strtrim(legend_idx)});
                props = get(leg);
                props.UserData.plotHandles = plotHandles;
                set(gca, 'xticklabel', 1:length(obj.rerp_plot_spec.ts_idx));
            end
            
            %Plot statistical significance markers
            sig_plot=[];
            for j=1:length(obj.rerp_plot_spec.ts_idx)
                if rsquare_significance(j)
                    hold all;
                    sig_plot = plot(j-1, vals(j) ,'s', 'LineWidth', 1, 'MarkerEdgeColor',line_props.Color,'MarkerSize', 10);
                end
            end
            
            %Put significance in legend
            if ~isempty(sig_plot)
                plotHandles = [props.UserData.plotHandles sig_plot];
                leg = legend(plotHandles, {props.UserData.lstrings{:} ['significant @ p < ' num2str(obj.rerp_plot_spec.significance_level)]});
                props = get(leg);
                props.UserData.plotHandles = plotHandles;
            end
            
            %Restore the information in the handle object
            set(leg, 'UserData', props.UserData);
            pr= get(gca,'UserData');
            
            %Display the new legend in the figure
            pr.legend=leg;
            set(gca, 'UserData',pr);
            
            hcmenu = uicontextmenu;
            uimenu(hcmenu, 'Label', 'Publish graph', 'Callback', @RerpResult.gui_publish);
            set(gca,'uicontextmenu', hcmenu);
            set(gca, 'ygrid', 'on');
            
            if obj.rerp_profile.settings.type_proc
                type= 'Channel';
            else
                type= 'Component';
            end
            
            xlabel([type ' - decreasing R ^2 order']);
            ylabel('R ^2');
            title('Rsquare performance by time series');
            
            a = get(gca);
            tmin = min(tmin, a.YLim(1));
            tmax = max(tmax, a.YLim(2));
            
            axis([-1 length(obj.rerp_plot_spec.ts_idx) tmin tmax]);
            line([-1 length(obj.rerp_plot_spec.ts_idx)],[0 0],'linewidth',1,'color',[0 0 0]);
            
        end
        
        function plotRerpEventRsquared(obj, h)
            %Plot the average rsquare taking into account only specific event types, with significance markings.
            %Time series are ranked on average rsquare, with ttest taken across xvalidation folds.
            %   Usage:
            %       rerp_result.plotRerpEventRsquared(rerp_plot_spec.ts_idx, rerp_plot_spec.significance_level, rerp_plot_spec.event_idx);
            %       rerp_result.plotRerpEventRsquared(rerp_plot_spec.ts_idx, rerp_plot_spec.significance_level, rerp_plot_spec.event_idx, h);
            import rerp_dependencies.*
            
            if ~exist('h','var')
                h=figure;
            end
            
            hold all;
            tags=obj.get_plotting_params;
            
            if isempty(obj.rerp_plot_spec.event_idx)
                obj.rerp_plot_spec.event_idx = 1:size(obj.average_event_rsquare, 1);
            end
            
            if ~isempty(obj.rerp_plot_spec.significance_level)
                rsquare_significance = obj.get_event_rsquare_significance;
            end
            
            datasetname = regexp(obj.rerp_profile.eeglab_dataset_name,'.*[\\\/](.*).set','tokens');
            datasetname = {{regexprep(datasetname{1}{1},'[\_]','\\\_')}};
            
            vals=obj.average_event_rsquare(obj.rerp_plot_spec.event_idx, obj.rerp_plot_spec.ts_idx_event_types(:));
            overall_max = max(max(vals));
            overall_min = min(min(min(vals), 0));
            
            m=1;
            for i=1:length(obj.rerp_plot_spec.event_idx)
                
                vals = obj.average_event_rsquare(obj.rerp_plot_spec.event_idx(i), obj.rerp_plot_spec.ts_idx_event_types(i,:));
                this_rsquare_significance = rsquare_significance(i, obj.rerp_plot_spec.ts_idx_event_types(i,:));
                
                tmax = max(vals);
                tmin = min(min(vals), 0);
                
                hold all;
                scrollsubplot(3,1,m,h);
                
                p=plot(0:(size(obj.rerp_plot_spec.ts_idx_event_types,2)-1), vals);
                line_props = get(p);
                set(gca,'xtickmode','manual');
                set(gca, 'xtick', 0:(size(obj.rerp_plot_spec.ts_idx_event_types,2)-1));
                props=get(findobj(h,'Tag', ['legend_' num2str(i)]));
                legend_idx=[datasetname{1}{1} ' - ' obj.analysis_name];
                
                if isempty(props)
                    leg = legend(legend_idx);
                    props = get(leg);
                    props.UserData.plotHandles = p;
                    props.UserData.lstrings={legend_idx};
                    set(leg,'UserData', props.UserData, 'Tag',['legend_' num2str(i)]);
                    set(gca, 'xticklabel', obj.rerp_plot_spec.ts_idx_event_types(i,:));
                else
                    plotHandles = [props.UserData.plotHandles p];
                    leg = legend(plotHandles, {props.UserData.lstrings{:} legend_idx});
                    props = get(leg);
                    props.UserData.plotHandles = plotHandles;
                    set(leg,'UserData', props.UserData);
                    set(gca, 'xticklabel', 1:length(obj.rerp_plot_spec.ts_idx_event_types(i,:)));
                end
                
                %Plot significance markers
                sig_plot = [];
                for j=1:length(this_rsquare_significance)
                    if this_rsquare_significance(j)
                        hold all;
                        sig_plot=plot(j-1, vals(j) , 's', 'LineWidth', 1, 'MarkerEdgeColor',line_props.Color,'MarkerSize', 10);
                    end
                end
                
                %Put significance in legend
                if ~isempty(sig_plot)
                    plotHandles = [props.UserData.plotHandles sig_plot];
                    leg = legend(plotHandles, {props.UserData.lstrings{:} ['significant @ p < ' num2str(obj.rerp_plot_spec.significance_level)]});
                    props = get(leg);
                    props.UserData.plotHandles = plotHandles;
                end
                
                set(leg,'UserData', props.UserData);
                pr= get(gca,'UserData');
                pr.legend=leg;
                set(gca, 'UserData',pr);
                set(gca, 'ygrid', 'on');
                
                hcmenu = uicontextmenu;
                uimenu(hcmenu, 'Label', 'Publish graph', 'Callback', @RerpResult.gui_publish);
                set(gca,'uicontextmenu', hcmenu);
                
                if obj.rerp_profile.settings.type_proc
                    type= 'Channel';
                else
                    type= 'Component';
                end
                
                xlabel([type ' - decreasing R ^2 order']);
                ylabel('R ^2');
                title(['Rsquare performance by time series: ' tags{obj.rerp_plot_spec.event_idx(i)}]);
                
                a = get(gca);
                if obj.rerp_plot_spec.constant_scale
                    tmin = overall_min - abs(overall_min)*.2;
                    tmax = overall_max + abs(overall_max)*.2;
                end
                
                tmin = min(tmin, a.YLim(1));
                tmax = max(tmax, a.YLim(2));
                axis_min(i)=tmin;
                axis_max(i)=tmax;
                
                axis([-1 size(obj.rerp_plot_spec.ts_idx_event_types,2) tmin tmax]);
                line([-1 size(obj.rerp_plot_spec.ts_idx_event_types,2)],[0 0],'linewidth',1,'color',[0 0 0]);
                m=m+1;
            end
            
            if obj.rerp_plot_spec.constant_scale
                ymin = min(axis_min);
                ymax = max(axis_max);
                m=1;
                for i=1:length(obj.rerp_plot_spec.event_idx)
                    scrollsubplot(3,1,m,h);
                    ylim([ymin ymax]);
                    m=m+1;
                end
            end
        end
        
        function plotRerpImage(obj, h)
            %Plot raw epochs, modeled epochs and difference (noise) epochs
            %Sorts epochs based on the time delay to another event type
            %   Usage:
            %     rerp_result.plotRerpImage;
            %
            %   Parameters:
            import rerp_dependencies.*
            
            if ~exist('h','var')
                h=figure;
            end
            
            boundary = obj.rerp_plot_spec.rerp_image_boundary;
            
            eeg_parts=regexp(obj.rerp_profile.eeglab_dataset_name, '(.*[\\\/])(.*.set)', 'tokens');
            
            try
                EEG=pop_loadset('filename', eeg_parts{1}{2}, 'filepath', eeg_parts{1}{1});
            catch
                uiwait(msgbox('RerpResult: could not find the dataset for this result, please locate','.set not found' ,'warn'));
                EEG=pop_loadset;
            end
            
            if (~obj.rerp_profile.settings.type_proc) && isempty(EEG.icaact)
                EEG.icaact=eeg_getica(EEG);
            end
            
            assert(obj.ersp_flag~=1, 'RerpResult: this profile was run on time-fequency data (rERSP); plotRerpImage is invalid');
            
            if isempty(obj.rerp_plot_spec.ts_idx)
                obj.rerp_plot_spec.ts_idx = 1:size(EEG.data,1);
            end
            
            if obj.rerp_profile.settings.type_proc==0
                assert(~isempty(EEG.icaact)&& size(EEG.icaact,3)==1,'RerpResult: this profile is set for ICA; populate EEG.icaact with continuous ICA activations');
                data = EEG.icaact(obj.rerp_profile.include_comps(obj.rerp_plot_spec.ts_idx),:)';
            else
                assert(~isempty(EEG.data) && size(EEG.data,3)==1,'RerpResult: EEG.data must be populated with continuous data');
                data = EEG.data(obj.rerp_profile.include_chans(obj.rerp_plot_spec.ts_idx),:)';
            end
            
            % Replace artifact indexes with data median for plotting
            if obj.rerp_profile.settings.artifact_rejection_enable
                if obj.rerp_profile.settings.artifact_variable_enable
                    data(obj.rerp_profile.variable_artifact_indexes,:)=repmat(median(data), [nnz(obj.rerp_profile.variable_artifact_indexes), 1]);
                else
                    data(obj.rerp_profile.computed_artifact_indexes,:)=repmat(median(data), [nnz(obj.rerp_profile.computed_artifact_indexes), 1]);
                end
            end
            
            [tags, estimates, xaxis_ms, epoch_boundaries] = obj.get_plotting_params;
            locking_tag = tags{obj.rerp_plot_spec.locking_idx};
            this_epoch_boundaries=epoch_boundaries{obj.rerp_plot_spec.locking_idx};
            locking_estimate = estimates(obj.rerp_plot_spec.locking_idx, obj.rerp_plot_spec.ts_idx);
            
            if ~isempty(obj.rerp_plot_spec.delay_idx)
                delay_tag = tags{obj.rerp_plot_spec.delay_idx};
            else
                delay_tag=[];
            end
            
            num_samples = ceil((boundary(2)-boundary(1))*obj.rerp_profile.sample_rate);
            %this_epoch_boundaries = epoch_boundaries{obj.rerp_plot_spec.locking_idx};
            this_xaxis_ms = ((0:(num_samples-1))'/obj.rerp_profile.sample_rate + boundary(1))*1000;
            
            disp('RerpResult: generating modeled data');
            [predictor, data_pad] = obj.rerp_profile.predictor;
            extra_pad = ceil(max(0, this_epoch_boundaries(1)-boundary(1))*obj.rerp_profile.sample_rate);
            front_pad = data_pad(1)+extra_pad;
            %data = [zeros(data_pad(1), size(data,2)); data; zeros(data_pad(2), size(data,2)); zeros(num_samples,size(data,2))];
            data = [zeros(front_pad, size(data,2)); data; zeros(data_pad(2) + num_samples, size(data,2))];
            modeled_data = [zeros(extra_pad, size(data,2)); predictor*obj.rerp_estimate(:,obj.rerp_plot_spec.ts_idx); zeros(num_samples,size(data,2))];
            noise = data-modeled_data;

            disp('RerpResult: getting epochs');
            %idx_start = max(max(ceil(this_epoch_boundaries(1)*obj.rerp_profile.sample_rate)+data_pad(1)), 1);
            [data_epochs, event_nums] = obj.get_rerp_epochs(data, locking_tag, front_pad, boundary);
            [modeled_epochs] = obj.get_rerp_epochs(modeled_data, locking_tag,  front_pad, boundary);
            [noise_epochs] = obj.get_rerp_epochs(noise, locking_tag, front_pad, boundary);
            
            %Threshold epochs to modeled data to get good color
            %range.
            max_model = max(max(modeled_epochs));
            min_model = min(min(modeled_epochs));
            data_epochs(data_epochs > max_model) = max_model;
            data_epochs(data_epochs < min_model) = min_model;
            noise_epochs(noise_epochs > max_model) = max_model;
            noise_epochs(noise_epochs < min_model) = min_model;
            
            if ~isempty(delay_tag)
                disp('RerpResult: calculating order of trials');
                sorting_var = (obj.get_delay_times(event_nums, delay_tag, num_samples)/obj.rerp_profile.sample_rate)*1000;
            else
                sorting_var=[];
            end
            
            m=1;
            for i=1:length(obj.rerp_plot_spec.ts_idx)
                this_ts = obj.rerp_plot_spec.ts_idx(i);
                
                if obj.rerp_profile.settings.type_proc
                    ts = 'Channel';
                    tsn = num2str(obj.rerp_profile.include_chans(this_ts));
                else
                    ts = 'Component';
                    tsn = num2str(obj.rerp_profile.include_comps(this_ts));
                end
                
                if obj.rerp_profile.settings.hed_enable
                    v = 'Tag';
                else
                    v = 'Event type';
                end
                
                % Plot the data epochs
                scrollsubplot(4,1,m,h);
                erpimage(data_epochs(:,:,i), sorting_var, this_xaxis_ms, ['Data epochs - ' v ': ' locking_tag ', ' ts ': ' tsn]);
                
                % Plot the modeled epochs
                scrollsubplot(4,1,m+1,h);
                erpimage(modeled_epochs(:,:,i), sorting_var, this_xaxis_ms, ['Modeled epochs - ' v ': ' locking_tag ', ' ts ': ' tsn]);
                
                % Plot the noise epochs
                scrollsubplot(4,1,m+2,h);
                erpimage(noise_epochs(:,:,i), sorting_var, this_xaxis_ms, ['Difference epochs - ' v ': ' locking_tag ', ' ts ': ' tsn]);
                
                % Plot the rerp estimates
                scrollsubplot(4,1,m+3,h);
                plot(xaxis_ms{obj.rerp_plot_spec.locking_idx}', locking_estimate{i});
                title('rERP estimates');
                xlabel('time (ms)');
                ylabel('epoch number');
                legend([v '(locking): ' locking_tag]);
                xlim([min(this_xaxis_ms) max(this_xaxis_ms)]);
                hcmenu = uicontextmenu;
                uimenu(hcmenu, 'Label', 'Publish graph', 'Callback', @RerpResult.gui_publish);
                set(gca,'uicontextmenu', hcmenu);
                m=m+4;
            end
            
        end
        
        function plotRersp(obj, h)
            %Plot regressed ERSP estimates
            %   Usage:
            %       rerp_result.plotRersp(rerp_plot_spec.event_idx, rerp_plot_spec.ts_idx);
            %       rerp_result.plotRersp(rerp_plot_spec.event_idx, rerp_plot_spec.ts_idx, h);
            import rerp_dependencies.*
            
            assert(obj.ersp_flag==1, 'RerpResult: plotRersp is invalid for this result, use plotRerpEventTypes or plotRerpTimeSeries instead');
            
            if ~exist('h','var')
                figure;
            end
            
            
            nbins = obj.rerp_profile.settings.nbins;
            sr = obj.rerp_profile.sample_rate;
            y_axis_hz = logspace(0,log10(floor((sr/2))), nbins);
            y_axis_hz = num2cell(num2str(y_axis_hz(:)), 2);
            [tags, estimates, xaxis_ms] = obj.get_plotting_params;
            
            m=1;
            for i=obj.rerp_plot_spec.ts_idx
                if obj.rerp_profile.settings.type_proc
                    ts = 'Channel';
                    tsn = num2str(obj.rerp_profile.include_chans(i));
                else
                    ts = 'Component';
                    tsn = num2str(obj.rerp_profile.include_comps(i));
                end
                
                if obj.rerp_profile.settings.hed_enable
                    v = 'Tag';
                else
                    v = 'Event type';
                end
                
                for j=obj.rerp_plot_spec.event_idx
                    % Plot the rERSP estimates
                    this_xaxis_ms = xaxis_ms{j};
                    this_estimate = estimates{j,i};
                    this_tag = tags{j};
                    
                    scrollsubplot(1,1,m,h);
                    a=imagesc(this_xaxis_ms, 1:nbins, this_estimate');
                    colormap('jet');
                    title(['rERSP (db Power), ' ts ': ' tsn ', ' v ': ' this_tag]);
                    set(gca, 'YTick', linspace(1,nbins,10));
                    set(gca, 'YTickLabel', y_axis_hz(linspace(1,nbins,10)));  
                    xlabel('time ( ms )');
                    ylabel('frequency ( Hz ) ');
                    axis xy;
                    %                     num_yticks=20;
                    %                     y_axis_hz_tick = (1:(num_yticks-1))*sr/(num_yticks*2);
                    %                     set(gca, 'ytick', 1:num_yticks);
                    %                     set(gca, 'yticklabel', flipdim(y_axis_hz_tick,2));
                    
                    hcmenu = uicontextmenu;
                    uimenu(hcmenu, 'Label', 'Publish graph', 'Callback', @RerpResult.gui_publish);
                    set(a,'uicontextmenu', hcmenu)
                    
                    m=m+1;
                end
            end
            
        end
        
        function plotGridSearch(obj, h)
            %Plot predictive surfaces and optimal values from regularization grid search
            %   Usage:
            %       rerp_result.plotGridSearch(rerp_plot_spec.ts_idx);
            %       rerp_result.plotGridSearch(rerp_plot_spec.ts_idx, h);
            import rerp_dependencies.*
            
            if ~exist('h', 'var')
                h=figure;
            end
            
            assert(obj.rerp_profile.settings.regularization_enable==1, 'RerpResult: regularization is not enabled for this profile');
            assert(obj.rerp_profile.settings.cross_validate_enable==1, 'RerpResult: cross validation is not enabled for this profile');
            
            this_level = obj.gridsearch;
            m=1;
            if obj.rerp_profile.settings.type_proc
                incl_ts=obj.rerp_profile.include_chans;
                ts = 'Channel';
            else
                incl_ts=obj.rerp_profile.include_comps;
                ts = 'Component';
            end
            
            while 1
                gr = this_level.grid_results;
                
                % Plot surfaces for indicated time series at this level
                for i=obj.rerp_plot_spec.ts_idx
                    this_lambda_range = this_level.lambda_range{i};
                    pred_surf = zeros(1,numel(gr));
                    
                    for j=1:numel(gr)
                        pred_surf(j) = gr{j}.average_total_rsquare(i);
                    end
                    
                    scrollsubplot(1,1,m,h);
                    
                    if strcmp(obj.rerp_profile.settings.penalty_func, 'Elastic net')
                        [lambda_grid_L1, lambda_grid_L2] = meshgrid(this_lambda_range(:,1), this_lambda_range(:, 2));
                        p=mesh(lambda_grid_L1, lambda_grid_L2, reshape(pred_surf, size(lambda_grid_L2)));
                        props=get(p);
                        grid on;
                        colormap('jet');
                        
                        title(['Predictive surface, Elastic net, '  ts ': ' num2str(incl_ts(i)) ', level ' num2str(m)]);
                        xlabel('\lambda  1');
                        ylabel('\lambda  2');
                        zlable('average rsquare');
                        zlim([max(0,min(min(pred_surf))) max(max(pred_surf))]);
                        
                        hold on;
                        opt_lambda1 = this_lambda_range(lambda_grid_L1==max(max(pred_surf)));
                        opt_lambda2 = this_lambda_range(lambda_grid_L2==max(max(pred_surf)));
                        line([opt_lambda1 opt_lambda1], [opt_lambda2 opt_lambda2], [0 max(max(pred_surf))],'Color',props.Color,'LineStyle',props.LineStyle, 'linewidth',1);
                        
                        zpos = (max(max(pred_surf)) + max(0,min(min(pred_surf))))/2;
                        text(opt_lambda1, opt_lambda2, zpos, ['   lambda  =  (' num2str(opt_lambda1) ', ' num2str(opt_lambda2) ')'], 'color',props.Color);
                    end
                    
                    if strcmp(obj.rerp_profile.settings.penalty_func, 'L1 norm')
                        p=plot(this_lambda_range, pred_surf');
                        props=get(p);
                        grid on;
                        title(['Predictive surface, L1 norm penalty, ' ts ': ' num2str(incl_ts(i)) ', level ' num2str(m)]);
                        ylim([max(0, min(pred_surf)), max(pred_surf)]);
                        xlabel('\lambda');
                        ylabel('R ^2');
                        
                        hold on;
                        opt_lambda = this_lambda_range(pred_surf==max(pred_surf));
                        line([opt_lambda opt_lambda],[0 max(pred_surf)*1.1],'Color',props.Color,'LineStyle',props.LineStyle, 'linewidth',1);
                        
                        ypos = (max(pred_surf) + max(0, min(pred_surf)))/2;
                        text(opt_lambda, ypos, ['   lambda  =  ' num2str(opt_lambda)], 'Color',props.Color);
                    end
                    
                    if strcmp(obj.rerp_profile.settings.penalty_func, 'L2 norm')
                        p=plot(this_lambda_range, pred_surf');
                        props=get(p);
                        grid on;
                        title(['Predictive surface, L2 norm penalty, ' ts ': ' num2str(incl_ts(i)) ', level ' num2str(m)]);
                        ylim([max(0, min(pred_surf)), max(pred_surf)]);
                        xlabel('\lambda');
                        ylabel('R ^2');
                        hold on;
                        opt_lambda = this_lambda_range(pred_surf==max(pred_surf));
                        line([opt_lambda opt_lambda],[0 max(pred_surf)],'Color',props.Color,'LineStyle',props.LineStyle, 'linewidth',1);
                        
                        ypos = (max(pred_surf) + max(0, min(pred_surf)))/2;
                        text(opt_lambda, ypos, ['   lambda  =  ' num2str(opt_lambda)],'Color',props.Color);
                    end
                    
                    
                    hcmenu = uicontextmenu;
                    uimenu(hcmenu, 'Label', 'Publish graph', 'Callback', @RerpResult.gui_publish);
                    set(gca,'uicontextmenu', hcmenu)
                    
                end
                
                % Proceed to next level in grid search
                try
                    this_level=this_level.next_stage.gridsearch;
                    m=m+1;
                catch
                    break;
                end
            end
        end
        
        function setPlotTimeSeries(obj, time_series_nums)
            %Set RerpResult to plot specific channels or components
            %   Usage:
            %       rerp_result.setPlotTimeSeries([17 2]);
            if ischar(time_series_nums)
                time_series_nums = str2double(time_series_nums);
            end
            
            if obj.rerp_profile.settings.type_proc
                [~, idx] = intersect(obj.rerp_profile.include_chans, time_series_nums);
                [~, not_idx] = setdiff(time_series_nums, obj.rerp_profile.include_chans);
                if ~isempty(not_idx)
                    disp('RerpResult: the following channel numbers were not found:');
                    disp(time_series_nums(not_idx));
                end
            else
                [~, idx] = intersect(obj.rerp_profile.include_comps, time_series_nums);
                [~, not_idx] = setdiff(time_series_nums, obj.rerp_profile.include_comps);
                if ~isempty(not_idx)
                    disp('RerpResult: the following component numbers were not found:');
                    disp(time_series_nums(not_idx));
                end
            end
            
            obj.rerp_plot_spec.ts_idx = idx;
        end
        
        function setPlotEventTypes(obj, event_types)
            %Set RerpResult to plot specific event types or HED tags
            %   Usage:
            %       rerp_result.setPlotEventTypes({'1' '2' '64'});
            %       rerp_result.setPlotEventTypes({'stimulus/visual/target' 'stimulus/expected'});
            if obj.rerp_profile.settings.hed_enable
                tags = obj.get_plotting_params;
                [~, idx] = intersect(tags, event_types);
                [~, not_idx] = setdiff(event_types, tags);
                if ~isempty(not_idx)
                    disp('RerpResult: the following HED tags were not found:');
                    disp(event_types(not_idx));
                end
            else
                [~, idx] = intersect(obj.rerp_profile.include_event_types, event_types);
                [~, not_idx] = setdiff(event_types, obj.rerp_profile.include_event_types);
                if ~isempty(not_idx)
                    disp('RerpResult: the following event types were not found:');
                    disp(event_types(not_idx));
                end
            end
            
            obj.rerp_plot_spec.event_idx = idx;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Utility
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [tags, estimates, xaxis_ms, epoch_boundaries] = get_plotting_params(obj)
            %Return a representation of the rerp estimate that is convenient for plotting
            %   Usage:
            %       [tags, estimates, xaxis_ms, epoch_boundaries] = rerp_result.get_plotting_params;
            %           tags: cell array of either hed tags or event types
            %           estimates: cell array of the ERP waveforms associated
            %               with each element in tags
            %           xaxis_ms: x axis for plotting each element in estimates
            %           epoch_boundaries: the epoch boundaries used for each element in tag
            import rerp_dependencies.*
            
            p=obj.rerp_profile;
            
            if obj.ersp_flag
                nbins = obj.rerp_profile.settings.nbins;
                raw_estimate = reshape(obj.rerp_estimate, [size(obj.rerp_estimate,1), size(obj.rerp_estimate,2)/nbins, nbins]);
            else
                raw_estimate = obj.rerp_estimate;
            end
            
            tags = obj.tags;
            idx = obj.parameter_idx_layout;
            stripped_tags=RerpTagList.strip_brackets(tags);
            estimates=cell(length(tags), size(raw_estimate, 2));
            xaxis_ms=cell(size(tags));
            epoch_boundaries=cell(size(tags));
            
            for i=1:length(stripped_tags)
                this_tag=stripped_tags{i};
                if ~isempty(intersect(this_tag, p.settings.continuous_tag))
                    epoch_boundaries{i}=p.settings.continuous_epoch_boundaries;
                else
                    epoch_boundaries{i}=p.settings.category_epoch_boundaries;
                end
                
                xaxis_ms{i} = RerpResult.get_xaxis_ms(epoch_boundaries{i}, p.sample_rate);
                for j=1:size(raw_estimate, 2)
                    estimates{i,j} = squeeze(raw_estimate(idx{i}, j, :));
                end
            end
        end
        
        function setLastResult(obj)
            %Save as results/last.rerp_result
            %   Usage:
            %       rerp_result.setLastResult;
            obj.saveRerpResult('path', fullfile(RerpProfile.rerp_path, 'results','last.rerp_result'));
        end
        
        function saveRerpResult(obj, varargin)
            %Save a result to disk
            %   Usage:
            %       rerp_result.saveRerpResult;
            %           opens a gui to choose the path to save result
            %
            %       rerp_result.saveRerpResult('rerp_path', '/data/projects/RSVP');
            %           opens gui starting at that path
            %
            %       rerp_result.saveRerpResult('path', '/data/projects/RSVP/exp_53.rerp_result');
            %           save this result to the specific path (will create the
            %           dir if does not exist)
            import rerp_dependencies.*
            
            p=inputParser;
            addOptional(p,'path',[]);
            addOptional(p,'rerp_path', fullfile(RerpProfile.rerp_path, 'results'));
            
            parse(p, varargin{:});
            temp = regexp(obj.rerp_profile.eeglab_dataset_name, '.*[\\\/](.*)\.set', 'tokens');
            
            if ~isempty(temp)
                fn = temp{1}{1};
            else
                fn='';
            end
            
            path = p.Results.path;
            rerp_path=p.Results.rerp_path;
            if isempty(path)
                %No path specified, launch GUI
                [filename, pathname] = uiputfile('*.rerp_result', 'Save rerp result as:', fullfile(rerp_path, fn));
                path = [pathname filename];
                
            else
                path2file = regexp(path, '(.*)[\\\/].*$','tokens');
                path2file = path2file{1}{1};
                if isempty(dir(path2file))
                    mkdir(path2file);
                end
                
                filename=1;
            end
            
            path2file = regexp(path, '(.*)(?:\.rerp_result)','tokens');
            path2file = path2file{1}{1};
            
            %Save result to disk
            if filename
                try
                    save([path2file '.rerp_result'], 'obj','-mat');
                    fprintf('RerpResult: saved result to disk, %s\n', path);
                catch e
                    fprintf('RerpResult: could not save the specified result to disk %s\n', path);
                    rethrow(e);
                end
            end
        end
        
        function modeled_data = getDataModel(obj)
            %Synthesize model esitmate of the original continuous data
            %   Usage: modeled_data = rerp_result.getDataModel;
            import rerp_dependencies.*
            if ~isempty(obj.rerp_estimate)
                [predictor, data_pad] = obj.rerp_profile.predictor;
                modeled_data = predictor*obj.rerp_estimate;
                modeled_data = modeled_data((data_pad(1)+1):(end-data_pad(2)),:);
            else
                modeled_data=[];
                disp('RerpResult: no anaysis results present, can not synthesize modeled data');
            end
        end
        
        function rsq = setSortIdx(obj)
            for i=1:length(obj)
                if obj(i).ersp_flag
                    nbins=obj(i).rerp_profile.settings.nbins;
                    rsq = max(reshape(obj(i).average_total_rsquare, [nbins, length(obj(i).average_total_rsquare)/nbins]));
                else
                    rsq = obj(i).average_total_rsquare;
                end
                
                %This sets the time series order to be sorted by rsquare,
                %if that is turned on
                if obj(i).rerp_plot_spec.sort_by_r2
                    [~, obj(i).rerp_plot_spec.sort_idx] = sort(rsq, 'descend');
                else
                    obj(i).rerp_plot_spec.sort_idx = 1:length(rsq);
                end
                
            end
        end
    end
    
    methods (Static=true)
        
        function rerp_result = loadRerpResult(varargin)
            %Load a RerpResult from disk
            %   Usage:
            %       rerp_result = RerpResult.loadRerpResult;
            %           Select .rerp_result file using GUI
            %
            %       rerp_result = RerpResult.loadRerpResult('rerp_path', '/data/projects/RSVP');
            %           Open GUI starting at that path
            %
            %       rerp_result = RerpResult.loadRerpResult('path', '/data/projects/RSVP/exp_53.rerp_result');
            %           Load rerp_result from that path, if present
            import rerp_dependencies.*
            
            rerp_result=[];
            p=inputParser;
            addOptional(p,'path',[]);
            addOptional(p,'rerp_path', fullfile(RerpProfile.rerp_path,'results'));
            
            parse(p, varargin{:});
            
            if isempty(p.Results.path)
                %No path specified, launch GUI
                [filename, pathname] = uigetfile({'*.rerp_result'}, 'Load rerp result:', p.Results.rerp_path, 'multiselect','on');
                if iscell(filename)
                    path = cellfun(@(x) [pathname x], filename, 'UniformOutput', false);
                elseif filename
                    path = cellstr(fullfile(pathname, filename));
                    filename=cellstr(filename);
                else
                    return;
                end
            else
                pathname = p.Results.path;
                %If path was a directoy, we load all the .rerp_result files
                %in that directory
                if isdir(pathname)
                    profdir=dir(fullfile(pathname, '*.rerp_result'));
                    filename={profdir(:).name};
                    path=cellfun(@(x) fullfile(pathname, x), filename, 'uniformoutput', false);
                else
                    filename=regexp(pathname, '.*[\\\/](.*\.rerp_result)','match');
                    path={pathname};
                end
            end
            
            %Read result from disk
            rerp_result=cell(1,length(path));
            for i=1:length(path)
                if path{i}
                    try
                        res = load(path{i}, '-mat');
                        rerp_result{i} = res.obj;
                        rerp_result{i}.name = filename{i};
                        
                    catch
                        fprintf('RerpResult: could not read the specified result, %s\n', path{i});
                    end
                end
            end
            
            rerp_result = rerp_result(cell2mat(cellfun(@(x) ~isempty(x), rerp_result, 'UniformOutput', false)));
            rerp_result = [rerp_result{:}];
        end
        
        function result = combineRerpResults(rerp_results)
            %Parallel process individual time courses for regression, then combine those results back into a single object
            %rerp_results parameter is a cell array of RerpResult objects, from
            %the same dataset, in increasing channel/IC number order.
            %   Usage:
            %       combined_rerp_result = RerpResult.combineRerpResults(rerp_results);
            %           Combine a cell array of RerpResult objects into a
            %           single object (parallelize computation on channels).
            rerp_results=rerp_results(~cellfun('isempty',rerp_results));
            
            result={};
            if ~isempty(rerp_results)
                result =  copy(rerp_results{1});
                
                m=length(result.lambda)+1;
                for i=2:length(rerp_results)
                    this_result = rerp_results{i};
                    nchans = size(this_result.rerp_estimate,2);
                    idx = m:(m+nchans-1);
                    
                    if ~isempty(result.rerp_profile)
                        if result.rerp_profile.settings.type_proc
                            result.rerp_profile.include_chans(idx) = this_result.rerp_profile.include_chans;
                        else
                            result.rerp_profile.include_comps(idx) = this_result.rerp_profile.include_comps;
                        end
                    end
                    
                    result.compute_time_seconds = result.compute_time_seconds + this_result.compute_time_seconds;
                    
                    try
                        result.lambda(idx) = this_result.lambda;
                    catch
                    end
                    
                    result.rerp_estimate(:, idx) = this_result.rerp_estimate;
                    
                    try
                        result.admm_residual(:, idx) = this_result.admm_residual;
                    catch
                    end
                    
                    if ~isempty(result.average_total_rsquare);
                        result.average_total_rsquare(idx) = this_result.average_total_rsquare;
                    end
                    
                    if ~isempty(result.average_event_rsquare);
                        result.average_event_rsquare(:,idx) = this_result.average_event_rsquare;
                    end
                    
                    for j=1:length(result.total_xval_folds)
                        result.total_xval_folds(j).noise_variance(idx)=this_result.total_xval_folds(j).noise_variance;
                        result.total_xval_folds(j).data_variance(idx)=this_result.total_xval_folds(j).data_variance;
                        result.total_xval_folds(j).num_samples(idx)=this_result.total_xval_folds(j).num_samples;
                    end
                    
                    for j=1:length(result.event_xval_folds)
                        result.event_xval_folds(j).noise_variance(:,idx)=this_result.event_xval_folds(j).noise_variance;
                        result.event_xval_folds(j).data_variance(:,idx)=this_result.event_xval_folds(j).data_variance;
                        result.event_xval_folds(j).num_samples(:,idx)=this_result.event_xval_folds(j).num_samples;
                    end
                    
                    result.gridsearch = RerpResult.mergeGridSearch(result, this_result, idx);
                    
                    m=m+nchans;
                end
            end
        end
        
        function xaxis_ms = get_xaxis_ms(epoch_boundaries, sample_rate)
            %Get x axis ticks for plotting erp waveform
            %   Usage: xaxis_ms = RerpResult.get_xaxis_ms([-1 2], 256)
            epoch_length = epoch_boundaries(2)-epoch_boundaries(1);
            ns=ceil(sample_rate*epoch_length);
            xaxis_ms=ceil(((0:(ns-1))+epoch_boundaries(1)*sample_rate)*1000.0/sample_rate);
        end
    end
    
    methods (Hidden=true)
        function rsquare_significance = get_total_rsquare_significance(obj)
            % Returns time series numbers where rsquare was statistically
            % different from zero mean at p<rerp_plot_spec.significance_level
            import rerp_dependencies.ttestw
            
            data_variance=zeros(length(obj.total_xval_folds), length(obj.total_xval_folds(1).data_variance));
            noise_variance=zeros(length(obj.total_xval_folds), length(obj.total_xval_folds(1).noise_variance));
            weight=zeros(length(obj.total_xval_folds), length(obj.total_xval_folds(1).num_samples));
            
            for i=1:length(obj.total_xval_folds)
                data_variance(i,:) = obj.total_xval_folds(i).data_variance;
                noise_variance(i,:) = obj.total_xval_folds(i).noise_variance;
                weight(i,:) = obj.total_xval_folds(i).num_samples;
            end
            
            rsquare = 1 - noise_variance./data_variance;
            rsquare_significance=zeros([1,size(rsquare,2)]);
            for i=1:size(rsquare,2)
                rsquare_significance(i) = squeeze(ttestw(rsquare(:,i), 0, weight(:,i), 'Alpha', obj.rerp_plot_spec.significance_level));
            end
        end
        
        function rsquare_significance = get_event_rsquare_significance(obj)
            % Returns time series numbers where rsquare was statistically
            % different from zero mean at p<rerp_plot_spec.significance_level.
            import rerp_dependencies.ttestw
            
            data_variance=zeros([length(obj.event_xval_folds) size(obj.event_xval_folds(1).data_variance)]);
            noise_variance=zeros([length(obj.event_xval_folds) size(obj.event_xval_folds(1).noise_variance)]);
            weight=zeros([length(obj.event_xval_folds) size(obj.event_xval_folds(1).num_samples)]);
            for i=1:length(obj.event_xval_folds)
                data_variance(i,:,:) = obj.event_xval_folds(i).data_variance;
                noise_variance(i,:,:) = obj.event_xval_folds(i).noise_variance;
                weight(i,:,:) = obj.event_xval_folds(i).num_samples;
            end
            
            rsquare = 1 - noise_variance./data_variance;
            rsquare_significance=zeros([size(rsquare,2), size(rsquare,3)]);
            for i=1:size(rsquare,2)
                for j=1:size(rsquare,3)
                    rsquare_significance(i, j) = squeeze(ttestw(rsquare(:,i, j), 0, weight(:,i, j), 'Alpha', obj.rerp_plot_spec.significance_level));
                end
            end
            
            rsquare_significance = reshape(rsquare_significance, size(obj.event_xval_folds(1).data_variance));
            rsquare_significance(isnan(rsquare_significance))=0;
        end
        
        function delay = get_delay_times(obj, event_nums, delay_var, num_samples)
            import rerp_dependencies.*
            
            regexp_str_in_parentheses = '.*\((.*)\).*';
            regexp_str_out_parentheses = '(.*)(?:\s*\(.*\))?';
            
            delay_context_tag = strtrim(regexp(delay_var, regexp_str_in_parentheses, 'tokens'));
            delay_locking_tag = strtrim(regexp(delay_var, regexp_str_out_parentheses, 'tokens'));
            
            delay=zeros(1,length(event_nums));
            events = obj.rerp_profile.these_events;
            
            m=1;
            if obj.rerp_profile.settings.hed_enable
                % Now find the delay for each delay_var event
                for i=1:length(event_nums)
                    this_evt_num=event_nums(i);
                    this_latency = events.latencyInFrame(this_evt_num);
                    if ~isempty(delay_var)
                        for j=this_evt_num:length(events.label)
                            
                            if events.latencyInFrame(j) > (this_latency + num_samples)
                                delay(m) = num_samples;
                                break;
                            end
                            
                            these_hed_tags = hedTree.hed_tag_count(regexp(events.hedTag{j}, '[,;]','split'), 0, 0);
                            if ~isempty(intersect(these_hed_tags, delay_locking_tag{1}{1}))
                                if isempty(delay_context_tag)
                                    delay_latency=events.latencyInFrame(j);
                                    delay(m) = delay_latency-this_latency;
                                    break;
                                else
                                    if ~isempty(intersect(these_hed_tags, delay_context_tag{1}{1}))
                                        delay_latency=events.latencyInFrame(j);
                                        delay(m) = delay_latency-this_latency;
                                        break;
                                    end
                                end
                            end
                        end
                    end
                    m=m+1;
                end
                
            else
                for i=event_nums
                    this_latency = events.latencyInFrame(i);
                    for j=i:length(events.label)
                        if events.latencyInFrame > (this_latency + num_samples)
                            delay(m) = num_samples;
                            break;
                        end
                        
                        this_event = events.label{j};
                        if strcmp(delay_var, this_event)
                            delay_latency=events.latencyInFrame(j);
                            delay(m) = delay_latency-this_latency;
                            break;
                        end
                    end
                    m=m+1;
                end
            end
        end
        
        function [rerp_epochs, event_nums] = get_rerp_epochs(obj, data, locking_var, front_pad, boundaries)
            import rerp_dependencies.RerpTagList;
            num_samples = ceil((boundaries(2)-boundaries(1))*obj.rerp_profile.sample_rate); 

            % Extract epochs corresponding to tags or event codes
            events = obj.rerp_profile.these_events;
            rerp_epochs = zeros([num_samples, length(events.label), size(data,2)]);
            m=1;
            
            regexp_str_in_parentheses = '.*\((.*)\).*';
            context_tag = regexp(locking_var, regexp_str_in_parentheses, 'tokens');
            locking_tag = RerpTagList.strip_label({locking_var});
            locking_tag=locking_tag{1};
            if obj.rerp_profile.settings.hed_enable
                % Get event numbers for locking_var (possibly in a context
                % group).
                if ~isempty(context_tag)
                    for i=1:length(obj.rerp_profile.context_group)
                        this_group = obj.rerp_profile.context_group(i);
                        for j=1:length(this_group.children)
                            this_child = this_group.children(i);
                            if strcmp(strtrim(context_tag{1}{1}), this_child.tag)
                                idx = find(strcmp(locking_tag, this_child.included_tag), 1);
                                event_nums = this_child.included_ids{idx};
                            end
                        end
                    end
                else
                    tag_idx = strcmp(locking_tag, obj.rerp_profile.hed_tree.uniqueTag);
                    event_nums = obj.rerp_profile.hed_tree.originalHedStringId{tag_idx};
                end
                
            else
                event_nums=zeros(0,1);
                for i=1:length(events.label)
                    this_event = events.label{i};
                    if strcmp(locking_tag, this_event)
                        this_latency=events.latencyInFrame(i);
                        rerp_epochs(:,m,:) = data(this_latency:(this_latency+num_samples-1),:);
                        event_nums(m) = i;
                        m=m+1;
                    end
                    
                end
            end
            
            m=1;
            %Extract the epochs for the event numbers selected
            for i=event_nums(:)'
                this_latency = round(events.latencyInFrame(i));
                this_latency = this_latency + front_pad + ceil(boundaries(1)*obj.rerp_profile.sample_rate);
                rerp_epochs(:,m,:) = data(this_latency:(this_latency+num_samples-1),:);
                m=m+1;
            end
            rerp_epochs=rerp_epochs(:,1:(m-1),:);
        end
    end
    
    methods (Static=true, Hidden=true)
        function main_gridsearch = mergeGridSearch(main_result, added_result, idx)
            % Recursively combine grid search results from two rerp_results
            main_gridsearch=main_result.gridsearch;
            added_gridsearch=added_result.gridsearch;
            
            if (~isempty(main_gridsearch)) && (~length(fieldnames(main_gridsearch))==0)
                main_gridsearch.lambda_range(idx) = added_gridsearch.lambda_range;
                
                for i=1:length(main_gridsearch.grid_results)
                    main_gridsearch.grid_results{i} = RerpResult.combineRerpResults({main_gridsearch.grid_results{i}, added_gridsearch.grid_results{i}});
                end
                
                try
                    main_gridsearch.next_stage.gridsearch = RerpResult.mergeGridSearch(main_gridsearch.next_stage, added_gridsearch.next_stage, idx);
                catch
                end
            end
        end
        
        function gui_publish(varargin)
            % Callback for axes context menu: publishes the figure
            ax=gca;
            h=figure;
            
            newax=copyobj(ax, h);
            pr=get(ax,'UserData');
            
            try
                lstr = get(pr.legend,'String');
                legend(lstr);
                lprops=get(pr.legend,'UserData');
                try
                    p=copyobj(lprops.plotHandles,newax);
                    legend(p,lstr);
                catch
                    legend(lstr);
                end
            catch
                
            end
            
            set(gca,  'ActivePositionProperty', 'OuterPosition',...
                'OuterPosition', [0 0 1 1]);
            
            rerp_publish_gui('hFig', h);
        end
    end
    
    methods (Access = protected)
        % Override copyElement method:
        function cpObj = copyElement(obj)
            import rerp_dependencies.*
            
            % Make a shallow copy of all properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            
            % Make a deep copy of these objects
            if isa(obj.rerp_profile,'RerpProfile')
                cpObj.rerp_profile = copy(obj.rerp_profile);
            end
            
            if isa(obj.rerp_plot_spec,'RerpPlotSpec')
                cpObj.rerp_plot_spec = copy(obj.rerp_plot_spec);
            end
            
            if isa(obj.total_xval_folds, 'RerpXvalFold')
                cpObj.total_xval_folds = copy(obj.total_xval_folds);
            end
            
            if isa(obj.event_xval_folds,'RerpXvalFold')
                cpObj.event_xval_folds = copy(obj.event_xval_folds);
            end
        end
    end
end
