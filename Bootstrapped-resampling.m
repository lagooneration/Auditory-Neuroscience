%% Code for bootstrap resampling

 

   % compute RMS across channels for a given condition in one subject.

   % RMS = subject x condition x time

  rms(s,jj,:) = sqrt(nanmean(data.^2,1));

 

 

%% Calculate bootstrap mean and SD for all conditions

   

    B = 1000;

    for ii = 1:numel(cond_list)

       

        % DIM: sub * cond * time

        x = squeeze(rms(:,ii,:));

       

        %  Force x to be time*repetitions

        if size(x,2) > size(x,1), x=x'; end

       

        [mn,sd] = fBootstrapRMS(x,B);

       

        bs_mean(ii,:) = mn';

       

        bs_std(ii,:) = sd';

    end

 

 

    h_fig = figure(2); clf

 

    for ii = to_plot

       a = bs_mean(ii,:)';

       plot(D.time, a);

       hold on;

    end

    legend(strrep(cond_list(to_plot),'_','-'));

 

    for ii = to_plot

        % 2*SE of boostrap resampling

        b = 2*bs_std(ii,:)'/sqrt(numel(subject_list));

        a = bs_mean(ii,:)';

       

        style = 'meanbased';

        if strcmp(style,'zerobased'); % plot bs sd on zero-line

            Y=[b;-flipud(b)]';

        elseif strcmp(style,'meanbased');

            Y=[b+a;flipud(-b+a)]';

        end

       

        abscissa = D.time(:);

        X = [abscissa; flipud(abscissa)];

        h = fill(X,Y,C(1,:),'edgecolor','none','facealpha',1); hold on;

        % plot(abscissa,a*0,'k'); % plot zero line

    end

    set(gca,'children',flipud(get(gca,'children'))); % Send shaded areas in background
