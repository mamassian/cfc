% plot function for confidence forced-choice
%
% 23-MAR-2019 - pascal mamassian

function cfc_plot_02(grouped_data, varargin)

% -> default optional arguments
dflt_human_model        = [];	% human against model choice
dflt_choi_by_resp       = false;    % choice by response
dflt_psychometric       = [];    % psychometric functions

% -> parse all arguments
ip = inputParser;
addRequired(ip, 'grouped_data', @isnumeric);
addParameter(ip, 'human_model', dflt_human_model, @isnumeric);
addParameter(ip, 'choi_by_resp', dflt_choi_by_resp, @islogical);
addParameter(ip, 'psychometric', dflt_psychometric, @isnumeric);
parse(ip, grouped_data, varargin{:});
plot_human_model = ip.Results.human_model;
plot_choi_by_resp = ip.Results.choi_by_resp;
plot_psychometric = ip.Results.psychometric;


stims2 = grouped_data(:, 1:2);
resps2 = grouped_data(:, 3:4);
choic2 = grouped_data(:, 5:6);
tasks2 = grouped_data(:, 7:8);

% -> nb tasks
tsk_vals = unique(tasks2);
nb_unique_tasks = size(tsk_vals, 1);

% -> nb kinds of trials (s1, s2, r1, r2, c2)
nb_knds = size(grouped_data, 1);

% -> nb trials per (s1, s2, r1, r2), duplicated in a vector
total1 = sum(choic2, 2);
total2 = repmat(total1, 1, 2);

% -> fraction of chosen interval 1
fraction_chosen = choic2(:, 1)./total1;

% -> list of unique stimuli across all tasks
stimtask = complex(stims2, tasks2);
[stm_vals, ~, stm_ic2] = unique(stimtask);
nb_unique_stims2 = size(stm_vals, 1);

% -> unique stimulus values for each task
stims_dim = NaN(nb_unique_tasks, nb_unique_stims2);
stims_ind = NaN(nb_unique_tasks, nb_unique_stims2);
stims_nb = zeros(nb_unique_tasks, 1);
for tt = 1:nb_unique_tasks
    tsk_inds = find(imag(stm_vals) == tt);
    nbs = length(tsk_inds);
    [stims_dim(tt, 1:nbs), ivals] = sort(real(stm_vals(tsk_inds)));
    stims_ind(tt, 1:nbs) = tsk_inds(ivals)';
    stims_nb(tt) = nbs;
end


if (plot_psychometric)
    unsorted_prob_lst = NaN(nb_unique_tasks, nb_unique_stims2);
    chosen_rsp_data_lst = NaN(nb_unique_tasks, nb_unique_stims2);
    declined_rsp_data_lst = NaN(nb_unique_tasks, nb_unique_stims2);
    chosen_rsp_model_lst = NaN(nb_unique_tasks, nb_unique_stims2);
    for tt = 1:nb_unique_tasks
        nbs = stims_nb(tt);
        tsk_inds = stims_ind(tt, 1:nbs);

        for ww = 1:nbs
            inds = (stm_ic2 == tsk_inds(ww));      % indices that have the same stim

            nb_resps = sum(total2(inds));       % nb. of responses for this particular stim
            resp_vals = resps2(inds) .* total2(inds);
            rsp_prob = sum(resp_vals) / nb_resps;   % unsorted responses
            unsorted_prob_lst(tt, ww) = rsp_prob;

            nb_resps = sum(choic2(inds));       % nb. of responses for this particular stim
            resp_vals = resps2(inds) .* choic2(inds);
            chosen_rsp_data = sum(resp_vals) / nb_resps;    % chosen responses in data
            chosen_rsp_data_lst(tt, ww) = chosen_rsp_data;

            decli2 = total2 - choic2;
            nb_resps = sum(decli2(inds));       % nb. of responses for this particular stim
            resp_vals = resps2(inds) .* decli2(inds);
            declined_rsp_data = sum(resp_vals) / nb_resps;    % declined responses in data
            declined_rsp_data_lst(tt, ww) = declined_rsp_data;

            nn1_plot_full = [plot_psychometric, 1-plot_psychometric] .* total2;
            nb_resps = sum(nn1_plot_full(inds));       % nb. of responses for this particular stim
            resp_vals = resps2(inds) .* nn1_plot_full(inds);
            chosen_rsp_model = sum(resp_vals) / nb_resps;    % chosen responses in model
            chosen_rsp_model_lst(tt, ww) = chosen_rsp_model;
        end

    end

    % -> plot psychometric function
    fig = figure;
    set(fig, 'Position', [65 600 720 340]);
    for tt = 1:nb_unique_tasks
        subplot(1, nb_unique_tasks, tt);
        hold on;
        set(gca, 'FontName', 'Arial'); set(gca, 'FontSize', 24);
        line([0 0], [0 1], 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
        line([-5 5], [0.5 0.5], 'Color', [0.5 0.5 0.5], 'LineWidth', 2);

        nbs = stims_nb(tt);
        xvals = stims_dim(tt, 1:nbs);
        unsorted = unsorted_prob_lst(tt, 1:stims_nb(tt));
        chosen = chosen_rsp_data_lst(tt, 1:stims_nb(tt));
        declined = declined_rsp_data_lst(tt, 1:stims_nb(tt));
        model = chosen_rsp_model_lst(tt, 1:stims_nb(tt));

        col = [0.5, 0.5, 0.5];
        plot(xvals, unsorted, 'o', 'MarkerSize', 16, 'MarkerFaceColor', col, 'Color', col);

        col = [0, 0.8, 0.5];
        plot(xvals, chosen, 'o', 'MarkerSize', 16, 'MarkerFaceColor', col, 'Color', col);
        plot(xvals, model, '-', 'Color', col, 'LineWidth', 3);
                
        col = [0.8, 0, 0.5];
        plot(xvals, declined, 'o--', 'MarkerSize', 16, 'MarkerFaceColor', col, 'Color', col, 'LineWidth', 3);

        ylim([0, 1]);
        xlabel('Signal Strength');
        ylabel('Proportion ''Right''');
    end
end


if (plot_human_model)
    fig = figure;
    set(fig, 'Position', [785 600 430 340]);
    hold on;
    set(gca, 'FontName', 'Arial'); set(gca, 'FontSize', 24);
    line([0 1], [0 1], 'Color', [0.5 0.5 0.5], 'LineWidth', 2);

    sizes = 100 * total1 ./ max(total1);
    col = [0.8, 0.0, 0.5];
    scatter(plot_human_model, fraction_chosen, sizes, col, 'filled', 'LineWidth', 2);

    axis([0, 1, 0, 1]);
    axis('square');
    xlabel('Model Choice');
    ylabel('Human Choice');

end


if (plot_choi_by_resp)
    choi11 = NaN(stims_nb(1), stims_nb(1), 4);
    count11 = zeros(stims_nb(1), stims_nb(1), 2);
    if (nb_unique_tasks > 1)
        choi22 = NaN(stims_nb(2), stims_nb(2), 4);
        choi12 = NaN(stims_nb(2), stims_nb(1), 4);
        count12 = zeros(stims_nb(2), stims_nb(1), 2);
    end

    % -> plot matrix of choices for pairs of stimulus intensities
    for kk = 1:nb_knds
        stim_pair = stims2(kk, :);
        resp_pair = resps2(kk, :);
        choi_pair = choic2(kk, :);
        task_pair = tasks2(kk, :);

        % -> convert vector {00, 01, 10, 11} into index {1, 2, 3, 4}
        resp_knd = bin2dec(char(resp_pair + '0')) + 1;

        choice_prob = choi_pair(1) / sum(choi_pair);

        if (task_pair(1) == 1)
            stim1_ind = find(stim_pair(1) == stims_dim(1,:));
            if (task_pair(2) == 1)
                stim2_ind = find(stim_pair(2) == stims_dim(1,:));
                choi11(stim2_ind, stim1_ind, resp_knd) = choice_prob;
                count11(stim2_ind, stim1_ind, 1) = count11(stim2_ind, stim1_ind, 1) + choi_pair(1);
                count11(stim2_ind, stim1_ind, 2) = count11(stim2_ind, stim1_ind, 2) + choi_pair(2);
            else
                stim2_ind = find(stim_pair(2) == stims_dim(2,:));
                choi12(stim2_ind, stim1_ind, resp_knd) = choice_prob;
                count12(stim2_ind, stim1_ind, 1) = count12(stim2_ind, stim1_ind, 1) + choi_pair(1);
                count12(stim2_ind, stim1_ind, 2) = count12(stim2_ind, stim1_ind, 2) + choi_pair(2);
            end
        else
            stim1_ind = find(stim_pair(1) == stims_dim(2,:));
            if (task_pair(2) == 1)
                stim2_ind = find(stim_pair(2) == stims_dim(1,:));
                resp_knd = bin2dec(char(fliplr(resp_pair) + '0')) + 1;
                choi12(stim1_ind, stim2_ind, resp_knd) = choice_prob;
                count12(stim1_ind, stim2_ind, 1) = count12(stim1_ind, stim2_ind, 1) + choi_pair(2);
                count12(stim1_ind, stim2_ind, 2) = count12(stim1_ind, stim2_ind, 2) + choi_pair(1);
            else
                stim2_ind = find(stim_pair(2) == stims_dim(2,:));
                choi22(stim2_ind, stim1_ind, resp_knd) = choice_prob;
            end
        end

    end

    fig = figure;
    set(fig, 'Position', [65 45 720 480]);
    knd2ind = [4, 1, 5, 2];     % convert response kind to index for subplot
    for resp_knd = 1:4
        plot_ind = knd2ind(resp_knd);
        subplot(2, 3, plot_ind);

        % -> assume for now that either there is a single task
        % -> or the two tasks are always present in a confidence pair
        if (nb_unique_tasks == 1)
            choi_im = choi11(:, :, resp_knd);
        else
            choi_im = choi12(:, :, resp_knd);
        end

        im = imagesc(flipud(choi_im), [0, 1]);
        set(im, 'AlphaData', ~isnan(flipud(choi_im)));
        axis('square');
        if (nb_unique_tasks == 1)
            xlabel('Interval 1 Stimulus');
            ylabel('Interval 2 Stimulus');   
        else
            xlabel('Task 1 Stimulus');
            ylabel('Task 2 Stimulus');   
        end
        set(gca, 'FontName', 'Arial'); set(gca, 'FontSize', 12);

        uu = 1:stims_nb(1);
        set(gca,'XTick', uu);
        set(gca,'XTickLabel', round(stims_dim(1,uu),2));
        vv = 1:stims_nb(nb_unique_tasks);
        set(gca,'YTick', vv);
        set(gca,'YTickLabel', fliplr(round(stims_dim(nb_unique_tasks,vv),2)));

        title(sprintf('Responses: %s', dec2bin(resp_knd-1,2)));
    end
    subplot(2, 3, 3);
    if (nb_unique_tasks == 1)
        choi_im = count11(:, :, 1) ./ sum(count11, 3);
    else
        choi_im = count12(:, :, 1) ./ (count12(:, :, 1) + count12(:, :, 2));
    end
    im = imagesc(flipud(choi_im), [0, 1]);
    set(im, 'AlphaData', ~isnan(flipud(choi_im)));
    axis('square');
    if (nb_unique_tasks == 1)
        xlabel('Interval 1 Stimulus');
        ylabel('Interval 2 Stimulus');   
    else
        xlabel('Task 1 Stimulus');
        ylabel('Task 2 Stimulus');   
    end
    set(gca, 'FontName', 'Arial'); set(gca, 'FontSize', 12);
    uu = 1:stims_nb(1);
    set(gca,'XTick', uu);
    set(gca,'XTickLabel', round(stims_dim(1,uu),2));
    vv = 1:stims_nb(nb_unique_tasks);
    set(gca,'YTick', vv);
    set(gca,'YTickLabel', fliplr(round(stims_dim(nb_unique_tasks,vv),2)));
    title('Any Response');

    colormap('cool');
    sbpl = subplot(2, 3, 6);
    sbpl.Position = sbpl.Position - [0.1, 0, 0.0, 0];
    axis('off')
    hcb=colorbar;
    set(gca, 'FontName', 'Arial'); set(gca, 'FontSize', 16);
    ylabel(hcb, 'Choice Probability', 'FontName', 'Arial', 'FontSize', 28);

end

end