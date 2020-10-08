% CONFIDENCE FORCED-CHOICE TOOLBOX  v0.3
%
% cfc_plot
%   This is a collection of plot functions for confidence forced-choice.
%
% INPUT:
%   'grouped_data': grouped data per (s1, s2, r1, r2):
%       1st col: stimulus intensity interval 1
%       2nd col: stimulus intensity interval 2
%       3rd col: perceptual decision interval 1
%       4th col: perceptual decision interval 2
%       5th col: nb of confidence choices for interval 1
%       6th col: nb of confidence choices for interval 2
%       7th col: stimulus task for interval 1
%       8th col: stimulus task for interval 2
%
% POSSIBLE PARAMETERS:
%   'human_model': plot fraction chosen for human against a model provided
%                  as argument
%
%   'choi_by_resp': plot matrix of choices for pairs of stimulus strengths
%
%   'zscore': plot z-scores instead of probabilities
%
%   'type1_psychometric': plot psychometric functions for percepts
%       if there is an argument 'arg', plot model data in 'arg'
%           format (n1, n0) for each 'knd'
%
%   'type1_psychometric_all': plot psychometric functions for percepts
%           for each interval and task
%
%   'type2_psychometric': plot psychometric functions for confidence choices
%       if there is an argument 'arg', plot model data in 'arg'
%           format 'choice for intrvl 1' for each (s1, s2, perc, tsk)
%
%   'type2_residuals': plot confidence choice residuals
%
%
% OUTPUT:
%   'plot_data' = struct
%        (vector of struct if number of unique tasks > 1)
%     plot_data(task_no).task               task number
%     plot_data(task_no).sensory_strength	list of stimuli
%     plot_data(task_no).unsorted_prob      Type1 responses (unsorted)
%     plot_data(task_no).chosen_prob        Type1 responses for confidence chosen
%     plot_data(task_no).declined_prob      Type1 responses for confidence declined
%     plot_data(task_no).unsorted_count     nb. of Type1 responses (unsorted)
%     plot_data(task_no).chosen_count       nb. of Type1 responses for confidence chosen
%     plot_data(task_no).declined_count     nb. of Type1 responses for confidence declined
%
%
%   'plot_data_intrvl' = struct
%        similar to 'plot_data', except that data are provided for the two intervals
%     plot_data_intrvl(task_ind).tasks      pair of task numbers for intervals 1 and 2
%     plot_data_intrvl(task_ind).unsorted_prob_intrvl1
%     plot_data_intrvl(task_ind).unsorted_prob_intrvl2
%     plot_data_intrvl(task_ind).chosen_prob_intrvl1
%     plot_data_intrvl(task_ind).chosen_prob_intrvl2
%     plot_data_intrvl(task_ind).declined_prob_intrvl1
%     plot_data_intrvl(task_ind).declined_prob_intrvl2
%     plot_data_intrvl(task_ind).unsorted_count_intrvl1
%     plot_data_intrvl(task_ind).unsorted_count_intrvl2
%     plot_data_intrvl(task_ind).chosen_count_intrvl1
%     plot_data_intrvl(task_ind).chosen_count_intrvl2
%     plot_data_intrvl(task_ind).declined_count_intrvl1
%     plot_data_intrvl(task_ind).declined_count_intrvl2
% 
%
%
% 21-SEP-2020 - pascal mamassian
 
% 29-FEB-2020 - pascal mamassian
% 29-JUN-2020 - generalize to nb. tasks > 2
% 21-SEP-2020 - cleaned up


function [plot_data, plot_data_intrvl] = cfc_plot(grouped_data, varargin)

    % -> default optional arguments
    dflt_human_model         = [];	   % human against model choice
    dflt_choi_by_resp        = false;  % choice by response
    dflt_zscore              = false;  % plot z-scores instead of probabilities
    dflt_t1_psychometric     = false;     % type-1 psychometric functions
    dflt_t1_psychometric_all = [];     % type-1 psychometric functions for all 8 cases
    dflt_t2_psychometric     = [];     % type-2 psychometric functions
    dflt_t2_residuals        = [];     % type-2 residuals (re. super-ideal)

    % -> parse all arguments
    ip = inputParser;
    addRequired(ip, 'grouped_data', @isnumeric);
    addParameter(ip, 'human_model', dflt_human_model, @isnumeric);
    addParameter(ip, 'choi_by_resp', dflt_choi_by_resp, @islogical);
    addParameter(ip, 'zscore', dflt_zscore, @islogical);
    addParameter(ip, 'type1_psychometric', dflt_t1_psychometric);
    addParameter(ip, 'type1_psychometric_all', dflt_t1_psychometric_all);
    addParameter(ip, 'type2_psychometric', dflt_t2_psychometric);
    addParameter(ip, 'type2_residuals', dflt_t2_residuals, @isnumeric);

    parse(ip, grouped_data, varargin{:});
    plot_human_model = ip.Results.human_model;
    plot_choi_by_resp = ip.Results.choi_by_resp;
    plot_zscore = ip.Results.zscore;
    plot_t1_psychometric = ip.Results.type1_psychometric;
    plot_t1_psychometric_all = ip.Results.type1_psychometric_all;
    plot_t2_psychometric = ip.Results.type2_psychometric;
    plot_t2_residuals = ip.Results.type2_residuals;



    % -> prepare output variables
    plot_data = struct;
    plot_data_intrvl = struct;

    % -> un-wrap data matrix
    stims2 = grouped_data(:, 1:2);
    resps2 = grouped_data(:, 3:4);
    choic2 = grouped_data(:, 5:6);
    tasks2 = grouped_data(:, 7:8);

    % -> check how many tasks were run
    tsk_vals = unique(tasks2);
    nb_unique_tasks = size(tsk_vals, 1);

    % -> how many different pairs of tasks
    tsk_pairs_vals = unique(tasks2, 'rows');
    nb_unique_tsk_pairs = size(tsk_pairs_vals, 1);

    % -> nb kinds of trials (s1, s2, r1, r2, c2)
    nb_knds = size(grouped_data, 1);

    % -> nb trials per (s1, s2, r1, r2)
    total1 = sum(choic2, 2);
    total2 = repmat(total1, 1, 2);  % duplicate in a (1,2) vector

    % -> fraction of chosen interval 1 for each trial kind
    fraction_chosen = choic2(:, 1)./total1;

    % -> list of unique stimuli across all tasks
    stimtask = complex(stims2, tasks2);
    [stm_vals, ~, stm_ic2] = unique(stimtask);
    nb_unique_stims2 = size(stm_vals, 1);

    % -> unique stimulus values for each task
    stims_dim = NaN(nb_unique_tasks, nb_unique_stims2);
    stims_ind = NaN(nb_unique_tasks, nb_unique_stims2);
    stims_nb = zeros(nb_unique_tasks, 1);
    for task_no = 1:nb_unique_tasks
        tsk_inds = find(imag(stm_vals) == task_no);
        nbs = length(tsk_inds);
        [stims_dim(task_no, 1:nbs), ivals] = sort(real(stm_vals(tsk_inds)));
        stims_ind(task_no, 1:nbs) = tsk_inds(ivals)';
        stims_nb(task_no) = nbs;
    end
    stims_nb_max = max(stims_nb);


    % -> internal function variables for 'plot_t1_psychometric'
    unsorted_prob_lst     = NaN(nb_unique_tasks, nb_unique_stims2);
    chosen_rsp_data_lst   = NaN(nb_unique_tasks, nb_unique_stims2);
    declined_rsp_data_lst = NaN(nb_unique_tasks, nb_unique_stims2);
    unsorted_count_lst    = NaN(nb_unique_tasks, nb_unique_stims2);
    chosen_count_lst      = NaN(nb_unique_tasks, nb_unique_stims2);
    declined_count_lst    = NaN(nb_unique_tasks, nb_unique_stims2);


    for task_no = 1:nb_unique_tasks
        nbs = stims_nb(task_no);
        tsk_inds = stims_ind(task_no, 1:nbs);

        for ww = 1:nbs
            inds = (stm_ic2 == tsk_inds(ww));      % indices that have the same stim

            nb_resps = sum(total2(inds));       % nb. of responses for this particular stim
            resp_vals = resps2(inds) .* total2(inds);
            rsp_prob = sum(resp_vals) / nb_resps;   % unsorted responses
            unsorted_prob_lst(task_no, ww) = rsp_prob;
            unsorted_count_lst(task_no, ww) = nb_resps;

            nb_resps = sum(choic2(inds));       % nb. of responses for this particular stim
            resp_vals = resps2(inds) .* choic2(inds);
            chosen_rsp_data = sum(resp_vals) / nb_resps;    % chosen responses in data
            chosen_rsp_data_lst(task_no, ww) = chosen_rsp_data;
            chosen_count_lst(task_no, ww) = nb_resps;

            decli2 = total2 - choic2;
            nb_resps = sum(decli2(inds));       % nb. of responses for this particular stim
            resp_vals = resps2(inds) .* decli2(inds);
            declined_rsp_data = sum(resp_vals) / nb_resps;    % declined responses in data
            declined_rsp_data_lst(task_no, ww) = declined_rsp_data;
            declined_count_lst(task_no, ww) = nb_resps;
        end

        plot_data(task_no).task = tsk_vals(task_no);
        plot_data(task_no).sensory_strength = stims_dim(task_no, 1:nbs);
        plot_data(task_no).unsorted_prob = unsorted_prob_lst(task_no, 1:nbs);
        plot_data(task_no).chosen_prob = chosen_rsp_data_lst(task_no, 1:nbs);
        plot_data(task_no).declined_prob = declined_rsp_data_lst(task_no, 1:nbs);
        plot_data(task_no).unsorted_count = unsorted_count_lst(task_no, 1:nbs);
        plot_data(task_no).chosen_count = chosen_count_lst(task_no, 1:nbs);
        plot_data(task_no).declined_count = declined_count_lst(task_no, 1:nbs);

    end



    % -> variables for Type 2 in matrix (s1, s2, perc_pairs, tsk_pairs)
    human_choices = NaN(stims_nb_max, stims_nb_max, 4, nb_unique_tsk_pairs);    % prob choose intrvl 1
    intrvl1_counts = zeros(stims_nb_max, stims_nb_max, 4, nb_unique_tsk_pairs); % counts for both intrvls
    intrvl1_freq1 = zeros(stims_nb_max, stims_nb_max, 4, nb_unique_tsk_pairs);  % counts for intrvl 1
    model_choices = NaN(stims_nb_max, stims_nb_max, 4, nb_unique_tsk_pairs);    % model
    residual_choices = NaN(stims_nb_max, stims_nb_max, 4, nb_unique_tsk_pairs); % model to plot residuals

    for kk = 1:nb_knds
        stim_pair = stims2(kk, :);
        resp_pair = resps2(kk, :);
        choi_pair = choic2(kk, :);
        task_pair = tasks2(kk, :);


        % -> convert vector {00, 01, 10, 11} into index {1, 2, 3, 4}
        resp_ind = bin2dec(char(resp_pair + '0')) + 1;

        % -> find index corresponding to the task pair
        [~, task_ind] = ismember(task_pair, tsk_pairs_vals, 'rows');

        % -> indices of stimuli presented in 1st and 2nd intervals
        stim1_ind = find(stim_pair(1) == stims_dim(task_pair(1),:));
        stim2_ind = find(stim_pair(2) == stims_dim(task_pair(2),:));

        % -> probability of choosing the 1st interval for humans
        human_choices(stim1_ind, stim2_ind, resp_ind, task_ind) = fraction_chosen(kk);

        % -> probability of choosing the 1st interval for model
        if (length(plot_t2_psychometric) > 1)
            model_choices(stim1_ind, stim2_ind, resp_ind, task_ind) = plot_t2_psychometric(kk);
        end

        % -> residuals of choice probabilities
        if (length(plot_t2_residuals) > 1)
            residual_choices(stim1_ind, stim2_ind, resp_ind, task_ind) = plot_t2_residuals(kk);
        end

        % -> frequency of occurence
        intrvl1_freq1(stim1_ind, stim2_ind, resp_ind, task_ind) = choi_pair(1);
        intrvl1_counts(stim1_ind, stim2_ind, resp_ind, task_ind) = total1(kk);
    end        

    counts_max = max(max(max(max(intrvl1_counts))));


    % -> data for psychometric functions split across tasks and intervals
    plot_tsk_intrvl = zeros(1, nb_unique_tsk_pairs);

    for task_ind = 1:nb_unique_tsk_pairs
        task_pair = tsk_pairs_vals(task_ind, :);

        task1 = task_pair(1);
        task2 = task_pair(2);

        nbs1 = stims_nb(task1);
        nbs2 = stims_nb(task2);

        % -> temporary variables
        unsorted_prob_intrvl1  = NaN(1, nbs1);
        unsorted_prob_intrvl2  = NaN(1, nbs2);
        chosen_prob_intrvl1    = NaN(1, nbs1);
        chosen_prob_intrvl2    = NaN(1, nbs2);
        declined_prob_intrvl1  = NaN(1, nbs1);
        declined_prob_intrvl2  = NaN(1, nbs2);
        unsorted_count_intrvl1 = NaN(1, nbs1);
        unsorted_count_intrvl2 = NaN(1, nbs2);
        chosen_count_intrvl1   = NaN(1, nbs1);
        chosen_count_intrvl2   = NaN(1, nbs2);
        declined_count_intrvl1 = NaN(1, nbs1);
        declined_count_intrvl2 = NaN(1, nbs2);

        for stim1_ind = 1:nbs1
            resp0 = sum(sum(intrvl1_counts(stim1_ind, :, [1,2], task_ind)));
            resp1 = sum(sum(intrvl1_counts(stim1_ind, :, [3,4], task_ind)));
            unsorted_count_intrvl1(stim1_ind) = resp0 + resp1;
            unsorted_prob_intrvl1(stim1_ind) = resp1 / (resp0 + resp1);
            chos0 = sum(sum(intrvl1_freq1(stim1_ind, :, [1,2], task_ind)));
            chos1 = sum(sum(intrvl1_freq1(stim1_ind, :, [3,4], task_ind)));
            chosen_count_intrvl1(stim1_ind) = chos0 + chos1;
            chosen_prob_intrvl1(stim1_ind) = chos1 / (chos0 + chos1);
            decl0 = resp0 - chos0;
            decl1 = resp1 - chos1;
            declined_count_intrvl1(stim1_ind) = decl0 + decl1;
            declined_prob_intrvl1(stim1_ind) = decl1 / (decl0 + decl1);
        end

        for stim2_ind = 1:nbs2
            resp0 = sum(sum(intrvl1_counts(:, stim2_ind, [1,3], task_ind)));
            resp1 = sum(sum(intrvl1_counts(:, stim2_ind, [2,4], task_ind)));
            unsorted_count_intrvl2(stim2_ind) = resp0 + resp1;
            unsorted_prob_intrvl2(stim2_ind) = resp1 / (resp0 + resp1);
            decl0 = sum(sum(intrvl1_freq1(:, stim2_ind, [1,3], task_ind)));
            decl1 = sum(sum(intrvl1_freq1(:, stim2_ind, [2,4], task_ind)));
            declined_count_intrvl2(stim2_ind) = decl0 + decl1;
            declined_prob_intrvl2(stim2_ind) = decl1 / (decl0 + decl1);
            chos0 = resp0 - decl0;
            chos1 = resp1 - decl1;
            chosen_count_intrvl2(stim2_ind) = chos0 + chos1;
            chosen_prob_intrvl2(stim2_ind) = chos1 / (chos0 + chos1);
        end


        % -> fill in data structure for output
        if (~any(isnan(unsorted_prob_intrvl1)))
            plot_tsk_intrvl(task_ind) = 1;

            plot_data_intrvl(task_ind).tasks = tsk_pairs_vals(task_ind, :);

            plot_data_intrvl(task_ind).unsorted_prob_intrvl1 = unsorted_prob_intrvl1;
            plot_data_intrvl(task_ind).unsorted_prob_intrvl2 = unsorted_prob_intrvl2;
            plot_data_intrvl(task_ind).chosen_prob_intrvl1   = chosen_prob_intrvl1;
            plot_data_intrvl(task_ind).chosen_prob_intrvl2   = chosen_prob_intrvl2;
            plot_data_intrvl(task_ind).declined_prob_intrvl1 = declined_prob_intrvl1;
            plot_data_intrvl(task_ind).declined_prob_intrvl2 = declined_prob_intrvl2;

            plot_data_intrvl(task_ind).unsorted_count_intrvl1 = unsorted_count_intrvl1;
            plot_data_intrvl(task_ind).unsorted_count_intrvl2 = unsorted_count_intrvl2;
            plot_data_intrvl(task_ind).chosen_count_intrvl1   = chosen_count_intrvl1;
            plot_data_intrvl(task_ind).chosen_count_intrvl2   = chosen_count_intrvl2;
            plot_data_intrvl(task_ind).declined_count_intrvl1 = declined_count_intrvl1;
            plot_data_intrvl(task_ind).declined_count_intrvl2 = declined_count_intrvl2;
        end

    end


    % -> ******************************************************************* <-
    % -> plot psychometric function for perceptual decisions
    goaheadandplot = 0;
    if (islogical(plot_t1_psychometric))
        if (plot_t1_psychometric)
            goaheadandplot = 1;
        end
    elseif (isstruct(plot_t1_psychometric))
        goaheadandplot = 2;
    end

    if (goaheadandplot)
        fig = figure;
        fig_pos = get(fig, 'Position');
        if (nb_unique_tasks == 1)
            set(fig, 'Position', [fig_pos(1) fig_pos(2) 430 340]);
        elseif (nb_unique_tasks == 2)
            set(fig, 'Position', [fig_pos(1) fig_pos(2) 720 340]);
        else
            set(fig, 'Position', [fig_pos(1) fig_pos(2) 720 480]);
        end

        for task_no = 1:nb_unique_tasks
            if (nb_unique_tasks <= 2)
                subplot(1, nb_unique_tasks, task_no);
            else
                subplot(2, ceil(nb_unique_tasks/2), task_no);
            end
            hold on;
            set(gca, 'FontName', 'Arial'); set(gca, 'FontSize', 24);

            nbs = stims_nb(task_no);
            xvals = stims_dim(task_no, 1:nbs);
            xrange = xvals(end) - xvals(1);
            xmin = xvals(1) - xrange/10;
            xmax = xvals(end) + xrange/10;

            unsorted_prob = unsorted_prob_lst(task_no, 1:stims_nb(task_no));
            chosen_prob = chosen_rsp_data_lst(task_no, 1:stims_nb(task_no));
            declined_prob = declined_rsp_data_lst(task_no, 1:stims_nb(task_no));
            unsorted_count = unsorted_count_lst(task_no, 1:stims_nb(task_no));
            chosen_count = chosen_count_lst(task_no, 1:stims_nb(task_no));
            declined_count = declined_count_lst(task_no, 1:stims_nb(task_no));

            line([0 0], [0 1], 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
            line([xmin xmax], [0.5 0.5], 'Color', [0.5 0.5 0.5], 'LineWidth', 2);

            max_size = max(unsorted_count) * 0.4;
            col = [0.5, 0.5, 0.5];
            sizes = 100 * unsorted_count ./ max_size + 0.01;
            sizes(sizes==0) = 0.1;
            scatter(xvals, unsorted_prob, sizes, col, 'filled', 'LineWidth', 2);
            plot(xvals, unsorted_prob, '-', 'Color', col, 'LineWidth', 3);

            max_size = max(max(chosen_count), max(declined_count)) * 0.4;
            col = [0, 0.8, 0.5];
            sizes = 100 * chosen_count ./ max_size + 0.01;
            sizes(sizes==0) = 0.1;
            scatter(xvals, chosen_prob, sizes, col, 'filled', 'LineWidth', 2);
            if (goaheadandplot == 2)
                model = plot_t1_psychometric(task_no).chosen_prob;
                plot(xvals, model, '-', 'Color', col, 'LineWidth', 3);
            else
                plot(xvals, chosen_prob, '-', 'Color', col, 'LineWidth', 3);
            end

            col = [0.8, 0, 0.5];
            sizes = 100 * declined_count ./ max_size + 0.01;
            sizes(sizes==0) = 0.1;
            scatter(xvals, declined_prob, sizes, col, 'filled', 'LineWidth', 2);
            if (goaheadandplot == 2)
                model = plot_t1_psychometric(task_no).declined_prob;
                plot(xvals, model, '--', 'Color', col, 'LineWidth', 3);
            else
                plot(xvals, declined_prob, '--', 'Color', col, 'LineWidth', 3);
            end

            xlim([xmin, xmax]);
            ylim([0, 1]);
            xlabel('Stimulus Strength');
            ylabel('Proportion ''Right''');

            if (nb_unique_tasks > 1)
                title(sprintf('Task %d', task_no));
            end
        end
    end


    % -> ******************************************************************* <-
    % -> plot all (task_pairs x intervals) psychometric functions for perceptual decisions
    goaheadandplot = 0;
    if (islogical(plot_t1_psychometric_all))
        if (plot_t1_psychometric_all)
            goaheadandplot = 1;
        end
    elseif (isstruct(plot_t1_psychometric_all))
        goaheadandplot = 2;
    end

    if (goaheadandplot)

        plot_tsk_nb = sum(plot_tsk_intrvl);
        plot_no = -1;

        fig = figure;
        fig_pos = get(fig, 'Position');

        switch (plot_tsk_nb)
            case 1
                set(fig, 'Position', [fig_pos(1) fig_pos(2) 800 500]);
            case 2
                set(fig, 'Position', [fig_pos(1) fig_pos(2) 800 800]);
            case 4
                set(fig, 'Position', [fig_pos(1) fig_pos(2) 800 1000]);
            otherwise
                set(fig, 'Position', [fig_pos(1) fig_pos(2) 800 800]);
        end

        nbs_vect = NaN(1, 2);
        xvals_mat = NaN(plot_tsk_nb*2, stims_nb_max);
        unsorted_prob_lst = NaN(plot_tsk_nb*2, stims_nb_max);
        chosen_prob_lst = NaN(plot_tsk_nb*2, stims_nb_max);
        declined_prob_lst = NaN(plot_tsk_nb*2, stims_nb_max);
        unsorted_count_lst = NaN(plot_tsk_nb*2, stims_nb_max);
        chosen_count_lst = NaN(plot_tsk_nb*2, stims_nb_max);
        declined_count_lst = NaN(plot_tsk_nb*2, stims_nb_max);
        task_vect = NaN(1, plot_tsk_nb*2);


        for task_ind = 1:nb_unique_tsk_pairs
            task_pair = tsk_pairs_vals(task_ind, :);

            task1 = task_pair(1);
            task2 = task_pair(2);
            nbs_vect(1) = stims_nb(task1);
            nbs_vect(2) = stims_nb(task2);


            if (plot_tsk_intrvl(task_ind) == 1)
                plot_no = plot_no + 2;

                task_vect(plot_no) = task1;
                xvals_mat(plot_no,1:nbs_vect(1)) = stims_dim(task1,1:nbs_vect(1));
                unsorted_prob_lst(plot_no,1:nbs_vect(1)) = plot_data_intrvl(task_ind).unsorted_prob_intrvl1;
                chosen_prob_lst(plot_no,1:nbs_vect(1)) = plot_data_intrvl(task_ind).chosen_prob_intrvl1;
                declined_prob_lst(plot_no,1:nbs_vect(1)) = plot_data_intrvl(task_ind).declined_prob_intrvl1;

                unsorted_count_lst(plot_no,1:nbs_vect(1)) = plot_data_intrvl(task_ind).unsorted_count_intrvl1;
                chosen_count_lst(plot_no,1:nbs_vect(1)) = plot_data_intrvl(task_ind).chosen_count_intrvl1;
                declined_count_lst(plot_no,1:nbs_vect(1)) = plot_data_intrvl(task_ind).declined_count_intrvl1;

                task_vect(plot_no+1) = task2;
                xvals_mat(plot_no+1,1:nbs_vect(2)) = stims_dim(task2,1:nbs_vect(2));
                unsorted_prob_lst(plot_no+1,1:nbs_vect(2)) = plot_data_intrvl(task_ind).unsorted_prob_intrvl2;
                chosen_prob_lst(plot_no+1,1:nbs_vect(2)) = plot_data_intrvl(task_ind).chosen_prob_intrvl2;
                declined_prob_lst(plot_no+1,1:nbs_vect(2)) = plot_data_intrvl(task_ind).declined_prob_intrvl2;

                unsorted_count_lst(plot_no+1,1:nbs_vect(2)) = plot_data_intrvl(task_ind).unsorted_count_intrvl2;
                chosen_count_lst(plot_no+1,1:nbs_vect(2)) = plot_data_intrvl(task_ind).chosen_count_intrvl2;
                declined_count_lst(plot_no+1,1:nbs_vect(2)) = plot_data_intrvl(task_ind).declined_count_intrvl2;

            end
        end


        for plot_no = 1:(plot_tsk_nb*2)

            intrvl = mod((plot_no+1), 2) + 1;
            task_ind = ceil(plot_no/2);

            curr_task = task_vect(plot_no);
            curr_nbs_vect = stims_nb(curr_task);

            xvals = xvals_mat(plot_no, 1:curr_nbs_vect);
            xrange = xvals(end) - xvals(1);
            xmin = xvals(1) - xrange/10;
            xmax = xvals(end) + xrange/10;

            unsorted_prob = unsorted_prob_lst(plot_no, 1:curr_nbs_vect);
            chosen_prob = chosen_prob_lst(plot_no, 1:curr_nbs_vect);
            declined_prob = declined_prob_lst(plot_no, 1:curr_nbs_vect);

            unsorted_count = unsorted_count_lst(plot_no, 1:curr_nbs_vect);
            chosen_count = chosen_count_lst(plot_no, 1:curr_nbs_vect);
            declined_count = declined_count_lst(plot_no, 1:curr_nbs_vect);

            subplot(plot_tsk_nb, 2, plot_no);
            hold on;
            set(gca, 'FontName', 'Arial'); set(gca, 'FontSize', 24);


            line([0 0], [0 1], 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
            line([xmin xmax], [0.5 0.5], 'Color', [0.5 0.5 0.5], 'LineWidth', 2);

            max_size = max(unsorted_count) * 0.4;
            col = [0.5, 0.5, 0.5];
            sizes = 100 * unsorted_count ./ max_size + 0.01;
            scatter(xvals, unsorted_prob, sizes, col, 'filled', 'LineWidth', 2);
            plot(xvals, unsorted_prob, '-', 'Color', col, 'LineWidth', 3);

            max_size = max(max(chosen_count), max(declined_count)) * 0.4;
            col = [0, 0.8, 0.5];
            sizes = 100 * chosen_count ./ max_size + 0.01;
            scatter(xvals, chosen_prob, sizes, col, 'filled', 'LineWidth', 2);
            if (goaheadandplot == 2)
                if (intrvl == 1)
                    model = plot_t1_psychometric_all(task_ind).chosen_prob_intrvl1;
                else
                    model = plot_t1_psychometric_all(task_ind).chosen_prob_intrvl2;
                end
                plot(xvals, model, '-', 'Color', col, 'LineWidth', 3);
            else
                plot(xvals, chosen_prob, '-', 'Color', col, 'LineWidth', 3);
            end

            col = [0.8, 0, 0.5];
            sizes = 100 * declined_count ./ max_size + 0.01;
            scatter(xvals, declined_prob, sizes, col, 'filled', 'LineWidth', 2);
            if (goaheadandplot == 2)
                if (intrvl == 1)
                    model = plot_t1_psychometric_all(task_ind).declined_prob_intrvl1;
                else
                    model = plot_t1_psychometric_all(task_ind).declined_prob_intrvl2;
                end
                plot(xvals, model, '--', 'Color', col, 'LineWidth', 3);
            else
                plot(xvals, declined_prob, '--', 'Color', col, 'LineWidth', 3);
            end

            xlim([xmin, xmax]);
            ylim([0, 1]);
            xlabel('Stimulus Strength');
            ylabel('Proportion ''Right''');

            title(sprintf('Task %d, Intrvl %d', task_vect(plot_no), intrvl));

        end

    end


    % -> ******************************************************************* <-
    % -> plot of fraction chosen for human against model
    if (plot_human_model)
        fig = figure;
        fig_pos = get(fig, 'Position');
        set(fig, 'Position', [fig_pos(1) fig_pos(2) 430 340]);
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



    % -> ******************************************************************* <-
    % -> plot matrix of choices for pairs of stimulus strengths

    ind2resp = ['L', 'R'];
    axes_sizx = 0.25;
    axes_sizy = 0.35;

    if (plot_choi_by_resp)

        for task_ind = 1:nb_unique_tsk_pairs
            choice_mat = squeeze(human_choices(:, :, :, task_ind));
            tasks = tsk_pairs_vals(task_ind, :);
            task_intrvl1 = tasks(1);
            task_intrvl2 = tasks(nb_unique_tasks);

            % -> generate a new figure only if there is something to plot
            if ~all(all(all(isnan(choice_mat))))

                fig = figure;
                fig_pos = get(fig, 'Position');
                set(fig, 'Position', [fig_pos(1) fig_pos(2) 720 480]);

                for resp_knd = 1:4
                    resp_ind = dec2bin(resp_knd-1,2);
                    resp1 = bin2dec(resp_ind(1));
                    resp2 = bin2dec(resp_ind(2));
                    axes_minx = 0.08 + resp1*0.3;
                    axes_miny = 0.1 + resp2*0.5;
                    axes('Position',[axes_minx axes_miny axes_sizx axes_sizy]);

                    nb_x = stims_nb(task_intrvl1);
                    nb_y = stims_nb(task_intrvl2);

                    choi_im = choice_mat(1:nb_x, 1:nb_y, resp_knd);
                    choi_im = flipud(choi_im');
                    im = imagesc(choi_im, [0, 1]);
                    set(im, 'AlphaData', ~isnan(choi_im));
                    axis('square');
                    set(gca, 'FontName', 'Arial'); set(gca, 'FontSize', 12);

                    uu = 1:nb_x;
                    set(gca,'XTick', uu);
                    set(gca,'XTickLabel', round(stims_dim(task_intrvl1,uu),2));
                    vv = 1:nb_y;
                    set(gca,'YTick', vv);
                    set(gca,'YTickLabel', fliplr(round(stims_dim(task_intrvl2,vv),2)));

                    if (resp_knd == 1)
                        text(0.41, -0.2, 'Interval 1 Stimulus Strength', ...
                            'FontSize', 20, 'HorizontalAlignment', 'left', ...
                            'Units', 'normalized');
                        text(-0.25, 0.46, 'Interval 2 Stimulus Strength', ...
                            'FontSize', 20, 'HorizontalAlignment', 'left', ...
                            'Units', 'normalized',...
                            'Rotation', 90);
                    end

                    title(sprintf('(D_1=%c, D_2=%c)', ind2resp(resp1+1), ind2resp(resp2+1)));
                end


                axes('Position',[0.7 0.6 axes_sizx axes_sizy]);

                choi_intrvl1 = sum(squeeze(intrvl1_freq1(:, :, :, task_ind)), 3);
                choi_counts = sum(squeeze(intrvl1_counts(:, :, :, task_ind)), 3);
                choi_im = choi_intrvl1 ./ choi_counts;
                choi_im = choi_im(1:nb_x, 1:nb_y);
                choi_im = flipud(choi_im');
                im = imagesc(choi_im, [0, 1]);
                set(im, 'AlphaData', ~isnan(choi_im));
                axis('square');

                xlabel('Interval 1 Stimulus Strength');
                ylabel('Interval 2 Stimulus Strength');   
                set(gca, 'FontName', 'Arial'); set(gca, 'FontSize', 12);
                uu = 1:nb_x;
                set(gca,'XTick', uu);
                set(gca,'XTickLabel', round(stims_dim(task_intrvl1,uu),2));
                vv = 1:nb_y;
                set(gca,'YTick', vv);
                set(gca,'YTickLabel', fliplr(round(stims_dim(task_intrvl2,vv),2)));
                title('Any Response');

                colormap('cool');
                axes('Position',[0.5 0.1 axes_sizx axes_sizy]);
                axis('off');
                hcb = colorbar;
                set(gca, 'FontName', 'Arial'); set(gca, 'FontSize', 16);
                ylabel(hcb, {'Probability to'; 'Choose Interval 1'}, ...
                    'FontName', 'Arial', 'FontSize', 20);

                if (nb_unique_tasks > 1)
                    text(2.25, 0.1, sprintf('Tasks:\n(%d, %d)', task_intrvl1, task_intrvl2), ...
                        'FontSize', 20, 'HorizontalAlignment', 'center', ...
                        'Units', 'normalized');
                end
            end
        end

    end


    % -> ******************************************************************* <-
    % -> plot Type-2 psychometric functions
    if (any(plot_t2_psychometric))

        for task_ind = 1:nb_unique_tsk_pairs
            choice_mat = squeeze(human_choices(:, :, :, task_ind));
            tasks = tsk_pairs_vals(task_ind, :);

            % -> generate a new figure only if there is something to plot
            if ~all(all(all(isnan(choice_mat))))

                fig = figure;
                fig_pos = get(fig, 'Position');
                set(fig, 'Position', [fig_pos(1) fig_pos(2) 720 480]);

                task_intrvl1 = tasks(1);
                task_intrvl2 = tasks(2);
                nbs1 = stims_nb(task_intrvl1);
                nbs2 = stims_nb(task_intrvl2);
                xvals = stims_dim(task_intrvl1, 1:nbs1);
                yvals = stims_dim(task_intrvl2, 1:nbs2);

                xrange = xvals(end) - xvals(1);
                xmin = xvals(1) - xrange/10;
                xmax = xvals(end) + xrange/10;

                if (plot_zscore)
                    zvals = [0.05, 0.2, 0.5, 0.8, 0.95];
                else
                    zvals = 0:0.2:1;
                end
                xvals_mat = repmat(xvals, nbs2, 1);
                col_mat = cool(nbs2);

                for resp_knd = 1:4
                    resp_ind = dec2bin(resp_knd-1,2);
                    resp1 = bin2dec(resp_ind(1));
                    resp2 = bin2dec(resp_ind(2));
                    axes_minx = 0.08 + resp1*0.3;
                    axes_miny = 0.1 + resp2*0.5;
                    axes('Position',[axes_minx axes_miny axes_sizx axes_sizy]);

                    hold on;
                    if (plot_zscore)
                        line([xmin, xmax], [0 ,0], ...
                            'Color', [0.5 0.5 0.5], 'LineWidth', 2);
                    else
                        line([xmin, xmax], [0.5 ,0.5], ...
                            'Color', [0.5 0.5 0.5], 'LineWidth', 2);
                    end

                    human_choices_knd = choice_mat(:, :, resp_knd);
                    human_choices_knd = human_choices_knd(1:nbs1, 1:nbs2);
                    human_choices_knd = human_choices_knd';
                    if (plot_zscore)
                        human_choices_knd = norminv(human_choices_knd);
                    end

                    plot(xvals_mat', human_choices_knd', ':', 'LineWidth', 1);

                    choi_counts = squeeze(intrvl1_counts(1:nbs1, 1:nbs2, resp_knd, task_ind));
                    choi_counts = choi_counts';

                    sizes = 100 * choi_counts / counts_max;
                    sizes(sizes==0) = NaN;

                    if (length(plot_t2_psychometric) > 1)
                        model_choices_knd = squeeze(model_choices(:, :, resp_knd, task_ind));
                        model_choices_knd = model_choices_knd(1:nbs1, 1:nbs2);
                        model_choices_knd = model_choices_knd';

                        if (plot_zscore)
                            model_choices_knd = norminv(model_choices_knd);
                        end
                    end

                    for ii=1:nbs2
                        col = col_mat(ii,:);
                        scatter(xvals, human_choices_knd(ii,:), sizes(ii,:), col, 'filled', 'LineWidth', 2);

                        if (length(plot_t2_psychometric) > 1)
                            plot(xvals, model_choices_knd(ii,:), '-', 'Color', col, 'LineWidth', 2);
                        end
                    end


                    set(gca, 'FontName', 'Arial'); set(gca, 'FontSize', 12);

                    set(gca,'XTick', xvals, 'XTickLabel', xvals);
                    xlim([xmin xmax]);
                    if (plot_zscore)
                        set(gca,'YTick', norminv(zvals), 'YTickLabel', zvals);
                        ylim([-2.5 2.5]);
                    else
                        set(gca,'YTick', zvals, 'YTickLabel', zvals);
                        ylim([0 1]);
                    end
                    axis('square');

                    if (resp_knd == 1)
                        text(0.41, -0.2, 'Interval 1 Stimulus Strength', ...
                            'FontSize', 20, 'HorizontalAlignment', 'left', ...
                            'Units', 'normalized');
                        text(-0.25, 0.46, 'Probability to Choose Interval 1', ...
                            'FontSize', 20, 'HorizontalAlignment', 'left', ...
                            'Units', 'normalized',...
                            'Rotation', 90);
                    end

                    title(sprintf('(D_1=%c, D_2=%c)', ind2resp(resp1+1), ind2resp(resp2+1)));
                end

                axes('Position',[0.57 0.1 axes_sizx axes_sizy]);
                hold on;
                set(gca, 'FontName', 'Arial'); set(gca, 'FontSize', 18);
                for ii=1:nbs2
                    col = col_mat(ii,:);
                    plot(0, 0, 'LineWidth', 3, 'Color', col);
                end
                axis('off');
                legend(string(yvals));

                text(0.75, 1.1, sprintf('Interval 2\nStimulus Strength'), ...
                    'FontSize', 20, 'HorizontalAlignment', 'center', ...
                    'Units', 'normalized');
                if (nb_unique_tasks > 1)
                    text(1.25, 0.1, sprintf('Tasks:\n(%d, %d)', task_intrvl1, task_intrvl2), ...
                        'FontSize', 20, 'HorizontalAlignment', 'center', ...
                        'Units', 'normalized');
                end

            end
        end
    end


    % -> ******************************************************************* <-
    % -> plot Type-2 residuals
    if (length(plot_t2_residuals) > 1)

%         for task_ind = 1:4
        for task_ind = 1:nb_unique_tsk_pairs
            choice_mat = squeeze(human_choices(:, :, :, task_ind));
            tasks = tsk_pairs_vals(task_ind, :);

            % -> generate a new figure only if there is something to plot
            if ~all(all(all(isnan(choice_mat))))

                fig = figure;
                fig_pos = get(fig, 'Position');
                set(fig, 'Position', [fig_pos(1) fig_pos(2) 720 480]);

                nbs1 = stims_nb(1);
                nbs2 = stims_nb(nb_unique_tasks);
                xvals = stims_dim(1, 1:nbs1);
                yvals = stims_dim(nb_unique_tasks, 1:nbs2);
                col_mat = cool(nbs2);

                diff_mat = NaN(stims_nb(1), stims_nb(nb_unique_tasks), 4);
                diff_cnt = NaN(stims_nb(1), stims_nb(nb_unique_tasks), 4);


                subplot(1, 2, 1);
                hold on;
                line([xvals(1), xvals(end)]*1.1, [0 ,0], 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
                for resp_knd = 1:4

                    human_choices_knd = choice_mat(:, :, resp_knd);
                    human_choices_knd = human_choices_knd';

                    choi_counts = squeeze(intrvl1_counts(:, :, resp_knd, task_ind));
                    choi_counts = choi_counts';
                    sizes = 100 * choi_counts / counts_max;

                    sizes(sizes==0) = NaN;

                    model_choices_knd = squeeze(residual_choices(:, :, resp_knd, task_ind));
                    model_choices_knd = model_choices_knd';

                    diff_vals = human_choices_knd - model_choices_knd;

                    for ii=1:nbs2
                        col = col_mat(ii,:);

                        if (resp_knd > 2)
                            scatter(xvals, diff_vals(ii,:), sizes(ii,:), col, 'filled', 'LineWidth', 2);
                            plot(xvals, diff_vals(ii,:), ':', 'LineWidth', 1);
                            diff_mat(:, ii, resp_knd) = diff_vals(ii,:);
                            diff_cnt(:, ii, resp_knd) = sizes(ii,:);
                        else
                            scatter(-xvals, diff_vals(ii,:), sizes(ii,:), col, 'filled', 'LineWidth', 2);
                            plot(-xvals, diff_vals(ii,:), ':', 'LineWidth', 1);
                            diff_mat(:, ii, resp_knd) = fliplr(diff_vals(ii,:));
                            diff_cnt(:, ii, resp_knd) = fliplr(sizes(ii,:));
                        end
                    end

                    zval_max = round(max(max(max(abs(diff_vals)))), 2);
                    zvals = linspace(-zval_max, zval_max, 9);

                    xlabel('Interval 1 Stimulus');
                    ylabel('Prob Interval 1 Choice');
                    set(gca, 'FontName', 'Arial'); set(gca, 'FontSize', 12);


                    set(gca,'XTick', xvals, 'XTickLabel', xvals);
                    set(gca,'YTick', zvals, 'YTickLabel', zvals);
                    ylim([-zval_max zval_max]);
                    axis('square');


                    if (nb_unique_tasks > 1)
                        caption = sprintf('Tasks: (%d, %d)', ...
                            tasks(1), tasks(2));
                        title(caption, 'FontSize', 14);   
                    end

                end

                count = nansum(nansum(diff_cnt,3), 2);
                diff_mean = nansum(nansum(diff_mat .* diff_cnt,3), 2) ./ count;
                plot(xvals, diff_mean, 'o-', 'MarkerSize', 12, 'MarkerFaceColor', [0,0,0], ...
                    'Color', [0,0,0], 'LineWidth', 3);
                err2 = diff_cnt.^2 .* diff_mat.^2;
                err2 = nansum(nansum(nansum(err2, 3), 2));
                fprintf('err2 = %7.3f\n', err2);


                sbpl = subplot(1, 2, 2);
                sbpl.Position = sbpl.Position - [0.1, 0.2, 0.0, 0];
                hold on;
                set(gca, 'FontName', 'Arial'); set(gca, 'FontSize', 18);
                for ii=1:nbs2
                    col = col_mat(ii,:);
                    plot(0, 0, 'LineWidth', 3, 'Color', col);
                end
                axis('off');
                legend(string(yvals));

                title('               Interval 2 Stimulus', 'FontSize', 14);
            end
        end
    end

end


% -> THE END <-
% ------------------------------------------------------------------------

