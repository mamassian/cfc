% function to group trials in a confidence forced-choice experiment
%
% INPUT:
%   raw_data: matrix of:
%       1st col: stimulus intensity interval 1
%       2nd col: stimulus intensity interval 2
%       3rd col: perceptual decision interval 1 (1 = 'A',  0 = 'B')
%       4th col: perceptual decision interval 2 (1 = 'A',  0 = 'B')
%       5th col: confidence choice for interval 1 (1 = chosen,  0 = declined)
%       6th col: confidence choice for interval 2 (1 = chosen,  0 = declined)
%       7th col: stimulus task for interval 1 (optional)
%       8th col: stimulus task for interval 2 (optional)
%  Possible parameter:
%       'bin': ask to bin stimuli for each task in (value) bins
%       'intrvl_sym': are the two intervals symmetric? (either
%                       'true' (default) or 'false')
%
% OUTPUT:
%   wrap_data: grouped data per (s1, s2, r1, r2):
%       1st col: stimulus intensity interval 1
%       2nd col: stimulus intensity interval 2
%       3rd col: perceptual decision interval 1
%       4th col: perceptual decision interval 2
%       5th col: nb of confidence choices for interval 1
%       6th col: nb of confidence choices for interval 2
%       7th col: stimulus task for interval 1
%       8th col: stimulus task for interval 2
%   type2_resp: Type 2 choice probability to select interval 1 (rather than 2)
%
% EXAMPLES:
%   wrap_data = cfc_group_04(raw_data);
%   [wrap_data, type2_resp] = cfc_group_04(raw_data);
%   wrap_data = cfc_group_04(raw_data, 'bins', 5);
%   
%
% 06-JAN-2019 - pascal mamassian
% 28-JAN-2019 - pm: changed to list for all rows of 'stims'
% 15-FEB-2019 - pm: generalize for different sensory tasks
% 17-MAR-2019 - pm: add optional argument to bin data
% 29-MAR-2019 - pm: generalize for staircase values

function [wrap_data, type2_resp] = cfc_group_04(raw_data, varargin)

    % -> default optional arguments
    dflt_nb_bins = 0;          % 0 is a message to not bin input data
    dflt_intrvl_sym = true;    % intervals are symmetric, so we can group trials across them
    
    % -> parse all arguments
    ip = inputParser;
    valid_nb_bins = @(xx) isnumeric(xx) && isscalar(xx) && (xx > 0);
    addRequired(ip, 'raw_data', @isnumeric);
    addParameter(ip, 'bins', dflt_nb_bins, valid_nb_bins);
    addParameter(ip, 'intrvl_sym', dflt_intrvl_sym, @islogical);
    parse(ip, raw_data, varargin{:});
    nb_bins = ip.Results.bins;
    intrvl_sym = ip.Results.intrvl_sym;

    nb_trials = size(raw_data, 1);

    % -> make sure the task is explicit for each interval
    if (size(raw_data, 2) < 7)
        % -> by default, assume that there is only one sensory task
        tasks_nn = repmat([1, 1], nb_trials, 1);
        raw_data = [raw_data, tasks_nn];
    end
    
    % -> check how many tasks were run
    tsk_vals = unique(raw_data(:, 7:8));
    nb_tasks = length(tsk_vals);

    if (nb_bins > 0)
        stm_data = [raw_data(:, [1, 7]); raw_data(:, [2, 8])];
        stm_mean = NaN(nb_tasks, nb_bins);
        inds_mat = NaN(nb_trials, 2);   % bin indices for the two intervals
        for tt = 1:nb_tasks
            inds3 = (stm_data(:,2) == tt);
            inds4 = find(inds3);
            stm_tsk_vals = stm_data(inds3);
            qq = quantile(stm_tsk_vals, nb_bins - 1);
            [~, ~, stm_tsk_bin] = histcounts(stm_tsk_vals, [-inf, qq, inf]);
            for nn = 1:nb_bins
                inds5 = stm_tsk_bin == nn;
                stm_mean(tt, nn) = mean(stm_tsk_vals(inds5));
                inds6 = inds4(inds5);
                inds1 = (inds6 <= nb_trials);    % indices in 1st interval
                inds2 = (inds6 > nb_trials);     % indices in 2nd interval
                inds_mat(inds6(inds1), 1) = nn;
                inds_mat(inds6(inds2)-nb_trials, 2) = nn;
            end
        end
                
        % -> replace raw stimulus value by mean of bin
        for rr = 1:nb_trials
            tsk1 = raw_data(rr, 7);
            tsk2 = raw_data(rr, 8);
            raw_data(rr, 1) = stm_mean(tsk1, inds_mat(rr, 1));
            raw_data(rr, 2) = stm_mean(tsk2, inds_mat(rr, 2));
        end
    end
    
    % -> extract the interval chosen
    type2_choice_prob = raw_data(:, 5) ./ sum(raw_data(:, 5:6), 2);
    
    [knd_vals, ~, knd_ic] = unique(raw_data(:, [1:4, 7:8]), 'rows');
    nb_knds = size(knd_vals, 1);
    type2_resp = NaN(nb_knds, 1);
    wrap_data = NaN(nb_knds, 8);
    wrap_data(:, [1:4, 7:8]) = knd_vals;

    for ww = 1:nb_knds
        inds = (knd_ic == ww);     % indices that have the same kind
        type2_choice_vals = type2_choice_prob(inds);

        nn1 = sum(type2_choice_vals);
        nn0 = sum(1 - type2_choice_vals);
        
        type2_resp(ww) = nn1 / (nn1 + nn0);
        wrap_data(ww, 5:6) = [nn1, nn0];
    end
    
    if (intrvl_sym)
        % -> combine symmetric pairs (s2, s1) -> (s1, s2)
        nb_knds2 = 0;
        wrap_data2 = NaN(size(wrap_data));
        knd_vals2 = NaN(size(knd_vals));
        for ww = 1:nb_knds
            uu = knd_vals(ww, :);
            vv = [uu(2), uu(1), uu(4), uu(3), uu(6), uu(5)];
            [aa, bb] = ismember(vv, knd_vals2, 'rows');
            if (aa)
                wrap_data2(bb, [5,6]) = wrap_data2(bb, [5,6]) + wrap_data(ww, [6,5]);
            else
                nb_knds2 = nb_knds2 + 1;
                knd_vals2(nb_knds2, :) = uu;
                wrap_data2(nb_knds2, :) = wrap_data(ww, :);
            end
        end
        wrap_data = wrap_data2(1:nb_knds2, :);
        type2_resp = wrap_data(:, 5) ./ sum(wrap_data(:, 5:6), 2);
    end
    
end
