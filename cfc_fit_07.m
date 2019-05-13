% wrap function for confidence model fit
%
% INPUT:
%   wrap_data: grouped data per (s1, s2, r1, r2):
%       1st col: stimulus intensity interval 1
%       2nd col: stimulus intensity interval 2
%       3rd col: perceptual decision interval 1
%       4th col: perceptual decision interval 2
%       5th col: nb of confidence choice for interval 1
%       6th col: nb of confidence choice for interval 2
%       7th col: stimulus task for interval 1
%       8th col: stimulus task for interval 2
%   eff_key: efficiency key:
%       0 = just compute Type 1 & ideal perf. & super-ideal perf.
%       1 = compute efficiency (default)
%       2 = compute full model assuming crit2 = crit1
%       3 = compute full model (noise2, crit2, boost2)
%
% OUTPUT: conf_st = struct
%   conf_st.param_noise1:       sdtev of noise for Type 1 decision (0 = perfectly sensitive)
%   conf_st.param_crit1:        criterion for Type 1 decision (criterion)
%   conf_st.fitted_perf1_uns:   unsorted perf (fitted)
%   conf_st.fitted_perf1_idl:   ideal perf
%   conf_st.fitted_choi2_idl:   ideal choice
%   conf_st.tmp_noise2eq:       equivalent noise 2 for ideal perf
%   conf_st.fitted_perf1_supidl:    super-ideal perf
%   conf_st.fitted_choi2_supidl:    super-ideal choice
%   conf_st.param_efficiency:       efficiency
%   conf_st.param_noise2:           sdtev of noise for Type 2 decision (0 = ideal)
%   conf_st.param_crit2:            criterion for Type 2 decision (criterion)
%   conf_st.param_boost2:           boost2 (fraction super-ideal)
%   conf_st.param_conf_bias:        conf_bias (bias in favour of task 2)
%   conf_st.param_intrvl1_bias:     intrvl_bias (bias in favour of interval 1)
%
% 04-JUN-2018: pascal mamassian
% 29-DEC-2018: pm: take into account type 2 criterion
% 21-JAN-2019: pm: adapted for cfc format
% 16-FEB-2019 - pm: generalize for different sensory tasks
% 13-MAR-2019 - pm: normalize confidence evidence by sensory sensitivity
% 21-MAR-2019 - pm: do not fit second noise2 parameter (undetermined)


function cfc_struct = cfc_fit_07(wrap_data, varargin)

    % -> start timer
    tic;
    
    debug1 = 0;
%     debug1 = 1;
    
    % -> default optional arguments
    dflt_efficiency       = true;	% compute efficiency
    dflt_conf_noise_boost = false;	% split confidence noise and boost
    dflt_conf_crit        = false;  % add confidence criterion
    dflt_conf_bias        = false;  % add confidence bias
    dflt_intrvl_bias      = false;  % add interval bias
    
    % -> parse all arguments
    ip = inputParser;
    addRequired(ip, 'wrap_data', @isnumeric);
    addParameter(ip, 'efficiency', dflt_efficiency, @islogical);
    addParameter(ip, 'conf_noise_boost', dflt_conf_noise_boost, @islogical);
    addParameter(ip, 'conf_crit', dflt_conf_crit, @islogical);
    addParameter(ip, 'conf_bias', dflt_conf_bias, @islogical);
    addParameter(ip, 'intrvl_bias', dflt_intrvl_bias, @islogical);
    parse(ip, wrap_data, varargin{:});
    compute_efficiency = ip.Results.efficiency;
    split_conf_noise_boost = ip.Results.conf_noise_boost;
    add_conf_crit = ip.Results.conf_crit;
    add_conf_bias = ip.Results.conf_bias;
    add_intrvl_bias = ip.Results.intrvl_bias;


    % -> order of parameters
    ind_1noise1 = 1;
    ind_1crit1  = 2;
    ind_2noise1 = 3;
    ind_2crit1  = 4;
    ind_1noise2 = 5;
    ind_1boost  = 6;
    ind_2noise2 = 7;
    ind_2boost  = 8;
    ind_1crit2  = 9;
    ind_2crit2 = 10;
    ind_bias2  = 11;
    ind_intrvl = 12;
%     ind_effic  = 13;
    params_nb = 12;

    stims = wrap_data(:, 1:2);
    type1_resp = wrap_data(:, 3:4);
    nn1_choices = wrap_data(:, 5:6);
    tasks = wrap_data(:, 7:8);
    total2 = repmat(sum(nn1_choices, 2), 1, 2);
    stimtask = complex(stims, tasks);
    [stm_vals, ~, stm_ic] = unique(stimtask);
    nb_unique_stims2 = size(stm_vals, 1);

    % -> starting values and bounds for parameters:
    % -> (1)1noise1, (2)1crit1, (3)2noise1, (4)2crit1
    % -> (5)1noise2, (6)1boost, (7)2noise2, (8)2boost, (9)1crit2, (10)2crit2 
    % -> (11)bias2, (12)intrvl
%     params0 = [1.0, 0.1, 1.0, 0.1, ...
%         1.7, 0.2, 1.7, 0.2, 0.1, 0.1, ...
%         0.1, 0.5];
    params0 = [5.0, 0.1, 1.0, 0.1, ...
        1.7, 0.8, 1.7, 0.2, 0.1, 0.1, ...
        0.1, 0.5];
    paramsLB = [0.0, -Inf, 0.0, -Inf, ...
        0.0, 0.0, 0.0, 0.0, -Inf, -Inf, ...
        -Inf, 0.0];
    paramsUB = [Inf, Inf, Inf, Inf, ...
        Inf, 1.0, Inf, 1.0, Inf, Inf, ...
        Inf, 1.0];
    params_reset = [1.0, 0.0, 1.0, 0.0, ...
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
        0.0, 0.5];

    
    % -> display on-going information
    verbose_on = 1;
%     verbose_on = 0;

    options1 = optimset;
    if (verbose_on)
    %     options2 = optimset('Display','iter');    % <<
        options2 = optimset('Display','iter', 'TolFun',1e-3, 'TolX',1e-3);
    else
        options2 = optimset('TolFun',1e-3, 'TolX',1e-3);
    end

    noise2_idl = NaN;
    noise2_dat = NaN;
    efficiency = NaN;
    chosen_eff = NaN;
    nn1_eff = NaN;
    chosen_full = NaN;
    nn1_full = NaN;
    
    
    model_params = params_reset;

    
    % -> step 0: check how many tasks were run
    [tsk_vals, ~, ~] = unique(tasks);
    nb_unique_tasks = size(tsk_vals, 1);

    % -> Type 1 parameters
    sens_noise_task = NaN(1, nb_unique_tasks);  % sensory noise, per task
    sens_crit_task = NaN(1, nb_unique_tasks);  % sensory criterion, per task

    % -> Type 2 parameters
    conf_noise_task = NaN(1, nb_unique_tasks);  % confidence noise, per task
    conf_crit_task = NaN(1, nb_unique_tasks);  % confidence criterion, per task
    conf_boost_task = NaN(1, nb_unique_tasks);  % confidence boost, per task
    
    conf2_bias = NaN;
    intrvl1_bias = NaN;


    % -> step 1: fit Type 1 performance
    % ->   collect all cases where stimulus is identical
    nn1_list = NaN(nb_unique_stims2, 2);     % (nn1, nn0)
    for ww = 1:nb_unique_stims2
        inds = (stm_ic == ww);     % indices that have the same stim
        resp_vals = type1_resp(inds);
        resp_nb = total2(inds);     % nb. of responses for each type
        nn1_list(ww, 1) = sum(resp_vals .* resp_nb);
        nn1_list(ww, 2) = sum((1 - resp_vals) .* resp_nb);
    end

    params0_0 = params0(1:2);       % noise1, crit1
    paramsLB_0 = paramsLB(1:2);     paramsUB_0 = paramsUB(1:2);
    for tt = 1:nb_unique_tasks
        inds = find(imag(stm_vals) == tt);
        stm_tsk_vals = real(stm_vals(inds, :));
        nn1_tsk_list = nn1_list(inds, :);
        
        my_fun_0 = @(pp) normcdf(stm_tsk_vals, pp(2), pp(1));
        [param_type1, loglike] = fitnll(my_fun_0, nn1_tsk_list, ...
            params0_0, paramsLB_0, paramsUB_0, options1);

        sens_noise_task(tt) = param_type1(1);
        sens_crit_task(tt)  = param_type1(2);
        fprintf('Unsorted: (crit1, noise1) = (%7.3f, %7.3f)\n', ...
            sens_crit_task(tt), sens_noise_task(tt));
        fprintf('Log-Likelihood = %7.3f\n', loglike);
        
        % -> store fitted params for step 2
        if (tt == 1)
            model_params(ind_1noise1) = param_type1(1);
            model_params(ind_1crit1)  = param_type1(2);
        else
            model_params(ind_2noise1) = param_type1(1);
            model_params(ind_2crit1)  = param_type1(2);
        end
    end
    

    % -> step 2: get simulated performance for ideal (noise2 = 0; boost2 = 0)
    model_params(ind_1noise2) = 0.0;
    model_params(ind_1boost) = 0.0;
    model_params(ind_1crit2) = model_params(ind_1crit1);
    if (nb_unique_tasks > 1)
        model_params(ind_2noise2) = 0.0;
        model_params(ind_2boost) = 0.0;
        model_params(ind_2crit2) = model_params(ind_2crit1);
    end
    [chosen_idl, nn1_idl] = cfc_core_08(wrap_data, model_params);
    
    % -> super-ideal
    model_params(ind_1boost) = 1.0;
    if (nb_unique_tasks > 1)
        model_params(ind_2boost) = 1.0;
    end
    [chosen_sdl, nn1_sdl] = cfc_core_08(wrap_data, model_params);

    
            
    % -> compute the confidence efficiency
    if (compute_efficiency)

        % -> step 3: get the equivalent noise 2 for this ideal performance
        % -> efficiency is defined relative to (boost2 = 1)
        wrap_idl = wrap_data;
        wrap_idl(:, 5:6) = nn1_idl;
        % -> (1)1noise1, (2)1crit1, (3)2noise1, (4)2crit1
        % -> (5)1noise2, (6)1boost, (7)2noise2, (8)2boost, (9)1crit2, (10)2crit2 
        % -> (11)bias2, (12)intrvl
        params_set = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0];
        switch nb_unique_tasks
            case 1
                fixed_vals = [sens_noise_task, sens_crit_task, NaN(1,2), ...
                    1.0, NaN(1,2), sens_crit_task, NaN, params_reset(11:12)];
            case 2
                fixed_vals = [sens_noise_task(1), sens_crit_task(1), ...
                    sens_noise_task(2), sens_crit_task(2), 1.0, NaN, 1.0, ...
                    sens_crit_task(1), sens_crit_task(2), params_reset(11:12)];
        end
%         switch nb_unique_tasks
%             case 1
%                 params_set = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0];
%                 fixed_vals = [sens_noise_task, sens_crit_task, NaN(1,2), ...
%                     1.0, NaN(1,2), sens_crit_task, NaN, params_reset(11:12)];
%             case 2
%                 params_set = [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0];
%                 fixed_vals = [sens_noise_task(1), sens_crit_task(1), ...
%                     sens_noise_task(2), sens_crit_task(2), 1.0, 1.0, ...
%                     sens_crit_task(1), sens_crit_task(2), params_reset(11:12)];
%         end
        params_set = logical(params_set);
        params0_idl = params0(params_set);
        paramsLB_idl = paramsLB(params_set);
        paramsUB_idl = paramsUB(params_set);
        
        
        my_fun_idl = @(pp) cfc_core_wrap(wrap_idl, pp, params_set, fixed_vals);
        
        
        [noise2_idl, loglike] = fitnll(my_fun_idl, nn1_idl, ...
            params0_idl, paramsLB_idl, paramsUB_idl, options2);
        fprintf('Ideal: noise2 = %7.3f\n', noise2_idl);
        fprintf('Log-Likelihood = %7.3f\n', loglike);
% aa = cfc_core_08(wrap_idl, define_params(noise2_idl));      
% sum(log(aa) .* wrap_idl(:,5) + log(1 - aa) .* wrap_idl(:,6))

if (debug1)
    bb_lst = 0.1:0.2:0.9;
    bb_nb = length(bb_lst);
    nse_lst1 = NaN(1, bb_nb);
    for bb_ind = 1:bb_nb
        bb_val = bb_lst(bb_ind);
        switch nb_unique_tasks
            case 1
                fixed_vals2 = [sens_noise_task, sens_crit_task, NaN(1,2), ...
                    bb_val, NaN(1,2), sens_crit_task, NaN, params_reset(11:12)];
            case 2
                fixed_vals2 = [sens_noise_task(1), sens_crit_task(1), ...
                    sens_noise_task(2), sens_crit_task(2), bb_val, NaN, bb_val, ...
                    sens_crit_task(1), sens_crit_task(2), params_reset(11:12)];
        end
        my_fun_idl = @(pp) cfc_core_wrap(wrap_idl, pp, params_set, fixed_vals2);

        [noise2_idl2, loglike] = fitnll(my_fun_idl, nn1_idl, ...
            params0_idl, paramsLB_idl, paramsUB_idl, options2);

        nse_lst1(bb_ind) = noise2_idl2;
        fprintf('Ideal (%7.3f): noise2 = %7.3f\n', bb_val, noise2_idl2);
        fprintf('Log-Likelihood = %7.3f\n', loglike);
    end

    fig = figure;
    set(fig, 'Position', [875 455 790 350]);
    subplot(1,2,1);
    hold on;
    set(gca, 'FontName', 'Arial'); set(gca, 'FontSize', 24);
    xvals = nse_lst1;
    yvals = bb_lst;
    plot(xvals, yvals, 'ko', 'MarkerSize', 12);
    plot(noise2_idl, 1.0, 'r+', 'MarkerSize', 12);
    xlabel('noise');
    ylabel('boost');
    xlim([0, 3]);
    yy = 0.0:.1:1.0;
    % xx = spline(yvals, xvals, yy);
    xx = pchip(yvals, xvals, yy);
    plot(xx, yy, 'b-', 'LineWidth', 2);
    title('Ideal');
end

        % -> step 4: efficiency computation
        % -> get the noise2 for the actual data, assuming (boost2 = 1)
        my_fun_eff = @(pp) cfc_core_wrap(wrap_data, pp, params_set, fixed_vals);

        [noise2_dat, loglike] = fitnll(my_fun_eff, nn1_choices, ...
            params0_idl, paramsLB_idl, paramsUB_idl, options2);

        % -> use ratio of variance (std squared) for definition of efficiency
        efficiency = (noise2_idl ./ noise2_dat).^2;
        
        fprintf('Data: noise2 = %7.3f\n', noise2_dat);
        fprintf('Efficiency = %7.3f\n', efficiency);
        fprintf('Log-Likelihood = %7.3f\n', loglike);

if (debug1)        
    nse_lst2 = NaN(1, bb_nb);
    for bb_ind = 1:bb_nb
        bb_val = bb_lst(bb_ind);
        switch nb_unique_tasks
            case 1
                fixed_vals2 = [sens_noise_task, sens_crit_task, NaN(1,2), ...
                    bb_val, NaN(1,2), sens_crit_task, NaN, params_reset(11:12)];
            case 2
                fixed_vals2 = [sens_noise_task(1), sens_crit_task(1), ...
                    sens_noise_task(2), sens_crit_task(2), bb_val, NaN, bb_val, ...
                    sens_crit_task(1), sens_crit_task(2), params_reset(11:12)];
        end
        my_fun_eff = @(pp) cfc_core_wrap(wrap_data, pp, params_set, fixed_vals2);

        [noise2_dat2, loglike] = fitnll(my_fun_eff, nn1_choices, ...
            params0_idl, paramsLB_idl, paramsUB_idl, options2);

        nse_lst2(bb_ind) = noise2_dat2;
        fprintf('Data (%7.3f): noise2 = %7.3f\n', bb_val, noise2_dat2);
        fprintf('Log-Likelihood = %7.3f\n', loglike);
    end

    subplot(1,2,2);
    hold on;
    set(gca, 'FontName', 'Arial'); set(gca, 'FontSize', 24);
    xvals = nse_lst2;
    yvals = bb_lst;
    plot(xvals, yvals, 'ko', 'MarkerSize', 12);
    plot(noise2_dat, 1.0, 'r+', 'MarkerSize', 12);
    xlabel('noise');
    ylabel('boost');
    yy = 0.0:.1:1.0;
    % xx = spline(yvals, xvals, yy);
    xx = pchip(yvals, xvals, yy);
    plot(xx, yy, 'b-', 'LineWidth', 2);
    title('Data');
    
    eff_lst = (nse_lst1 ./ nse_lst2).^2;
    fprintf('%7.3f \n', eff_lst);

end        
        
        
        [chosen_eff, nn1_eff] = my_fun_eff(noise2_dat);
        
        % -> by default, in case we don't run the full model, take this
        chosen_full = chosen_eff;
        
    end
    
    
    % -> compute full model
    if (split_conf_noise_boost || add_conf_crit || add_conf_bias || add_intrvl_bias)
        
        % -> (1)1noise1, (2)1crit1, (3)2noise1, (4)2crit1
        % -> (5)1noise2, (6)1boost, (7)2noise2, (8)2boost, (9)1crit2, (10)2crit2 
        % -> (11)bias2, (12)intrvl
        params_set = zeros(1, params_nb);
        fixed_vals = NaN(1, params_nb);
        
        % -> check which parameters we should fit
        switch nb_unique_tasks
            case 1
                fixed_vals(ind_1noise1) = sens_noise_task;
                fixed_vals(ind_1crit1) = sens_crit_task;
                fixed_vals(ind_2noise1) = NaN;
                fixed_vals(ind_2crit1) = NaN;
                params_set(ind_1noise2) = 1;
                fixed_vals(ind_2noise2) = NaN;
                if (split_conf_noise_boost)
                    params_set(ind_1boost) = 1;
                else
                    fixed_vals(ind_1boost) = 1.0;
                end
                fixed_vals(ind_2boost) = NaN;
                if (add_conf_crit)
                    params_set(ind_1crit2) = 1;
                else
                    fixed_vals(ind_1crit2) = sens_crit_task;
                end
                fixed_vals(ind_2crit2) = NaN;
            case 2
                fixed_vals(ind_1noise1) = sens_noise_task(1);
                fixed_vals(ind_1crit1) = sens_crit_task(1);
                fixed_vals(ind_2noise1) = sens_noise_task(2);
                fixed_vals(ind_2crit1) = sens_crit_task(2);
                params_set(ind_1noise2) = 1;
                fixed_vals(ind_2noise2) = NaN;
                if (split_conf_noise_boost)
                    params_set(ind_1boost) = 1;
                    params_set(ind_2boost) = 1;
                else
                    fixed_vals(ind_1boost) = 1.0;
                    fixed_vals(ind_2boost) = 1.0;
                end
                if (add_conf_crit)
                    params_set(ind_1crit2) = 1;
                    params_set(ind_2crit2) = 1;
                else
                    fixed_vals(ind_1crit2) = sens_crit_task(1);
                    fixed_vals(ind_2crit2) = sens_crit_task(2);
                end
        end
        if (add_conf_bias)
            params_set(ind_bias2) = 1;
        else
            fixed_vals(ind_bias2) = 0.0;
        end
        if (add_intrvl_bias)
            params_set(ind_intrvl) = 1;
        else
            fixed_vals(ind_intrvl) = 0.5;
        end
%         if ((compute_efficiency) && (~split_conf_noise_boost))
%             params_set(ind_effic) = 1;
%             params_set(6:7) = zeros(1,2);
%             fixed_vals(6:7) = NaN(1,2);
%         else
%             fixed_vals(ind_effic) = NaN;
%         end
        params_set = logical(params_set);
        fixed_vals = fixed_vals(~params_set);
        params0full = params0(params_set);
        paramsLBfull = paramsLB(params_set);
        paramsUBfull = paramsUB(params_set);


        my_fun_full = @(pp) cfc_core_wrap(wrap_data, pp, params_set, fixed_vals);


        % -> run the fit by maximum-likelihood
        [paramBest, loglike] = fitnll(my_fun_full, nn1_choices, ...
            params0full, paramsLBfull, paramsUBfull, options2);

        
        % -> (1)1noise1, (2)1crit1, (3)2noise1, (4)2crit1
        % -> (5)1noise2, (6)1boost, (7)2noise2, (8)2boost, (9)1crit2, (10)2crit2 
        % -> (11)bias2, (12)intrvl
        params_ind = find(params_set == 1);
        fprintf('Best fit of full model\n');
        switch nb_unique_tasks
            case 1
                [pp_on, pp_ind] = ismember([ind_1noise2, ind_1boost], params_ind); 
                if (pp_on)
                    conf_noise_task = paramBest(pp_ind(1));
                    conf_boost_task = paramBest(pp_ind(2));
                    fprintf(' (noise2, boost2) = (%7.3f, %7.3f)\n', ...
                            conf_noise_task, conf_boost_task);
                end
                    
                [pp_on, pp_ind] = ismember(ind_1crit2, params_ind); 
                if (pp_on)
                    conf_crit_task = paramBest(pp_ind);
                    fprintf(' crit2 = %7.3f\n', conf_crit_task);
                end
            case 2
                [pp_on, pp_ind] = ismember([ind_1noise2, ind_1boost, ...
                    ind_2boost], params_ind); 
                if (pp_on)
                    conf_noise_task = paramBest(pp_ind(1));
                    conf_boost_task = paramBest([pp_ind(2), pp_ind(3)]);
                    fprintf(' noise2 (shared between tasks 1 & 2) = %7.3f\n', ...
                        conf_noise_task);
                    fprintf(' boost2 for (task 1, task 2) = (%7.3f, %7.3f)\n', ...
                        conf_boost_task(1), conf_boost_task(2));
                end
                
                [pp_on, pp_ind] = ismember([ind_1crit2, ind_2crit2], params_ind); 
                if (pp_on)
                    conf_crit_task = paramBest(pp_ind);
                    fprintf(' crit2 = %7.3f (task 1) and %7.3f (task 2)\n', ...
                        conf_crit_task(1), conf_crit_task(2));
                end
                
                [pp_on, pp_ind] = ismember(ind_bias2, params_ind); 
                if (pp_on)
                    conf2_bias = paramBest(pp_ind);
                    fprintf(' confidence bias (overcfd for task 2 rel. to task 1) = %7.3f\n', ...
                        conf2_bias);
                end
        end
        [pp_on, pp_ind] = ismember(ind_intrvl, params_ind); 
        if (pp_on)
            intrvl1_bias = paramBest(pp_ind);
            fprintf(' bias to choose 1st interval (0.5 is no bias) = %7.3f\n', ...
                intrvl1_bias);
        end
%         [pp_on, pp_ind] = ismember(ind_effic, params_ind); 
%         if (pp_on)
%             noise2_dat = paramBest(pp_ind);
%             efficiency = noise2_idl ./ noise2_dat;
%             fprintf(' efficiency = %7.3f\n', efficiency);
%         end


        fprintf('Log-Likelihood = %7.3f\n', loglike);

        [chosen_full, nn1_full] = my_fun_full(paramBest);
        
    end

    
    % -> stop timer
    fprintf('\n');
    toc;
    fprintf('\n********************\n\n');

        
    % -> fill in struct
    cfc_struct = struct;
    cfc_struct.param_noise1 = sens_noise_task;      % sensory noise
    cfc_struct.param_crit1 = sens_crit_task;        % sensory criterion
    cfc_struct.chosen_ideal = chosen_idl;           % ideal choice
    cfc_struct.nn1_ideal = nn1_idl;                 % ideal responses
    cfc_struct.chosen_super_ideal = chosen_sdl;     % super-ideal choice
    cfc_struct.nn1_super_ideal = nn1_sdl;           % super-ideal responses
    
    cfc_struct.tmp_noise2eq = noise2_idl;           % equivalent noise 2 for ideal perf
    cfc_struct.tmp_noise2dat = noise2_dat;          % equivalent noise 2 for human perf
    cfc_struct.param_efficiency = efficiency;       % efficiency
    cfc_struct.chosen_eff = chosen_eff;             % efficiency choice
    cfc_struct.nn1_eff = nn1_eff;                   % efficiency responses
    
    cfc_struct.param_noise2 = conf_noise_task;      % sdtev of noise for Type 2 decision (0 = ideal)
    cfc_struct.param_crit2 = conf_crit_task;        % criterion for Type 2 decision (criterion)
    cfc_struct.param_boost2 = conf_boost_task;      % fraction super-ideal (1 - fraction ideal)
    cfc_struct.chosen_full = chosen_full;           % full-model choice
    cfc_struct.nn1_full = nn1_full;                 % full-model responses

    cfc_struct.param_conf_bias = conf2_bias;        % conf_bias (bias in favour of task 2)
    cfc_struct.param_intrvl1_bias = intrvl1_bias;    % intrvl_bias (bias in favour of interval 1)

    
    % ---------------------------------------------------------------------
    % -> define nested functions
    
    % -> wrapper around cfc_core function to flexibly add parameters
    function [pred_resp, pred_nn1] = cfc_core_wrap(wrap_data, variable_prms, params_set, fixed_vals)
      params = zeros(size(params_set));
      params(~params_set) = fixed_vals;
      params(params_set) = variable_prms;
      [pred_resp, pred_nn1] = cfc_core_08(wrap_data, params);
    end

    
    % ->   fit by maximizing the log-likelihood
    function [params_best, loglike] = fitnll(fit_fcn, nn1_lst, params_0, params_LB, params_UB, fit_options)

        fun = @(xx) loglikefcn(xx, fit_fcn, nn1_lst);
        
        [params_best, nll_best] = fminsearchbnd(fun, params_0, params_LB, params_UB, fit_options);
        loglike = - nll_best;
        
        % ->   negative summed log-likelihood
        function nglglk = loglikefcn(pp, ff, nn)
            ypred = ff(pp);
% if ((any(ypred == 0)) || (any(ypred == 1)))
%     fprintf('err: %7.3f\n', ypred');
% end
            ypred(ypred == 0) = 1e-6;
            ypred(ypred == 1) = 1 - 1e-6;
%             ypred(ypred == 0) = 1e-2;
%             ypred(ypred == 1) = 1 - 1e-2;

            ll1 = log(ypred);
            ll0 = log(1.0 - ypred);
            ll_vect = nn(:,1) .* ll1 + nn(:,2) .* ll0;  % vector of log-likelihoods
            nglglk = - sum(ll_vect);      % minimize (-log) likelihood
            
%             ydata = nn(:,1) ./ (nn(:,1) + nn(:,2));
%             sqrerr = (ydata - ypred).^2;
%             nglglk = sum(sqrerr);
        end
        
    end

end

