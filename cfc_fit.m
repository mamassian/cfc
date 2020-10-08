% CONFIDENCE FORCED-CHOICE TOOLBOX  v0.3
%
% cfc_fit
%   This is the wrap function for confidence model fit.
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
% OPTIONAL PARAMETERS:
%
%   'model_parameters': model free parameters, as a structure
%                       increment numbers for desired parameters
%                       (use vectors if more than 1 task,
%                       use NaN for padding):
%       'sens_noise' : sensory (Type 1) sdtev of noise (0 = perfectly sensitive)
%       'sens_crit'  : sensory (Type 1) criterion
%       'conf_noise' : confidence (Type 2) sdtev of noise (0 = ideal)
%       'conf_boost' : fraction super-ideal (0 = ideal, 1 = super-ideal)
%       'conf_crit'  : confidence (Type 2) criterion
%       'intrvl_bias': bias in favour of interval 2 (over interval 1)
%       'conf_bias'  : overconfidence relative to one of the tasks
%
%   'model_fixed_values': model fixed values parameters, as a structure
%                       (use vectors if more than 1 task):
%       'tasks_list' : vector of tasks, e.g. '1' or '[1, 2]'
%       'sens_noise' : sensory (Type 1) sdtev of noise (0 = perfectly sensitive)
%       'sens_crit'  : sensory (Type 1) criterion
%       'conf_noise' : confidence (Type 2) sdtev of noise (0 = ideal)
%       'conf_boost' : fraction super-ideal (0 = ideal, 1 = super-ideal)
%       'conf_crit'  : confidence (Type 2) criterion
%       'intrvl_bias': bias in favour of interval 2 (over interval 1)
%       'conf_bias'  : overconfidence relative to one of the tasks
%
%   'boost_init': redo the fit from multiple boost starting points
%       (provide list here as a vector) to avoid local minima
%
%   'skip_efficiency': skip computing efficiency (for faster computations)
%
%   'verbose': verbose flag:
%       0: remove all online reports
%       1: print parameter estimates
%       2: print parameter estimates and fitting progress
%
%
% OBSOLETE:
%
%   'efficiency': compute confidence efficiency
%
%   'conf_noise_boost': attempt to split the efficiency into its
%       "confidence noise" and "confidence boost" components
%
%   'conf_crit': estimate confidence criterion (default assume no bias)
%
%   'conf_bias': estimate confidence bias (default assume no bias)
%
%   'intrvl_bias': estimate interval bias (default assume no bias)
%
%   'set_conf_boost': sometimes, it is interesting to test a model where the
%       confidence boost is fixed at some value (typically 0 or 1)
%
%
% OUTPUT: cfc_struct = struct
%   cfc_struct.sens_noise        sdtev of sensory noise (0 = perfect sensitivity)
%   cfc_struct.sens_crit         sensory criterion
%   cfc_struct.conf_noise        sdtev of confidence noise (0 = ideal)
%   cfc_struct.conf_boost        confidence boost (fraction super-ideal)
%   cfc_struct.conf_crit         confidence criterion
%   cfc_struct.intrvl_bias       intrvl_bias (bias in favour of interval 1)
%   cfc_struct.conf_bias         confidence bias (relative to a task set to 1.0)
%   cfc_struct.efficiency        efficiency
%   cfc_struct.choice_prob_ideal            ideal choice probabilities
%   cfc_struct.choice_prob_super_ideal      super-ideal choice probabilities
%   cfc_struct.choice_prob_eff              efficiency choice probabilities
%   cfc_struct.choice_prob_model            full-model choice probabilities
%   cfc_struct.loglike                      log-likelihood of best fit
%
% -> obsolete
%   cfc_struct.chosen_noise1       sensory noise for chosen
%   cfc_struct.chosen_crit1        sensory criterion for chosen
%   cfc_struct.nn1_ideal           ideal responses
%   cfc_struct.nn1_super_ideal     super-ideal responses
%   cfc_struct.tmp_noise2eq        equivalent noise 2 for ideal perf
%   cfc_struct.tmp_noise2dat       equivalent noise 2 for human perf
%   cfc_struct.nn1_eff             efficiency responses
%   cfc_struct.nn1_full            full-model responses
% 
% EXAMPLES OF USE:
%   cfc_struct = cfc_fit(grouped_data)
%   cfc_struct = cfc_fit(grouped_data, 'boost2_init', [0.2, 0.5, 0.8])
% 
%
% 21-SEP-2020 - pascal mamassian

% 04-JUN-2018: pascal mamassian
% 29-DEC-2018: pm: take into account type 2 criterion
% 21-JAN-2019: pm: adapted for cfc format
% 16-FEB-2019 - pm: generalize for different sensory tasks
% 13-MAR-2019 - pm: normalize confidence evidence by sensory sensitivity
% 21-MAR-2019 - pm: do not fit second noise2 parameter (undetermined)
% 17-AUG-2019 - pm: fixed type 2 criterion
% 03-FEB-2020 - pm: changed definition of confidence criterion
% 29-FEB-2020 - pm: cleaned up
% 10-JUN-2020 - pm: fast normcdf
% 29-JUN-2020 - pm: added verbose flag
% 04-JUL-2020 - pm: updated arguments with model parameters
% 21-SEP-2020 - pm: cleaned up



function cfc_struct = cfc_fit(grouped_data, varargin)

    % -> default optional arguments
    dflt_model            = struct; % pass on specific model
    dflt_fixed            = struct; % pass on fixed values
    dflt_boost2_init      = [];     % initial value (or list of values) for boost2
    dflt_skip_efficiency  = false;  % skip computing efficiency
    dflt_verbose          = 1;      % verbose flag

    
    % -> parse all arguments
    ip = inputParser;
    ip.StructExpand = false;
    addRequired(ip, 'grouped_data', @isnumeric);
    addParameter(ip, 'model_parameters', dflt_model, @isstruct);
    addParameter(ip, 'model_fixed_values', dflt_fixed, @isstruct);
    addParameter(ip, 'boost_init', dflt_boost2_init, @isnumeric);
    addParameter(ip, 'skip_efficiency', dflt_skip_efficiency, @islogical);
    addParameter(ip, 'verbose', dflt_verbose, @isnumeric);
    
    
    parse(ip, grouped_data, varargin{:});
    my_model_params = ip.Results.model_parameters;
    my_fixed_params = ip.Results.model_fixed_values;
    boost_init_lst = ip.Results.boost_init;
    skip_efficiency = ip.Results.skip_efficiency;
    verbose_flag = ip.Results.verbose;
    
    
    compute_efficiency = ~skip_efficiency;	% compute efficiency


    options1 = optimset;
    if (verbose_flag >= 2)
        options2 = optimset('Display','iter', 'TolFun',1e-3, 'TolX',1e-3);
        tstart = tic;   % start timer
    else
        options2 = optimset('TolFun',1e-3, 'TolX',1e-3);
    end

    % -> default values of the parameters
    default_sens_noise  = 1.0;
    default_sens_crit   = 0.0;
    default_conf_noise  = 0.0;
    default_conf_boost  = 0.0;
    default_conf_crit   = 0.0;
    default_intrvl_bias = 0.0;
    default_conf_bias   = 1.0;
    
    % -> initial values of the parameters
%     initial_sens_noise  = 1.2;
%     initial_sens_crit   = 0.1;
    initial_sens_noise  = 3.2;
    initial_sens_crit   = 4.0;
    initial_conf_noise  = 0.5;
    initial_conf_boost  = 0.2;
    initial_conf_crit   = 0.4;
    initial_intrvl_bias = 0.6;
    initial_conf_bias   = 1.5;
    
    % -> lower and upper bound values of the parameters
    lo_bnd_sens_noise   = 0.0;      hi_bnd_sens_noise   = Inf;
    lo_bnd_sens_crit    = -Inf;     hi_bnd_sens_crit    = Inf;
    lo_bnd_conf_noise   = 0.0;      hi_bnd_conf_noise   = Inf;
    lo_bnd_conf_boost   = 0.0;      hi_bnd_conf_boost   = 1.0;
    lo_bnd_conf_crit    = -Inf;     hi_bnd_conf_crit    = Inf;
    lo_bnd_intrvl_bias  = -Inf;     hi_bnd_intrvl_bias  = Inf;
    lo_bnd_conf_bias    = -Inf;     hi_bnd_conf_bias    = Inf;
    

% -> ******************************************************************* <-
% -> step 0: extract data and initialize internal variables
    
    % -> un-wrap data matrix
    stims2 = grouped_data(:, 1:2);
    resps2 = grouped_data(:, 3:4);
    choic2 = grouped_data(:, 5:6);
    tasks2 = grouped_data(:, 7:8);
    
    % -> nb trials per (s1, s2, r1, r2)
    total1 = sum(choic2, 2);
    total2 = repmat(total1, 1, 2);  % duplicate in a (1,2) vector
    
    % -> list of unique stimuli across all tasks
    stimtask = complex(stims2, tasks2);
    [stm_vals, ~, stm_ic] = unique(stimtask);
    nb_unique_stims2 = size(stm_vals, 1);

    % -> check how many tasks were run
    tsk_vals = unique(tasks2)';
    nb_unique_tasks = length(tsk_vals);
    
    % -> Type 1 parameters
    sens_noise_task = NaN(1, nb_unique_tasks);  % sensory noise, per task
    sens_crit_task = NaN(1, nb_unique_tasks);  % sensory criterion, per task

    % -> Type 2 parameters
    conf_noise_task = NaN(1, nb_unique_tasks);  % confidence noise, per task
    conf_crit_task = NaN(1, nb_unique_tasks);  % confidence criterion, per task
    conf_boost_task = NaN(1, nb_unique_tasks);  % confidence boost, per task
    
    conf2_bias = NaN;
    intrvl1_bias = NaN;

    chosen_full = NaN;
    nn1_full = NaN;
    
    
    
    % -> 'params_set' is a struct that defines the free parameters, and 
    %    their order. Put '0' if the variable is not a free parameter.
    default_params_set = struct;
    default_params_set.tasks_list = zeros(1, nb_unique_tasks);
    default_params_set.sens_noise = zeros(1, nb_unique_tasks);
    default_params_set.sens_crit  = zeros(1, nb_unique_tasks);
    default_params_set.conf_noise = zeros(1, nb_unique_tasks);
    default_params_set.conf_boost = zeros(1, nb_unique_tasks);
    default_params_set.conf_crit  = zeros(1, nb_unique_tasks);
    default_params_set.intrvl_bias = 0;
    default_params_set.conf_bias  = zeros(1, nb_unique_tasks);
    
    % -> default values of the parameters
    default_params_val = struct;
    default_params_val.tasks_list = 1:nb_unique_tasks;
    default_params_val.sens_noise = ones(1, nb_unique_tasks) * default_sens_noise;
    default_params_val.sens_crit  = ones(1, nb_unique_tasks) * default_sens_crit;
    default_params_val.conf_noise = ones(1, nb_unique_tasks) * default_conf_noise;
    default_params_val.conf_boost = ones(1, nb_unique_tasks) * default_conf_boost;
    default_params_val.conf_crit  = ones(1, nb_unique_tasks) * default_conf_crit;
    default_params_val.intrvl_bias = default_intrvl_bias;
    default_params_val.conf_bias  = ones(1, nb_unique_tasks) * default_conf_bias;
    
    % -> initial values of the parameters
    initial_params = struct;
    initial_params.tasks_list = 1:nb_unique_tasks;
    initial_params.sens_noise = ones(1, nb_unique_tasks) * initial_sens_noise;
    initial_params.sens_crit  = ones(1, nb_unique_tasks) * initial_sens_crit;
    initial_params.conf_noise = ones(1, nb_unique_tasks) * initial_conf_noise;
    initial_params.conf_boost = ones(1, nb_unique_tasks) * initial_conf_boost;
    initial_params.conf_crit  = ones(1, nb_unique_tasks) * initial_conf_crit;
    initial_params.intrvl_bias = initial_intrvl_bias;
    initial_params.conf_bias  = ones(1, nb_unique_tasks) * initial_conf_bias;
    
    % -> lower and upper bound values of the parameters
    lo_bnd_params = struct;
    lo_bnd_params.tasks_list = 1:nb_unique_tasks;
    lo_bnd_params.sens_noise = ones(1, nb_unique_tasks) * lo_bnd_sens_noise;
    lo_bnd_params.sens_crit  = ones(1, nb_unique_tasks) * lo_bnd_sens_crit;
    lo_bnd_params.conf_noise = ones(1, nb_unique_tasks) * lo_bnd_conf_noise;
    lo_bnd_params.conf_boost = ones(1, nb_unique_tasks) * lo_bnd_conf_boost;
    lo_bnd_params.conf_crit  = ones(1, nb_unique_tasks) * lo_bnd_conf_crit;
    lo_bnd_params.intrvl_bias = lo_bnd_intrvl_bias;
    lo_bnd_params.conf_bias  = ones(1, nb_unique_tasks) * lo_bnd_conf_bias;
    
    hi_bnd_params = struct;
    hi_bnd_params.tasks_list = 1:nb_unique_tasks;
    hi_bnd_params.sens_noise = ones(1, nb_unique_tasks) * hi_bnd_sens_noise;
    hi_bnd_params.sens_crit  = ones(1, nb_unique_tasks) * hi_bnd_sens_crit;
    hi_bnd_params.conf_noise = ones(1, nb_unique_tasks) * hi_bnd_conf_noise;
    hi_bnd_params.conf_boost = ones(1, nb_unique_tasks) * hi_bnd_conf_boost;
    hi_bnd_params.conf_crit  = ones(1, nb_unique_tasks) * hi_bnd_conf_crit;
    hi_bnd_params.intrvl_bias = hi_bnd_intrvl_bias;
    hi_bnd_params.conf_bias  = ones(1, nb_unique_tasks) * hi_bnd_conf_bias;
    
    
    % -> check how many parameters we should fit
    my_params_cell = struct2cell(my_model_params);
    my_params_mat = cell2mat(my_params_cell');
    param_free_nb = max(my_params_mat);

    % -> check which parameters we should fit
    params_set = default_params_set;    % reset params_set
    params_set1 = default_params_set;   % params_set for only Type 1
    params_set2 = default_params_set;   % params_set for only Type 2
    param_free_nb1 = 0;     % nb free parameters for Type 1
    param_free_nb2 = 0;     % nb free parameters for Type 2

    my_fld_nms = fieldnames(my_model_params);
    for my_kk = 1:param_free_nb
        my_fld_ind = cellfun(@(xx) find(xx == my_kk), my_params_cell, 'UniformOutput', false);
        my_uu = find(~cellfun(@isempty, my_fld_ind));

        % -> allow for multiple use of a variable across model parameters
        for my_pp = 1:length(my_uu)
            my_vv = my_fld_ind{my_uu(my_pp)};

            % -> allow for multiple use of a variable within model parameters
            if (strcmp('sens_noise', my_fld_nms{my_uu(my_pp)}) || ...
                        strcmp('sens_crit', my_fld_nms{my_uu(my_pp)}))
                param_free_nb1 = param_free_nb1 + 1;
            else
                param_free_nb2 = param_free_nb2 + 1;
            end
            
            for my_qq = 1:length(my_vv)
                params_set.(my_fld_nms{my_uu(my_pp)})(my_vv(my_qq)) = my_kk;
                if (strcmp('sens_noise', my_fld_nms{my_uu(my_pp)}) || ...
                        strcmp('sens_crit', my_fld_nms{my_uu(my_pp)}))
                    params_set1.(my_fld_nms{my_uu(my_pp)})(my_vv(my_qq)) = ...
                        my_kk;
                else
                    params_set2.(my_fld_nms{my_uu(my_pp)})(my_vv(my_qq)) = ...
                        my_kk - param_free_nb1;
                end
            end
        end
    end


    % -> if (part of) a model was set, pick the corresponding parameters
    % -> check which values are fixed
    fixed_set = default_params_set;

    % -> use the best fit of Type 1 to estimate Type 2 parameters
    fixed_vals = default_params_val;

    my_params_cell = struct2cell(my_fixed_params);

    my_fld_nms = fieldnames(my_fixed_params);
    my_fld_ind = cellfun(@(xx) ~isnan(xx), my_params_cell, 'UniformOutput', false);
    my_uu = find(~cellfun(@isempty, my_fld_ind));

    % -> allow for multiple use of a variable across model parameters
    for my_pp = 1:length(my_uu)
        my_ww = find(my_fld_ind{my_uu(my_pp)});

        % -> allow for multiple use of a variable within model parameters
        for my_qq = my_ww

            % -> remember that this variable was set by user
            fixed_set.(my_fld_nms{my_uu(my_pp)})(my_qq) = 1;

            % -> remember the value of this variable
            fixed_vals.(my_fld_nms{my_uu(my_pp)})(my_qq) = ...
                my_fixed_params.(my_fld_nms{my_uu(my_pp)})(my_qq);
        end
    end



% -> ******************************************************************* <-
% -> step 1: fit Type 1 performance
    
    % -> collect all cases where stimulus is identical (within a task)
    % -> and split in two according to perceptual response
    nn1_list = NaN(nb_unique_stims2, 2);    % (nn1, nn0)
    for ww = 1:nb_unique_stims2
        inds = (stm_ic == ww);     % indices that have the same stim
        resp_vals = resps2(inds);
        resp_nb = total2(inds);     % nb. of responses for each type
        
        nn1_list(ww, 1) = sum(resp_vals .* resp_nb);
        nn1_list(ww, 2) = sum((1 - resp_vals) .* resp_nb);
    end
    
    % -> check that Type 1 is indeed fitted
    params_set_type1 = default_params_set;
    if (isfield(my_model_params, 'sens_noise'))
        params_set_type1.sens_noise = my_model_params.sens_noise;
    else
        params_set_type1.sens_noise = 1:nb_unique_tasks;
    end
    if (isfield(my_model_params, 'sens_crit'))
        params_set_type1.sens_crit = my_model_params.sens_crit;
    else
        params_set_type1.sens_crit = (nb_unique_tasks+1):(2*nb_unique_tasks);
    end
    
    % -> fit cumulative Gaussian
    params_set_cumul = params_set_type1;
    params_set_cumul.sens_noise = 1;
    params_set_cumul.sens_crit = 2;
    params0_cumul  = extract_params_from_struct(params_set_cumul, initial_params);
    paramsLB_cumul = extract_params_from_struct(params_set_cumul, lo_bnd_params);
    paramsUB_cumul = extract_params_from_struct(params_set_cumul, hi_bnd_params);

    for tt = 1:nb_unique_tasks
        tsk = tsk_vals(tt);
        
        inds = find(imag(stm_vals) == tsk);
        stm_tsk_vals = real(stm_vals(inds, :));
        nn1_tsk_list = nn1_list(inds, :);

        my_fun_0 = @(pp) basic_normcdf(stm_tsk_vals, pp(2), pp(1));
        [param_type1, loglike] = fitnll(my_fun_0, nn1_tsk_list, ...
            params0_cumul, paramsLB_cumul, paramsUB_cumul, options1);

        sens_noise_task(tt) = param_type1(1);
        sens_crit_task(tt)  = param_type1(2);
        if (verbose_flag >= 2)
            fprintf('Unsorted trials for task %d: (crit1, noise1) = (%7.3f, %7.3f)\n', ...
                tsk, sens_crit_task(tt), sens_noise_task(tt));
            fprintf('  Log-Likelihood = %7.3f\n\n', loglike);
        end

    end
    
    % -> store fitted params for step 2
    model_params_type1 = default_params_val;
    model_params_type1.tasks_list = tsk_vals;
    model_params_type1.sens_noise = sens_noise_task;
    model_params_type1.sens_crit  = sens_crit_task;

    % -> set initial values obtained from single task fit
    params0_0 = extract_params_from_struct(params_set_type1, model_params_type1);
    paramsLB_0 = extract_params_from_struct(params_set_type1, lo_bnd_params);
    paramsUB_0 = extract_params_from_struct(params_set_type1, hi_bnd_params);

    % -> redo fit of Type 1 using the model passed on
    xx_vals = NaN(nb_unique_stims2, 2);
    xx_ind = 0;
    nn1_tsk_list = NaN(nb_unique_stims2, 2);
    for tt = 1:nb_unique_tasks
        tsk = tsk_vals(tt);
        
        inds = find(imag(stm_vals) == tsk);
        stm_tsk_vals = real(stm_vals(inds, :));
        
        xx_ind_new = xx_ind + length(stm_tsk_vals);
        xx_vals((xx_ind+1):xx_ind_new, 1) = tsk;
        xx_vals((xx_ind+1):xx_ind_new, 2) = stm_tsk_vals;
        
        nn1_tsk_list((xx_ind+1):xx_ind_new, :) = nn1_list(inds, :);
        
        xx_ind = xx_ind_new;
    end
    my_fun_0 = @(pp) cfc_type1(xx_vals, pp, params_set_type1, fixed_vals);
    [param_type1, loglike] = fitnll(my_fun_0, nn1_tsk_list, ...
        params0_0, paramsLB_0, paramsUB_0, options1);

    
    % -> store fitted params for step 2
    model_params_type2 = pack_params_in_struct(param_type1, params_set_type1, fixed_vals);

    % -> fill in fitted values
    sens_noise_task = model_params_type2.sens_noise;
    sens_crit_task  = model_params_type2.sens_crit;

    
% -> ******************************************************************* <-
% -> step 2: get simulated performance for ideal (noise2 = 0; boost2 = 0)
    
    [chosen_idl, nn1_idl] = cfc_core(grouped_data, model_params_type2);
    
    % -> super-ideal
    model_params_type2.conf_boost = ones(1, nb_unique_tasks);

    [chosen_sdl, nn1_sdl] = cfc_core(grouped_data, model_params_type2);

    % -> put back in ideal state
    model_params_type2.conf_boost = zeros(1, nb_unique_tasks);
    
    % -> extract parameters of Type 1
    paramBest1 = extract_params_from_struct(params_set1, model_params_type2);

        
% -> ******************************************************************* <-
% -> step 3: get the equivalent confidence noise for this ideal performance
% -> efficiency is defined relative to (boost = 1)
    
    if (compute_efficiency)
        wrap_idl = grouped_data;
        wrap_idl(:, 5:6) = nn1_idl;

        params_set_idl = default_params_set;
        fixed_vals_idl = model_params_type2;
        
        if (param_free_nb2 > 0)
            params_set_idl.conf_noise = params_set2.conf_noise;
        else
            params_set_idl.conf_noise = 1:nb_unique_tasks;
        end
        
        % -> set confidence boost to 1.0
        fixed_vals_idl.conf_boost = ones(1, nb_unique_tasks);
                
        params0_idl  = extract_params_from_struct(params_set_idl, initial_params);
        paramsLB_idl = extract_params_from_struct(params_set_idl, lo_bnd_params);
        paramsUB_idl = extract_params_from_struct(params_set_idl, hi_bnd_params);
        
        my_fun_idl = @(pp) cfc_core_wrap(wrap_idl, pp, params_set_idl, fixed_vals_idl);
        
        [noise2_idl, loglike] = fitnll(my_fun_idl, nn1_idl, ...
            params0_idl, paramsLB_idl, paramsUB_idl, options2);
        if (verbose_flag >= 1)
            for tt = 1:length(noise2_idl)
                tsk = tsk_vals(tt);
                fprintf('Ideal for task %d:', tsk);
                fprintf(' Equivalent confidence noise = %7.3f\n', noise2_idl(tt));
            end
            fprintf('  Log-Likelihood = %7.3f\n', loglike);
        end


% -> ******************************************************************* <-
% -> step 4: efficiency computation

        % -> get the equivalent confidence noise for the actual data, assuming (boost = 1)
        my_fun_eff = @(pp) cfc_core_wrap(grouped_data, pp, params_set_idl, fixed_vals_idl);

        [noise2_dat, loglike] = fitnll(my_fun_eff, choic2, ...
            params0_idl, paramsLB_idl, paramsUB_idl, options2);

        % -> use ratio of variance (std dev squared) for definition of efficiency
        efficiency = (noise2_idl ./ noise2_dat).^2;
        
        if (verbose_flag >= 1)
            for tt = 1:length(noise2_dat)
                tsk = tsk_vals(tt);
                fprintf('Data for task %d:', tsk);
                fprintf(' Equivalent confidence noise = %7.3f\n', noise2_dat(tt));
            end
            fprintf('  Log-Likelihood = %7.3f\n', loglike);
            for tt = 1:length(efficiency)
                tsk = tsk_vals(tt);
                fprintf('Efficiency for task %d = %7.3f\n', tsk, efficiency(tt));
            end
            fprintf('\n');
        end

        [chosen_eff, nn1_eff] = my_fun_eff(noise2_dat);
        
    else
        noise2_idl = NaN;
        noise2_dat = NaN;
        efficiency = NaN;
        chosen_eff = NaN;
        nn1_eff = NaN;
    end
        
    
% -> ******************************************************************* <-
    % -> step 5: compute full model
    
    if (param_free_nb2 > 0)

        % -> use the best fit of Type 1 to estimate Type 2 parameters
        fixed_vals.sens_noise = model_params_type2.sens_noise;
        fixed_vals.sens_crit  = model_params_type2.sens_crit;

        
        % -> redo the fit from multiple boost2 starting points to avoid local minima
        if (isempty(boost_init_lst))
            boost_init_lst = initial_conf_boost;
        end
        boost2_init_nb = length(boost_init_lst);
        paramBest_mat = NaN(boost2_init_nb, param_free_nb2);
        loglike_lst = NaN(boost2_init_nb, 1);

        % -> do the fit for multiple starting values of 'boost2'
        for bb = 1:boost2_init_nb
            boost2_init_val = boost_init_lst(bb);

            params0_bb = initial_params;
            params0_bb.conf_boost = ones(1, nb_unique_tasks) * boost2_init_val;

            params0full  = extract_params_from_struct(params_set2, params0_bb);
            paramsLBfull = extract_params_from_struct(params_set2, lo_bnd_params);
            paramsUBfull = extract_params_from_struct(params_set2, hi_bnd_params);

            my_fun_full = @(pp) cfc_core_wrap(grouped_data, pp, params_set2, fixed_vals);


            % -> run the fit by maximum-likelihood
            [paramBest, loglike] = fitnll(my_fun_full, choic2, ...
                params0full, paramsLBfull, paramsUBfull, options2);

            paramBest_mat(bb, :) = paramBest;
            loglike_lst(bb) = loglike;

            if ((boost2_init_nb > 1) && (verbose_flag >= 1))
                fprintf('For boost2_init = %5.3f, loglike = %7.3f\n', ...
                    boost2_init_val, loglike);
            end
            
            if ((boost2_init_nb > 1) && (verbose_flag >= 1))
                paramBest_struct = pack_params_in_struct(paramBest, params_set2, fixed_vals);
                conf_noise_vals = paramBest_struct.conf_noise;
                conf_boost_vals = paramBest_struct.conf_boost;
                for tt = 1:nb_unique_tasks
                    fprintf('  Task%d: (noise2, boost2) = (%7.3f, %7.3f)\n', ...
                        tt, conf_noise_vals(tt), conf_boost_vals(tt));
                end
            end

        end     % loop on param boost init

        % -> pick the initial boost value that led to maximum likelihood
        [~, ii] = max(loglike_lst);
        paramBest2 = paramBest_mat(ii, :);
        paramBest_struct = pack_params_in_struct(paramBest2, params_set2, fixed_vals);
        paramBest = [paramBest1, paramBest2];
        loglike = loglike_lst(ii);


        % -> **************************************************************
        % -> fill in best fitted parameters (and print in command window)
        if (verbose_flag >= 1)
            fprintf('Best fit of full model:\n');
            
            my_params_cell = struct2cell(params_set);
            my_fld_nms = fieldnames(params_set);

            for my_kk = 1:param_free_nb
                my_fld_ind = cellfun(@(xx) find(xx == my_kk), my_params_cell, 'UniformOutput', false);
                my_uu = find(~cellfun(@isempty, my_fld_ind));

                % -> allow for multiple use of a variable across model parameters
                for my_pp = 1:length(my_uu)
                    my_vv = my_fld_ind{my_uu(my_pp)};

                    % -> allow for multiple use of a variable within model parameters
                    for my_qq = 1:length(my_vv)
                        fprintf('  Parameter #%d: %11s(%d) = %7.3f\n', my_kk, ...
                            my_fld_nms{my_uu(my_pp)}, my_vv(my_qq), paramBest(my_kk));
                    end
                end
            end
            
            fprintf('  Full Model Log-Likelihood = %7.3f\n', loglike);
        end
        
        conf_noise_task = paramBest_struct.conf_noise;
        conf_boost_task = paramBest_struct.conf_boost;
        conf_crit_task  = paramBest_struct.conf_crit;
        intrvl1_bias    = paramBest_struct.intrvl_bias;
        conf2_bias      = paramBest_struct.conf_bias;

        [chosen_full, nn1_full] = my_fun_full(paramBest2);
        
    end

    
    % -> stop timer
    if (verbose_flag >= 2)
        elapsedTime = toc(tstart);
        fprintf('\nElapsed time to perform fit = %7.3f sec\n', elapsedTime);
    end
    if (verbose_flag >= 1)
        fprintf('\n********************\n\n');
    end

        
    % -> fill in struct
    cfc_struct = struct;
    cfc_struct.sens_noise = sens_noise_task;    % sensory noise
    cfc_struct.sens_crit = sens_crit_task;      % sensory criterion
    cfc_struct.conf_noise = conf_noise_task;    % sdtev of noise for Type 2 decision (0 = ideal)
    cfc_struct.conf_boost = conf_boost_task;    % fraction super-ideal (1 - fraction ideal)
    cfc_struct.conf_crit = conf_crit_task;      % criterion for Type 2 decision (criterion)
    cfc_struct.intrvl_bias = intrvl1_bias;      % intrvl_bias (bias in favour of interval 1)
    cfc_struct.conf_bias = conf2_bias;          % conf_bias (bias in favour of task 2)
    
    cfc_struct.efficiency = efficiency;                % efficiency
    cfc_struct.equiv_conf_noise_ideal = noise2_idl;    % equivalent noise 2 for ideal perf
    cfc_struct.equiv_conf_noise_human = noise2_dat;    % equivalent noise 2 for human perf
    
    cfc_struct.choice_prob_ideal = chosen_idl;          % ideal choice probabilities
    cfc_struct.choice_prob_super_ideal = chosen_sdl;    % super-ideal choice probabilities
    cfc_struct.choice_prob_eff = chosen_eff;            % efficiency choice probabilities
    cfc_struct.choice_prob_model = chosen_full;         % full-model choice probabilities
    
    cfc_struct.choice_pair_ideal = nn1_idl;             % ideal choice responses
    cfc_struct.choice_pair_super_ideal = nn1_sdl;       % super-ideal choice responses
    cfc_struct.choice_pair_eff = nn1_eff;               % efficiency choice responses
    cfc_struct.choice_pair_model = nn1_full;            % full-model choice responses

    cfc_struct.loglike = loglike;               % log-likelihood of best fit

    

    
    % ---------------------------------------------------------------------
    % -> define nested functions
    
    % -> This is the Type 1 function for confidence model fit
    function yy_vals = cfc_type1(xx_vals, variable_prms, params_set, fixed_vals)

        noise_set_lst = params_set.sens_noise;
        crit_set_lst = params_set.sens_crit;
        tsk_lst = fixed_vals.tasks_list;
        noise_fix_lst = fixed_vals.sens_noise;
        crit_fix_lst = fixed_vals.sens_crit;
        nb_tasks = length(noise_set_lst);
        noise_val_lst = NaN(1, nb_tasks);
        crit_val_lst = NaN(1, nb_tasks);

        for t1_tt = 1:nb_tasks
            if (noise_set_lst(t1_tt))
                noise_val_lst(t1_tt) = variable_prms(noise_set_lst(t1_tt));
            else
                noise_val_lst(t1_tt) = noise_fix_lst(t1_tt);
            end

            if (crit_set_lst(t1_tt))
                crit_val_lst(t1_tt) = variable_prms(crit_set_lst(t1_tt));
            else
                crit_val_lst(t1_tt) = crit_fix_lst(t1_tt);
            end
        end


        nb_stims = size(xx_vals, 1);
        [t1_tsk_vals, ~, tsk_ic] = unique(xx_vals(:,1));
        nb_tsks = size(t1_tsk_vals, 1);

        yy_vals = NaN(nb_stims, 1);

        for t1_tt = 1:nb_tsks
            t1_inds = (tsk_ic == t1_tt);     % indices that have the same task

            xx = xx_vals(t1_inds, 2);

            my_tsk = t1_tsk_vals(t1_tt);
            tsk_ind = find(tsk_lst == my_tsk, 1);
            noise = noise_val_lst(tsk_ind);
            crit = crit_val_lst(tsk_ind);

            yy_vals(t1_inds) = basic_normcdf(xx, crit, noise);
        end
    end



    % -> wrapper around cfc_core function to flexibly add parameters
    function [pred_resp, pred_nn1] = cfc_core_wrap(wrap_data, variable_prms, params_set, fixed_vals)

        params_vals = pack_params_in_struct(variable_prms, params_set, fixed_vals);
        [pred_resp, pred_nn1] = cfc_core(wrap_data, params_vals);
    end


    function params_struct = pack_params_in_struct(variable_prms, params_set, fixed_vals)
    
        params_struct = fixed_vals;
        params_cell = struct2cell(params_set);
        fld_nms = fieldnames(params_set);
        
        variable_nb = length(variable_prms);

        for kk = 1:variable_nb
            fld_ind = cellfun(@(xx) find(xx == kk), params_cell, 'UniformOutput', false);
            uu = find(~cellfun(@isempty, fld_ind));
            
            % -> allow for multiple use of a variable across model parameters
            for pp = 1:length(uu)
                vv = fld_ind{uu(pp)};
                
                % -> allow for multiple use of a variable within model parameters
                for qq = 1:length(vv)
                    params_struct.(fld_nms{uu(pp)})(vv(qq)) = variable_prms(kk);
                end
            end
        end
    end


    % -> transform parameter structure into a vector of free parameters
    function params_vals = extract_params_from_struct(params_set, params_struct)
        
        params_cell = struct2cell(params_set);
        fld_nms = fieldnames(params_set);
        
        params_mat = cell2mat(params_cell');
        variable_nb = max(params_mat);
        
        params_vals = NaN(1, variable_nb);

        for kk = 1:variable_nb
            fld_ind = cellfun(@(xx) find(xx == kk), params_cell, 'UniformOutput', false);
            uu = find(~cellfun(@isempty, fld_ind));
            vv = fld_ind{uu(1)};
            params_vals(kk) = params_struct.(fld_nms{uu(1)})(vv(1));
        end
    end

    
    % ->   fit by maximizing the log-likelihood
    function [params_best, loglike] = fitnll(fit_fcn, nn1_lst, params_0, params_LB, params_UB, fit_options)

        fun = @(xx) loglikefcn(xx, fit_fcn, nn1_lst);
        [params_best, nll_best] = fminsearchbnd(fun, params_0, params_LB, params_UB, fit_options);
        loglike = - nll_best;
    end


    % ->   negative summed log-likelihood
    function nglglk = loglikefcn(pp, ff, nn)
        
        ypred = ff(pp);
        ypred(ypred == 0) = 1e-6;
        ypred(ypred == 1) = 1 - 1e-6;

        ll1 = log(ypred);
        ll0 = log(1.0 - ypred);

        ll_vect = nn(:,1) .* ll1 + nn(:,2) .* ll0;  % vector of log-likelihoods
        nglglk = - sum(ll_vect);      % minimize (-log) likelihood
    end


end

% -> THE END <-
% ------------------------------------------------------------------------
