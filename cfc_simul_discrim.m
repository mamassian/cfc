% CONFIDENCE FORCED-CHOICE TOOLBOX  v0.2
%
% cfc_simul_discrim
%   This function generates simulation of data of confidence forced-choice
%   for a Type 1 discrimination task.
%
%
% INPUT:
%   'simul_params': structure of simulation parameters:
%                   (vector of struct if number of unique tasks > 1)
%       simul_params(task_no).sens_intens: difficulty levels (constant stimuli)
%       simul_params(task_no).sens_intens_min: minimum of range for stimuli (uniform sampling)
%       simul_params(task_no).sens_intens_max: maximum of range for stimuli (uniform sampling)
%       simul_params(task_no).method: use of method of constant stimuli ('1' or '0')
%       simul_params(1).nb_trials: number of confidence pairs to simulate
%
%   'model_params': model values parameters, as a structure
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
% OUTPUT:
%   'raw_data': matrix of:
%       1st col: stimulus intensity interval 1
%       2nd col: stimulus intensity interval 2
%       3rd col: perceptual decision interval 1 (1 = 'A',  0 = 'B')
%       4th col: perceptual decision interval 2 (1 = 'A',  0 = 'B')
%       5th col: confidence choice for interval 1 (1 = chosen,  0 = declined)
%       6th col: confidence choice for interval 2 (1 = chosen,  0 = declined)
%       7th col: stimulus task for interval 1
%       8th col: stimulus task for interval 2
% 
%
% 06-MAR-2020 - pascal mamassian
% 30-AUG-2020 - pm: updated decision rule for interval bias
% 27-SEP-2020 - pm: changed format of simul_params and model_params


function raw_data = cfc_simul_discrim(simul_params, model_params)


% *************************************************************************
% -> PARAMETERS <-

% -> Type 1 parameters
tasks_list      = model_params.tasks_list;
sens_noise_task = model_params.sens_noise;  % sensory noise, per task
sens_crit_task  = model_params.sens_crit;   % sensory criterion, per task
nb_tasks = length(tasks_list);

% -> Type 2 parameters (compulsory)
conf_noise_task = model_params.conf_noise;  % confidence noise, per task
conf_boost_task = model_params.conf_boost;  % confidence boost, per task

% -> Type 2 parameters (optional)
% -> confidence criterion, per task
if (any(strcmp(fieldnames(model_params), 'conf_crit')))
    conf_crit_task  = model_params.conf_crit;   
else
    conf_crit_task = zeros(1, nb_tasks);
end

% -> interval bias, in favour of interval 1
if (any(strcmp(fieldnames(model_params), 'intrvl_bias')))
    intrvl_bias = model_params.intrvl_bias;
else
    intrvl_bias = 0.0;
end

% -> confidence bias: overconfidence (>1) relative to one of the tasks (=1)
if (any(strcmp(fieldnames(model_params), 'conf_bias')))
    conf_scale = model_params.conf_bias;
else
    conf_scale = ones(1, nb_tasks);
end


% -> these are actually (confidence) pairs of trials
nb_trials = simul_params(1).nb_trials;


% *************************************************************************

tasks_nn = NaN(nb_trials, 2);   % stimulus tasks for intervals 1, 2
stims_nn = NaN(nb_trials, 2);   % stimulus intensities for intervals 1, 2
percs_nn = NaN(nb_trials, 2);   % percepts for intervals 1, 2
choic_nn = NaN(nb_trials, 2);   % confidence choice: [1, 0] if interval 1, [0, 1] if interval 2
sens_smpl = NaN(1, 2);  % the 2 sensory samples in a confidence pair used for Type 1
conf_evd = NaN(1, 2);       % confidence evidence
conf_mag = NaN(1, 2);       % confidence magnitude
resp1 = NaN(1, 2);  % decision based on sensory evidence (Type 1)
resp2 = NaN(1, 2);  % decision based on confidence evidence (Type 2)


for tt = 1:nb_trials
    
    % -> interval kk of confidence pair
    for intrv = 1:2
        
        % -> choose randomly one of the tasks
        ind_task = randi(nb_tasks);  % index of task for interval kk

        tasks_nn(tt, intrv) = tasks_list(ind_task); % which task

        % -> choose stimulus and make response

        % -> choose randomly a stimulus in the set of possible stimuli
        if (simul_params(ind_task).method)
            sens_intens = simul_params(ind_task).sens_intens;   % range of stimuli for that task
            nb_sens_intens = length(sens_intens);   % nb of stimuli for that task

            ind_intens = randi(nb_sens_intens);   % index of stimulus for interval 'intrv'
            stims_nn(tt, intrv) = sens_intens(ind_intens);  % (noise-less) stimulus intensity
        else
            sens_intens_min = simul_params(task_no).sens_intens_min;
            sens_intens_max = simul_params(task_no).sens_intens_max;

            stim_val = rand*(sens_intens_max - sens_intens_min) + sens_intens_min;
            stims_nn(tt, intrv) = stim_val;  % (noise-less) stimulus intensity
        end

        % -> take independent noisy samples of the stimuli
        sens_smpl(intrv) = stims_nn(tt, intrv) + randn * sens_noise_task(ind_task);

        % -> sensory decision based on side of sensory criterion
        resp1(intrv) = (sens_smpl(intrv) > sens_crit_task(ind_task));
        percs_nn(tt, intrv) = resp1(intrv);

        % -> add fraction super-ideal
        boosted_sens_smpl = (1 - conf_boost_task(ind_task))*sens_smpl(intrv) + ...
            conf_boost_task(ind_task)*stims_nn(tt, intrv);

        % -> confidence sample
        conf_smpl = boosted_sens_smpl - sens_crit_task(ind_task) - conf_crit_task(ind_task);

        % -> scale to sensory sensitivity
        scaled_smpl = conf_smpl / sens_noise_task(ind_task) * conf_scale(ind_task);

        % -> confidence evidence: corrupt by Type 2 noise
        conf_evd(intrv) = scaled_smpl + randn * conf_noise_task(ind_task);

        % -> confidence magnitude
        conf_mag(intrv) = abs(conf_evd(intrv));

        % -> sensory decision based on confidence sample
        resp2(intrv) = (conf_evd(intrv) > 0);    

        % -> if pseudo-perceptual decision is inconsistent with
        %    perceptual decision, that perceptual decision might well be wrong
        if (resp2(intrv) ~= resp1(intrv))
            conf_mag(intrv) = - conf_mag(intrv);
        end

    end

    
    % -> confidence decision rule, adding interval bias
    pref_intrvl1 = conf_mag(1) - conf_mag(2) + intrvl_bias;
    
    % -> are you more confident on 1st interval?
    choice = (pref_intrvl1 > 0);   % 1 if intrv1, 0 if intrv2
    choice_intrvl = 2 - choice;    % 1 if intrv1, 2 if intrv2
        
    choic_nn(tt, choice_intrvl) = 1;
    choic_nn(tt, 3 - choice_intrvl) = 0;
    
end

% -> group data by identical stimuli and responses
raw_data = [stims_nn, percs_nn, choic_nn, tasks_nn];

end

% -> THE END
