% CONFIDENCE FORCED-CHOICE TOOLBOX  v0.3
%
% cfc_core
%   This function is the core function for confidence forced-choice
%   predictions. It computes the probability of choosing interval 1 for a 
%   given trial of a confidence forced-choice experiment.
%
% INPUT:
%   'grouped_data': grouped data per (s1, s2, r1, r2):
%       1st col: stimulus intensity interval 1
%       2nd col: stimulus intensity interval 2
%       3rd col: perceptual decision interval 1
%       4th col: perceptual decision interval 2
%       5th col: nb of confidence choices for interval 1
%       6th col: nb of confidence choices for interval 2
%       7th col: stimulus task for interval 1 (optional)
%       8th col: stimulus task for interval 2 (optional)
%
%   'model_params_vals': model parameters, as a structure (use vectors if more than 1 task):
%       'tasks_list' : vector of tasks, e.g. '1' or '[1, 2]'
%       'sens_noise' : sensory (Type 1) sdtev of noise (0 = perfectly sensitive)
%       'sens_crit'  : sensory (Type 1) criterion
%       'conf_noise' : confidence (Type 2) sdtev of noise (0 = ideal)
%       'conf_boost' : fraction super-ideal (0 = ideal, 1 = super-ideal)
%       'conf_crit'  : confidence (Type 2) criterion
%       'intrvl_bias': bias in favour of interval 2 (over interval 1)
%       'conf_bias'  : overconfidence relative to one of the tasks
%
%
% OUTPUT:
%   conf_choice_prob: Type 2 choice probability to select interval 1 (rather than 2)
%
%   conf_choice_freq: nb of confidence choices for intervals 1 and 2 provided so
%       that their sum across the two intervals matches the sum in
%       'grouped_data' (columns 5 + 6)
%
%   choose2_resp: P(ch=1 | r1, r2)
%
%
% 20-SEP-2020 - pascal mamassian

% 05-NOV-2017 - pascal mamassian
% 24-APR-2018 - pm: split integral space into subspaces to place discontinuities at boundaries
% 29-DEC-2018 - pm: take into account type 2 criterion
% 06-JAN-2019 - pm: computes all 4 type 2 choice probabilities
% 07-FEB-2019 - pm: expand to all stimuli of an experiment
% 10-FEB-2019 - pm: introduce 'small triangle' for type 2 criterion
% 15-FEB-2019 - pm: generalize for different sensory tasks
% 13-MAR-2019 - pm: normalize confidence evidence by sensory sensitivity
% 21-MAR-2019 - pm: do not fit second noise2 parameter (undetermined)
% 29-MAR-2019 - pm: take care of extreme probability values
% 17-AUG-2019 - pm: fixed type 2 criterion
% 21-AUG-2019 - pm: cleaned up
% 22-SEP-2019 - pm: updated interval bias
% 03-FEB-2020 - pm: changed definition of confidence criterion
% 10-JUN-2020 - pm: fast normcdf from Tarryn
% 01-JUL-2020 - pm: changed format of model_params
% 30-AUG-2020 - pm: updated decision rule for interval bias
% 20-SEP-2020 - pm: cleaned up


function [conf_choice_prob, conf_choice_freq, choose2_resp] = ...
    cfc_core(grouped_data, model_params_vals)


% -> nb. of cases of (s1, s2, r1, r2, t1, t2)
nb_knds = size(grouped_data, 1);
conf_choice_prob = NaN(nb_knds, 1);
choose2_resp = NaN(nb_knds, 1);

if (size(grouped_data, 2) < 7)
    % by default, assume that there is only one sensory task
    tasks_nn = repmat([1, 1], nb_knds, 1);
    grouped_data = [grouped_data, tasks_nn];
end


% -> tolerances for Matlab integral functions (smaller is better)
atol = 1e-6;
rtol = 1e-2;

% -> display on-going information
% verbose_on = 1;
verbose_on = 0;


% -> Type 1 parameters
tasks_list      = model_params_vals.tasks_list;
sens_noise_task = model_params_vals.sens_noise;  % sensory noise, per task
sens_crit_task  = model_params_vals.sens_crit;   % sensory criterion, per task
nb_tasks = length(tasks_list);

% -> Type 2 parameters (compulsory)
conf_noise_task = model_params_vals.conf_noise;  % confidence noise, per task
conf_boost_task = model_params_vals.conf_boost;  % confidence boost, per task

% -> Type 2 parameters (optional)
% -> confidence criterion, per task
if (any(strcmp(fieldnames(model_params_vals), 'conf_crit')))
    conf_crit_task  = model_params_vals.conf_crit;   
else
    conf_crit_task = zeros(1, nb_tasks);
end

% -> interval bias, in favour of interval 1
if (any(strcmp(fieldnames(model_params_vals), 'intrvl_bias')))
    intrvl_bias = model_params_vals.intrvl_bias;
else
    intrvl_bias = 0.0;
end

% -> confidence bias: overconfidence (>1) relative to one of the tasks (=1)
if (any(strcmp(fieldnames(model_params_vals), 'conf_bias')))
    conf_scale = model_params_vals.conf_bias;
else
    conf_scale = ones(1, nb_tasks);
end


% -> rescale by sensory sensitivity
conf_scale = conf_scale ./ sens_noise_task;

% -> constrain (2noise2 = 1noise2) if there are 2 tasks but 2noise2 is NaN
noise2_inds = isnan(conf_noise_task);
conf_noise_task(noise2_inds) = conf_noise_task(1);

% -> avoid delta functions in integrals due to zero confidence noise
conf_noise_task(conf_noise_task == 0) = 1e-6;

% -> approximate infinity by restricting bounds in integrals, per task
intsz_task = 10.0 .* sens_noise_task;


% -> go through all the cases of (s1, s2, r1, r2, t1, t2)
for kk = 1:nb_knds

    % -> Type 1 tasks (usually '1' or '2', but can be any integer)
    tsk1 = grouped_data(kk, 7);    % task for interval 1
    tsk2 = grouped_data(kk, 8);    % task for interval 2
    
    % -> match task number to the list of tasks to extract relevant params
    [~, tsk1_ind] = ismember(tsk1, tasks_list);  % task index for intrvl 1
    [~, tsk2_ind] = ismember(tsk2, tasks_list);  % task index for intrvl 2
    
    % -> Stimulus intensities
    mu1 = grouped_data(kk, 1);
    x1_min = mu1 - intsz_task(tsk1_ind);
    x1_max = mu1 + intsz_task(tsk1_ind);

    mu2 = grouped_data(kk, 2);
    x2_min = mu2 - intsz_task(tsk2_ind);
    x2_max = mu2 + intsz_task(tsk2_ind);

    % -> Type 1 responses (1 = 'A', 0 = 'B')
    resp1 = grouped_data(kk, 3);
    resp2 = grouped_data(kk, 4);

    % -> parameters for 1st interval
    sens_noise1 = sens_noise_task(tsk1_ind);
    sens_crit1 = sens_crit_task(tsk1_ind);
    conf_noise1 = conf_noise_task(tsk1_ind);
    conf_crit1 = conf_crit_task(tsk1_ind);
    conf_boost1 = conf_boost_task(tsk1_ind);
    conf_scale1 = conf_scale(tsk1_ind);
    
    % -> parameters for 2nd interval
    sens_noise2 = sens_noise_task(tsk2_ind);
    sens_crit2 = sens_crit_task(tsk2_ind);
    conf_noise2 = conf_noise_task(tsk2_ind);
    conf_crit2 = conf_crit_task(tsk2_ind);
    conf_boost2 = conf_boost_task(tsk2_ind);
    conf_scale2 = conf_scale(tsk2_ind);
    
    
    % -> probability density function of joint sensory evidence:
    %    P(s1, s2)
    noisy1_hdl = @(xx, yy) ...
        noisy1_fun(xx, yy, mu1, mu2, sens_noise1, sens_noise2);


    % -> probability to choose interval 1 given a pair of (smpl1, smpl2):
    %    P(ch=1 | s1, s2)
    choice2_hdl = @(xx, yy) ...
        choice2_fun(xx, yy, mu1, mu2, ...
                sens_crit1, conf_noise1, conf_crit1, conf_boost1, conf_scale1, ...
                sens_crit2, conf_noise2, conf_crit2, conf_boost2, conf_scale2, ...
                intrvl_bias);

    % -> probability to choose interval 1 for a pair of (smpl1, smpl2):
    %    P(ch=1, s1, s2)
    choose_intrvl_1_hdl = @(xx, yy) ...
        noisy1_hdl(xx, yy) .* choice2_hdl(xx, yy);



    % -> go through all 4 possible scenarii of response pairs
    if (resp1 == 1)
        if (resp2 == 1)
            % -> Type 1 responses were 'A' and 'A'
            % -> P(r1=1, r2=1)
            joint_prob = (1 - basic_normcdf(sens_crit1, mu1, sens_noise1)) * ...
                (1 - basic_normcdf(sens_crit2, mu2, sens_noise2));
            % -> P(ch=1 | r1=1, r2=1)
            choose_intrvl1 = integral2(choose_intrvl_1_hdl, sens_crit1, x1_max, ...
                sens_crit2, x2_max,'AbsTol',atol,'RelTol',rtol);

        else
            % -> Type 1 responses were 'A' and 'B'
            joint_prob = (1 - basic_normcdf(sens_crit1, mu1, sens_noise1)) * ...
                basic_normcdf(sens_crit2, mu2, sens_noise2);
            choose_intrvl1 = integral2(choose_intrvl_1_hdl, sens_crit1, x1_max, ...
                x2_min, sens_crit2,'AbsTol',atol,'RelTol',rtol);
            
        end
    else
        
        if (resp2 == 1)
            % -> Type 1 responses were 'B' and 'A'
            joint_prob = basic_normcdf(sens_crit1, mu1, sens_noise1) * ...
                (1 - basic_normcdf(sens_crit2, mu2, sens_noise2));
            choose_intrvl1 = integral2(choose_intrvl_1_hdl, x1_min, sens_crit1, ...
                sens_crit2, x2_max,'AbsTol',atol,'RelTol',rtol);
            
        else
            % -> Type 1 responses were 'B' and 'B'
            joint_prob = basic_normcdf(sens_crit1, mu1, sens_noise1) * ...
                basic_normcdf(sens_crit2, mu2, sens_noise2);
            choose_intrvl1 = integral2(choose_intrvl_1_hdl, x1_min, sens_crit1, ...
                x2_min, sens_crit2,'AbsTol',atol,'RelTol',rtol);
            
        end
    end

    % -> probability to choose interval 1
    conf_choice_prob(kk) = choose_intrvl1 / joint_prob;
    
    choose2_resp(kk) = choose_intrvl1;

    if (verbose_on)
        if (mod(kk, 50) == 0)
            fprintf('.\n');
        elseif (mod(kk, 10) == 0)
            fprintf('. ');
        else
            fprintf('.');
        end
    end

end

% -> take care of very small joint_prob that are rounded off to 0.0
conf_choice_prob(isnan(conf_choice_prob)) = 0.5;
conf_choice_prob(conf_choice_prob < 0.0) = 0.0;
conf_choice_prob(conf_choice_prob > 1.0) = 1.0;

if (verbose_on)
    fprintf('\n');
end

nn1_choices = grouped_data(:, 5:6);
nn1_tmp1 = sum(nn1_choices, 2);
nn1_tmp2 = round(nn1_tmp1 .* conf_choice_prob);
conf_choice_freq = [nn1_tmp2, nn1_tmp1 - nn1_tmp2];


% ------------------------------------------------------------------------
% -> internal functions for model

% -> probability density function of (smpl1, smpl2): 
%    P(s1, s2)
function prob = noisy1_fun(smpl1, smpl2, mu1, mu2, sd1, sd2) 
    prob = normpdf(smpl1, mu1, sd1) .* normpdf(smpl2, mu2, sd2);
end


% -> probability to choose interval 1 for a pair of (smpl1, smpl2):
%    P(ch=1 | s1, s2)
function prob = choice2_fun(smpl1, smpl2, mu1, mu2, ...
        s_crit1, c_noise1, c_crit1, c_boost1, c_beta1, ...
        s_crit2, c_noise2, c_crit2, c_boost2, c_beta2, c_gamma) 
    
    % -> this function is often called within a 2D integral, i.e. with
    %    a matrix of samples, so extract dimension of matrix
    [smplx_nb, smply_nb] = size(smpl1);
    prob = NaN(smplx_nb, smply_nb);
                
    for smplx_ii=1:smplx_nb
        for smply_ii=1:smply_nb
            smpl_val1 = smpl1(smplx_ii, smply_ii);
            smpl_val2 = smpl2(smplx_ii, smply_ii);
            
            % -> perceptual decision (1 = 'A', 0 = 'B') for intervals 1 & 2
            intrvl1_dec = (smpl_val1 > s_crit1);
            intrvl2_dec = (smpl_val2 > s_crit2);

            % -> sliced mean for intervals 1 and 2
            c_Ox = (smpl_val1 + (mu1 - smpl_val1)*c_boost1 - s_crit1 - c_crit1) * c_beta1;
            c_Oy = (smpl_val2 + (mu2 - smpl_val2)*c_boost2 - s_crit2 - c_crit2) * c_beta2;

            % -> rescale the axes so that the confidence noise distribution
            %    is bivariate Normal with unit variance (and no covariance)
            c_Oxr = c_Ox / c_noise1;
            c_Oyr = c_Oy / c_noise2;

            % -> rescale interval bias
            c_g = c_gamma / sqrt(c_noise1*c_noise1 + c_noise2*c_noise2);

            % -> rotate coordinates Oxy clockwise to Ouv so that Ou
            %    is parallel to the separating line for AB and BA
            slope_yx1 = - c_noise1 / c_noise2;     % slope of that line
            c_Ov1 = cfc_rotate2(c_Oxr, c_Oyr, slope_yx1);

            % -> rotate coordinates Oxy counter-clockwise to Ou'v' so that Ou'
            %    is parallel to the separating line for AA and BB
            slope_yx2 = c_noise1 / c_noise2;     % slope of that line
            c_Ov2 = cfc_rotate2(c_Oxr, c_Oyr, slope_yx2);

            % -> go through all 4 possible scenarii
            if (intrvl1_dec == 1)
                if (intrvl2_dec == 1)
                    % -> Type 1 responses were 'A' and 'A'
                    half_space = basic_normcdf(c_g, c_Ov2, 1.0);
                else
                    % -> Type 1 responses were 'A' and 'B'
                    half_space = 1 - basic_normcdf(-c_g, c_Ov1, 1.0);
                end
            else
                if (intrvl2_dec == 1)
                    % -> Type 1 responses were 'B' and 'A'
                    half_space = basic_normcdf(c_g, c_Ov1, 1.0);
                else
                    % -> Type 1 responses were 'B' and 'B'
                    half_space = 1 - basic_normcdf(-c_g, c_Ov2, 1.0);
                end
            end

            prob(smplx_ii, smply_ii) = half_space;
        end
    end
    

end


% -> rotate coordinates Oxy clockwise to Ouv so that Ou
%    is parallel to the line that has slope 'aa'
function mv = cfc_rotate2(mx, my, aa)
    cc = sqrt(1 + aa^2);
    mv = (-aa* mx + my) / cc;
end


end

% -> THE END <-
% ------------------------------------------------------------------------