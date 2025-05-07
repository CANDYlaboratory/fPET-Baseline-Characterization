% This script is meant to demonstrate the use of the various functions in
% this github repository.

%% Load data, simulate TAC, and plot

% Load example resting-state constant-infusion fPET-FDG data
dataset = load("resting_CI_sample_subjects.mat");
subjects = ["sub0" + num2str((1:5)')]'; % define subject names
n_subs = 5;

% Set simulation parameters
Tmax = 338; % total time points (90 minutes given 16-second timestep)
CI_AIF = load("cont_infusion_group_AIF.mat").cont_infusion_group_AIF; % arterial input function

% Simulate a time activity curve using 2-compartment model
simTAC = fPETsimulateTAC([0.102 0.13 0.062 0.0068], CI_AIF, Tmax, 16/60);

% Define time in minutes
time = ((1:Tmax)*16/60 - 8/60)'; % subtracting 8 seconds places the timepoint in the middle of the PET frame

% Plot each subjects' mean TAC and calculate the mean across subjects

mean_TAC = zeros(Tmax, 1);

figure;
hold on;

for sub = subjects
    plot(time, mean(dataset.(sub)(1:Tmax,:), 2), 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5])
    mean_TAC = mean_TAC + mean(dataset.(sub)(1:Tmax,:), 2);
end

mean_TAC = mean_TAC/n_subs;

plot(time, mean_TAC, 'r', 'LineWidth', 2);

hold off;

% Scale the simulated TAC to the data and plot it
scale = pinv(simTAC)*mean_TAC;
simTAC = simTAC*scale;

hold on;
plot(time, simTAC, 'b--', 'LineWidth', 2);
hold off;

xlabel("Time (minutes)");
xlim([0 90]);

ylabel("Activity (kBq/cc)");
ylim([0 55]);



%% Detrend data

figure;
hold on;

% Initialize an array to store detrended data for each subject
detrended_data = zeros(Tmax,n_subs);

% Loop through each subject to detrend their data
for isub = 1:n_subs
    sub = subjects(isub);
    tmp_data = dataset.(sub)(1:Tmax,:); % truncate data to 90 minutes

    % Generate baseline model
    baseline = fPETgetBaseline(tmp_data, "EXP2", [], 16/60);

    % Detrend using baseline model
    tmp_data_detrended = tmp_data - baseline*(pinv(baseline)*tmp_data);

    % Plot mean detrended TAC for each subject
    plot(time, mean(tmp_data_detrended, 2), 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5])

    % Save subject's mean detrended TAC
    detrended_data(:,isub) = mean(tmp_data_detrended, 2);
end

% Plot mean detrended data across all subjects
plot(time, mean(detrended_data, 2), 'LineWidth', 2, 'Color', 'r');

% Detrend and plot simulated TAC
baseline_sim = fPETgetBaseline(simTAC, "EXP2", [], 16/60);
simTAC_detrended = simTAC - baseline_sim*(pinv(baseline_sim)*simTAC);
plot(time, simTAC_detrended, 'LineWidth', 2, 'Color', 'b')

%% GLM Analaysis of Data

% Create a sham task regressor to test on the resting-state data
sham_task_minutes = [zeros(20, 1); (1:10)'; 10*ones(10, 1); (11:20)'; 20*ones(10,1); (21:30)'; 30*ones(10,1); (31:40)'; 40*ones(10, 1)]; % ramp task defined in minutes
sham_task = interp1((1:100)*60 - 30, sham_task_minutes', (1:338)*16 - 8); % translated to the 16-second bins of the data
sham_task = sham_task(:); % make sure it is a column vector

% Exclude the first 10 minutes from GLM fitting
n_exclude = round(10*60/16);

% Initialize arrays for storing GLM results across subjects
gammas = zeros(100, n_subs);
t_scores = zeros(100, n_subs);
perc_sig_change = zeros(100, n_subs);

% Run GLM analysis for each subject
for isub = 1:n_subs
    sub = subjects(isub);
    tmp_data = dataset.(sub)(1:Tmax,:);

    % Generate baseline model
    baseline = fPETgetBaseline(tmp_data, "EXP2", sham_task, 16/60, 'exclude', n_exclude);

    % Create design and contrast matrices
    design_matrix = [baseline sham_task];
    contrast_matrix = [zeros(1,size(baseline,2)) ones(1, size(sham_task, 2))];

    % Perform GLM regression and save results
    [t_scores(:,isub), gammas(:,isub), ~, ~, ~, perc_sig_change(:,isub)] = fPETregress(design_matrix, tmp_data, contrast_matrix, 'exclude', n_exclude, 'robust', 'signal change');
end

% Calculate group-level random-effect t-scores
rfx_t_scores = mean(gammas, 2)./(std(gammas, [], 2)/sqrt(n_subs));

% Return index of any p < 0.05 significant (de)activation, FDR corrected using the Benjamini-Hochburg procedure
significant_index = benjaminiHochburg(rfx_t_scores, 0.05, 'tvals', n_subs); % no significance with only 5 subjects, this is just an example of how to use the function