% fPETregress - Perform ordinary least squares (OLS) or weighted least
% squares (WLS) regression on fPET-FDG data.
%
% Sean Coursey, 2024.11.03
% Jingyuan Chen's CANDY Lab
% The MGH/HST Martinos Center for Biomedical Imaging
%
% This function performs GLM regression on the input data, providing
% regression coefficients (B), residuals (epsilon), t-values (tvals),
% and percentage signal changes (perc_sig_change, if requested).
%
% The function can handle either standard design and contrast matrix inputs
% or cell-arrays to apply a different matrix/contrast to each voxel/ROI.
% It also offers several optional features through additional arguments:
%
%   - 'ortho': Orthogonalizes the design matrix (X) before regression.
%              (not recommended)
%
%   - 'robust': Computes robust standard errors using the HC3 method, 
%               which accounts for heteroscedasticity in the residuals.
%               (recommended)
%
%   - 'exclude': Excludes a specified number of timepoints from both
%                the design matrix (X) and the response variable (Y).
%
%   - 'signal change': Calculates percentage signal change based on
%                      Godbersen et al. 2024.
%
%   - 'WLS': Performs WLS instead of OLS regression. (Mutually exclusive
%            with 'robust'.) Requires the next varargin input after the
%            'WLS' flag to be a vector of variances for each timepoint, or
%            a matrix with a column for each voxel/ROI.
%
% Input arguments:
%   X_in: Design matrix (predictor variables). Can be either a matrix or a 
%         cell array of matrices.
%   Y_in: Response variable matrix.
%   C_in: Contrast matrix. Can be either matrix or struct with matrices 
%         C_in.c1, C_in.c2, &tc, for testing multiple contrasts
%         simultaneously.
%   varargin: Optional flags for orthogonalization, robust errors, exclusion, 
%             and signal change calculation.
%
% Output:
%   tvals: t-statistics for each condition's coefficients.
%   gammas: Coefficients for the specified conditions.
%   gamvars: Variance of the gamma coefficients.
%   B: Regression coefficients for the model.
%   epsilon: Residuals (errors) between predicted and actual values.
%   perc_sig_change: Percentage signal change (if 'signal change' option is enabled).

function [tvals, gammas, gamvars, B, epsilon, perc_sig_change] = fPETregress(X_in, Y_in, C_in, varargin)
    % Initialize flags for optional parameters
    ORTHO = false;   % If true, orthogonalize X before regression
    ROBUST = false;  % If true, use robust standard errors (HC3 method)
    WLS = false;     % If true, perform weighted least squares regression
    SIGNAL_CHANGE = false; % If true, calculate percentage signal change
    numArgsIn = length(varargin); % Number of optional arguments
    cell_bool = false; % Flag to check if X_in is a cell array

    % Parse optional input arguments
    for i = 1:numArgsIn
        if strcmpi(varargin{i}, 'ortho')
            ORTHO = true; % Set orthogonalize flag
        elseif strcmpi(varargin{i}, 'robust')
            ROBUST = true; % Set robust standard errors flag
        elseif strcmpi(varargin{i}, 'exclude')
            % Exclude data up to a certain index
            EXCLUDE = varargin{i+1};
            if ~cell_bool
                X_in = X_in((EXCLUDE+1):end, :); % Exclude rows from X
            else
                for k = 1:num_cell
                    X_in{k} = X_in{k}((EXCLUDE+1):end, :); % Exclude rows for each ROI
                end
            end
            Y_in = Y_in((EXCLUDE+1):end, :); % Exclude rows from Y
            i = i+1; % Skip next argument since it's the exclude index
        elseif strcmpi(varargin{i}, 'signal change')
            SIGNAL_CHANGE = true; % Set signal change flag
        elseif strcmpi(varargin{i}, 'WLS')
            WLS = true;
            variances_in = varargin{i+1};
            i = i+1;
        end
    end

    % Display warning if both WLS and ROBUST are true
    if WLS && ROBUST
        disp("Warning: Overriding ROBUST to false because WLS is true.")
        ROBUST = false;
    end

    % Check if X_in is a cell array
    if iscell(X_in)
        cell_bool = true;
        num_cell = length(X_in); % Number of cells (voxels/ROIs)
    end

    % If WLS, check if variances_in is a matrix, if it is, make X_in a cell array
    if WLS
        if size(variances_in, 2) > 1
            if ~cell_bool
                cell_bool = true;
                num_cell = size(variances_in, 2);
                X_in = repmat({X_in}, 1, num_cell);
            end
        end
    end

    % Orthogonalize X if needed
    if ORTHO
        if ~cell_bool
            X = orthogonalize(X_in);
        else
            X = cell(1,num_cell); % Initialize cell array for orthogonalized X
            for k = 1:num_cell
                X{k} = orthogonalize(X_in{k}); % Orthogonalize the design matrix for each ROI
            end
        end
    else
        X = X_in; % Use X_in as-is
    end
    
    % Initialize variables
    Y = Y_in;
    C = C_in;

    % Check if num_cell == number of ROIs
    if cell_bool
        assert(num_cell == size(Y_in, 2), "The number of cells does not match the number of ROIs.");
    end
    
    % If C is not a structure, convert it into one
    if ~isstruct(C)
        tmp = C;
        clear C
        C.c1 = tmp;
        clear tmp
    end
    

    if ~WLS
        % Perform ordinary least squares regression
        if ~cell_bool
            B = pinv(X)*Y; % Calculate beta coefficients (pseudo-inverse method)
            epsilon = Y - X*B; % Calculate residuals
        else
            B = cell(1,num_cell); % Initialize cell array for beta coefficients
            epsilon = zeros(size(Y)); % Initialize residuals matrix
            for k = 1:num_cell
                B{k} = pinv(X{k})*Y(:,k); % Calculate beta coefficients for each ROI
                epsilon(:,k) = Y(:,k) - X{k}*B{k}; % Calculate residuals for each ROI
            end
        end
    else
        % Perform weighted least squares regression
        if ~cell_bool
            % Calculate whitened X and Y matrices
            variances = variances_in(:); % Insure variances_in is a column vector
            X_w = X.*(variances.^(-0.5));
            Y_w = Y.*(variances.^(-0.5));
            % Calculate beta coefficients
            B = pinv(X_w)*Y_w;
            epsilon = Y - X*B; % Calculate residuals
            epsilon_w = Y_w - X_w*B; % Calculate whitened residuals
        else
            B = cell(1,num_cell); % Initialize cell array for beta coefficients
            epsilon = zeros(size(Y)); % Initialize residuals matrix
            epsilon_w = zeros(size(Y)); % Initialize whitened residuals matrix
            X_w = cell(1, num_cell); % Initialize whitened design matrix cell array
            for k = 1:num_cell
                % Calculate whitened X matrix and Y vector
                variances = variances_in(:,k); % Use variances from the proper ROI
                X_w{k} = X{k}.*(variances.^(-0.5));
                Y_w = Y(:,k).*(variances.^(-0.5));
                % Calculate beta coefficients for each ROI
                B{k} = pinv(X_w{k})*Y_w;
                epsilon(:,k) = Y(:,k) - X{k}*B{k}; % Calculate residuals for each ROI
                epsilon_w(:,k) = Y_w - X_w{k}*B{k}; % Calculate whitened residuals for each ROI
            end
        end
    end

    % Initialize percentage signal change array
    perc_sig_change = zeros(length(fieldnames(C)), size(Y, 2));
    
    % If signal change flag is set, calculate percentage signal change
    if SIGNAL_CHANGE
        for i = 1:length(fieldnames(C))
            act_slope = zeros(1, size(Y, 2)); % Initialize activation slope array
            weight = zeros(1, size(Y, 2)); % Initialize weight array
            if ~cell_bool
                tmpC = C.("c"+num2str(i)) == 1; % Identify the active condition
                for j = 1:length(tmpC)
                    if tmpC(j)
                        w = sum(diff(X(:,j)) == max(diff(X(:,j)))); % Calculate weight
                        act_slope = act_slope + w*max(diff(X(:,j)))*B(j, :); % Calculate activation slope
                        weight = weight + w; % Update weight
                    end                
                end
            else
                for k = 1:num_cell
                    tmpC = C.("c"+num2str(i)){k} == 1; % Identify active condition for each voxel/region
                    for j = 1:length(tmpC)
                        if tmpC(j)
                            w = sum(diff(X{k}(:,j)) == max(diff(X{k}(:,j)))); % Calculate weight
                            act_slope(k) = act_slope(k) + w*max(diff(X{k}(:,j)))*B{k}(j); % Calculate activation slope
                            weight(k) = weight(k) + w; % Update weight
                        end
                    end
                end
            end
            act_slope = act_slope ./ weight; % Normalize actual slope by weight
            lin_mdl = [ones(size(Y, 1), 1) (1:size(Y, 1))']; % Linear model
            fit_inds = (2*floor(size(Y, 1)/3):size(Y, 1))'; % Fit indices for baseline (last third of data)
            if ~cell_bool
                base_slope = [0 1]*(pinv(lin_mdl(fit_inds, :))*(X(fit_inds, ~tmpC)*B(~tmpC,:))); % Calculate baseline slope
            else
                base_slope = zeros(size(act_slope)); % Initialize baseline slope for each voxel/region
                for k = 1:num_cell
                    base_slope(k) = [0 1]*(pinv(lin_mdl(fit_inds, :))*(X{k}(fit_inds, ~(C.("c"+num2str(i)){k} == 1))*B{k}(~(C.("c"+num2str(i)){k} == 1)))); % Calculate baseline slope for each voxel/region
                end
            end
            perc_sig_change(i, :) = 100 * act_slope ./ base_slope; % Calculate percentage signal change
        end
    end

    % Initialize gamma and variance arrays
    gammas = zeros(length(fieldnames(C)), size(Y, 2));
    if ROBUST
        gamvars_robust = zeros(size(gammas)); % Initialize robust variance array

        if ~cell_bool
            all_gamvars = zeros(size(B)); % Initialize robust variance for all data
            parfor i = 1:size(epsilon, 2)
                all_gamvars(:,i) = diag(HC3(X,epsilon(:,i))); % Compute robust variances using HC3
            end
        else
            all_gamvars = cell(size(B)); % Initialize robust variance for each ROI
            for i = 1:size(epsilon, 2)
                all_gamvars{i} = diag(HC3(X{i},epsilon(:,i))); % Compute robust variances for each ROI
            end
        end
    else
        gamvars = zeros(size(gammas)); % Initialize non-robust variance array
    end

    % Calculate gamma coefficients and variances
    if ~cell_bool
        for i = 1:length(fieldnames(C))
            gammas(i,:) = C.("c"+num2str(i))*B; % Calculate gamma coefficients
            if ROBUST
                gamvars_robust(i,:) = C.("c"+num2str(i))*all_gamvars; % Calculate robust variances
            elseif WLS
                gamvars(i, :) = C.("c"+num2str(i))*((X_w')*X_w)^(-1)*(C.("c"+num2str(i))').*(sum(epsilon_w.^2, 1)./(size(Y, 1) - size(X, 2))); % Calculate WLS variances
            else
                gamvars(i, :) = C.("c"+num2str(i))*((X')*X)^(-1)*(C.("c"+num2str(i))').*(sum(epsilon.^2, 1)./(size(Y, 1) - size(X, 2))); % Calculate OLS variances
            end
        end
    else
        for i = 1:length(fieldnames(C))
            for k = 1:num_cell
                gammas(i, k) = C.("c"+num2str(i)){k}*B{k}; % Calculate gamma coefficients for each ROI
                if ROBUST
                    gamvars_robust(i, k) = C.("c"+num2str(i)){k}*all_gamvars{k}; % Calculate robust variances for each ROI
                elseif WLS
                    gamvars(i, k) = C.("c"+num2str(i))*((X_w{k}')*X_W{k})^(-1)*(C.("c"+num2str(i))').*(sum(epsilon_w(:,k).^2, 1)./(size(Y, 1) - size(X{k}, 2))); % Calculate WLS variances for each ROI
                else
                    gamvars(i, k) = C.("c"+num2str(i))*((X{k}')*X{k})^(-1)*(C.("c"+num2str(i))').*(sum(epsilon(:,k).^2, 1)./(size(Y, 1) - size(X{k}, 2))); % Calculate OLS variances for each ROI
                end
            end
        end
    end

    % Calculate t-values using either robust or non-robust variances
    if ROBUST
        tvals = gammas ./ sqrt(gamvars_robust);
        gamvars = gamvars_robust;
    else
        tvals = gammas ./ sqrt(gamvars);
    end
end

function covarianceMtrx = HC3(X,e)
    n = numel(e);
    e = reshape(e, n, 1);
    X = reshape(X, n, numel(X)/n);

    H = X * (X'*X)^-1 * X';
    h = diag(H);

    u = (e./(1 - h)).^2;

    covarianceMtrx = (X'*X)^-1 * X'*diag(u)*X * (X'*X)^-1;
end