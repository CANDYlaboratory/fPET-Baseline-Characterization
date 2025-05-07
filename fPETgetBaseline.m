% fPETgetBaseline - Computes the baseline for functional PET data based on the 
% specified detrending type.
%
% Sean Coursey, 2024.11.03
% Jingyuan Chen's CANDY Lab
% The MGH/HST Martinos Center for Biomedical Imaging
%
%   baseline = fPETgetBaseline(data_in, type, nuisance, varargin) returns 
%   a baseline signal, where the baseline model varies based on the specified 
%   'type' of detrending (e.g., polynomial fit, GSR, biexponential, etc.)
%
% Parameters:
%   data_in  - Matrix of time-series data (time points x number of signals)
%   type     - String indicating the detrending method to apply. Options include:
%              'Exp2', 'SA', 'P(n)MT' (where n is an integer), 'GSR', 'Poly(n)', etc.
%   nuisance - Matrix of nuisance variables to include in the baseline fit.
%   varargin - Additional parameters depending on the specified detrending type:
%              - 'dt' (time step size)
%              - 'AIF' (arterial input function, required for 'SA' and 'LocalSA' types)
%              - 'b0_off' (optional, excludes the intercept term if true)
%              - 'plot' (optional, if true enables debug plotting)
%              - 'exclude' (time index to exclude data before)
%              - 'lin_fit_start' (start time for linear fitting, in 'Exp2' type)
%
% Returns:
%   baseline - Baseline model (matrix or cell array depending on the 'type')
    
function baseline = fPETgetBaseline(data_in, type, nuisance, varargin)
    
    % Set default flags and parameters
    b0_bool = true;
    plot_bool = false;
    start = 1;
    T = size(data_in, 1);
    fit_idx = ones(T,1) == 1;

    % Parse parameters specific to certain detrending types
    if ~isempty(regexp(type, "[Ee][Xx][Pp]2", 'once'))
        dt = varargin{1};
        lin_fit_start = 45;
        start = 2;
    elseif ~isempty(regexp(type, "[Ss][Aa]", 'once'))
        dt = varargin{1};
        AIF = varargin{2};
        start = 3;
    end

    % Parse optional arguments in varargin
    for i = start:length(varargin)
        if strcmpi(varargin{i}, 'b0_off')
            b0_bool = false;
        elseif strcmpi(varargin{i}, 'plot')
            plot_bool = varargin{i+1}; i = i+1;
        elseif strcmpi(varargin{i}, 'exclude')
            fit_idx = (1:T) > varargin{i+1}; i = i+1;
        elseif strcmpi(varargin{i}, 'lin_fit_start')
            lin_fit_start = varargin{i+1}; i = i+1;
        end
    end

    % Generate baseline model given detrending type
    if ~isempty(regexp(type, "[Pp]\d+MT", 'once'))
        % Polynomial approximation of mean TAC with nuisance
        n = str2num(regexp(type, "\d*", 'match', 'once'));
        polyn = repmat((1:T)', 1, n+1);
        for i = 0:n
            polyn(:,i+1) = polyn(:,i+1).^i;
        end
        B = pinv([polyn(fit_idx,:) nuisance(fit_idx,:)])*mean(data_in(fit_idx, :), 2);
        baseline = polyn*B(1:size(polyn, 2));
        if b0_bool
            baseline = [ones(T, 1) baseline];
        end

    elseif ~isempty(regexp(type, "^[Pp]\d+", 'once'))
        % Simple polynomial detrending
        n = str2num(regexp(type, "\d*", 'match', 'once'));
        baseline = repmat((1:T)', 1, n+1);
        for i = 0:n
            baseline(:,i+1) = baseline(:,i+1).^i;
        end
        if ~b0_bool
            baseline = baseline(:,2:end);
        end
        
        if plot_bool
            mt = mean(data_in, 2);
            figure;
            subplot(2,1,1); plot(mt)
            hold on;
            plot(baseline*(pinv(baseline)*mt))
            hold off;
            legend(["Mean TAC" "Poly3 Fit"])
            subplot(2,1,2); plot(mt - baseline*(pinv(baseline)*mt))
        end

    elseif ~isempty(regexp(type, "^[Mm][Tt]", 'once'))
        % Mean TAC detrending
        baseline = mean(data_in, 2);
        if b0_bool
            baseline = [ones(size(data_in, 1), 1) baseline];
        end

    elseif ~isempty(regexp(type, "[Ee][Xx][Pp]2", 'once'))
        % Linear plus biexponential detrending
        mt = mean(data_in, 2);
        
        if plot_bool
            figure;
            ax1 = subplot(3, 1, 1); plot(mt)
        end

        t = dt*(1:T)';
        lineMdl = [ones(T,1) t];
        initLineFit = pinv(lineMdl(floor(lin_fit_start/dt):T, :))*mt(floor(lin_fit_start/dt):T);
        
        if plot_bool
            hold on;
            plot(lineMdl*initLineFit);
            hold off;
            title("Mean TAC")
        end

        mtLinRem = mt - lineMdl(:,2)*initLineFit(2);
        
        [a, b, c, p, q] = fit_biexp(t, mtLinRem);
        
        if plot_bool
            subplot(3, 1, 2); plot(mtLinRem);
            hold on;
            plot(a + b*exp(p*t) + c*exp(q*t));
            hold off;
            title("Linear Trend Removed, p="+num2str(p)+", q="+num2str(q))
        end

        if p ~= q
            baseline = [ones(T, 1) t exp(p*t) exp(q*t)];
        else
            baseline = [ones(T, 1) t exp(p*t)];
        end
        
        if ~b0_bool
            baseline = baseline(:,2:end);
        end
        
        if plot_bool
            axes(ax1);
            hold on;
            plot(baseline*(pinv(baseline)*mt));
            hold off;
            title("Mean TAC and model fit")
    
            subplot(3,1,3); plot(mt - baseline*(pinv(baseline)*mt));
            title("Residuals")
        end
        
    elseif ~isempty(regexp(type, "^[Ss][Aa]", 'once'))
        % Spectral analysis detrending with nuisance
        mt = mean(data_in, 2);
        baseline = spectral_baseline(mt, AIF, nuisance, T, dt);
        if b0_bool
            baseline = [ones(T,1) baseline];
        end

        if plot_bool
            figure;
            subplot(3,1,1);
            plot(mt); hold on; plot(baseline*(pinv(baseline)*mt(1:T))); hold off;
            subplot(3,1,2);
            plot(baseline);
            subplot(3,1,3);
            plot(mt(1:T) - baseline*(pinv(baseline)*mt(1:T)));
        end

    elseif ~isempty(regexp(type, "^[Ll][Oo][Cc][Aa][Ll][Ss][Aa]", 'once'))
        % Local SA detrending - separate baseline model per ROI, returns
        % a cell array of a matrix per ROI
        baseline = cell(1,size(data_in, 2));
        for k = 1:size(data_in, 2)
            if b0_bool
                baseline{k} = [ones(T, 1) spectral_baseline(data_in(:,k), AIF, nuisance, T, dt)];
            else
                baseline{k} = spectral_baseline(data_in(:,k), AIF, nuisance, T, dt);
            end
        end
    else
        ME = MException("Error creating baseline: did not recognize detrending type");
        throw(ME);
    end
end


%% Helper Functions
function [a b c p q] = fit_biexp(x, y)
    n = length(x);
    S = 0;
    SS = 0;
    for k = 2:n
        Sk = S(k-1) + 0.5*(y(k) + y(k-1))*(x(k) - x(k-1));
        S = [S Sk];
        SSk = SS(k-1) + 0.5*(S(k) + S(k-1))*(x(k) - x(k-1));
        SS = [SS SSk];
    end
    S = S'; SS = SS';

    mat1 = [sum(SS.^2)      sum(SS.*S)     sum(SS.*(x.^2)) sum(SS.*x) sum(SS); ...
            sum(SS.*S)      sum(S.^2)      sum(S.*(x.^2))  sum(S.*x)  sum(S); ...
            sum(SS.*(x.^2)) sum(S.*(x.^2)) sum(x.^4)       sum(x.^3)  sum(x.^2); ...
            sum(SS.*x)      sum(S.*x)      sum(x.^3)       sum(x.^2)  sum(x); ...
            sum(SS)         sum(S)         sum(x.^2)       sum(x)     n];

    mat2 = [sum(SS.*y) sum(S.*y) sum((x.^2).*y) sum(x.*y) sum(y)]';

    fit1 = pinv(mat1)*mat2;
    p = real(0.5*(fit1(2) + sqrt(fit1(2)^2 + 4*fit1(1))));
    q = real(0.5*(fit1(2) - sqrt(fit1(2)^2 + 4*fit1(1))));
    beta = exp(p*x);
    eta = exp(q*x);

    mat3 = [n         sum(beta)      sum(eta); ...
            sum(beta) sum(beta.^2)   sum(beta.*eta); ...
            sum(eta)  sum(beta.*eta) sum(eta.^2)];

    mat4 = [sum(y) sum(beta.*y) sum(eta.*y)]';

    fit2 = pinv(mat3)*mat4;
    a = fit2(1); b = fit2(2); c = fit2(3);
end
