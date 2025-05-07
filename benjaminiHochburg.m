function posBool = benjaminiHochburg(pVals, alpha, varargin)

    for i = 1:length(varargin)
        if strcmpi(varargin{i}, 'tvals')
            tVals = pVals;
            df = varargin{i+1}; i = i+1;
            pVals = (1 - tcdf(abs(tVals),df))*2;
        end
    end

    [pValsSorted, ordering] = sort(pVals);

    m = length(pVals);

    k = find(pValsSorted <= alpha*(1:m)'/m, 1, 'last');

    posIdx = ordering(1:k);

    posBool = zeros(length(pVals), 1);

    posBool(posIdx) = 1;

    posBool = (posBool == 1);

end