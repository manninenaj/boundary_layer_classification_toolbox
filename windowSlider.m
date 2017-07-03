function [ Xo,nsamples,n_nonnan_samples ] = windowSlider( Xi, win, f, edges, prc_val )
%WINDOWSLIDER calculates mean, median, standard deviation, variance, min,
%max, skewness, kurtosis, or fraction on ones within a sliding window.
%
% Usage:
% Xo = windowSlider(Xi,win)
% Xo = windowSlider(Xi,win,f)
% Xo = windowSlider(Xi,win,f,edges)
%
% Inputs (default):
% - Xi        Input array, either a vector or a matrix
% - win       Size of sliding window, either scalar (3) or vector ([3,3])
% - f         handle: (@nanmean), @nanmedian, @nanvar, @nanstd, @nanmax, @nanmin
% - edges     Edges padded with ('nans'), or with 'zeros'
%
% Outputs:
%  - Xo        Calculated quantity in the same dimensions as the input array
%
% 15 June 2017
% Antti Manninen
% antti.j.manninen(at)helsinki.fi
% Department of Physics, University of Helsinki, Finland


% Check the inputs and defaults
if nargin < 1
    error 'Not enough input arguments!'
elseif nargin < 2
    switch isvector(Xi)
        case 1
            win = 3;
        case 0
            win = [3,3];
    end
    f = @nanmean;
    edges = 'nans';
elseif nargin < 3
    f = @nanmean;
    edges = 'nans';
elseif nargin < 4
    edges = 'nans';
elseif nargin > 4 && isempty(edges)
    edges = 'nans';
elseif nargin > 5
    error 'Too many input arguments!'
end

% Check that the inputs meet the requirements
switch isvector(Xi) || ismatrix(Xi)
    case 1
        switch isnumeric(win) && length(win) < 3 && ...
                (isscalar(win) || ...
                isvector(win))
            case 1
                if length(win) > 1
                    switch abs(win(1)) <= intmax% && ...
                            %mod(win(1),2)
                        case 1
                            switch isa(f, 'function_handle')% && ...
%                                     (strcmp(func2str(f),'nanmean') || ...
%                                     strcmp(func2str(f),'nanmedian') || ...
%                                     strcmp(func2str(f),'nanstd') || ...
%                                     strcmp(func2str(f),'nanvar') || ...
%                                     strcmp(func2str(f),'nanmax') || ...
%                                     strcmp(func2str(f),'nanmin') || ...
%                                     strcmp(func2str(f),'skewness') || ...
%                                     strcmp(func2str(f),'kurtosis') || ...
%                                     strcmp(func2str(f),'fraction'))

                                case 1
                                    switch ischar(edges) && ...
                                            (strcmp(edges,'nans') || ...
                                            strcmp(edges,'zeros'))
                                        case 0
                                            error 'The fourth input has to be a string specifying how the edges are handled: ''nans'' (default), or ''zeros''!'
                                    end
                                case 0
                                    error 'The third input has to be a function handle: @nanmean (default), @nanmedian, @nanvar, @nanstd, @skewness, or @kurtosis!'
                            end
                        case 0
                            error 'The second input, size of the sliding window, has to be less than 'intmax'!'
                    end
                end
            case 0
                error 'The second input specifies the length (if input array is a vector) or size (if input array is a matrix) of the sliding window, and has to be scalar or vector with max length of 2!'
        end
    case 0
        error 'The first input has to be a 1-D or 2-D array!'
end

% Determine padding
switch edges
    case 'nans'
        Xi_p = [nan(floor(win(1)/2),size(Xi,2)); Xi;...
            nan(floor(win(1)/2),size(Xi,2))];
        Xi_p = [nan(size(Xi_p,1),floor(win(2)/2)), Xi_p,...
            nan(size(Xi_p,1),floor(win(2)/2))];
    case 'zeros'
        Xi_p = [zeros(floor(win(1)/2),size(Xi,2)); Xi;...
            zeros(floor(win(1)/2),size(Xi,2))];
        Xi_p = [zeros(size(Xi_p,1),floor(win(2)/2)), Xi_p,...
            zeros(size(Xi_p,1),floor(win(2)/2))];
end

if strcmp(func2str(f),'fraction')
    warning('With ''@fraction'' the input array has to be ones and zeros!')
end

% Calculate the quantity
Xo_p = Xi_p .* 0;
n_nonnan_samples = Xi_p .* 0;
nsamples = Xi_p .* 0;
for i = 1+floor(win(1)/2):size(Xi_p,1)-floor(win(1)/2)
    for j = 1+floor(win(2)/2):size(Xi_p,2)-floor(win(2)/2)
        % Extract sub sections and calculate the quantity
        tmp = Xi_p(i-floor(win(1)/2):i+floor(win(1)/2),...
            j-floor(win(2)/2):j+floor(win(2)/2));
        switch func2str(f)
            case 'fraction'
                Xo_p(i,j) = sum(tmp(:))./numel(tmp(:));
            case 'prctile'
                if isempty(prc_val)
                    error 'Specify percentile as the 5th input!'
                else
                    Xo_p(i,j) = prctile(tmp(:),prc_val);
                end
            case 'skewness'
                Xo_p(i,j) = my_skewness(tmp(:));
                n_nonnan_samples(i,j) = sum(~isnan(tmp(:)));
                nsamples(i,j) = numel(tmp(:));
            otherwise
                Xo_p(i,j) = f(tmp(:));
                n_nonnan_samples(i,j) = sum(~isnan(tmp(:)));
                nsamples(i,j) = numel(tmp(:));
        end
    end
end

% Remove padded edges
Xo = Xo_p(1+floor(win(1)/2):end-floor(win(1)/2),...
    1+floor(win(2)/2):end-floor(win(2)/2));

end

