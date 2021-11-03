function [vX] = spm_vec(X,varargin)
% Vectorise a numeric, cell or structure array - a compiled routine
% FORMAT [vX] = spm_vec(X)
% X  - numeric, cell or structure array[s]
% vX - vec(X)
%
% See spm_unvec
%__________________________________________________________________________
%
% e.g.:
% spm_vec({eye(2) 3}) = [1 0 0 1 3]'
%__________________________________________________________________________
% Copyright (C) 2005-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_vec.m 6110 2014-07-21 09:36:13Z karl $


%error('spm_vec.c not compiled - see Makefile')

% initialise X and vX
%--------------------------------------------------------------------------
if nargin > 1
    X = [{X},varargin];
end


% vectorise numerical arrays
%--------------------------------------------------------------------------
if isnumeric(X)
    vX = X(:);

% vectorise logical arrays
%--------------------------------------------------------------------------
elseif islogical(X)
    vX = X(:);

% vectorise structure into cell arrays
%--------------------------------------------------------------------------
elseif isstruct(X)
    vX = [];
    f   = fieldnames(X);
    X    = X(:);
    for i = 1:numel(f)
        vX = cat(1,vX,spm_vec({X.(f{i})}));
    end

% vectorise cells into numerical arrays
%--------------------------------------------------------------------------
elseif iscell(X)
    vX   = [];
    for i = 1:numel(X)
        vX = cat(1,vX,spm_vec(X{i}));
    end
else
    vX = [];
end
