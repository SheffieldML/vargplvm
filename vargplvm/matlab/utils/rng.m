% A Simple utility that mimics rng.m that is included in newer matlab
% distributions but not in older ones.

% Andreas C. Damianou, 2013

function R = rng(varargin)

DEFAULT_SEED = 1e4;

if exist('rng', 'builtin')
    if nargin == 1
        rng(varargin{:})
    elseif nargin < 1
        R = rng;
    end
else
    if nargin == 1
        givenSeed = varargin{1};
        if ischar(givenSeed) && strcmp(givenSeed, 'default')
            rand('seed', DEFAULT_SEED);
        else
            rand('seed', givenSeed);
        end
    elseif nargin < 1 % get seed
        R = rand('seed');
    elseif nargin > 1
        warning('Arguments 2, ... ignored!')
    end
end