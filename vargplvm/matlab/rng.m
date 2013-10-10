% A Simple utility that mimics rng.m that is included in newer matlab
% distributions but not in older ones.

% Andreas C. Damianou, 2013

function R = rng(givenSeed)

DEFAULT_SEED = 1e4;

if exist('rand', 'builtin')
    if nargin == 1
        rng(givenSeed)
    else
        R = rng;
    end
else
    if nargin == 1
        if ischar(givenSeed) && strcmp(givenSeed, 'default')
            rand('seed', DEFAULT_SEED);
        else
            rand('seed', givenSeed);
        end
    else % get seed
        R = rand('seed');
    end
end