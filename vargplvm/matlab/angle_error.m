% ANGLE_ERROR 
% VARGPLVM

function [m,m_vec] = angle_error(Y1,Y2)

% usage [m,e] = angle_error(Y1,Y2); 
% where Y1 and Y2 are matrices of same size
%
% e is a 54 vector corr to rms error in each angle across all data points
% m is the mean of e
%
% ERRORS GIVEN ARE IN DEGREES

[d,n] = size(Y1);

mat_error = acos(cos((pi/180)*(Y1-Y2)));
mat_error = (180/pi)*mat_error;
% cos and acos to handle wrap around effects in angle space

sq_mat_error = mat_error.*mat_error;

m_vec = zeros(d,1);

%m_vec = mean(mat_error,2); 
%L1 norm
m_vec = sqrt(mean(sq_mat_error,2));

m = mean(m_vec);