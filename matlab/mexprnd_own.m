function r = mexprnd_own(Leff, L, m, n)
%MEXP_OWN Function to compute a modified exponential random variable
%with pdf % (L+x)/(L+Leff)*1/Leff*exp(-x/Leff)
%
%Default call:
%r = mexprnd_own(Leff, L, m, m)
%
%------
%Input:
%------
%Leff: BS-VR length
%L: Array length
%m: Number of output rows 
%n: Number of output columns
%------
%Output:
%------
%r: The random variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This file is a part of the COST2100 channel model.
%
%This program, the COST2100 channel model, is free software: you can 
%redistribute it and/or modify it under the terms of the GNU General Public 
%License as published by the Free Software Foundation, either version 3 of 
%the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful, but 
%WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
%or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
%for more details.
%
%If you use it for scientific purposes, please consider citing it in line 
%with the description in the Readme-file, where you also can find the 
%contributors.
%
%You should have received a copy of the GNU General Public License along 
%with this program. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sanity check.
if (nargin < 3)
    m = 1;
end

if (nargin < 4)
    n = 1;
end

if ((L<0) || (Leff<0))
    % Illegal parameter values
    r = zeros(m,n);
elseif (Leff==Inf)
    % Set all to infinite.
    r = inf*ones(m,n);
else
    % Generate random variables with the desired PDF by the inverse CDF
    % method (see Bishop's book "Pattern Recognition and Machine Learning."
    t1 = rand(m,n);
    r = -Leff*lambertw(-1,(t1-1)*exp(-(L+Leff)/Leff)*(L+Leff)/Leff) - L - Leff;
end