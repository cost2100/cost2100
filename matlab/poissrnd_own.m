function r = poissrnd_own(lambda)
%POISSRND_OWN Draw values from a Poisson-distributed random variable
%r = poissrnd_own(lambda) output a Poisson random generated variable

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

% Reference:
% https://www.johndcook.com/blog/2010/06/14/generating-poisson-random-values/

if (lambda > 100)
    error('Only lambda <= 100 has been statistically tested (see test_poissrnd_own.m). To run the routine, comment out this error message at your own risk.');
end

L = exp(-1*lambda); % Interval between event occurences
p = 1; % Unit interval
k = 0; % Occurence count

while 1
    k = k+1; % Occurence count
    p = p*rand(1); % Rest of the interval
    % p <= L, unit interval finishes
    if (p<=L)
        break;
    end
end

r = k-1;
end