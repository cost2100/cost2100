function C = eadf(phi, theta, eadf)
%EADF Generate antenna response at (phi,theta) in [0,1)x[0,1) from EADF.
%eadf is an M-port array
%
%Default call:
%C = eadf(phi, theta, eadf)

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

[M, U, V] = size(eadf); % M: # of antenna ports
                        % U: # of elevation angles
                        % V: # of azimuth angles

theta = theta(:);
u = floor((U-1)/2)+1-U:floor((U-1)/2);
u = ifftshift(u).';
    
phi = phi(:);
v = floor((V-1)/2)+1-V:floor((V-1)/2);
v = ifftshift(v).';
    
C = zeros(M, length(theta));

ff = reshape(eadf, M*U, V);

for z = 1:length(phi)
    vv = exp(1j*2*pi*v*phi(z))/V;
    zz_tmp = ff*vv;
    uu = exp(1j*2*pi*u*theta(z))/U;
    zz = reshape(zz_tmp, M, U)*uu;
    C(:, z) = zz;
end
