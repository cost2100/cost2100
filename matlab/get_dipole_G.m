function [Gtx, Grx] = get_dipole_G(Ntx,Nrx)
%get_dipole_G Get the dipole array antenna responses
%Default call 
%-------
%Input:
%-------
%
%Ntx: Number of transmit antenna
%Nrx: Number of receive antenna
%
%------
%Output:
%------
%Gtx: Antenna response at Tx side
%Grx: Antenna response at Rx side

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (C)2008 Meifang Zhu, Lund University, Sweden
%This file is part of cost2100.
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% therotical dipole antenna gain
Theta = (-90:1:90)/180*pi;
Phi = (0:1:360)/180*pi;
Np = length(Phi);
Nt = length(Theta);
G = zeros(Np,Nt);
for nn = 1:Np
    for mm = 1:Nt
        theta = Theta(mm);
        if theta == -90/180*pi || theta == 90 /180*pi
            G(nn,mm) = 0;
        else
            G(nn,mm) = cos(pi/2*sin(theta))/cos(theta);
        end
    end
end


% generate the response for linear antenna array, when they are placed with
% distance half wavelength
Anttx = zeros(Ntx,Np,Nt);
for nn = 1:Ntx
    for mm = 1:size(G,1)
        phi = Phi(mm);
        Anttx(nn,mm,:) = G(mm,:)*exp(-1j*pi*cos(phi)*(nn-1));
    end
end

Gtx.antennaResponse = Anttx;
Gtx.minResponse = min(min(min(abs(Anttx))));
Gtx.azimuthRange = Phi;
Gtx.elevationRange = Theta; 
Gtx.dangle = 1/180*pi;

Antrx = zeros(Nrx,Np,Nt);
for nn = 1:Nrx
    for mm = 1:size(G,1)
        phi = Phi(mm);
        Antrx(nn,mm,:) = G(mm,:)*exp(-1j*pi*cos(phi)*(nn-1));
    end
end

Grx.antennaResponse = Anttx;
Grx.minResponse = min(min(min(abs(Anttx))));
Grx.azimuthRange = Phi;
Grx.elevationRange = Theta; 
Grx.dangle = 1/180*pi;