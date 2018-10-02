function h_omni_MIMO_VLA = create_IR_omni_MIMO_VLA(link, freq, delta_f, Band)
%CREATE_IR_OMNI_VLA Get the channel impulse reponse matrix for a very large
%array
%
%Default call:
%h_omni_MIMO_VLA = create_IR_omni_MIMO_VLA(link, freq, delta_f, Band)
%
%------
%Input:
%------
%link: simulated links from the COST 2100 channel model
%freq: start/end frequencies [Hz]
%delta_f: the difference between two frequency bins
%Band: 'Wideband' or 'Narrowband'
%------
%Output:
%------
%h_omni_MIMO_VLA: channel impulse response matrix for a very large array

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

switch Band
    case 'Wideband'
        Nfreq = freq(1):delta_f:freq(2);
    case 'Narrowband'
        Nfreq = .5*(freq(1) + freq(2));
end

Nr = 1; % Assuming single-antenna terminals
first_snapshot = 1;
last_snapshot = size(link.channel,2); 
nBS_pos = size(link.channel,1); % number of positions along large array
H = zeros(last_snapshot,length(Nfreq),Nr,nBS_pos);
for jj = first_snapshot:last_snapshot      
    for kk = 1:nBS_pos
        % h = [phi_mpc_BS theta_mpc_BS phi_mpc_MS theta_mpc_MS tau_mpc amp_h]
        % unfiletered impulse response, angle in radii
        channel1 = link.channel{kk,jj}.h;
        channel2 = link.channel{kk,jj}.h_los;
        channel = [channel1; channel2];
        MPC_delay = channel(:,5); % Delay
        MPC_amp = channel(:,6); % Amplitude
        H11 = MPC_amp.' * exp(-1i*2*pi*MPC_delay*Nfreq);
        H(jj,:,1,kk) = H11;
    end
end

h_omni_MIMO_VLA = ifft(H,[],2);
