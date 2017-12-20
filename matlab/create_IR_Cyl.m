function ir_Cyl_patch = create_IR_Cyl(link, freq, delta_f, BSantPattern, MSantPattern)
%GET_IR_CYL Create impulse reponse produced by a compact cylindrical array
%containing 128 antenna ports or 64 dual-polarized patch antennas 
%
%Default call:
%ir_Cyl_patch = create_IR_Cyl(link, freq, delta_f, BSantPattern, MSantPattern)
%
%------
%Input:
%------
%link: Simulated links from the COST 2100 channel model
%freq: Start/end frequencies [Hz]
%delta_f: The difference between two frequency bins
%BSantPattern: Antenna pattern at BS side
%MSantPattern: Antenna pattern at MS side
%------
%Output:
%------
%ir_Cyl_patch: Impulse response for a cylindrical array

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

% Allow for some speed improvement by igonoring MPCs from directions not 
% from the main lobe, slow version means considering all MPCs
% (FAST_VERSION=0)
FAST_VERSION = 1; 

freq_vec = freq(1):delta_f:freq(2); % Vector of frequency points

numBS = size(link, 1); % Number of BSs
numMS = size(link, 2); % Number of MSs
numSS = size(link(1,1).channel, 2); % Number of snapshots
numFreq = length(freq_vec); % Number of frequency points

farfield_BS = BSantPattern.farfield_BS_Rep;
field_MS = MSantPattern.field;

ant_phi_reso = farfield_BS(1).phi(2)-farfield_BS(1).phi(1); % BS antenna pattern azimuth resolution [rad]
ant_theta_reso = farfield_BS(1).theta(2)-farfield_BS(1).theta(1); % BS antenna pattern elevation resolution [rad]
user_ant_phi_reso = field_MS.phi(2)-field_MS.phi(1); % MS antenna pattern azimuth resolution [rad]

% Antenna pattern at BS side, cylindrical array
Nant = 128;
Ntheta = pi/ant_theta_reso+1; % every 5 degrees, elevation
Nphi = 2*pi/ant_phi_reso+1; % every 5 degrees, azimuth
ant_gain_theta = zeros(Nant, Ntheta, Nphi); % vertical-polarized
ant_gain_phi = zeros(Nant, Ntheta, Nphi); % horizontal-polarized
angle_theta = zeros(Nant, Ntheta);
angle_phi = zeros(Nant, Nphi);

for idx_ant = 1:Nant % Loop over 128 ports
    % Vertical-polarized part
    ant_gain_theta(idx_ant, :, :) = double(reshape(farfield_BS(idx_ant).E_theta, length(farfield_BS(idx_ant).theta), length(farfield_BS(idx_ant).phi)));
    % Horizontal-polarized part
    ant_gain_phi(idx_ant, :, :) = double(reshape(farfield_BS(idx_ant).E_phi, length(farfield_BS(idx_ant).theta), length(farfield_BS(idx_ant).phi)));
    
    angle_theta(idx_ant, :) = double(farfield_BS(idx_ant).theta); % 0:5:180 in rad
    angle_phi(idx_ant, :) = double(farfield_BS(idx_ant).phi); % 0:5:360 in rad
    
    % NOTE!
    % phi is measured from positive x-axis [0, 360]
    % theta is measured from positsive z-axis [0, 180]
end

% Antenna pattern at user side 
IDX_f = find(field_MS.f==(freq(1)+(freq(2)-freq(1))/2)); % Find the index of the center frequency
user_ant_gain_v_temp = reshape(field_MS.E_theta(:, IDX_f), length(field_MS.theta), length(field_MS.phi));
user_ant_gain_h_temp = reshape(field_MS.E_phi(:, IDX_f), length(field_MS.theta), length(field_MS.phi));

% Convert phi to [0, 360]
user_ant_gain_v = user_ant_gain_v_temp(:, [(length(field_MS.phi)/2+1):length(field_MS.phi), 1:length(field_MS.phi)/2]); % v-polarized part
user_ant_gain_h = user_ant_gain_h_temp(:, [(length(field_MS.phi)/2+1):length(field_MS.phi), 1:length(field_MS.phi)/2]); % h-polarized part

% Initial rotation
user_ant_rot_h = -1*pi+2*pi*rand(1, numMS); % User initial phase, randomly generated between -pi and pi
user_ant_rot_all_ss = -1*pi+2*pi*rand(1, numMS); % [-pi, pi]
user_ant_rot_per_ss = user_ant_rot_all_ss/numSS;

               
H_Cyl_patch = zeros(numSS, numFreq, numMS, Nant); % [snapshot, freq, user, BS ant]

for idx_BS = 1:numBS
    for idx_MS = 1:numMS
        
        % MS antenna initial rotation
        N_rot = round(user_ant_rot_h(idx_MS)/user_ant_phi_reso)+1; % every phi_reso [rad] (1, numMS)
        if N_rot>=0 
            user_ant_gain_v = user_ant_gain_v(:, [(length(field_MS.phi)-N_rot):length(field_MS.phi), 1:(length(field_MS.phi)-N_rot-1)]);
            user_ant_gain_h = user_ant_gain_h(:, [(length(field_MS.phi)-N_rot):length(field_MS.phi), 1:(length(field_MS.phi)-N_rot-1)]);
        else
            user_ant_gain_v = user_ant_gain_v(:, [abs(N_rot):length(field_MS.phi), 1:abs(N_rot)-1]);
            user_ant_gain_h = user_ant_gain_h(:, [abs(N_rot):length(field_MS.phi), 1:abs(N_rot)-1]);
        end
        
        for idx_ss = 1:numSS
            channel_MPC = link(idx_BS, idx_MS).channel{1, idx_ss}.h; % MPCs
            channel_LOS = link(idx_BS, idx_MS).channel{1, idx_ss}.h_los; % LOS component

            channel = [channel_MPC; channel_LOS];

            MPC_delay = channel(:, 5); % Delay
            MPC_amp_vv = channel(:, 6).*channel(:, 7); % Amplitude vv
            MPC_amp_vh = channel(:, 6).*channel(:, 8); % Amplitude vh
            MPC_amp_hh = channel(:, 6).*channel(:, 9); % Amplitude hh
            MPC_amp_hv = channel(:, 6).*channel(:, 10); % Amplitude hv

            % Convert the MPC angles to fit angles of antenna pattern
            MPC_azi_MS = angle(exp(1j*channel(:, 3))); % Angle in azimuth at MS
            MPC_azi_MS(MPC_azi_MS<0) = MPC_azi_MS(MPC_azi_MS<0)+2*pi;
            MPC_ele_MS = abs(angle(exp(1j*channel(:, 4)))); % Angle in elvation at MS

            MPC_azi_BS = angle(exp(1j*channel(:, 1))); % Angle in azimuth at BS
            MPC_azi_BS(MPC_azi_BS<0) = MPC_azi_BS(MPC_azi_BS<0)+2*pi;
            MPC_ele_BS = abs(angle(exp(1j*channel(:, 2)))-pi/2); % Angle in elevation at BS

            angle_phi_idx_MS = round(MPC_azi_MS/user_ant_phi_reso)+1; % Every user_ant_phi_reso rad
            if ~isempty(find(angle_phi_idx_MS>size(user_ant_gain_v, 2), 1))
                angle_phi_idx_MS(angle_phi_idx_MS>size(user_ant_gain_v, 2)) = 1;
            end
            
            angle_theta_idx_MS = round(MPC_ele_MS/user_ant_phi_reso)+1; % Every user_ant_phi_reso degree
            if ~isempty(find(angle_theta_idx_MS>size(user_ant_gain_v, 1), 1))
                angle_theta_idx_MS(angle_theta_idx_MS>size(user_ant_gain_v, 1)) = 1;
            end

            angle_phi_idx_BS = round(MPC_azi_BS/ant_phi_reso)+1; % Every 5 degree at BS ant pattern
            angle_theta_idx_BS = round(MPC_ele_BS/ant_theta_reso)+1; % Every 5 degree at BS ant pattern
            
            % MS antenna rotation per snapshot
            N_rot_ss = round(user_ant_rot_per_ss(idx_MS)/user_ant_phi_reso)+1; % Every 2.8125 degree
            if N_rot_ss>=0 
                user_ant_gain_v = user_ant_gain_v(:, [(length(field_MS.phi)-N_rot_ss):length(field_MS.phi), 1:(length(field_MS.phi)-N_rot_ss-1)]);
                user_ant_gain_h = user_ant_gain_h(:, [(length(field_MS.phi)-N_rot_ss):length(field_MS.phi), 1:(length(field_MS.phi)-N_rot_ss-1)]);
            else
                user_ant_gain_v = user_ant_gain_v(:, [abs(N_rot_ss):length(field_MS.phi), 1:abs(N_rot_ss)-1]);
                user_ant_gain_h = user_ant_gain_h(:, [abs(N_rot_ss):length(field_MS.phi), 1:abs(N_rot_ss)-1]);
            end

            mpc_user_gain_v = diag(user_ant_gain_v(angle_theta_idx_MS, angle_phi_idx_MS));
            mpc_user_gain_h = diag(user_ant_gain_h(angle_theta_idx_MS, angle_phi_idx_MS));

            % Combine with antenna pattern at BS
            for idx_ant = 1:Nant

                if (FAST_VERSION == 0)
                    keep_IDX = 1:length(MPC_azi_BS);
                else
                    % Allow for some speed improvement by igonoring MPCs
                    % from directions not from the main lobe
                    [~, max_IDX] = max(squeeze(ant_gain_theta(idx_ant, round(Ntheta/2), :)).^2+squeeze(ant_gain_phi(idx_ant, round(Ntheta/2), :)).^2);
                    max_angle = angle_phi(idx_ant, max_IDX)+pi+pi/3; % +60 degrees
                    if max_angle>=2*pi
                        max_angle = max_angle-2*pi;
                    end
                    min_angle = angle_phi(idx_ant, max_IDX)+pi-pi/3; % -60 degrees
                    if min_angle>=2*pi
                        min_angle = min_angle-2*pi;
                    end

                    if max_angle > min_angle
                        keep_IDX = find(MPC_azi_BS<min_angle | MPC_azi_BS>max_angle);
                    else
                        keep_IDX = find(MPC_azi_BS<min_angle & MPC_azi_BS>max_angle);
                    end
                end
                              
                if ~isempty(keep_IDX)       
                    % Vertical polarization
                    complex_gain_from_v = MPC_amp_vv(keep_IDX).*diag(squeeze(ant_gain_theta(idx_ant, angle_theta_idx_BS(keep_IDX), angle_phi_idx_BS(keep_IDX))))+...
                                          MPC_amp_vh(keep_IDX).*diag(squeeze(ant_gain_phi(idx_ant, angle_theta_idx_BS(keep_IDX), angle_phi_idx_BS(keep_IDX))));

                    complex_gain_from_h = MPC_amp_hv(keep_IDX).*diag(squeeze(ant_gain_theta(idx_ant, angle_theta_idx_BS(keep_IDX), angle_phi_idx_BS(keep_IDX))))+...
                                          MPC_amp_hh(keep_IDX).*diag(squeeze(ant_gain_phi(idx_ant, angle_theta_idx_BS(keep_IDX), angle_phi_idx_BS(keep_IDX))));

                    % Add antenna pattern at the MS side                  
                    complex_gain = mpc_user_gain_v(keep_IDX).*complex_gain_from_v+...
                                   mpc_user_gain_h(keep_IDX).*complex_gain_from_h;

                    % Get transfer function for each frequency
                    for idx_freq = 1:numFreq
                        H_Cyl_patch(idx_ss, idx_freq, idx_MS, idx_ant) = sum(complex_gain.*exp(-1j*2*pi*freq_vec(idx_freq)*MPC_delay(keep_IDX)));
                    end  
                else 
                    % Empty MPCs
                    for idx_freq = 1:numFreq
                        H_Cyl_patch(idx_ss, idx_freq, idx_MS, idx_ant) = 0;
                    end  
                end
            end
        end
    end
end

% Convert channel transfer function to impulse response
ir_Cyl_patch = ifft(H_Cyl_patch, [], 2); 