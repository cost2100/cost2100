function ir_Cyl_patch = create_IR_Cyl_EADF(link, freq, delta_f, BSantEADF, MSantPattern, user_ant_rot_h, user_ant_rot_all_ss)
%GET_IR_CYL Create impulse reponse produced by a compact cylindrical EADF array
%containing 128 antenna ports or 64 dual-polarized patch antennas 
%
%Default call:
%ir_Cyl_patch = create_IR_Cyl_EADF(link, freq, delta_f, BSantEADF, MSantPattern)
%
%------
%Input:
%------
%link: Simulated links from the COST 2100 channel model
%freq: Start/end frequencies [Hz]
%delta_f: The difference between two frequency bins
%BSantEADF: EADF at BS side
%MSantPattern: Antenna pattern at MS side
%user_ant_rot_h: Users' initial orientation
%user_ant_rot_all_ss: Users' overall rotation
%------
%Output:
%------
%ir_Cyl_patch: Impulse response for a cylindrical EADF array

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

freq_vec = freq(1):delta_f:freq(2); % Vector of frequency points

numBS = size(link, 1); % Number of BSs
numMS = size(link, 2); % Number of MSs
numSS = size(link(1,1).channel, 2); % Number of snapshots
numFreq = length(freq_vec); % Number of frequency points

field_MS = MSantPattern.field;

user_ant_phi_reso = field_MS.phi(2)-field_MS.phi(1); % MS antenna pattern azimuth resolution [rad]

% Extract the EADF at BS side - cylindrical array
Nant = 128;
for i = 1:Nant
    eadf_BS_theta(i,:,:) = BSantEADF(i).E_theta;
    eadf_BS_phi(i,:,:) = BSantEADF(i).E_phi;
end

% Antenna pattern at user side 
IDX_f = find(field_MS.f==(freq(1)+(freq(2)-freq(1))/2)); % find the index of the center frequency
user_ant_gain_v_temp = reshape(field_MS.E_theta(:, IDX_f), length(field_MS.theta), length(field_MS.phi));
user_ant_gain_h_temp = reshape(field_MS.E_phi(:, IDX_f), length(field_MS.theta), length(field_MS.phi));

% Convert phi to [0, 360]
user_ant_gain_v_temp = user_ant_gain_v_temp(:, [(length(field_MS.phi)/2+1):length(field_MS.phi), 1:length(field_MS.phi)/2]); % v-polarized part
user_ant_gain_h_temp = user_ant_gain_h_temp(:, [(length(field_MS.phi)/2+1):length(field_MS.phi), 1:length(field_MS.phi)/2]); % h-polarized part

% Initial rotation
if (nargin<6)
    user_ant_rot_h = -1*pi+2*pi*rand(1, numMS); % User initial phase, ramdonly generate between -pi and pi
end
if (nargin<7)
    user_ant_rot_all_ss = -1*pi+2*pi*rand(1, numMS); % [-pi, pi]
end
user_ant_rot_per_ss = user_ant_rot_all_ss/numSS;

H_Cyl_patch = zeros(numSS, numFreq, numMS, Nant); % [snapshot, freq, user, BS ant]

for idx_BS = 1:numBS
    for idx_MS = 1:numMS
        
        % MS antenna initial rotation
        N_rot = round(user_ant_rot_h(idx_MS)/user_ant_phi_reso)+1; % Every phi_reso [rad] (1, numMS)
        if N_rot>=0 
            user_ant_gain_v = user_ant_gain_v_temp(:, [(length(field_MS.phi)-N_rot):length(field_MS.phi), 1:(length(field_MS.phi)-N_rot-1)]);
            user_ant_gain_h = user_ant_gain_h_temp(:, [(length(field_MS.phi)-N_rot):length(field_MS.phi), 1:(length(field_MS.phi)-N_rot-1)]);
        else
            user_ant_gain_v = user_ant_gain_v_temp(:, [abs(N_rot):length(field_MS.phi), 1:abs(N_rot)-1]);
            user_ant_gain_h = user_ant_gain_h_temp(:, [abs(N_rot):length(field_MS.phi), 1:abs(N_rot)-1]);
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

            angle_phi_idx_MS = round(MPC_azi_MS/user_ant_phi_reso)+1; % Every user_ant_phi_reso rad
            if ~isempty(find(angle_phi_idx_MS>size(user_ant_gain_v, 2), 1))
                angle_phi_idx_MS(angle_phi_idx_MS>size(user_ant_gain_v, 2)) = 1;
            end
            
            angle_theta_idx_MS = round(MPC_ele_MS/user_ant_phi_reso)+1; % Every user_ant_phi_reso degree
            if ~isempty(find(angle_theta_idx_MS>size(user_ant_gain_v, 1), 1))
                angle_theta_idx_MS(angle_theta_idx_MS>size(user_ant_gain_v, 1)) = 1;
            end

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
            MPC_azi_BS = mod(channel(:, 1),2*pi); % Angle in azimuth at BS
            MPC_zen_BS = elevation2zenith(channel(:, 2)); % Angle in elevation at BS

            if (1)
                % The following must hold, so that the EADF can be called
                % correctly
                if (sum(MPC_azi_BS<0 | MPC_azi_BS>2*pi)>0)
                    error('MPC BS azimuth angles out of range.');
                end
                if (sum(MPC_zen_BS<0 | MPC_zen_BS>pi)>0)
                    error('MPC BS azimuth angles out of range.');
                end
            end
                
            % Obtain the BS antenna response from the EADF
            mpc_BS_gain_v = eadf(MPC_azi_BS/2/pi, MPC_zen_BS/2/pi, eadf_BS_theta);
            mpc_BS_gain_h = eadf(MPC_azi_BS/2/pi, MPC_zen_BS/2/pi, eadf_BS_phi);            
            
            % Activate for visualization/debugging
            if (0)
                figure;
                plot(MPC_azi_BS/2/pi); 
                hold on; 
                plot(MPC_zen_BS/2/pi,'r');
                hold on;
                axis([1 length(MPC_azi_BS) 0 1]);
                grid on;

                figure;
                for i = 1:32
                    subplot(4,8,mod(i-1,32)+1);
                    plot(20*log10(abs(mpc_BS_gain_v(i,:))));
                    hold on;
                    plot(20*log10(abs(mpc_BS_gain_h(i,:))),'r');
                    hold on;
                    axis([1 length(MPC_azi_BS) -50 0]);
                    grid on;
                end
            end            
            
            % Vertical polarization
            complex_gain_from_v = mpc_BS_gain_v*diag(MPC_amp_vv) + mpc_BS_gain_h*diag(MPC_amp_vh);

            % Horizontal polarization
            complex_gain_from_h = mpc_BS_gain_v*diag(MPC_amp_hv) + mpc_BS_gain_h*diag(MPC_amp_hh);
            
            % Activate for visualization/debugging
            if (0)
                figure;
                max_val_v = max(20*log10(abs(vec(complex_gain_from_v))));
                max_val_h = max(20*log10(abs(vec(complex_gain_from_h))));
                max_val = max(max_val_v,max_val_h);
                min_val_v = min(20*log10(abs(vec(complex_gain_from_v))));
                min_val_h = min(20*log10(abs(vec(complex_gain_from_h))));
                min_val = min(min_val_v,min_val_h);
                for i = 1:32
                    subplot(4,8,mod(i-1,32)+1);
                    plot(20*log10(abs(complex_gain_from_v(i,:))));
                    hold on;
                    plot(20*log10(abs(complex_gain_from_h(i,:))),'r');
                    hold on;
                    axis([1 length(MPC_azi_BS) max_val-50 max_val]);
                    grid on;
                end
            end
            
            % Add antenna pattern at the MS side
            complex_gain = complex_gain_from_v*diag(mpc_user_gain_v) +...
                           complex_gain_from_h*diag(mpc_user_gain_h);
            
            % Activate for visualization/debugging
            if (0)
                figure;
                max_val = max(20*log10(abs(vec(complex_gain))));
                min_val = min(20*log10(abs(vec(complex_gain))));
                for i = 97:128
                    subplot(4,8,mod(i-1,32)+1);
                    plot(20*log10(abs(complex_gain(i,:))));
                    hold on;
                    axis([1 length(MPC_azi_BS) max_val-50 max_val]);
                    grid on;
                end
            end

            % Get transfer function for each frequency
            for idx_freq = 1:numFreq
                H_Cyl_patch(idx_ss, idx_freq, idx_MS, :) = complex_gain*exp(-1j*2*pi*freq_vec(idx_freq)*MPC_delay);
            end      
        end
    end
end

% Convert channel transfer function to impulse response
ir_Cyl_patch = ifft(H_Cyl_patch, [], 2); 

function theta_zenith = elevation2zenith(theta_elevation)
theta_zenith = pi/2 - theta_elevation;
