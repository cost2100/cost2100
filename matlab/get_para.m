function [paraEx, paraSt] = get_para(network, scenario, freq, snapRate, snapNum, posBS, spacingBS, posBSnum, posMS, veloMS)
%GET_PARA Generate the external and stochastic parameters of the scenario
%
%Default call: 
%[paraEx paraSt] = get_para(network, scenario, freq, snapRate, snapNum, posBS, spacingBS, posBSnum, posMS, veloMS)
%
%------
%Input:
%------
%network :
%'IndoorHall_5GHz','SemiUrban_300MHz','Indoor_CloselySpacedUser_2_6GHz','SemiUrban_CloselySpacedUser_2_6GHz', or 'SemiUrban_VLA_2_6GHz'.
%scenario: 'LOS' or 'NLOS'
%freq: Frequency range, [freq_start freq_stop] [Hz]
%snapRate: Number of snapshots per s [Hz]
%snapNum: Number of simulated snapshots
%posBS: Position matrix of BSs, size = (number of BS, [x y z]) [m]
%spacingBS: Inter-position spacing, for large arrays [m]
%posBSnum: Number of positions at each BS site, for large arrays
%posMS: Initial position matrix of MSs, size = (number of MS, [x y z]) [m]
%veloMS: Velocity vector of MSs, size = (number of MS, [x y z]) [m/s]
%
%------
%Output:
%------
%paraEx: External parameters
% .network: Network type
% .freq_start: Start frequency [Hz]
% .freq_stop: Stop frequency [Hz]
% .center_freq: Carrier frequency [Hz]
% .net_radii: Network cell radius [m]
% .sample_rate: Delay sampling rate [s]
% .c0: Wave speed [m/s]
% .pos_BS: BS position [m]
% .pos_MS: MS position [m]
% .velo_MS: MS moving speed vector [m/s]
% .snap_rate: Channel snapshot rate [s]
% .snap_num: Channel snapshot number [no unit]
% .bandwidth: Bandwidth [Hz]
% .delay_max: Maximum delay, 5 times of net_radii [s]
% .nfloor: Average number of floors between BS and MS [no unit] (picocell)
% .h_BS: BS height [m] (macrocell)
% .h_MS: MS height [m] (macrocell)
% .h_rooftop: Rooftop height [m] (macrocell)
% .w_r: Road width [m] (macrocell)
% .w_b: Street width [m] (macrocell)
% .phi_road: Road orientation [deg] (macrocell)

%paraSt: Stochastic parameters
% .r_c: Visibility region radius for clusters [m]
% .l_c: Transition region radius for clusters [m]
% .k_tau: Cluster power attenuation coefficient [dB/us]
% .tau_b: Cut-off delay for cluster power attenuation, equation 6.74 [s]
% .d_co: LOS cut-off distance [m]
% .r_l: Visibility region radius for LOS [m]
% .l_l: Transition region radius for LOS [m]
% .mu_k: LOS power factor mean  [dB]
% .sigma_k: LOS power factor dB scale std [dB]
% .BSLocal: Activity of BS local cluster
% .MSLocal: Activity of MS local cluster
% .mu_n_c_far: Average number of far clusters, mean
% .n_c_far: Average number of far clusters in the scenario
% .k_sel: Proportion of twin/single clusters
% .n_mpc: Number of MPCs per cluster
% .mu_tau: Delay spread mean [s]
% .sigma_tau: Delay spread dB scale std [dB]
% .mu_phi_BS: AoD mean [deg]
% .sigma_phi_BS: AoD dB std [dB]
% .mu_theta_BS: EoD mean [deg]
% .sigma_theta_BS: EoD dB scale std [dB]
% .mu_phi_MS: AoA mean [deg]
% .sigma_phi_MS: AoA dB scale std [dB]
% .mu_theta_MS: EoA mean [deg]
% .sigma_theta_MS: EoA dB scale std [dB]
% .sigma_sh: Std of shadowing [dB]
% .rho: Correlation matrix
% .corr_mat: Cholesky factorized correlation matrix
% .phi_c: Mean/std azimuth of cluster to VR, mean [deg]  
% .pdf_phi_c: phi_c distribution
% .theta_c: Mean/std elevation of cluster to VR, mean [deg] 
% .pdf_theta_c: theta_c distribution
% .para_r_c: Statistics distance cluster-BS/VR [m]
% .pdf_r_c: r_c distribution
% .mu__tauCLink: Cluster link delay extra, mean [s]
% .sigma__tauCLink: Cluster link delay extra, dB scale std
% .sigma_sf: Shadow fading, dB std, [dB]
% .BS_common: Ratio of BS-common cluster to total cluster ([2,3,...]BS-CC)
% .MS_common: Average number of VRs in VR group
% .dmc_*: Parameters of DMC

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

paraEx.network = network; % Network type
paraEx.scenario = scenario; % Scenario type
paraEx.overSample = 4; % Oversampling factor for the impulse response shaping

switch network
    case 'IndoorHall_5GHz'
        paraEx.freq_start = freq(1); % Reference start frequency [Hz]
        paraEx.freq_stop = freq(2); % Reference stop frequency [Hz]
        paraEx.net_radii = 30; % Cell radius [m]    
        paraEx.n_floor = 0;
    case 'SemiUrban_300MHz'
        paraEx.freq_start = freq(1); % Reference start frequency [Hz]
        paraEx.freq_stop = freq(2); % Reference stop frequency [Hz]
        paraEx.net_radii = 500; % Cell radius [m]    
        paraEx.n_floor = 0; % Number of floors between BS and MS
    case 'Indoor_CloselySpacedUser_2_6GHz'
        paraEx.freq_start = freq(1); % Reference start frequency [Hz]
        paraEx.freq_stop = freq(2); % Reference stop frequency [Hz]
        paraEx.net_radii = 4.5; % Cell radius [m]
        paraEx.n_floor = 0; % Number of floors between BS and MS
    case 'SemiUrban_CloselySpacedUser_2_6GHz'
        paraEx.freq_start = freq(1); % Reference start frequency [Hz]
        paraEx.freq_stop = freq(2); % Reference stop frequency [Hz]
        paraEx.net_radii = 75; % Cell radius [m]
        paraEx.n_floor = 0; % Number of floors between BS and MS
    case 'SemiUrban_VLA_2_6GHz'
        paraEx.freq_start = freq(1); % Reference start frequency [Hz]
        paraEx.freq_stop = freq(2); % Reference stop frequency [Hz]
        paraEx.net_radii = 200; % Cell radius [m]
        paraEx.n_floor = 0; % Number of floors between BS and MS
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common external parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paraEx.snap_rate = snapRate; % Channel snapshot rate [s]
paraEx.snap_num = snapNum; % Channel snapshot number
paraEx.c0 = 3e8; % Wave speed [m/s]
paraEx.center_freq = (freq(1)+freq(2))/2; % Carrier frequency [Hz]
paraEx.bandwidth = freq(2)-freq(1); % Bandwidth [Hz]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BS/MS position and movement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paraEx.num_BS = length(posBS(:,1)); % Number of BS
paraEx.num_MS = length(posMS(:,1)); % Number of MS
paraEx.pos_BS = posBS; % BS position
paraEx.pos_MS = posMS; % MS position
paraEx.velo_MS = veloMS;
paraEx.BS_spacing = spacingBS;
paraEx.BS_spacing_num = posBSnum;
% Base station spatial range [m], for VLA
if (strcmp(network,'SemiUrban_VLA_2_6GHz'))
    paraEx.BS_range = norm(paraEx.BS_spacing)*paraEx.BS_spacing_num;
else
    paraEx.BS_range = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delay resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paraEx.sample_rate = 1/paraEx.bandwidth; % Delay sample rate [s]
paraEx.delay_max = paraEx.net_radii*5/paraEx.c0; % Maximum delay 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get ST parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paraSt = get_paraSt(paraEx);

end

function paraSt = get_paraSt(paraEx)
% Sub function to get the stochastic parameters

switch paraEx.network           
    %%%%%%%%%%%%%%%%%%%%%%%%
    case 'IndoorHall_5GHz' % Aalto Obstructive LOS measurements
    %%%%%%%%%%%%%%%%%%%%%%%%
        paraSt.bs_vr_mu = Inf; % BS-VR length, Inf means BS-VR feature has no effect
        paraSt.bs_vr_sigma = 0; % BS-VR length
        paraSt.bs_vr_slope_mu = 0; % BS-VR power slope, 0 means BS-VR feature has no effect
        paraSt.bs_vr_slope_sigma = 0; % BS-VR power slope, normal, sigma

        paraSt.r_c = 2.72; % Visibility region size [m]
        paraSt.l_c = 1; % Transition region size [m]
        paraSt.k_tau = 16; % Cluster power decaying factor lin. [/us]
        paraSt.tau_b = 30/3e8*1e6; % Cluster power decaying cut-off delay [us] (cell size)
        
        paraSt.d_co = 30; % LOS cut-off distance [m]
        paraSt.r_l = 0; % LOS VR size [m]
        paraSt.l_l = 0; % LOS TR size [m]
        paraSt.mu_k = 10^(-10/10); % LOS power factor, mean linear
        paraSt.sigma_k = 7.2; % LOS power factor, std dB
        paraSt.power_factor = paraSt.mu_k*10^(randn(1)*paraSt.sigma_k/10); % LOS power factor        
                
        paraSt.BSLocal = true; % BS Local cluster activity true/false
        paraSt.MSLocal = true; % MS local cluster activity true/false
        
        paraSt.mu_n_c_far = 3.8; % Average number of far clusters, mean
        n_c_far = poissrnd_own(paraSt.mu_n_c_far); % Possion distribution
        if n_c_far==0 n_c_far = 1; end % Guarantee min 1 far cluster
        paraSt.n_c_far = n_c_far; % Average number of far clusters
        
        paraSt.k_sel = 0; % K-selection factor (ratio of single cluster)
        
        E = 3e-9; % Delay spread mean [s]
        V = 0.6e-9; % Delay spread std [s]
        paraSt.mu_tau = 10^(log10(exp(1))*(log(E)-1/2*log(1+V^2/E^2) )); % Delay spread mean [s]
        paraSt.sigma_tau = 10*log10(exp(1))*sqrt(log(1+V^2/E^2)); % Delay spread dB scale std        
        E=3; % AoD spread mean [deg]
        V=2; % AoD spread std [deg]		
        paraSt.mu_phi_BS = 10^(log10(exp(1))*(log(E)-1/2*log(1+V^2/E^2) )); % AoD spread mean [deg]
        paraSt.sigma_phi_BS = 10*log10(exp(1))*sqrt(log(1+V^2/E^2)); % AoD spread dB scale std        
        E=1.5; % EoD spread mean [deg]
        V=1.5; % EoD spread std [deg]
        paraSt.mu_theta_BS = 10^(log10(exp(1))*(log(E)-1/2*log(1+V^2/E^2) )); % EoD spread mean [deg]
        paraSt.sigma_theta_BS = 10*log10(exp(1))*sqrt(log(1+V^2/E^2)); % EoD dB scale spread
        E=4; % AoA spread mean [deg]
        V=2; % AoA spread std [deg]		
        paraSt.mu_phi_MS = 10^(log10(exp(1))*(log(E)-1/2*log(1+V^2/E^2) )); % AoA spread mean [deg]      
        paraSt.sigma_phi_MS = 10*log10(exp(1))*sqrt(log(1+V^2/E^2)); % AoA dB scale spread std 
        E=1; % EoA spread mean [deg]
        V=1; % EoA spread mean [deg]
        paraSt.mu_theta_MS = 10^(log10(exp(1))*(log(E)-1/2*log(1+V^2/E^2) )); % EoA spread mean [deg]
        paraSt.sigma_theta_MS = 10*log10(exp(1))*sqrt(log(1+V^2/E^2)); % EoA dB scale std       
        paraSt.sigma_sf = 0; % Shadow fading std dB
        
        paraSt.mpc_gain_mu = Inf; % MPC-VR radius, lognormal, mu, Inf. deactivates feature
        paraSt.mpc_gain_sigma = 0; % MPC-VR radius, lognormal, sigma
        		
        paraSt.n_mpc = 3; % Number of MPCs per cluster
        
        % Cross-correlation matrix
        paraSt.rho =[1  0  0  0  0  0;...  % Shadow fading
                     0  1  0  0  0  0;...  % Delay spread
                     0  0  1  0  0  0;...  % theta BS
                     0  0  0  1  0  0;...  % phi BS
                     0  0  0  0  1  0;...  % theta MS
                     0  0  0  0  0  1];    % phi MS   
        paraSt.corr_mat=chol(paraSt.rho); % Cholesky factorization
               
        % Cluster distribution and cluster link delay
		paraSt.phi_c = [43 39]; % Mean/std azimuth of cluster to VR, mean [deg]        
		paraSt.pdf_phi_c = 'norm'; % phi_c normal distribution
		paraSt.theta_c = [0 0]; % Mean/std elevation of cluster to visibility region [deg]    
		paraSt.pdf_theta_c = 'norm';  % theta_c normal distribution
        paraSt.para_r_c = [0 30]; % Min/max distance from cluster to BS/MS        
		paraSt.pdf_r_c = 'unif'; % r_c uniform distribution        
        paraSt.mu_tauCLink = 1.07e-7; % Cluster link delay mean [s]
        paraSt.min_tauCLink = 3.36e-8; % Cluster link delay dB scale std        
        
        % Polarization
        paraSt.mu_xpdv = 0;
        paraSt.sigma_xpdv = 0; % Mean and std for XPDV
        paraSt.mu_xpdh = 0; 
        paraSt.sigma_xpdh = 0; % Mean and std for XPDH
        paraSt.mu_cpr =0;
        paraSt.sigma_cpr = 0; % Mean and std for CPR
        
        % For multi-link extension
        paraSt.BS_common = [0.35 0.65]; % BS-common cluster ratio to total number of clusters (1,2,3,...)
        paraSt.MS_common = 4; % Average number of VRs in VR group for one MS-common cluster
        
        % DMC parameters
        paraSt.dmc_spread_mean = 0; % DMC spatial spread [m]
        paraSt.dmc_spread_sigma = 0; % DMC spatial spread [m]
        paraSt.dmc_beta = 0; % DMC delay power decaying [dB]	
        
    %%%%%%%%%%%%%%%%%%%%%%%%%
    case 'SemiUrban_300MHz' % Parameters for two scenarios were extracted
    %%%%%%%%%%%%%%%%%%%%%%%%%
        switch paraEx.scenario
            %%%%%%%%%%%
            case 'LOS'
            %%%%%%%%%%%             
                paraSt.bs_vr_mu = Inf; % BS-VR length, Inf means BS-VR feature has no effect
                paraSt.bs_vr_sigma = 0; % BS-VR length
                paraSt.bs_vr_slope_mu = 0; % BS-VR power slope, 0 means BS-VR feature has no effect
                paraSt.bs_vr_slope_sigma = 0; % BS-VR power slope, normal, sigma

                paraSt.r_c = 32.8; % Visibility region size [m]
                paraSt.l_c = 16.8; % Transition region size [m]
                paraSt.k_tau = log(10^(12.1/10)); % Cluster power decaying factor lin. [/us]
                paraSt.tau_b = 2.4e-6; % Cluster power decaying cut-off delay [us] (cell size)

                paraSt.d_co = 350; % LOS cut-off distance [m]
                paraSt.r_l = 343; % LOS VR size [m]
                paraSt.l_l = 93; % LOS TR size [m]
                paraSt.mu_k = 10^(-4.7/10); % LOS power factor, mean linear
                paraSt.sigma_k = 2.0; % LOS power factor, std dB
                paraSt.power_factor = paraSt.mu_k*10^(randn(1)*paraSt.sigma_k/10); % LOS power factor        

                paraSt.BSLocal = false; % BS Local cluster activity true/false
                paraSt.MSLocal = true; % MS local cluster activity true/false
               
                n_c_far = 6;
                paraSt.n_c_far = n_c_far; % Average number of far clusters                

                paraSt.k_sel = 0.1; % K-selection factor (ratio of single cluster)

                paraSt.mu_tau = 1.39e-7; % Delay spread mean [s]
                paraSt.sigma_tau = 3.66; % Delay spread dB scale std        

                paraSt.mu_phi_BS =  0.255/pi*180; % AoD spread mean [deg]
                paraSt.sigma_phi_BS = 2.43; %1.9567; % AoD spread dB scale std        

                paraSt.mu_theta_BS = 1.578/pi*180; % EoD spread mean [deg]
                paraSt.sigma_theta_BS = 0; % EoD dB scale spread

                paraSt.mu_phi_MS =  0.2575/pi*180; % AoA spread mean [deg]      
                paraSt.sigma_phi_MS = 2.68; % AoA dB scale spread std 

                paraSt.mu_theta_MS = 1.578; % EoA spread mean [deg]
                paraSt.sigma_theta_MS = 0; % EoA dB scale std       
                paraSt.sigma_sf = 2.05; % Shadow fading std dB
                
                paraSt.mpc_gain_mu = Inf; % MPC-VR radius, lognormal, mu, Inf. deactivates feature
                paraSt.mpc_gain_sigma = 0; % MPC-VR radius, lognormal, sigma

                paraSt.n_mpc = 27; % Number of MPCs per cluster

                % Cross-correlation matrix
                paraSt.rho =[1    0.0  0  0.0  0  0.0;...  % Shadow fading
                             0.0  1    0  0.9  0  0.9;...  % Delay spread
                             0    0    1  0    0  0;...    % theta BS
                             0.0  0.9  0  1    0  0.9;...  % phi BS
                             0    0    0  0    1  0;...    % theta MS
                             0.0  0.9  0  0.9  0  1];      % phi MS   
                paraSt.corr_mat=chol(paraSt.rho); % Cholesky factorization

                % Cluster distribution and cluster link delay
                paraSt.phi_c = [43 39]; % Mean/std azimuth of cluster to VR, mean [deg]        
                paraSt.pdf_phi_c = 'norm'; % phi_c normal distribution
                paraSt.theta_c = [0 0]; % Mean/std elevation of cluster to visibility region [deg]    
                paraSt.pdf_theta_c = 'norm'; % theta_c normal distribution
                paraSt.para_r_c = [0 paraEx.net_radii]; % Min/max distance from cluster to BS/MS        
                paraSt.pdf_r_c = 'unif'; % r_c uniform distribution        
                paraSt.mu_tauCLink = 8.54e-7; % Cluster link delay mean [s]
                paraSt.min_tauCLink = 4.79e-8; % Cluster link delay dB scale std        

                % Polarization
                paraSt.mu_xpdv = 0;
                paraSt.sigma_xpdv = 0; % Mean and std for XPDV
                paraSt.mu_xpdh = 0; 
                paraSt.sigma_xpdh = 0; % Mean and std for XPDH
                paraSt.mu_cpr =0;
                paraSt.sigma_cpr = 0; % Mean and std for CPR

                % For multi-link extension
                if size(paraEx.pos_BS,1)>1
                    error('Parameters for more than 1 BS are not available')                        
                end
                if size(paraEx.pos_MS,1)> 1  
                    if size(paraEx.pos_MS,1) == 2 
                        paraSt.BS_common = 0; % BS-common cluster ratio to total number of clusters (1,2,3,...)
                        dist = sqrt(sum((paraEx.pos_MS(1,:) - paraEx.pos_MS(2,:)).^2));
                        % The modeling of common cluster decay ratio is
                        % not finished yet, so only round the dist to a
                        % certain range for now and get the
                        % corresponding common cluster ratio from the
                        % Eucap paper
                        if dist == 0
                            paraSt.MS_common = 1; % Average number of VRs in VR group for one MS-common cluster
                        elseif dist > 0 && dist <= 5
                            paraSt.MS_common = 0.52; 
                        elseif dist > 5 && dist <= 10
                            paraSt.MS_common = 0.50;
                        elseif dist > 10 && dist <= 15
                            paraSt.MS_common = 0.46;
                        elseif dist > 10 && dist <= 20
                            paraSt.MS_common = 0.44;
                        elseif dist > 20 && dist <= 25
                            paraSt.MS_common = 0.42;
                        elseif dist > 25 && dist <= 30
                            paraSt.MS_common = 0.42;
                        elseif dist > 30 && dist <= 35
                            paraSt.MS_common = 0.38;
                        elseif dist > 35 && dist <= 40
                            paraSt.MS_common = 0.38;
                        else
                            paraSt.MS_common = 0.1; % When the distance is larger than 40 meters, we just give a small common ratio
                        end
                    else
                        error('Parameters for more than 2 MSs are not available')                            
                    end
                else
                    % Single link set up
                    paraSt.BS_common = 0; % BS-common cluster ratio to total number of clusters (1,2,3,...)
                    paraSt.MS_common = 0; % Average number of VRs in VR group for one MS-common cluster
                end

                % DMC parameters
                paraSt.dmc_spread_mean = 0; % DMC spatial spread [m]
                paraSt.dmc_spread_sigma = 0; % DMC spatial spread [m]
                paraSt.dmc_beta = 0; % DMC delay power decaying [dB]
            %%%%%%%%%%%%%%%%
            case 'NLOS'
            %%%%%%%%%%%%%%%%    
                paraSt.bs_vr_mu = Inf; % BS-VR length, Inf means BS-VR feature has no effect
                paraSt.bs_vr_sigma = 0; % BS-VR length
                paraSt.bs_vr_slope_mu = 0; % BS-VR power slope, 0 means BS-VR feature has no effect
                paraSt.bs_vr_slope_sigma = 0; % BS-VR power slope, normal, sigma

                paraSt.r_c = 24.5; % Visibility region size [m]
                paraSt.l_c = 12.2; % Transition region size [m]
                paraSt.k_tau = log(10^(7.2/10)); % Cluster power decaying factor lin. [/us]
                paraSt.tau_b = 4.2e-6; % Cluster power decaying cut-off delay [us] (cell size)

                paraSt.d_co = 0; % LOS cut-off distance [m]
                paraSt.r_l = 0; % LOS VR size [m]
                paraSt.l_l = 0; % LOS TR size [m]
                paraSt.mu_k = 10^(0/10); % LOS power factor, mean linear
                paraSt.sigma_k = 0; % LOS power factor, std dB
                paraSt.power_factor = paraSt.mu_k*10^(randn(1)*paraSt.sigma_k/10); % LOS power factor        

                paraSt.BSLocal = false; % BS Local cluster activity true/false
                paraSt.MSLocal = true; % MS local cluster activity true/false
               
                n_c_far = 6;
                paraSt.n_c_far = n_c_far; % Average number of far clusters
                
                paraSt.k_sel = 0.2; % K-selection factor (ratio of single cluster)

                paraSt.mu_tau = 3.17e-7; % Delay spread mean [s]
                paraSt.sigma_tau = 2.5; % Delay spread dB scale std        

                paraSt.mu_phi_BS =  0.3254/pi*180; % AoD spread mean [deg]
                paraSt.sigma_phi_BS = 2.02; % 1.9567; % AoD spread dB scale std        

                paraSt.mu_theta_BS = 1.578/pi*180; % EoD spread mean [deg]
                paraSt.sigma_theta_BS = 0; % EoD dB scale spread

                paraSt.mu_phi_MS = 0.3318/pi*180; % AoA spread mean [deg]      
                paraSt.sigma_phi_MS = 2.03; % 1.8599; % AoA dB scale spread std 

                paraSt.mu_theta_MS = 1.578; % EoA spread mean [deg]
                paraSt.sigma_theta_MS = 0; % EoA dB scale std       
                paraSt.sigma_sf = 2.27; % Shadow fading std dB

                paraSt.mpc_gain_mu = Inf; % MPC-VR radius, lognormal, mu, Inf. deactivates feature
                paraSt.mpc_gain_sigma = 0; % MPC-VR radius, lognormal, sigma

                paraSt.n_mpc = 48; % Number of MPCs per cluster

                % Cross-correlation matrix
                paraSt.rho =[1    -0.1  0  -0.1  0  -0.1;...  % Shadow fading
                            -0.1   1    0   0.9  0   0.9;...  % Delay spread
                             0     0    1   0    0   0;...    % theta BS
                            -0.1   0.9  0   1    0   0.9;...  % phi BS
                             0     0    0   0    1   0;...    % theta MS
                            -0.1   0.9  0   0.9  0   1];      % phi MS   
                paraSt.corr_mat=chol(paraSt.rho); % Cholesky factorization

                % Cluster distribution and cluster link delay
                paraSt.phi_c = [43 39]; % Mean/std azimuth of cluster to VR, mean [deg]        
                paraSt.pdf_phi_c = 'norm'; % phi_c normal distribution
                paraSt.theta_c = [0 0]; % Mean/std elevation of cluster to visibility region [deg]    
                paraSt.pdf_theta_c = 'norm'; % theta_c normal distribution
                paraSt.para_r_c = [0 paraEx.net_radii]; % Min/max distance from cluster to BS/MS        
                paraSt.pdf_r_c = 'unif'; % r_c uniform distribution        
                paraSt.mu_tauCLink = 1.02e-6; % Cluster link delay mean [s]
                paraSt.min_tauCLink = 5.24e-8; % Cluster link delay dB scale std        

                % Polarization
                paraSt.mu_xpdv = 0;
                paraSt.sigma_xpdv = 0; % Mean and std for XPDV
                paraSt.mu_xpdh = 0; 
                paraSt.sigma_xpdh = 0; % Mean and std for XPDH
                paraSt.mu_cpr =0;
                paraSt.sigma_cpr = 0; % Mean and std for CPR

                
                % For multi-link extension
                paraSt.BS_common = 0; % BS-common cluster ratio to total number of clusters (1,2,3,...)
                paraSt.MS_common = 0; % Average number of VRs in VR group for one MS-common cluster
              

                % DMC parameters
                paraSt.dmc_spread_mean = 0; % DMC spatial spread [m]
                paraSt.dmc_spread_sigma = 0; % DMC spatial spread [m]
                paraSt.dmc_beta = 0; % DMC delay power decaying [dB]
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Indoor_CloselySpacedUser_2_6GHz' % Massive MIMO with a compact array.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch paraEx.scenario
            %%%%%%%%%%%%%%%%
            case 'LOS'
                paraSt.bs_vr_mu = Inf; % BS-VR length, Inf means BS-VR feature has no effect
                paraSt.bs_vr_sigma = 0; % BS-VR length
                paraSt.bs_vr_slope_mu = 0; % BS-VR power slope, 0 means BS-VR feature has no effect.
                paraSt.bs_vr_slope_sigma = 0; % BS-VR power slope, normal, sigma
                
                paraSt.r_c = 5; % Visibility region size at MS [m]
                paraSt.l_c = .5; % Transition region size at MS [m]
 
                paraSt.k_tau = 31; % Cluster power decaying factor lin. /us
                paraSt.tau_b = .25e-6; % Cluster power decaying cut-off delay [s]

                paraSt.d_co = paraEx.net_radii; %LOS cut-off distance [m]
                paraSt.r_l = paraEx.net_radii; %LOS VR size [m]
                paraSt.l_l = 0; %LOS TR size [m]
                paraSt.mu_k = 10^(-5.2/10); %LOS power factor, mean linear
                paraSt.sigma_k = 2.9; %LOS power factor, std dB
                paraSt.power_factor = paraSt.mu_k*10^(randn(1)*paraSt.sigma_k/10); %LOS power factor

                paraSt.BSLocal = false; % BS Local cluster activity true/false
                paraSt.MSLocal = false; % MS local cluster activity true/false
               
                n_c_far = 15;
                paraSt.n_c_far = n_c_far; % Average number of far clusters                

                paraSt.k_sel = 1; % K-selection factor (ratio of the number of single bounce clusters relative to twin clusters)
                
                E = 5.5e-9; %Delay spread mean [s]
                V = 2.0e-9; %Delay spread std [s]
                paraSt.mu_tau = 10^(log10(exp(1))*(log(E)-1/2*log(1+V^2/E^2) )); %Delay spread mean [s]
                paraSt.sigma_tau = 10*log10(exp(1))*sqrt(log(1+V^2/E^2)); %Delay spread dB scale std

                E = 5.2; %AoD spread mean [deg]
                V = 2.7; %AoD spread std [deg]
                paraSt.mu_phi_BS = 10^(log10(exp(1))*(log(E)-1/2*log(1+V^2/E^2) )); %AoD spread mean [deg]
                paraSt.sigma_phi_BS = 10*log10(exp(1))*sqrt(log(1+V^2/E^2)); %AoD spread dB scale std

                E = 4.5; %EoD spread mean [deg]
                V = 3.0; %EoD spread std [deg]
                paraSt.mu_theta_BS = 10^(log10(exp(1))*(log(E)-1/2*log(1+V^2/E^2) )); %EoD spread mean [deg]
                paraSt.sigma_theta_BS = 10*log10(exp(1))*sqrt(log(1+V^2/E^2)); %EoD dB scale spread

		        E=4; % AoA spread mean [deg]
        		V=2; % AoA spread std [deg]		
		        paraSt.mu_phi_MS = 10^(log10(exp(1))*(log(E)-1/2*log(1+V^2/E^2) )); % AoA spread mean [deg]      
		        paraSt.sigma_phi_MS = 10*log10(exp(1))*sqrt(log(1+V^2/E^2)); % AoA dB scale spread std 

		        E=1; % EoA spread mean [deg]
		        V=1; % EoA spread mean [deg]
		        paraSt.mu_theta_MS = 10^(log10(exp(1))*(log(E)-1/2*log(1+V^2/E^2) )); % EoA spread mean [deg]
		        paraSt.sigma_theta_MS = 10*log10(exp(1))*sqrt(log(1+V^2/E^2)); % EoA dB scale std       

                paraSt.sigma_sf = 2.7; % shadow fading std dB
                
                r_mpc_gain = .5; % 3 dB power decay in meters
                mpc_gain_width = sqrt(-1*r_mpc_gain^2/(2*log(0.7))); % width of MPC Gaussian gain function 
                paraSt.mpc_gain_mu = 10*log10(mpc_gain_width); % MPC-VR radius, lognormal, mu, Inf. deactivates feature.
                paraSt.mpc_gain_sigma = 0; % MPC-VR radius, lognormal, sigma.

                n_mpc_eff = 10; % # of effective MPCs
                paraSt.n_mpc = round(n_mpc_eff*paraSt.r_c^2/r_mpc_gain^2); % # of MPCs per cluster
              
                %Cross-correlation matrix
                paraSt.rho = [1.00   -0.45   -0.50   -0.56    0.00     0.00;... % shadow fading
                             -0.45    1.00    0.34    0.70    0.00     0.00;... % delay spread 
                             -0.50    0.34    1.00    0.50    0.00     0.00;... % theta BS
                             -0.55    0.70    0.50    1.00    0.00     0.00;... % phi BS
                              0.00    0.00    0.00    0.00    1.00     0.00;... % theta MS
                              0.00    0.00    0.00    0.00    0.00     1.00];   % phi MS
                paraSt.corr_mat = chol(paraSt.rho);    % Cholesky factorization

                % Cluster distribution & cluster link delay
                paraSt.phi_c = [0 106]; % mean/std azimuth of cluster seen from VR at MS side [deg]        
                paraSt.pdf_phi_c = 'norm'; % phi_c follows normal distribution
                paraSt.theta_c = [4 77]; % mean/std elevation of cluster seen from VR at MS side [deg]    
                paraSt.pdf_theta_c = 'norm';  % theta_c follows normal distribution
                paraSt.para_r_c = [0 16]; % Min/max distance from cluster to BS/MS
                paraSt.pdf_r_c = 'unif'; % r_c uniform distribution        
                paraSt.mu_tauCLink = 8.54e-7; % Cluster link delay mean [s]
                paraSt.min_tauCLink = 4.79e-8; % Cluster link delay dB scale std    
                
                % Polarization
                paraSt.mu_xpdv = 9; % [dB]
                paraSt.sigma_xpdv = 3; % Mean and std for XPDV
                paraSt.mu_xpdh = 9; % [dB]
                paraSt.sigma_xpdh = 3; % Mean and std for XPDH
                
                paraSt.mu_cpr = 0; % [dB]
                paraSt.sigma_cpr = 0; % Mean and std for CPR

                % For multi-link extension.
                paraSt.BS_common = 0; %BS-common cluster ratio to total number of clusters (1,2,3,...)
                paraSt.MS_common = 0; %Average number of VRs in VR group for one MS-common cluster

                %DMC parameters
                paraSt.dmc_spread_mean = 0; %dmc spatial spread [m]
                paraSt.dmc_spread_sigma = 0; %dmc spatial spread [m]
                paraSt.dmc_beta = 0; %dmc delay power decaying [dB]
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'SemiUrban_CloselySpacedUser_2_6GHz' % Massive MIMO with a compact array.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch paraEx.scenario
            case 'LOS'
                paraSt.bs_vr_mu = Inf; % BS-VR length, Inf means BS-VR feature has no effect.
                paraSt.bs_vr_sigma = 0; % BS-VR length.
                paraSt.bs_vr_slope_mu = 0; % BS-VR power slope, 0 means BS-VR feature has no effect.
                paraSt.bs_vr_slope_sigma = 0; % BS-VR power slope, normal, sigma
                
                paraSt.r_c = 10; % Visibility region size at MS [m]
                paraSt.l_c = 2; % Transition region size at MS [m]
 
                paraSt.k_tau = log(10^(79.6/10)); % Cluster power decaying factor lin. [/us]
                paraSt.tau_b = 1.7e-6; % Cluster power decaying cut-off delay [s] (cell size)

                paraSt.d_co = 350; % LOS cut-off distance [m]
                paraSt.r_l = 343; % LOS VR size at MS [m]
                paraSt.l_l = 93; % LOS TR size at MS [m]
                paraSt.r_l_bs = 200; % LOS VR size at BS [m]
                paraSt.mu_k = 10^(2.8/10); % LOS power factor, mean linear
                paraSt.sigma_k = 0.8; % LOS power factor, std dB
                paraSt.power_factor = paraSt.mu_k*10^(randn(1)*paraSt.sigma_k/10); % LOS power factor   

                paraSt.BSLocal = false; % BS Local cluster activity true/false
                paraSt.MSLocal = false; % MS local cluster activity true/false
               
                n_c_far = 15;
                paraSt.n_c_far = n_c_far; % Average number of far clusters                

                paraSt.k_sel = 1; % K-selection factor (ratio of the number of single bounce clusters relative to twin clusters)

                paraSt.mu_tau = 0.02e-6; % Delay spread mean [s]
                paraSt.sigma_tau = 0.01; % Delay spread dB scale std        

                paraSt.mu_phi_BS = 8.5; % AoD spread mean [deg]
                paraSt.sigma_phi_BS = 1.9; % AoD spread dB scale std        

                paraSt.mu_theta_BS = 7; % EoD spread mean [deg]
                paraSt.sigma_theta_BS = 1.9; % EoD dB scale spread

                paraSt.mu_phi_MS = 14.8; % AoA spread mean [deg]
                paraSt.sigma_phi_MS = 2.68; % AoA dB scale spread std 

                paraSt.mu_theta_MS = 10^0.6; % EoA spread mean [deg] - borrow from 3GPP 3D model
                paraSt.sigma_theta_MS = 0.16*10; % EoA dB scale std - borrow from 3GPP 3D model
                
                paraSt.sigma_sf = 5.8; % Shadow fading std dB
                
                r_mpc_gain = 3; % 3 dB power decay in meters
                mpc_gain_width = sqrt(-1*r_mpc_gain^2/(2*log(0.7))); % Width of MPC Gaussian gain function 
                paraSt.mpc_gain_mu = 10*log10(mpc_gain_width); % MPC-VR radius, lognormal, mu, Inf. deactivates feature.
                paraSt.mpc_gain_sigma = 0; % MPC-VR radius, lognormal, sigma.
                
                n_mpc_eff = 5; % Number of effective MPCs
                paraSt.n_mpc = round(n_mpc_eff*paraSt.r_c^2/r_mpc_gain^2); % Number of MPCs per cluster

              
                % Cross-correlation matrix
                paraSt.rho = [ 1    -0.5  -0.8  -0.8  0  0;...  % Shadow fading
                              -0.5   1     0.4   0.6  0  0;...  % Delay spread
                              -0.8   0.4   1     0.7  0  0;...  % theta BS
                              -0.8   0.6   0.7   1    0  0;...  % phi BS
                               0     0     0     0    1  0;...  % theta MS
                               0     0     0     0    0  1];    % phi MS   
                paraSt.corr_mat = chol(paraSt.rho); % Cholesky factorization

                % Cluster distribution and cluster link delay
                paraSt.phi_c = [43 39]; % Mean/std azimuth of cluster seen from VR at MS side [deg]        
                paraSt.pdf_phi_c = 'norm'; % phi_c follows normal distribution
                paraSt.theta_c = [10 5]; % Mean/std elevation of cluster seen from VR at MS side [deg]    
                paraSt.pdf_theta_c = 'norm';  % theta_c follows normal distribution
                paraSt.para_r_c = [0 paraEx.net_radii]; % Min/max distance from cluster to BS/MS        
                paraSt.pdf_r_c = 'unif'; % r_c uniform distribution        
                paraSt.mu_tauCLink = 8.54e-7; % Cluster link delay mean [s]
                paraSt.min_tauCLink = 4.79e-8; % Cluster link delay dB scale std    
                
                % Polarization
                paraSt.mu_xpdv = 9;
                paraSt.sigma_xpdv = 3; % Mean and std for XPDV
                paraSt.mu_xpdh = 9; 
                paraSt.sigma_xpdh = 3; % Mean and std for XPDH
                paraSt.mu_cpr = 0;
                paraSt.sigma_cpr = 0; % Mean and std for CPR

                % For multi-link extension
                paraSt.BS_common = 0; % BS-common cluster ratio to total number of clusters (1,2,3,...)
                paraSt.MS_common = 0; % Average number of VRs in VR group for one MS-common cluster

                % DMC parameters
                paraSt.dmc_spread_mean = 0; % DMC spatial spread [m]
                paraSt.dmc_spread_sigma = 0; % DMC spatial spread [m]
                paraSt.dmc_beta = 0; % DMC delay power decaying [dB]
            %%%%%%%%%%%%%%%%
            case 'NLOS'
            %%%%%%%%%%%%%%%%
                paraSt.bs_vr_mu = Inf; % BS-VR length, lognormal, mu
                paraSt.bs_vr_sigma = 0; % BS-VR length, lognormal, sigma                
                paraSt.bs_vr_slope_mu = 0; % BS-VR power slope, normal, mu
                paraSt.bs_vr_slope_sigma = 0; % BS-VR power slope, normal, sigma
                
                paraSt.r_c = 24.5; % Visibility region size at MS [m]
                paraSt.l_c = 12.2; % Transition region size at MS [m]
 
                paraSt.k_tau = log(10^(20/10)); % Cluster power decaying factor lin. [/us]
                paraSt.tau_b = 1.7e-6; % Cluster power decaying cut-off delay [s] (cell size)
               
                paraSt.d_co = 0; % LOS cut-off distance [m]
                paraSt.r_l = 0; % LOS VR size at MS [m]
                paraSt.l_l = 0; % LOS TR size at MS [m]
                paraSt.r_l_bs = 0; % LOS VR size at BS [m]
                paraSt.mu_k = 0; % LOS power factor, mean linear
                paraSt.sigma_k = 0; % LOS power factor, std dB
                paraSt.power_factor = paraSt.mu_k*10^(randn(1)*paraSt.sigma_k/10); % LOS power factor        

                paraSt.BSLocal = false; % BS Local cluster activity true/false
                paraSt.MSLocal = true; % MS local cluster activity true/false
               
                n_c_far = 14;
                paraSt.n_c_far = n_c_far; % Average number of far clusters                

                paraSt.k_sel = 0.8; % K-selection factor (ratio of the number of single bounce clusters relative to twin clusters)

                paraSt.mu_tau = 0.06e-6; % Delay spread mean [s]
                paraSt.sigma_tau = 0.01; % Delay spread dB scale std        

                paraSt.mu_phi_BS = 9.8; % AoD spread mean [deg]
                paraSt.sigma_phi_BS = 2.2; % AoD spread dB scale std        

                paraSt.mu_theta_BS = 8.9; % EoD spread mean [deg]
                paraSt.sigma_theta_BS = 1.9; % EoD dB scale spread

                paraSt.mu_phi_MS = 19; % AoA spread mean [deg]
                paraSt.sigma_phi_MS = 2.03; % AoA dB scale spread std 

                paraSt.mu_theta_MS = 10^0.88; % EoA spread mean [deg] - borrow from 3GPP 3D model
                paraSt.sigma_theta_MS = 0.16*10; % EoA dB scale std - borrow from 3GPP 3D model
                
                paraSt.sigma_sf = 5; % Shadow fading std dB

                paraSt.n_mpc = 58; % Number of MPCs per cluster
                mpc_gain_width = 1.2; % Width of MPC gain function 
                paraSt.mpc_gain_mu = 10*log10(mpc_gain_width); % MPC-VR radius, lognormal, mu, Inf. deactivates feature.
                paraSt.mpc_gain_sigma = 0; % MPC-VR radius, lognormal, sigma. Inf. deactivates feature.

                % Cross-correlation matrix
                paraSt.rho = [ 1    -0.4  -0.7  -0.8  0  0;...  % Shadow fading
                              -0.4   1     0.2   0.4  0  0;...  % Delay spread
                              -0.7   0.2   1     0.7  0  0;...  % theta BS
                              -0.8   0.4   0.7   1    0  0;...  % phi BS
                               0     0     0     0    1  0;...  % theta MS
                               0     0     0     0    0  1];    % phi MS   
                paraSt.corr_mat = chol(paraSt.rho); % Cholesky factorization

                % Cluster distribution and cluster link delay
                paraSt.phi_c = [43 39]; % Mean/std azimuth of cluster seen from VR at MS side [deg]        
                paraSt.pdf_phi_c = 'norm'; % phi_c follows normal distribution
                paraSt.theta_c = [10 5]; % Mean/std elevation of cluster seen from VR at MS side [deg]    
                paraSt.pdf_theta_c = 'norm';  % theta_c follows normal distribution
                paraSt.para_r_c = [0 paraEx.net_radii]; % Min/max distance from cluster to BS/MS        
                paraSt.pdf_r_c = 'unif'; % r_c uniform distribution        
                paraSt.mu_tauCLink = 1.02e-6; % Cluster link delay mean [s]
                paraSt.min_tauCLink = 5.24e-8; % Cluster link delay dB scale std     
                
                % Polarization
                paraSt.mu_xpdv = 8;
                paraSt.sigma_xpdv = 3; % Mean and std for XPDV
                paraSt.mu_xpdh = 8; 
                paraSt.sigma_xpdh = 3; % Mean and std for XPDH       
                paraSt.mu_cpr = 0;
                paraSt.sigma_cpr = 0; % Mean and std for CPR

                % For multi-link extension
                paraSt.BS_common = 0; % BS-common cluster ratio to total number of clusters (1,2,3,...)
                paraSt.MS_common = 0; % Average number of VRs in VR group for one MS-common cluster

                % DMC parameters
                paraSt.dmc_spread_mean = 0; % DMC spatial spread [m]
                paraSt.dmc_spread_sigma = 0; % DMC spatial spread [m]
                paraSt.dmc_beta = 0; % DMC delay power decaying [dB]
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'SemiUrban_VLA_2_6GHz' % Massive MIMO with physically large array at 2.6 GHz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         switch paraEx.scenario
            %%%%%%%%%%%%%%%%
            case 'LOS'
            %%%%%%%%%%%%%%%%
                paraSt.bs_vr_mu = 0.7; % BS-VR length, lognormal, mu
                paraSt.bs_vr_sigma = 2; % BS-VR length, lognormal, sigma                
                paraSt.bs_vr_slope_mu = 0; % BS-VR power slope, normal, mu
                paraSt.bs_vr_slope_sigma = 0.9; % BS-VR power slope, normal, sigma                
                
                paraSt.r_c = 24.5; % Visibility region size at MS [m]
                paraSt.l_c = 12.2; % Transition region size at MS [m]
                paraSt.k_tau = log(10^(27.83/10)); % Cluster power decaying factor lin. [/us]
                paraSt.tau_b = 0.87e-6; % Cluster power decaying cut-off delay [us] (cell size)

                paraSt.d_co = 350; % LOS cut-off distance [m]
                paraSt.r_l = 343; % LOS VR size at MS [m]
                paraSt.l_l = 93; % LOS TR size at MS [m]
                paraSt.r_l_bs = 7.35; % LOS VR size at BS [m]
                paraSt.mu_k = 10^(-5.19/10); % LOS power factor, mean linear
                paraSt.sigma_k = 3.47; % LOS power factor, std dB
                paraSt.power_factor = paraSt.mu_k*10^(randn(1)*paraSt.sigma_k/10); % LOS power factor        

                paraSt.BSLocal = false; % BS Local cluster activity true/false
                paraSt.MSLocal = false; % MS local cluster activity true/false
               
				n_c_far = 12.5;
				paraSt.n_c_far = n_c_far; % Average number of far clusters                

                paraSt.k_sel = 0.1; % K-selection factor (ratio of the number of single bounce clusters relative to twin clusters)

                paraSt.mu_tau = 0.15e-6; % Delay spread mean [s]
                paraSt.sigma_tau = 3.20; % Delay spread dB scale std        

                paraSt.mu_phi_BS =  11.04; % AoD spread mean [deg]
                paraSt.sigma_phi_BS = 2.93; % AoD spread dB scale std        

                paraSt.mu_theta_BS = 0; % EoD spread mean [deg]
                paraSt.sigma_theta_BS = 0; % EoD dB scale spread

                paraSt.mu_phi_MS =  0.2575/pi*180; % AoA spread mean [deg]      
                paraSt.sigma_phi_MS = 2.68; % AoA dB scale spread std 

                paraSt.mu_theta_MS = 0; % EoA spread mean [deg]
                paraSt.sigma_theta_MS = 0; % EoA dB scale std       
                
                paraSt.sigma_sf = 5.84; % Shadow fading std dB

                paraSt.n_mpc = 30; % Number of MPCs per cluster                
                paraSt.mpc_gain_mu = Inf; % MPC-VR radius, lognormal, mu, Inf. deactivates feature.
                paraSt.mpc_gain_sigma = 0; % MPC-VR radius, lognormal, sigma.

                % Cross-correlation matrix
                paraSt.rho = [1     0.35  0  0.09  0  0;...  % Shadow fading
                              0.35  1     0  0.27  0  0;...  % Delay spread
                              0     0     1  0     0  0;...  % theta BS
                              0.09  0.27  0  1     0  0;...  % phi BS
                              0     0     0  0     1  0;...  % theta MS
                              0     0     0  0     0  1];    % phi MS   
                paraSt.corr_mat = chol(paraSt.rho); % Cholesky factorization

                % Cluster distribution and cluster link delay
                paraSt.phi_c = [43 39]; % Mean/std azimuth of cluster seen from VR at MS side [deg]        
                paraSt.pdf_phi_c = 'norm'; % phi_c follows normal distribution
                paraSt.theta_c = [0 0]; % Mean/std elevation of cluster seen from VR at MS side [deg]    
                paraSt.pdf_theta_c = 'norm';  % theta_c follows normal distribution
                paraSt.para_r_c = [0 paraEx.net_radii]; % Min/max distance from cluster to BS/MS        
                paraSt.pdf_r_c = 'unif'; % r_c uniform distribution        
                paraSt.mu_tauCLink = 8.54e-7; % Cluster link delay mean [s]
                paraSt.min_tauCLink = 4.79e-8; % Cluster link delay dB scale std    
                
                % Polarization
                paraSt.mu_xpdv = 0;
                paraSt.sigma_xpdv = 0; % Mean and std for XPDV
                paraSt.mu_xpdh = 0; 
                paraSt.sigma_xpdh = 0; % Mean and std for XPDH
                paraSt.mu_cpr =0;
                paraSt.sigma_cpr = 0; % Mean and std for CPR

                % For multi-link extension
                if size(paraEx.pos_BS,1)>1
                    error('Several BSs with a very large array are not yet supported!');
                end
                paraSt.BS_common = 0; % BS-common cluster ratio to total number of clusters (1,2,3,...)
                paraSt.MS_common = 0; % Average number of VRs in VR group for one MS-common cluster

                % DMC parameters
                paraSt.dmc_spread_mean = 0; % DMC spatial spread [m]
                paraSt.dmc_spread_sigma = 0; % DMC spatial spread [m]
                paraSt.dmc_beta = 0; % DMC delay power decaying [dB]
            %%%%%%%%%%%%%%%%
            case 'NLOS'
            %%%%%%%%%%%%%%%%
                paraSt.bs_vr_mu = 0.7; % BS-VR length, lognormal, mu is mean of associated normal distribution
                paraSt.bs_vr_sigma = 2; % BS-VR length, lognormal, sigma is standard deviation of normal distribution                
                paraSt.bs_vr_slope_mu = 0; % BS-VR power slope, normal, mu
                paraSt.bs_vr_slope_sigma = 0.9; % BS-VR power slope, normal, sigma
                
                paraSt.r_c = 24.5; % Visibility region radius [m] - borrow from 300 MHz measurements
                paraSt.l_c = 12.2; % Transition region size [m] - borrow from 300 MHz measurements
                
                paraSt.k_tau = log(10^(43/10)); % Cluster power decaying factor lin. [/us]
                paraSt.tau_b = 0.91e-6; % Cluster power decaying cut-off delay [s] (cell size)

                paraSt.d_co = 0; % LOS cut-off distance [m]
                paraSt.r_l = 0; % LOS VR size [m]
                paraSt.l_l = 0; % LOS TR size [m]
                paraSt.mu_k = 10^(0/10); % LOS power factor, mean linear
                paraSt.sigma_k = 0; % LOS power factor, std dB
                paraSt.power_factor = paraSt.mu_k*10^(randn(1)*paraSt.sigma_k/10); % LOS power factor: 0 for NLOS       

                paraSt.BSLocal = false; % BS Local cluster activity true/false
                paraSt.MSLocal = false; % MS local cluster activity true/false
               
                % We have
                %
                %   n_c_far(L,Leff) = L*lambda + Leff*lambda,
                %
				% where L is the array length, Leff the mean BS-VR length, and lambda the BS-VR density.
                lambda = 2.3; % [BS-VR/m]
                paraSt.n_c_far = (paraEx.BS_range + paraSt.bs_vr_mu)*lambda; % Average number of far clusters
                
                paraSt.k_sel = 1; % K-selection factor (ratio of single cluster)

                paraSt.mu_tau = 0.14e-6; % Delay spread mean [s]
                paraSt.sigma_tau = 2.85; % Delay spread dB scale std        

                paraSt.mu_phi_BS = 6.96; % AoD spread mean [deg]
                paraSt.sigma_phi_BS = 2.39; % AoD spread dB scale std        

                paraSt.mu_theta_BS = 1.578/pi*180; % EoD spread mean [deg] - borrow from 300 MHz measurements
                paraSt.sigma_theta_BS = 0; % EoD dB scale spread - borrow from 300 MHz measurements

                paraSt.mu_phi_MS =  0.3318/pi*180; % AoA spread mean [deg] - borrow from 300 MHz measurements
                paraSt.sigma_phi_MS = 2.03; % 1.8599; % AoA dB scale spread std - borrow from 300 MHz measurements

                paraSt.mu_theta_MS = 1.578; % EoA spread mean [deg] - borrow from 300 MHz measurements
                paraSt.sigma_theta_MS = 0; % EoA dB scale std - borrow from 300 MHz measurements
                
                paraSt.sigma_sf = 7.55; % Shadow fading std dB

                paraSt.n_mpc = 31; % Number of MPCs per cluster
                paraSt.mpc_gain_mu = Inf; % MPC-VR radius, lognormal, mu, Inf. deactivates feature.
                paraSt.mpc_gain_sigma = 0; % MPC-VR radius, lognormal, sigma.

                % Cross-correlation matrix
                paraSt.rho =[ 1     -0.09  0  0.04  0  0;...  % Shadow fading
                             -0.09   1     0  0.42  0  0;...  % Delay spread
                              0      0     1  0     0  0;...  % theta BS
                              0.04   0.42  0  1     0  0;...  % phi BS
                              0      0     0  0     1  0;...  % theta MS
                              0      0     0  0     0  1];    % phi MS   
                paraSt.corr_mat=chol(paraSt.rho); % Cholesky factorization

                % Cluster distribution and cluster link delay
                paraSt.phi_c = [43 39]; % Mean/std azimuth of cluster to VR, mean [deg]        
                paraSt.pdf_phi_c = 'norm'; % phi_c normal distribution
                paraSt.theta_c = [0 0]; % Mean/std elevation of cluster to visibility region [deg]    
                paraSt.pdf_theta_c = 'norm'; % theta_c normal distribution
                paraSt.para_r_c = [0 paraEx.net_radii]; % Min/max distance from cluster to BS/MS        
                paraSt.pdf_r_c = 'unif'; % r_c uniform distribution        
                paraSt.mu_tauCLink = 1.02e-6; % Cluster link delay mean [s], for twin clusters - borrow from 300 MHz measurements
                paraSt.min_tauCLink = 5.24e-8; % Cluster link delay dB scale std, for twin clusters - borrow from 300 MHz measurements   

                % Polarization
                paraSt.mu_xpdv = 0;
                paraSt.sigma_xpdv = 0; % Mean and std for XPDV
                paraSt.mu_xpdh = 0; 
                paraSt.sigma_xpdh = 0; % Mean and std for XPDH
                paraSt.mu_cpr =0;
                paraSt.sigma_cpr = 0; % Mean and std for CPR
                
                % For multi-link extension
                if size(paraEx.pos_BS,1)>1
                    error('Several BS with a very-large array not yet supported!');
                end
                paraSt.BS_common = 0; % BS-common cluster ratio to total number of clusters (1,2,3,...)
                paraSt.MS_common = 0; % Average number of VRs in VR group for one MS-common cluster
              
                % DMC parameters
                paraSt.dmc_spread_mean = 0; % DMC spatial spread [m]
                paraSt.dmc_spread_sigma = 0; % DMC spatial spread [m]
                paraSt.dmc_beta = 0; % DMC delay power decaying [dB]
         end  
    otherwise
        disp('Specific network cannot be found!');
end
end