function [paraEx paraSt] = get_para(network,scenario,Nlink,Band,freq,snapRate, snapNum, posBS,posMS,veloMS)
%GET_PARA Generate the external and stochastic parameters of the scenario
%Default call: [paraEx paraSt] = get_para(network,scenario,Nlink,freq,snapRate, snapNum,
%posBS,posMS,veloMS)
%
%------
%Input:
%------
%network: network type, fours scenarios are available: 
%'macrocell', 'microcell', 'picocell', 'aalto'
%freq: frequency range, 2x1 vector, [freq_start freq_stop] [Hz]
%snapRate: channel snapshot rate, [s]
%snapeNum: channel snapshot number
%BSPos: position matrix of BSs, size = (number of BS, [x y z])  [m]
%MSPos: initial position matrix of MSs, size = (number of MS, [x y z])  [m]
%MSVelo: velocity vector of MSs, size = (number of MS, [x y z])  [m/s]
%
%------
%Output:
%------
%paraEx: external parameters
% .network: network type
% .freq_start: start frequency [Hz]
% .freq_stop: stop frequency [Hz]
% .freq: carrier frequency [Hz]
% .net_radii: network cell radius [m]
% .sample_rate: delay sampling rate [s]
% .c0: wave speed [m/s]
% .pos_BS: BS position [m]
% .pos_MS: MS position [m]
% .velo_MS: MS moving speed vector [m]
% .nfloor: average number of floors between BS and MS( picocell, aalto)
% .snap_rate: channel snapshot rate [s]
% .snap_num: channel snapshot number 
% .bandwidth: bandwidth [Hz]
% .delay_max: maximum delay, 5 times of net_radii [s]
%paraSt: stochastic parameters
% .r_c: visibility region radius for clusters [m]
% .l_c: transition region radius for clusters [m]
% .k_tau: cluster power attenuation coefficient [dB/mu_s]
% .tau_b: cut-off delay for cluster power attenuation [s]
% .r_l: visibility region radius for LOS [m]
% .l_l: transition region radius for LOS [m]
% .mu_k: LOS power factor mean 
% .sigma_k: LOS power factor dB scale std
% .BSLocal: Activity of BS local cluster
% .MSLocal: Activity of MS local cluster
% .mu_n_c_far: average number of far clusters, mean
% .n_c_far: average number of far clusters in the scenario
% .k_sel: proportion of twin/single clusters
% .n_mpc: number of MPCs per cluster
% .mu_tau: delay spread mean [s]
% .sigma_stau: delay spread dB scale std
% .mu_phi_BS: AoD mean [deg]
% .sigma_phi_BS: AoD dB std
% .mu_theta_BS: EoD mean [deg]
% .sigma_theta_BS: EoD dB scale std
% .mu_phi_MS: AoA mean [deg]
% .sigma_phi_MS: AoA dB scale std
% .mu_theta_MS: EoA mean [deg]
% .sigma_theta_MS: EoA dB scale std
% .sigma_sh: std of shadowing [dB]
% .rho: correlation matrix
% .corr_mat: Cholesky factorized correlation matrix
% .phi_c: mean/std azimuth of cluster to VR, mean [deg]  
% .pdf_phi_c: phi_c distribution
% .theta_c: mean/std elevation of cluster to VR, mean [deg] 
% .pdf_theta_c: theta_c distribution
% .r_c: statistics distance cluster-BS/VR [m]
% .pdf_r_c: r_c distribution
% .mu__tauCLink: cluster link delay extra, mean [s]
% .sigma__tauCLink: cluster link delay extra, dB scale std
% .BS_common: Ratio of BS-common cluster to total cluster ([2,3,...]BS-CC)
% .MS_common: average number of VRs in VR group
% .dmc_*: parameters of dmc
%
%See also: cost2100

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (C)2008 LIU Ling-Feng, ICTEAM, UCL, Belgium 
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

paraEx.network = network; %Network type
paraEx.scenario = scenario; %Scenario type
paraEx.Nlink = Nlink;% Multiple link or single link
paraEx.band = Band; % Bandwidth
paraEx.overSample = 4; % Oversampling factor for the impulse response shaping
switch network
    case 'SemiUrban_300MHz'
        paraEx.freq_start = 2.75e8;%[Hz] %Reference start frequency
        paraEx.freq_stop = 2.95e8;%[Hz] %Reference stop frequency
        paraEx.net_radii = 500; %cell radius [m]    
        paraEx.n_floor = 0; %Number of floors between BS and MS
    case 'IndoorHall_5GHz' %OLOS scenario
        paraEx.freq_start = 5.3e9-60e6;%[Hz]
        paraEx.freq_stop = 5.3e9+60e6;%[Hz]
        paraEx.net_radii = 30; %cell radius [m]    
        paraEx.n_floor = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Common external parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paraEx.snap_rate = snapRate; %Channel snapshot rate [s]
paraEx.snap_num = snapNum; %Channel snapshot number
paraEx.c0 = 3e8; %Wave speed [m/s]

paraEx.freq = (freq(1)+freq(2))/2; %Carrier frequency [Hz]
paraEx.bandwidth = freq(2)-freq(1); %Bandwidth [Hz]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BS/MS position & movement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paraEx.num_BS = length(posBS(:,1));
paraEx.num_MS = length(posMS(:,1));
paraEx.pos_BS = posBS;
paraEx.pos_MS = posMS;
paraEx.velo_MS = veloMS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Delay resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paraEx.sample_rate = 1/paraEx.bandwidth; %Delay sample rate [s]
%paraEx.delay_max = paraEx.net_radii/paraEx.c0; %Maximum delay 
paraEx.delay_max = paraEx.net_radii*5/paraEx.c0; %Maximum delay 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get ST parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paraSt = get_paraSt(paraEx);

end

function paraSt = get_paraSt(paraEx)
%Sub function to get the stochastic parameters

switch paraEx.network
   %% in the 300 MHz case two senario were extracted
    %%%%%%%%%%%%%%%%        
    case 'SemiUrban_300MHz' 
    %%%%%%%%%%%%%%%%   
        switch paraEx.scenario
            %%%%%%%%%%%%%%%%
            case 'LOS'
            %%%%%%%%%%%%%%%%             
                paraSt.r_c = 32.8; %Visibility region size [m]
                paraSt.l_c = 16.8; %Transition region size [m]
                paraSt.k_tau = log(10^(12.1/10));%Cluster power decaying factor lin. /us
                paraSt.tau_b = 2.4e-6;%Cluster power decaying cut-off delay us (cell size)

                paraSt.d_co = 350; %LOS cut-off distance [m]
                paraSt.r_l = 343; %LOS VR size [m]
                paraSt.l_l = 93; %LOS TR size [m]
                paraSt.mu_k = 10^(-4.7/10);%LOS power factor, mean linear
                paraSt.sigma_k = 2.0; %LOS power factor, std dB
                paraSt.power_factor = paraSt.mu_k*10^(randn(1)*paraSt.sigma_k/10); %LOS power factor        

                paraSt.BSLocal = false; %BS Local cluster activity true/false
                paraSt.MSLocal = true; %MS local cluster activity true/false
               
                n_c_far = 6;
                paraSt.n_c_far = n_c_far; %Average number of far clusters                

                paraSt.k_sel = 0.1; %K-selection factor (ratio of single cluster)

                paraSt.mu_tau = 1.39e-7; %Delay spread mean [s]
                paraSt.sigma_tau = 3.66; %Delay spread dB scale std        

                paraSt.mu_phi_BS =  0.255/pi*180; %AoD spread mean [deg]
                paraSt.sigma_phi_BS = 2.43;%1.9567; %AoD spread dB scale std        

                paraSt.mu_theta_BS = 1.578/pi*180; %EoD spread mean [deg]
                paraSt.sigma_theta_BS = 0; %EoD dB scale spread

                paraSt.mu_phi_MS =  0.2575/pi*180; %AoA spread mean [deg]      
                paraSt.sigma_phi_MS = 2.68; %AoA dB scale spread std 

                paraSt.mu_theta_MS = 1.578; %EoA spread mean [deg]
                paraSt.sigma_theta_MS = 0; %EoA dB scale std       
                paraSt.sigma_sf = 2.05; %shadow fading std dB

                paraSt.n_mpc = 27; %Number of MPCs per cluster

                %Cross-correlation matrix
                paraSt.rho =[1    0.0   0    0.0    0    0.0;...        % shadow fading
                            0.0   1    0    0.9    0    0.9;...        % delay spread
                             0    0    1    0    0    0;...        % theta BS
                            0.0    0.9  0    1    0    0.9;...        % phi BS
                             0    0    0    0    1    0;...        % theta MS
                            0.0    0.9    0    0.9    0    1];          % phi MS   
                paraSt.corr_mat=chol(paraSt.rho);    %Cholesky factorization

                %Cluster distribution & cluster link delay
                paraSt.phi_c = [43 39]; %mean/std azimuth of cluster to VR, mean [deg]        
                paraSt.pdf_phi_c = 'norm'; %phi_c normal distribution
                paraSt.theta_c = [0 0]; %mean/std elevation of cluster to visibility region [deg]    
                paraSt.pdf_theta_c = 'norm';  %theta_c normal distribution
                paraSt.para_r_c = [0 paraEx.net_radii]; %Min/max distance from cluster to BS/MS        
                paraSt.pdf_r_c = 'unif'; %r_c uniform distribution        
                paraSt.mu_tauCLink = 8.54e-7;%Cluster link delay mean [s]
                paraSt.min_tauCLink = 4.79e-8;% %Cluster link delay dB scale std        

                %Polarization
                paraSt.mu_xpdv = 0;
                paraSt.sigma_xpdv = 0; %Mean and std for XPDV
                paraSt.mu_xpdh = 0; 
                paraSt.sigma_xpdh = 0; %Mean and std for XPDH
                paraSt.mu_cpr =0;
                paraSt.sigma_cpr = 0; %Mean and std for CPR

                %For multi-link extension,
                if strcmp(paraEx.Nlink,'Multiple')
                    % more than 1 BS, parameters not extrated
                    if size(paraEx.pos_BS,1)>1
                        error('Parameters for more than 1 BS are not available')                        
                    end
                    if size(paraEx.pos_MS,1)> 1  
                        if size(paraEx.pos_MS,1) == 2 
                            paraSt.BS_common = 0; %BS-common cluster ratio to total number of clusters (1,2,3,...)
                            dist = sqrt(sum((paraEx.pos_MS(1,:) - paraEx.pos_MS(2,:)).^2));
                            % the modeling of common cluster decay ratio is
                            % not finish yet, so only round the dist to a
                            % certain range for now and get the
                            % corresponding common cluster ratio from the
                            % Eucap paper
                            if dist == 0
                                paraSt.MS_common = 1;%Average number of VRs in VR group for one MS-common cluster
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
                                paraSt.MS_common = 0.1; % when the distance is large than 40 meters, we just give a small common ratio.
                            end
                        else
                            error('Parameters for more than 2 MSs are not available')                            
                        end
                    end
                else
                    % signle link set up
                    paraSt.BS_common = 0;
                    paraSt.MS_common = 0;
                end

                %DMC parameters
                paraSt.dmc_spread_mean = 0; %dmc spatial spread [m]
                paraSt.dmc_spread_sigma = 0; %dmc spatial spread [m]
                paraSt.dmc_beta = 0; %dmc delay power decaying [dB]
            %%%%%%%%%%%%%%%%
            case 'NLOS'
            %%%%%%%%%%%%%%%%    
                paraSt.r_c = 24.5; %Visibility region size [m]
                paraSt.l_c = 12.2; %Transition region size [m]
                paraSt.k_tau = log(10^(7.2/10));%Cluster power decaying factor lin. /us
                paraSt.tau_b = 4.2e-6;%Cluster power decaying cut-off delay us (cell size)

                paraSt.d_co = 0; %LOS cut-off distance [m]
                paraSt.r_l = 0; %LOS VR size [m]
                paraSt.l_l = 0; %LOS TR size [m]
                paraSt.mu_k = 10^(0/10);%LOS power factor, mean linear
                paraSt.sigma_k = 0; %LOS power factor, std dB
                paraSt.power_factor = paraSt.mu_k*10^(randn(1)*paraSt.sigma_k/10); %LOS power factor        

                paraSt.BSLocal = false; %BS Local cluster activity true/false
                paraSt.MSLocal = true; %MS local cluster activity true/false
               
                n_c_far = 6;
                paraSt.n_c_far = n_c_far; %Average number of far clusters
                
                paraSt.k_sel = 0.2; %K-selection factor (ratio of single cluster)

                paraSt.mu_tau = 3.17e-7; %Delay spread mean [s]
                paraSt.sigma_tau = 2.5; %Delay spread dB scale std        

                paraSt.mu_phi_BS =  0.3254/pi*180; %AoD spread mean [deg]
                paraSt.sigma_phi_BS = 2.02;%1.9567; %AoD spread dB scale std        

                paraSt.mu_theta_BS = 1.578/pi*180; %EoD spread mean [deg]
                paraSt.sigma_theta_BS = 0; %EoD dB scale spread

                paraSt.mu_phi_MS =  0.3318/pi*180;%AoA spread mean [deg]      
                paraSt.sigma_phi_MS = 2.03;%1.8599; %AoA dB scale spread std 

                paraSt.mu_theta_MS = 1.578; %EoA spread mean [deg]
                paraSt.sigma_theta_MS = 0; %EoA dB scale std       
                paraSt.sigma_sf = 2.27; %shadow fading std dB

                paraSt.n_mpc = 48; %Number of MPCs per cluster

                %Cross-correlation matrix
                paraSt.rho =[1    -0.1   0    -0.1    0    -0.1;...        % shadow fading
                            -0.1     1    0    0.9    0    0.9;...        % delay spread
                             0    0    1    0    0    0;...        % theta BS
                            -0.1     0.9  0    1    0    0.9;...        % phi BS
                             0    0    0    0    1    0;...        % theta MS
                            -0.1    0.9    0    0.9    0    1];          % phi MS   
                paraSt.corr_mat=chol(paraSt.rho);    %Cholesky factorization

                %Cluster distribution & cluster link delay
                paraSt.phi_c = [43 39]; %mean/std azimuth of cluster to VR, mean [deg]        
                paraSt.pdf_phi_c = 'norm'; %phi_c normal distribution
                paraSt.theta_c = [0 0]; %mean/std elevation of cluster to visibility region [deg]    
                paraSt.pdf_theta_c = 'norm';  %theta_c normal distribution
                paraSt.para_r_c = [0 paraEx.net_radii]; %Min/max distance from cluster to BS/MS        
                paraSt.pdf_r_c = 'unif'; %r_c uniform distribution        
                paraSt.mu_tauCLink = 1.02e-6;%Cluster link delay mean [s]
                paraSt.min_tauCLink = 5.24e-8;% %Cluster link delay dB scale std        

                %Polarization
                paraSt.mu_xpdv = 0;
                paraSt.sigma_xpdv = 0; %Mean and std for XPDV
                paraSt.mu_xpdh = 0; 
                paraSt.sigma_xpdh = 0; %Mean and std for XPDH
                paraSt.mu_cpr =0;
                paraSt.sigma_cpr = 0; %Mean and std for CPR

                
                %For multi-link extension
                paraSt.BS_common = 0; %BS-common cluster ratio to total number of clusters (1,2,3,...)
                paraSt.MS_common = 0; %Average number of VRs in VR group for one MS-common cluster
              

                %DMC parameters
                paraSt.dmc_spread_mean = 0; %dmc spatial spread [m]
                paraSt.dmc_spread_sigma = 0; %dmc spatial spread [m]
                paraSt.dmc_beta = 0; %dmc delay power decaying [dB]
        end
    %%%%%%%%%%%%%%%%        
    case 'IndoorHall_5GHz' %Aalto Obstructive LOS measurements
    %%%%%%%%%%%%%%%%        
        paraSt.r_c = 2.72; %Visibility region size [m]
        paraSt.l_c = 1; %Transition region size [m]
        paraSt.k_tau = 16; %Cluster power decaying factor lin. /us
        paraSt.tau_b = 30/3e8*1e6; %Cluster power decaying cut-off delay us (cell size)
        
        paraSt.d_co = 30; %LOS cut-off distance [m]
        paraSt.r_l = 0; %LOS VR size [m]
        paraSt.l_l = 0; %LOS TR size [m]
        paraSt.mu_k = 10^(-10/10); %LOS power factor, mean linear
        paraSt.sigma_k = 7.2; %LOS power factor, std dB
        paraSt.power_factor = paraSt.mu_k*10^(randn(1)*paraSt.sigma_k/10); %LOS power factor        
                
        paraSt.BSLocal = true; %BS Local cluster activity true/false
        paraSt.MSLocal = true; %MS local cluster activity true/false
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        paraSt.mu_n_c_far = 3.8; %Average number of far clusters, mean
        n_c_far = poissrnd(paraSt.mu_n_c_far); %Possion distribution
        if n_c_far==0 n_c_far = 1; end %Guarantee min 1 far cluster
        paraSt.n_c_far = n_c_far; %Average number of far clusters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        paraSt.k_sel = 0; %K-selection factor (ratio of single cluster)
        
        E = 3e-9; %Delay spread mean [s]
        V = 0.6e-9; %Delay spread std [s]
        paraSt.mu_tau = 10^(log10(exp(1))*(log(E)-1/2*log(1+V^2/E^2) )); %Delay spread mean [s]
        paraSt.sigma_tau = 10*log10(exp(1))*sqrt(log(1+V^2/E^2)); %Delay spread dB scale std        
        E=3; %AoD spread mean [deg]
        V=2; %AoD spread std [deg]		
        paraSt.mu_phi_BS = 10^(log10(exp(1))*(log(E)-1/2*log(1+V^2/E^2) )); %AoD spread mean [deg]
        paraSt.sigma_phi_BS = 10*log10(exp(1))*sqrt(log(1+V^2/E^2)); %AoD spread dB scale std        
        E=1.5; %EoD spread mean [deg]
        V=1.5; %EoD spread std [deg]
        paraSt.mu_theta_BS = 10^(log10(exp(1))*(log(E)-1/2*log(1+V^2/E^2) )); %EoD spread mean [deg]
        paraSt.sigma_theta_BS = 10*log10(exp(1))*sqrt(log(1+V^2/E^2)); %EoD dB scale spread
        E=4; %AoA spread mean [deg]
        V=2; %AoA spread std [deg]		
        paraSt.mu_phi_MS = 10^(log10(exp(1))*(log(E)-1/2*log(1+V^2/E^2) )); %AoA spread mean [deg]      
        paraSt.sigma_phi_MS = 10*log10(exp(1))*sqrt(log(1+V^2/E^2)); %AoA dB scale spread std 
        E=1; %EoA spread mean [deg]
        V=1; %EoA spread mean [deg]
        paraSt.mu_theta_MS = 10^(log10(exp(1))*(log(E)-1/2*log(1+V^2/E^2) )); %EoA spread mean [deg]
        paraSt.sigma_theta_MS = 10*log10(exp(1))*sqrt(log(1+V^2/E^2)); %EoA dB scale std       
        paraSt.sigma_sf = 0; %shadow fading std dB
        		
        paraSt.n_mpc = 3; %Number of MPCs per cluster
        
        %Cross-correlation matrix
        paraSt.rho =[1    0    0    0    0    0;...        % shadow fading
                     0    1    0    0    0    0;...        % delay spread
                     0    0    1    0    0    0;...        % theta BS
                     0    0    0    1    0    0;...        % phi BS
                     0    0    0    0    1    0;...        % theta MS
                     0    0    0    0    0    1];          % phi MS   
        paraSt.corr_mat=chol(paraSt.rho);    %Cholesky factorization
               
        %Cluster distribution & cluster link delay
		paraSt.phi_c = [43 39]; %mean/std azimuth of cluster to VR, mean [deg]        
		paraSt.pdf_phi_c = 'norm'; %phi_c normal distribution
		paraSt.theta_c = [0 0]; %mean/std elevation of cluster to visibility region [deg]    
		paraSt.pdf_theta_c = 'norm';  %theta_c normal distribution
        paraSt.para_r_c = [0 30]; %Min/max distance from cluster to BS/MS        
		paraSt.pdf_r_c = 'unif'; %r_c uniform distribution        
        paraSt.mu_tauCLink = 1.07e-7; %Cluster link delay mean [s]
        paraSt.min_tauCLink = 3.36e-8; %Cluster link delay dB scale std        
        
        %Polarization
        paraSt.mu_xpdv = 0;
        paraSt.sigma_xpdv = 0; %Mean and std for XPDV
        paraSt.mu_xpdh = 0; 
        paraSt.sigma_xpdh = 0; %Mean and std for XPDH
        paraSt.mu_cpr =0;
        paraSt.sigma_cpr = 0; %Mean and std for CPR
        
        %For multi-link extension
        paraSt.BS_common = [0.35 0.65]; %BS-common cluster ratio to total number of clusters (1,2,3,...)
        paraSt.MS_common = 4; %Average number of VRs in VR group for one MS-common cluster
        
        %DMC parameters
        paraSt.dmc_spread_mean = 0; %dmc spatial spread [m]
        paraSt.dmc_spread_sigma = 0; %dmc spatial spread [m]
        paraSt.dmc_beta = 0; %dmc delay power decaying [dB]	
end
end