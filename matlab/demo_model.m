%% Demo model to run the COST 2100 channel model
% This is an example file of how to run the model. The input parameters, as
% stated below, are first chosen before the main function cost2100 run. The
% result from cost2100 is combined with different antenna patterns, with
% the available options stated as Output below.
%
% For more information about the available networks, references and version
% history, see the file 'Readme.txt'. If you use the COST 2100 channel model
% for publications, please refer to the stated publications.
%
%------
%Input:
%------
% Network : 'IndoorHall_5GHz','SemiUrban_300MHz','Indoor_CloselySpacedUser_2_6GHz','SemiUrban_CloselySpacedUser_2_6GHz', or 'SemiUrban_VLA_2_6GHz'
% Band : 'Wideband' or 'Narrowband'
% Link: 'Multiple' or 'Single'
% Antenna: 'SISO_omni', 'MIMO_omni', 'MIMO_dipole', 'MIMO_measured', 'MIMO_Cyl_patch', 'MIMO_VLA_omni'
% scenario: 'LOS' or 'NLOS'        
% freq: Frequency band [Hz]
% snapRate: Number of snapshots per s
% snapNum: Number of simulated snapshots         
% BSPosCenter: Center position of BS array [x, y, z] [m]
% BSPosSpacing: Inter-position spacing [m], for large arrays
% BSPosNum: Number of positions at each BS site, for large arrays
% MSPos: Position of MSs [m]
% MSVelo: Velocity of MSs [m/s]
%
%------
%Output:
%------ 
% 1) SISO_omni: Transfer function for SISO omni-directional antenna
% create_IR_omni: users have to set up the frequency separation, delta_f
% 
% 2) MIMO_omni: Transfer function for MIMO omini-directional antenna
% create_IR_omni_MIMO: users have to set up the frequency separation, delta_f.
% Only 2 by 2 MIMO system is implemented.
% 
% 3) MIMO_dipole: Transfer function for a theoretical antenna response for 
% any size of lambda/2-spaced linear dipole antenna arrays. An Ntx-by-Nrx theoretical 
% antenna array response is generated and the correponding 
% channel transfer function is simulated.
% 
% 4) MIMO_measured: Transfer function for any measured MIMO antenna response
% get_H: users have to provide the full antenna response at the BS and 
% MS sides, and also the rotation of the antenna arrays. The antenna 
% response mat file have to be the same format as 'antSample.mat' file.
%
% 5) MIMO_Cyl_patch: Transfer function for a synthetic pattern of a cylindrical array 
% with 128 antennas.
% get_IR_Cyl_patch: users have to set up the frequency separation, delta_f, as well as
% provide the full antennas response at the BS and MS sides. The antennas response mat 
% file have to be the same format as the 'BS_Cyl_AntPattern.mat' and 'MS_AntPattern_User.mat' 
% files.

% 6) MIMO_VLA_omni: Transfer function for a physically large array with 128 omni-directional 
% antennas, with lambda/2 inter-element separation, and MS with omni-directional antenna.
% create_IR_omni_MIMO_VLA: users have to set up the frequency separation, delta_f

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

% Choose a Network type out of 
% {'IndoorHall_5GHz','SemiUrban_300MHz','Indoor_CloselySpacedUser_2_6GHz','SemiUrban_CloselySpacedUser_2_6GHz','SemiUrban_VLA_2_6GHz'}
% to parameterize the COST2100 model
Network = 'Indoor_CloselySpacedUser_2_6GHz';
% In COST2100, # links = # BSs x # MSs
% Set Link type to `Multiple' if you work with more than one link
% Set Link type to `Single' otherwise
Link = 'Multiple';
% Choose an Antenna type out of
% {'SISO_omni', 'MIMO_omni', 'MIMO_dipole', 'MIMO_measured', 'MIMO_Cyl_patch', 'MIMO_VLA_omni'}
Antenna = 'MIMO_Cyl_patch';
% ...and type of channel: {'Wideband','Narrowband'}.
Band = 'Wideband';

% Here are some tested combinations of the above variables:
% 'IndoorHall_5GHz', 'Single', 'SISO_omni', 'Wideband'
% 'SemiUrban_300MHz', 'Single', 'SISO_omni', 'Wideband'
% 'SemiUrban_300MHz', 'Multiple', 'MIMO_omni', 'Wideband'
% 'Indoor_CloselySpacedUser_2_6GHz', 'Multiple', 'MIMO_Cyl_patch', 'Wideband'
% 'SemiUrban_CloselySpacedUser_2_6GHz', 'Multiple', 'MIMO_Cyl_patch', 'Wideband'
% 'SemiUrban_VLA_2_6GHz', 'Single', 'MIMO_VLA_omni', 'Wideband'
% 'SemiUrban_VLA_2_6GHz', 'Multiple', 'MIMO_VLA_omni', 'Wideband'

switch Network
    %%%%%%%%%%%%%%%%%%%%%%
    case 'IndoorHall_5GHz'
    %%%%%%%%%%%%%%%%%%%%%%
        switch Link
            case 'Single'
                scenario = 'LOS'; % {'LOS'} only LOS is available
                freq = [-10e6 10e6]+5.3e9; % [Hz}
                snapRate = 1; % Number of snapshots per s
                snapNum = 100; % Number of snapshots
                MSPos  = [10 5  0]; % [m]
                MSVelo = -[0 .1 0]; % [m/s]
                BSPosCenter  = [10 10 0]; % Center position of BS array [x, y, z] [m]
                BSPosSpacing = [0 0 0]; % Inter-position spacing (m), for large arrays
                BSPosNum = 1; % Number of positions at each BS site, for large arrays
            case 'Multiple'
                error('IndoorHall_5GHz does not support multiple links.');
        end
    %%%%%%%%%%%%%%%%%%%%%%%
    case 'SemiUrban_300MHz'
    %%%%%%%%%%%%%%%%%%%%%%%
        switch Link
           case 'Single'
                scenario = 'LOS'; % {'LOS', 'NLOS'} is available
                freq = [2.75e8 2.95e8]; %  [Hz]
                snapRate = 1; % Number of snapshots per s
                snapNum = 100; % Number of snapshots        
                MSPos  = [100 -200 0]; % [m]
                MSVelo = [-0.2 0.9 0]; % [m/s]
                BSPosCenter  = [0 0 0]; % Center position of BS array [x, y, z] [m]
                BSPosSpacing = [0 0 0]; % Inter-position spacing (m), for large arrays
                BSPosNum = 1; % Number of positions at each BS site, for large arrays
           case 'Multiple'
                scenario = 'LOS'; % {'LOS'} only LOS is available
                freq = [2.75e8 2.95e8]; % [Hz]
                snapRate = 1; % Number of snapshots per s
                snapNum = 100; % Number of snapshots        
                MSPos  = [100 -200 0;
                          120 -200 0]; % [m]
                MSVelo = [-0.2 0.9 0;
                          -0.2 0.9 0]; % [m/s]
                BSPosCenter  = [0 0 0]; % Cnter position of BS array [x, y, z] [m]
                BSPosSpacing = [0 0 0]; % Inter-position spacing (m), for large arrays
                BSPosNum = 1; % Number of positions at each BS site, for large arrays
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Indoor_CloselySpacedUser_2_6GHz'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        scenario = 'LOS';
        freq = [2.58e9 2.62e9]; % starting freq. - ending freq. [Hz]
        snapNum = 50; % number of snapshots (given MSVelo, cover about .25 m)
        snapRate = 50; % number of snapshots per second (sample at 0.05 m/snapshot)

        % closely-spaced users
        MSPos  = [   -2.5600    1.7300    2.2300;...
                     -3.0800    1.7300    2.2300;...
                     -2.5600    2.6200    2.5800;...
                     -4.6400    1.7300    2.2300;...
                     -2.5600    4.4000    3.3000;...
                     -3.0800    3.5100    2.9400;...
                     -3.6000    4.4000    3.3000;...
                     -4.1200    4.4000    3.3000;...
                     -4.1200    2.6200    2.5800]; % [x, y, z] (m)

        MSVelo = repmat([-.25,0,0],9,1); % [x, y, z] (m/s)

        BSPosCenter  = [0.30 -4.37 3.20]; % center position of BS array [x, y, z] (m)
        BSPosSpacing = [0 0 0]; % inter-position spacing (m), for large arrays.
        BSPosNum = 1; % number of positions at each BS site, for large arrays.
        
        BSPosCenter = BSPosCenter - mean(MSPos); % center users a origo
        MSPos = MSPos - repmat(mean(MSPos),size(MSPos,1),1); % center users a origo
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'SemiUrban_CloselySpacedUser_2_6GHz'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        scenario = 'LOS'; % {'LOS', 'NLOS'}
        freq = [2.58e9 2.62e9]; % Starting freq. - ending freq. [Hz]
        snapRate = 10; % Number of snapshots per second (sample at 0.05 m/snapshot)
        snapNum = 5; % Number of snapshots (given MSVelo, cover about .25 m)

        % Closely-spaced users
        MSPos  = [-27, -10, 0;...
                  -27-1.5, -10+1.5, 0;...
                  -27, -10+1.5, 0;...
                  -27+1.5, -10+1.5, 0;...
                  -27-1.5, -10, 0;...
                  -27+1.5, -10, 0;...
                  -27-1.5, -10-1.5, 0;...
                  -27, -10-1.5, 0;...
                  -27+1.5, -10-1.5, 0]; % [x, y, z] [m]

%         % Well-separated users
%         MSPos  = [27, 10, 0;...
%                   -27, 10, 0;...
%                   -27, -10, 0;...
%                   27, -10, 0;...
%                   10, 27, 0;...
%                   -10, 27, 0;...
%                   -10, -27, 0;...
%                   10, -27, 0;...
%                   5, 20, 0]; % [x, y, z] (m)

        MSVelo = [0.4, 0.3, 0;...
                  0.3, -0.4, 0;...
                  -0.5, 0.1, 0;...
                  -0.3, -0.4, 0;...
                  -0.4, -0.2, 0;...
                  0.3,  0.4, 0;...
                  0.3, -0.3, 0;...
                  -0.45,  0.1, 0;...
                  0.4, -0.3, 0]; % [x, y, z] [m/s]

        BSPosCenter  = [0 0 8]; % Center position of BS array [x, y, z] [m]
        BSPosSpacing = [0 0 0]; % Inter-position spacing (m), for large arrays
        BSPosNum = 1; % Number of positions at each BS site, for large arrays     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'SemiUrban_VLA_2_6GHz'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
        scenario = 'LOS'; % {'LOS', 'NLOS'}
        freq = [2.57e9 2.62e9]; % [Hz]
        snapRate = 1; % Number of snapshots per second
        snapNum = 1; % Number of snapshots
        MSPos  = [30 30 0; -10 -10 0]; % [x, y, z] [m]
        MSVelo = [0 0 0; 0 0 0]; % [m/s]
        BSPosCenter  = [0 0 0]; % Center position of BS array [x, y, z] [m]
        BSPosSpacing = [0.0577 0 0]; % Inter-position spacing (m), for large arrays
        BSPosNum = 128; % Number of positions at each BS site, for large arrays
end

  
tic
%% Get the MPCs from the COST 2100 channel model
[...
    paraEx,...       % External parameters
    paraSt,...       % Stochastic parameters
    link,...         % Simulated propagation data for all links [nBs,nMs]
    env...           % Simulated environment (clusters, clusters' VRs, etc.)
] = cost2100...
(...
    Network,...      % Model environment
    scenario,...     % LOS or NLOS
    freq,...         % [starting freq., ending freq.]
    snapRate,...     % Number of snapshots per second
    snapNum,...      % Total # of snapshots
    BSPosCenter,...  % Center position of each BS
    BSPosSpacing,... % Position spacing for each BS (parameter for physically very-large arrays)
    BSPosNum,...     % Number of positions on each BS (parameter for physically very-large arrays)
    MSPos,...        % Position of each MS
    MSVelo...        % Velocity of MS movements
    );         
toc

%% Visualize the generated environment
if 1  
    switch Network
        case {'IndoorHall_5GHz','SemiUrban_300MHz'}   
             visual_channel(paraEx, paraSt, link, env);
        case {'SemiUrban_VLA_2_6GHz','SemiUrban_CloselySpacedUser_2_6GHz','Indoor_CloselySpacedUser_2_6GHz'}   
             visualize_channel_env(paraEx, paraSt, link, env); axis equal; view(2);
    end   
end

%% Ccombine propagation data with antenna patterns
% Construct the channel data
% The following is example code
% End users can write their own code

switch Antenna
    %%%%%%%%%%%%%%%%%%
    case 'SISO_omni' % SISO between one BS and one MS
    %%%%%%%%%%%%%%%%%%
        switch Link
            %%%%%%%%%%%%%%
            case 'Single'
            %%%%%%%%%%%%%%
                delta_f = (freq(2)-freq(1))/256;
                h_omni = create_IR_omni(link,freq,delta_f,Band);
                switch Band
                    case 'Wideband'
                        H_omni = fft(h_omni,[],2);                        
                        figure;
                        mesh((freq(1):delta_f:freq(2))*1e-6,1:size(H_omni,1),10*log10(abs(H_omni)));
                        xlabel('Frequency [MHz]')
                        ylabel('Snapshots')
                        zlabel('Power [dB]')
                        title('Frequency response for the SISO channel')

                        figure;
                        plot(1:snapNum, pow2db(mean(abs(H_omni).^2, 2)));
                        xlabel('Snapshots')
                        ylabel('Power [dB]')
                        title('Frequency response for the SISO channel')
                    case 'Narrowband'
                        figure,plot(1:size(h_omni,2),10*log10(abs(h_omni)))                        
                        xlabel('Snapshots')
                        ylabel('Power [dB]')
                        title('Impulse response for the SISO channel')
                end
            %%%%%%%%%%%%%%%%
            case 'Multiple'
            %%%%%%%%%%%%%%%%
                delta_f = 7.8125e4;
                h_omni_1 = create_IR_omni(link(1),freq,delta_f,Band);
                h_omni_2 = create_IR_omni(link(2),freq,delta_f,Band);
                
                switch Band
                    case 'Wideband'
                        H_omni_1 = fft(h_omni_1,[],2);
                        H_omni_2 = fft(h_omni_2,[],2);
                        figure,mesh((freq(1):delta_f:freq(2))*1e-6,1:size(H_omni_1,1),10*log10(abs(H_omni_1)))
                        xlabel('Frequency [MHz]')
                        ylabel('Snapshots')
                        zlabel('Power [dB]')       
                        title('Frequency response for the SISO channel in simulated link 1')
                        figure,mesh((freq(1):delta_f:freq(2))*1e-6,1:size(H_omni_2,1),10*log10(abs(H_omni_2)))
                        xlabel('Frequency [MHz]')
                        ylabel('Snapshots')
                        zlabel('Power [dB]')
                        title('Frequency response for the SISO channel in simulated link 2') 
                    case 'Narrowband'
                        figure,plot(1:size(h_omni_1,2),10*log10(abs(h_omni_1)))                        
                        xlabel('Snapshots')
                        ylabel('Power [dB]')       
                        title('Impulse response for the SISO channel in simulated link 1')
                        figure,plot(1:size(h_omni_2,2),10*log10(abs(h_omni_2)))                        
                        xlabel('Snapshots')
                        ylabel('Power [dB]')       
                        title('Impulse response for the SISO channel in simulated link 2')
                end
        end 
        
    %%%%%%%%%%%%%%%%%%
    case 'MIMO_omni' % MIMO omni-directional between one BS and one MS
    %%%%%%%%%%%%%%%%%%
        switch Link
            %%%%%%%%%%%%%
            case 'Single'
            %%%%%%%%%%%%%%%%%
                % Channel transfer funtion with two omni-directional 
                % antennas at both BS and MS, distance between the two antennas is lambda/2
                delta_f = 7.8125e4;
                h_omni_MIMO = create_IR_omni_MIMO(link,freq,delta_f,Band);
                switch Band
                    case 'Wideband'
                        H_omni_MIMO = fft(h_omni_MIMO,[],2);
                        figure,mesh((freq(1):delta_f:freq(2))*1e-6,1:size(H_omni_MIMO(:,:,1,2),1),10*log10(abs(H_omni_MIMO(:,:,1,2))))
                        xlabel('Frequency [MHz]')
                        ylabel('Snapshots')
                        zlabel('Power [dB]')               
                        title('Frequency response for the channel between antenna 1 at Rx side and antenna 2 at Tx side ')
                    case 'Narrowband'
                        figure,plot(1:size(h_omni_MIMO,1),10*log10(abs(h_omni_MIMO(:,1,2))))                        
                        xlabel('Snapshots')
                        ylabel('Power [dB]')
                        title('Impulse response for the SISO channel')
                end
            %%%%%%%%%%%%%%%%
            case 'Multiple'
            %%%%%%%%%%%%%%%%
                % Channel transfer funtion with two omni-directional 
                % antennas at both BS and MS, distance between the two antennas is lambda/2
                delta_f = 7.8125e4;
                h_omni_MIMO_Link1 = create_IR_omni_MIMO(link(1),freq,delta_f,Band);
                h_omni_MIMO_Link2 = create_IR_omni_MIMO(link(2),freq,delta_f,Band);
                switch Band
                    case 'Wideband'
                        H_omni_MIMO_Link1 = fft(h_omni_MIMO_Link1,[],2);
                        H_omni_MIMO_Link2 = fft(h_omni_MIMO_Link2,[],2);
                        figure,mesh((freq(1):delta_f:freq(2))*1e-6,1:size(H_omni_MIMO_Link1(:,:,1,2),1),10*log10(abs(H_omni_MIMO_Link1(:,:,1,2))))
                        xlabel('Frequency [MHz]')
                        ylabel('Snapshots')
                        zlabel('Power [dB]')       
                        title('\fontsize{8} Frequency response for the channel between antenna 1 at Rx side and antenna 2 at Tx side in simulated link 1')
                        figure,mesh((freq(1):delta_f:freq(2))*1e-6,1:size(H_omni_MIMO_Link2(:,:,1,2),1),10*log10(abs(H_omni_MIMO_Link2(:,:,1,2))))
                        xlabel('Frequency [MHz]')
                        ylabel('Snapshots')
                        zlabel('Power [dB]')       
                        title('\fontsize{8} Frequency response for the channel between antenna 1 at Rx side and antenna 2 at Tx side in simulated link 2')
                    case 'Narrowband'
                        figure,plot(1:size(h_omni_MIMO_Link1(:,1,2),1),10*log10(abs(h_omni_MIMO_Link1(:,1,2))))                        
                        xlabel('Snapshots')
                        ylabel('Power [dB]')       
                        title('\fontsize{8} Impulse response for the channel between antenna 1 at Rx side and antenna 2 at Tx side in simulated link 1')
                        figure,plot(1:size(h_omni_MIMO_Link2(:,1,2),1),10*log10(abs(h_omni_MIMO_Link2(:,1,2))))                       
                        xlabel('Snapshots')
                        ylabel('Power [dB]')       
                        title('\fontsize{8} Impulse response for the channel between antenna 1 at Rx side and antenna 2 at Tx side in simulated link 2')
        
                end
        end 
    %%%%%%%%%%%%%%%%%%%%%
    case 'MIMO_dipole'  % MIMO theoretical dipole between one BS and one MS
    %%%%%%%%%%%%%%%%%%%%%
        switch Link
            %%%%%%%%%%%%%%
            case 'Single'
            %%%%%%%%%%%%%%
                % Channel transfer function for a theoretical dipole 2-by-2 linear antenna array
                % Generate the theoretical dipole antenna array, number of
                % antennas is changeable
                Ntx = 2;
                Nrx = 2;
                [Gtx, Grx] = get_dipole_G(Ntx,Nrx);        
                txRot = [0,0]; % Rotation of the antenna array [azi,ele]
                rxRot = [0,0];
                for Nsnap = 1:size(link.channel,2)
                    channel = link.channel{Nsnap};
                    % Get impulse response
                    h = get_H(channel, Gtx, Grx, txRot,rxRot,paraEx); % h[delay, rx, tx]                    
                    h_MIMO_dipole(Nsnap,:,:,:) = h;
                end
                switch Band
                    case 'Wideband'
                        H_MIMO_dipole = fft(h_MIMO_dipole,[],2);
                        figure,mesh(linspace(freq(1),freq(2),size(H_MIMO_dipole,2))*1e-6,1:size(H_MIMO_dipole(:,:,1,2),1),10*log10(abs(H_MIMO_dipole(:,:,1,2))))
                        xlabel('Frequency [MHz]')
                        ylabel('Snapshots')
                        zlabel('Power [dB]')               
                        title('Frequency response for the channel between antenna 1 at Rx side and antenna 2 at Tx side ')
                    case 'Narrowband'
                        figure,plot(1:size(h_MIMO_dipole(:,1,2),1),10*log10(abs(h_MIMO_dipole(:,1,2))))                        
                        xlabel('Snapshots')
                        ylabel('Power [dB]')               
                        title('Impulse response for the channel between antenna 1 at Rx side and antenna 2 at Tx side ')
                end
            %%%%%%%%%%%%%%%%
            case 'Multiple'
            %%%%%%%%%%%%%%%%
                % Channel transfer function for a theoretical dipole 2-by-2 linear antenna array
                % Generate the theoretical dipole antenna array, number of
                % antennas is changeable
                Ntx = 2;
                Nrx = 2;
                [Gtx, Grx] = get_dipole_G(Ntx,Nrx); 
                txRot = [0,0]; % Rotation of the antenna array [azi,ele]
                rxRot1 = [0,0];
                for Nsnap = 1:size(link(1).channel,2)
                    channel = link(1).channel{Nsnap};                    
                    h = get_H(channel, Gtx, Grx, txRot,rxRot1,paraEx); % h[delay, rx, tx]                   
                    h_MIMO_dipole_Link1(Nsnap,:,:,:) = h;
                end

                Ntx = 2;
                Nrx = 2;
                [Gtx, Grx] = get_dipole_G(Ntx,Nrx); 
                txRot = [0,0]; % Rotation of the antenna array [azi,ele]
                rxRot2 = [0,0];
                for Nsnap = 1:size(link(2).channel,2)
                    channel = link(2).channel{Nsnap};                    
                    h = get_H(channel, Gtx, Grx, txRot,rxRot2,paraEx); % h[delay, rx, tx]                    
                    h_MIMO_dipole_Link2(Nsnap,:,:,:) = h;
                end
                switch Band
                    case 'Wideband'
                        H_MIMO_dipole_Link1 = fft(h_MIMO_dipole_Link1,[],2);
                        H_MIMO_dipole_Link2 = fft(h_MIMO_dipole_Link2,[],2);
                        figure,mesh(linspace(freq(1),freq(2),size(H_MIMO_dipole_Link1,2))*1e-6,1:size(H_MIMO_dipole_Link1(:,:,1,2),1),10*log10(abs(H_MIMO_dipole_Link1(:,:,1,2))))
                        xlabel('Frequency [MHz]')
                        ylabel('Snapshots')
                        zlabel('Power [dB]')       
                        title('\fontsize{8} Frequency response for the channel between antenna 1 at Rx side and antenna 2 at Tx side in simulated link 1')
                        figure,mesh(linspace(freq(1),freq(2),size(H_MIMO_dipole_Link2,2))*1e-6,1:size(H_MIMO_dipole_Link2(:,:,1,2),1),10*log10(abs(H_MIMO_dipole_Link2(:,:,1,2))))
                        xlabel('Frequency [MHz]')
                        ylabel('Snapshots')
                        zlabel('Power [dB]')       
                        title('\fontsize{8} Frequency response for the channel between antenna 1 at Rx side and antenna 2 at Tx side in simulated link 2')
                    case 'Narrowband'
                        figure,plot(1:size(h_MIMO_dipole_Link1(:,1,2),1),10*log10(abs(h_MIMO_dipole_Link1(:,1,2))))                        
                        xlabel('Snapshots')
                        ylabel('Power [dB]')       
                        title('\fontsize{8} Impulse response for the channel between antenna 1 at Rx side and antenna 2 at Tx side in simulated link 1')
                        figure,plot(1:size(h_MIMO_dipole_Link2(:,1,2),1),10*log10(abs(h_MIMO_dipole_Link2(:,1,2))))                        
                        xlabel('Snapshots')
                        ylabel('Power [dB]')       
                        title('\fontsize{8} Impulse response for the channel between antenna 1 at Rx side and antenna 2 at Tx side in simulated link 2')
                end
        end 
    %%%%%%%%%%%%%%%%%%%%%%
    case 'MIMO_measured' % Measured antenna patterns
    %%%%%%%%%%%%%%%%%%%%%%
        switch Link
            %%%%%%%%%%%%%%
            case 'Single'
            %%%%%%%%%%%%%%
                % Channel transfer function for any measured MIMO antenna
                % response with same format as 'antSample.mat' is possible
                load antSample.mat
                txFull = antFull;
                rxFull = antFull;
                txRot = [0,0]; % Rotation of the antenna array [azi,ele]
                rxRot = [0,0];
                for Nsnap = 1:size(link.channel,2)
                    channel = link.channel{Nsnap};
                    % Get impulse response
                    h = get_H(channel, txFull, rxFull, txRot,rxRot,paraEx); % h[delay, rx, tx]                    
                    h_MIMO(Nsnap,:,:,:) = h;
                end
                switch Band
                    case 'Wideband'
                        H_MIMO = fft(h_MIMO,[],2);
                        figure,mesh(linspace(freq(1),freq(2),size(H_MIMO,2))*1e-6,1:size(H_MIMO(:,:,1,2),1),10*log10(abs(H_MIMO(:,:,1,2))))
                        xlabel('Frequency [MHz]')
                        ylabel('Snapshots')
                        zlabel('Power [dB]')               
                        title('Frequency response for the channel between antenna 1 at Rx side and antenna 2 at Tx side ')
                    case 'Narrowband'
                        figure,plot(1:size(h_MIMO(:,1,2),1),10*log10(abs(h_MIMO(:,1,2))))                        
                        xlabel('Snapshots')
                        ylabel('Power [dB]')               
                        title('Impulse response for the channel between antenna 1 at Rx side and antenna 2 at Tx side ')
                        
                end
            %%%%%%%%%%%%%%%%
            case 'Multiple'
            %%%%%%%%%%%%%%%%
                % Channel transfer function for any measured MIMO antenna
                % response with same format as 'antSample.mat' is possible
                load antSample.mat
                txFull = antFull;
                rxFull = antFull;       
                txRot = [0,0]; % Rotation of the antenna array [azi,ele]
                rxRot1 = [0,0];
                for Nsnap = 1:size(link(1).channel,2)
                    channel = link(1).channel{Nsnap};
                    % Get impulse response
                    h = get_H(channel, txFull, rxFull, txRot,rxRot1,paraEx); % h[delay, rx, tx]                    
                    h_MIMO_Link1(Nsnap,:,:,:) = h;
                end

                load antSample.mat
                txFull = antFull;
                rxFull = antFull;  
                txRot = [0,0]; % Rotation of the antenna array [azi,ele]
                rxRot2 = [0,0];
                for Nsnap = 1:size(link(2).channel,2)
                    channel = link(2).channel{Nsnap};
                    % Get impulse response
                    h = get_H(channel, txFull, rxFull, txRot,rxRot2,paraEx); % h[delay, rx, tx]                    
                    h_MIMO_Link2(Nsnap,:,:,:) = h;
                end
                switch Band
                    case 'Wideband'  
                        H_MIMO_Link1 = fft(h_MIMO_Link1,[],2);
                        H_MIMO_Link2 = fft(h_MIMO_Link2,[],2);
                        figure,mesh(linspace(freq(1),freq(2),size(H_MIMO_Link1,2))*1e-6,1:size(H_MIMO_Link1(:,:,1,2),1),10*log10(abs(H_MIMO_Link1(:,:,1,2))))
                        xlabel('Frequency [MHz]')
                        ylabel('Snapshots')
                        zlabel('Power [dB]')       
                        title('\fontsize{8} Frequency response for the channel between antenna 1 at Rx side and antenna 2 at Tx side in simulated link 1')

                        figure,mesh(linspace(freq(1),freq(2),size(H_MIMO_Link2,2))*1e-6,1:size(H_MIMO_Link2(:,:,1,2),1),10*log10(abs(H_MIMO_Link2(:,:,1,2))))
                        xlabel('Frequency [MHz]')
                        ylabel('Snapshots')
                        zlabel('Power [dB]')       
                        title('\fontsize{8} Frequency response for the channel between antenna 1 at Rx side and antenna 2 at Tx side in simulated link 2')
                    case 'Narrowband'
                        figure,plot(1:size(h_MIMO_Link1(:,1,2),1),10*log10(abs(h_MIMO_Link1(:,1,2))))
                        xlabel('Snapshots')
                        ylabel('Power [dB]')       
                        title('\fontsize{8} Impulse response for the channel between antenna 1 at Rx side and antenna 2 at Tx side in simulated link 1')

                        figure,plot(1:size(h_MIMO_Link2(:,1,2),1),10*log10(abs(h_MIMO_Link2(:,1,2))))
                        xlabel('Snapshots')
                        ylabel('Power [dB]')       
                        title('\fontsize{8} Impulse response for the channel between antenna 1 at Rx side and antenna 2 at Tx side in simulated link 2')

                end
        end 
    %%%%%%%%%%%%%%%%%%%%%%
    case 'MIMO_VLA_omni' % MIMO omni-directional for very-large arrays
    %%%%%%%%%%%%%%%%%%%%%%
        switch Link
            %%%%%%%%%%%%%%
            case 'Single'
            %%%%%%%%%%%%%%
                % Channel transfer function with 128 omni-directional at the
                % BS, with lambda/2 inter-element separation, and one MS
                % with one omni-directional antenna
                delta_f = (freq(2)-freq(1))/256;
                h_omni_MIMO = create_IR_omni_MIMO_VLA(link,freq,delta_f,Band);
                switch Band
                    case 'Wideband'
                        H_omni_MIMO = fft(h_omni_MIMO,[],2);
                        figure,mesh((freq(1):delta_f:freq(2))*1e-6,1:size(H_omni_MIMO(:,:,1,2),1),10*log10(abs(H_omni_MIMO(:,:,1,2))))
                        xlabel('Frequency [MHz]')
                        ylabel('Snapshots')
                        zlabel('Power [dB]')               
                        title('Frequency response for the channel between antenna 1 at Rx side and antenna 2 at Tx side ')
                    case 'Narrowband'
                        figure,plot(1:size(h_omni_MIMO,1),10*log10(abs(h_omni_MIMO(:,1,2))))                        
                        xlabel('Snapshots')
                        ylabel('Power [dB]')
                        title('Impulse response for the SISO channel')
                end
            %%%%%%%%%%%%%%%%
            case 'Multiple'
            %%%%%%%%%%%%%%%%
                % Channel transfer function with 128 omni-directional at the
                % BS, with lambda/2 inter-element separation, and two MSs
                % with one omni-directional antenna each
                delta_f = (freq(2)-freq(1))/256;
                h_omni_MIMO_Link1 = create_IR_omni_MIMO_VLA(link(1),freq,delta_f,Band);
                h_omni_MIMO_Link2 = create_IR_omni_MIMO_VLA(link(2),freq,delta_f,Band);
                switch Band
                    case 'Wideband'
                        H_omni_MIMO_Link1 = fft(h_omni_MIMO_Link1,[],2);
                        H_omni_MIMO_Link2 = fft(h_omni_MIMO_Link2,[],2);
                        
                        figure;
                        subplot(1,2,1);
                        mesh(1:BSPosNum, (freq(1):delta_f:freq(2))*1e-6, log10(abs(squeeze(H_omni_MIMO_Link1(1, :, 1, :))).^2));
                        xlabel('Base station antennas')
                        ylabel('Frequency [MHz]')
                        zlabel('Power [dB]')  
                        title('User #1, Snapshot #1');
                        axis square;
                        subplot(1,2,2);
                        mesh(1:BSPosNum, (freq(1):delta_f:freq(2))*1e-6, log10(abs(squeeze(H_omni_MIMO_Link2(1, :, 1, :))).^2));
                        xlabel('Base station antennas')
                        ylabel('Frequency [MHz]')
                        zlabel('Power [dB]')  
                        title('User #2, Snapshot #1');
                        axis square;

                        figure;
                        subplot(1,2,1);
                        y1 = pow2db(squeeze(mean(mean(abs(H_omni_MIMO_Link1).^2, 1), 2)));
                        plot(1:BSPosNum, y1);
                        xlabel('Base station antennas')
                        ylabel('Power [dB]')
                        title('User #1, Snapshot #1');
                        axis tight;
                        subplot(1,2,2);
                        y2 = pow2db(squeeze(mean(mean(abs(H_omni_MIMO_Link2).^2, 1), 2)));
                        plot(1:BSPosNum, y2);
                        xlabel('Base station antennas')
                        ylabel('Power [dB]')
                        title('User #2, Snapshot #1');
                        axis tight;

                        figure;
                        subplot(1,2,1);
                        y1 = pow2db(squeeze(mean(mean(abs(H_omni_MIMO_Link1).^2, 1), 4)));
                        plot((freq(1):delta_f:freq(2))*1e-6, y1);
                        xlabel('Frequency [MHz]')
                        ylabel('Power [dB]')
                        title('User #1, Snapshot #1');
                        axis tight;
                        subplot(1,2,2);
                        y2 = pow2db(squeeze(mean(mean(abs(H_omni_MIMO_Link2).^2, 1), 4)));
                        plot((freq(1):delta_f:freq(2))*1e-6, y2);
                        xlabel('Frequency [MHz]')
                        ylabel('Power [dB]')
                        title('User #2, Snapshot #1');
                        axis tight;
                    case 'Narrowband'
                        H_omni_MIMO_Link1 = fft(h_omni_MIMO_Link1,[],2);
                        H_omni_MIMO_Link2 = fft(h_omni_MIMO_Link2,[],2);
                        
                        figure;
                        subplot(1,2,1);
                        y1 = pow2db(squeeze(mean(mean(abs(H_omni_MIMO_Link1).^2, 1), 2)));
                        plot(1:BSPosNum, y1);
                        xlabel('Base station antennas')
                        ylabel('Power [dB]')
                        title('User #1, Snapshot #1');
                        axis tight;
                        subplot(1,2,2);
                        y2 = pow2db(squeeze(mean(mean(abs(H_omni_MIMO_Link2).^2, 1), 2)));
                        plot(1:BSPosNum, y2);
                        xlabel('Base station antennas')
                        ylabel('Power [dB]')
                        title('User #2, Snapshot #1');
                        axis tight;
                end
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%
    case 'MIMO_Cyl_patch' % 128-element cylindrical array at BS
    %%%%%%%%%%%%%%%%%%%%%%%   
        USE_EADF = 1; % Use EADF or antenna pattern
    
        if (USE_EADF)
            BSantEADF = load('BS_Cyl_EADF.mat','F'); % BS antenna EADF
        else
            BSantPattern = load('BS_Cyl_AntPattern.mat'); % BS antenna array pattern
        end
        
        Nbs_ant = 128; % Number of BS antennas
        Nms = size(MSPos, 1); % Number of MS
   
        MSantPattern = load('MS_AntPattern_User.mat'); % MS antenna pattern with user effect
        
        delta_f = (freq(2)-freq(1))/256;

        tic
        % Get impulse response [snapshot, delay, ms, bs ant]
        if (USE_EADF)
            ir_Cyl_Patch = create_IR_Cyl_EADF(link, freq, delta_f, BSantEADF.F, MSantPattern);
        else
            ir_Cyl_Patch = create_IR_Cyl(link, freq, delta_f, BSantPattern, MSantPattern);
        end

        % Get transfer function [snapshot, freq, ms, bs ant]
        H_transfer = fft(ir_Cyl_Patch, [], 2);
        
        toc

        % Test plot transfer function [snapshot, freq, ms ant, bs ant]
        figure, mesh(1:Nbs_ant, (freq(1):delta_f:freq(2))*1e-6, log10(abs(squeeze(H_transfer(1, :, 1, :))).^2));
        xlabel('Base station antennas')
        ylabel('Frequency [MHz]')
        zlabel('Power [dB]')  
        title('User #1, Snapshot # 1');

        figure, plot(1:Nbs_ant, pow2db(squeeze(mean(mean(abs(H_transfer).^2, 1), 2))));
        xlabel('Base station antennas')
        ylabel('Power [dB]')
        title('User channel power over the base station array');

        % Condition number
        H_norm = zeros(size(H_transfer));
        for idx_user = 1:size(H_transfer, 3)
            norm_pow = sum(sum(sum(abs(H_transfer(:, :, idx_user, :)).^2, 1), 2), 4)/(size(H_transfer, 1)*size(H_transfer, 2)*size(H_transfer, 4));
            H_norm(:, :, idx_user, :) = H_transfer(:, :, idx_user, :)./sqrt(norm_pow);
        end
        cond_num = zeros(size(H_transfer, 1), size(H_transfer, 2));
        sig_val_all = zeros(size(H_transfer, 1), size(H_transfer, 2), size(H_transfer, 3));
        cond_num_all = zeros(size(H_transfer, 1), size(H_transfer, 2), size(H_transfer, 3));
        for idx_ss = 1:size(H_transfer, 1)
            for idx_freq = 1:size(H_transfer, 2)
                sig_val = svd(squeeze(H_norm(idx_ss, idx_freq, :, :)));
                cond_num(idx_ss, idx_freq) = max(sig_val)./min(sig_val);

                sig_val_all(idx_ss, idx_freq, :) = sig_val;
                cond_num_all(idx_ss, idx_freq, :) = max(sig_val)./sig_val;
            end
        end
        figure, cdfplot(20*log10(cond_num(:)));
        xlabel('Condition number [dB]'), ylabel('CDF');

        figure,
        for idx_user = 1:size(H_transfer, 3)
            cdfplot(20*log10(reshape(cond_num_all(:, :, idx_user), [], 1)));
            hold on
        end
        hold off
        xlabel('Ratio of maximum eigenvalue and other eigenvalues [dB]'), ylabel('CDF');        
end

return;