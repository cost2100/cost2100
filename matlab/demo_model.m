% This version 2.3.3 is implemented by Meifang Zhu, Lund University, Sweden.
% SemiUrban_300MHz: both NLOS and LOS single link MIMO simulations 
% are supported. Outdoor LOS multiple  MIMO link simulation is supported. 
% IndoorHall_5GHz: OLOS single link MIMO simulation is supported.

% If you use the COST 2100 channel model for publications, please refer to: 
% L. Liu, J. Poutanen, F. Quitin, K. Haneda, F. Tufvesson, P. De Doncker,
% P. Vainikainen and C. Oestges, “The COST 2100 MIMO channel model,”
% IEEE Wireless Commun., vol 19, issue 6, pp 92-99, Dec. 2012.

% Further details about the COST 2100 channel model can be found in:
% Roberto Verdone (Editor), Alberto Zanella (Editor)
% Pervasive Mobile and Ambient Wireless Communications Pervasive Mobile 
% and Ambient Wireless Communications, ISBN 978-1-4471-2315-6,
% Springer, 2012. 

% If you use the SemiUrban_300MHz scenario, 
% further information can be found in:
% 1. Meifang Zhu, Gunnar Eriksson, and Fredrik Tufvesson, 
% "The COST 2100 Channel Model: Parameterization and Validation 
% Based on Outdoor MIMO Measurements at 300 MHz", 
% IEEE Transactions on Wireless Commun..
% 2. Meifang Zhu and Fredrik Tufvesson, "Virtual Multi-link Propagation 
% Investigation of an Outdoor Scenario At 300 MHz," Proceedings of the 
% 7th European Conference on Antennas and Propagation (EUCAP), Gothenburg,
% Sweden, April 2013.

% If you use the IndoorHall_5GHz scenario, 
% further information can be found in:
% 1. V. M. Kolmonen, P. Almers, J. Salmi, J. Koivunen, A. Richter,
% F. Tufvesson, A. Molisch, P. Vainikainen, "A dynamic dual-link 
% ideband MIMO channel sounder for 5.3 GHz," IEEE Transactions on 
% Instrumentation and Measurement, Vol. 59, No. 4, pp. 873-883, 2010.
% 2. J. Poutanen, K. Haneda, L. Lin, C. Oestges, F. Tufvesson , 
% P. Vainikainen, "Parameterization of the COST 2100 MIMO channel 
% modeling in indoor scenarios," Proceedings of the 5th European 
% Conference on Antennas and Propagation (EUCAP), Rome, Italy, 
% pp. 3606-3610, April 2011.


%------
%Input:
%------
% network : 'SemiUrban_300MHz' or 'IndoorHall_5GHz'
% Band : 'Wideband' or 'Narrowband'
% Link: 'Multiple' or 'Single'
% Antenna: 'SISO_omni', 'MIMO_omni', 'MIMO_dipole', 'MIMO_measured'
% scenario: 'LOS' or 'NLOS'        
% freq: Frequency band in Hz
% snapRate: Number of snapshots per s
% snapNum: Number of simulated snapshots         
% BSPos: Position of BS, unit m
% MSPos: Position of MSs, unit m
% MSVelo: Velocity of MSs,  unit m/s

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (C)2013 Meifang Zhu, Lund University, Sweden 
%This work is created within the COST 2100/COST IC1004 framework.
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

clear all; clc;
close all
network = 'IndoorHall_5GHz'   %{'IndoorHall_5GHz', 'SemiUrban_300MHz'}
Band = 'Wideband'              %{'Wideband','Narrowband'}
Link = 'Single'              %{'Single', 'Multiple'}
Antenna = 'SISO_omni'      %{'SISO_omni', 'MIMO_omni', 'MIMO_dipole', 'MIMO_measured'}

switch network
    %%%%%%%%%%%%%%%
    case 'SemiUrban_300MHz'
        %%%%%%%%%%%%%%%  
        switch Link
           case 'Single'
                scenario = 'LOS' % {'LOS', 'NLOS'} is possible
                % initialize frequency and snapshot rate
                freq = [2.75e8 2.95e8]; %Hz
                snapRate = 1; % number of snapshots per s
                snapNum = 100;         
                BSPos  = [0 0 0];
                MSPos  = [100 -200 0];%m
                MSVelo = [-0.2 0.9 0];%m/s
           case 'Multiple'
                scenario = 'LOS' % {'LOS'} only LOS is available
                % initialize frequency and snapshot rate
                freq = [2.75e8 2.95e8]; %Hz
                snapRate = 1; % number of snapshots per s
                snapNum = 100;         
                BSPos  = [0 0 0];           
                MSPos  = [100 -200 0;
                          120 -200 0];%m
                MSVelo = [-0.2 0.9 0;
                          -0.2 0.9 0];%m/s
        end
    %%%%%%%%%%%%%%%
    case 'IndoorHall_5GHz'
    %%%%%%%%%%%%%%%
        switch Link
            case 'Single'
                scenario = 'LOS' % {'LOS'} only LOS is available
                freq = [-10e6 10e6]+5.3e9; %Hz
                snapRate = 1; % number of snapshots per s
                snapNum = 100;
                BSPos  = [10   10  0;];
                MSPos  = [10   5   0];
                MSVelo = [0    0.001   0];
        end
end


%% get the MPCs from the COST 2100 channel model
[paraEx paraSt link env BS MS ] = cost2100(network, scenario, Link, Band, freq, snapRate, snapNum, BSPos,MSPos,MSVelo);

%% reconstruct the channel function
switch Link
    %% single link simulation
    case 'Single'        
        switch Antenna
            case 'SISO_omni'
                delta_f = 7.8125e4;
                h_omni = create_IR_omni(link,freq,delta_f,Band);                                         
                switch Band
                    case 'Wideband'
                        H_omni = fft(h_omni,[],2);                        
                        figure,mesh((freq(1):delta_f:freq(2))*1e-6,1:size(H_omni,1),10*log10(abs(H_omni)))
                        xlabel('Frequency [MHz]')
                        ylabel('Snapshots')
                        zlabel('Power [dB]')
                        title('Frequency response for the SISO channel')
                    case 'Narrowband'
                        % channel transfer function for SISO omni-directional
                        % antenna                       
                        h_omni = create_IR_omni(link,freq,delta_f,Band);
                        figure,plot(1:size(h_omni,2),10*log10(abs(h_omni)))                        
                        xlabel('Snapshots')
                        ylabel('Power [dB]')
                        title('Impulse response for the SISO channel')
                end
            case 'MIMO_omni'                
                % channel transfer funtion with two omni-directional 
                % antennas at both BS and MS, distance between the two antennas is \lambda/2
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
            case 'MIMO_dipole'   
                % channel transfer function for a theoretical dipole 2-by-2 linear antenna array
                % generate the theoretical dipole antenna array, number of antenna
                % is changeable
                Ntx = 2;
                Nrx = 2;
                [Gtx, Grx] = get_dipole_G(Ntx,Nrx);        
                txRot = [0,0];% rotation of the antenna array [azi,ele]
                rxRot = [0,0];
                for Nsnap = 1:size(link.channel,2)
                    channel = link.channel{Nsnap};
                    % get impulse response
                    h = get_H(channel, Gtx, Grx, txRot,rxRot,paraEx);% h[delay, rx, tx]                    
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
            case 'MIMO_measured'
                % channel transfer function for any measured MIMO antenna response with same format of 'antSample.mat'
                load antSample.mat
                txFull = antFull;
                rxFull = antFull;
                txRot = [0,0];% rotation of the antenna array [azi,ele]
                rxRot = [0,0];
                for Nsnap = 1:size(link.channel,2)
                    channel = link.channel{Nsnap};
                    % get impulse response
                    h = get_H(channel, txFull, rxFull, txRot,rxRot,paraEx);% h[delay, rx, tx]                    
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
        end    
    %% multiple link simulation
    case 'Multiple'
        switch Antenna
            case 'SISO_omni'                       
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
            case 'MIMO_omni'                                 
                % channel transfer funtion with two omni-directional 
                % antennas at both BS and MS, distance between the two antennas is \lambda/2
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
            case 'MIMO_dipole' 
                % channel transfer function for a theoretical dipole 2-by-2 linear antenna array
                % generate the theoretical dipole antenna array, number of antenna
                % is changeable
                Ntx = 2;
                Nrx = 2;
                [Gtx, Grx] = get_dipole_G(Ntx,Nrx); 
                txRot = [0,0];% rotation of the antenna array [azi,ele]
                rxRot1 = [0,0];
                for Nsnap = 1:size(link(1).channel,2)
                    channel = link(1).channel{Nsnap};                    
                    h = get_H(channel, Gtx, Grx, txRot,rxRot1,paraEx);% h[delay, rx, tx]                   
                    h_MIMO_dipole_Link1(Nsnap,:,:,:) = h;
                end

                Ntx = 2;
                Nrx = 2;
                [Gtx, Grx] = get_dipole_G(Ntx,Nrx); 
                txRot = [0,0];% rotation of the antenna array [azi,ele]
                rxRot2 = [0,0];
                for Nsnap = 1:size(link(2).channel,2)
                    channel = link(2).channel{Nsnap};                    
                    h = get_H(channel, Gtx, Grx, txRot,rxRot2,paraEx);% h[delay, rx, tx]                    
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

            case 'MIMO_measured'
                % channel transfer function for any measured MIMO antenna response with same format of 'antSample.mat'
                load antSample.mat
                txFull = antFull;
                rxFull = antFull;       
                txRot = [0,0];% rotation of the antenna array [azi,ele]
                rxRot1 = [0,0];
                for Nsnap = 1:size(link(1).channel,2)
                    channel = link(1).channel{Nsnap};
                    % get impulse response
                    h = get_H(channel, txFull, rxFull, txRot,rxRot1,paraEx);% h[delay, rx, tx]                    
                    h_MIMO_Link1(Nsnap,:,:,:) = h;
                end

                load antSample.mat
                txFull = antFull;
                rxFull = antFull;  
                txRot = [0,0];% rotation of the antenna array [azi,ele]
                rxRot2 = [0,0];
                for Nsnap = 1:size(link(2).channel,2)
                    channel = link(2).channel{Nsnap};
                    % get impulse response
                    h = get_H(channel, txFull, rxFull, txRot,rxRot2,paraEx);% h[delay, rx, tx]                    
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
end


