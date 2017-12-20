function h_omni_MIMO = create_IR_omni_MIMO(link, freq, delta_f, Band)
%CREATE_IR_OMNI_MIMO Get the MIMO channel impulse reponse matrix
%
%Default call:
%h_omni_MIMO = create_IR_omni_MIMO(link, freq, delta_f, Band)
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
%h_omni: MIMO channel impulse response matrix 

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
        Nt = 2;
        Nr = 2;
        first_snapshot = 1;
        last_snapshot = size(link.channel,2);
        Nfreq = freq(1):delta_f:freq(2);
        H = zeros(last_snapshot,length(Nfreq),Nr,Nt);
        
        for jj = first_snapshot:last_snapshot              
            channel1 = link.channel{jj}.h; 
            channel2 = link.channel{jj}.h_los;
            channel = [channel1 ;channel2];
            MPC_delay = channel(:,5);    
            MPC_amp = channel(:,6); 
            MPC_AOA = channel(:,3);
            MPC_AOD = channel(:,1);         
            H11 = zeros(1,length(Nfreq));
            H21 = zeros(1,length(Nfreq));
            H12 = zeros(1,length(Nfreq));
            H22 = zeros(1,length(Nfreq));
            
            for I = 1:length(Nfreq)
                H11(I) = sum((MPC_amp).*exp(-1j*2*pi*Nfreq(I)*MPC_delay));  
                H21(I) = sum((MPC_amp).*exp(-1j*2*pi*Nfreq(I)*MPC_delay).*exp(-1j*pi*cos(MPC_AOA)));
                H12(I) = sum((MPC_amp).*exp(-1j*2*pi*Nfreq(I)*MPC_delay).*exp(-1j*pi*cos(MPC_AOD)));
                H22(I) = sum((MPC_amp).*exp(-1j*2*pi*Nfreq(I)*MPC_delay).*exp(-1j*pi*cos(MPC_AOD)).*exp(-1j*pi*cos(MPC_AOA)));        
            end  
            
            H(jj,:,1,1) = H11;
            H(jj,:,2,1) = H21;
            H(jj,:,1,2) = H12;
            H(jj,:,2,2) = H22;        
        end

        h_omni_MIMO = ifft(H,[],2);
    case 'Narrowband'
        Nt = 2;
        Nr = 2;
        first_snapshot = 1;
        last_snapshot = size(link.channel,2);
        Nfreq = freq(1):delta_f:freq(2);
        H = zeros(last_snapshot,length(Nfreq),Nr,Nt);
        
        for jj = first_snapshot:last_snapshot              
            channel1 = link.channel{jj}.h; 
            channel2 = link.channel{jj}.h_los;
            channel = [channel1 ;channel2];              
            MPC_amp = channel(:,6); 
            MPC_AOA = channel(:,3);
            MPC_AOD = channel(:,1);            

            h11 = sum(MPC_amp);  
            h21 = sum((MPC_amp).*exp(-1j*pi*cos(MPC_AOA)));
            h12 = sum((MPC_amp).*exp(-1j*pi*cos(MPC_AOD)));
            h22 = sum((MPC_amp).*exp(-1j*pi*cos(MPC_AOD)).*exp(-1j*pi*cos(MPC_AOA)));        

            h(jj,1,1) = h11;
            h(jj,2,1) = h21;
            h(jj,1,2) = h12;
            h(jj,2,2) = h22;        
        end

        h_omni_MIMO = h;
end