function h_omni = create_IR_omni(link, freq, delta_f, Band)
%CREATE_IR_OMNI Get the SISO channel impulse reponse matrix
%
%Default call:
%h_omni = create_IR_omni(link, freq, delta_f, Band)
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
%h_omni: SISO channel impulse response matrix 

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
        first_snapshot = 1;
        last_snapshot = size(link.channel,2);
        bandwidth = freq(2)-freq(1);
        H = [];
        
        for jj = first_snapshot:last_snapshot  
            channel1 = link.channel{jj}.h; 
            channel2 = link.channel{jj}.h_los;
            channel = [channel1; channel2];
            MPC_delay = channel(:,5);    
            MPC_amp = channel(:,6); 
            Hf=[];
            Nfreq = freq(1):delta_f:freq(2);    
            
            for I = 1:length(Nfreq)
                Hf(I) = sum((MPC_amp).*exp(-1j*2*pi*Nfreq(I)*MPC_delay));        
            end
            
            H(jj,:) = Hf;
        end
        
        h_omni = ifft(H,[],2);
    case 'Narrowband'
        first_snapshot = 1;
        last_snapshot = size(link.channel,2);        
        h = zeros(1,last_snapshot);
        
        for jj = first_snapshot:last_snapshot 
            channel1 = link.channel{jj}.h; 
            channel2 = link.channel{jj}.h_los;
            channel = [channel1; channel2];
            MPC_amp = channel(:,6); 
            h(jj) = sum(MPC_amp);            
        end
        
        h_omni = h;
end