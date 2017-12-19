function h_omni_MIMO = create_IR_omni_MIMO(link,freq,delta_f,Band)
% link: simulated results from the COST 2100 channel model
% freq: start frequency and end frequency
% delta_f: the difference between two frequency bins
% Band: 'Wideband' or 'Narrowband'
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