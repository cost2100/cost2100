function h_omni = create_IR_omni(link,freq,delta_f,Band)
% link: simulated results from the COST 2100 channel model
% freq: start frequency and end frequency
% delta_f: the difference between two frequency bins
% Band: 'Wideband' or 'Narrowband'
switch Band
    case 'Wideband'
        first_snapshot = 1;
        last_snapshot = size(link.channel,2);
        bandwidth = freq(2)-freq(1);
        H = [];
        for jj = first_snapshot:last_snapshot  
            channel1 = link.channel{jj}.h; 
            channel2 = link.channel{jj}.h_los;
            channel = [channel1 ;channel2];
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
            channel = [channel1 ;channel2];
            MPC_amp = channel(:,6); 
            h(jj) = sum(MPC_amp);            
        end
        h_omni = h;;
end