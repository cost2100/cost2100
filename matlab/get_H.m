function H = get_H(channel, txFull, rxFull, txRot, rxRot, paraEx)
%GET_H Get the MIMO channel impulse reponse matrix
%
%Default call:
%H = get_H(channel, txFull, rxFull, txRot, rxRot, paraEx)
%-------
%Input:
%-------
%channel: Channel information
%txFull: Tx antenna array radiation pattern
% .antennaResponse[numAntenna,azimuth,elevation]
% .azimuthRange
% .elevationRange
% .dAngle: Angle resolution
%rxFull: Rx antenna array radiation pattern
%see above structure of txFull
%txRot: Tx antenna array rotation [azi ele]
%rxRot: Rx antenna array rotation [azi ele]
%paraEx: External parameters
%------
%Output:
%------
%H: MIMO channel impulse response matrix

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

if paraEx.freq_start==paraEx.freq_stop
    type = 'Narrowband'; 
else
    type = 'Wideband'; 
end

sample_rate = paraEx.sample_rate; % IR sampling rate
tauIndMax = ceil(paraEx.delay_max/sample_rate); % Maximum delay index
M=length(txFull.antennaResponse(:,1,1));
N=length(rxFull.antennaResponse(:,1,1));
txAziRange=txFull.azimuthRange;
txEleRange=txFull.elevationRange;
rxAziRange=rxFull.azimuthRange;
rxEleRange=rxFull.elevationRange;

switch type
    case 'Narrowband'
        H=zeros(M,N);
    case 'Wideband'
        H=zeros(tauIndMax,M,N);
end

for m = 1:length(channel.h(:,1))
    tauInd = ceil(channel.h(m,5)/sample_rate);
    if tauInd > tauIndMax % If out of the range of delay
        continue;
    elseif tauInd <= 0
        continue;
    else
        aod=channel.h(m,1);
        eod=channel.h(m,2);
        bsv=[sin(eod)*cos(aod) sin(eod)*sin(aod) cos(eod)];
        bsv=bsv*rotate_matrix(txRot(1),txRot(2));
        aod=acos(bsv(3));
        eod=atan2(bsv(2),bsv(1));
        
        aoa=channel.h(m,3);
        eoa=channel.h(m,4);
        msv=[sin(eoa)*cos(aoa) sin(eoa)*sin(aoa) cos(eoa)];
        msv=msv*rotate_matrix(rxRot(1),rxRot(2));
        aoa=acos(msv(3));
        eoa=atan2(msv(2),msv(1));
        
        [tmp txAziInd]= min(abs(txAziRange-aod*180/pi));
        [tmp txEleInd]= min(abs(txEleRange-eod*180/pi));
        [tmp rxAziInd]= min(abs(rxAziRange-aoa*180/pi));
        [tmp rxEleInd]= min(abs(rxEleRange-eoa*180/pi));
        
        Tx=squeeze(txFull.antennaResponse(:,txAziInd,txEleInd));
        Rx=squeeze(rxFull.antennaResponse(:,rxAziInd,rxEleInd));
        
        switch type
            case 'Narrowband' 
                H = H+channel.h(m,6)*sqrt(Tx)*sqrt(conj(Rx'));
            case 'Wideband' 
                H(tauInd,:,:) = squeeze(H(tauInd,:,:))+channel.h(m,6)*sqrt(Tx)*sqrt(conj(Rx'));
        end
    end
end
        
% LOS
tauInd=ceil(channel.h_los(5)/sample_rate);
if tauInd > tauIndMax % If out of the range of delay    
elseif tauInd <= 0
else
    aod=channel.h_los(1);
    eod=channel.h_los(2);
    bsv=[sin(eod)*cos(aod) sin(eod)*sin(aod) cos(eod)];
    bsv=bsv*rotate_matrix(txRot(1),txRot(2));
    aod=acos(bsv(3));
    eod=atan2(bsv(2),bsv(1));

    aoa=channel.h_los(3);
    eoa=channel.h_los(4);
    msv=[sin(eoa)*cos(aoa) sin(eoa)*sin(aoa) cos(eoa)];
    msv=msv*rotate_matrix(rxRot(1),rxRot(2));
    aoa=acos(msv(3));
    eoa=atan2(msv(2),msv(1));
        
    [tmp txAziInd]= min(abs(txAziRange-aod*180/pi));
    [tmp txEleInd]= min(abs(txEleRange-eod*180/pi));
    [tmp rxAziInd]= min(abs(rxAziRange-aoa*180/pi));
    [tmp rxEleInd]= min(abs(rxEleRange-eoa*180/pi));
    Tx=squeeze(txFull.antennaResponse(:,txAziInd,txEleInd));
    Rx=squeeze(rxFull.antennaResponse(:,rxAziInd,rxEleInd));
    switch type
        case 'Narrowband'
            H = H+channel.h_los(6)*sqrt(Tx)*sqrt(conj(Rx'));
        case 'Wideband'
            H(tauInd,:,:) = squeeze(H(tauInd,:,:))+channel.h_los(6)*sqrt(Tx)*sqrt(conj(Rx'));
    end
end

% Pulse shaping of H
% [num den]=rcosine(1,paraEx.overSample,'fir/normal'); %Raised cosine filter para.
% HOld=H;
% clear H
% switch type
%     case 'Narrowband'
%         for m=1:M
%             for n=1:N
%                h=HOld(m,n);
%                H(:,m,n)=rcosflt(h,1,paraEx.overSample,'filter/Fs',num,den);
%             end
%         end
%     case 'Wideband'
%         for m=1:M
%             for n=1:N
%                 h=HOld(:,m,n);
%                 H(:,m,n)=rcosflt(h,1,paraEx.overSample,'filter/Fs',num,den);
%             end
%         end
% end