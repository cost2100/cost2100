function IR = get_IR(channel, paraEx, paraSt)
%GET_IR Function to compute the impulse response of the SISO channel
%Default call IR = get_IR(channel, paraEx, paraSt)
%-------
%Input:
%-------
%
%channel: channel 
%paraEx, paraSt: external and stochastic parameters
%
%------
%Output:
%------
%IR: unfiltered Impulse Response

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (C)2008 LIU Ling-Feng, ICTEAM, UCL, Belgium 
%This file is part of cost2100.
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

sample_rate = paraEx.sample_rate; %Sampling rate
delay_max = paraEx.delay_max; %Maxium delay is 5 times of the net radius
irMax = ceil(delay_max/sample_rate); %The maximum IR delay length
IR = zeros(irMax,1);

for m = 1:length(channel.h(:,1))
        mInd = mod(ceil(channel.h(m,5)/sample_rate),irMax);
        if mInd==0 mInd = length(irMax); end
        IR(mInd) = IR(mInd)+channel.h(m,6);
%         IR_v(mInd) = IR(mInd)+channel.h(n_c,n_mpc,7);        
%         IR_h(mInd) = IR(mInd)+channel.h(n_c,n_mpc,8);        
end