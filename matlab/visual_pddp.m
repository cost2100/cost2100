function visual_pddp(channel, paraEx, paraSt)
%VISUAL_PDDP Visualization of power delay direction profile 
%Default call: visual_pddp(channel, paraEx, paraSt)

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

tau_max_stack = [];
caxis_stack = [];

for m = 1:length(channel)
    if numel(channel{m}.h)~=0
        tau_max_stack = [tau_max_stack max(round(channel{m}.h(:,5)/paraEx.sample_rate))];
    end
end

tau_max = max(tau_max_stack);
tau_min = 0;
% We assume isotropic tranceivers
Ptau = zeros(length(channel),tau_max+21); % Power-delay
Paod = zeros(length(channel),361);
Peod = zeros(length(channel),181);
Paoa = zeros(length(channel),361);
Peoa = zeros(length(channel),181);
for m = 1:length(channel)    
    h_tmp = channel{m}.h(:,:);  
    if numel(h_tmp)~=0
        ind_tau = round(h_tmp(:,5)/paraEx.sample_rate)+1;
        ind_aod = round(h_tmp(:,1)*180/pi)+181;
        ind_eod = round(h_tmp(:,2)*180/pi)+91;
        ind_aoa = round(h_tmp(:,3)*180/pi)+181;
        ind_eoa = round(h_tmp(:,4)*180/pi)+91;        
        
        amp_h = h_tmp(:,6);
            
        if length(ind_aod)>0
            for n = 1:length(ind_tau)
                if ind_tau(n)>0
                    Ptau(m,ind_tau(n)) = Ptau(m,ind_tau(n)) + amp_h(n);
                end
            end
        end
        if length(ind_aod)>0
            for n = 1:length(ind_aod)
                Paod(m,ind_aod(n)) = Paod(m,ind_aod(n)) + amp_h(n);
            end
        end
        if length(ind_eod)>0
            for n = 1:length(ind_eod)
                Peod(m,ind_eod(n)) = Peod(m,ind_eod(n)) + amp_h(n);
            end
        end
        if length(ind_aoa)>0
            for n = 1:length(ind_aoa)
                Paoa(m,ind_aoa(n)) = Paoa(m,ind_aoa(n)) + amp_h(n);
            end
        end
        if length(ind_eoa)>0
            for n = 1:length(ind_eoa)
                Peoa(m,ind_eoa(n)) = Peoa(m,ind_eoa(n)) + amp_h(n);
            end
        end
    end
    % LOS
    tau_0 = round(channel{m}.h_los(5)/paraEx.sample_rate)+1;
    aod_0 = round(channel{m}.h_los(1)*180/pi)+181;
    eod_0 = round(channel{m}.h_los(2)*180/pi)+91;
    aoa_0 = round(channel{m}.h_los(3)*180/pi)+181;
    eoa_0 = round(channel{m}.h_los(4)*180/pi)+91;
    if tau_0>0
    Ptau(m,tau_0) = Ptau(m,tau_0)+channel{m}.h_los(6);
    end
    Paod(m,aod_0) = Paod(m,aod_0)+channel{m}.h_los(6);
    Peod(m,eod_0) = Peod(m,eod_0)+channel{m}.h_los(6);
    Paoa(m,aoa_0) = Paoa(m,aoa_0)+channel{m}.h_los(6);
    Peoa(m,eoa_0) = Peoa(m,eoa_0)+channel{m}.h_los(6);
end    

fig1 = figure;
surf([0:tau_max+20]*paraEx.sample_rate,[0:length(channel)-1]*paraEx.snap_rate, 20*log10(abs(Ptau)+eps) )
axis([paraEx.sample_rate (tau_max+20)*paraEx.sample_rate 0 (length(channel)-1)*paraEx.snap_rate])
shading interp
view([90 90])
colorbar 
caxis_stack = [caxis_stack caxis];
xlabel('Delay[S]')
ylabel('Time line[S]')
setFontsize
%--------------------------------------------------
fig2 = figure;
surf([0:360],[0:length(channel)-1]*paraEx.snap_rate,20*log10(abs(Paod)+eps))
axis([0 360 0 (length(channel)-1)*paraEx.snap_rate])
shading interp
view([90 90])
colorbar 
caxis_stack = [caxis_stack caxis];
xlabel('Azimuth of Departure[\circ]')
ylabel('Time line[S]')
setFontsize
%--------------------------------------------------
fig3 = figure;
surf([0:180],[0:length(channel)-1]*paraEx.snap_rate,20*log10(abs(Peod)+eps))
axis([1 180 0 (length(channel)-1)*paraEx.snap_rate])
shading interp
view([90 90])
colorbar 
caxis_stack = [caxis_stack caxis];
xlabel('Elevation of Departure[\circ]')
ylabel('Time line[S]')
setFontsize
%--------------------------------------------------
fig4 = figure;
surf([0:360],[0:length(channel)-1]*paraEx.snap_rate,20*log10(abs(Paoa)+eps))
axis([0 360 0 (length(channel)-1)*paraEx.snap_rate])
shading interp
view([90 90])
colorbar 
caxis_stack = [caxis_stack caxis];
xlabel('Azimuth of Arrival[\circ]')
ylabel('Time line[S]')
setFontsize
%--------------------------------------------------
fig5 = figure;
surf([0:180],[0:length(channel)-1]*paraEx.snap_rate,20*log10(abs(Peoa)+eps))
axis([1 180 0 (length(channel)-1)*paraEx.snap_rate])
shading interp
view([90 90])
colorbar 
caxis_stack = [caxis_stack caxis];
xlabel('Elevation of Arrival[\circ]')
ylabel('Time line[S]')
setFontsize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unify the color axis
val = (max(caxis_stack)-min(caxis_stack))*0.6;
figure(fig1)
caxis([min(caxis_stack)+val max(caxis_stack)])
figure(fig2)
caxis([min(caxis_stack)+val max(caxis_stack)])
figure(fig3)
caxis([min(caxis_stack)+val max(caxis_stack)])
figure(fig4)
caxis([min(caxis_stack)+val max(caxis_stack)])
figure(fig5)
caxis([min(caxis_stack)+val max(caxis_stack)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
