function VRtable = get_VRtable(paraEx, paraSt)
%GET_VRTABLE Function to generate the VR assignment table, which includes
%the commonness of the VR to the BS and VR group assignment
%
%Default call: 
%VRtable = getVRtable(paraEx, paraSt)
%------
%Input:
%------
%paraEx: External parameters
%paraSt: Stochastic parameters
%------
%Output:
%------
%VRtable: the VR assignment table (numBS, numVR, [BSassign VRgroupAssign])

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

tic % For generating state of rand state
numBS = paraEx.num_BS; % BS number
rho_VR = (paraSt.n_c_far)/(pi*(paraSt.r_c-paraSt.l_c)^2); % VR density
n_VR = round(rho_VR*pi*paraEx.net_radii^2); % Total number of VRs in the cell of one BS
numVRlogi = n_VR*numBS; % Total number of logical VRs for all BSs
numVRtrue = round(numVRlogi/sum([1:numBS].*paraSt.BS_common)); % Total number of physical VRs for all BSs

VRtable = zeros(numBS,numVRlogi,2);

% Assign the k-BS-common VRs
for m = 1:numBS
	num_BSCC(m) = round(numVRtrue*paraSt.BS_common(m)); % Total number of each kBS-CC
end
for m = numBS:-1:2	
	for n = sum(num_BSCC(end:-1:m+1))+1:sum(num_BSCC(end:-1:m))
		
		% Check if BSs VR density is reached
		for nBS = 1:numBS
			BSFull(nBS) = (sum(VRtable(nBS,:,1))>=n_VR);
        end
        
        if sum(~BSFull)<m % If unfilled BSs are less than m to assign a mBS-CC
            continue
        else
            % The index of unfilled BS
            BSEmp = find(~BSFull);
            numBSEmp = length(BSEmp);
            
            S = toc*1e5; rand('twister',S); % Initialize a rand state            
            idx = BSEmp(ceil(rand(m,1)*numBSEmp+eps/1e80)); % Generate un-repeated random integers
            idxChk = diff(sort(idx,'ascend')); % idx check
            while any(idxChk==0) % If idx contain duplicated numbers
                S = toc*1e5; rand('twister',S); % Initialize a rand state                
                idx = BSEmp(ceil(rand(m,1)*numBSEmp+eps/1e80)); % Generate un-repeated random integers
                idxChk = diff(sort(idx,'ascend')); % idx check
            end            
            VRtable(idx,n,1)=1; % Assign the VR to BS(idx)
        end
	end	
end

% Assign the non-common VRs
VRCount = length(find(sum(squeeze(VRtable(:,:,1)),1)>0));
for m = 1:numBS
	if sum(VRtable(m,:,1))<n_VR
		ncVR = n_VR-sum(VRtable(m,:,1)); % Remaining number of non-common VR
		VRtable(m,VRCount+1:VRCount+ncVR,1) = 1;
		VRCount = VRCount+ncVR;		
	end
end

% Truncate the VRtable
VRtable = VRtable(:,1:VRCount,:);

% Assign the VR groups
clusterCount = 0;
for m = 1:length(VRtable(1,:,1))
	if VRtable(1,m,2)==0 % If the VR doesn't have a cluster
		clusterCount = clusterCount+1; 
		VRtable(:,m,2) = clusterCount; % Assign a new cluster	
        if isempty(find(VRtable(1, :, 2)==0))
            continue; % After attributing the last VR
        end
        
		num_VRgroup = poissrnd_own(paraSt.MS_common); % Number of extra VRs in the VR group
		if num_VRgroup==0
			continue % If single VR assigned to the cluster
		else
			VREmp = find(VRtable(1, :, 2)==0); % VRs that don't have a cluster
			numVREmp = length(VREmp);
                        
            while num_VRgroup > numVREmp || num_VRgroup ==0 % If there are not enough VRs to group
                num_VRgroup = poissrnd_own(paraSt.MS_common);
            end
            
            S = toc*1e5; rand('twister',S); % Initialize the rand state            
			idx = VREmp(ceil(rand(num_VRgroup,1)*numVREmp+eps/1e80)); % Generate un-repeated random integers
			idxChk = diff(sort(idx,'ascend')); % Check idx		            
			if num_VRgroup >1			
				while any(idxChk==0)	
                    S = toc*1e5; rand('twister',S); % Initialize the rand state                    
					idx = VREmp(ceil(rand(num_VRgroup,1)*numVREmp+eps/1e80));
					idxChk = diff(sort(idx,'ascend'));                    
				end
            end
            
			VRtable(:,idx,2) = clusterCount; % Group the VR(idx)
		end
	else
		continue
	end
end