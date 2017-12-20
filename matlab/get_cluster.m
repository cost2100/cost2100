function cluster = get_cluster(VR, VRtable, paraEx, paraSt)
%GET_CLUSTER Function to generate the cluster
%
%Default call: 
%cluster = get_cluster(VR, VRtable, paraEx, paraSt)
%------
%Input:
%------
%paraEx: External parameters
%paraSt: Stochastic parameters
%VRtable: VR assignment table
%VR: VR distribution
%------
%Output:
%------
%cluster: Cluster information
% .idx: Label of cluster
% .refBS: Reference BS for cluster distribution
% .refVR: Reference VR for cluster distribution
% .type: Cluster type: 0 local, 1 single, 2 twin
% .pos_c_BS: Position of cluster at BS side
% .pos_c_MS: Position of cluster at MS side
% .a_c_BS: Delay spatial spread at BS side
% .b_c_BS: AoD spatial spread at BS side
% .h_c_BS: EoD spatial spread at BS side
% .a_c_MS: Delay spatial spread at MS side
% .b_c_MS: AoA spatial spread at MS side
% .h_c_MS: EoA spatial spread at MS side
% .shadow_f: Shadowing fading
% .tau_c_link: Cluster link delay
% .Phi_c_BS: Azimuth angle viewed from BS side, for computation
% .Theta_c_BS: Elevation angle viewed from BS side, for computation
% .Phi_c_MS: Azimuth angle viewed from MS side, for computation
% .Theta_c_MS: Elevation angle viewed from MS side, for computation

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

nCluster = max(VRtable(1, :, 2)); % Number of clusters
nBS = paraEx.num_BS; % Number of BSs
posBS = paraEx.pos_BS; % Positions at BS

% Get the random delay angular spread for each cluster
corr_randn = randn(nCluster, 6)*paraSt.corr_mat; % Generate cross-correlation matrix for shadow fading, cluster spreads, etc
shadow_f = 10.^(0.1*corr_randn(:,1)*paraSt.sigma_sf); % Shadow fading
tau_c = paraSt.mu_tau.*10.^(0.1*paraSt.sigma_tau*corr_randn(:, 2)); % Delay spread
theta_c_BS = paraSt.mu_theta_BS*10.^(0.1*paraSt.sigma_theta_BS*corr_randn(:, 3)); % Elevation spread BS
phi_c_BS = paraSt.mu_phi_BS*10.^(0.1*paraSt.sigma_phi_BS*corr_randn(:, 4)); % Azimtuh spread BS
theta_c_MS = paraSt.mu_theta_MS*10.^(0.1*paraSt.sigma_theta_MS*corr_randn(:, 5)); % Elevation spread MS
phi_c_MS = paraSt.mu_phi_MS*10.^(0.1*paraSt.sigma_phi_MS*corr_randn(:, 6)); % Azimuth spread MS
d_tau = tau_c*paraEx.c0/2; % Spatial delay spread

for m = 1:nCluster    
    cluster(m).idx = m; % Label the cluster

    VRGrp = find(VRtable(1,:,2)==m); % Find the VR group
    
    % Find the reference BS
    BSWeight = sum(VRtable(:,VRGrp,1),2);
    BSIdx = find( BSWeight==max(BSWeight) );    
    numBSIdx = length(BSIdx);    
    cluster(m).refBS = BSIdx(ceil(rand(1)*numBSIdx+eps/1e80));

    % Find the reference VR
    VRBS = VRGrp(find(VRtable(cluster(m).refBS,VRGrp,1)==1));
    BSCom = sum(VRtable(:,VRBS,1));
    if numel(VRBS)==0 keyboard; end
    if any(BSCom>1)
        VRCom = VRBS(BSCom>1);
        numVRCom = length(VRCom);
        cluster(m).refVR= VRCom(ceil(rand(1)*numVRCom+eps/1e80));
    else
        numVRGrp = length(VRBS);
        cluster(m).refVR= VRBS(ceil(rand(1)*numVRGrp+eps/1e80));
    end   
    
    % Determine the cluster type: 1 single, 2 twin cluster
    cluster(m).type = (rand>paraSt.k_sel)+1;
end

% Parameterize the clusters
for m = 1:nCluster
    switch cluster(m).type
        case 1 % Single cluster
            % BS side
            pos_VR = VR(cluster(m).refVR, :); % VR position [x y]
            pos_BS = posBS(cluster(m).refBS, :); % BS position [x y z]
            VR_direct = pos_VR-pos_BS(1:2); % VR to BS vector
            VR_direct(3) = 0;
            
            [phi, theta, r] = cart2sph(VR_direct(1), VR_direct(2), VR_direct(3));
            phi_new = phi+deg2rad(get_random(paraSt.pdf_phi_c, paraSt.phi_c)); % Determine the cluster to BS-VR angle
            theta_new = theta+deg2rad(get_random(paraSt.pdf_theta_c, paraSt.theta_c)); % Determine the cluster to BS-VR angle
            r_new = get_random(paraSt.pdf_r_c, paraSt.para_r_c); % Determine the cluster to BS distance
            [c_direct(1), c_direct(2), c_direct(3)] = sph2cart(phi_new, theta_new, r_new); % Cluster to VR vector
            
            cluster(m).pos_c_BS = pos_BS+c_direct; % Determine cluster at BS side 
            cluster(m).pos_c_MS = cluster(m).pos_c_BS; % Determine cluster at MS side
            
            d_c_BS = calc_dist(cluster(m).pos_c_BS, pos_BS);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Cluster spread, shadowing
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cluster(m).a_c_BS = d_tau(m); % Spatial delay spread at BS side (radius) [m]
            cluster(m).b_c_BS = d_c_BS*tan(phi_c_BS(m)/180*pi/2); % Width of cluster (spatial AoD spread at BS side, radius) [m]
            cluster(m).h_c_BS = d_c_BS*tan(theta_c_BS(m)/180*pi/2);% Height of cluster (spatial EoD spread at BS side, radius) [m]
            cluster(m).a_c_MS = cluster(m).a_c_BS; % Spatial delay spread at MS side [m]
            cluster(m).b_c_MS = cluster(m).b_c_BS; % Spatial AoA spread at MS side
            cluster(m).h_c_MS = cluster(m).h_c_BS; % Spatial EoA spread at MS side
            cluster(m).shadow_f = shadow_f(m); % Cluster shadow fading
            cluster(m).tau_c_link = 0; % No link delay
            
        case 2 % Twin cluster
            % BS side
            pos_VR = VR(cluster(m).refVR, :); % VR position [x y]
            pos_BS = posBS(cluster(m).refBS, :); % BS position [x y z]
            VR_direct = pos_VR-pos_BS(1:2); % VR to BS vector
            VR_direct(3) = 0;
            
            [phi, theta, r] = cart2sph(VR_direct(1), VR_direct(2), VR_direct(3));
            phi_new = phi+deg2rad(get_random(paraSt.pdf_phi_c, paraSt.phi_c)); % Determine the cluster to BS-VR angle
            theta_new = theta+deg2rad(get_random(paraSt.pdf_theta_c, paraSt.theta_c)); % Determine the cluster to BS-VR angle
            r_new = get_random(paraSt.pdf_r_c,paraSt.para_r_c); % Determine the cluster to BS distance            
            [c_direct(1), c_direct(2), c_direct(3)] = sph2cart(phi_new, theta_new, r_new); % Cluster to VR vector
            
            cluster(m).pos_c_BS = pos_BS+c_direct; % Determine cluster at BS side
            
            d_c_BS = calc_dist(cluster(m).pos_c_BS, pos_BS); % Distance between cluster and BS
            
            % MS side
            d_c_MS = d_c_BS*tan(phi_c_BS(m)/180*pi/2)/tan(phi_c_MS(m)/180*pi/2);
            phi_new = phi+deg2rad(get_random(paraSt.pdf_phi_c, paraSt.phi_c));
            theta_new = theta+deg2rad(get_random(paraSt.pdf_theta_c, paraSt.theta_c));
            [c_direct(1), c_direct(2), c_direct(3)] = sph2cart(phi_new, theta_new, d_c_MS);
            pos_c_MS = c_direct+[pos_VR 0]; % Cluster position from the VR
            cluster(m).pos_c_MS = pos_c_MS; % Cluster position at MS side
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Cluster spread, shadowing, link delay
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cluster(m).a_c_BS = d_tau(m); % Delay spread [m]
            cluster(m).b_c_BS = d_c_BS*tan(phi_c_BS(m)/180*pi/2); % Azimuth spread seen from BS side
            cluster(m).h_c_BS = d_c_BS*tan(theta_c_BS(m)/180*pi/2); % Elevation spread seen from BS side
            
            cluster(m).a_c_MS = cluster(m).a_c_BS; % Delay spread seen from the MS side
            cluster(m).b_c_MS = d_c_MS*tan(phi_c_MS(m)/180*pi/2); % Azimuth spread seen from MS side
            cluster(m).h_c_MS = d_c_MS*tan(theta_c_MS(m)/180*pi/2); % Elevation spread seen from MS side
            
            cluster(m).shadow_f = shadow_f(m);
%             % generate a randn extra delay factor, modified by meifang
%             para = [paraSt.mu_tauCLink paraSt.sigma_tauCLink];
%             var =  para(1)+randn(1)*para(2);
%             while var < 0
%                 var =  para(1)+randn(1)*para(2);
%             end
%             LOS_delay = r/3e8; % the distance between BS and VR
%             cluster(m).tau_c_link = LOS_delay + var-(d_c_MS+d_c_BS)/3e38;
            
            %cluster(m).tau_c_link = ...%calc_dist(cluster(m).pos_c_BS,cluster(m).pos_c_MS)/paraEx.c0+...
            %paraSt.mu_tauCLink*10.^(0.1*paraSt.sigma_tauCLink*randn(1)); %Cluster link delay            
            % change to exponential distribution by Meifang 2012-02-03
            % first generate a uniform U(0,1) then go to expoinentail distribution            
            while (1)
                tempvar = rand(1);
                tau_c_link = - paraSt.mu_tauCLink * log(1-tempvar) ;
                if tau_c_link >= paraSt.min_tauCLink
                    cluster(m).tau_c_link = tau_c_link;
                    break;
                end
            end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For computation purpose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = 1:nCluster
    pos_BS = posBS(cluster(m).refBS,:); % BS position [x y z]
    pos_VR = VR(cluster(m).refVR,:); % VR position [x y]    
    [Phi_c_BS, Theta_c_BS, tmp] = cart2sph(pos_BS(1)-cluster(m).pos_c_BS(1),pos_BS(2)-cluster(m).pos_c_BS(2),pos_BS(3)-cluster(m).pos_c_BS(3));
    [Phi_c_MS, Theta_c_MS, tmp] = cart2sph(pos_VR(1)-cluster(m).pos_c_MS(1),pos_VR(2)-cluster(m).pos_c_MS(2),0-cluster(m).pos_c_MS(3));
    cluster(m).Phi_c_BS = Phi_c_BS; % Azimuth seen from BS side
    cluster(m).Theta_c_BS = Theta_c_BS; % Elevation seen from BS side
    cluster(m).Phi_c_MS = Phi_c_MS; % Azimuth spread seen MS side
    cluster(m).Theta_c_MS = Theta_c_MS; % Elevation seen from MS side
end

end 

function var = get_random(pdf,para)
    switch pdf
        case 'norm'
            var = para(1)+randn(1)*para(2);
        case 'unif'
            var = rand(1)*(para(2)-para(1))+para(1);
        case 'mycase'
            var = 200;
    end
end                                 