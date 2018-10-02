function visualize_channel_env(paraEx, paraSt, link, env)
%VISUALIZE_CHANNEL_ENV Visualization function to display the generated
%channel environment
%Default calling: visual_channel(paraEx, paraSt, link, env)

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

figure, 
% Visualize simulation area (cell) at MS side
draw_circ([0, 0], paraEx.net_radii, 'k', 2);
axis equal
hold on 

% Visualize cluster visibility regions
for idx_MS_VR = 1:size(env.MS_VR, 1)
    draw_circ(env.MS_VR(idx_MS_VR, :), paraSt.r_c, '-.b', 1);
end
% BS positions
for idx_BS = 1:paraEx.num_BS
    plot(paraEx.pos_BS(idx_BS, 1), paraEx.pos_BS(idx_BS, 2), 'rx', 'MarkerSize', 10, 'LineWidth', 3);
end

% MS starting positions
for idx_MS = 1:paraEx.num_MS
    plot(paraEx.pos_MS(idx_MS, 1), paraEx.pos_MS(idx_MS, 2), 'kx', 'MarkerSize', 10, 'LineWidth', 3);
end

% MS ending positions due to linear movement
for idx_MS = 1:paraEx.num_MS
    pos_MS_end = paraEx.pos_MS(idx_MS, :)+paraEx.velo_MS(idx_MS, :)/paraEx.snap_rate*(paraEx.snap_num-1);
    plot(pos_MS_end(1), pos_MS_end(2), 'gx', 'MarkerSize', 10, 'LineWidth', 3);
end

% Active VR for 1st BS, 1st MS and 1st snapshot
active_cluster_idx = link(1, 1).channel{1, 1}.active_c; % Active cluster
for idx_active_c = 1:length(active_cluster_idx)
    active_MS_VR_idx = find(env.VRtable(1, :, 2)==active_cluster_idx(idx_active_c));
    draw_circ(env.MS_VR(active_MS_VR_idx, :), paraSt.r_c, '-.r', 1);
end
grid on
hold off
title('Simulation area at MS side, cluster visibility regions (active ones for 1st user), BS and MS positions');
xlabel('x [m]'), ylabel('y [m]');

figure,
% Visualize simulation area at BS side (the 1st BS)
draw_line([paraEx.pos_BS(1, 1)-paraEx.BS_range/2, paraEx.pos_BS(1, 2)], [paraEx.pos_BS(1, 1)+paraEx.BS_range/2, paraEx.pos_BS(1, 2)], 'k');
hold on
for idx_BS_pos = 1:paraEx.BS_spacing_num(1)
    pos_vec = (-1*paraEx.BS_spacing(1, 1)/2-paraEx.BS_spacing(1, 1)*(paraEx.BS_spacing_num(1)/2-1)):paraEx.BS_spacing(1, 1):(paraEx.BS_spacing(1, 1)/2+paraEx.BS_spacing(1, 1)*(paraEx.BS_spacing_num(1)/2-1));
    if ~isempty(pos_vec)
        BS_simu_pos(:, 1) = paraEx.pos_BS(1, 1)+pos_vec; % [x]
    else
        BS_simu_pos(:, 1) = paraEx.pos_BS(1, 1);
    end
    BS_simu_pos(:, 2) = paraEx.pos_BS(1, 2); % [y]
    BS_simu_pos(:, 3) = paraEx.pos_BS(1, 3); % [z]
    plot(BS_simu_pos(:, 1), BS_simu_pos(:, 2), 'rx', 'MarkerSize', 10, 'LineWidth', 3);
end

% Active cluster VR at BS side
for idx_active_c = 1:length(active_cluster_idx)
    actvie_BS_VR_idx = active_cluster_idx(idx_active_c);
    if (env.BS_VR_len(actvie_BS_VR_idx) < Inf)
        draw_circ(env.BS_VR(actvie_BS_VR_idx, :), env.BS_VR_len(actvie_BS_VR_idx)/2, '-.r', 1);
    end
end

hold off
grid on
axis equal
if (paraEx.BS_range>0)
    xlim([paraEx.pos_BS(1, 1)-paraEx.BS_range*10, paraEx.pos_BS(1, 1)+paraEx.BS_range*10]);
end
title('Simulation area at BS side, cluster visibility regions (active ones for 1st user)');
xlabel('x [m]'), ylabel('y [m]');

% Visualize clusters
figure,
for idx_BS = 1:paraEx.num_BS
    plot3(paraEx.pos_BS(idx_BS, 1), paraEx.pos_BS(idx_BS, 2), paraEx.pos_BS(idx_BS, 3), 'rx', 'MarkerSize', 10, 'LineWidth', 3);
end
hold on

% MS starting positions
for idx_MS = 1:paraEx.num_MS
    plot3(paraEx.pos_MS(idx_MS, 1), paraEx.pos_MS(idx_MS, 2), paraEx.pos_MS(idx_MS, 3), 'kx', 'MarkerSize', 10, 'LineWidth', 3);
end

% MS ending positions due to linear movement
for idx_MS = 1:paraEx.num_MS
    pos_MS_end = paraEx.pos_MS(idx_MS, :)+paraEx.velo_MS(idx_MS, :)/paraEx.snap_rate*(paraEx.snap_num-1);
    plot3(pos_MS_end(1), pos_MS_end(2), pos_MS_end(3), 'gx', 'MarkerSize', 10, 'LineWidth', 3);
end

% Active clusters (1st BS, 1 user, 1st snapshot)
for idx_active_c = 1:length(active_cluster_idx)
    switch env.cluster(active_cluster_idx(idx_active_c)).type
        case 1 % Single cluster
            pos_c_BS = env.cluster(active_cluster_idx(idx_active_c)).pos_c_BS;
            radi_c_BS = [env.cluster(active_cluster_idx(idx_active_c)).a_c_BS,...
                         env.cluster(active_cluster_idx(idx_active_c)).b_c_BS,...
                         env.cluster(active_cluster_idx(idx_active_c)).h_c_BS];
            rot_c_BS = [env.cluster(active_cluster_idx(idx_active_c)).Phi_c_BS,...
                        env.cluster(active_cluster_idx(idx_active_c)).Theta_c_BS];
            draw_ellipsoid(pos_c_BS, radi_c_BS, rot_c_BS, 1, 0.3, [], []);
            
            % MPC
            mpc_peak_pos(:, 1) = env.mpc(active_cluster_idx(idx_active_c)).gain_center(:, 1)+...
                                 squeeze(env.MS_VR(env.VRtable(1, :, 2)==active_cluster_idx(idx_active_c), 1));
            mpc_peak_pos(:, 2) = env.mpc(active_cluster_idx(idx_active_c)).gain_center(:, 2)+...
                                 squeeze(env.MS_VR(env.VRtable(1, :, 2)==active_cluster_idx(idx_active_c), 2));
            mpc_gain_sigma = env.mpc(active_cluster_idx(idx_active_c)).gain_radius;
            a_mpc_gain = exp(-1*((mpc_peak_pos(:, 1)-paraEx.pos_MS(1, 1)).^2./(2*mpc_gain_sigma.^2)+(mpc_peak_pos(:, 2)-paraEx.pos_MS(1,2)).^2./(2*mpc_gain_sigma.^2)));
            [sort_val, mpc_idx_sort] = sort(a_mpc_gain, 'ascend');
            scatter3(env.mpc(active_cluster_idx(idx_active_c)).pos_BS(mpc_idx_sort, 1),...
                     env.mpc(active_cluster_idx(idx_active_c)).pos_BS(mpc_idx_sort, 2),...
                     env.mpc(active_cluster_idx(idx_active_c)).pos_BS(mpc_idx_sort, 3),...
                     10, pow2db(abs(a_mpc_gain(mpc_idx_sort)).^2), 'filled');
            colormap(jet), colorbar, caxis([-60, 0]);
            
            % MPC-VR 3-dB 
            for mpc_vr_idx = find(a_mpc_gain.^2 > .5).'
                draw_circ(mpc_peak_pos(mpc_vr_idx,1:2), mpc_gain_sigma(mpc_vr_idx)*sqrt(log(2)), 'b--', 1);
                plot(mpc_peak_pos(mpc_vr_idx,1), mpc_peak_pos(mpc_vr_idx,2), 'k.', 'MarkerSize', 10);
                draw_line(mpc_peak_pos(mpc_vr_idx,1:2), env.mpc(active_cluster_idx(idx_active_c)).pos_BS(mpc_vr_idx, 1:2), 'k--');
            end
        case 2 % Twin cluster
            % From BS side
            pos_c_BS = env.cluster(active_cluster_idx(idx_active_c)).pos_c_BS;
            radi_c_BS = [env.cluster(active_cluster_idx(idx_active_c)).a_c_BS,...
                         env.cluster(active_cluster_idx(idx_active_c)).b_c_BS,...
                         env.cluster(active_cluster_idx(idx_active_c)).h_c_BS];
            rot_c_BS = [env.cluster(active_cluster_idx(idx_active_c)).Phi_c_BS,...
                        env.cluster(active_cluster_idx(idx_active_c)).Theta_c_BS];
            draw_ellipsoid(pos_c_BS, radi_c_BS, rot_c_BS, 1, 0.3, [], []);
            
            % From MS side
            pos_c_MS = env.cluster(active_cluster_idx(idx_active_c)).pos_c_MS;
            radi_c_MS = [env.cluster(active_cluster_idx(idx_active_c)).a_c_MS,...
                         env.cluster(active_cluster_idx(idx_active_c)).b_c_MS,...
                         env.cluster(active_cluster_idx(idx_active_c)).h_c_MS];
            rot_c_MS = [env.cluster(active_cluster_idx(idx_active_c)).Phi_c_MS,...
                        env.cluster(active_cluster_idx(idx_active_c)).Theta_c_MS];
            draw_ellipsoid(pos_c_MS, radi_c_MS, rot_c_MS, 1, 0.3, [], []);
            
            % Link the twin clusters
            draw_line(pos_c_BS, pos_c_MS, 'k');
            
            % MPC
            mpc_peak_pos(:, 1) = env.mpc(active_cluster_idx(idx_active_c)).gain_center(:, 1)+...
                                 squeeze(env.MS_VR(env.VRtable(1, :, 2)==active_cluster_idx(idx_active_c), 1));
            mpc_peak_pos(:, 2) = env.mpc(active_cluster_idx(idx_active_c)).gain_center(:, 2)+...
                                 squeeze(env.MS_VR(env.VRtable(1, :, 2)==active_cluster_idx(idx_active_c), 2));
            mpc_gain_sigma = env.mpc(active_cluster_idx(idx_active_c)).gain_radius;
            a_mpc_gain = exp(-1*((mpc_peak_pos(:, 1)-paraEx.pos_MS(1, 1)).^2./(2*mpc_gain_sigma.^2)+(mpc_peak_pos(:, 2)-paraEx.pos_MS(1,2)).^2./(2*mpc_gain_sigma.^2)));
            [sort_val, mpc_idx_sort] = sort(a_mpc_gain, 'ascend');
            scatter3(env.mpc(active_cluster_idx(idx_active_c)).pos_MS(mpc_idx_sort, 1),...
                     env.mpc(active_cluster_idx(idx_active_c)).pos_MS(mpc_idx_sort, 2),...
                     env.mpc(active_cluster_idx(idx_active_c)).pos_MS(mpc_idx_sort, 3),...
                     10, pow2db(abs(a_mpc_gain(mpc_idx_sort)).^2), 'filled');
            colormap(jet), colorbar, caxis([-60, 0]);
            scatter3(env.mpc(active_cluster_idx(idx_active_c)).pos_BS(mpc_idx_sort, 1),...
                     env.mpc(active_cluster_idx(idx_active_c)).pos_BS(mpc_idx_sort, 2),...
                     env.mpc(active_cluster_idx(idx_active_c)).pos_BS(mpc_idx_sort, 3),...
                     10, pow2db(abs(a_mpc_gain(mpc_idx_sort)).^2), 'filled');
            colormap(jet), colorbar, caxis([-60, 0]);
            
            % MPC-VR 3-dB 
            for mpc_vr_idx = find(a_mpc_gain.^2 > .5).'
                draw_circ(mpc_peak_pos(mpc_vr_idx,1:2), mpc_gain_sigma(mpc_vr_idx)*sqrt(log(2)), 'b--', 1);
                plot(mpc_peak_pos(mpc_vr_idx,1), mpc_peak_pos(mpc_vr_idx,2), 'k.', 'MarkerSize', 10);
                draw_line(mpc_peak_pos(mpc_vr_idx,1:2), env.mpc(active_cluster_idx(idx_active_c)).pos_MS(mpc_vr_idx, 1:2), 'k--');
                draw_line(env.mpc(active_cluster_idx(idx_active_c)).pos_BS(mpc_vr_idx, 1:2), paraEx.pos_BS(1, 1:2), 'k--');
            end
    end
end
% Plot MPCs

hold off
grid on
title('Simulated environment: active clusters for 1st user');
xlabel('x [m]'), ylabel('y [m]'), zlabel('z [m]');