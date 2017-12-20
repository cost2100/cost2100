function visual_channel(paraEx, paraSt, link, env)
%VISUAL_CHANNEL Visualization function to display the sampled channel
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

channel = link.channel;
BS = link.BS;
MS = link.MS;
cluster = env.cluster;
mpc = env.mpc;

cluster_a_color = 1e3; % Color value
alpha_val  = 0.4; % Transparency
figure_h = gcf;

for ind1 = 1:length(channel) % For each snapshot
    clf     
    
    plot_cluster = channel{ind1}.active_c;
    
    for ind2 = 1:length(plot_cluster)                
        figure(figure_h)        
        hold on
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Draw BS and MS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        plot3(BS(ind1).pos(1),BS(ind1).pos(2),BS(ind1).pos(3),'*','markerSize',14,'lineWidth',3) % Plot the BS position
        plot3(MS(ind1).pos(1),MS(ind1).pos(2),MS(ind1).pos(3),'o','markerSize',14,'lineWidth',3) % Plot the MS position
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Draw single clusters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if cluster(plot_cluster(ind2)).type ==1 

            x = cluster(plot_cluster(ind2)).pos_c_BS(1);
            y = cluster(plot_cluster(ind2)).pos_c_BS(2);
            z = cluster(plot_cluster(ind2)).pos_c_BS(3);
            rotate(1) = cluster(plot_cluster(ind2)).Phi_c_BS;
            rotate(2) = cluster(plot_cluster(ind2)).Theta_c_BS;

            center = [x y z];
            r_x = cluster(plot_cluster(ind2)).a_c_BS;
            r_y = cluster(plot_cluster(ind2)).b_c_BS;
            r_z = cluster(plot_cluster(ind2)).h_c_BS;
            radius = [r_x r_y r_z];            

            weight = channel{ind1}.a_VR(ind2)*cluster_a_color; % The weight is according to attenuation

            draw_ellipsoid(center,radius,rotate,weight,alpha_val,'',figure_h)
            axis ([-2 2 -2 2 -0.1 0.1]*paraEx.net_radii)
            hold on
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Draw MPCs for cluster
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            for ind3 = 1:paraSt.n_mpc

                center = squeeze(mpc(plot_cluster(ind2)).pos_BS(ind3,:));
                center = center(:)';
                radius = paraEx.net_radii/500*[1 1 1]; % Drawing radius for MPCs
                rotate = [0 0];

%Obs!                weight = abs(channel{ind1}.h((ind2-1)*paraSt.n_mpc+ind3,6));
                plot3(center(1),center(2),center(3),'k.','markerSize',14,'lineWidth',3);
%Obs!                draw_ellipsoid(center,radius,rotate,weight,alpha_val,'',figure_h);

%Obs!                draw_line(MS(ind1).pos,center,'--r',figure_h); % Connect the MS to MPCs
%Obs!                draw_line(BS(ind1).pos,center,'--r',figure_h); % Connect the BS to MPCs
            end

        elseif cluster(plot_cluster(ind2)).type ==2 % Twin cluster
            x_BS = cluster(plot_cluster(ind2)).pos_c_BS(1);
            y_BS = cluster(plot_cluster(ind2)).pos_c_BS(2);
            z_BS = cluster(plot_cluster(ind2)).pos_c_BS(3);
            center_BS = [x_BS y_BS z_BS];
            r_x = cluster(plot_cluster(ind2)).a_c_BS;
            r_y = cluster(plot_cluster(ind2)).b_c_BS;
            r_z = cluster(plot_cluster(ind2)).h_c_BS;
            radius = [r_x r_y r_z];
            rotate(1) = cluster(plot_cluster(ind2)).Phi_c_BS;
            rotate(2) = cluster(plot_cluster(ind2)).Theta_c_BS;

            weight = channel{ind1}.a_VR(ind2)*cluster_a_color;
            draw_ellipsoid(center_BS,radius,rotate,weight,alpha_val,'',figure_h)
            axis ([-2 2 -2 2 -0.1 0.1]*paraEx.net_radii)
            hold on

            x_MS = cluster(plot_cluster(ind2)).pos_c_MS(1);
            y_MS = cluster(plot_cluster(ind2)).pos_c_MS(2);
            z_MS = cluster(plot_cluster(ind2)).pos_c_MS(3);
            center_MS = [x_MS y_MS z_MS];
            r_x = cluster(plot_cluster(ind2)).a_c_MS;
            r_y = cluster(plot_cluster(ind2)).b_c_MS;
            r_z = cluster(plot_cluster(ind2)).h_c_MS;
            radius = [r_x r_y r_z];
            rotate(1) = cluster(plot_cluster(ind2)).Phi_c_MS;
            rotate(2) = cluster(plot_cluster(ind2)).Theta_c_MS;
            
            weight = channel{ind1}.a_VR(ind2)*cluster_a_color; % The weight is according to attenuation
            draw_ellipsoid(center_MS,radius,rotate,weight,alpha_val,'',figure_h)

            draw_line(center_BS,center_MS,'.-b',figure_h);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Draw MPCs for clusters
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            for ind3 = 1:paraSt.n_mpc

                center_BS = squeeze(mpc(plot_cluster(ind2)).pos_BS(ind3,:));
                center_BS = center_BS(:)';                
                center_MS = squeeze(mpc(plot_cluster(ind2)).pos_MS(ind3,:));
                center_MS = center_MS(:)';
                
                radius = paraEx.net_radii/500*[1 1 1]; % Drawing radius for MPCs
                rotate = [0 0];

                weight = abs(channel{ind1}.h((ind2-1)*paraSt.n_mpc+ind3,6));

                draw_ellipsoid(center_BS,radius,rotate,weight,alpha_val,'',figure_h);
                draw_ellipsoid(center_MS,radius,rotate,weight,alpha_val,'',figure_h);
                
                draw_line(BS(ind1).pos,center_BS,'--r',figure_h); % Connect the BS to MPCs                
                draw_line(MS(ind1).pos,center_MS,'--r',figure_h); % Connect the MS to MPCs
            end
        end
    end
    
    if paraSt.BSLocal        
        x_BS = BS(ind1).cluster_local.pos_c_BS(1);
        y_BS = BS(ind1).cluster_local.pos_c_BS(2);
        z_BS = BS(ind1).cluster_local.pos_c_BS(3);
        center_BS = [x_BS y_BS z_BS];
        r_x = BS(ind1).cluster_local.a_c_BS;
        r_y = BS(ind1).cluster_local.b_c_BS;
        r_z = BS(ind1).cluster_local.h_c_BS;
        radius = [r_x r_y r_z];
        rotate(1) = 0;
        rotate(2) = 0;
        
        weight = 1*cluster_a_color;
        %draw_ellipsoid(center_BS,radius,rotate,weight,alpha_val,'',figure_h)
        axis ([-2 2 -2 2 -0.1 0.1]*paraEx.net_radii)
        hold on
        
        for ind3 = 1:paraSt.n_mpc
                center = squeeze(BS(ind1).mpc_local.pos_BS(ind3,:));
                center = center(:)';

                radius = paraEx.net_radii/500*[1 1 1]; % Drawing radius for MPCs
                rotate = [0 0]; % Depicted as balls

                weight = abs(channel{ind1}.h((ind2-1)*paraSt.n_mpc+ind3,6));

                draw_ellipsoid(center,radius,rotate,weight,alpha_val,'',figure_h);

                draw_line(MS(ind1).pos,center,'--r',figure_h); % Connect the MS to MPCs
                draw_line(BS(ind1).pos,center,'--r',figure_h); % Connect the BS to MPCs
        end
    end
    
    if paraSt.MSLocal
        x_MS = MS(ind1).cluster_local.pos_c_MS(1);
        y_MS = MS(ind1).cluster_local.pos_c_MS(2);
        z_MS = MS(ind1).cluster_local.pos_c_MS(3);
        center_MS = [x_MS y_MS z_MS];
        r_x = MS(ind1).cluster_local.a_c_BS;
        r_y = MS(ind1).cluster_local.b_c_BS;
        r_z = MS(ind1).cluster_local.h_c_BS;
        radius = [r_x r_y r_z];        
        rotate(1) = 0;
        rotate(2) = 0;
        
        weight = 1*cluster_a_color;
        %draw_ellipsoid(center_MS,radius,rotate,weight,alpha_val,'',figure_h)
        axis ([-2 2 -2 2 -0.1 0.1]*paraEx.net_radii)
        hold on
        
        for ind3 = 1:paraSt.n_mpc
                center = squeeze(MS(ind1).mpc_local.pos_MS(ind3,:));
                center = center(:)';

                radius = paraEx.net_radii/500*[1 1 1]; % Drawing radius for MPCs
                rotate = [0 0]; % Depicted as balls

                weight = abs(channel{ind1}.h((ind2-1)*paraSt.n_mpc+ind3,6));

                draw_ellipsoid(center,radius,rotate,weight,alpha_val,'',figure_h);

                draw_line(MS(ind1).pos,center,'--r',figure_h); % Connect the MS to MPCs
                draw_line(BS(ind1).pos,center,'--r',figure_h); % Connect the BS to MPCs
        end
    end
         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Draw cell radius
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    degr = [1:360]/180*pi;
    xx = paraEx.net_radii*cos(degr)+BS(1).pos(1);
    yy = paraEx.net_radii*sin(degr)+BS(1).pos(2);
    zz = zeros(size(xx));
    figure(figure_h)
    plot3(xx,yy,zz,'--k','LineWidth',2)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Text and view
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     viewP = [   0.9061    0.4232    0.0000   -0.6646;
%                -0.3549   0.7599    0.5446   -0.4748;
%                -0.2305   0.4935   -0.8387    8.9481;
%                 0         0         0         1.0000];
%     view(viewP)
%     text(3996 ,-8672, 2890,'BS','FontSize',18);
%     figure(figure_h)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % View in azimuth plane
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    view([0 90])
    text(3996 ,-8672, 2890,'BS','FontSize',18);
    axis equal
%    axis square
    figure(figure_h)
    legend('BS','MS')
    setFontsize
    pause(0.5)
    %         movframe(ind1) = getframe(gcf);  %Uncomment if you need to
    %         generate the figures of each snapshot
    %         saveas(gcf,num2str(ind1),'jpg');
end

%Uncomment if you need to generate a movie
% save tryMovie
% 
% movie2avi(movframe,'test.avi', 'fps', 10, 'quality', 100)