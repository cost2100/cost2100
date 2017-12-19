function draw_circ(center, r, control, fig_h)
%DRAW_CIRC Drawing function to draw 2D circles
%draw_circ(center, r, control, fig_h) draw a circle at position center with
%radius r in fig_h. control is the control option

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (C)2008 LIU Ling-Feng ,Université catholique de Louvain, Belgium
%This file is part of cost2100_model.

%cost2100_model is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.

%cost2100_model is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.

%You should have received a copy of the GNU General Public License
%along with cost2100_model.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = [0:360]/180*pi;

x = center(1)+r*cos(phi);
y = center(2)+r*sin(phi);

if nargin ==4
figure(fig_h)
end
plot(x,y,control,'LineWidth',2);