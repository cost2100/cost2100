function draw_ellipsoid(center, radius, rotate, weight, alpha_val, control, figure_h)
%DRAW_ELLIPSOID Drawing function to draw 3D ellipsoid
%draw_ellipsoid(center, radius, rotate, weight, alpha_val, control, figure_h)
%draw a ellipsoid in figure_h at position center with radius and tilt angle 
%rotate. The color is determined by weight and the transparency is 
%determined by alpha_val. control is the control option

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

[x y z] = sphere;

x = reshape(x,1,length(x(:)));
y = reshape(y,1,length(y(:)));
z = reshape(z,1,length(z(:)));
tmp = [x' y' z'];
tmp = tmp*diag(radius)*rotate_matrix(rotate(1),rotate(2));
tmp = tmp + repmat(center,length(x(:)),1);

x = tmp(:,1);
y = tmp(:,2);
z = tmp(:,3);
x = reshape(x,sqrt(length(x(:))),sqrt(length(x(:))));
y = reshape(y,sqrt(length(y(:))),sqrt(length(y(:))));
z = reshape(z,sqrt(length(z(:))),sqrt(length(z(:))));
c = weight.*ones(size(x));

surface(x,y,z,c);
alpha(alpha_val);
shading flat
lighting gouraud