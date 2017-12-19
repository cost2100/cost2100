function  draw_line( a,b,control,fig_h )
%DRAW_LINE Drawing function of draw a line between point a and b
%draw_line( a,b,control,fig_h ) plot a line between point a and b in fig_h
%drawing control control

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

if nargin == 4
    figure(fig_h);
end
if length(a) == 3
    plot3([a(1) b(1)],[a(2) b(2)],[a(3) b(3)],control);
else
    plot([a(1) b(1)],[a(2) b(2)],control);
end
