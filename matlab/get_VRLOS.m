function VRLOS = get_VRLOS( BS, paraEx, paraSt)
%GET_VRLOS Generate the LOS VR for each BS
%Default call: VRLOS = get_VRLOS( BS, paraEx, paraSt)
%------
%Input:
%------
%BS: BS information
% .pos: BS position
%paraEx, paraSt: external and stochastic parameters
%------
%Output:
%------
%VRLOS: LOS VR position, VR([x y])
%
%See also: cost2100, get_para, get_VR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (C)2008 LIU Ling-Feng, ICTEAM, UCL, Belgium 
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
x = rand*(paraSt.d_co-paraSt.r_l);
y = rand*(paraSt.d_co-paraSt.r_l);

while (sqrt(x^2+y^2)>paraSt.d_co)
    x = rand*(paraSt.d_co-paraSt.r_l);
    y = rand*(paraSt.d_co-paraSt.r_l);
end

VRLOS = [x+BS.pos(1) y+BS.pos(2)];