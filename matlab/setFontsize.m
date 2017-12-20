% A script to change the fontsize of every notation in a figure

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

FontSize = 24;

h = gcf;
c = findobj(gcf);
for n = 1:length(c),
    try
        set(c(n), 'FontSize', FontSize)
    catch
    end
%     if strcmp( get(c(n),'Type'), 'axes' )
%         set(get(c(2),'XLabel'), 'FontSize', FontSize)
%         set(get(c(2),'YLabel'), 'FontSize', FontSize)
%         set(get(c(2),'ZLabel'), 'FontSize', FontSize)
%         set(get(c(2),'Title'), 'FontSize', FontSize)
%     end
end

try
    h = get(gca, 'title');
    set(h, 'FontSize', FontSize)
    h = get(gca, 'xlabel');
    set(h, 'FontSize', FontSize)
    h = get(gca, 'ylabel');
    set(h, 'FontSize', FontSize)
    h = get(gca, 'zlabel');
    set(h, 'FontSize', FontSize)
catch
end

% set(fighandle,'paperpositionmode','auto')
% switch(filename(end-2:end))
%     case 'jpg'
%         print(h, '-djpeg', filename)
%     case 'eps'
%         print(h, '-depsc2', filename)
%     otherwise
%         error('unknown file format')
% end