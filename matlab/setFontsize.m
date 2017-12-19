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

%A script to change the fontsize of every notation in a figure
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