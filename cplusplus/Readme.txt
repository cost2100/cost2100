This is a C++ implementation for the COST2100 channel model. Its support single and multiple link simulations. It supports indoor OLOS scenario simulation for now. The code is developed by Weiwei Jiang under GPL. The code is currently maintained by Meifang Zhu, Lund University. For any problem and question, please directly contact us via:

email: meifang.zhu@eit.lth.se
tel: +46 46 2229246
address: Institutionen för Elektro- och informationsteknik, Lunds universitet, Box 118, 221 00 Lund, Sweden.

This C++ implementation have been developed based on the current MATLAB framework and TD(10)12092.

1. function list
COST2100_Channel : interface for the system level simulation. You can integrate this function into any of your C++ project. The return variable of this function is the corresponding channel matrix.

COST2100_Specification:  channel specification including all the main majors of the channel model. In the specification, you can modify you single or multiple link parameters, such as common cluster ratio. Also you can develop a scenario you wanted, and now only indoor OLOS scenario has been developed. Antenna pattern can be modified as well to generate proper channel matrix.

2. How to run the code
a. modify the function COST2100_Spectification if you have a different scenario besides indoor OLOS scenario.

b. Integrate the function COST2100_Channel into your simulation chain, specify the input for the function and channel matrix can be generated directly.

3. Copyright
Copyright (C)2008 Weiwei Jiang, Lund University, Sweden.
This program, cost2100, is a free software: you can redistribute it and/or 
modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or(at your 
option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
