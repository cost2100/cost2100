\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
The updated version 2.3.5 is implemented by Meifang Zhu, Lund University, Sweden.
The get_H and get_dipole_G have been updated

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
The updated version 2.3.4 is implemented by Meifang Zhu, Lund University, Sweden.
The update_chan has been updated

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
The updated version 2.3.3 is implemented by Meifang Zhu, Lund University, Sweden.
Two scenarios are implemented so far
SemiUrban_300MHz: both NLOS and LOS single link MIMO simulations 
are supported. For outdoor LOS multiple link MIMO simulation is supported as well. 
IndoorHall_5GHz: OLOS single link MIMO simulation is supported.

If you use the COST 2100 channel model for publications, please refer to: 
L. Liu, J. Poutanen, F. Quitin, K. Haneda, F. Tufvesson, P. De Doncker,
P. Vainikainen and C. Oestges, “The COST 2100 MIMO channel model,”
IEEE Wireless Commun., vol 19, issue 6, pp 92-99, Dec. 2012.

Further details about the COST 2100 channel model can be found in:
Roberto Verdone (Editor), Alberto Zanella (Editor)
Pervasive Mobile and Ambient Wireless Communications Pervasive Mobile 
and Ambient Wireless Communications, ISBN 978-1-4471-2315-6,
Springer, 2012. 

If you use the SemiUrban_300MHz scenario, 
further information can be found in:
1. Meifang Zhu, Gunnar Eriksson, and Fredrik Tufvesson, 
"The COST 2100 Channel Model: Parameterization and Validation 
Based on Outdoor MIMO Measurements at 300 MHz", 
IEEE Transactions on Wireless Commun..
2. Meifang Zhu and Fredrik Tufvesson, "Virtual Multi-link Propagation 
Investigation of an Outdoor Scenario At 300 MHz," Proceedings of the 
7th European Conference on Antennas and Propagation (EUCAP), Gothenburg,
Sweden, April 2013.

If you use the IndoorHall_5GHz scenario, 
further information can be found in:
1. V. M. Kolmonen, P. Almers, J. Salmi, J. Koivunen, A. Richter,
F. Tufvesson, A. Molisch, P. Vainikainen, "A dynamic dual-link 
wideband MIMO channel sounder for 5.3 GHz," IEEE Transactions on 
Instrumentation and Measurement, Vol. 59, No. 4, pp. 873-883, 2010.
2. J. Poutanen, K. Haneda, L. Lin, C. Oestges, F. Tufvesson , 
P. Vainikainen, "Parameterization of the COST 2100 MIMO channel 
modeling in indoor scenarios," Proceedings of the 5th European 
Conference on Antennas and Propagation (EUCAP), Rome, Italy, 
pp. 3606-3610, April 2011.



The script "demo_model.m" provides an example for testing the COST 2100 model. It selects scenario, link type, and initiates the simulation.
The output of demo_model.m is channel functions, with size dependent on the parameter choice and setup: 

1) SISO_omni: Transfer function for SISO omni-directional antenna
create_IR_omni: users have to set up the frequency separation, delta_f
 
2) MIMO_omni: Transfer function for MIMO omini-directional antenna
 create_IR_omni_MIMO: users have to set up the frequency separation, delta_f.
 Only 2 by 2 MIMO system is implemented.
 
3) MIMO_dipole: Transfer function for a theoretical antenna response for 
 any size of lambda/2-spaced linear dipole antenna arrays. An Ntx-by-Nrx theoretical 
 antenna array response is generated and the correponding 
 channel transfer function is simulated.
 
4) MIMO_measured: Transfer function for any measured MIMO antenna response
 get_H: users have to provide the full antenna response at the BS and 
 MS sides, and also the rotation of the antenna arrays. The antenna 
 response mat file have to be the same format as 'antSample.mat' file.


The modifications are listed in document changelog.

Meifang Zhu, 2013.02.06


\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 The code is developed under GPL. The COST2100 code was originally written by Lingfeng Liu, Universtié catholique de Louvain (UCL). 

email: lingfeng.liu@uclouvain.be
tel: +32 (0) 10 47 81 05
address: Bâtiment Maxwell, Place du Levent 3, 1348 Louvain-la-Neuve, Belgium


1,Function list
calc_dist
calc_pathloss
cost2100
demo_model
draw_circ
draw_ellpsoid
get_channel
get_channel_los
get_cluster_local
get_common
get_cluster
get_dmc
get_H
get_IR
get_mpc
get_para
get_VR
get_VRLOS
get_VRtable
rotate_matrix
setFontsize
update_chan
visual_pddp
visual_channel

2,Other files list
unit.txt: list of parameter unit defined in the model
antSample.mat: antenna array radiation pattern sample file

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
3, How to run the code

The main routine of the code is cost2100

The visualization routines are named as 'visual_*'

To understand each function,  type help function_name in matlab

The test script is demo_model.m
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
4, How to customize the channel implementation

Please read get_para for more instructions

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
5, Copyright
Copyright (C)2008 LIU Ling-Feng, Université catholique de Louvain, Belgium
This program, cost2100, is free software: you can redistribute it and/or 
modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or(at your 
option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

6, Acknowledgements
Thanks for Helmut Hofstetter for his previous code on the COST 273 
channel model which inspires me a lot.

Thanks for Claude Oestges and Nicolai Czink for their feedbacks and theory 
supports.

Thanks for Katsuyuki Haneda, Juho Puotanen, Fredrik Tufvesson, and Meifang 
Zhu for their active participations and code testing.

Also thanks for my wife Qin, always be patient during my coding time.

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
7, Some further references

[1] L.M. Correia, Mobile Broadband Multimedia Networks. Academic Press, 2006.
[2] N. Czink and C. Oestges, “The COST 273 MIMO channel model: Three kinds of clusters,” IEEE 10th Int. Sym., ISSSTA’08, pp. 282–286, 2008.
[3] L. Liu, N. Czink, and C. Oestges, Implementing the COST 273 MIMO channel model, NEWCOM 2009