\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
The updated version 3.1 is implemented by Sara Gunnarsson and Jose 
Flordelis, Lund University, Lund, Sweden. Massive MIMO scenario 
Indoor_CloselySpacedUser_2_6GHz has been added, together with some 
modifications to the VLA and MPC gain function extensions for Massive MIMO. 
Additionally, bug fixes and speed optimizations are also included.

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
The updated version 3.0 is implemented by Sara Gunnarsson and Jose 
Flordelis, Lund University, Lund, Sweden. Massive MIMO extensions are now 
fully backwards compatible, and the four scenarios (IndoorHall_5GHz, 
SemiUrban_300MHz, SemiUrban_VLA_2_6GHz, SemiUrban_CloselySpacedUser_2_6GHz) 
can be run. Some bug fixes and code restructuring have also been done.

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
The Massive MIMO-extended version was implemented by Xiang Gao, 
Lund University, Sweden, and was partly supported by the Seventh Framework 
Programme (FP7) of the European Union under grant agreement no. 619086 
(MAMMOET). Two new scenarios are implemented:
SemiUrban_CloselySpacedUser_2_6GHz: both NLOS and LOS as well as both single and 
multiple link MIMO simulations are supported.
SemiUrban_VLA_2_6GHz: both NLOS and LOS as well as both single and multiple link 
MIMO simulations are supported.

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
The updated version 2.3.5 is implemented by Meifang Zhu, Lund University, Sweden.
The get_H and get_dipole_G have been updated

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
The updated version 2.3.4 is implemented by Meifang Zhu, Lund University, Sweden.
The update_chan has been updated

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
The updated version 2.3.3 is implemented by Meifang Zhu, Lund University, Sweden.
Two scenarios are implemented so far:
SemiUrban_300MHz: both NLOS and LOS single link MIMO simulations 
are supported. For outdoor LOS multiple link MIMO simulation is supported as well. 
IndoorHall_5GHz: OLOS single link MIMO simulation is supported.

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
Previous version history:
1.00   25/05/2005 by HH (ftw)                    
1.10   08/07/2006 by HH (Eurecom)
1.2.0  22/09/2008 by LIU (UCL) 
1.2.1  29/10/2008 by LIU (UCL)
1.3    19/01/2010 by LIU (UCL)
2.1    04/02/2010 by LIU (UCL)   
2.2    08/04/2010 by LIU (UCL) 
2.3    15/03/2013 by Meifang Zhu (Lund University)

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
If you use the COST 2100 channel model for publications, please refer to: 
L. Liu, J. Poutanen, F. Quitin, K. Haneda, F. Tufvesson, P. De Doncker,
P. Vainikainen and C. Oestges, “The COST 2100 MIMO channel model,”
IEEE Wireless Commun., vol 19, issue 6, pp 92-99, Dec. 2012.

Further details about the COST 2100 channel model can be found in:
Roberto Verdone (Editor), Alberto Zanella (Editor)
Pervasive Mobile and Ambient Wireless Communications Pervasive Mobile 
and Ambient Wireless Communications, ISBN 978-1-4471-2315-6,
Springer, 2012. 

If you use the SemiUrban_300MHz scenario, further information can be found in:
1. Meifang Zhu, Gunnar Eriksson, and Fredrik Tufvesson, 
"The COST 2100 Channel Model: Parameterization and Validation 
Based on Outdoor MIMO Measurements at 300 MHz", 
IEEE Transactions on Wireless Commun..
2. Meifang Zhu and Fredrik Tufvesson, "Virtual Multi-link Propagation 
Investigation of an Outdoor Scenario At 300 MHz," Proceedings of the 
7th European Conference on Antennas and Propagation (EUCAP), Gothenburg,
Sweden, April 2013.

If you use the IndoorHall_5GHz scenario, further information can be found in:
1. V. M. Kolmonen, P. Almers, J. Salmi, J. Koivunen, A. Richter,
F. Tufvesson, A. Molisch, P. Vainikainen, "A dynamic dual-link 
wideband MIMO channel sounder for 5.3 GHz," IEEE Transactions on 
Instrumentation and Measurement, Vol. 59, No. 4, pp. 873-883, 2010.
2. J. Poutanen, K. Haneda, L. Lin, C. Oestges, F. Tufvesson , 
P. Vainikainen, "Parameterization of the COST 2100 MIMO channel 
modeling in indoor scenarios," Proceedings of the 5th European 
Conference on Antennas and Propagation (EUCAP), Rome, Italy, 
pp. 3606-3610, April 2011.

If you use the SemiUrban_CloselySpacedUser_2_6GHz scenario, further information can be found in:
1. Flordelis, J., Gao, X., Dahman, G., Rusek, F., Edfors, O., & Tufvesson, F. (2015). 
Spatial Separation of CloselySpaced Users in Measured Massive Multi-User MIMO Channels. 
In 2015 IEEE International Conference on Communications (ICC) (pp. 1441-1446). 
IEEE--Institute of Electrical and Electronics Engineers Inc..
2. Gao, X., Flordelis, J., Dahman, G., Tufvesson, F., & Edfors, O. (2015). 
Massive MIMO Channel Modeling - Extension of the COST 2100 Model. Paper presented at 
Joint NEWCOM/COST Workshop on Wireless Communications (JNCW), Barcelona, Spain.
3. Bourdoux, A., Desset, C., van der Perre, L., Dahman, G., Edfors, O., Flordelis, J., 
Medbo, J. (2015). D1.2 MaMi Channel Characteristics: Measurement Results. MAMMOET.

If you use the SemiUrban_VLA_2_6GHz scenario, further information can be found in:
1. Gao, X., Flordelis, J., Dahman, G., Tufvesson, F., & Edfors, O. (2015). 
Massive MIMO Channel Modeling - Extension of the COST 2100 Model. Paper presented at 
Joint NEWCOM/COST Workshop on Wireless Communications (JNCW), Barcelona, Spain.
2. Bourdoux, A., Desset, C., van der Perre, L., Dahman, G., Edfors, O., Flordelis, J., 
Medbo, J. (2015). D1.2 MaMi Channel Characteristics: Measurement Results. MAMMOET.

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
The script "demo_model.m" provides an example for testing the COST 2100 model. It selects network, 
link, antenna, band and scenario and initiates the simulation.

The output of demo_model.m is a channel function, with size dependent on the parameter choice and setup: 

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

5) MIMO_Cyl_patch: Transfer function for a synthetic pattern of a cylindrical array 
 with 128 antennas.
 get_IR_Cyl_patch: users have to set up the frequency separation, delta_f, as well as
 provide the full antennas response at the BS and MS sides. The antennas response mat 
 file have to be the same format as the 'BS_Cyl_AntPattern.mat' and 'MS_AntPattern_User.mat' 
 files.

6) MIMO_VLA_omni: Transfer function for a physically large array with 128 omni-directional 
 antennas, with lambda/2 inter-element separation, and MS with omni-directional antenna.
 create_IR_omni_MIMO_VLA: users have to set up the frequency separation, delta_f
 
Modifications are listed in the document changelog.

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
The code is developed under GPL. The COST2100 code was originally written by Lingfeng Liu, Universtié catholique de Louvain (UCL). 

Email: lingfeng.liu@uclouvain.be
Tel: +32 (0) 10 47 81 05
Address: Bâtiment Maxwell, Place du Levent 3, 1348 Louvain-la-Neuve, Belgium

1, Function list
calc_dist
calc_pathloss
cost2100
create_IR_Cyl
create_IR_Cyl_EADF
create_IR_omni
create_IR_omni_MIMO
create_IR_omni_MIMO_VLA
demo_model
draw_circ
draw_ellpsoid
draw_line
eadf
get_BS_VR_para
get_channel
get_channel_los
get_cluster
get_cluster_local
get_dipole_G
get_dmc
get_H
get_mpc
get_MS_VR
get_para
get_VRLOS
get_VRtable
lognrnd_own
mexprnd_own
poissrnd_own
rotate_matrix
setFontsize
update_chan
visual_pddp
visual_channel
visualize_channel_env

2, Other files list
antSample.mat: antenna array radiation pattern sample file
BS_Cyl_AntPattern.mat: synthetic antenna array radiation pattern sample file
BS_Cyl_EADF: EADF file to be used as antenna pattern for the cylindrical array
changelog.txt: list of modifications
GPL.txt: GNU general public license
MS_AntPattern_User: antenna radiation pattern sample file
unit.txt: list of parameter unit defined in the model

3, misc
poissrnd_own_parallel
test_poissrnd_own

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
3, How to run the code

The main routine of the code is cost2100

The visualization routines are named as 'visual_*' or 'visualize.*'

To understand each function, type help function_name in Matlab

The test script is demo_model.m

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
4, How to customize the channel implementation

Please read get_para for more instructions

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
5, Copyright
This file is a part of the COST2100 channel model.

This program, the COST2100 channel model, is free software: you can 
redistribute it and/or modify it under the terms of the GNU General Public 
License as published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
for more details.

If you use it for scientific purposes, please consider citing it in line 
with the description in the Readme-file, where you also can find the 
contributors.

You should have received a copy of the GNU General Public License along 
with this program. If not, see <http://www.gnu.org/licenses/>.

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