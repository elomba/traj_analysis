V34 :0x24 gpucells
9 Cells.f90 S624 0
12/16/2023  00:42:25
use g_types public 0 indirect
use linkcell_d public 0 direct
use iso_c_binding public 0 indirect
use nvf_acc_common public 0 indirect
use cudafor_lib_la public 0 indirect
use cudafor_la public 0 direct
use gpcodes public 0 direct
enduse
D 58 26 648 8 647 7
D 67 26 651 8 650 7
D 76 26 648 8 647 7
D 97 26 745 8 744 7
D 5147 23 9 2 15286 15293 0 0 1 0 0
 0 15288 11 11 15289 15289
 0 15291 15289 11 15292 15292
D 5150 23 9 1 11 15289 0 0 1 0 0
 0 15288 11 11 15289 15289
D 5153 23 9 2 15294 15301 0 0 1 0 0
 0 15296 11 11 15297 15297
 0 15299 15297 11 15300 15300
D 5156 23 9 1 11 15297 0 0 1 0 0
 0 15296 11 11 15297 15297
S 624 24 0 0 0 9 1 0 5013 10005 0 A 0 0 0 0 B 0 211 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 211 0 0 0 0 0 0 gpucells
R 647 25 7 iso_c_binding c_ptr
R 648 5 8 iso_c_binding val c_ptr
R 650 25 10 iso_c_binding c_funptr
R 651 5 11 iso_c_binding val c_funptr
R 685 6 45 iso_c_binding c_null_ptr$ac
R 687 6 47 iso_c_binding c_null_funptr$ac
R 688 26 48 iso_c_binding ==
R 690 26 50 iso_c_binding !=
R 744 25 6 nvf_acc_common c_devptr
R 745 5 7 nvf_acc_common cptr c_devptr
R 751 6 13 nvf_acc_common c_null_devptr$ac
R 789 26 51 nvf_acc_common =
S 20927 23 5 0 4 0 20932 624 159303 0 0 A 0 0 0 0 B 0 219 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 gpu_graph_cell
S 20928 7 3 1 0 5147 1 20927 156370 808204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 r
S 20929 6 1 1 0 6 1 20927 158330 808004 3000 A 0 0 0 0 B 0 219 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nmol
S 20930 6 1 1 0 6 1 20927 2375 808004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dim
S 20931 7 3 1 0 5150 1 20927 155649 808204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 sidel
S 20932 14 5 0 4 0 1 20927 159303 200 400000 A 0 0 0 0 B 0 219 0 0 0 0 0 6212 6 0 0 0 0 0 0 0 0 0 0 0 0 219 0 624 0 0 0 0 gpu_graph_cell gpu_graph_cell 
F 20932 6 20928 20937 20938 20939 20931 20940
S 20933 6 1 0 0 7 1 20927 159318 40808006 3000 A 0 0 0 0 B 0 223 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_15552
S 20934 6 1 0 0 7 1 20927 72959 40808006 3000 A 0 0 0 0 B 0 223 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_15554
S 20935 6 1 0 0 7 1 20927 72969 40808006 3000 A 0 0 0 0 B 0 223 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_15557
S 20936 6 1 0 0 7 1 20927 159328 40808006 3000 A 0 0 0 0 B 0 223 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_15559
S 20937 6 3 1 0 6 1 20927 158384 800004 7000 A 0 0 0 0 B 0 219 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _V_nmol
S 20938 6 3 1 0 6 1 20927 158392 800004 7000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _V_dim
S 20939 1 3 0 0 9 1 20927 158399 4 7000 A 0 0 0 0 B 0 219 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _V_rcl2
S 20940 1 3 1 0 6 1 20927 158407 4 7000 A 0 0 0 0 B 0 219 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _V_jmin
S 20941 23 5 0 4 0 20946 624 159338 0 0 A 0 0 0 0 B 0 264 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 gpu_adj_cell
S 20942 7 3 1 0 5153 1 20941 156370 808204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 r
S 20943 6 1 1 0 6 1 20941 158330 808004 3000 A 0 0 0 0 B 0 264 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nmol
S 20944 6 1 1 0 6 1 20941 2375 808004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dim
S 20945 7 3 1 0 5156 1 20941 155649 808204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 sidel
S 20946 14 5 0 4 0 1 20941 159338 200 400000 A 0 0 0 0 B 0 264 0 0 0 0 0 6219 5 0 0 0 0 0 0 0 0 0 0 0 0 264 0 624 0 0 0 0 gpu_adj_cell gpu_adj_cell 
F 20946 5 20942 20951 20952 20953 20945
S 20947 6 1 0 0 7 1 20941 72979 40808006 3000 A 0 0 0 0 B 0 268 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_15560
S 20948 6 1 0 0 7 1 20941 159351 40808006 3000 A 0 0 0 0 B 0 268 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_15562
S 20949 6 1 0 0 7 1 20941 159361 40808006 3000 A 0 0 0 0 B 0 268 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_15565
S 20950 6 1 0 0 7 1 20941 159371 40808006 3000 A 0 0 0 0 B 0 268 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_15567
S 20951 6 3 1 0 6 1 20941 158384 800004 7000 A 0 0 0 0 B 0 264 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _V_nmol
S 20952 6 3 1 0 6 1 20941 158392 800004 7000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _V_dim
S 20953 1 3 0 0 9 1 20941 158399 4 7000 A 0 0 0 0 B 0 264 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _V_rcl2
A 68 1 0 0 0 58 685 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 71 1 0 0 0 67 687 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 141 1 0 0 0 97 751 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 15286 1 0 0 11861 7 20936 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 15287 1 0 0 7186 6 20930 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 15288 7 0 0 13668 7 15287 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 15289 1 0 0 14885 7 20933 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 15290 1 0 0 7183 6 20929 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 15291 7 0 0 13671 7 15290 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 15292 1 0 0 14870 7 20934 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 15293 1 0 0 2509 7 20935 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 15294 1 0 0 7224 7 20950 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 15295 1 0 0 11871 6 20944 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 15296 7 0 0 13677 7 15295 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 15297 1 0 0 11867 7 20947 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 15298 1 0 0 11868 6 20943 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 15299 7 0 0 13680 7 15298 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 15300 1 0 0 11869 7 20948 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 15301 1 0 0 11872 7 20949 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
J 133 1 1
V 68 58 7 0
S 0 58 0 0 0
A 0 6 0 0 1 2 0
J 134 1 1
V 71 67 7 0
S 0 67 0 0 0
A 0 6 0 0 1 2 0
J 36 1 1
V 141 97 7 0
S 0 97 0 0 0
A 0 76 0 0 1 68 0
Z
