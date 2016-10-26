#!/bin/tcsh -f

#@ s = 0

#foreach i (`seq 0 9999`)
foreach i (`seq 1 10000`)

setenv size 0


if ((-e /mss/clas/e1e/production/simulation_2pi/sim_conv_w_145_155_gleb/t21_${i}.root)) then
setenv size `cat /mss/clas/e1e/production/simulation_2pi/sim_conv_w_145_155_gleb/t21_${i}.root | grep size | sed -e 's/size=//g'`
endif
if ($size < 19000000) then
#jremove /mss/clas/e1e/production/simulation_2pi/sim_conv_w_145_155_gleb/t21_${i}.root
echo $i'  '$size
#echo t21_${i}.root
endif



#if ($size < 200000000) then
#@ s++
#echo $i

#endif

#if (!(-e /mss/home/skorodum/e1e/simulation/sim_conv_w_15_17_gleb/t21_${i}.root)) then
#if (!(-e /mss/home/gleb/e1e/sim_2014_w_1275_1525/recsis/recsis${i}.bos)) then
#if (!(-e /mss/home/gleb/e1e/sim_2015_w_15_17/gsim/gsim${i}.bos)) then
#echo $i
#@ s++
#endif

end

#echo $s
