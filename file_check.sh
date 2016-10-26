#!/bin/tcsh -f

#@ s = 0

#foreach i (`seq 0 9999`)
foreach i (`seq 1 10000`)
setenv size 0
setenv size `ls -all /cache/clas/e1e/production/simulation_2pi/gleb_2016/low_w/low_w${i}.root  | awk '{print $5}'`
#echo $size
#if ((-e /mss/home/skorodum/e1e/simulation/sim_conv_w_165_185_gleb/t21_${i}.root)) then
#setenv size `cat /mss/home/skorodum/e1e/simulation/sim_conv_w_165_185_gleb/t21_${i}.root | grep size | sed -e 's/size=//g'`
if ($size < 230000000) then
jremove /mss/home/skorodum/e1e/simulation/sim_conv_w_165_185_gleb/t21_${i}.root
rm -f /cache/clas/e1e/production/simulation_2pi/gleb_2016/low_w/low_w${i}.root
jcache get /mss/clas/e1e/production/simulation_2pi/gleb_2016/low_w/low_w${i}.root
echo $i
#echo t21_${i}.root
endif

#endif
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
