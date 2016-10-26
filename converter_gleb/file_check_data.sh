#!/bin/tcsh -f

@ s = 0

foreach i (`seq 1 185`)
#foreach i (`seq 1 50000`)

#setenv size `cat /mss/home/skorodum/e1e/simulation/sim_conv_w_15_17_gleb/t21_${i}.root | grep size | sed -e 's/size=//g'`
#echo $size



if ((-e /mss/home/skorodum/e1e/GLEB_DATA_E1E_CONVERTED/t21_data${i}.root)) then
setenv size `cat  /mss/home/skorodum/e1e/GLEB_DATA_E1E_CONVERTED/t21_data${i}.root | grep size | sed -e 's/size=//g'`
echo $i $size
endif
#if ($size < 200000000) then


#@ s++
#echo $i

#if ($size < 500000000) then
#@ s++
#echo $i
#endif

#if (!(-e /mss/home/skorodum/e1e/simulation/sim_conv_w_15_17_gleb/t21_${i}.root)) then

#if (!(-e /mss/home/gleb/e1e/sim_2015_w_15_17/recsis/recsis${i}.bos)) then
if (!(-e /mss/home/skorodum/e1e/GLEB_DATA_E1E_CONVERTED/t21_data${i}.root)) then
echo $i
@ s++
endif

end

echo $s
