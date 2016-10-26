#!/bin/tcsh -f


foreach i (`seq 0 10000`)

#foreach file ( /cache/mss/home/gleb/e1e/sim2014_w_1675_1825/goa/goa_out*.hbook )

setenv size 0

#if ((-e /mss/home/skorodum/e1e/simulation/sim_w_1275_1525_gleb_Sep2015/nt10/nt10_${i}.root)) then
setenv size `cat  /mss/home/skorodum/e1e/simulation/sim_conv_w_1275_1525_gleb/t21_${i}.root | grep size | sed -e 's/size=//g'`

#endif

if ($size < 203000000) then
echo $size" "$i
#jremove /mss/home/skorodum/e1e/simulation/sim_conv_w_1275_1525_gleb/t21_${i}.root
#jcache remove /mss/home/skorodum/e1e/simulation/sim_conv_w_1275_1525_gleb/t21_${i}.root
endif



end
