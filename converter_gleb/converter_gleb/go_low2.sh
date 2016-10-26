#!/bin/tcsh -f

@ ll = 1
foreach k (`seq  10895 5 10895`)

@ z = $k / 5 + 1
#echo $z

echo "PROJECT: e1e" > jsub_new
echo "TRACK: simulation" >>jsub_new
echo "JOBNAME: convert_low2_${k}" >>jsub_new
echo "MAIL: gleb@jlab.org" >>jsub_new
echo "SINGLE_JOB: TRUE" >>jsub_new
echo "MEMORY: 3000 MB" >>jsub_new
echo "DISK_SPACE: 15 GB" >>jsub_new
echo "OS:  centos65" >>jsub_new
@ k1 = $k + 1
@ k2 = $k + 2
@ k3 = $k + 3
@ k4 = $k + 4
set  a = `sed -n ${ll}p qqq`
@ ll++

echo "INPUT_FILES:/mss/home/skorodum/e1e/simulation/sim_w_145_155_gleb_Sep2015/nt10/nt10_${k}.root /mss/home/skorodum/e1e/simulation/sim_w_145_155_gleb_Sep2015/nt10/nt10_${k1}.root /mss/home/skorodum/e1e/simulation/sim_w_145_155_gleb_Sep2015/nt10/nt10_${k2}.root /mss/home/skorodum/e1e/simulation/sim_w_145_155_gleb_Sep2015/nt10/nt10_${k3}.root /mss/home/skorodum/e1e/simulation/sim_w_145_155_gleb_Sep2015/nt10/nt10_${k4}.root /volatile/clas/clase1/gleb/test/2pi_analysis_e1e/converter_gleb/converter_gleb/h10tot21" >>jsub_new
echo "COMMAND: /volatile/clas/clase1/gleb/test/2pi_analysis_e1e/converter_gleb/converter_gleb/execut_low2.sh ${k}" >>jsub_new
echo "OUTPUT_DATA: out.root" >>jsub_new
echo "OUTPUT_TEMPLATE: /mss/clas/e1e/production/simulation_2pi/gleb_2016/low_w_145_155/low_w_145_155${z}.root" >>jsub_new

#set  a = `sed -n ${z}p qqq`
#sed -n ${z}p qqq
#echo $a
#if (!(-e /mss/clas/e1e/production/simulation_2pi/gleb_2016/low_w/low_w${z}.root)) then
#/site/bin/jsub jsub_new
echo ${z}
#endif

#setenv size 0
#if ((-e /mss/clas/e1e/production/simulation_2pi/gleb_2016/low_w/low_w${z}.root)) then
#setenv size `sed -n '3p' /mss/clas/e1e/production/simulation_2pi/gleb_2016/low_w/low_w${z}.root | sed -e 's/size=//g'`
#endif

#if ($size < 230000000) then
#echo ${size}'    '${z}
#endif

rm jsub_new

end
