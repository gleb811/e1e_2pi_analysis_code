#!/bin/tcsh -f

set arr = ( 1112 ) 

foreach k (`seq  340 5 345`)

@ z = $k / 5 + 1 - 68
echo "med_w$arr[$z].root"

echo "PROJECT: e1e" > jsub_new
echo "TRACK: simulation" >>jsub_new
echo "JOBNAME: convert_med${k}" >>jsub_new
echo "MAIL: gleb@jlab.org" >>jsub_new
echo "SINGLE_JOB: TRUE" >>jsub_new
echo "MEMORY: 3000 MB" >>jsub_new
echo "DISK_SPACE: 15 GB" >>jsub_new
echo "OS:  centos65" >>jsub_new
@ k1 = $k + 1
@ k2 = $k + 2
@ k3 = $k + 3
@ k4 = $k + 4
echo "INPUT_FILES: /mss/clas/e1e/production/simulation_2pi/gleb_2016/w_15_17_raw/nt10/nt10_${k}.root /mss/clas/e1e/production/simulation_2pi/gleb_2016/w_15_17_raw/nt10/nt10_${k1}.root /mss/clas/e1e/production/simulation_2pi/gleb_2016/w_15_17_raw/nt10/nt10_${k2}.root /mss/clas/e1e/production/simulation_2pi/gleb_2016/w_15_17_raw/nt10/nt10_${k3}.root /mss/clas/e1e/production/simulation_2pi/gleb_2016/w_15_17_raw/nt10/nt10_${k4}.root /volatile/clas/clase1-6/gleb/2pi_analysis_e1e/converter_gleb/converter_gleb/h10tot21" >>jsub_new
echo "COMMAND: /volatile/clas/clase1-6/gleb/2pi_analysis_e1e/converter_gleb/converter_gleb/execut_med.sh ${k}" >>jsub_new
echo "OUTPUT_DATA: out.root" >>jsub_new
echo "OUTPUT_TEMPLATE: /mss/clas/e1e/production/simulation_2pi/gleb_2016/w_15_17/med_w$arr[$z].root" >>jsub_new


#if (!(-e /mss/clas/e1e/production/simulation_2pi/gleb_2016/w_15_17/med_w${z}.root)) then
/site/bin/jsub jsub_new
echo ${z}
#endif

rm jsub_new

end
