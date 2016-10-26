#!/bin/tcsh -f

foreach k (`seq  11 5 50000`)

@ z = $k / 5 + 1
echo $z

echo "PROJECT: e1e" > jsub_new
echo "TRACK: simulation" >>jsub_new
echo "JOBNAME: convert_high${k}" >>jsub_new
echo "MAIL: gleb@jlab.org" >>jsub_new
echo "SINGLE_JOB: TRUE" >>jsub_new
echo "MEMORY: 3000 MB" >>jsub_new
echo "DISK_SPACE: 15 GB" >>jsub_new
echo "OS:  centos65" >>jsub_new
@ k1 = $k + 1
@ k2 = $k + 2
@ k3 = $k + 3
@ k4 = $k + 4
echo "INPUT_FILES: /mss/home/skorodum/e1e/simulation/sim_w_165_185_gleb_Sep2015/nt10/nt10_${k}.root /mss/home/skorodum/e1e/simulation/sim_w_165_185_gleb_Sep2015/nt10/nt10_${k1}.root /mss/home/skorodum/e1e/simulation/sim_w_165_185_gleb_Sep2015/nt10/nt10_${k2}.root /mss/home/skorodum/e1e/simulation/sim_w_165_185_gleb_Sep2015/nt10/nt10_${k3}.root /mss/home/skorodum/e1e/simulation/sim_w_165_185_gleb_Sep2015/nt10/nt10_${k4}.root /volatile/clas/clase1/gleb/test/2pi_analysis_e1e/converter_gleb/converter_gleb/h10tot21" >>jsub_new
echo "COMMAND: /volatile/clas/clase1/gleb/test/2pi_analysis_e1e/converter_gleb/converter_gleb/execut_high.sh ${k}" >>jsub_new
echo "OUTPUT_DATA: out.root" >>jsub_new
echo "OUTPUT_TEMPLATE: /mss/clas/e1e/production/simulation_2pi/gleb_2016/high_w/high_w${z}.root" >>jsub_new

/site/bin/jsub jsub_new

rm jsub_new

end
