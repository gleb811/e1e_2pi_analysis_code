#!/bin/tcsh -f


foreach i (`seq 0 100`)




echo "PROJECT: e1e" > jsub_new
echo "TRACK: simulation" >>jsub_new
echo "JOBNAME: run_w_145_155_${i}" >>jsub_new
echo "MAIL: gleb@jlab.org" >>jsub_new
echo "SINGLE_JOB: TRUE" >>jsub_new
echo "MEMORY: 8000 MB" >>jsub_new
echo "DISK_SPACE: 15 GB" >>jsub_new
echo "TIME: 4320" >>jsub_new
echo "OS:  centos65" >>jsub_new
echo "INPUT_FILES: /volatile/clas/clase1-6/gleb/2pi_analysis_e1e_new_bin/h10tot21 /volatile/clas/clase1-6/gleb/2pi_analysis_e1e_new_bin/new_ratio.root /volatile/clas/clase1-6/gleb/2pi_analysis_e1e_new_bin/phel_integr_fract.txt" >>jsub_new
echo "COMMAND: /volatile/clas/clase1-6/gleb/2pi_analysis_e1e_new_bin/execute.sh ${i}" >>jsub_new
echo "OUTPUT_DATA: out.root" >>jsub_new
echo "OUTPUT_TEMPLATE: /volatile/clas/clase1-6/gleb/2pi_analysis_e1e_new_bin/w_145_155_${i}.root" >>jsub_new





if (!(-e  /volatile/clas/clase1-6/gleb/2pi_analysis_e1e_new_bin/w_145_155_${i}.root)) then
/site/bin/jsub jsub_new
echo ${i}

endif

rm jsub_new



end
