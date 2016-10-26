#!/bin/tcsh -f


foreach k (`seq  145 145`)

echo "PROJECT: e1e" > jsub_new
echo "TRACK: simulation" >>jsub_new
echo "JOBNAME: run${k}" >>jsub_new
echo "MAIL: gleb@jlab.org" >>jsub_new
echo "SINGLE_JOB: TRUE" >>jsub_new
echo "MEMORY: 3500 MB" >>jsub_new
echo "DISK_SPACE: 15 GB" >>jsub_new
echo "OS:  centos65" >>jsub_new
echo "INPUT_FILES: /volatile/clas/clase1/gleb/test/2pi_analysis_e1e/converter_gleb/converter_gleb/h10tot21" >>jsub_new
echo "COMMAND: /volatile/clas/clase1/gleb/test/2pi_analysis_e1e/converter_gleb/converter_gleb/run_data${k}.sh" >>jsub_new
echo "OUTPUT_DATA: out.root" >>jsub_new
echo "OUTPUT_TEMPLATE: /mss/home/gleb/e1e/data_converted_2016/t21_data${k}.root" >>jsub_new



#if (!(-e /mss/home/skorodum/e1e/GLEB_DATA_E1E_CONVERTED/t21_data${k}.root)) then
#echo $k
/site/bin/jsub jsub_new
#endif

rm jsub_new

end
