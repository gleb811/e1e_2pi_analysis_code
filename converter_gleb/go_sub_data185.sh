#!/bin/tcsh -f




foreach k (`seq  185 185`)

echo "PROJECT: e1e" > jsub_data_185
echo "TRACK: simulation" >>jsub_data_185
echo "JOBNAME: run${k}" >>jsub_data_185
echo "MAIL: skorodum@jlab.org" >>jsub_data_185
echo "SINGLE_JOB: TRUE" >>jsub_data_185
echo "MEMORY: 3500 MB" >>jsub_data_185
echo "DISK_SPACE: 15 GB" >>jsub_data_185
echo "OS:  centos62" >>jsub_data_185
echo "INPUT_FILES:/volatile/clas/clase1-6/skorodum/converter_gleb/h10tot21_data" >>jsub_data_185
echo "COMMAND: /volatile/clas/clase1-6/skorodum/converter_gleb/run_data${k}.sh" >>jsub_data_185
echo "OUTPUT_DATA: out.root" >>jsub_data_185
echo "OUTPUT_TEMPLATE: /mss/home/skorodum/e1e/GLEB_DATA_E1E_CONVERTED/t21_data${k}.root" >>jsub_data_185



if (!(-e /mss/home/skorodum/e1e/GLEB_DATA_E1E_CONVERTED/t21_data${k}.root)) then
echo $k
/site/bin/jsub jsub_data_185
endif

rm jsub_data_185

end
