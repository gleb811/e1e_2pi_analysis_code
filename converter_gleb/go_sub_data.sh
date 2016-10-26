#!/bin/tcsh -f





@ x = 0
@ y = 1
while ($x <= 1843)
echo '#\!/bin/tcsh -f' > run_data$y.sh
echo ' ' >> run_data$y.sh
@ z = $x + 1
sed -e "1,$z d" gleb_file_list |  sed  -e "10q" | sed -e 's$/mss$jget /mss$g' | sed -e "s/.root/.root ./2">> run_data$y.sh
echo ' ' >> run_data$y.sh
echo 'echo "1" > inp' >> run_data$y.sh
echo 'echo "2.039" >> inp' >> run_data$y.sh
echo 'echo "10" >> inp' >> run_data$y.sh
sed -e "1,$z d" gleb_file_list |  sed  -e "10q" | sed -e 's$/mss/home/gleb/e1e/cooked_root/h10/$$g' | sed -e 's$clas$echo "clas$' | sed -e 's$root$root" >> inp$'>> run_data$y.sh
echo 'echo "out.root" >> inp' >> run_data$y.sh
echo './h10tot21_data<inp' >> run_data$y.sh
chmod +x run_data$y.sh
@ y += 1

@ x += 10
end

echo $y

@ s = $y - 2
echo $s

foreach k (`seq  1 $s`)

echo "PROJECT: e1e" > jsub_new
echo "TRACK: simulation" >>jsub_new
echo "JOBNAME: run${k}" >>jsub_new
echo "MAIL: skorodum@jlab.org" >>jsub_new
echo "SINGLE_JOB: TRUE" >>jsub_new
echo "MEMORY: 3500 MB" >>jsub_new
echo "DISK_SPACE: 15 GB" >>jsub_new
echo "OS:  centos62" >>jsub_new
echo "INPUT_FILES:/volatile/clas/clase1-6/skorodum/converter_gleb/h10tot21_data" >>jsub_new
echo "COMMAND: /volatile/clas/clase1-6/skorodum/converter_gleb/run_data${k}.sh" >>jsub_new
echo "OUTPUT_DATA: out.root" >>jsub_new
echo "OUTPUT_TEMPLATE: /mss/home/skorodum/e1e/GLEB_DATA_E1E_CONVERTED/t21_data${k}.root" >>jsub_new



if (!(-e /mss/home/skorodum/e1e/GLEB_DATA_E1E_CONVERTED/t21_data${k}.root)) then
echo $k
/site/bin/jsub jsub_new
endif

rm jsub_new

end
