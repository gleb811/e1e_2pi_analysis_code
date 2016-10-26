#!/bin/tcsh -f

set j=1821

foreach i (`seq 0 10000`)

if (!(-e /mss/home/skorodum/e1e/simulation/sim_conv_w_1275_1525_gleb/t21_${i}.root)) then

if (-e run${i}.sh) then
rm run${i}.sh
endif
touch run${i}.sh
chmod +x run${i}.sh
echo '#\!/bin/tcsh -f' >> run${i}.sh
echo 'source /group/clas/builds/environment.csh' >> run${i}.sh
echo 'echo "2" > inp' >> run${i}.sh
echo 'echo "2.039" >> inp' >> run${i}.sh
echo 'echo "5" >> inp' >> run${i}.sh

foreach k (`seq 1 5`)
setenv size 0
while ($size < 120000000)
@ j++
if ((-e /mss/home/skorodum/e1e/simulation/sim_w_1275_1525_gleb_Sep2015/nt10/nt10_${j}.root)) then
setenv size `cat  /mss/home/skorodum/e1e/simulation/sim_w_1275_1525_gleb_Sep2015/nt10/nt10_${j}.root | grep size | sed -e 's/size=//g'`
endif
end
echo "jget /mss/home/skorodum/e1e/simulation/sim_w_1275_1525_gleb_Sep2015/nt10/nt10_${j}.root ." >> run${i}.sh
echo "echo 'nt10_${j}.root' >> inp" >> run${i}.sh
end


echo 'echo "out.root" >> inp' >> run${i}.sh
echo './h10tot21_sim<inp' >> run${i}.sh

echo "PROJECT: e1e" > jsub_new
echo "TRACK: simulation" >>jsub_new
echo "JOBNAME: run_sim_low_w${i}" >>jsub_new
echo "MAIL: skorodum@jlab.org" >>jsub_new
echo "SINGLE_JOB: TRUE" >>jsub_new
echo "MEMORY: 3500 MB" >>jsub_new
echo "DISK_SPACE: 15 GB" >>jsub_new
echo "TIME: 4320" >>jsub_new
echo "OS:  centos65" >>jsub_new
#echo "INPUT_FILES:/volatile/clas/clase1-6/skorodum/converter_gleb/h10tot21_sim /volatile/clas/clase1-6/skorodum/converter_gleb/nt10maker_mctk_new" >>jsub_new
echo "INPUT_FILES:/volatile/clas/clase1-6/skorodum/converter_gleb/h10tot21_sim" >>jsub_new
echo "COMMAND: /volatile/clas/clase1-6/skorodum/converter_gleb/run${i}.sh" >>jsub_new
echo "OUTPUT_DATA: out.root" >>jsub_new
echo "OUTPUT_TEMPLATE: /mss/home/skorodum/e1e/simulation/sim_conv_w_1275_1525_gleb/t21_${i}.root" >>jsub_new







/site/bin/jsub jsub_new


rm jsub_new






endif

end
