#!/bin/tcsh -f

#set j=0
#foreach i (`seq 1 1000`)

#if (!(-e /cache/mss/home/gleb/e1e/sim2014/recsis_sim${i}.bos)) then
#if (!(-e /cache/mss/home/gleb/e1e/sim2014/ceb/out_ceb${i}.hbook)) then

#echo $i
#sed -e "s/test1/sim${i}/g" jsub_test > jsub${i}

#jsub jsub${i}

#rm jsub${i}

#@ j++

#endif

#end

#echo $j


foreach i (`seq 1 50000`)

@ j = ($i  - 1) / 5

if (-e run${j}.sh) then
rm run${j}.sh
endif
touch run${j}.sh
chmod +x run${j}.sh
echo '#\!/bin/tcsh -f' >> run${j}.sh
echo 'source /group/clas/builds/environment.csh' >> run${j}.sh

#echo "jget /mss/home/skorodum/e1e/simulation/sim_w_145_155_gleb_Sep2015/nt10/nt10_${i}.root ." >> run${j}.sh


if (($i % 5) == 0) then

@ s = $i - 4

#foreach k (`seq  $s $i`)
#echo "./nt10maker_mctk_new -t3 -ont10_${k}.hbook /cache/mss/home/gleb/e1e/sim_2014_w_1675_1825/recsis/recsis${k}.bos" >> run${j}.sh
#echo "h2root nt10_${k}.hbook" >> run${j}.sh
#end

echo 'echo "2" > inp' >> run${j}.sh
echo 'echo "2.039" >> inp' >> run${j}.sh
echo 'echo "5" >> inp' >> run${j}.sh

foreach k (`seq  $s $i`)
echo "echo 'nt10_${k}.root' >> inp" >> run${j}.sh
end

echo 'echo "out.root" >> inp' >> run${j}.sh
echo './h10tot21_sim<inp' >> run${j}.sh

#foreach k (`seq  $s $i`)
#echo "if (!(-e /mss/home/gleb/e1e/sim_2015_w_15_17/nt10/nt10_${k}.root)) then" >> run${j}.sh
#echo "jput nt10_${k}.root /mss/home/gleb/e1e/sim_2015_w_15_17/nt10/" >> run${j}.sh
#echo "endif" >> run${j}.sh
#end
@ i1 = $i - 1
@ i2 = $i - 2
@ i3 = $i - 3
@ i4 = $i - 4


echo "PROJECT: e1e" > jsub_new
echo "TRACK: simulation" >>jsub_new
echo "JOBNAME: run_sim_hiw_w${j}" >>jsub_new
echo "MAIL: skorodum@jlab.org" >>jsub_new
echo "SINGLE_JOB: TRUE" >>jsub_new
echo "MEMORY: 3000 MB" >>jsub_new
echo "DISK_SPACE: 15 GB" >>jsub_new
echo "OS:  centos65" >>jsub_new
#echo "INPUT_FILES:/volatile/clas/clase1-6/skorodum/converter_gleb/h10tot21_sim /volatile/clas/clase1-6/skorodum/converter_gleb/nt10maker_mctk_new" >>jsub_new
echo "INPUT_FILES:/volatile/clas/clase1-6/skorodum/converter_gleb/h10tot21_sim /mss/home/skorodum/e1e/simulation/sim_w_145_155_gleb_Sep2015/nt10/nt10_${i1}.root /mss/home/skorodum/e1e/simulation/sim_w_145_155_gleb_Sep2015/nt10/nt10_${i2}.root /mss/home/skorodum/e1e/simulation/sim_w_145_155_gleb_Sep2015/nt10/nt10_${i3}.root /mss/home/skorodum/e1e/simulation/sim_w_145_155_gleb_Sep2015/nt10/nt10_${i4}.root /mss/home/skorodum/e1e/simulation/sim_w_145_155_gleb_Sep2015/nt10/nt10_${i}.root" >>jsub_new
echo "COMMAND: /volatile/clas/clase1-6/skorodum/converter_gleb/run${j}.sh" >>jsub_new
echo "OUTPUT_DATA: out.root" >>jsub_new
echo "OUTPUT_TEMPLATE: /mss/clas/e1e/production/simulation_2pi/sim_conv_w_145_155_gleb/t21_${j}.root" >>jsub_new





if (!(-e /mss/clas/e1e/production/simulation_2pi/sim_conv_w_145_155_gleb/t21_${j}.root)) then
/site/bin/jsub jsub_new
endif

rm jsub_new





endif


end
