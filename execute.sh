#!/bin/tcsh -f
echo "1" > inp
echo "2.039" >> inp
echo "1" >> inp
echo "/cache/clas/e1e/production/simulation_2pi/gleb_2016/low_w_145_155/low_w_145_1551.root" >> inp
echo "1" >> inp
echo "/cache/clas/e1e/production/simulation_2pi/gleb_2016/low_w_145_155/low_w_145_1551.root" >> inp
echo "100" >> inp
foreach i (`seq 1 100`)
@ num =  $argv[1] * 100 + $i
echo "/cache/clas/e1e/production/simulation_2pi/gleb_2016/low_w_145_155/low_w_145_155${num}.root" >> inp
end
echo "out.root" >> inp
setenv dir `pwd`
cd /apps/root/PRO/root/
source bin/thisroot.csh
cd $dir
./h10tot21<inp
