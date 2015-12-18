#!/bin/sh

source ~/env32.sh
cd ~/Compact/utility/batch

list=`echo $@ | tr " " ","`
echo $list
sed -i "s/vector<int> incSq = {.*};/vector<int> incSq = {$list};/" combine_batch.cpp

cd ../build
make -j60

dirName=`date +"%s"`
mkdir $dirName
mv combine $dirName
cd $dirName


echo "mcfiles= /afs/cern.ch/user/n/nlurkin/work/L2_L3/dflt/mu.root /afs/cern.ch/user/n/nlurkin/work/L2_L3/dflt/pi.root
mcIndex=0 1
mcout=mu.root pi.root 
mccolors=330 440 550 630 740 850
#mccolors=330 330 330 330 330 330
brs=3.353E-2 2.066E-1
#brs=2.066E-1 2.066E-1 2.066E-1 2.066E-1 2.066E-1 2.066E-1
mclegends=K^{+}->#pi^{0}_{d}#mu^{+}#nu K^{+}->#pi^{+}#pi^{0}_{d}
#mclegends=K^{+}->#pi^{+}#pi^{0}_{d} K^{+}->#pi^{+}#pi^{0}_{d} K^{+}->#pi^{+}#pi^{0}_{d} K^{+}->#pi^{+}#pi^{0}_{d} K^{+}->#pi^{+}#pi^{0}_{d} K^{+}->#pi^{+}#pi^{0}_{d}

datafiles= /afs/cern.ch/user/n/nlurkin/work/L2_L3/dflt/data.root
dataout=data.root
datacolors=360
datalegends=Data

scanid=1" > listFile.spec.comb

bsub -q 8nh "source ~/env32.sh;
cd ~/Compact/utility/build/$dirName
./combine listFile.spec.comb
./combine listFile.spec.comb -b"

