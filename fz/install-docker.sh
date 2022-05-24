#!/bin/bash
#echo -e "This script is used to install the solver in the docker container used by the minizinc challenge. It should be placed in the directory ~/entry_data in the orignal minizinc chalenge container"
#read -p "Press [Enter] key to start installation of Mistral-2.0"
#apt-get update
#apt-get -y install nano
apt-get -y install wget
cd /entry_data
apt-get -y install software-properties-common
apt-add-repository universe
echo "deb http://dk.archive.ubuntu.com/ubuntu/ xenial main" >> /etc/apt/sources.list 
echo "deb http://dk.archive.ubuntu.com/ubuntu/ xenial universe" >> /etc/apt/sources.list
apt-get -y update
apt-get -y install g++-5
export CXX=g++-5
export CCC=g++-5
#apt-get -y install g++
apt-get -y install python pip
apt-get -y install libxml2-dev
apt-get -y install make 
apt-get -y install git
apt-get -y install bison 
apt-get -y install flex 
export CXX=g++-5
export CCC=g++-5
git clone https://github.com/ehebrard/Mistral-2.0.git
cd Mistral-2.0
wget https://boostorg.jfrog.io/artifactory/main/release/1.79.0/source/boost_1_79_0.tar.bz2
tar -xvf boost_1_79_0.tar.bz2
mv boost_1_79_0 boost
rm boost_1_79_0.tar.bz2
git clone https://github.com/xcsp3team/XCSP3-CPP-Parser.git 
git clone https://github.com/xcsp3team/XCSP3-Java-Tools.git
cd fz
cp template-minizinc-docker template.mk
#update template.mk to change the location of boost 
make
echo "test if it works \n \n \n \n \n \n \n \n \n \n \n "
./mistral-fz ../data/zinc/amaze3.fzn
./mistral-mzn ../data/zinc/amaze3.mzn ../data/zinc/2012-03-29.dzn 
#To build the parallel version
mv ./mistral-fzn ../..
cd ..
make clean 
cd fz 
make clean 
make parallel 
mv ../../mistral-fzn ./
#Here  mistral-fzn and mistral-fzn_parallel are in fz/
cd /entry_data
cp Mistral-2.0/fz/fzn-exec fzn-exec
cp Mistral-2.0/fz/mistral-fzn* ./
cp Mistral-2.0/fz/exec ./exec
#echo "PATH=/entry_data/Mistral-2.0/fz:$PATH" >> ~/.bashrc
#source ~/.bashrc
chmod 777 *
cp Mistral-2.0/fz/mznlib/* mzn-lib/
echo "Mistral-2.0 Installed."
cd 
#echo "Should I edit fzn-exec and exec to use /entry_data/mistral-fzn and /entry_data/fzn-exec ? "
#test :
./fzn-exec --print_sta Mistral-2.0/data/zinc/amaze3.fzn > fzn-test.txt
./exec --print_sta  Mistral-2.0/data/zinc/amaze3.mzn  Mistral-2.0/data/zinc/2012-03-29.dzn > mzn-test.txt
cat fzn-test.txt
cat mzn-test.txt
/minizinc/mzn-exec-free /minizinc/test.mzn /minizinc/2.dzn
/minizinc/mzn-exec-par -p 2 /minizinc/test.mzn /minizinc/2.dzn
/minizinc/mzn-exec-free Mistral-2.0/data/zinc/amaze3.fzn 
/minizinc/mzn-exec-free Mistral-2.0/data/zinc/amaze3.mzn  Mistral-2.0/data/zinc/2012-03-29.dzn
