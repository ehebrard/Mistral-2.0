echo -e "This script is used to install the solver in the docker container used by the minizinc challenge. It should be placed in the directory ~/entry_data in the orignal minizinc chalenge container\n"
read -p "Press [Enter] key to start installation of Mistral-2.0"


apt-get update
apt-get -y install wget
cd /entry_data

apt-get -y install nano
apt-get -y install software-properties-common
apt-add-repository universe
apt-get -y update
apt-get -y install g++-5
apt-get -y install python-pip
apt-get -y install wget
apt-get -y install libxml2-dev
apt-get -y install make 
apt-get -y install git
apt-get install bison 
apt-get install flex 

export CXX=g++-5
export CCC=g++-5


git clone https://github.com/ehebrard/Mistral-2.0.git
cd Mistral-2.0
git clone https://github.com/xcsp3team/XCSP3-CPP-Parser.git 
git clone https://github.com/xcsp3team/XCSP3-Java-Tools.git
#make clean
cd fz
read -p "I should edit template.mk and ../template.mk to use g++-5"
make clean
make

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
cp Mistral-2.0/fz/mistral-fz fzn-exec
cp Mistral-2.0/fz/mistral-fzn* ./
cp Mistral-2.0/fz/mistral-mzn ./exec
#export PATH=$PATH:/entry_data/Mistral-2.0/fz/
#echo "PATH=/entry_data/Mistral-2.0/fz:$PATH" >> ~/.bashrc
chmod 777 *
cp Mistral-2.0/fz/mznlib/* mzn-lib/
echo "Mistral-2.0 Installed."

echo "I should edit fzn-exec and exec to use /entry_data/mistral-fzn and /entry_data/fzn-exec "

#test :

./fzn-exec Mistral-2.0/data/zinc/amaze3.fzn 
./exec Mistral-2.0/data/zinc/amaze3.mzn  Mistral-2.0/data/zinc/2012-03-29.dzn 

/minizinc/mzn-exec-free /minizinc/test.mzn /minizinc/2.dzn
/minizinc/mzn-exec-par -p 2 /minizinc/test.mzn /minizinc/2.dzn
 
