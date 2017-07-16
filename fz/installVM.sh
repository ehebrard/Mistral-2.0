echo "This script should be placed in the directory ~/entry_data in the Virtual Machine provided by the challenge organizers"

echo "You should probably change the keyboard layout with \n dpkg-reconfigure keyboard-configuration (then restart) \n" 
echo "github, g++, python, and make should be installed before this step"

read -p "Press [Enter] key to start installation of Mistral-2.0"

git clone --recursive https://github.com/ehebrard/Mistral-2.0.git
cd Mistral-2.0
make clean
cd fz
make clean
make 

#To build the parallel version, uncomment the following: 
#mv fz/mistral-fzn ../
#make clean 
#cd fz 
#make clean 
#make parallel 
#mv ../../mistral-fzn ./

#Here  mistral-fzn and mistral-fzn_parallel are in fz/

cd /home/user/entry_data
cp Mistral-2.0/fz/mistral-fz fzn-exec
cp Mistral-2.0/fz/mistral-mzn ./exec

chmod 777 *exec
#echo "export PATH=\$DIR/entry_data/Mistral-2.0/fz:\$PATH" >> ../bin/challenge_env.sh
echo "PATH=/home/user/entry_data/Mistral-2.0/fz:$PATH" >> /home/user/.bashrc 

cp Mistral-2.0/fz/mznlib/* mzn-lib/
echo "Mistral-2.0 Installed. You need to logout then login to complete installation"
echo "Please edit ~/bin/challenge_env.sh to support only free (and parallel search?)”

#read -p "Press [Enter] key to start test"
#nano fzn-exec

#
#OPTIONAL :  to install X
#apt-get install xfce4
#apt-get install gedit
#apt-get install Midori

#Test if mznlib is ok! 

#cd

#echo "exec-free black-hole 12 black-hole.mzn 12.dzn
#exec-free black-hole 6 black-hole.mzn 6.dzn
#exec-free on-call-rostering 10s-50d oc-roster.mzn 10s-50d.dzn" > entry_data/listfile


#echo "Testing Global Constraints... [only cumulative]"

#cd
#cp /home/user/entry_data/Mistral-2.0/fz/test--mzn-lib.sh ./
#wget homepages.laas.fr/msiala/minizinc2014/test--mzn-lib.sh 
#chmod +x test--mzn-lib.sh 
#./test--mzn-lib.sh 
#rm test--mzn-lib.sh 
#ls 
#echo "End"
