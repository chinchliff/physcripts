# setup gui
sudo apt-get install ubuntu-desktop

# get rid of ubuntu one
sudo apt-get purge ubuntuone-client* python-ubuntuone-*

# install text editors
sudo apt-get install emacs
sudo apt-get install geany

# install phylogenetics software
sudo apt-get install mafft
sudo apt-get install muscle

# install dev libraries
sudo apt-get install libnlopt-dev
sudo apt-get install liblapack-dev
sudo apt-get install libblas-dev 
sudo apt-get install libboost-dev
sudo apt-get install libatlas-dev
#sudo apt-get install libarmadillo-dev

# install dev tools
sudo apt-get install cmake

# setup git
sudo apt-get install git
sudo apt-get install meld
sudo apt-get install xclip

# setup scripts repo
cd ~/
git clone https://github.com/chinchliff/physcripts.git && mv physcripts scripts

# setup java
sudo apt-get purge openjdk*
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:webupd8team/java && sudo apt-get update && sudo apt-get install oracle-java7-installer

# setup chrome
cd ~/Downloads
wget https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb && sudo dpkg -i google-chrome-stable_current_amd64.deb

# setup eclipse
! wget http://ftp.osuosl.org/pub/eclipse/technology/epp/downloads/release/kepler/R/eclipse-java-kepler-R-linux-gtk-x86_64.tar.gz &&
tar -xvf eclipse-java-kepler-R-linux-gtk-x86_64.tar.gz &&
sudo mv eclipse /opt/
!
sudo apt-get install maven



# setup network share
sudo apt-get install nfs-common
sudo echo "blackrimlab:x:2071:cody: # added" >> /etc/group
sudo echo "eebsmithnfs1.value.storage.umich.edu:/eebsmithnfs1 /media/eebsmithnfs1   	nfs  rsize=32768,wsize=32768,timeo=14,intr,rw,users,exec 0 0 # added for remote storage at umich" >> /etc/fstab
sudo mkdir /media/eebsmithnfs1
sudo chown cody:blackrimlab /media/eebsmithnfs1
sudo chmod +rwx /media/eebsmithnfs1

# done
echo "Need to restart to mount network disk"



# this goes in the .bashrc file
alias ls="ls -lGA --color"
