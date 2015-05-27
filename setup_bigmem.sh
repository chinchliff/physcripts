RAXMLTYPE=SSE3

cd ~

# setup paths
mkdir $HOME/bin
mkdir -p $HOME/local/lib/python2.6/site-packages
echo 'export PATH=~/scripts:~/bin:$PATH' >> ~/.bashrc
echo 'export PYTHONPATH=~/scripts:~/local/lib/python2.6/site-packages:$PYTHONPATH' >> ~/.bashrc
source ~/.bashrc

# raxml
git clone http://github.com/stamatak/standard-raxml.git
cd standard-raxml
make -f Makefile.$RAXMLTYPE.PTHREADS.gcc
make -f Makefile.$RAXMLTYPE.gcc
mv raxmlHPC-$RAXMLTYPE-PTHREADS ~/bin
mv raxmlHPC-$RAXMLTYPE ~/bin
cd ../

# phyx
git clone git@github.com:FePhyFoFum/phyx.git
cd phyx/src
make -f Makefile.in pxbp
mv pxbp ~/bin

# physcripts
git clone http://github.com/chinchliff/physcripts.git
mv physcripts scripts

# indelible
scp hinchlif@smithbigmem1.eeb.lsa.umich.edu:~/bin/indelible ~/bin

# python libraries
easy_install --prefix=~/local -U dendropy
easy_install --prefix=~/local -U argparse

# simulation script
git clone git@github.com:chinchliff/taxonjackknife.git
cd taxonjackknife
git checkout py2.6
mkdir data_products

# temp dir
mkdir ~/temp
