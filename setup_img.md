```
wget https://raw.githubusercontent.com/cjprybol/reproducibility-via-singularity/master/ubuntu.def
sudo singularity create --size 500 demo.img
sudo singularity bootstrap demo.img ubuntu.def
sudo singularity shell -w --contain demo.img
mkdir /scratch /share
apt-get update && apt-get install -y build-essential curl wget && apt-get clean
wget https://github.com/pachterlab/kallisto/releases/download/v0.43.0/kallisto_linux-v0.43.0.tar.gz
tar -zxvf kallisto_linux-v0.43.0.tar.gz
rm kallisto_linux-v0.43.0.tar.gz
ln -s /kallisto_linux-v0.43.0/kallisto /usr/local/bin
exit
```
