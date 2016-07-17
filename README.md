# reproducibility-via-singularity
This repository provides a demo for how to use git and [singularity](https://github.com/gmkurtzer/singularity) to make your research more portable and reproducible.

# TL;DR

run this code on a linux-based operating system with singularity and git installed and in the `$PATH`
```bash
git clone {this repo}
cd {this repo}
wget {image}
./reproduce_project.sh
```

# Introduction

Reproducibility is core value of all scientific work, but getting it right is difficult. Using git, mercurial, SVN, or other kinds of version control software is a great start. Not only do these tools encourage useful practices such as versioning, development branches when exploring new features, and documentation via commit logs, they also make collaborative coding easier because multiple people can work on their own instances of the code independantly, and merge changes when they are finished. Web services such as github, gitlab, and bitbucket provide great tools for hosting your code repositories, which not only means you have a remote backup in case you lose your local copy, but these services also provide a very easy method by which you can share your code with collaborators or other researchers who would like to reproduce your work. All it takes is a `git clone` or, if they don't have git installed, they can download a compressed .zip file. However, even if someone has all of your code, there are many reasons they may not be able to reproduce your results. Some reasons are straightforward. A common issue for many researchers in academia or industry is that they do not have the admin-level permissions needed to install the software on their system. Another common issue is that the original author and the end-user may be running different versions of the same software, or the compute cluster the end user has access to many be running a different distribution of linux that requires a different set of installation instructions. One way around this is to run your analysis within a virtual machine (VM) that has all of the necessary software installed. This VM can then be distributed. VM's are great, but they are isolated from the machine they are run on, have a great computational cost and software overhead, and lead to much larger files that need to be managed and maintained. A new computing paradigm coined "containers" are essentially lighweight VMs that strip out much of the software bloat and computional overhead of full VMs, while also allowing a more natural flow between the isolated run-time environment of the container and the "host" system that the container is run on. Several instances of containers have started to gain traction, with the most widely-known implementation being [Docker](https://www.docker.com/). While Docker is a wonderful environment and has a strong user base with enterprise support (e.g. runs on AWS and Google Cloud), there are security issues that have prevented Docker from being the de-facto container adopted by researchers, who are more likely to have access to a shared computing cluster than they are to have access to a private cloud computing environment. [singularity](https://github.com/gmkurtzer/singularity) was developed with this audience in mind, and provides many features that system admins who maintain these computing resources require, such as encapsulation of the environment, no user contextual changes or root escalation allowed, no root owned daemon processes. Together, this means that you can ship a fully-functional container with all of your required software, exactly as you used it, and anyone else with singularity installed can run your code using that software. This project will serve as an example of how to setup a singularity container with several pieces of software commonly used by life-science researchers. I will demonstrate how to modify your code to utilize your 'container-ized' computing environment, and I will provide a set of scripts and an image so you can try out an example!

# As an end user looking to reproduce an analysis, what do I need to do?

The advantages of using singularity really shine for the end user. Once you have downloaded a singularity image, you can run the software installed on the container as if it were installed on the host computer. A nice way to conceptualize a container is as extension of your computer, rather than some isolated item on your computer. It's as if you could plug in a USB drive that came preloaded with all of the software you needed, and you can call the software programs and have it work immediately without needing to worry about installing each program 1 by 1 or worrying about dependency issues (software that your software requires to run; like how java-based applications require you to already have java installed). That's because unlike a USB drive or other drive that depends on the operating system of your computer, singularity containers have all of the operating system components they need to run so you can avoid all of the headaches of installing everything one by one. You also don't have to worry about accidentally breaking anything in the container, because the critical components inside of the container are protected by read-write permissions the same way that the core operating system components on your host computer are. For example, a general user on a shared computing cluster at a university is unlikely to have permissions to edit the `/bin` directory. They might not even be able to look to see what is inside. Those same restrictions and rules apply to whatever is inside of the container. Anything you run that exists inside of the container has the same limits on the host computer as the user does. If you are not able to edit a folder, then neither is the container. But if you have a folder in which you would like to run an analysis, and you generate files in the process of running that analysis, whatever you execute inside of the container can create, remove, and edit those files with the same permissions as yourself.

example
```bash
# if you would normally run a command as
CMD --flag1 --flag2 argument1 argument2 ...
# to run that command from within the singularity container
# you can simply pre-prend the appropriate singularity command
singularity exec my_singularity_container.img CMD --flag1 --flag2 argument1 argument2 ...
```

show me a real example
```
# python3 is not installed on the local system
vagrant@jessie:~$ python3
-bash: python3: command not found
# but it is installed inside of the container
vagrant@jessie:~$ singularity exec test.img python3
Python 3.5.2 |Anaconda 4.1.1 (64-bit)| (default, Jul  2 2016, 17:53:06)
[GCC 4.4.7 20120313 (Red Hat 4.4.7-1)] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>>
```

# OK, I'm interested. How do I build a singularity image to use with my own project?

Because containers interact so closely with the host system, installing singularity and building singularity images require root and/or sudo level permissions. While many users of shared computing clusters are unlikely to have these permission levels, they are likely to have those permissions on whatever local computer (laptop/desktop) they use to ssh into the shared computing cluster. If you happen to have linux-based OS with sudo permissions, go ahead and jump to the next section. If you have a Windows OS or an Apple OS, you'll need to get a linux distribution running. If you're already familiar with the process of setting up virtual machines, go ahead with whatever method you know. If installing a virtual machine is a new thing for you, I recommend checking out [VirtualBox](https://www.virtualbox.org/). If you don't have root access on your local machine (windows mac or linux), check with the system admin to see if they'd be willing to set up a virtual machine for you.

How I did it: starting from a mac with homebrew installed
```bash
# get necessary software
brew tap caskroom/cask

# install virtualbox 5.0.24

brew cask install vagrant
# make a folder to host the vagrant VM information
mkdir singularity-vm
# move into it
cd singularity-vm
# generate a virtual machine
vagrant init ubuntu/trusty64; vagrant up --provider virtualbox
# edit vagrant file to have 4gb memory
# uncomment the following lines
# config.vm.provider "virtualbox" do |vb|
# vb.memory = "1024"
# end
# & change memory to = "4096"

vagrant halt && vagrant up && vagrant ssh
sudo apt-get install -y build-essential git vim autoconf libtool curl debootstrap
# clone the singularity software from github
git clone https://github.com/gmkurtzer/singularity.git && cd singularity && ./autogen.sh && ./configure --prefix=/usr/local && make && sudo make install
```



# I've got a linux machine with singularity installed. How do I build an image

```bash
# This will create an image with a MAXIMUM size of 8Gb
# adjust to your needs
# if you run out of space, you can update the maximum via
# singularity expand --size {Mib to expand by} your_image.img
sudo singularity create -s 24000 test.img
# This will set up a basic ubuntu image inside of the container
sudo singularity bootstrap test.img $HOME/singularity/examples/ubuntu.def
# this will enter into a bash shell in the container
# -w means allow user write permissions on the container
# -w requires sudo privileges, and the container is not writable
# by default
# --contain means restrict the filesystem to only include
# the filesystem of the container. i.e. don't accidentally write to host
sudo singularity shell -w --contain test.img
```

Now you should be inside of the image. Feel free to jump around the file system to learn your way around.






```bash
cd /root
# singularity containers
mkdir /scratch /share
apt-get update && apt-get install -y build-essential curl git python-setuptools ruby nettle-dev
rm -r /usr/local && git clone https://github.com/Linuxbrew/brew.git /usr/local && chmod -R 777 /usr/local
# add user to install brew things with, because it complains if you install as root



useradd -m user && su user && cd /home/user
brew install --force-bottle openssl open-mpi
brew install curl automake cmake curl git libtool parallel pigz wget
brew tap homebrew/science
brew install abyss art bamtools bcftools beagle bedtools bowtie bowtie2 blat bwa exonerate \
fastq-tools fastqc gmap-gsnap hmmer2 htslib jellyfish kallisto last lighter novoalign openblas picard-tools \
plink r samtools snap-aligner snpeff soapdenovo tophat trimmomatic varscan vcflib vcfanno vcftools velvet
exit
cd / && rm /environment && wget --no-check-certificate https://gist.githubusercontent.com/cjprybol/e3baaabf9b95e65e765b9231d1594325/raw/9d8391b29fac7d0ed7b84442e1b3ebe2d3df3a36/environment
chmod -R 775 /usr/local
# blat dependency mysql takes 30 minutes to build
# consider dropping


# back to root
mkdir /Software
cd /Software
git clone git://github.com/JuliaLang/julia.git
cd julia
git checkout release-0.4
touch Make.user
echo "USE_SYSTEM_GMP=1" >> Make.user
echo "USE_SYSTEM_MPFR=1" >> Make.user
make
ln -s /Software/julia/julia /usr/local/bin

# wget --no-check-certificate https://julialang.s3.amazonaws.com/bin/linux/x64/0.4/julia-0.4.6-linux-x86_64.tar.gz
# tar xvzf julia-0.4.6-linux-x86_64.tar.gz && rm julia-0.4.6-linux-x86_64.tar.gz && ln -s /Software/julia-*/bin/julia /usr/local/bin
# rmdir /home/vagrant/.julia

# in future releases, I'll be able to just add these to the $PATH
# rather than obliterating the system defaults
cd /Software
wget http://repo.continuum.io/archive/Anaconda3-4.1.0-Linux-x86_64.sh && \
bash Anaconda3-4.1.0-Linux-x86_64.sh -b -p /Software/anaconda3 && \
rm Anaconda3-4.1.0-Linux-x86_64.sh
PATH="/Software/anaconda3/bin:/usr/local/sbin:/usr/local/bin:/bin:/sbin:/usr/sbin:/usr/bin"

conda config --add channels r
conda config --add channels bioconda
conda install -y pyaml pybedtools pyfasta pysam python-igraph pyvcf theano
conda install --channel https://conda.anaconda.org/conda-forge tensorflow
pip install keras
conda install -y -c r r
conda install -y --channel bioconda cufflinks cutadapt freebayes rsem rtg-tools sailfish salmon sambamba star plink2 trinity

exit
exit?
sudo singularity shell -w test.img
gatk-register /home/vagrant/GenomeAnalysisTK-3.6.tar.bz2
exit
```
```bash
R
install.packages( c("dplyr", "tidyr", "stringr", "lubridate", "ggplot2 ", "Hmisc", "caret", "randomForest", "survival", "parallel", "shiny", "glmnet", "datatable", "devtools", "Rcpp", "reshape2", "colorspace", "RColorBrewer", "plyr"))
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("affy", "affyio", "annotate", "biobase", "biocgenerics", "biocinstaller", "biocparallel", "biocstyle", "biomart", "biostrings", "deseq", "deseq2", "edger", "iranges", "limma", "preprocesscore", "rhdf5", "rsamtools", "s4vectors", "variantannotation", "pheatmap", "ChIPpeakAnno"))
quit()
```

```bash
julia
# install julia packages outside of the container
# is there a way to install if not installed? require?
pkgs = [ "DataFrames", "FreqTables", "Distributions", "GLM", "HypothesisTests", "Nettle", "IJulia", "RCall", "NormalizeQuantiles", "Plots", "PyPlot", "GR", "BenchmarkTools", "MLBase", "NullableArrays" ]
map(Pkg.add, pkgs)
map(Pkg.build, pkgs)
using DataFrames, FreqTables, Distributions, GLM, HypothesisTests, Nettle, IJulia, RCall, NormalizeQuantiles, Plots, GR, BenchmarkTools, MLBase, NullableArrays, PyPlot
```

todo
poretools marginalign

```
conda list | tail -n+3 | awk '{print $1, $2}' > temp.info
brew list --versions >> temp.info

julia --version >> temp.info
anything else I add via git

column -t temp.info | sort -u
```
