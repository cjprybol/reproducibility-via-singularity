# reproducibility-via-singularity
Repositories like Github, Bitbucket, and Gitlab have made writing and distributing research code a possibility for everyone. This advancement in democratizing code has spread much more quickly among researchers than the partner advancment required for true reproducability, which is a ability to distribute the computing environment necessary to faithfully execute the instructions. Given the diversity of projects that are conducted on shared HPC environments at many universities and laboratories, it's an impossible task for the managers of these clusters to manage and install the software necessary to run them all. The dependency libraries of software necessary to simulate models for the engineering teams may conflict with the dependencies of the sequencing analysis pipelines run by life scientists, and these both may conflict with the software needed by the computer science teams, yet they all work on the university or laboratory cluster. Package managers like Home/Linuxbrew and Anaconda have made it far easier for researchers to access software, but as most package managers tend to prioritize current releases over past releases, getting the proper configuration to reproduce a code execution remains a non-trivial practice. Singularity is a runtime container developed specifically for the shared-computing clusters of universities and laboratories that enables users to share completely mobile computing environments as easily as they share documents and image files. Singularity containers are stripped-down linux operating systems that can be loaded with software and configured exactly as desired. As long as another researcher as Singularity installed, they can run your code against the same version of software in the same computing environment as you ran the code against. I herein present an example of how to use Singularity to build an environment. I also present sample code and data such that you can reproduce a simple analysis, exactly as I performed it. If you find the exercise of reproducing my analysis from just a few lines of code at the command line, then I hope you consider using Singularity in your own research to make your work reproducable to others.

# TL;DR

run this code on a linux-based operating system with singularity and git installed and in the `$PATH`
```bash
git clone {this repo}
cd {this repo}
wget {image}
./reproduce_project.sh
```

# What is Singularity?

Singularity is a feature complete linux operating system, along with any files and software you want, wrapped into a single file that can be easily moved from computer to computer and reliably executed and interacted with via the Singularity runtime engine. Because it only has the aspects of the operating system required to execute code via the command line, they are leaner than full-fledged virtual machines, that can also replicate things like a disk drive and graphical user interfaces.

# What is Singularity not?

Singularity requires some understanding of how operating systems work in order to effectively use it and understand how it works. Singularity will generate an operating system inside of the image file, also refered to as a "container", and Singularity will ensure that the same code will be executed the same way when executed against that container. However, just as you may experience installation issues on your local computer, you may experience installations inside of the container. This is not necessarily a problem with the container, and is more likely an issue that the software you are attempting to install may not have been tested and troubleshooted against your particular linux flavor.

# How do I interact with a container?

**example**

if you would normally run a command as:
```bash
CMD --flag1 --flag2 argument1 argument2 ...
```
to run that command from within the singularity container you can simply pre-prend the appropriate singularity command
```bash
singularity exec my_singularity_container.img CMD --flag1 --flag2 argument1 argument2 ...
```

**show me a real example**

Here, python3 is not installed on the host system, and cannot be run
```bash
vagrant@jessie:~$ python3
-bash: python3: command not found
```
but this computer has singularity installed, and a container that has python3 installed
```bash
vagrant@jessie:~$ singularity exec test.img python3
Python 3.5.2 |Anaconda 4.1.1 (64-bit)| (default, Jul  2 2016, 17:53:06)
[GCC 4.4.7 20120313 (Red Hat 4.4.7-1)] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>>
```

# How do I build a singularity image to use with my own project?

Because containers interact so closely with the host system, installing singularity and building singularity images require root and/or sudo level permissions. While many users of shared computing clusters are unlikely to have these permission levels, they are likely to have those permissions on whatever local computer (laptop/desktop) they use to ssh into the shared computing cluster. If you happen to have linux-based OS with sudo permissions, go ahead and jump to the appropriately named **start here if linux user with root/admin/sudo priviliges**. If you have a Windows OS or an Apple OS, or even a linux OS without root/admin/sudo priviliges, you'll need to get a linux virtual machine running where you do have root permissions. If you're already familiar with the process of setting up virtual machines, go ahead with whatever method you know. If installing a virtual machine is a new thing for you, I recommend checking out [VirtualBox](https://www.virtualbox.org/). If you don't have root access on your local machine, ask a system admin for assistance.

**starting from a mac with homebrew installed**
```bash
brew tap caskroom/cask
brew cask install vagrant
# brew cask install virtualbox
```
Quick aside here: this implementation doesn't work with the new VirtualBox 5.1, and thus I've had to bail on the current (as of July 16, 2016) virtualbox cask install. Go to [VirtualBox's old build page](https://www.virtualbox.org/wiki/Download_Old_Builds_5_0) and install the appropriate version for your system
```bash
mkdir singularity-vm && cd singularity-vm
# generate a virtual machine
vagrant init ubuntu/trusty64; vagrant up --provider virtualbox && vagrant halt
```
Now we've initalized a lightweight virtual machine via a combination of vagrant and virtualbox. The virtual machine has some preset defaults that, in my experience, I've had to adjust. First and foremost the initial memory allocation. I'm unsure if the default is 1Gb RAM or less, but I only found success when bumping the memory allocation up to 4Gb.

```bash
vi Vagrantfile
```

Find the code block that looks like this
```
# config.vm.provider "virtualbox" do |vb|
#   # Display the VirtualBox GUI when booting the machine
#   vb.gui = true
#
#   # Customize the amount of memory on the VM:
#   vb.memory = "1024"
# end
```

and uncomment the code block and the memory line, setting the memory to 4Gb
```
config.vm.provider "virtualbox" do |vb|
#   # Display the VirtualBox GUI when booting the machine
#   vb.gui = true
#
#   # Customize the amount of memory on the VM:
  vb.memory = "4096"
end
```

now, restart the virtual machine and ssh into it
```bash
vagrant up && vagrant ssh
```
You are now a linux user with root/admin/sudo priviliges


**start here if linux user with root/admin/sudo priviliges**

Here we will install Singularity starting from an Ubuntu 14.04 LTS "Trusty" 64-bit base installation
```bash
sudo apt-get install -y build-essential git vim autoconf libtool curl debootstrap
git clone https://github.com/gmkurtzer/singularity.git && cd singularity && ./autogen.sh && ./configure --prefix=/usr/local && make && sudo make install
```

now we want to create a singularity container, load it with an operating system, and install our reproducable computing environment onto it.

first, we need to allocate the file. Here `--size` represents the maximum size that the container is allowed to take on.
```bash
sudo singularity create --size 24000 test.img
```
now we will preload a Ubuntu 14.04 LTS "Trusty" 64-bit base install
```bash
sudo singularity bootstrap test.img $HOME/singularity/examples/ubuntu.def
```
and now, we have a container! but it doesn't have much loaded onto it. let's enter into a bash shell that is running inside of the container. by default, containers have the ability to read and write data to the filesystem of the host, yet containers themselves are immutable. here, we will launch our shell session with the `--contain` flag, which blocks our ability to interact with the hosting computer, and the `--writable` flag such that we can modify the contents of the container
```bash
sudo singularity shell -w --contain test.img
```

Now you should be inside of the image. Feel free to jump around the file system to learn your way around.


As mentioned in the introduction, Linuxbrew and Anaconda are two great package managers that greatly simplify the process of installing and managing software. The default software available via apt-get is often out of date compared to those available via either Linuxbrew and Anaconda, and additionally, both of these package managers have recipes to install a wide-variety of software useful for my particular domain.

By default, when entering a container, you enter into it and stay in the same directory on the host where you started from. Run `pwd` and notice the `(unreachable)` prepension to the path. Let's change to the user directory of our current user inside of the container, which is `root`
```bash
cd /root
```

First, I will add two directories that I founded I need to add to actually get my singularity container to work on a remote system. On the cluster I work on, all of the labs have environments that are setup for their group as so
```
/scratch/PI/{PI_name}/"All lab data and projects go here"
/share/PI/{PI_name}/"All lab software and executables, and some reference files, go here"
```
The base directories for these paths, `/scratch` & `/share` are not standard linux directories, and are not found on the container. I had to add them to the container in order to get the container to resolve paths when I tried to use the container on the cluster
```bash
mkdir /scratch /share
```

Now we install required system dependencies for other software. Software required for your case may be different
```bash
apt-get update && apt-get install -y build-essential curl wget git python-setuptools ruby nettle-dev
```
and now we install linuxbrew, and open it for editing by all users
```bash
rm -r /usr/local && git clone https://github.com/Linuxbrew/brew.git /usr/local && chmod -R 777 /usr/local
```
and install my other favorite package manager and it's full python3 scientific computing environment, Anaconda
```bash
mkdir /Software && cd /Software
# anaconda and linuxbrew expect non-root usage, so make temporary user
useradd -m user && su user
wget http://repo.continuum.io/archive/Anaconda3-4.1.0-Linux-x86_64.sh && \
bash Anaconda3-4.1.0-Linux-x86_64.sh -b -p /Software/anaconda3 && \
rm Anaconda3-4.1.0-Linux-x86_64.sh
```
and make these new package managers available via our executable path
```bash
PATH="/Software/anaconda3/bin:/usr/local/sbin:/usr/local/bin:/bin:/sbin:/usr/sbin:/usr/bin"
```
now our recipe packed package maangers are available for us to install the latest and greatest bioinformatics software tools. If you have never reviewed their offerings, I recommend you check out their offerings here -> [homebrew-science](https://github.com/Homebrew/homebrew-science) [bioconda channel of anaconda](https://github.com/bioconda/bioconda-recipes/tree/master/recipes)
```bash
cd /home/user
brew install --force-bottle openssl open-mpi
brew install curl automake cmake curl git libtool parallel pigz wget
brew tap homebrew/science
brew install abyss art bamtools bcftools beagle bedtools bowtie bowtie2 blat bwa exonerate \
fastq-tools fastqc gmap-gsnap hmmer2 htslib jellyfish kallisto last lighter novoalign openblas picard-tools \
plink r samtools snap-aligner snpeff soapdenovo tophat trimmomatic varscan vcflib vcfanno vcftools velvet

conda config --add channels r
conda config --add channels bioconda
conda install -y pyaml pybedtools pyfasta pysam python-igraph pyvcf theano
conda install --channel https://conda.anaconda.org/conda-forge tensorflow
pip install keras
conda install -y -c r r
conda install -y --channel bioconda cufflinks cutadapt freebayes rsem rtg-tools sailfish salmon sambamba star plink2 trinity
```
now, I'll install Julia from source via git and github, as an example of how to manually install software not available via the package managers. Anything you need that isn't covered, go ahead and add it in here.
```bash
cd /Software
git clone git://github.com/JuliaLang/julia.git
cd julia
git checkout release-0.4
touch Make.user
echo "USE_SYSTEM_GMP=1" >> Make.user
echo "USE_SYSTEM_MPFR=1" >> Make.user
make
ln -s /Software/julia/julia /usr/local/bin
```
And now that we've got our system fully loaded with the software we want, let's make our extensions to `$PATH` to include the software installed by linuxbrew and anaconda permanent. You can do this too by modifying the `/environment` file in your container (**version >= 2.1**). I've got one saved on github in a gist that has the path updates needed for this container.
```bash
exit # back to root
cd / && rm /environment && wget --no-check-certificate https://gist.githubusercontent.com/cjprybol/e3baaabf9b95e65e765b9231d1594325/raw/9d8391b29fac7d0ed7b84442e1b3ebe2d3df3a36/environment
```
Now that linuxbrew is all set with software, let's remove the global ability to modify it's contents
```bash
chmod -R 775 /usr/local
```
exit the container to the host linux
```bash
exit
```
another tool I use that cannot be installed via either package manager due to licensing restrictions is gatk. load the compressed download of gatk into the folder where you will launch the container shell. then we will jump back into the container and link it to our anaconda path, and hop back out of the container.
```bash
sudo singularity shell -w test.img
gatk-register /home/vagrant/GenomeAnalysisTK-3.6.tar.bz2
exit
```
you're all done, you've built a great base-image for computational genomics! adjust these installation steps and software installed to your needs. you may have noticed that I installed lots of bioinformatics command line tools, and all of the packages I want for python, but not for R or Julia. Why? I found that because R and Julia both cooperate very nicely with installing their packages outside of the container, you can ship around a base language installation and have the first line of your analysis code download the necessary libraries. Anaconda has a method to do this with the python packages, however it did not work when I tried it.
```bash
singularity exec test.img R
install.packages( c("dplyr", "tidyr", "stringr", "lubridate", "ggplot2 ", "Hmisc", "caret", "randomForest", "survival", "parallel", "shiny", "glmnet", "datatable", "devtools", "Rcpp", "reshape2", "colorspace", "RColorBrewer", "plyr"))
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("affy", "affyio", "annotate", "biobase", "biocgenerics", "biocinstaller", "biocparallel", "biocstyle", "biomart", "biostrings", "deseq", "deseq2", "edger", "iranges", "limma", "preprocesscore", "rhdf5", "rsamtools", "s4vectors", "variantannotation", "pheatmap", "ChIPpeakAnno"))
quit()
```

```bash
singularity exec test.img julia
pkgs = [ "DataFrames", "FreqTables", "Distributions", "GLM", "HypothesisTests", "Nettle", "IJulia", "RCall", "NormalizeQuantiles", "Plots", "PyPlot", "GR", "BenchmarkTools", "MLBase", "NullableArrays" ]
map(Pkg.add, pkgs)
map(Pkg.build, pkgs)
using DataFrames, FreqTables, Distributions, GLM, HypothesisTests, Nettle, IJulia, RCall, NormalizeQuantiles, Plots, GR, BenchmarkTools, MLBase, NullableArrays, PyPlot
```

and, in progress, a method to return the version number of every piece of software installed on the container.
```bash
conda list | tail -n+3 | awk '{print $1, $2}' > temp.info
brew list --versions >> temp.info

julia --version >> temp.info
anything else I add via git

column -t temp.info | sort -u
```
