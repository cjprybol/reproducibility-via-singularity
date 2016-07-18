# Make your research project more reproducable with [Singularity](http://singularity.lbl.gov)
Services like [GitHub](https://github.com/), [Bitbucket](https://bitbucket.org), and [GitLab](https://about.gitlab.com/) have democratized access to affordable (or free!) tools that reinforce reproducability in research, via the code repositories that these services offer. These services make it easier to backup, version control, collaborate on, and distribute code (or any other text-based file). These features make it easier for researchers to write and maintain high-quality code. These features also increase the chance that someone else will review the code for errors or bugs. These reviews can be done by some combination of reading and executing the code. Unfortunately, unless the code is run on the exact same computer, while logged in as the same user, getting the same code to run the same way can be a research project in and of itself.


# Why Singularity?
Universities and research laboratories often conduct their work on shared HPC clusters. These clusters are all set up on slightly different configurations of hardware and operating system software. Trying to recreate an environment to re-run code exactly as it was executed on another cluster, an environment where all of the C pointers align in the same way and the Java versions sync up, is beyond the abilities, free time, and account priviliges of most researchers who may want to try to reproduce your results. Singularity is an implementation of a "container", and an engine to run those containers, that allows researchers to isolate the environment needed to produce a result away from the resources and code that is already available. This means that your colleague at University X can run the analysis exactly the same way on their cluster as you are running it on your cluster at University Y, and all it requires is sharing a git repository and a container image.

# TL;DR

Run this code on a linux-based operating system with singularity installed and in the `$PATH`
```bash
git clone git@github.com:cjprybol/reproducibility-via-singularity.git
cd reproducibility-via-singularity
wget https://stanfordmedicine.box.com/shared/static/i108j365ewur9cudcbehmd5aij9df2gu.img -O demo.img
./reproduce_project.sh
```

# What does that code do?

I encourage you to read and execute the `reproduce_project.sh` file and see for yourself! But in summary, it will download 1 sample of paired-end RNA sequencing reads, and quantify transcript isoform abundances in the sample. Run `less Kallisto/abundance.tsv` to see the results!

This is a contrived example, designed to run quickly and efficiently to demonstrate a point that reproducable research can be, and should be, be as simple as running 4 code instructions. However, by extending the software libraries in the container, and extending the analytical complexity in the analysis scripts, these same 4 lines of code can be sufficient to reproduce entire thesis projects. 

**Software Requirements**
  - RAM ~ 8Gb
  - CPU = 1
  - available disk space: ~3 Gb

time to run > 10 minutes on a modern laptop with decent download speeds

Requiring little more effort for others than the wait required to perform the computation, this simplicity of reproducability lowers the barrier for others to investigate your work. Together, this encourages responsible research conduct and increases the rate of knowledge transfer.

# What is Singularity?

A Singularity container is a complete linux kernel capable of executing code and running software. Because the container is a singe file, it can be easily moved to another computer with Singularity and run in exactly the same way. The Singularity engine on each cluster handles the translation to the unique hardware and OS configurations of each cluster.

# What is Singularity not?

A magic, drag and drop, point and click solution. Singularity requires some understanding of how operating systems work in order to effectively use it and understand it. That being said, it's also very well designed, and should be intuitive for anyone who has set up a computing environment on their personal computer, a cluster, a VM, or a cloud service. Singularity will generate an operating system inside of the image file, or "container". Singularity will ensure that whatever you install in the container will execute the same on any computer with the Singularity installed. You simply load the container with software and configure it as you would set up any other basic linux installation.

Containers require sudo/root/admin priviliges to modify and edit. This design feature makes Singularity containers incompatible with software that is programmed to download files to, or edit files within, the same directory (or sub-directories) where the the software is installed, as those directories will be within the immutable container. This design strategy is used by many annotation libraries that download data from public databases. However, if these libraries allow the option to write to a directory on the host, such as the users `$HOME` directory, or a directory specified by the user, than this will not be a problem. Another approach would be to simply pre-download data and pre-generate files, although this may result in container files that are too large to effectively share.

# How do I interact with a container?

**Example**

If you would normally run a command as:
```bash
CMD --flag1 --flag2 argument1 argument2 ...
```

That command can be executed by the software inside of the container, instead of by the host software, by prepending the appropriate command
```bash
singularity exec my_singularity_container.img CMD --flag1 --flag2 argument1 argument2 ...
```

**Show me a real example**

Here, python3 is not installed on the host system, and cannot be run
```bash
vagrant@jessie:~$ python3
-bash: python3: command not found
```

But python3 is installed inside of a container present in the current directory
```bash
vagrant@jessie:~$ singularity exec test.img python3
Python 3.5.2 |Anaconda 4.1.1 (64-bit)| (default, Jul  2 2016, 17:53:06)
[GCC 4.4.7 20120313 (Red Hat 4.4.7-1)] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>>
```

# How do I build a singularity image to use with my own project?

Because containers interact so closely with the host system, installing Singularity, and building Singularity images requires root/sudo/admin level permissions. While many users of shared computing clusters are unlikely to have these permission level on the clusters where they intend to run the analysis, they may have those permissions on whatever local computer (laptop/desktop) they use to ssh into the shared computing cluster. If you happen to have linux-based OS with sudo permissions, go ahead and jump to the section **Start here if linux user with root/sudo/admin priviliges**. If you have a Windows OS or an Apple OS (with our without root/sudo/admin), or a linux OS without root/sudo/admin, you'll need to get a linux virtual machine running where you do have root/sudo/admin permissions. If you're already familiar with the process of setting up virtual machines, go ahead with whatever method you know. If installing a virtual machine is a new thing for you, I recommend checking out [VirtualBox](https://www.virtualbox.org/). If you don't have root access on your local machine, ask a system admin for assistance.

# Starting from a mac with [Homebrew](http://brew.sh/) installed
```bash
brew tap caskroom/cask
brew cask install vagrant
# brew cask install virtualbox
```
**Quick Fix**: These instructions don't work with the new VirtualBox 5.1 `brew cask` install, and thus I've had to bail on the current (as of July 16, 2016) cask. Go to [VirtualBox's old build page](https://www.virtualbox.org/wiki/Download_Old_Builds_5_0) and install the appropriate version for your system. Or try the `brew cask` method and let me know that it works again and I'll update these instructions.
```bash
mkdir singularity-vm && cd singularity-vm
# generate a virtual machine
vagrant init ubuntu/trusty64; vagrant up --provider virtualbox && vagrant halt
```
Now we have a lightweight virtual machine implemented with vagrant and virtualbox. The virtual machine has some preset defaults that were too restrictive for running this demo. The default memory allocation is was unable to execute the code, however I was able to run the demo after setting the memory allocation to 8Gb.

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

Uncomment the code block and the memory line, setting the memory to 4Gb
```
config.vm.provider "virtualbox" do |vb|
#   # Display the VirtualBox GUI when booting the machine
#   vb.gui = true
#
#   # Customize the amount of memory on the VM:
  vb.memory = "8192"
  vb.cpus = 2
end
```

Restart the virtual machine and ssh into it
```bash
vagrant up && vagrant ssh
```
You are now a linux user with root/sudo/admin priviliges


# Start here if linux user with root/sudo/admin priviliges

Here we will install Singularity starting from an Ubuntu 14.04 LTS "Trusty" 64-bit base installation
```bash
sudo apt-get install -y build-essential git vim autoconf libtool curl debootstrap
git clone https://github.com/gmkurtzer/singularity.git && cd singularity && ./autogen.sh && ./configure --prefix=/usr/local && make && sudo make install
```

We want to create a singularity container, load it with an operating system, and install our reproducable computing environment onto it.

First, we need to allocate the file. Here `--size` represents the maximum size (in Mb) that the container is allowed to take on.
```bash
sudo singularity create --size 24000 test.img
```

We will preload a Ubuntu 14.04 LTS "Trusty" 64-bit base install
```bash
sudo singularity bootstrap test.img $HOME/singularity/examples/ubuntu.def
```

We have a container! But it doesn't have much loaded onto it. Let's enter into a bash shell that is running inside of the container. By default, containers have the ability to read and write data to the filesystem of the host, yet containers themselves are immutable. Here, we will launch our shell session with the `--contain` flag, which blocks our ability to interact with the hosting computer, and the `--writable` flag such that we can modify the contents of the container
```bash
sudo singularity shell -w --contain test.img
```

Now you should be inside of the image. Feel free to jump around the file system to learn your way around.

[Linuxbrew](http://linuxbrew.sh/) and [Anaconda](https://www.continuum.io/downloads) are two great package managers that greatly simplify the process of installing and managing software. The default software available via apt-get is often out of date compared to those available via either Linuxbrew and Anaconda, and additionally, both of these package managers have recipes to install a wide-variety of software useful for researchers.

By default, when entering a container, you enter into it and stay in the same directory on the host where you started from. Run `pwd` and notice the `(unreachable)` prepension to the path. Let's change to the user directory of our current user inside of the container, which is `root`
```bash
cd /root
```

First, I will add two directories that I founded I need to add to actually get my singularity container to work on a remote system. On the cluster I work on, all of the labs have environments that are setup for their group as so
```
/scratch/PI/{PI_name}/"All lab data and projects go here"
/share/PI/{PI_name}/"All lab software and executables, and some reference files, go here"
```
The base directories for these paths, `/scratch` & `/share` are not standard linux directories, and are not found on the container. I had to add them to the container in order to get the container to resolve paths when I tried to use the container on the cluster. **If you experience a similar issue, let me know. If there are other common directories used on other HPC clusters, I'll consider adding them here**
```bash
mkdir /scratch /share
```

Here we install required system dependencies for other software. Software required for your case may be different
```bash
apt-get update && apt-get install -y build-essential curl wget git python-setuptools ruby nettle-dev
```

Install Linuxbrew, and open it for editing by all users (temporarily)
```bash
rm -r /usr/local && git clone https://github.com/Linuxbrew/brew.git /usr/local && chmod -R 777 /usr/local
```

Make a new directory for installing additional software
```bash
mkdir /Software && cd /Software
```

Linuxbrew expects non-root usage, so make temporary user. This will ensure any software installed is executable without root/sudo/admin permissions outside of the container as well.
```bash
useradd -m user && su user
```

Install Anaconda
```bash
wget http://repo.continuum.io/archive/Anaconda3-4.1.0-Linux-x86_64.sh && \
bash Anaconda3-4.1.0-Linux-x86_64.sh -b -p /Software/anaconda3 && \
rm Anaconda3-4.1.0-Linux-x86_64.sh
```

Add the new software to the `$PATH`
```bash
PATH="/Software/anaconda3/bin:/usr/local/sbin:/usr/local/bin:/bin:/sbin:/usr/sbin:/usr/bin"
```

Now let's use our package managers to quickly and easily install and configure software into the container. You can review their full offerings here -> [Homebrew-science](https://github.com/Homebrew/homebrew-science) & [Bioconda channel of anaconda](https://github.com/bioconda/bioconda-recipes/tree/master/recipes)

**Linuxbrew**
```bash
cd /home/user
brew install --force-bottle openssl open-mpi
brew install curl automake cmake git libtool parallel pigz wget
brew tap homebrew/science
brew install abyss art bamtools bcftools beagle bedtools bowtie bowtie2 blat bwa exonerate fastq-tools fastqc gmap-gsnap hmmer2 htslib jellyfish kallisto last lighter novoalign openblas picard-tools plink r samtools snap-aligner snpeff soapdenovo tophat trimmomatic varscan vcflib vcfanno vcftools velvet

```

**Anaconda**
```bash
conda config --add channels r
conda config --add channels bioconda
conda install -y pyaml pybedtools pyfasta pysam python-igraph pyvcf theano
conda install --channel https://conda.anaconda.org/conda-forge tensorflow
pip install keras
conda install -y -c r r
conda install -y --channel bioconda cufflinks cutadapt freebayes rsem rtg-tools sailfish salmon sambamba star plink2 trinity
conda clean -y --all
```

I'll install Julia from source via GitHub, as an example of how to manually install software not available via the package managers, as well as a plug for the language (which I recommend you try out!). Feel free to add your own custom recipes for installing software here, to configure a container that meets your needs.
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

We've got our system fully loaded with the software we want, but our `$PATH` update was only for this session. We'll need to make our `$PATH` updates permanent to make the software installed inside of the `/Software` directory available by name alone. Alternatively, you can also specify the full path when calling executables inside of the container. In Singularity version >= 2.1, you can update the `$PATH` by modifying the `/environment` file, which is loaded each time you interact with the container. I have a pre-prepared `/environment` file that I saved to a gist for easy access.
```bash
exit # back to root
cd / && rm /environment && wget --no-check-certificate https://gist.githubusercontent.com/cjprybol/e3baaabf9b95e65e765b9231d1594325/raw/9d8391b29fac7d0ed7b84442e1b3ebe2d3df3a36/environment
```

Now that linuxbrew is all set with software, let's remove the global ability to modify it's contents
```bash
chmod -R 775 /usr/local
```

Exit the container to the host linux
```bash
exit
```

Another tool I use that cannot be installed via either package manager due to licensing restrictions is [Genome Analysis Toolkit](https://www.broadinstitute.org/gatk/). I've downloaded the GATK installer to the same directory where my container is. By entering the container again without the `--contain` flag, we allow the container to interact with the host system, and copy GATK into Anaconda
```bash
sudo singularity shell -w test.img
gatk-register /home/vagrant/GenomeAnalysisTK-3.6.tar.bz2
exit
```

You're all done, you've built a great base-image for computational genomics! Adjust these installation steps to your needs.

You may have noticed that I installed lots of bioinformatics command line tools, and all of the packages I want for Python, but none for R or Julia. Why? and How can we install packages to R and Julia if we can't modify the contents of the container? By trial and error I discovered that R is happy to install packages **either** inside of the container **or** outside of the container, but it has trouble when the package library is split between both the container and the host. I think it's easier to just let R install packages to the user's `$HOME` directory, and this can be included in the source code that is shipped along-with the container. Julia, on the other hand, always wants to install libraries outside of the container, so the initial steps of installing Julia packages must be part of the source code. Anaconda provides a method to copy it's package library onto the host, however it did not work for me. Because Anaconda's package library is so complete, I was less worried about locking the packages inside of the container, and just decided to preload the packages and lock the library inside of the container.

**Load R packages**
```bash
singularity exec test.img R
install.packages( c("dplyr", "tidyr", "stringr", "lubridate", "ggplot2 ", "Hmisc", "caret", "randomForest", "survival", "parallel", "shiny", "glmnet", "datatable", "devtools", "Rcpp", "reshape2", "colorspace", "RColorBrewer", "plyr"))
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("affy", "affyio", "annotate", "biobase", "biocgenerics", "biocinstaller", "biocparallel", "biocstyle", "biomart", "biostrings", "deseq", "deseq2", "edger", "iranges", "limma", "preprocesscore", "rhdf5", "rsamtools", "s4vectors", "variantannotation", "pheatmap", "ChIPpeakAnno"))
quit(save="no")
```

**Load Julia packages**
```bash
singularity exec test.img julia
pkgs = [ "DataFrames", "FreqTables", "Distributions", "GLM", "HypothesisTests", "Nettle", "IJulia", "RCall", "NormalizeQuantiles", "Plots", "PyPlot", "GR", "BenchmarkTools", "MLBase", "NullableArrays" ]
map(Pkg.add, pkgs)
map(Pkg.build, pkgs)
using DataFrames, FreqTables, Distributions, GLM, HypothesisTests, Nettle, IJulia, RCall, NormalizeQuantiles, Plots, GR, BenchmarkTools, MLBase, NullableArrays, PyPlot
```

# How do I get the version numbers of installed software?
```bash
conda list | tail -n+3 | awk '{print $1, $2}' > temp.info
brew list --versions >> temp.info
julia --version | awk '{print $1, $3}' >> temp.info
column -t temp.info | sort -u
R
installed.packages(fields = c("Package", "Version"))[,c("Package", "Version")]
quit(save="no")
julia
Pkg.installed()
exit()
```

# Do I have to manually build a container by hand each time I want to do this?
No! But I do recommend the interactive method when trying to install software for the first time, since it'll be easier to troubleshoot any issues. But once you have a reliable installation recipe ready, checkout [the documentation](http://singularity.lbl.gov/#bootstrap) for how to create a definition file for bootstrapping ready-to-go containers!
