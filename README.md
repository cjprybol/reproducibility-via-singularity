# Reproducibility via [Singularity](https://github.com/gmkurtzer/singularity)

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

## Make your research more reproducible with Singularity

Services like [GitHub](https://github.com/), [Bitbucket](https://bitbucket.org), and [GitLab](https://about.gitlab.com/) have democratized access to affordable (or free!) tools that reinforce reproducibility in research, via the code repositories that these services offer. These services make it easier to backup, version control, collaborate on, and distribute code (or any other text-based file). These features make it easier for researchers to write and maintain high-quality code. These features also increase the chance that someone else will review the code for errors or bugs. These reviews can be done by some combination of reading and executing the code. Unfortunately, unless the code is run on the exact same computer, while logged in as the same user, getting the same code to run the same way can be a research project in and of itself.


## Why Singularity?
Universities and research laboratories often conduct their work on shared HPC clusters. These clusters are all set up on slightly different configurations of hardware and operating system software. Trying to recreate an environment to re-run code exactly as it was executed on another cluster, an environment where all of the C pointers align the Java versions sync, up is beyond the abilities, free time, and account privileges of most researchers who may want to try to reproduce your results. Singularity is an implementation of a "container", and an engine to run those containers. Containers allow researchers to isolate the software environment needed to produce a result away from the configuration and operating system of the computer that the analysis will be run on. This means that your colleague at University X can run the analysis exactly the same way on their cluster as you are running it on your cluster at University Y, and all it requires is sharing a git repository (code) and a container image (software).

## TL;DR

**Hardware requirements for this example**
  - RAM ~ 4 GB
  - CPU = 1
  - available disk space ~ 1 GB

Run this code on a linux-based operating system with Singularity installed and in the `$PATH`
```bash
git clone https://github.com/cjprybol/reproducibility-via-singularity.git && \
cd reproducibility-via-singularity && \
wget --no-check-certificate https://stanfordmedicine.box.com/shared/static/fbuap57hb3ywl7f741t1x4mcq9akznl4.img -O demo.img && \
singularity exec demo.img bash reproduce_project.sh
```

**Time to run** < 5 minutes on a modern laptop with ~ 10 MB/s download rates

## What does that code do?

I encourage you to read and execute the `reproduce_project.sh` file and see for yourself! But in summary, it will download 1 sample of paired-end RNA sequencing reads, and quantify transcript isoform abundances in the sample. You can examine the results by running the command `less Kallisto/abundance.tsv` from inside the project directory!

This is a contrived example, designed to run quickly and efficiently to demonstrate a point that reproducing research analyses can be, and should be, be as simple as running 4 code instructions. You could condense the instructions further, but I wanted to keep the steps of acquiring the code and the acquiring the software environment to run that code separate, to highlight the idea of making each modular. By extending the software libraries in the container, and extending the analytical complexity in the analysis scripts, these same 4 lines of code can be sufficient to reproduce entire thesis projects.

Requiring little more effort for others than the wait required to perform the computation, this simplicity of reproducibility lowers the barrier for others to investigate your work. This encourages responsible research conduct and promotes increased rates of knowledge transfer.

## What is Singularity?

A Singularity container is a complete linux kernel capable of executing code and running software. Because the container is a single file, it can be easily moved between computers with Singularity and run in exactly the same way. The Singularity engine on each cluster handles the translation to cluster-specific machine code, ensuring that your code produces the same output, regardless of the underlying configuration. You simply need to install software into the container, and configure it as you would any other linux-based operating system.

## What is Singularity not?

A magic, drag and drop, point and click solution. Singularity requires some understanding of how operating systems work in order to effectively use it and understand it. That being said, it's also very well designed, and should be intuitive for anyone who has set up a computing environment on their personal computer, a cluster, a virtual machine, or a cloud service.

Containers require sudo/root privileges to modify and edit. This means that when creating containers, users must take care to note whether or not the software installed inside of the container will try to download or edit files that exist inside of the directory structure of the container. If a container is shipped without those files pre-downloaded and pre-initialized, end users without sudo/root permissions will experience errors as the software fails to modify files inside of the container. For most situations, the strategy is simply to preload the data! But for some use cases, such as when using software that downloads large databases (> 100GB) for data annotation purposes, you may wish to explore the configuration options of that software. If the software can download data to a directory that will exist *OUTSIDE* of the container, then users only need to worry about sharing containers that may range in size from a few hundred MB to several GB (the disk size of the software), rather than massive images that could easily blossom to sizes larger than a TB.

## How do I interact with a container?

**Example**

If you would normally run a command as:
```bash
CMD --flag1 --flag2 argument1 argument2 ...
```

That command can be executed by the environment inside of the container, instead of by the host, by prepending the `singularity exec` command, followed by the container you wish to use.
```bash
singularity exec <your container>.img CMD --flag1 --flag2 argument1 argument2 ...
```

You can also chain commands together using standard shell syntax you are already used to, such as `|` pipes and `&&` operators

If you compose your analysis as a series of shell scripts that then call other processes, you can just use the shell interpreter within the container.

will run on the host, as usual
```bash
bash reproduce_analysis.sh
```

will run using the bash environment and software within the container
```bash
singularity exec <your container>.img bash reproduce_analysis.sh
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

## How do I build a Singularity image to use with my own project?

Because containers interact so closely with the host system, installing Singularity, and building Singularity images requires sudo/root level permissions. While many users of shared computing clusters are unlikely to have these permission levels on the clusters where they intend to run the analysis, they may have those permissions on whatever local computer (laptop/desktop) they use to ssh into the shared computing cluster. If you happen to have linux-based OS with sudo permissions, go ahead and jump to the section [Start here if linux user with sudo/root privileges](https://github.com/cjprybol/reproducibility-via-singularity#start-here-if-linux-user-with-sudoroot-privileges). If you have a Windows OS or an Apple OS (with OR without sudo/root), or a linux OS without sudo/root, you'll need to get access to a linux machine where you do have sudo/root permissions. For most, the easiest way to do so is to launch a virtual machine that runs inside of your personal computer. If you're already familiar with the process of setting up virtual machines, go ahead with whatever method you know. If installing a virtual machine is a new thing for you, I recommend checking out [VirtualBox](https://www.virtualbox.org/). If you don't have the admin permissions to setup a virtual machine on your computer, ask a system administrator for assistance.

## Starting from a mac with [Homebrew](http://brew.sh/) installed
```bash
brew tap caskroom/cask
brew cask install vagrant
# brew cask install virtualbox
```
**Quick fix**: These instructions don't work with the new VirtualBox 5.1 `brew cask` install, and thus I've had to bail on the current (as of July 16, 2016) cask. Go to [VirtualBox's old build page](https://www.virtualbox.org/wiki/Download_Old_Builds_5_0) and install the appropriate version for your system. Or try the `brew cask` method and let me know that it works again and I'll update these instructions.
```bash
mkdir singularity-vm && cd singularity-vm
# generate a virtual machine
vagrant init ubuntu/trusty64; vagrant up --provider virtualbox && vagrant halt
```
Now we have a virtual machine implemented with vagrant and virtualbox. The virtual machine has some preset defaults that were too restrictive for running this demo. The default memory allocation is was unable to execute the code, however I was able to run the demo after setting the memory allocation to 4GB.

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

Uncomment the code block and the memory line, setting the memory to 4GB.
```
config.vm.provider "virtualbox" do |vb|
#   # Display the VirtualBox GUI when booting the machine
#   vb.gui = true
#
#   # Customize the amount of memory on the VM:
  vb.memory = "4096"
end
```

 You are welcome to increase this further, and it's also possible allocate additional CPUs (if you have them) by adding `vb.cpus = <number of cpus you want to use>` below the `vb.memory` line, like so.
 ```
 config.vm.provider "virtualbox" do |vb|
#   # Display the VirtualBox GUI when booting the machine
#   vb.gui = true
#
#   # Customize the amount of memory on the VM:
  vb.memory = "8192"
  vb.cpus = "2"
end
 ```

Restart the virtual machine and ssh into it
```bash
vagrant up && vagrant ssh
```
You are now a linux user with sudo/root privileges!


## Start here if linux user with sudo/root privileges

Here we will install Singularity starting from an Ubuntu 14.04 LTS "Trusty" 64-bit base installation
```bash
sudo apt-get update && \
sudo apt-get install -y build-essential git vim autoconf libtool curl debootstrap && \
git clone https://github.com/gmkurtzer/singularity.git && \
cd singularity && \
./autogen.sh && \
./configure --prefix=/usr/local && \
make && \
sudo make install && \
cd ..
```

We want to create a Singularity container, load it with an operating system, and install and configure the software necessary to run our analysis onto it.

First, we need to allocate the file. Here `--size` represents the maximum size (in MiB) that the container is allowed to take on. This container is allowed to take on a maximum of 20GiB. **Interesting aside** Containers are initialized as sparse images. If you evaluate the allocated space for the image with `ls -lah test.img`, and compare that disk size to what is returned by `du -h test.img`, you'll see that the files are only keeping track of the informative content installed, rather than the total possible disk space they could use.
```bash
sudo singularity create --size 20000 test.img
```

We will preload a Ubuntu 14.04 LTS "Trusty" 64-bit base install
```bash
wget https://raw.githubusercontent.com/cjprybol/reproducibility-via-singularity/master/ubuntu.def && \
sudo singularity bootstrap test.img ubuntu.def && \
rm ubuntu.def
```
Congratulations, you've made a container!

We have a container, but it doesn't have much software loaded onto it. Let's enter into a bash shell that is running inside of the container. By default, containers have the ability to read and write data to the filesystem of the host, yet containers themselves are immutable. Here, we will launch our shell session with the `--contain` flag, which blocks our ability to interact with the hosting computer, and the `--writable` flag such that we can modify the contents of the container
```bash
sudo singularity shell --writable --contain test.img
```

Now you should be inside of the image. Feel free to jump around the file system to learn your way around.

[Linuxbrew](http://linuxbrew.sh/) and [Anaconda](https://www.continuum.io/downloads) are two great package managers that simplify the process of installing and managing software. The default software available via apt-get is often out of date compared to those available via either Linuxbrew and Anaconda, and additionally, both of these package managers have recipes to install a wide-variety of software useful for researchers.

First, I will add two directories that I found I needed to add to get my Singularity container to work on my University's cluster. On this cluster, the lab spaces are setup as follows:
```
/scratch/PI/{PI_name}/"All lab data and projects go here"
/share/PI/{PI_name}/"All lab software and executables, and some reference files, go here"
```
The base directories for these paths, `/scratch` & `/share` are not standard linux directories, and are not found on the container. I had to add them to the container in order to get the container to resolve paths when I tried to use the container on the cluster. `/local-scratch` was added to suppress a warning `WARNING: Non existant 'bind point' in container: '/local-scratch'` **If you experience a similar issue, let me know. If there are other common directories used on other HPC clusters, I'll add them here**
```bash
mkdir /scratch /share /local-scratch
```

Here we install required system dependencies for other software. The dependencies necessary to install the software you require may be different.
```bash
cd / && \
apt-get update -y && \
apt-get install -y alien build-essential cmake curl ed gdebi-core git libsm6 libxrender1 libfontconfig1 lsb-release nettle-dev python-setuptools ruby software-properties-common vim wget zlib1g-dev
```

Make a new directory for installing additional software
```bash
mkdir /Software && \
cd /Software
```

Install Linuxbrew
```bash
git clone https://github.com/Linuxbrew/brew.git /Software/.linuxbrew
```

Install Anaconda
```bash
wget http://repo.continuum.io/archive/Anaconda3-4.1.1-Linux-x86_64.sh && \
bash Anaconda3-4.1.1-Linux-x86_64.sh -b -p /Software/anaconda3 && \
rm Anaconda3-4.1.1-Linux-x86_64.sh
```

Add the new software to the `$PATH`
```bash
PATH="/Software/.linuxbrew/bin:/Software/anaconda3/bin:$PATH"
```

With each package manager, we will install a few essential and commonly used libraries, and then extend the library of available software by loading additional channels. We will tap into the [Homebrew-science](https://github.com/Homebrew/homebrew-science) channel as well as the [Bioconda channel of anaconda](https://github.com/bioconda/bioconda-recipes/tree/master/recipes). There are other channels available for both package managers as well, with installation recipes for software specific to other domains of research.

Linuxbrew

Temp fix for open-mpi until the linuxbrew bottle comes back
```bash
rm /Software/.linuxbrew/Library/Taps/homebrew/homebrew-core/Formula/open-mpi.rb && \
wget -P /Software/.linuxbrew/Library/Taps/homebrew/homebrew-core/Formula https://raw.githubusercontent.com/Linuxbrew/homebrew-core/4e682fe09ae928bd468a421fc1e8067a54799a3a/Formula/open-mpi.rb && \
cd /Software && \
wget http://linuxbrew.bintray.com/bottles/open-mpi-1.10.2_1.x86_64_linux.bottle.tar.gz && \
cd $(brew --prefix)/Cellar && tar -zxf /Software/open-mpi-1.10.2_1.x86_64_linux.bottle.tar.gz && brew link open-mpi && \
cd /Software && \
rm open-mpi-1.10.2_1.x86_64_linux.bottle.tar.gz
```

```bash
brew tap homebrew/dupes && \
brew install --force-bottle automake bash binutils cmake coreutils curl file-formula findutils gawk gcc git gnu-sed gnu-tar gnu-which grep libtool libgit2 make open-mpi parallel pigz util-linux wget && \
ln -sf /Software/.linuxbrew/bin/bash /bin/bash && \
brew tap homebrew/science && \
brew install abyss art bamtools bcftools beagle bedops bedtools bowtie bowtie2 blat bwa clustal-omega clustal-w exonerate fastq-tools fastqc hisat hmmer htslib igv jellyfish last novoalign openblas picard-tools plink repeatmasker samtools snap-aligner soapdenovo sratoolkit tophat trimmomatic varscan vcflib vcftools velvet && \
rm -r $(brew --cache)
```

Anaconda
```bash
conda update -y conda && \
conda update -y anaconda && \
conda config --add channels r && \
conda config --add channels bioconda && \
conda install -y pyaml pybedtools pyfasta pysam python-igraph pyvcf theano && \
conda install -y --channel https://conda.anaconda.org/conda-forge tensorflow && /Software/anaconda3/bin/pip install keras && \
conda install -y --channel r r && \
conda install -y --channel bioconda cramtools cufflinks cutadapt freebayes gatk impute2 kallisto lighter pindel plink2 rsem sailfish salmon sambamba star trinity && \
conda clean -y --all
```

Install the latest R from CRAN, install and configure packages
```bash
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9 && \
add-apt-repository -y "deb http://cran.rstudio.com/bin/linux/ubuntu $(lsb_release -s -c)/" && \
apt-get update -y && \
apt-get install -y r-base r-base-dev && \
apt-get clean && \
wget --no-check-certificate https://raw.githubusercontent.com/cjprybol/reproducibility-via-singularity/master/install_packages.R && \
chmod 775 install_packages.R && \
R --vanilla --slave < install_packages.R && \
rm install_packages.R
```

Install [Julia](http://julialang.org/). Note: the `MARCH=x86-64` instructs Julia to compile to run on a generic 64-bit CPU, enabling the executable to run on different hardware on different clusters. `USE_SYSTEM_LIBGIT2=1` fixes a bug with Julia not finding the libcurl libraries in linuxbrew.
```bash
cd /Software && \
git clone git://github.com/JuliaLang/julia.git && \
cd julia && \
git checkout release-0.5 && \
make MARCH=x86-64 USE_SYSTEM_LIBGIT2=1 && \
ln -s /Software/julia/julia /usr/local/bin/julia
```

Install [RTG core](http://realtimegenomics.com/products/rtg-core/). **This software is restricted to non-commercial use**. If you intend to use it commercially, you'll have to buy a license through the website, or alternatively just skip this installation.
```bash
cd /Software && \
wget --no-check-certificate https://github.com/RealTimeGenomics/rtg-core/releases/download/3.6.2/rtg-core-non-commercial-3.6.2-linux-x64.zip && \
unzip rtg-core-non-commercial-3.6.2-linux-x64.zip && \
rm rtg-core-non-commercial-3.6.2-linux-x64.zip && \
ln -s /Software/rtg-core-non-commercial-3.6.2/rtg /usr/local/bin && \
echo "n" | rtg
```

Install [Cell Ranger](http://support.10xgenomics.com/single-cell/software/pipelines/latest/what-is-cell-ranger)
```bash
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/bcl2fastq/bcl2fastq2-v2.17.1.14-Linux-x86_64.zip && \
unzip bcl2fastq2-v2.17.1.14-Linux-x86_64.zip && \
rm bcl2fastq2-v2.17.1.14-Linux-x86_64.zip && \
alien -i bcl2fastq2-v2.17.1.14-Linux-x86_64.rpm && \
rm bcl2fastq2-v2.17.1.14-Linux-x86_64.rpm && \
wget --no-check-certificate -O cellranger-1.1.0.tar.gz <get the url by requesting access on the 10x website> && \
tar -xzvf cellranger-1.1.0.tar.gz && \
rm cellranger-1.1.0.tar.gz && \
ln -s /Software/cellranger-1.1.0/cellranger /usr/local/bin
```

Here we will install the [Genome Analysis Toolkit](https://www.broadinstitute.org/gatk/). GATK is license restricted, and you can acquire a copy by going to website and accepting the terms of agreement (and purchase a license, if you're working commercially). If you have the option to host your copy on a private FTP server, you can save yourself a few steps by downloading your copy directly into the container with `wget`.

```bash
wget <ftp address for your copy> && \
gatk-register GenomeAnalysisTK-3.6.tar.bz2 && \
rm GenomeAnalysisTK-3.6.tar.bz2
```

We've got our system fully loaded with the software we want, but our `$PATH` update was only for this session. We'll need to make our `$PATH` updates permanent to make the software installed inside of the `/Software` directory available by name alone (e.g. calling `python3`, rather than `/Software/anaconda3/bin/python3`). In Singularity version >= 2.1, you can update the `$PATH` by modifying the `/environment` file, which is loaded each time you interact with the container. This functionality is not present in earlier versions of Singularity.
```bash
cd / && \
rm /environment && \
wget --no-check-certificate https://raw.githubusercontent.com/cjprybol/reproducibility-via-singularity/master/environment
```

We'll download another script that will list all software installed with version numbers. Note that if you install any software manually from source that is not included in this example, you'll need to update the script to include that software in the list.
```bash
wget --no-check-certificate https://raw.githubusercontent.com/cjprybol/reproducibility-via-singularity/master/singularity && \
chmod 775 singularity
```

That's it! I've tried to cover a wide array of software, but this list won't cover the needs of every user, so go ahead and adjust these installation steps to your needs.

## How do I see the full list of available software with version numbers?
```bash
singularity run test.img
```

And it should print something that looks like this (but not exactly like this. `...` indicate that I'm skipping lines in the output, and these versions are not up to date). The Python packages and the command line tools installed via Anaconda are both listed under `Anaconda`
```bash
Using Anaconda Cloud api site https://api.anaconda.org
abyss                   1.9.0_1       Homebrew
alabaster               0.7.8         Anaconda
anaconda-client         1.4.0         Anaconda
anaconda                custom        Anaconda
anaconda-navigator      1.2.1         Anaconda
argcomplete             1.0.0         Anaconda
argtable                2.13          Homebrew
art                     031915        Homebrew
astropy                 1.1.2         Anaconda
autoconf                2.69          Homebrew
automake                1.15          Homebrew
babel                   2.3.3         Anaconda
backports               1.0           Anaconda
bamtools                2.4.0_1       Homebrew
bash                    4.3.46        Homebrew
bcftools                1.3.1         Anaconda
bcftools                1.3.1         Homebrew
bcl2fastq               v2.17.1.14    User_Install
...                     ...           ...
julia                   0.5.0-rc2     GitHub
...                     ...           ...
R                       3.3.1         CRAN
...                     ...           ...
affy                    1.50.0      R_package
...                     ...           ...
```

## Do I have to manually build a container each time I want to do this?
No! But I do recommend the interactive method when trying to install software for the first time, since it'll be easier to troubleshoot any issues. But once you have a reliable installation recipe ready, checkout [the documentation](http://singularity.lbl.gov/#bootstrap) for how to create a definition file for bootstrapping ready-to-go containers!

## Can I run Jupyter notebooks on remote servers with this?

This will start a notebook on the remote server
```bash
singularity exec test.img jupyter notebook --no-browser --ip \*
```

And on your local computer
```bash
ssh -NL localhost:9999:${remote-node}:8888 your_username@your_domain.com
```

Then just go to the url `localhost:9999` in your web browser.

## I built a container inside of a vagrant VM. How do I get the container out of the VM and onto the server where I work?

Exit out of the container and vagrant VM back to your host computer. From within the same directory where your Vagrantfile is (the directory where you initialized your VM)
```bash
vagrant plugin install vagrant-scp
vagrant scp default:/home/vagrant/test.img .
```

Now the container will be available on your local machine, and you can use `scp` or whatever your preferred method is for moving the container to where you need it!

## Troubleshooting

Just as you might experience problems configuring your user environment on a new cluster, you are likely to experience some issues trying to configure your container for the first time. Even if you install all of the software successfully, there is a chance that when you transfer the container to the cluster on which you work, you may find some new errors when trying to execute your code. It may take a few iterative feedback cycles of testing and debugging before you have it all worked out.

## Naming conventions

Pick what works for you.

If you need an idea, what I'm doing is moving up major versions (i.e. `1.0 -> 2.0 -> 3.0 -> ...`) every time I compile from scratch, where `1.0` is the first working container that passed all of my testing. With every edit to add additional software, fix bigs, or update to more recent versions, I'm bumping up the minor version (i.e. `1.0 -> 1.1 -> 1.2 -> ...`). Upon project completion, versions can be tagged and archived along with the source code and raw data.

## Contributions
Thank you to everyone who has contributed!
- [gmkurtzer](https://github.com/gmkurtzer): building Singularity, and feedback on examples and document contents
- [dwaggott](https://github.com/dwaggott): providing a great list of R packages, and an example of how to install R packages from the command line via an R script
- [kprybol](https://github.com/kprybol): another extensive list of useful R packages, and an explanation of R's package system (and it's use of Suggests, Depends, Imports, LinkingTo, and Enhances)
