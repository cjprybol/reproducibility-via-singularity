# Make your research more reproducible with [Singularity](http://singularity.lbl.gov)
Services like [GitHub](https://github.com/), [Bitbucket](https://bitbucket.org), and [GitLab](https://about.gitlab.com/) have democratized access to affordable (or free!) tools that reinforce reproducibility in research, via the code repositories that these services offer. These services make it easier to backup, version control, collaborate on, and distribute code (or any other text-based file). These features make it easier for researchers to write and maintain high-quality code. These features also increase the chance that someone else will review the code for errors or bugs. These reviews can be done by some combination of reading and executing the code. Unfortunately, unless the code is run on the exact same computer, while logged in as the same user, getting the same code to run the same way can be a research project in and of itself.


# Why Singularity?
Universities and research laboratories often conduct their work on shared HPC clusters. These clusters are all set up on slightly different configurations of hardware and operating system software. Trying to recreate an environment to re-run code exactly as it was executed on another cluster, an environment where all of the C pointers align the Java versions sync, up is beyond the abilities, free time, and account privileges of most researchers who may want to try to reproduce your results. Singularity is an implementation of a "container", and an engine to run those containers. Containers allow researchers to isolate the software environment needed to produce a result away from the configuration and operating system of the computer that the analysis will be run on. This means that your colleague at University X can run the analysis exactly the same way on their cluster as you are running it on your cluster at University Y, and all it requires is sharing a git repository (code) and a container image (software).

# TL;DR

**Hardware requirements for this example**
  - RAM ~ 4 GB
  - CPU = 1
  - available disk space ~ 1 GB

Run this code on a linux-based operating system with Singularity installed and in the `$PATH`
```bash
git clone https://github.com/cjprybol/reproducibility-via-singularity.git
cd reproducibility-via-singularity
wget --no-check-certificate https://stanfordmedicine.box.com/shared/static/fbuap57hb3ywl7f741t1x4mcq9akznl4.img -O demo.img
./reproduce_project.sh
```

**Time to run** < 5 minutes on a modern laptop with ~ 10 MB/s download rates

# What does that code do?

I encourage you to read and execute the `reproduce_project.sh` file and see for yourself! But in summary, it will download 1 sample of paired-end RNA sequencing reads, and quantify transcript isoform abundances in the sample. You can examine the results by running the command `less Kallisto/abundance.tsv` from inside the project directory!

This is a contrived example, designed to run quickly and efficiently to demonstrate a point that reproducible research can be, and should be, be as simple as running 4 code instructions. However, by extending the software libraries in the container, and extending the analytical complexity in the analysis scripts, these same 4 lines of code can be sufficient to reproduce entire thesis projects.

Requiring little more effort for others than the wait required to perform the computation, this simplicity of reproducibility lowers the barrier for others to investigate your work. This encourages responsible research conduct and promotes increased rates of knowledge transfer.

# What is Singularity?

A Singularity container is a complete linux kernel capable of executing code and running software. Because the container is a single file, it can be easily moved between computers with Singularity and run in exactly the same way. The Singularity engine on each cluster handles the translation to cluster-specific machine code, ensuring that your code produces the same output, regardless of the underlying configuration. You simply need to install software into the container, and configure it as you would any other linux-based operating system.

# What is Singularity not?

A magic, drag and drop, point and click solution. Singularity requires some understanding of how operating systems work in order to effectively use it and understand it. That being said, it's also very well designed, and should be intuitive for anyone who has set up a computing environment on their personal computer, a cluster, a virtual machine, or a cloud service.

Containers require sudo/root privileges to modify and edit. This means that when creating containers, users must take care to note whether or not the software installed inside of the container will try to download or edit files that exist inside of the directory structure of the container. If a container is shipped without those files pre-downloaded and pre-initialized, end users without sudo/root permissions will experience errors as the software fails to modify files inside of the container. For most situations, the strategy is simply to preload the data! But for some use cases, such as when using software that downloads large databases (> 100GB) for data annotation purposes, you may wish to explore the configuration options of that software. If the software can download data to a directory that will exist *OUTSIDE* of the container, then users only need to worry about sharing containers that may range in size from a few hundred MB to several GB (the disk size of the software), rather than massive images that could easily blossom to sizes larger than a TB.

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

You can also chain commands together using standard shell syntax you are already used to, such as `|` pipes and `&&` operators

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

# How do I build a Singularity image to use with my own project?

Because containers interact so closely with the host system, installing Singularity, and building Singularity images requires sudo/root level permissions. While many users of shared computing clusters are unlikely to have these permission levels on the clusters where they intend to run the analysis, they may have those permissions on whatever local computer (laptop/desktop) they use to ssh into the shared computing cluster. If you happen to have linux-based OS with sudo permissions, go ahead and jump to the section [Start here if linux user with sudo/root privileges](https://github.com/cjprybol/reproducibility-via-singularity#start-here-if-linux-user-with-sudoroot-privileges). If you have a Windows OS or an Apple OS (with OR without sudo/root), or a linux OS without sudo/root, you'll need to get access to a linux machine where you do have sudo/root permissions. For most, the easiest way to do so is to launch a virtual machine that runs inside of your personal computer. If you're already familiar with the process of setting up virtual machines, go ahead with whatever method you know. If installing a virtual machine is a new thing for you, I recommend checking out [VirtualBox](https://www.virtualbox.org/). If you don't have the admin permissions to setup a virtual machine on your computer, ask a system administrator for assistance.

# Starting from a mac with [Homebrew](http://brew.sh/) installed
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

Uncomment the code block and the memory line, setting the memory to 4GB
```
config.vm.provider "virtualbox" do |vb|
#   # Display the VirtualBox GUI when booting the machine
#   vb.gui = true
#
#   # Customize the amount of memory on the VM:
  vb.memory = "4096"
end
```

Restart the virtual machine and ssh into it
```bash
vagrant up && vagrant ssh
```
You are now a linux user with sudo/root privileges!


# Start here if linux user with sudo/root privileges

Here we will install Singularity starting from an Ubuntu 14.04 LTS "Trusty" 64-bit base installation
```bash
sudo apt-get install -y build-essential git vim autoconf libtool curl debootstrap
git clone https://github.com/gmkurtzer/singularity.git && cd singularity && ./autogen.sh && ./configure --prefix=/usr/local && make && sudo make install
```

We want to create a Singularity container, load it with an operating system, and install and configure the software necessary to run our analysis onto it.

First, we need to allocate the file. Here `--size` represents the maximum size (in MiB) that the container is allowed to take on. This container is allowed to take on a maximum of 15GiB. **Interesting aside** Containers are initialized as sparse images. If you evaluate the allocated space for the image with `ls -lah my_container.img`, and compare that disk size to what is returned by `du -h my_container.img`, you'll see that the files are only keeping track of the informative content installed, rather than the total possible disk space they could use.
```bash
sudo singularity create --size 15000 test.img
```

We will preload a Ubuntu 14.04 LTS "Trusty" 64-bit base install
```bash
sudo singularity bootstrap test.img $HOME/singularity/examples/ubuntu.def
```
Congratulations, you've made a container!

We have a container, but it doesn't have much software loaded onto it. Let's enter into a bash shell that is running inside of the container. By default, containers have the ability to read and write data to the filesystem of the host, yet containers themselves are immutable. Here, we will launch our shell session with the `--contain` flag, which blocks our ability to interact with the hosting computer, and the `--writable` flag such that we can modify the contents of the container
```bash
sudo singularity shell --writable --contain test.img
```

Now you should be inside of the image. Feel free to jump around the file system to learn your way around.

[Linuxbrew](http://linuxbrew.sh/) and [Anaconda](https://www.continuum.io/downloads) are two great package managers that simplify the process of installing and managing software. The default software available via apt-get is often out of date compared to those available via either Linuxbrew and Anaconda, and additionally, both of these package managers have recipes to install a wide-variety of software useful for researchers.

By default, when entering into a shell session within a container, you stay in the same directory on the host where you started from. Run `pwd` and notice the `(unreachable)` pre-pension to the path. Our current directory is immutable because we have both remained in a directory on the host computer, but also specified that we want to `--contain` our actions to only interact with the contents of the container.

Let's change to the user directory of our current user, which is `root`, inside of the container so we have a mutable space that allows for temporary file generation.
```bash
cd /root
```

First, I will add two directories that I found I needed to add to get my Singularity container to work on my University's cluster. On this cluster, the lab spaces are setup as follows:
```
/scratch/PI/{PI_name}/"All lab data and projects go here"
/share/PI/{PI_name}/"All lab software and executables, and some reference files, go here"
```
The base directories for these paths, `/scratch` & `/share` are not standard linux directories, and are not found on the container. I had to add them to the container in order to get the container to resolve paths when I tried to use the container on the cluster. **If you experience a similar issue, let me know. If there are other common directories used on other HPC clusters, I'll add them here**
```bash
mkdir /scratch /share
```

Here we install required system dependencies for other software. The dependencies necessary to install the software you require may be different.
```bash
apt-get update && apt-get install -y build-essential cmake curl wget git python-setuptools ruby nettle-dev ed && apt-get clean
```

Make a new directory for installing additional software
```bash
mkdir /Software && cd /Software
```

Install Linuxbrew
```bash
git clone https://github.com/Linuxbrew/brew.git /Software/.linuxbrew
```

Install Anaconda
```bash
wget http://repo.continuum.io/archive/Anaconda3-4.1.0-Linux-x86_64.sh && bash Anaconda3-4.1.0-Linux-x86_64.sh -b -p /Software/anaconda3 && rm Anaconda3-4.1.0-Linux-x86_64.sh
```

Add the new software to the `$PATH`
```bash
PATH="/Software/.linuxbrew/bin:/Software/anaconda3/bin:$PATH"
```

Now let's use our package managers to quickly and easily install and configure software into the container. You can review many additional software offerings here -> [Homebrew-science](https://github.com/Homebrew/homebrew-science) & [Bioconda channel of anaconda](https://github.com/bioconda/bioconda-recipes/tree/master/recipes), although these channels are just the tip of the iceberg for what's available!

Linuxbrew
```bash
brew install --force-bottle open-mpi && brew install automake cmake curl git libtool parallel pigz wget
brew tap homebrew/science && brew install abyss art bamtools bcftools beagle bedops bedtools bowtie bowtie2 blat bwa clustal-omega clustal-w exonerate fastq-tools fastqc gmap-gsnap hisat hmmer htslib igv jellyfish kallisto last lighter novoalign openblas picard-tools plink r repeatmasker samtools snap-aligner snpeff soapdenovo sratoolkit tophat trimmomatic varscan vcflib vcfanno vcftools velvet
rm -r $(brew --cache)
```

Anaconda
```bash
conda update -y conda && conda update -y anaconda
conda config --add channels r && conda config --add channels bioconda
conda install -y pyaml pybedtools pyfasta pysam python-igraph pyvcf theano
conda install -y --channel https://conda.anaconda.org/conda-forge tensorflow
# linuxbrew is in path BEFORE anaconda, so must call pip for anaconda python3 by full path
/Software/anaconda3/bin/pip install keras
conda install -y --channel r r && conda install -y --channel bioconda cramtools cufflinks cutadapt freebayes gatk impute2 pindel plink2 rsem sailfish salmon sambamba star trinity && conda clean -y --all
```

Setup R packages
```bash
wget --no-check-certificate https://raw.githubusercontent.com/cjprybol/reproducibility-via-singularity/master/install_packages.R && chmod 775 install_packages.R && ./install_packages.R
rm install_packages.R
```

Install [Julia](http://julialang.org/)
```bash
wget https://julialang.s3.amazonaws.com/bin/linux/x64/0.4/julia-0.4.6-linux-x86_64.tar.gz && tar xfz julia-0.4.6-linux-x86_64.tar.gz && rm julia-0.4.6-linux-x86_64.tar.gz && ln -s /Software/julia-*/bin/julia /usr/local/bin/julia
```

Install [RTG core](http://realtimegenomics.com/products/rtg-core/). **NOTE** This software is license restricted. It's free for non-commercial academic use, but if you intend to use it commercially you'll have to buy a license (alternatively, just skip this installation).
```bash
cd /Software && wget --no-check-certificate https://github.com/RealTimeGenomics/rtg-core/releases/download/3.6.2/rtg-core-non-commercial-3.6.2-linux-x64.zip && unzip rtg-core-non-commercial-3.6.2-linux-x64.zip && rm rtg-core-non-commercial-3.6.2-linux-x64.zip
ln -s /Software/rtg-core-non-commercial-3.6.2/rtg /usr/local/bin
```

Here is our first example of a configuration step that needs to be performed **BEFORE** trying to use the container without sudo/root permissions on the cluster. The first time you run RTG, it will ask whether or not it can perform logging to help the developers improve the software. Because we intend to run the container without sudo and not in `--writable` mode, any attempts RTG makes to save log files to disk will fail, so say no. RTG will save your answer to a config file inside of the directory where RTG is installed, so if you try this without using the `--writable` flag, it will fail.
```bash
Singularity.test.img> rtg
RTG has a facility to automatically send basic usage information to Real
Time Genomics. This does not contain confidential information such as
command-line parameters or dataset contents.

Would you like to enable automatic usage logging (y/n)? n
```

We've got our system fully loaded with the software we want, but our `$PATH` update was only for this session. We'll need to make our `$PATH` updates permanent to make the software installed inside of the `/Software` directory available by name alone (e.g. calling `python3`, rather than `/Software/anaconda3/bin/python3`). In Singularity version >= 2.1, you can update the `$PATH` by modifying the `/environment` file, which is loaded each time you interact with the container. This functionality is not present in earlier versions of Singularity.
```bash
cd / && rm /environment && wget --no-check-certificate https://raw.githubusercontent.com/cjprybol/reproducibility-via-singularity/master/environment
```

We'll download one more pre-written script that will list all software installed with version numbers. Note that if you install any software manually from source that is not included in this example, you'll need to update the script to include that software in the list.
```bash
wget --no-check-certificate https://raw.githubusercontent.com/cjprybol/reproducibility-via-singularity/master/singularity && chmod 775 singularity
```

Exit the container to the host linux
```bash
exit
```

Another tool I use that cannot be installed programmatically due to licensing restrictions is [Genome Analysis Toolkit](https://www.broadinstitute.org/gatk/). You'll need to make an account and accept some terms of agreement before you can access the software. I've downloaded the GATK installer to the same directory on the host computer where my container is. By entering the container again without the `--contain` flag, we allow the container to interact with the host system, and copy GATK into Anaconda. However, if you download the file and then host it on a personal server, you can automate this step with a `wget` or `curl` command too! This step will install gatk into our Anaconda software library, so you can download and link your file immediately after the last `conda install` step.
```bash
sudo singularity shell --writable test.img
gatk-register /home/vagrant/GenomeAnalysisTK-3.6.tar.bz2
rm GenomeAnalysisTK-3.6.tar.bz2
exit
```

You're all done, you've built a great base-image for computational genomics! Adjust these installation steps to your needs.

# How do I see the full list of available software with version numbers?
```bash
singularity run test.img
```

Output
```bash
vagrant@vagrant-ubuntu-trusty-64:~$ singularity run test.img
Using Anaconda Cloud api site https://api.anaconda.org
abyss                   1.9.0         Homebrew
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
bcftools                1.3.1         Anaconda
bcftools                1.3.1         Homebrew
beagle                  2.1.2         Homebrew
beautifulsoup4          4.4.1         Anaconda
bedops                  2.4.19        Homebrew
bedtools                2.24.0        Homebrew
bigreqsproto            1.1.2         Homebrew
bitarray                0.8.1         Anaconda
blast                   2.4.0         Homebrew
blat                    36            Homebrew
blaze                   0.10.1        Anaconda
bokeh                   0.12.0        Anaconda
boost                   1.60.0_2      Homebrew
boost                   1.60.0        Anaconda
boto                    2.40.0        Anaconda
bottleneck              1.0.0         Anaconda
bowtie                  1.1.2_1       Homebrew
bowtie                  1.1.2         Anaconda
bowtie2                 2.2.9         Homebrew
bwa                     0.7.13        Homebrew
bzip2                   1.0.6_1       Homebrew
bzip2                   1.0.6         Anaconda
cairo                   1.12.18       Anaconda
cffi                    1.6.0         Anaconda
chest                   0.2.3         Anaconda
click                   6.6           Anaconda
cloudpickle             0.2.1         Anaconda
clustal-omega           1.2.1         Homebrew
clustal-w               2.1           Homebrew
clyent                  1.2.2         Anaconda
cmake                   3.6.0         Homebrew
collectl                4.0.4         Anaconda
colorama                0.3.7         Anaconda
compositeproto          0.4.2         Homebrew
conda                   4.1.9         Anaconda
conda-build             1.21.2        Anaconda
conda-env               2.5.2         Anaconda
configobj               5.0.6         Anaconda
contextlib2             0.5.3         Anaconda
cramtools               3.0.b127      Anaconda
cryptography            1.4           Anaconda
cufflinks               2.2.1         Anaconda
curl                    7.49.0        Anaconda
curl                    7.49.1        Homebrew
cutadapt                1.10          Anaconda
cycler                  0.10.0        Anaconda
cython                  0.24          Anaconda
cytoolz                 0.8.0         Anaconda
damageproto             1.2.1         Homebrew
dask                    0.10.0        Anaconda
datashape               0.5.2         Anaconda
decorator               4.0.10        Anaconda
dill                    0.2.5         Anaconda
dmxproto                2.3.1         Homebrew
docutils                0.12          Anaconda
dri2proto               2.8           Homebrew
dri3proto               1.0           Homebrew
dynd-python             0.7.2         Anaconda
entrypoints             0.2.2         Anaconda
et_xmlfile              1.0.1         Anaconda
exonerate               2.2.0         Homebrew
expat                   2.2.0         Homebrew
fastcache               1.0.2         Anaconda
fastool                 0.1.4         Anaconda
fastqc                  0.11.3        Homebrew
fastq-tools             0.7           Homebrew
fixesproto              5.0           Homebrew
flask                   0.11.1        Anaconda
flask-cors              2.1.2         Anaconda
fontconfig              2.11.1        Anaconda
fontsproto              2.1.3         Homebrew
font-util               1.3.1         Homebrew
freebayes               1.0.2.29      Anaconda
freetype                2.5.5         Anaconda
freetype                2.6.5         Homebrew
gatk                    3.6           Anaconda
gcc                     5.3.0         Homebrew
gdbm                    1.12          Homebrew
get_terminal_size       1.0.0         Anaconda
gettext                 0.19.8.1      Homebrew
gevent                  1.1.1         Anaconda
git                     2.9.2         Homebrew
glib                    2.43.0        Anaconda
glib                    2.48.1        Homebrew
glproto                 1.4.17        Homebrew
gmap-gsnap              2016-05-01    Homebrew
gmp                     6.1.1         Homebrew
go                      1.6.3         Homebrew
gpatch                  2.7.5         Homebrew
greenlet                0.4.10        Anaconda
gsl                     1.16          Anaconda
gsl                     1.16          Homebrew
h5py                    2.6.0         Anaconda
harfbuzz                0.9.35        Anaconda
hdf5                    1.8.16_1      Homebrew
hdf5                    1.8.16        Anaconda
heapdict                1.0.0         Anaconda
hisat                   0.1.6b        Homebrew
hmmer                   3.1b2_1       Homebrew
htslib                  1.3.1         Anaconda
htslib                  1.3.1         Homebrew
icu                     54.1          Anaconda
idna                    2.1           Anaconda
igv                     2.3.68        Homebrew
imagesize               0.7.1         Anaconda
impute2                 2.3.2         Anaconda
inputproto              2.3.1         Homebrew
ipykernel               4.3.1         Anaconda
ipython                 4.2.0         Anaconda
ipython_genutils        0.1.0         Anaconda
ipywidgets              4.1.1         Anaconda
isl                     0.15          Homebrew
itsdangerous            0.24          Anaconda
java-jdk                8.0.45        Anaconda
jbig                    2.1           Anaconda
jdcal                   1.2           Anaconda
jdk                     1.8.0-60      Homebrew
jedi                    0.9.0         Anaconda
jellyfish               2.2.3         Anaconda
jellyfish               2.2.4         Homebrew
jinja2                  2.8           Anaconda
jpeg                    8d            Anaconda
jpeg                    8d            Homebrew
jsonschema              2.5.1         Anaconda
julia                   0.4.7-pre     User_Install
jupyter                 1.0.0         Anaconda
jupyter_client          4.3.0         Anaconda
jupyter_console         4.1.1         Anaconda
jupyter_core            4.1.0         Anaconda
kallisto                0.42.5        Homebrew
kbproto                 1.0.7         Homebrew
Keras                   1.0.6         Anaconda
last                    581           Homebrew
libdmx                  1.1.3         Homebrew
libdynd                 0.7.2         Anaconda
libedit                 20150325-3.1  Homebrew
libevent                2.0.22        Homebrew
libffi                  3.0.13        Homebrew
libffi                  3.2.1         Anaconda
libfontenc              1.1.3         Homebrew
libfs                   1.0.7         Homebrew
libgcc                  4.8.5         Anaconda
libgfortran             3.0.0         Anaconda
libice                  1.0.9         Homebrew
libidn                  1.32          Homebrew
libmagic                5.28          Homebrew
libmpc                  1.0.3         Homebrew
libpciaccess            0.13.4        Homebrew
libpng                  1.6.22        Anaconda
libpng                  1.6.23        Homebrew
libsm                   1.2.2         Homebrew
libsodium               1.0.10        Anaconda
libtiff                 4.0.6_1       Homebrew
libtiff                 4.0.6         Anaconda
libtool                 2.4.6_1       Homebrew
libx11                  1.6.3         Homebrew
libxau                  1.0.8         Homebrew
libxaw                  1.0.13        Homebrew
libxcb                  1.11.1        Homebrew
libxcomposite           0.4.4         Homebrew
libxcursor              1.1.14        Homebrew
libxdamage              1.1.4         Homebrew
libxdmcp                1.1.2         Homebrew
libxext                 1.3.3         Homebrew
libxfixes               5.0.1         Homebrew
libxfont                1.5.1         Homebrew
libxft                  2.3.2         Homebrew
libxi                   1.7.6         Homebrew
libxinerama             1.1.3         Homebrew
libxkbfile              1.0.9         Homebrew
libxml2                 2.9.2         Anaconda
libxml2                 2.9.4         Homebrew
libxmu                  1.1.2         Homebrew
libxpm                  3.5.11        Homebrew
libxrandr               1.5.0         Homebrew
libxrender              0.9.9         Homebrew
libxres                 1.0.7         Homebrew
libxscrnsaver           1.2.2         Homebrew
libxshmfence            1.2           Homebrew
libxslt                 1.1.28        Anaconda
libxt                   1.1.5         Homebrew
libxtst                 1.2.2         Homebrew
libxv                   1.0.10        Homebrew
libxvmc                 1.0.9         Homebrew
libxxf86dga             1.1.4         Homebrew
libxxf86vm              1.1.4         Homebrew
lighter                 1.0.7         Homebrew
llvmlite                0.11.0        Anaconda
locket                  0.2.0         Anaconda
lxml                    3.6.0         Anaconda
m4                      1.4.17        Homebrew
makedepend              1.0.5         Homebrew
markupsafe              0.23          Anaconda
matplotlib              1.5.1         Anaconda
mistune                 0.7.2         Anaconda
mkl                     11.3.3        Anaconda
mkl-service             1.1.2         Anaconda
mpfr                    3.1.4         Homebrew
mpmath                  0.19          Anaconda
multipledispatch        0.4.8         Anaconda
mysql                   5.7.13        Homebrew
nb_anacondacloud        1.1.0         Anaconda
nb_conda                1.1.0         Anaconda
nb_conda_kernels        1.0.3         Anaconda
nbconvert               4.2.0         Anaconda
_nb_ext_conf            0.2.0         Anaconda
nbformat                4.0.1         Anaconda
nbpresent               3.0.2         Anaconda
ncurses                 5.9           Anaconda
ncurses                 6.0_1         Homebrew
networkx                1.11          Anaconda
nltk                    3.2.1         Anaconda
nose                    1.3.7         Anaconda
notebook                4.2.1         Anaconda
novoalign               3.02.13       Homebrew
numba                   0.26.0        Anaconda
numexpr                 2.5.2         Anaconda
numpy                   1.10.4        Anaconda
odo                     0.5.0         Anaconda
openblas                0.2.14        Anaconda
openblas                0.2.18_2      Homebrew
open-mpi                1.10.2_1      Homebrew
openpyxl                2.3.2         Anaconda
openssl                 1.0.2h_1      Homebrew
openssl                 1.0.2h        Anaconda
pandas                  0.18.1        Anaconda
pango                   1.36.8        Anaconda
parafly                 r2013_01_21   Anaconda
parallel                20160622      Homebrew
partd                   0.3.4         Anaconda
patchelf                0.9_1         Homebrew
patchelf                0.9           Anaconda
pathlib2                2.1.0         Anaconda
path.py                 8.2.1         Anaconda
patsy                   0.4.1         Anaconda
pcre                    8.39          Anaconda
pcre                    8.39          Homebrew
pep8                    1.7.0         Anaconda
perl-app-cpanminus      1.7039        Anaconda
perl-file-find-rule     0.34          Anaconda
perl-module-build       0.4214        Anaconda
perl-number-compare     0.03          Anaconda
perl-pathtools          3.40          Anaconda
perl-scalar-list-utils  1.42          Anaconda
perl-text-glob          0.09          Anaconda
perl-threaded           5.22.0        Anaconda
pexpect                 4.0.1         Anaconda
picard-tools            2.2.4         Homebrew
pickleshare             0.7.2         Anaconda
pigz                    2.3.3         Homebrew
pillow                  3.2.0         Anaconda
pindel                  0.2.5b8       Anaconda
pip                     8.1.2         Anaconda
pixman                  0.32.6        Anaconda
pkg-config              0.29.1_1      Homebrew
plink                   1.07          Homebrew
plink2                  1.90b3.35     Anaconda
ply                     3.8           Anaconda
presentproto            1.0           Homebrew
psutil                  4.3.0         Anaconda
ptyprocess              0.5.1         Anaconda
py                      1.4.31        Anaconda
pyaml                   15.8.2        Anaconda
pyasn1                  0.1.9         Anaconda
pybedtools              0.7.8         Anaconda
pycosat                 0.6.1         Anaconda
pycparser               2.14          Anaconda
pycrypto                2.6.1         Anaconda
pycurl                  7.43.0        Anaconda
pyfasta                 0.5.2         Anaconda
pyflakes                1.2.3         Anaconda
pygments                2.1.3         Anaconda
pyopenssl               0.16.0        Anaconda
pyparsing               2.1.4         Anaconda
pyqt                    4.11.4        Anaconda
pysam                   0.9.0         Anaconda
pytables                3.2.2         Anaconda
pytest                  2.9.2         Anaconda
python                  2.7.12        Homebrew
python                  3.5.2         Anaconda
python-dateutil         2.5.3         Anaconda
python-igraph           0.7.1.post6   Anaconda
pytz                    2016.4        Anaconda
pyvcf                   0.6.8         Anaconda
pyyaml                  3.11          Anaconda
pyzmq                   15.2.0        Anaconda
qt                      4.8.7         Anaconda
qtconsole               4.2.1         Anaconda
qtpy                    1.0.2         Anaconda
r                       3.3.1_1       Homebrew
r                       3.3.1         Anaconda
randrproto              1.5.0         Homebrew
r-base                  3.3.1         Anaconda
r-boot                  1.3_18        Anaconda
r-class                 7.3_14        Anaconda
r-cluster               2.0.4         Anaconda
r-codetools             0.2_14        Anaconda
readline                6.2           Anaconda
readline                6.3.8_1       Homebrew
recordproto             1.14.2        Homebrew
redis                   3.2.0         Anaconda
redis-py                2.10.5        Anaconda
renderproto             0.11.1        Homebrew
repeatmasker            4.0.5         Homebrew
requests                2.10.0        Anaconda
resourceproto           1.2.0         Homebrew
r-foreign               0.8_66        Anaconda
r-kernsmooth            2.23_15       Anaconda
r-lattice               0.20_33       Anaconda
r-mass                  7.3_45        Anaconda
r-matrix                1.2_6         Anaconda
rmblast                 2.2.28_1      Homebrew
r-mgcv                  1.8_12        Anaconda
r-nlme                  3.1_128       Anaconda
r-nnet                  7.3_12        Anaconda
rope                    0.9.4         Anaconda
r-recommended           3.3.1         Anaconda
r-rpart                 4.1_10        Anaconda
rsem                    1.2.28        Anaconda
r-spatial               7.3_11        Anaconda
r-survival              2.39_4        Anaconda
RTG                     3.6.2         User_Install
ruamel_yaml             0.11.7        Anaconda
sailfish                0.10.1        Anaconda
salmon                  0.6.0         Anaconda
sambamba                0.6.3         Anaconda
samtools                0.1.19        Anaconda
samtools                1.3.1         Homebrew
scikit-image            0.12.3        Anaconda
scikit-learn            0.17.1        Anaconda
scipy                   0.17.1        Anaconda
scrnsaverproto          1.2.2         Homebrew
setuptools              23.0.0        Anaconda
simplegeneric           0.8.1         Anaconda
singledispatch          3.4.0.3       Anaconda
sip                     4.16.9        Anaconda
six                     1.10.0        Anaconda
slclust                 02022010      Anaconda
snap-aligner            0.15          Homebrew
snowballstemmer         1.2.1         Anaconda
snpeff                  4.2           Homebrew
soapdenovo              2.04.r240     Homebrew
sockjs-tornado          1.0.3         Anaconda
sphinx                  1.4.1         Anaconda
sphinx_rtd_theme        0.1.9         Anaconda
spyder                  2.3.9         Anaconda
sqlalchemy              1.0.13        Anaconda
sqlite                  3.13.0        Anaconda
sqlite                  3.13.0        Homebrew
sratoolkit              2.5.4         Homebrew
star                    2.5.2a        Anaconda
statsmodels             0.6.1         Anaconda
sympy                   1.0           Anaconda
szip                    2.1           Homebrew
tbb                     4.4_20150728  Anaconda
tbb                     4.4-20160128  Homebrew
terminado               0.6           Anaconda
texinfo                 6.1           Homebrew
theano                  0.7.0         Anaconda
tk                      8.5.18        Anaconda
toolz                   0.8.0         Anaconda
tophat                  2.1.1         Homebrew
tornado                 4.3           Anaconda
traitlets               4.2.1         Anaconda
trf                     4.07b         Homebrew
trimmomatic             0.35          Homebrew
trimmomatic             0.36          Anaconda
trinity                 2.2.0         Anaconda
unicodecsv              0.14.1        Anaconda
unzip                   6.0_2         Homebrew
util-linux              2.27.1        Homebrew
util-macros             1.19.0        Homebrew
varscan                 2.3.7         Homebrew
vcfanno                 0.0.2         Homebrew
vcflib                  1.0.0_1       Homebrew
vcftools                0.1.13        Homebrew
velvet                  1.2.10        Homebrew
videoproto              2.3.3         Homebrew
werkzeug                0.11.10       Anaconda
wget                    1.18          Homebrew
wheel                   0.29.0        Anaconda
xcmiscproto             1.2.2         Homebrew
xextproto               7.3.0         Homebrew
xf86bigfontproto        1.2.0         Homebrew
xf86dgaproto            2.1           Homebrew
xf86driproto            2.1.1         Homebrew
xf86vidmodeproto        2.3.1         Homebrew
xineramaproto           1.2.1         Homebrew
xlrd                    1.0.0         Anaconda
xlsxwriter              0.9.2         Anaconda
xlwt                    1.1.2         Anaconda
xorg                    20160314      Homebrew
xorg-protocols          latest        Homebrew
xproto                  7.0.28        Homebrew
xtrans                  1.3.5         Homebrew
xz                      5.2.2         Anaconda
xz                      5.2.2         Homebrew
yaml                    0.1.6         Anaconda
zeromq                  4.1.4         Anaconda
zlib                    1.2.8         Anaconda
zlib                    1.2.8         Homebrew

abind                   1.4-3        R_package
acepack                 1.3-3.3      R_package
aCGH                    1.50.0       R_package
adabag                  4.1          R_package
ade4                    1.7-4        R_package
affxparser              1.44.0       R_package
affy                    1.50.0       R_package
affyio                  1.42.0       R_package
ALL                     1.14.0       R_package
annotate                1.50.0       R_package
AnnotationDbi           1.34.4       R_package
AnnotationForge         1.14.2       R_package
AnnotationHub           2.4.2        R_package
argparse                1.0.1        R_package
assertthat              0.1          R_package
AssotesteR              0.1-10       R_package
AUC                     0.3.0        R_package
base                    3.3.1        R_package
base64enc               0.1-3        R_package
BayesX                  0.2-9        R_package
BB                      2014.10-1    R_package
bbmle                   1.0.18       R_package
bdsmatrix               1.3-2        R_package
BH                      1.60.0-2     R_package
bibtex                  0.4.0        R_package
biganalytics            1.1.14       R_package
biglm                   0.9-1        R_package
bigmemory               4.5.19       R_package
bigmemory.sri           0.1.3        R_package
bigRR                   1.3-10       R_package
bigtabulate             1.1.5        R_package
binGroup                1.1-0        R_package
Biobase                 2.32.0       R_package
BiocGenerics            0.18.0       R_package
BiocInstaller           1.22.3       R_package
BiocParallel            1.6.3        R_package
BiocStyle               2.0.2        R_package
biomaRt                 2.28.0       R_package
Biostrings              2.40.2       R_package
bit                     1.1-12       R_package
bit64                   0.9-5        R_package
bitops                  1.0-6        R_package
boot                    1.3-18       R_package
bootstrap               2015.2       R_package
Boruta                  5.0.0        R_package
BradleyTerry2           1.0-6        R_package
brew                    1.0-6        R_package
brglm                   0.5-9        R_package
broom                   0.4.1        R_package
BSgenome                1.40.1       R_package
bst                     0.3-13       R_package
btergm                  1.7.6        R_package
C50                     0.1.0-24     R_package
car                     2.1-2        R_package
caret                   6.0-70       R_package
caretEnsemble           2.0.0        R_package
caTools                 1.17.1       R_package
cgdsr                   1.2.5        R_package
changepoint             2.2.1        R_package
ChIPpeakAnno            3.6.5        R_package
chron                   2.3-47       R_package
Ckmeans.1d.dp           3.4.6        R_package
class                   7.3-14       R_package
cluster                 2.0.4        R_package
clValid                 0.6-6        R_package
coda                    0.18-1       R_package
codetools               0.2-14       R_package
coin                    1.1-2        R_package
colorspace              1.2-6        R_package
combinat                0.0-8        R_package
compiler                3.3.1        R_package
ConsensusClusterPlus    1.36.0       R_package
corpcor                 1.6.8        R_package
covr                    2.1.0        R_package
crayon                  1.3.2        R_package
Cubist                  0.0.18       R_package
curl                    0.9.7        R_package
DAAG                    1.22         R_package
DatABEL                 0.9-6        R_package
datasets                3.3.1        R_package
data.table              1.9.6        R_package
DBI                     0.4-1        R_package
deepSNV                 1.18.1       R_package
deldir                  0.1-12       R_package
DEoptimR                1.0-6        R_package
DESeq                   1.24.0       R_package
DESeq2                  1.12.3       R_package
deSolve                 1.13         R_package
devtools                1.12.0       R_package
dfoptim                 2016.7-1     R_package
DiagrammeR              0.8.4        R_package
dichromat               2.0-0        R_package
digest                  0.6.9        R_package
diptest                 0.75-7       R_package
dlm                     1.1-4        R_package
DNAcopy                 1.46.0       R_package
doMC                    1.3.4        R_package
doParallel              1.0.10       R_package
dplyr                   0.5.0        R_package
DT                      0.1          R_package
dtplyr                  0.0.1        R_package
dygraphs                0.9          R_package
e1071                   1.6-7        R_package
earth                   4.4.4        R_package
edgeR                   3.14.0       R_package
elasticnet              1.1          R_package
ellipse                 0.3-8        R_package
ensembldb               1.4.7        R_package
epitools                0.5-7        R_package
ergm                    3.6.0        R_package
ergm.count              3.2.2        R_package
estimability            1.1-1        R_package
evaluate                0.9          R_package
expm                    0.999-0      R_package
fastICA                 1.2-0        R_package
fastmatch               1.0-4        R_package
fasttime                1.0-1        R_package
fBasics                 3011.87      R_package
feather                 0.0.1        R_package
ff                      2.2-13       R_package
ffbase                  0.12.3       R_package
fGarch                  3010.82      R_package
fields                  8.4-1        R_package
filehash                2.3          R_package
findpython              1.0.1        R_package
flexmix                 2.3-13       R_package
flexsurv                1.0.0        R_package
foreach                 1.4.3        R_package
forecast                7.1          R_package
foreign                 0.8-66       R_package
formatR                 1.4          R_package
Formula                 1.2-1        R_package
fracdiff                1.4-2        R_package
fts                     0.9.9        R_package
futile.logger           1.4.3        R_package
futile.options          1.0.0        R_package
gam                     1.12         R_package
gamlss                  4.4-0        R_package
gamlss.data             4.3-4        R_package
gamlss.dist             4.3-6        R_package
gamm4                   0.2-3        R_package
gapminder               0.2.0        R_package
gbm                     2.1.1        R_package
gclus                   1.3.1        R_package
gdata                   2.17.0       R_package
geepack                 1.2-0.2      R_package
GenABEL                 1.8-0        R_package
GenABEL.data            1.0.0        R_package
genefilter              1.54.2       R_package
geneplotter             1.50.0       R_package
genetics                1.3.8.1      R_package
GenomeInfoDb            1.8.3        R_package
GenomicAlignments       1.8.4        R_package
GenomicFeatures         1.24.4       R_package
GenomicRanges           1.24.2       R_package
geosphere               1.5-5        R_package
GEOsubmission           1.24.0       R_package
GERGM                   0.7.4        R_package
getopt                  1.20.0       R_package
GGally                  1.2.0        R_package
ggfortify               0.2.0        R_package
ggplot2                 2.1.0        R_package
ggplot2movies           0.0.1        R_package
ggrepel                 0.5          R_package
ggthemes                3.2.0        R_package
ggvis                   0.4.2        R_package
git2r                   0.15.0       R_package
glmnet                  2.0-5        R_package
gmailr                  0.7.1        R_package
gmm                     1.5-2        R_package
gmodels                 2.16.2       R_package
GO.db                   3.3.0        R_package
googleVis               0.6.0        R_package
gplots                  3.0.1        R_package
graph                   1.50.0       R_package
graphics                3.3.1        R_package
grDevices               3.3.1        R_package
grid                    3.3.1        R_package
gridBase                0.4-7        R_package
gridExtra               2.2.1        R_package
gss                     2.1-5        R_package
gsubfn                  0.6-6        R_package
gtable                  0.2.0        R_package
gtools                  3.5.0        R_package
haplo.stats             1.7.7        R_package
haven                   0.2.1        R_package
hdi                     0.1-6        R_package
hexbin                  1.27.1       R_package
hglm                    2.1-1        R_package
hglm.data               1.0-0        R_package
highlight               0.4.7        R_package
highr                   0.6          R_package
Hmisc                   3.17-4       R_package
HMMcopy                 1.14.0       R_package
HSAUR2                  1.1-14       R_package
htmltools               0.3.5        R_package
htmlwidgets             0.6          R_package
httpuv                  1.3.3        R_package
httr                    1.2.1        R_package
hunspell                1.4.2        R_package
ibdreg                  0.2.5        R_package
idr                     1.2          R_package
igraph                  1.0.1        R_package
influenceR              0.1.0        R_package
inline                  0.3.14       R_package
interactiveDisplayBase  1.10.3       R_package
intergraph              2.0-2        R_package
ipred                   0.9-5        R_package
IRanges                 2.6.1        R_package
irlba                   2.0.0        R_package
Iso                     0.0-17       R_package
iterators               1.0.8        R_package
itertools               0.1-3        R_package
its                     1.1.8        R_package
jpeg                    0.1-8        R_package
jsonlite                1.0          R_package
kernlab                 0.9-24       R_package
KernSmooth              2.23-15      R_package
KFAS                    1.2.3        R_package
kknn                    1.3.1        R_package
klaR                    0.6-12       R_package
km.ci                   0.5-2        R_package
KMsurv                  0.1-5        R_package
knitr                   1.13         R_package
kohonen                 2.0.19       R_package
labeling                0.3          R_package
LaF                     0.6.2        R_package
Lahman                  4.0-1        R_package
lambda.r                1.1.9        R_package
lars                    1.2          R_package
lattice                 0.20-33      R_package
latticeExtra            0.6-28       R_package
lava                    1.4.3        R_package
lazyeval                0.2.0        R_package
LearnBayes              2.15         R_package
lfda                    1.1.1        R_package
lfe                     2.5-1998     R_package
limma                   3.28.16      R_package
linprog                 0.9-2        R_package
lintr                   1.0.0        R_package
lme4                    1.1-12       R_package
lmodel2                 1.7-2        R_package
lmtest                  0.9-34       R_package
locfit                  1.5-9.1      R_package
loo                     0.1.6        R_package
lpSolve                 5.6.13       R_package
lsmeans                 2.23-5       R_package
lubridate               1.5.6        R_package
magrittr                1.5          R_package
mail                    1.0          R_package
mapdata                 2.2-6        R_package
mapproj                 1.2-4        R_package
maps                    3.1.0        R_package
maptools                0.8-39       R_package
markdown                0.7.7        R_package
MASS                    7.3-45       R_package
Matrix                  1.2-6        R_package
MatrixModels            0.4-1        R_package
matrixStats             0.50.2       R_package
maxLik                  1.3-4        R_package
mboost                  2.6-0        R_package
mclust                  5.2          R_package
mcmc                    0.9-4        R_package
mda                     0.4-8        R_package
memoise                 1.0.0        R_package
MEMSS                   0.9-2        R_package
MetABEL                 0.2-0        R_package
methods                 3.3.1        R_package
mgcv                    1.8-12       R_package
mhsmm                   0.4.14       R_package
mice                    2.25         R_package
microbenchmark          1.4-2.1      R_package
mime                    0.5          R_package
miniUI                  0.1.1        R_package
minqa                   1.2.4        R_package
misc3d                  0.8-4        R_package
miscTools               0.6-16       R_package
missForest              1.4          R_package
mix                     1.0-9        R_package
mlbench                 2.1-1        R_package
mlmRev                  1.0-6        R_package
mlogit                  0.2-4        R_package
mnormt                  1.5-4        R_package
modeltools              0.2-21       R_package
mondate                 0.10.01.02   R_package
mstate                  0.2.9        R_package
MSwM                    1.2          R_package
muhaz                   1.2.6        R_package
multcomp                1.4-6        R_package
multicool               0.1-9        R_package
multtest                2.28.0       R_package
munsell                 0.4.3        R_package
mvtnorm                 1.0-5        R_package
network                 1.13.0       R_package
networkDynamic          0.9.0        R_package
nlme                    3.1-128      R_package
nloptr                  1.0.4        R_package
NMF                     0.20.6       R_package
NMFN                    2.0          R_package
nnet                    7.3-12       R_package
nnls                    1.4          R_package
numDeriv                2014.2-1     R_package
nws                     1.7.0.1      R_package
nycflights13            0.2.0        R_package
oligo                   1.36.1       R_package
oligoClasses            1.34.0       R_package
openssl                 0.9.4        R_package
openxlsx                3.0.0        R_package
optextras               2013-10.28   R_package
OptimalCutpoints        1.1-3        R_package
optimx                  2013.8.7     R_package
orcutt                  1.1          R_package
ordinal                 2015.6-28    R_package
packrat                 0.4.7-1      R_package
pamr                    1.55         R_package
parallel                3.3.1        R_package
party                   1.0-25       R_package
partykit                1.1-0        R_package
pbapply                 1.2-1        R_package
pbkrtest                0.4-6        R_package
permute                 0.9-0        R_package
pheatmap                1.0.8        R_package
pkgKitten               0.1.3        R_package
pkgmaker                0.22         R_package
PKI                     0.1-3        R_package
PKPDmodels              0.3.2        R_package
plm                     1.5-12       R_package
plotly                  3.6.0        R_package
plotmo                  3.1.4        R_package
plotrix                 3.6-2        R_package
pls                     2.5-0        R_package
plyr                    1.8.4        R_package
png                     0.1-7        R_package
PoissonSeq              1.1.2        R_package
poLCA                   1.4.1        R_package
polspline               1.1.12       R_package
praise                  1.0.0        R_package
preprocessCore          1.34.0       R_package
pROC                    1.8          R_package
prodlim                 1.5.7        R_package
profileModel            0.5-9        R_package
proto                   0.3-10       R_package
proxy                   0.4-16       R_package
psych                   1.6.6        R_package
purrr                   0.2.2        R_package
pwr                     1.1-3        R_package
quadprog                1.5-5        R_package
quantreg                5.26         R_package
R2HTML                  2.3.2        R_package
R2OpenBUGS              3.2-3.1      R_package
R6                      2.1.2        R_package
randomForest            4.6-12       R_package
ranger                  0.5.0        R_package
RankAggreg              0.5          R_package
RANN                    2.5          R_package
rARPACK                 0.11-0       R_package
rbenchmark              1.0.0        R_package
RBGL                    1.48.1       R_package
R.cache                 0.12.0       R_package
Rcgmin                  2013-2.21    R_package
RColorBrewer            1.1-2        R_package
Rcpp                    0.12.6       R_package
RcppArmadillo           0.7.100.3.1  R_package
RcppEigen               0.3.2.8.1    R_package
RCurl                   1.95-4.8     R_package
regioneR                1.4.2        R_package
registry                0.3          R_package
rem                     1.1.2        R_package
reporttools             1.1.2        R_package
reshape                 0.8.5        R_package
reshape2                1.4.1        R_package
rex                     1.1.1        R_package
rFerns                  2.0.1        R_package
RH2                     0.2.3        R_package
rhdf5                   2.16.0       R_package
Rhtslib                 1.4.3        R_package
rJava                   0.9-8        R_package
RJDBC                   0.2-5        R_package
rjson                   0.2.15       R_package
RJSONIO                 1.3-0        R_package
rlecuyer                0.3-4        R_package
rmarkdown               1.0          R_package
rmeta                   2.16         R_package
R.methodsS3             1.7.1        R_package
rms                     4.5-0        R_package
RMySQL                  0.10.9       R_package
rngtools                1.2.4        R_package
robustbase              0.92-6       R_package
ROCR                    1.0-7        R_package
R.oo                    1.20.0       R_package
roxygen2                5.0.1        R_package
rpart                   4.1-10       R_package
rPython                 0.0-6        R_package
R.rsp                   0.30.0       R_package
Rsamtools               1.24.0       R_package
RSclient                0.7-3        R_package
rsconnect               0.4.3        R_package
Rserve                  1.7-3        R_package
RSiena                  1.1-232      R_package
RSpectra                0.12-0       R_package
RSQLite                 1.0.0        R_package
rstan                   2.10.1       R_package
rstanarm                2.10.1       R_package
rstudioapi              0.6          R_package
rtracklayer             1.32.1       R_package
RUnit                   0.4.31       R_package
R.utils                 2.3.0        R_package
rversions               1.0.2        R_package
Rvmmin                  2013-11.12   R_package
S4Vectors               0.10.2       R_package
sandwich                2.3-4        R_package
scagnostics             0.2-4        R_package
scales                  0.4.0        R_package
scalreg                 1.0          R_package
scatterplot3d           0.3-37       R_package
segmented               0.5-1.4      R_package
sendmailR               1.2-1        R_package
seqinr                  3.3-0        R_package
seqminer                5.3          R_package
setRNG                  2013.9-1     R_package
setwidth                1.0-4        R_package
shapefiles              0.7          R_package
shiny                   0.13.2       R_package
shinyjs                 0.6          R_package
shinystan               2.2.0        R_package
shinythemes             1.0.1        R_package
SKAT                    1.2.1        R_package
sna                     2.3-2        R_package
snow                    0.4-1        R_package
sp                      1.2-3        R_package
spam                    1.3-0        R_package
SparseM                 1.7          R_package
spatial                 7.3-11       R_package
spdep                   0.6-5        R_package
speedglm                0.3-1        R_package
sphet                   1.6          R_package
splines                 3.3.1        R_package
splm                    1.3-7        R_package
spls                    2.2-1        R_package
sqldf                   0.4-10       R_package
stabledist              0.7-0        R_package
stabs                   0.5-1        R_package
StanHeaders             2.10.0-2     R_package
statmod                 1.4.24       R_package
statnet                 2016.4       R_package
statnet.common          3.3.0        R_package
stats                   3.3.1        R_package
stats4                  3.3.1        R_package
stringdist              0.9.4.1      R_package
stringi                 1.1.1        R_package
stringr                 1.0.0        R_package
strucchange             1.5-1        R_package
subselect               0.12-5       R_package
SummarizedExperiment    1.2.3        R_package
superpc                 1.09         R_package
SuppDists               1.1-9.2      R_package
survcomp                1.22.0       R_package
survival                2.39-4       R_package
survival                2.39-5       R_package
survivalROC             1.0.3        R_package
survMisc                0.5.3        R_package
svMisc                  0.9-70       R_package
svTools                 0.9-4        R_package
svUnit                  0.7-12       R_package
synchronicity           1.1.9.1      R_package
tables                  0.7.79       R_package
tcltk                   3.3.1        R_package
TeachingDemos           2.10         R_package
tergm                   3.4.0        R_package
testit                  0.5          R_package
testthat                1.0.2        R_package
texreg                  1.36.7       R_package
TH.data                 1.0-7        R_package
threejs                 0.2.2        R_package
tibble                  1.1          R_package
tidyr                   0.5.1        R_package
tikzDevice              0.10-1       R_package
timeDate                3012.100     R_package
timeSeries              3022.101.2   R_package
tis                     1.30         R_package
tnam                    1.6.2        R_package
tools                   3.3.1        R_package
tripack                 1.3-7        R_package
trust                   0.1-7        R_package
tseries                 0.10-35      R_package
ucminf                  1.1-3        R_package
urca                    1.2-9        R_package
utils                   3.3.1        R_package
VariABEL                0.9-2.1      R_package
VariantAnnotation       1.18.5       R_package
vars                    1.5-2        R_package
vcd                     1.4-1        R_package
vegan                   2.4-0        R_package
VennDiagram             1.6.17       R_package
VGAM                    1.0-2        R_package
viridis                 0.3.4        R_package
visNetwork              1.0.1        R_package
wbstats                 0.1          R_package
webshot                 0.3.2        R_package
whisker                 0.3-2        R_package
withr                   1.0.2        R_package
xergm                   1.7.3        R_package
xergm.common            1.7.4        R_package
xgboost                 0.4-4        R_package
xml2                    1.0.0        R_package
XML                     3.98-1.4     R_package
xtable                  1.8-2        R_package
xts                     0.9-7        R_package
XVector                 0.12.0       R_package
yaml                    2.1.13       R_package
zlibbioc                1.18.0       R_package
zoo                     1.7-13       R_package
```

# Do I have to manually build a container each time I want to do this?
No! But I do recommend the interactive method when trying to install software for the first time, since it'll be easier to troubleshoot any issues. But once you have a reliable installation recipe ready, checkout [the documentation](http://singularity.lbl.gov/#bootstrap) for how to create a definition file for bootstrapping ready-to-go containers!

# Contributions
Thank you to everyone who has contributed!
- [gmkurtzer](https://github.com/gmkurtzer): building Singularity, and feedback on examples and document contents
- [dwaggott](https://github.com/dwaggott): providing a great list of R packages, and an example of how to install R packages from the command line via an R script
- [kprybol](https://github.com/kprybol): another extensive list of useful R packages, and an explanation of R's package system (and it's use of Suggests, Depends, Imports, LinkingTo, and Enhances)
