# Make your research more reproducible with [Singularity](http://singularity.lbl.gov)
Services like [GitHub](https://github.com/), [Bitbucket](https://bitbucket.org), and [GitLab](https://about.gitlab.com/) have democratized access to affordable (or free!) tools that reinforce reproducibility in research, via the code repositories that these services offer. These services make it easier to backup, version control, collaborate on, and distribute code (or any other text-based file). These features make it easier for researchers to write and maintain high-quality code. These features also increase the chance that someone else will review the code for errors or bugs. These reviews can be done by some combination of reading and executing the code. Unfortunately, unless the code is run on the exact same computer, while logged in as the same user, getting the same code to run the same way can be a research project in and of itself.


# Why Singularity?
Universities and research laboratories often conduct their work on shared HPC clusters. These clusters are all set up on slightly different configurations of hardware and operating system software. Trying to recreate an environment to re-run code exactly as it was executed on another cluster, an environment where all of the C pointers align in the same way and the Java versions sync up, is beyond the abilities, free time, and account privileges of most researchers who may want to try to reproduce your results. Singularity is an implementation of a "container", and an engine to run those containers, that allows researchers to isolate the environment needed to produce a result away from the resources and code that is already available. This means that your colleague at University X can run the analysis exactly the same way on their cluster as you are running it on your cluster at University Y, and all it requires is sharing a git repository and a container image.

# TL;DR

**Hardware requirements for this example**
  - RAM ~ 4Gb
  - CPU = 1
  - available disk space ~ 1 Gb

Run this code on a linux-based operating system with singularity installed and in the `$PATH`
```bash
git clone https://github.com/cjprybol/reproducibility-via-singularity.git
cd reproducibility-via-singularity
wget --no-check-certificate https://stanfordmedicine.box.com/shared/static/fbuap57hb3ywl7f741t1x4mcq9akznl4.img -O demo.img
./reproduce_project.sh
```

time to run < 5 minutes on a modern laptop with decent download rates (10 Mb/s?)

# What does that code do?

I encourage you to read and execute the `reproduce_project.sh` file and see for yourself! But in summary, it will download 1 sample of paired-end RNA sequencing reads, and quantify transcript isoform abundances in the sample. Run `less Kallisto/abundance.tsv` to see the results!

This is a contrived example, designed to run quickly and efficiently to demonstrate a point that reproducible research can be, and should be, be as simple as running 4 code instructions. However, by extending the software libraries in the container, and extending the analytical complexity in the analysis scripts, these same 4 lines of code can be sufficient to reproduce entire thesis projects.

Requiring little more effort for others than the wait required to perform the computation, this simplicity of reproducibility lowers the barrier for others to investigate your work. This encourages responsible research conduct and increases the rate of knowledge transfer.

# What is Singularity?

A Singularity container is a complete linux kernel capable of executing code and running software. Because the container is a singe file, it can be easily moved to another computer with Singularity and run in exactly the same way. The Singularity engine on each cluster handles the translation to the unique hardware and OS configurations of each cluster.

# What is Singularity not?

A magic, drag and drop, point and click solution. Singularity requires some understanding of how operating systems work in order to effectively use it and understand it. That being said, it's also very well designed, and should be intuitive for anyone who has set up a computing environment on their personal computer, a cluster, a VM, or a cloud service. Singularity will generate an operating system inside of the image file, or "container". Singularity will ensure that whatever you install in the container will execute the same on any computer with Singularity installed. You simply load the container with software and configure it as you would set up any other basic linux installation.

Containers require sudo/root/admin privileges to modify and edit. This means that when creating containers, users must take care to note whether or not the software installed inside of the container will try to download or edit files that exist inside of the directory structure of the container during use. If a container is shipped without those files pre-downloaded and pre-initialized, end users without sudo permission will experience errors as the software fails to write to disk. For most situations, the strategy is simply to do the work of preloading the data! But for some use cases, such as when the software works by first downloading large databases that can range anywhere from several dozen Gb to several dozen Tb, preloading that much information and sharing it via containers may be an impractical strategy. Many pieces of software that work in this way are configurable to write, edit, and save files to a user-specified directory. If that is the case, you should be able to get around this issue by configuring the software to use a directory outside of the container.

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

Because containers interact so closely with the host system, installing Singularity, and building Singularity images requires root/sudo/admin level permissions. While many users of shared computing clusters are unlikely to have these permission levels on the clusters where they intend to run the analysis, they may have those permissions on whatever local computer (laptop/desktop) they use to ssh into the shared computing cluster. If you happen to have linux-based OS with sudo permissions, go ahead and jump to the section [Start here if linux user with root/sudo/admin privileges](https://github.com/cjprybol/reproducibility-via-singularity/blob/master/README.md#start-here-if-linux-user-with-rootsudoadmin-privileges). If you have a Windows OS or an Apple OS (with or without root/sudo/admin), or a linux OS without root/sudo/admin, you'll need to get a linux virtual machine running where you do have root/sudo/admin permissions. If you're already familiar with the process of setting up virtual machines, go ahead with whatever method you know. If installing a virtual machine is a new thing for you, I recommend checking out [VirtualBox](https://www.virtualbox.org/). If you don't have root access on your local machine, ask a system admin for assistance.

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
  vb.memory = "4096"
end
```

Restart the virtual machine and ssh into it
```bash
vagrant up && vagrant ssh
```
You are now a linux user with root/sudo/admin privileges


# Start here if linux user with root/sudo/admin privileges

Here we will install Singularity starting from an Ubuntu 14.04 LTS "Trusty" 64-bit base installation
```bash
sudo apt-get install -y build-essential git vim autoconf libtool curl debootstrap
git clone https://github.com/gmkurtzer/singularity.git && cd singularity && ./autogen.sh && ./configure --prefix=/usr/local && make && sudo make install
```

We want to create a singularity container, load it with an operating system, and install and configure the software necessary to run our analysis onto it.

First, we need to allocate the file. Here `--size` represents the maximum size (in Mb) that the container is allowed to take on. This container is allowed to take on a maximum of 15Gb. **Interesting aside** Containers are initialized as sparse images. If you evaluate the allocated space for the image with `ls -lah my_container.img`, and compare that disk size to what is returned by `du -h my_container.img`, you'll see that the files are only keeping track of the informative content installed, rather than the total possible disk space they could use.
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

[Linuxbrew](http://linuxbrew.sh/) and [Anaconda](https://www.continuum.io/downloads) are two great package managers that greatly simplify the process of installing and managing software. The default software available via apt-get is often out of date compared to those available via either Linuxbrew and Anaconda, and additionally, both of these package managers have recipes to install a wide-variety of software useful for researchers.

By default, when entering a container, you enter into it and stay in the same directory on the host where you started from. Run `pwd` and notice the `(unreachable)` pre-pension to the path. Our current directory is immutable because we have both remained in a directory on the host computer, but also specified that we want to `--contain` our actions to only the contents of the container.

Let's change to the user directory of our current user inside of the container, which is `root`, so we have a mutable space that can generate temporary files needed to install software.
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

Now let's use our package managers to quickly and easily install and configure software into the container. You can review their full offerings here -> [Homebrew-science](https://github.com/Homebrew/homebrew-science) & [Bioconda channel of anaconda](https://github.com/bioconda/bioconda-recipes/tree/master/recipes)

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

I'll install [Julia](http://julialang.org/) from source via GitHub, as an example of how to manually install software not available via the package managers, as well as a plug for the language (which I recommend you try out!). Add your own custom recipes for installing software here and configure a system that meets your needs.
```bash
git clone git://github.com/JuliaLang/julia.git && cd julia && git checkout release-0.4 && touch Make.user && echo "USE_SYSTEM_GMP=1" >> Make.user && echo "USE_SYSTEM_MPFR=1" >> Make.user && make
ln -s /Software/julia/julia /usr/local/bin
```

Install [RTG core](http://realtimegenomics.com/products/rtg-core/). **NOTE** This software is license restricted. It's free for non-commercial academic use, but if you intend to use it commercially you'll have to buy a license (alternatively, just skip this installation).
```bash
cd /Software && wget --no-check-certificate https://github.com/RealTimeGenomics/rtg-core/releases/download/3.6.2/rtg-core-non-commercial-3.6.2-linux-x64.zip
unzip rtg-core-non-commercial-3.6.2-linux-x64.zip && rm rtg-core-non-commercial-3.6.2-linux-x64.zip
ln -s /Software/rtg-core-non-commercial-3.6.2/rtg /usr/local/bin
```

The first time you run RTG, it will ask whether or not it can perform logging to help the developers improve the software. Because we intend to run the container without sudo and not in `--writable` mode, any attempts RTG makes to save log files to disk will probably fail, so just say no.
```bash
Singularity.test.img> rtg
RTG has a facility to automatically send basic usage information to Real
Time Genomics. This does not contain confidential information such as
command-line parameters or dataset contents.

Would you like to enable automatic usage logging (y/n)? n
```

We've got our system fully loaded with the software we want, but our `$PATH` update was only for this session. We'll need to make our `$PATH` updates permanent to make the software installed inside of the `/Software` directory available by name alone. Alternatively, you can also specify the full path when calling executables inside of the container. In Singularity version >= 2.1, you can update the `$PATH` by modifying the `/environment` file, which is loaded each time you interact with the container.
```bash
cd / && rm /environment && wget --no-check-certificate https://raw.githubusercontent.com/cjprybol/reproducibility-via-singularity/master/environment
```

We'll download one more pre-written script that will list all software installed with Linuxbrew and Anaconda, as well as the Julia version. It will sort the list and return to us an A-Z list of installed software with the version number for everything! Go [here](https://github.com/cjprybol/reproducibility-via-singularity/blob/master/README.md#how-do-i-get-the-version-numbers-of-installed-software) to see how to call it
```bash
wget --no-check-certificate https://raw.githubusercontent.com/cjprybol/reproducibility-via-singularity/master/singularity
chmod 775 singularity
```

Exit the container to the host linux
```bash
exit
```

Another tool I use that cannot be installed via either package manager due to licensing restrictions is [Genome Analysis Toolkit](https://www.broadinstitute.org/gatk/). I've downloaded the GATK installer to the same directory where my container is. By entering the container again without the `--contain` flag, we allow the container to interact with the host system, and copy GATK into Anaconda
```bash
sudo singularity shell --writable test.img
gatk-register /home/vagrant/GenomeAnalysisTK-3.6.tar.bz2
rm GenomeAnalysisTK-3.6.tar.bz2
exit
```

You're all done, you've built a great base-image for computational genomics! Adjust these installation steps to your needs.

# How do I get the version numbers of installed software?
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
```

# What about R and Julia Packages?

**R**
```bash
singularity exec test.img R
installed.packages(fields = c("Package", "Version"))[,c("Package", "Version")]
quit(save="no")
```

**Julia**
```bash
singularity exec test.img julia
Pkg.installed()
exit()
```

# Do I have to manually build a container each time I want to do this?
No! But I do recommend the interactive method when trying to install software for the first time, since it'll be easier to troubleshoot any issues. But once you have a reliable installation recipe ready, checkout [the documentation](http://singularity.lbl.gov/#bootstrap) for how to create a definition file for bootstrapping ready-to-go containers!

# Contributions
Thank you to everyone who has contributed!
- @gmkurtzer: building Singularity, and feedback on examples and document contents
- @dwaggott: providing a great list of R packages, and an example of how to install R packages from the command line via an R script
- @kprybol: another extensive list of useful R packages, and an explanation of R's package system (and it's use of Suggests, Depends, Imports, LinkingTo, and Enhances)
