singularity exec demo.img wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR308/003/SRR3084933/SRR3084933_1.fastq.gz
singularity exec demo.img wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR308/003/SRR3084933/SRR3084933_2.fastq.gz
singularity exec demo.img wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.transcripts.fa.gz
singularity exec demo.img kallisto index --make-unique --index=$(pwd)/index gencode.v24.transcripts.fa.gz
singularity exec demo.img kallisto quant --index=$(pwd)/index --output-dir=$(pwd)/Kallisto --threads=1 --plaintext SRR3084933_1.fastq.gz SRR3084933_2.fastq.gz
