singularity exec demo.img wget --no-check-certificate https://stanfordmedicine.box.com/shared/static/7ibiqg9j0xejkxjd25tvu16jb4ko1szl.gz -O reads.1.fq.gz
singularity exec demo.img wget --no-check-certificate https://stanfordmedicine.box.com/shared/static/gnyhj0zbc168emwk1bjzzmhodcnawn65.gz -O reads.2.fq.gz
singularity exec demo.img wget http://bio.math.berkeley.edu/kallisto/transcriptomes/Homo_sapiens.GRCh38.rel79.cdna.all.fa.gz
singularity exec demo.img kallisto index --make-unique --index=index Homo_sapiens.GRCh38.rel79.cdna.all.fa.gz
singularity exec demo.img kallisto quant --index=index --output-dir=Kallisto --threads=1 --plaintext reads.1.fq.gz reads.2.fq.gz
