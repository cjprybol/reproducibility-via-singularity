wget --no-check-certificate https://stanfordmedicine.box.com/shared/static/qf02pfpv1e87sxz2g1yssufnpe0szbo3.gz -O reads_1.fq.gz
wget --no-check-certificate https://stanfordmedicine.box.com/shared/static/xu1snsf8st8sbx3x1cejc2f4ltkn0qsc.gz -O reads_2.fq.gz
wget http://bio.math.berkeley.edu/kallisto/transcriptomes/Saccharomyces_cerevisiae.R64-1-1.rel81.cdna.all.fa.gz
kallisto index --make-unique --index=index Saccharomyces_cerevisiae.R64-1-1.rel81.cdna.all.fa.gz
kallisto quant --index=index --output-dir=Kallisto --threads=1 --plaintext reads_1.fq.gz reads_2.fq.gz
