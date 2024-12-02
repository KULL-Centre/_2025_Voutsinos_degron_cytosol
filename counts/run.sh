module load eautils

ls ../fastq_230313/*R1_001.fastq.gz  |  parallel -P 64 ./call_zerotol_paired.sh {};
ls ../fastq_230327/*R1_001.fastq.gz  |  parallel -P 64 ./call_zerotol_paired.sh {};
ls ../fastq_230406/*R1_001.fastq.gz  |  parallel -P 64 ./call_zerotol_paired.sh {};
ls ../fastq_221202/*R1_001.fastq.gz  |  parallel -P 64 ./call_zerotol_paired.sh {};

grep LABEL nohup.out | head -n1 > paring_stat.txt
grep RESULT nohup.out >> paring_stat.txt
