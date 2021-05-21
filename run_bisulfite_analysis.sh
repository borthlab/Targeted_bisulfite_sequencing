#---------------------------------------------------------
### Prepare reference genome ###
#---------------------------------------------------------
/usr/local/bioinf/bismark/bismark_v0.22.1/bismark_genome_preparation  $ref

### Preprocessing and alignment of raw reads to extract methylation levels
./nick_bisulfite_analysis.sh -sample 97036 -ref CMV
./nick_bisulfite_analysis.sh -sample 97037 -ref CMV
./nick_bisulfite_analysis.sh -sample 97038 -ref CMV
./nick_bisulfite_analysis.sh -sample 97039 -ref Fut8
./nick_bisulfite_analysis.sh -sample 97040 -ref Fut8
./nick_bisulfite_analysis.sh -sample 97041 -ref Fut8

