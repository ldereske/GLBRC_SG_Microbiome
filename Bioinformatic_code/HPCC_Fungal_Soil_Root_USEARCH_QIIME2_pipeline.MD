## Let  extract the reads from 20200302_Amplicon_PE300 
### GLBRC_Lib11_ITS_Soil_2018
```
cd /mnt/ufs18/rs-001/EvansLab/MMPRNT_2018/20200302_Amplicon_PE300
gunzip *.fastq.gz
```



####Let's demultiplex this dataset 
I created a bar code file in a text editor
```
cd /mnt/research/EvansLab/MMPRNT_2018/20200302_Amplicon_PE300/ 
mkdir demux_files

```
# Now let's try the QIIME2 versions
https://forum.qiime2.org/t/demultiplexing-and-trimming-adapters-from-reads-with-q2-cutadapt/2313

###### Import the multiplexed sequences

```
mkdir QIIME2_demux
cd QIIME2_demux/
mkdir compressed_seqs
cd ..

gzip GLBRC_Lib11_ITS_Evans_soil_2018_S1_L001_R1_001.fastq  QIIME2_demux/compressed_seqs/forward.fastq.gz
gzip GLBRC_Lib11_ITS_Evans_soil_2018_S1_L001_R2_001.fastq  QIIME2_demux/compressed_seqs/reverse.fastq.gz
gzip GLBRC_Lib11_ITS_Evans_soil_2018_S1_L001_I1_001.fastq  QIIME2_demux/compressed_seqs/barcodes.fastq.gz


export PATH="/mnt/research/EvansLab/Software/anaconda3/bin:$PATH"
source activate qiime2-2019.4
qiime tools import --show-importable-types
cd QIIME2_demux
qiime tools import --type EMPPairedEndSequences --input-path compressed_seqs --output-path multiplexed-seqs_paired_end.qza
#Imported compressed_seqs as EMPPairedEndDirFmt to multiplexed-seqs_paired_end.qza

nano Q2_demux_lib_MMPRNT018_soil_fungi.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=48:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name Q2_demux_lib_MMPRNT018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/MMPRNT_2018/20200302_Amplicon_PE300/QIIME2_demux
export PATH="/mnt/research/EvansLab/Software/anaconda3/bin:$PATH"
source activate qiime2-2019.4



qiime demux emp-paired --m-barcodes-file MMPRNT_fungal_QIIME_mapping_for.txt --m-barcodes-column BarcodeSequence --p-rev-comp-mapping-barcodes --p-no-golay-error-correction --i-seqs multiplexed-seqs_paired_end.qza --o-per-sample-sequences demultiplexed-seqs.qza --o-error-correction-details demux-details.qza
###End

sbatch Q2_demux_lib_MMPRNT018_soil_fungi.sbatch
#Submitted batch job 57000649
```

Let's try to visualize the demux sequences
```
qiime demux summarize  --i-data demultiplexed-seqs.qza  --o-visualization demux-seq_summary.qzv
#Saved Visualization to: demux-seq_summary.qzv


```

---
summary
Minimum:	87
Median:	143104.5
Mean:	140706.50520833334
Maximum:	320018
Total:	27015649
---

Let's export the demux files
```
qiime tools export  --input-path demultiplexed-seqs.qza  --output-path output/demultiplexed-seqs
#Exported demultiplexed-seqs.qza as SingleLanePerSamplePairedEndFastqDirFmt to directory output/demultiplexed-seqs
```
Let's unzip the resulting demultiplex files
```
cd /mnt/research/EvansLab/MMPRNT_2018/20200302_Amplicon_PE300/QIIME2_demux/output/demultiplexed-seqs/
gunzip *.fastq.gz
```

## 1) Quality checking
### 1a) First look at the quality of  raw unmerged seqs for run1
#https://www.drive5.com/usearch/manual/pipe_readprep_understand.html

```
mkdir fastq_info
```

#### make a forloop to run fastx_info on every file
```
nano fasta_info_fq.sh
#!/bin/bash --login
for fq in *.fastq
do
usearch10.0.240 -fastx_info $fq -output fastq_info/$fq
done
```

#### make file executable and run the for loop create in your `fasta_info_fq.sh` file
```
chmod +x fasta_info_fq.sh
./fasta_info_fq.sh
```

#### move to the fastq_info directory
```
cd fastq_info/
```

#### Now run the below code to summarize the fastq info for all of the forward and reverse reads


```

grep "^File" * > 20200302_Amplicon_PE300_fastq_lengths.txt
grep "^EE" * > 20200302_Amplicon_PE300_fastq_EE.txt

```

Look for any forward and reverse reads that look especially bad in terms of quality (high E is bad quality). This info will also be really helpful for troubleshooting later on (e.g. why some samples have extremely low read numbers)


##  2) Merge the forward and reverse sequences and trim adapters (for each run individually)
### 2a) Merge Pairs
#### Make sure you are in the folder with the extracted forward and reverse reads from Run 1
https://www.drive5.com/usearch/manual/merge_options.html
#### -alnout gives you a human readable text file of the alignment and misalignments for each pair merged.
#### -tabbedout give you extensive information on the quality of the merge
#### This step takes approximately 10 minutes

```

cd /mnt/research/EvansLab/Lukas/MMPRNT018/
#mkdir Fungal_libs
mkdir Fungal_libs/20200302_Amplicon_PE300

#Crap this was killed the first time I ran it so now I am going to submit a job

nano merg_seq_lib11_MMPRNT018_soil_fungi.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=3:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name merge_lib11_MMPRNT018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/MMPRNT_2018/20200302_Amplicon_PE300/QIIME2_demux/output/demultiplexed-seqs/

usearch11.0.667 -fastq_mergepairs *_R1*.fastq -relabel @ -fastqout /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200302_Amplicon_PE300/lib11_20200302_Amplicon_PE300_ITS1.fastq  -tabbedout /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200302_Amplicon_PE300/lib11_20200302_Amplicon_PE300_ITS1_pair_report.txt -alnout /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200302_Amplicon_PE300/lib11_20200302_Amplicon_PE300_ITS1_pair_aln.txt

###End
sbatch merg_seq_lib11_MMPRNT018_soil_fungi.sbatch
#Submitted batch job 58239587

```

Merge stats

---

#Run time 00:24:24
Totals:
  27015649  Pairs (27.0M)
  24135681  Merged (24.1M, 89.34%)
  17068562  Alignments with zero diffs (63.18%)
   2576142  Too many diffs (> 5) (9.54%)
    303826  No alignment found (1.12%)
         0  Alignment too short (< 16) (0.00%)
  13430327  Staggered pairs (49.71%) merged & trimmed
    263.53  Mean alignment length
    299.96  Mean merged length
      0.59  Mean fwd expected errors
      1.28  Mean rev expected errors
      0.05  Mean merged expected errors

---


###Check sample names to make sure all of the samples are present


```

usearch11.0.667 -fastx_get_sample_names /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200302_Amplicon_PE300/lib11_20200302_Amplicon_PE300_ITS1.fastq -output /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200302_Amplicon_PE300/lib11_20200302_Amplicon_PE300_ITS1_samples.txt

```

---
#Should be 192 samples

01:33 38Mb    100.0% 192 samples found
---


## Let  demultiplex the reads from 20200306_M02808_Amplicon_PE300 
### GLBRC_Lib9_ITS_Evans_soil_2018

```
cd /mnt/ufs18/rs-001/EvansLab/MMPRNT_2018/20200306_M02808_Amplicon_PE300

```
I created a bar code file in a text editor

### Let's use the QIIME2 demux
https://docs.qiime2.org/2020.2/plugins/available/demux/emp-paired/?highlight=demux

###### Import the multiplexed sequences
We need to rename and move the zipped sequence files
```
mkdir QIIME2_demux
cd QIIME2_demux/
mkdir compressed_seqs
cd ..

cp GLBRC_Lib9_ITS_Evans_soil_2018_S1_L001_R1_001.fastq.gz  QIIME2_demux/compressed_seqs/forward.fastq.gz
cp GLBRC_Lib9_ITS_Evans_soil_2018_S1_L001_R2_001.fastq.gz  QIIME2_demux/compressed_seqs/reverse.fastq.gz
cp GLBRC_Lib9_ITS_Evans_soil_2018_S1_L001_I1_001.fastq.gz  QIIME2_demux/compressed_seqs/barcodes.fastq.gz

#Now we need to convert the files to Q2 format
export PATH="/mnt/research/EvansLab/Software/anaconda3/bin:$PATH"
source activate qiime2-2019.4
#qiime tools import --show-importable-types
cd QIIME2_demux
qiime tools import --type EMPPairedEndSequences --input-path compressed_seqs --output-path multiplexed-seqs_paired_end.qza
#Imported compressed_seqs as EMPPairedEndDirFmt to multiplexed-seqs_paired_end.qza

nano Q2_demux_lib9_MMPRNT018_soil_fungi.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=48:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name Q2_demux_lib9_MMPRNT018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/MMPRNT_2018/20200306_M02808_Amplicon_PE300/QIIME2_demux
export PATH="/mnt/research/EvansLab/Software/anaconda3/bin:$PATH"
source activate qiime2-2019.4



qiime demux emp-paired --m-barcodes-file MMPRNT_QIIME_mapping_lib9.txt --m-barcodes-column BarcodeSequence --p-rev-comp-mapping-barcodes --p-no-golay-error-correction --i-seqs multiplexed-seqs_paired_end.qza --o-per-sample-sequences demultiplexed-seqs.qza --o-error-correction-details demux-details.qza

#Let's also export the demux files

qiime tools export  --input-path demultiplexed-seqs.qza  --output-path output/demultiplexed-seqs

###End

sbatch Q2_demux_lib9_MMPRNT018_soil_fungi.sbatch
#Run time 04:35:58
#Submitted batch job 58236758
#Saved SampleData[PairedEndSequencesWithQuality] to: demultiplexed-seqs.qza
#Saved ErrorCorrectionDetails to: demux-details.qza
#Exported demultiplexed-seqs.qza as SingleLanePerSamplePairedEndFastqDirFmt to directory output/demultiplexed-seqs
```

Let's try to visualize the demux sequences
```
export PATH="/mnt/research/EvansLab/Software/anaconda3/bin:$PATH"
source activate qiime2-2019.4
qiime demux summarize  --i-data demultiplexed-seqs.qza  --o-visualization demux-seq_summary.qzv
#Saved Visualization to: demux-seq_summary.qzv
#qiime tools view demux-seq_summary.qzv did not work well on the hpcc
#I use https://view.qiime2.org/ to visualize the files since I do not have QIIME2 on my computer 

```

---

#summary
Minimum:	1
Median:	112436.0
Mean:	116219.32984293194
Maximum:	402029
Total:	22197892
---

Let's unzip the resulting demultiplex files
```
cd /mnt/research/EvansLab/MMPRNT_2018/20200306_M02808_Amplicon_PE300/QIIME2_demux/output/demultiplexed-seqs/
gunzip *.fastq.gz
```

## 1) Quality checking
### 1a) First look at the quality of  raw unmerged seqs for run1
#https://www.drive5.com/usearch/manual/pipe_readprep_understand.html

```
mkdir fastq_info
```

#### make a forloop to run fastx_info on every file

```
nano fasta_info_fq.sh
#!/bin/bash --login
for fq in *.fastq
do
usearch10.0.240 -fastx_info $fq -output fastq_info/$fq
done

```

#### make file executable and run the for loop create in your `fasta_info_fq.sh` file
```
chmod +x fasta_info_fq.sh
./fasta_info_fq.sh
```

#### move to the fastq_info directory

```
cd fastq_info/

```

#### Now run the below code to summarize the fastq info for all of the forward and reverse reads

```

grep "^File" * > 20200306_M02808_Amplicon_PE300_fastq_lengths.txt
grep "^EE" * > 20200306_M02808_Amplicon_PE300_fastq_EE.txt

```

Look for any forward and reverse reads that look especially bad in terms of quality (high E is bad quality). This info will also be really helpful for troubleshooting later on (e.g. why some samples have extremely low read numbers)


##  2) Merge the forward and reverse sequences and trim adapters (for each run individually)
### 2a) Merge Pairs
#### Make sure you are in the folder with the extracted forward and reverse reads from Run 1
https://www.drive5.com/usearch/manual/merge_options.html
#### -alnout gives you a human readable text file of the alignment and misalignments for each pair merged.
#### -tabbedout give you extensive information on the quality of the merge
#### This step takes approximately 10 minutes

```

cd /mnt/research/EvansLab/Lukas/MMPRNT018/
#mkdir Fungal_libs
mkdir Fungal_libs/20200306_M02808_Amplicon_PE300
cd Fungal_libs/20200306_M02808_Amplicon_PE300
#I am going to submit a job

nano merg_seq_lib9_MMPRNT018_soil_fungi.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=3:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name merge_lib9_MMPRNT018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/MMPRNT_2018/20200306_M02808_Amplicon_PE300/QIIME2_demux/output/demultiplexed-seqs/

usearch11.0.667 -fastq_mergepairs *_R1*.fastq -relabel @ -fastqout /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200306_M02808_Amplicon_PE300/lib09_20200302_Amplicon_PE300_ITS1.fastq  -tabbedout /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200306_M02808_Amplicon_PE300/lib09_20200302_Amplicon_PE300_ITS1_pair_report.txt 

###End
sbatch merg_seq_lib9_MMPRNT018_soil_fungi.sbatch
#Submitted batch job 58336913

```

Merge stats

---
Totals:
  22197892  Pairs (22.2M)
  18539541  Merged (18.5M, 83.52%)
  12582620  Alignments with zero diffs (56.68%)
   3387923  Too many diffs (> 5) (15.26%)
    270428  No alignment found (1.22%)
         0  Alignment too short (< 16) (0.00%)
  12210769  Staggered pairs (55.01%) merged & trimmed
    265.46  Mean alignment length
    303.18  Mean merged length
      0.60  Mean fwd expected errors
      1.57  Mean rev expected errors
      0.05  Mean merged expected errors

---

###Check sample names to make sure all of the samples are present


```
usearch11.0.667 -fastx_get_sample_names /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200306_M02808_Amplicon_PE300/lib09_20200302_Amplicon_PE300_ITS1.fastq -output /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200306_M02808_Amplicon_PE300/lib09_20200302_Amplicon_PE300_ITS1_samples.txt
```
---
# Should be 192 samples but one sample lost during demux (MMPRNT1151)
01:09 38Mb    100.0% 190 samples found
# Sample MMPRNT1130 was lost in merging do to only having one read post demux.
---

## Let  demultiplex the reads from 20200306_M03127_Amplicon_PE300

### GLBRC_Lib10_ITS_Evans_soil_2018

```
cd /mnt/research/EvansLab/MMPRNT_2018/20200306_M03127_Amplicon_PE300

```
I created a bar code file in a text editor

### Let's use the QIIME2 demux
https://docs.qiime2.org/2020.2/plugins/available/demux/emp-paired/?highlight=demux

###### Import the multiplexed sequences
We need to rename and move the zipped sequence files
```
mkdir QIIME2_demux
cd QIIME2_demux/
mkdir compressed_seqs
cd ..

cp GLBRC_Lib10_ITS_Evans_soil_2018_S1_L001_R1_001.fastq.gz  QIIME2_demux/compressed_seqs/forward.fastq.gz
cp GLBRC_Lib10_ITS_Evans_soil_2018_S1_L001_R2_001.fastq.gz  QIIME2_demux/compressed_seqs/reverse.fastq.gz
cp GLBRC_Lib10_ITS_Evans_soil_2018_S1_L001_I1_001.fastq.gz  QIIME2_demux/compressed_seqs/barcodes.fastq.gz

#Now we need to convert the files to Q2 format and demux 

#qiime tools import --show-importable-types


cd QIIME2_demux

nano Q2_demux_lib10_MMPRNT018_soil_fungi.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=48:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name Q2_demux_lib10_MMPRNT018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/MMPRNT_2018/20200306_M03127_Amplicon_PE300/QIIME2_demux
export PATH="/mnt/research/EvansLab/Software/anaconda3/bin:$PATH"
source activate qiime2-2019.4

qiime tools import --type EMPPairedEndSequences --input-path compressed_seqs --output-path multiplexed-seqs_paired_end.qza

qiime demux emp-paired --m-barcodes-file MMPRNT_QIIME_mapping_lib10.txt --m-barcodes-column BarcodeSequence --p-rev-comp-mapping-barcodes --p-no-golay-error-correction --i-seqs multiplexed-seqs_paired_end.qza --o-per-sample-sequences demultiplexed-seqs.qza --o-error-correction-details demux-details.qza

#Let's also export the demux files

qiime tools export  --input-path demultiplexed-seqs.qza  --output-path output/demultiplexed-seqs

#Let's try to visualize the demux sequences too

qiime demux summarize  --i-data demultiplexed-seqs.qza  --o-visualization demux-seq_summary_lib10.qzv

###End

sbatch Q2_demux_lib10_MMPRNT018_soil_fungi.sbatch
#Submitted batch job 58240804
#Run time 05:04:47
```

Demux stats

---
Minimum:	5
Median:	162595.0
Mean:	163024.96858638743
Maximum:	326518
Total:	31137769

---
Let's unzip the resulting demultiplex files

```
cd  /mnt/research/EvansLab/MMPRNT_2018/20200306_M03127_Amplicon_PE300/QIIME2_demux/output/demultiplexed-seqs/
gunzip *.fastq.gz
```

## 1) Quality checking
### 1a) First look at the quality of  raw unmerged seqs for run1
#https://www.drive5.com/usearch/manual/pipe_readprep_understand.html

```
mkdir fastq_info
```

#### make a forloop to run fastx_info on every file
```
nano fasta_info_fq.sh
#!/bin/bash --login
for fq in *.fastq
do
usearch10.0.240 -fastx_info $fq -output fastq_info/$fq
done
```

#### make file executable and run the for loop create in your `fasta_info_fq.sh` file
```
chmod +x fasta_info_fq.sh
./fasta_info_fq.sh
```

#### move to the fastq_info directory
```
cd fastq_info/
```

#### Now run the below code to summarize the fastq info for all of the forward and reverse reads

```
grep "^File" * > 20200306_M03127_Amplicon_PE300_fastq_lengths.txt
grep "^EE" * > 20200306_M03127_Amplicon_PE300_fastq_EE.txt
```

Look for any forward and reverse reads that look especially bad in terms of quality (high E is bad quality). This info will also be really helpful for troubleshooting later on (e.g. why some samples have extremely low read numbers)


##  2) Merge the forward and reverse sequences and trim adapters (for each run individually)
### 2a) Merge Pairs
#### Make sure you are in the folder with the extracted forward and reverse reads from Run 1
https://www.drive5.com/usearch/manual/merge_options.html
#### -alnout gives you a human readable text file of the alignment and misalignments for each pair merged.
#### -tabbedout give you extensive information on the quality of the merge
#### This step takes approximately 10 minutes
```

cd /mnt/research/EvansLab/Lukas/MMPRNT018/
#mkdir Fungal_libs
mkdir Fungal_libs/20200306_M03127_Amplicon_PE300
cd Fungal_libs/20200306_M03127_Amplicon_PE300
#I am going to submit a job

nano merg_seq_lib10_MMPRNT018_soil_fungi.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=3:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name merge_lib10_MMPRNT018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/MMPRNT_2018/20200306_M03127_Amplicon_PE300/QIIME2_demux/output/demultiplexed-seqs/

usearch11.0.667 -fastq_mergepairs *_R1*.fastq -relabel @ -fastqout /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200306_M03127_Amplicon_PE300/lib10_20200306_Amplicon_PE300_ITS1.fastq  -tabbedout /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200306_M03127_Amplicon_PE300/lib10_20200306_Amplicon_PE300_ITS1_pair_report.txt 

###End
sbatch merg_seq_lib10_MMPRNT018_soil_fungi.sbatch
#Submitted batch job 58339234

```

Merge stats

---
Totals:
  31137769  Pairs (31.1M)
  26702553  Merged (26.7M, 85.76%)
  18781464  Alignments with zero diffs (60.32%)
   4155551  Too many diffs (> 5) (13.35%)
    279665  No alignment found (0.90%)
         0  Alignment too short (< 16) (0.00%)
  17346841  Staggered pairs (55.71%) merged & trimmed
    266.77  Mean alignment length
    301.58  Mean merged length
      0.48  Mean fwd expected errors
      1.13  Mean rev expected errors
      0.06  Mean merged expected errors



---

###Check sample names to make sure all of the samples are present


```
usearch11.0.667 -fastx_get_sample_names /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200306_M03127_Amplicon_PE300/lib10_20200306_Amplicon_PE300_ITS1.fastq -output /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200306_M03127_Amplicon_PE300/lib10_20200306_Amplicon_PE300_ITS1_samples.txt
```
---
# Should be 192 samples but one sample lost during demux (MMPRNT1434)
01:38 38Mb    100.0% 191 samples found
---


## Let  demultiplex the reads from 20200309_M02808_Amplicon_PE300

### GLBRC_Lib12_ITS_Evans_soil_2018

```
cd /mnt/research/EvansLab/MMPRNT_2018/20200309_M02808_Amplicon_PE300

```
I created a bar code file in a text editor

### Let's use the QIIME2 demux
https://docs.qiime2.org/2020.2/plugins/available/demux/emp-paired/?highlight=demux

###### Import the multiplexed sequences
We need to rename and move the zipped sequence files
```
mkdir QIIME2_demux
cd QIIME2_demux/
mkdir compressed_seqs
cd ..

cp GLBRC_Lib12_ITS_Evans_soil_2018_S1_L001_R1_001.fastq.gz  QIIME2_demux/compressed_seqs/forward.fastq.gz
cp GLBRC_Lib12_ITS_Evans_soil_2018_S1_L001_R2_001.fastq.gz  QIIME2_demux/compressed_seqs/reverse.fastq.gz
cp GLBRC_Lib12_ITS_Evans_soil_2018_S1_L001_I1_001.fastq.gz  QIIME2_demux/compressed_seqs/barcodes.fastq.gz

#Now we need to convert the files to Q2 format and demux 

#qiime tools import --show-importable-types


cd QIIME2_demux

nano Q2_demux_lib12_MMPRNT018_soil_fungi.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=20:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name Q2_demux_lib12_MMPRNT018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/MMPRNT_2018/20200309_M02808_Amplicon_PE300/QIIME2_demux
export PATH="/mnt/research/EvansLab/Software/anaconda3/bin:$PATH"
source activate qiime2-2019.4

#convert the files to the qza format
qiime tools import --type EMPPairedEndSequences --input-path compressed_seqs --output-path multiplexed-seqs_paired_end.qza


#demux the sequences
qiime demux emp-paired --m-barcodes-file MMPRNT_QIIME_mapping_lib12.txt --m-barcodes-column BarcodeSequence --p-rev-comp-mapping-barcodes --p-no-golay-error-correction --i-seqs multiplexed-seqs_paired_end.qza --o-per-sample-sequences demultiplexed-seqs.qza --o-error-correction-details demux-details.qza

#Let's also export the demux files

qiime tools export  --input-path demultiplexed-seqs.qza  --output-path output/demultiplexed-seqs

#Let's try to visualize the demux sequences too

qiime demux summarize  --i-data demultiplexed-seqs.qza  --o-visualization demux-seq_summary_lib12.qzv

###End

sbatch Q2_demux_lib12_MMPRNT018_soil_fungi.sbatch
#Submitted batch job 58241079
#03:45:18
```

Demux stats

---
Minimum:	328
Median:	113550.5
Mean:	115265.58421052631
Maximum:	258760
Total:	21900461

---
Let's unzip the resulting demultiplex files
```
cd  /mnt/research/EvansLab/MMPRNT_2018/20200309_M02808_Amplicon_PE300/QIIME2_demux/output/demultiplexed-seqs/
gunzip *.fastq.gz
```

## 1) Quality checking
### 1a) First look at the quality of  raw unmerged seqs for run1
#https://www.drive5.com/usearch/manual/pipe_readprep_understand.html

```
mkdir fastq_info
```

#### make a forloop to run fastx_info on every file
```
nano fasta_info_fq.sh
#!/bin/bash --login
for fq in *.fastq
do
usearch10.0.240 -fastx_info $fq -output fastq_info/$fq
done
```

#### make file executable and run the for loop create in your `fasta_info_fq.sh` file
```
chmod +x fasta_info_fq.sh
./fasta_info_fq.sh
```

#### move to the fastq_info directory
```
cd fastq_info/
```

#### Now run the below code to summarize the fastq info for all of the forward and reverse reads

```
grep "^File" * > 20200309_M02808_Amplicon_PE300_fastq_lengths.txt
grep "^EE" * > 20200309_M02808_Amplicon_PE300_fastq_EE.txt
```

Look for any forward and reverse reads that look especially bad in terms of quality (high E is bad quality). This info will also be really helpful for troubleshooting later on (e.g. why some samples have extremely low read numbers)


##  2) Merge the forward and reverse sequences and trim adapters (for each run individually)
### 2a) Merge Pairs
#### Make sure you are in the folder with the extracted forward and reverse reads from Run 1
https://www.drive5.com/usearch/manual/merge_options.html
#### -alnout gives you a human readable text file of the alignment and misalignments for each pair merged.
#### -tabbedout give you extensive information on the quality of the merge
#### This step takes approximately 10 minutes
```

cd /mnt/research/EvansLab/Lukas/MMPRNT018/
#mkdir Fungal_libs
mkdir Fungal_libs/20200309_M02808_Amplicon_PE300
cd Fungal_libs/20200309_M02808_Amplicon_PE300
#I am going to submit a job

nano merg_seq_lib12_MMPRNT018_soil_fungi.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=3:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name merge_lib12_MMPRNT018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/MMPRNT_2018/20200309_M02808_Amplicon_PE300/QIIME2_demux/output/demultiplexed-seqs/

usearch11.0.667 -fastq_mergepairs *_R1*.fastq -relabel @ -fastqout /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200309_M02808_Amplicon_PE300/lib12_20200309_Amplicon_PE300_ITS1.fastq  -tabbedout /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200309_M02808_Amplicon_PE300/lib12_20200309_Amplicon_PE300_ITS1_pair_report.txt 

###End
sbatch merg_seq_lib12_MMPRNT018_soil_fungi.sbatch
#Submitted batch job 58339876

```

Merge stats

---
Totals:
  21900461  Pairs (21.9M)
  18164594  Merged (18.2M, 82.94%)
  11979272  Alignments with zero diffs (54.70%)
   3484535  Too many diffs (> 5) (15.91%)
    251332  No alignment found (1.15%)
         0  Alignment too short (< 16) (0.00%)
  12234281  Staggered pairs (55.86%) merged & trimmed
    266.33  Mean alignment length
    299.84  Mean merged length
      0.69  Mean fwd expected errors
      1.69  Mean rev expected errors
      0.05  Mean merged expected errors



---

###Check sample names to make sure all of the samples are present


```
usearch11.0.667 -fastx_get_sample_names /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200309_M02808_Amplicon_PE300/lib12_20200309_Amplicon_PE300_ITS1.fastq -output /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200309_M02808_Amplicon_PE300/lib12_20200309_Amplicon_PE300_ITS1_samples.txt
```
```
#Should be 192 samples but there were two missing after demux (MMPRNT1842 and MMPRNT1862)

01:08 38Mb    100.0% 190 samples found
```

## Let  demultiplex the reads from 20200309_M03127_Amplicon_PE300

### GLBRC_Lib13_ITS_Evans_soil_2018

```
cd /mnt/research/EvansLab/MMPRNT_2018/20200309_M03127_Amplicon_PE300

#Look at the quality of the index file 
usearch11.0.667 -fastq_eestats2 GLBRC_Lib13_ITS_Evans_soil_2018_S1_L001_I1_001.fastq -output  eestats2_GLBRC_Lib13_ITS_Evans_soil_2018_S1_L001_I1_001.txt

```
I created a bar code file in a text editor

### Let's use the QIIME2 demux
https://docs.qiime2.org/2020.2/plugins/available/demux/emp-paired/?highlight=demux

###### Import the multiplexed sequences
We need to rename and move the zipped sequence files
```
mkdir QIIME2_demux
cd QIIME2_demux/
mkdir compressed_seqs
cd ..

cp GLBRC_Lib13_ITS_Evans_soil_2018_S1_L001_R1_001.fastq.gz  QIIME2_demux/compressed_seqs/forward.fastq.gz
cp GLBRC_Lib13_ITS_Evans_soil_2018_S1_L001_R2_001.fastq.gz  QIIME2_demux/compressed_seqs/reverse.fastq.gz
cp GLBRC_Lib13_ITS_Evans_soil_2018_S1_L001_I1_001.fastq.gz  QIIME2_demux/compressed_seqs/barcodes.fastq.gz

#Now we need to convert the files to Q2 format and demux 

#qiime tools import --show-importable-types


cd QIIME2_demux

nano Q2_demux_lib13_MMPRNT018_soil_fungi.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=20:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name Q2_demux_lib13_MMPRNT018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/MMPRNT_2018/20200309_M03127_Amplicon_PE300/QIIME2_demux
export PATH="/mnt/research/EvansLab/Software/anaconda3/bin:$PATH"
source activate qiime2-2019.4

#convert the files to the qza format
qiime tools import --type EMPPairedEndSequences --input-path compressed_seqs --output-path multiplexed-seqs_paired_end.qza


#demux the sequences
qiime demux emp-paired --m-barcodes-file MMPRNT_QIIME_mapping_lib13.txt --m-barcodes-column BarcodeSequence --p-rev-comp-mapping-barcodes --p-no-golay-error-correction --i-seqs multiplexed-seqs_paired_end.qza --o-per-sample-sequences demultiplexed-seqs.qza --o-error-correction-details demux-details.qza

#Let's also export the demux files

qiime tools export  --input-path demultiplexed-seqs.qza  --output-path output/demultiplexed-seqs

#Let's try to visualize the demux sequences too

qiime demux summarize  --i-data demultiplexed-seqs.qza  --o-visualization demux-seq_summary_lib13.qzv

###End

sbatch Q2_demux_lib13_MMPRNT018_soil_fungi.sbatch
#Submitted batch job 58711226
#Run time 05:45:42 first time (rerun qith correct barcode)
```

Demux stats

---

Minimum:	2
Median:	179349.0
Mean:	169704.16753926701
Maximum:	449166
Total:	32413496
---

Unzip the demultiplexed sequences

```

cd /mnt/research/EvansLab/MMPRNT_2018/20200309_M03127_Amplicon_PE300/QIIME2_demux/output/demultiplexed-seqs


```

## 1) Quality checking
### 1a) First look at the quality of  raw unmerged seqs for run1
#https://www.drive5.com/usearch/manual/pipe_readprep_understand.html

```
mkdir fastq_info
```

#### make a forloop to run fastx_info on every file
```
nano fasta_info_fq.sh
#!/bin/bash --login
for fq in *.fastq
do
usearch10.0.240 -fastx_info $fq -output fastq_info/$fq
done
```

#### make file executable and run the for loop create in your `fasta_info_fq.sh` file
```
chmod +x fasta_info_fq.sh
./fasta_info_fq.sh
```

#### move to the fastq_info directory
```
cd fastq_info/
```

#### Now run the below code to summarize the fastq info for all of the forward and reverse reads

```
grep "^File" * > 20200309_M03127_Amplicon_PE300_fastq_lengths.txt
grep "^EE" * > 20200309_M03127_Amplicon_PE300_fastq_EE.txt
```

Look for any forward and reverse reads that look especially bad in terms of quality (high E is bad quality). This info will also be really helpful for troubleshooting later on (e.g. why some samples have extremely low read numbers)

##  2) Merge the forward and reverse sequences and trim adapters (for each run individually)
### 2a) Merge Pairs
#### Make sure you are in the folder with the extracted forward and reverse reads from Run 1
https://www.drive5.com/usearch/manual/merge_options.html
#### -alnout gives you a human readable text file of the alignment and misalignments for each pair merged.
#### -tabbedout give you extensive information on the quality of the merge
#### This step takes approximately 10 minutes
```

cd /mnt/research/EvansLab/Lukas/MMPRNT018/
#mkdir Fungal_libs
mkdir Fungal_libs/20200309_M03127_Amplicon_PE300
cd Fungal_libs/20200309_M03127_Amplicon_PE300
#I am going to submit a job

nano merg_seq_lib13_MMPRNT018_soil_fungi.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=1:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name merge_cut_lib13_MMPRNT018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/MMPRNT_2018/20200309_M03127_Amplicon_PE300/QIIME2_demux/output/demultiplexed-seqs/

usearch11.0.667 -fastq_mergepairs *_R1*.fastq -relabel @ -fastqout /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200309_M03127_Amplicon_PE300/lib13_demux_20200309_Amplicon_PE300_ITS1-2.fastq  -tabbedout /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200309_M03127_Amplicon_PE300/lib13_demux_20200309_Amplicon_PE300_ITS1-2_pair_report.txt -fastqout_notmerged_fwd /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200309_M03127_Amplicon_PE300/NO_Match_lib13_demux_20200309_Amplicon_PE300_ITS1-2_for.fastq -fastqout_notmerged_rev /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200309_M03127_Amplicon_PE300/NO_Match_lib13_demux_20200309_Amplicon_PE300_ITS1-2_rev.fastq

###End

sbatch merg_seq_lib13_MMPRNT018_soil_fungi.sbatch

#Submitted batch job 58776582

```

Merge stats

---
#Run time 00:09:48
Totals:
  32413496  Pairs (32.4M)
  27393287  Merged (27.4M, 84.51%)
  18303690  Alignments with zero diffs (56.47%)
   4665259  Too many diffs (> 5) (14.39%)
    354950  No alignment found (1.10%)
         0  Alignment too short (< 16) (0.00%)
  15473712  Staggered pairs (47.74%) merged & trimmed
    261.80  Mean alignment length
    305.66  Mean merged length
      0.66  Mean fwd expected errors
      1.37  Mean rev expected errors
      0.07  Mean merged expected errors
---

###Check sample names to make sure all of the samples are present

```
usearch11.0.667 -fastx_get_sample_names /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200309_M03127_Amplicon_PE300/lib13_demux_20200309_Amplicon_PE300_ITS1-2.fastq -output /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200309_M03127_Amplicon_PE300/lib13_demux_20200309_Amplicon_PE300_ITS1-2_samples.txt
```
```
#Should be 191 samples
01:46 38Mb    100.0% 191 samples found

```

#Bonito Lab Root Field Libraries

## Let  demultiplex the reads from Lib_8

```
#20326995 seqs
cd /mnt/research/EvansLab/Lukas
mkdir GLBRC018_Root_Soil
```


+ /mnt/research/EvansLab/GLBRC_SG_Roots/Field2018_Fun/Lib_8_backup/

Pedro created a bar code file in a text editor

### Let's use the QIIME2 demux
https://docs.qiime2.org/2020.2/plugins/available/demux/emp-paired/?highlight=demux

###### Import the multiplexed sequences
We need to rename and move the zipped sequence files
```
cd /mnt/research/EvansLab/GLBRC_SG_Roots/Field2018_Fun/Lib_8_backup/
mkdir QIIME2_demux
cd QIIME2_demux/
mkdir compressed_seqs
cd ..

cp GLBRC_Lib8_ITS_soil_2018_S1_L001_R1_001.fastq.gz  QIIME2_demux/compressed_seqs/forward.fastq.gz
cp GLBRC_Lib8_ITS_soil_2018_S1_L001_R2_001.fastq.gz  QIIME2_demux/compressed_seqs/reverse.fastq.gz
cp GLBRC_Lib8_ITS_soil_2018_S1_L001_I1_001.fastq.gz  QIIME2_demux/compressed_seqs/barcodes.fastq.gz

#Now we need to convert the files to Q2 format and demux 

#qiime tools import --show-importable-types


cd QIIME2_demux

nano Q2_demux_lib8_GLBRC018_root_fungi.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=10:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name Q2_demux_lib8_GLBRC018_root      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd ${SLURM_SUBMIT_DIR}
export PATH="/mnt/research/EvansLab/Software/anaconda3/bin:$PATH"
source activate qiime2-2019.4

#convert the files to the qza format 
#qiime tools import --type EMPPairedEndSequences --input-path compressed_seqs --output-path multiplexed-seqs_paired_end.qza this was run previously
#deleted compressed file after this was run

#demux the sequences
qiime demux emp-paired --m-barcodes-file Lib_8_mapping.txt --m-barcodes-column BarcodeSequence --p-rev-comp-mapping-barcodes --p-no-golay-error-correction --i-seqs multiplexed-seqs_paired_end.qza --o-per-sample-sequences lib8_demultiplexed-seqs.qza --o-error-correction-details lib8_demux-details.qza

#Let's also export the demux files

qiime tools export  --input-path lib8_demultiplexed-seqs.qza  --output-path output/demultiplexed-seqs

#Let's try to visualize the demux sequences too

qiime demux summarize  --i-data lib8_demultiplexed-seqs.qza  --o-visualization demux-seq_summary_lib8.qzv

###End

sbatch Q2_demux_lib8_GLBRC018_root_fungi.sbatch
#Submitted batch job 20549039
#Run time 02:50:47
```

Demux stats

```
Demultiplexed sequence counts summary
Minimum:	1
Median:	42418.5
Mean:	72744.43233082707
Maximum:	350491
Total:	19350019

```

Unzip the demultiplexed sequences

```

cd /mnt/research/EvansLab/GLBRC_SG_Roots/Field2018_Fun/Lib_8_backup/QIIME2_demux/output/demultiplexed-seqs
gunzip *.fastq

```

## 1) Quality checking
### 1a) First look at the quality of  raw unmerged seqs for run1
#https://www.drive5.com/usearch/manual/pipe_readprep_understand.html

```
mkdir fastq_info
```

#### make a forloop to run fastx_info on every file
```
nano fasta_info_fq.sh
#!/bin/bash --login
for fq in *.fastq
do
usearch10.0.240 -fastx_info $fq -output fastq_info/$fq
done
```

#### make file executable and run the for loop create in your `fasta_info_fq.sh` file
```
chmod +x fasta_info_fq.sh
./fasta_info_fq.sh
```

#### move to the fastq_info directory
```
cd fastq_info/
```

#### Now run the below code to summarize the fastq info for all of the forward and reverse reads

```
grep "^File" * > Lib8_Amplicon_PE300_fastq_lengths.txt
grep "^EE" * > Lib8_Amplicon_PE300_fastq_EE.txt
```

Look for any forward and reverse reads that look especially bad in terms of quality (high E is bad quality). This info will also be really helpful for troubleshooting later on (e.g. why some samples have extremely low read numbers)

##  2) Merge the forward and reverse sequences and trim adapters (for each run individually)
### 2a) Merge Pairs
#### Make sure you are in the folder with the extracted forward and reverse reads from Run 1
https://www.drive5.com/usearch/manual/merge_options.html
#### -alnout gives you a human readable text file of the alignment and misalignments for each pair merged.
#### -tabbedout give you extensive information on the quality of the merge
#### This step takes approximately 10 minutes
```

cd /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil
mkdir Fungal_libs
mkdir Fungal_libs/LIb8_fun
cd Fungal_libs/LIb8_fun

#I am going to submit a job

nano merg_seq_lib8_GLBRC018_soil_fungi.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=1:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name merge_lib8_GLBRC018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/GLBRC_SG_Roots/Field2018_Fun/Lib_8_backup/QIIME2_demux/output/demultiplexed-seqs

usearch11.0.667 -fastq_mergepairs *_R1*.fastq -relabel @ -fastqout /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/LIb8_fun/merged_lib8_ITS1.fastq  -tabbedout /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/LIb8_fun/merged_lib8_ITS1_pair_report.txt -fastqout_notmerged_fwd /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/LIb8_fun/NO_Match_lib8_Amplicon_PE300_ITS1_for.fastq -fastqout_notmerged_rev /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/LIb8_fun/NO_Match_lib8_Amplicon_PE300_ITS1_rev.fastq

###End

sbatch merg_seq_lib8_GLBRC018_soil_fungi.sbatch

#Submitted batch job 20558017

```

Merge stats

```
Totals:
  19350019  Pairs (19.4M)
  16093657  Merged (16.1M, 83.17%)
   8283746  Alignments with zero diffs (42.81%)
   2895146  Too many diffs (> 5) (14.96%)
    361216  No alignment found (1.87%)
         0  Alignment too short (< 16) (0.00%)
   8555442  Staggered pairs (44.21%) merged & trimmed
    249.74  Mean alignment length
    304.23  Mean merged length
      0.85  Mean fwd expected errors
      1.88  Mean rev expected errors
      0.08  Mean merged expected errors
```

###Check sample names to make sure all of the samples are present

```
usearch11.0.667 -fastx_get_sample_names /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/LIb8_fun/merged_lib8_ITS1.fastq -output /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/LIb8_fun/merged_lib8_ITS1_samples.txt
```
```
#Should be  samples 268
#100.0% 265 samples found
#Amp2037 Amp2025 Amp1952 are missing from the seqs because they did not amplify or sequence

```
#I want to check the quality

```
usearch11.0.667 -fastq_eestats2 /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/LIb8_fun/merged_lib8_ITS1.fastq -output /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/LIb8_fun/merged_lib8_ITS1_eestats2.txt
```
```
03:25 661Mb   100.0% Reading reads

16093657 reads, max len 584, avg 304.2

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50   15746815( 97.8%)   15764316( 98.0%)   15765085( 98.0%)
   100   15465345( 96.1%)   15520222( 96.4%)   15524027( 96.5%)
   150   15408007( 95.7%)   15512197( 96.4%)   15521657( 96.4%)
   200   15363238( 95.5%)   15503548( 96.3%)   15517783( 96.4%)
   250   13494449( 83.8%)   13648603( 84.8%)   13667055( 84.9%)
   300    8162085( 50.7%)    8290121( 51.5%)    8312062( 51.6%)
   350    2856972( 17.8%)    3005080( 18.7%)    3064286( 19.0%)
   400    1712082( 10.6%)    1850713( 11.5%)    1924619( 12.0%)
   450     274110(  1.7%)     318656(  2.0%)     352592(  2.2%)
   500      64600(  0.4%)      88075(  0.5%)     107983(  0.7%)
   550      28208(  0.2%)      43657(  0.3%)      57028(  0.4%)

```

## Let  demultiplex the reads from Lib_5

```
#33877753 seq
#cd /mnt/research/EvansLab/Lukas
#mkdir GLBRC018_Root_Soil
```


+ /mnt/research/EvansLab/GLBRC_SG_Roots/Field2018_Fun/Lib_8_backup/

Pedro created a bar code file in a text editor

### Let's use the QIIME2 demux
https://docs.qiime2.org/2020.2/plugins/available/demux/emp-paired/?highlight=demux

###### Import the multiplexed sequences
We need to rename and move the zipped sequence files
```
cd /mnt/research/EvansLab/GLBRC_SG_Roots/Field2018_Fun/Lib_5_backup/
mkdir QIIME2_demux
cd QIIME2_demux/
mkdir compressed_seqs
cd ..

cp GLBRC_Lib5_ITS_root_2018_S1_L001_R1_001.fastq.gz  QIIME2_demux/compressed_seqs/forward.fastq.gz
cp GLBRC_Lib5_ITS_root_2018_S1_L001_R2_001.fastq.gz  QIIME2_demux/compressed_seqs/reverse.fastq.gz
cp GLBRC_Lib5_ITS_root_2018_S1_L001_I1_001.fastq.gz  QIIME2_demux/compressed_seqs/barcodes.fastq.gz

#Now we need to convert the files to Q2 format and demux 

#qiime tools import --show-importable-types


cd QIIME2_demux

nano Q2_demux_lib5_GLBRC018_root_fungi.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=10:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name Q2_demux_lib5_GLBRC018_root      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd ${SLURM_SUBMIT_DIR}
export PATH="/mnt/research/EvansLab/Software/anaconda3/bin:$PATH"
source activate qiime2-2019.4

#convert the files to the qza format
qiime tools import --type EMPPairedEndSequences --input-path compressed_seqs --output-path lib5_multiplexed-seqs_paired_end.qza
#deleted compressed file after this was run

#demux the sequences
qiime demux emp-paired --m-barcodes-file Lib_5_mapping.txt --m-barcodes-column BarcodeSequence --p-rev-comp-mapping-barcodes --p-no-golay-error-correction --i-seqs lib5_multiplexed-seqs_paired_end.qza --o-per-sample-sequences lib5_demultiplexed-seqs.qza --o-error-correction-details lib5_demux-details.qza

#Let's also export the demux files

qiime tools export  --input-path lib5_demultiplexed-seqs.qza  --output-path output/demultiplexed-seqs

#Let's try to visualize the demux sequences too

qiime demux summarize  --i-data lib5_demultiplexed-seqs.qza  --o-visualization demux-seq_summary_lib5.qzv

###End

sbatch Q2_demux_lib5_GLBRC018_root_fungi.sbatch
#Submitted batch job 20549052
#Run time 
```

Demux stats

```
Demultiplexed sequence counts summary
Minimum:	2001
Median:	100965.0
Mean:	96414.60597014925
Maximum:	179774
Total:	32298893

```

Unzip the demultiplexed sequences

```

cd /mnt/research/EvansLab/GLBRC_SG_Roots/Field2018_Fun/Lib_5_backup/QIIME2_demux/output/demultiplexed-seqs
gunzip *.fastq.gz

```

## 1) Quality checking
### 1a) First look at the quality of  raw unmerged seqs for run1
#https://www.drive5.com/usearch/manual/pipe_readprep_understand.html

```
mkdir fastq_info
```

#### make a forloop to run fastx_info on every file
```
nano fasta_info_fq.sh
#!/bin/bash --login
for fq in *.fastq
do
usearch10.0.240 -fastx_info $fq -output fastq_info/$fq
done
```

#### make file executable and run the for loop create in your `fasta_info_fq.sh` file
```
chmod +x fasta_info_fq.sh
./fasta_info_fq.sh
```

#### move to the fastq_info directory
```
cd fastq_info/
```

#### Now run the below code to summarize the fastq info for all of the forward and reverse reads

```
grep "^File" * > lib5_Amplicon_PE300_fastq_lengths.txt
grep "^EE" * > lib5_Amplicon_PE300_fastq_EE.txt
```

Look for any forward and reverse reads that look especially bad in terms of quality (high E is bad quality). This info will also be really helpful for troubleshooting later on (e.g. why some samples have extremely low read numbers)

##  2) Merge the forward and reverse sequences and trim adapters (for each run individually)
### 2a) Merge Pairs
#### Make sure you are in the folder with the extracted forward and reverse reads from Run 1
https://www.drive5.com/usearch/manual/merge_options.html
#### -alnout gives you a human readable text file of the alignment and misalignments for each pair merged.
#### -tabbedout give you extensive information on the quality of the merge
#### This step takes approximately 10 minutes
```

cd /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil
#mkdir Fungal_libs
mkdir Fungal_libs/lib5_fun
cd Fungal_libs/lib5_fun

#I am going to submit a job

nano merg_seq_lib5_GLBRC018_soil_fungi.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=1:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name merge_lib5_GLBRC018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/GLBRC_SG_Roots/Field2018_Fun/Lib_5_backup/QIIME2_demux/output/demultiplexed-seqs

usearch11.0.667 -fastq_mergepairs *_R1*.fastq -relabel @ -fastqout /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/lib5_fun/merged_lib5_ITS1.fastq  -tabbedout /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/lib5_fun/merged_lib5_ITS1_pair_report.txt -fastqout_notmerged_fwd /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/lib5_fun/NO_Match_lib5_Amplicon_PE300_ITS1_for.fastq -fastqout_notmerged_rev /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/lib5_fun/NO_Match_lib5_Amplicon_PE300_ITS1_rev.fastq

###End

sbatch merg_seq_lib5_GLBRC018_soil_fungi.sbatch

#Submitted batch job 20578491

```

Merge stats

```
1:53 1.1Gb   100.0% 85.5% merged

Totals:
  32298893  Pairs (32.3M)
  27630451  Merged (27.6M, 85.55%)
  19015773  Alignments with zero diffs (58.87%)
   4312799  Too many diffs (> 5) (13.35%)
    355643  No alignment found (1.10%)
         0  Alignment too short (< 16) (0.00%)
  19735989  Staggered pairs (61.10%) merged & trimmed
    237.15  Mean alignment length
    263.49  Mean merged length
      2.09  Mean fwd expected errors
      3.37  Mean rev expected errors
      0.04  Mean merged expected errors
```

###Check sample names to make sure all of the samples are present

```
usearch11.0.667 -fastx_get_sample_names /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/lib5_fun/merged_lib5_ITS1.fastq -output /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/lib5_fun/merged_lib5_ITS1_samples.txt
```
```
#Should be 336 samples
#01:31 38Mb    100.0% 335 samples found
#Amp1853 is missing

```

#I want to check the quality

```
usearch11.0.667 -fastq_eestats2 /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/lib5_fun/merged_lib5_ITS1.fastq -output /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/lib5_fun/merged_lib5_ITS1_eestats2.txt
```
```
02:46 661Mb   100.0% Reading reads

27630451 reads, max len 584, avg 263.5

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50   27134006( 98.2%)   27151178( 98.3%)   27151827( 98.3%)
   100   23285446( 84.3%)   23316192( 84.4%)   23317287( 84.4%)
   150   23242647( 84.1%)   23307095( 84.4%)   23309668( 84.4%)
   200   23205931( 84.0%)   23297023( 84.3%)   23301760( 84.3%)
   250   19704043( 71.3%)   19800277( 71.7%)   19806119( 71.7%)
   300    9408334( 34.1%)    9468566( 34.3%)    9474834( 34.3%)
   350    2331476(  8.4%)    2401364(  8.7%)    2427559(  8.8%)
   400     307514(  1.1%)     333243(  1.2%)     349569(  1.3%)
   450      83675(  0.3%)      96851(  0.4%)     107862(  0.4%)
   500      26700(  0.1%)      33684(  0.1%)      39993(  0.1%)
   550       1057(  0.0%)       1596(  0.0%)       2169(  0.0%)


```

## 3) Combine the seven merged sequence files 
### Make sure you do not use replicate sample names between the runs
### This step is only necessary if more than one Illumina run are being analyzed together.

```

cat /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200302_Amplicon_PE300/lib11_20200302_Amplicon_PE300_ITS1.fastq /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200306_M02808_Amplicon_PE300/lib09_20200302_Amplicon_PE300_ITS1.fastq /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200306_M03127_Amplicon_PE300/lib10_20200306_Amplicon_PE300_ITS1.fastq /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200309_M02808_Amplicon_PE300/lib12_20200309_Amplicon_PE300_ITS1.fastq /mnt/research/EvansLab/Lukas/MMPRNT018/Fungal_libs/20200309_M03127_Amplicon_PE300/lib13_demux_20200309_Amplicon_PE300_ITS1-2.fastq /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/lib5_fun/merged_lib5_ITS1.fastq /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/LIb8_fun/merged_lib8_ITS1.fastq > /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/combined_merged_ITS1_2_GLBRC018.fastq


```

Before we continue, you may want to check if the sample names are formatted correctly. USEARCH does some funny cutting during the merging step. Any hyphens or underscores can be problematic and you need to remove these (use sed command and merged_cut files)

Additionally, this is a good opportunity to double check that all of your samples merged and have unique IDs using [fastx_get_sample_names](https://www.drive5.com/usearch/manual/cmd_fastx_get_sample_names.html)

```
usearch11.0.667 -fastx_get_sample_names /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/combined_merged_ITS1_2_GLBRC018.fastq -output /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/combined_merged_ITS1_2_GLBRC018_samples.txt
```
```
#Should be 192 + 191 + 190 + 190 + 191 + 265 + 335 = 1554 samples
#09:49 38Mb    100.0% 1554 samples found

```
Looks good


## You can run the USEARCH version of phix removal
#https://www.drive5.com/usearch/manual/cmd_filter_phix.html
I ran a test about only had 43 hits within 87% of seqs so I am goign to skip this step


## Let's check quality and remove the primers and left over adapt and barcodes and such

```
#I am going to submit a job
cd /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/

nano qual_prefilter_GLBRC018_soil_root_fungi.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=1:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name qual_check_perfil_GLBRC018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/
usearch11.0.667 -fastq_eestats2 combined_merged_ITS1_2_GLBRC018.fastq -output combined_merged_ITS1_2_GLBRC018_eestats2.txt

#end

sbatch qual_prefilter_GLBRC018_soil_root_fungi.sbatch
#Submitted batch job 20586770
```
```
07:13 661Mb   100.0% Reading reads

158659764 reads, max len 584, avg 295.7

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50  156621481( 98.7%)  156754683( 98.8%)  156757709( 98.8%)
   100  151172438( 95.3%)  151507375( 95.5%)  151520630( 95.5%)
   150  150846326( 95.1%)  151459561( 95.5%)  151493500( 95.5%)
   200  150532083( 94.9%)  151394962( 95.4%)  151455103( 95.5%)
   250  137859263( 86.9%)  138866050( 87.5%)  138953307( 87.6%)
   300   68304166( 43.1%)   69011441( 43.5%)   69113788( 43.6%)
   350   19391825( 12.2%)   20246279( 12.8%)   20592944( 13.0%)
   400    6638491(  4.2%)    7241630(  4.6%)    7616774(  4.8%)
   450    1949315(  1.2%)    2256598(  1.4%)    2502337(  1.6%)
   500     535301(  0.3%)     676654(  0.4%)     802618(  0.5%)
   550      65723(  0.0%)      97331(  0.1%)     127910(  0.1%)
```



## 4) Filtering and Truncate the merged seqs  to MaxEE and set length using [fastq_filter](https://www.drive5.com/usearch/manual/cmd_fastq_filter.html)

#### 100 bp is within the expected amplification length and >94.9% of sequences meet this length cut off
My processes have been killed frequently... I think big HPCC is on to me!
```

nano qual_primer_filter_GLBRC018_soil_root_fungi.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=1:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name qual_prim_filt_perfil_GLBRC018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/


usearch11.0.667 -fastq_filter combined_merged_ITS1_2_GLBRC018.fastq -fastq_maxee 0.1 -fastq_minlen 100 -fastaout combined_merged_ITS1_2_GLBRC018_filtered.fa
####END
sbatch qual_primer_filter_GLBRC018_soil_root_fungi.sbatch
#Submitted batch job 20587265
```
```
25:53 282Mb   100.0% Filtering, 87.7% passed
 158659764  Reads (159M)                    
  12381854  Discarded reads with expected errs > 0.10
 139139092  Filtered reads (139M, 87.7%)

```

## [USEARCH phix removal](https://www.drive5.com/usearch/manual/cmd_filter_phix.html)

```

nano phix_filter_GLBRC018_soil_root_fungi.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=1:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name phix_filt_perfil_GLBRC018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu
cd /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/

usearch11.0.667 -filter_phix combined_merged_ITS1_2_GLBRC018_filtered.fa -output combined_merged_ITS1_2_GLBRC018_filtered_no_phix.fa -alnout combined_merged_ITS1_2_GLBRC018_filtered_phix_hits.txt
###END
.
sbatch phix_filter_GLBRC018_soil_root_fungi.sbatch
#Submitted batch job 20588874
#25:45 526Mb   100.0% Filtering for phix, 35 hits (0.0%)
```

### Now let's remove the primers from the entire fastq dataset
```
#make a fasta file with primers
nano EMP_AM__ITS1-2_fungal_primers.fa
>ITS1f 
CTTGGTCATTTAGAGGAAGTAA
>ITS4R
TCCTCCGCTTATTGATATGC
>ITS2r
GCTGCGTTCTTCATCGATGC
>fITS7o
GTGAATCATCRAATYTTT
>fITS7
GTGARTCATCGAATCTTTG

nano primer_filter_GLBRC018_soil_root_fungi.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=3:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name primer_filt_perfil_GLBRC018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/

usearch11.0.667 -fastx_trim_primer combined_merged_ITS1_2_GLBRC018_filtered.fa -db EMP_AM__ITS1-2_fungal_primers.fa -strand both  -maxdiffs 4 -width 40 -fastaout combined_merged_ITS1_2_GLBRC018_primer_filtered.fa -tabbedout combined_merged_ITS1_2_GLBRC018_primer_filtered_primer_re.txt


###

sbatch primer_filter_GLBRC018_soil_root_fungi.sbatch

#Submitted batch job 20606257
#33:51 44Mb    100.0% Processing
grep -c "^>" combined_merged_ITS1_2_GLBRC018_primer_filtered.fa 
#combined phix and primer filtering
3139139092-139116802=22290

```


## 5) Filter so we only have unique sequences with [fastx_uniques](https://www.drive5.com/usearch/manual/cmd_fastx_uniques.html)

## 6) Cluster into OTUS and filter out singletons
There are two options here. **(A)** uses the traditional approach and clusters sequences into 0.97 identity cutoff OTUs. **(B)** uses unoise3 to identify ZOTUs.

### 6A) Cluster into 0.97 OTUs using UPARSE and [cluster_otus](https://www.drive5.com/usearch/manual/cmd_cluster_otus.html)
This step will also denovo chimera check and filter out singletons.You can remove single sequences prior to clustering but singletons are also removed at the OTU clustering step (cluster_otus filters out OTUs <2 and unoise3 filters ZOTUs <8)

This step takes approximately 30 minutes with ~330,000 unique sequences


```

nano OTU_clustering_fungi_GLBRC018.sbatch

#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=50:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=3-4                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=3                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   clust_OTU_ZOTU_Fung_GLBRC018  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/

usearch11.0.667 -fastx_uniques combined_merged_ITS1_2_GLBRC018_primer_filtered.fa -fastaout uniques_combined_merged_ITS1_2_GLBRC018_primer_filtered.fa -sizeout

usearch11.0.667 -cluster_otus uniques_combined_merged_ITS1_2_GLBRC018_primer_filtered.fa -otus rep_set_combined_merged_ITS1_2_GLBRC_otus.fa -relabel OTU


###end of .sbatch

sbatch OTU_clustering_fungi_GLBRC018.sbatch
#Submitted batch job 20606856
```

```
#Derep
13:41 69.1Gb 139116802 seqs, 131671388 uniques, 126553588 singletons (96.1%)
13:41 69.1Gb Min size 1, median 1, max 238, avg 1.06

01:46:19 57.3Gb  100.0% Writing uniques_combined_merged_ITS1_2_GLBRC018_primer_filtered.fa


usearch11.0.667 -fastx_uniques combined_merged_ITS1_2_GLBRC018_primer_filtered.fa -fastaout uniques_combined_merged_ITS1_2_GLBRC018_primer_filtered.fa -sizeout

---Fatal error---
fclose(0x232d690)=-1


#OTUS
01:00:27 85Mb     23.9% 8039 OTUs, 7837 chimeras
01:00:27 85Mb    100.0% 8039 OTUs, 7837 chimeras

```




## 7) Map reads back to OTUs at a 97% similarity score using [otutab](https://www.drive5.com/usearch/manual/cmd_otutab.html)
**-id 0.97 -strand plus are defaults**
### 7A) Mapping reads to traditional 0.97 OTUS

This step takes approximately 1 hour 30 minutes with ~30,000 OTUs




```

nano F_OTU_mapping_GLBRC018.sbatch

#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########

#SBATCH --time=165:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=3-5                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=5                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   F_mapping_OTU_GLBRC018  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/

usearch11.0.667 -otutab combined_merged_ITS1_2_GLBRC018.fastq -otus rep_set_combined_merged_ITS1_2_GLBRC_otus.fa  -otutabout combined_merged_ITS1_2_GLBRC_OTU_table.txt -notmatchedfq combined_merged_ITS1_2_GLBRC_otu_unmapped.fq


###end of .sbatch
sbatch F_OTU_mapping_GLBRC018.sbatch
#Submitted batch job 20619372
```
```
#Output from OTUs 
128045995 / 158659764 mapped to OTUs (80.7%)     
146:36:51 344Mb  Writing combined_merged_ITS1_2_GLBRC_OTU_table.txt
146:36:51 344Mb  Writing combined_merged_ITS1_2_GLBRC_OTU_table.txt ...done.
```



## 8) Classifying taxa against the reference database using [sintax](https://www.drive5.com/usearch/manual/cmd_sintax.html)


### 8A) Classifying the traditional 0.97 OTUs
### 8B) Classifying the ZOTUs

```
nano taxa_class_GLBRC018_OTU.sbatch

#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=10:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   F_taxa_class_GLBRC018_OTU  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu
cd /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/

usearch11.0.667 -sintax rep_set_combined_merged_ITS1_2_GLBRC_otus.fa -db /mnt/ufs18/home-087/belldere/Databases/UNITe8.2/sh_general_release_all_04.02.2020/utax_reference_dataset_all_04.02.2020_fix.udb -tabbedout taxonomy_combined_merged_ITS1_2_GLBRC_otus.sintax -strand both

scontrol show job $SLURM_JOB_ID     ### write job information to output file

###end of .sbatch

sbatch taxa_class_GLBRC018_OTU.sbatch
#Submitted batch job 20674315

 
```

```
#OTU Classification 
01:53 1.3Gb   100.0% Processing


```

### 9a) Classifiying the taxa with the native install and no Bonito lab paths and the utax_reference_dataset_all_04.02.2020_fix.udb [CONSTAX](https://github.com/Gian77/CONSTAXv2) against the eukaryote UNITe 8.2
```

cd /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/

#CONSTAX requires ".fasta" files

cp rep_set_combined_merged_ITS1_2_GLBRC_otus.fa CONSTAX_out/REP_set_combined_merged_ITS1_2_GLBRC_otus.fasta


nano CONSTAX_v2_GLBRC018_OTU_full.sbatch

#!/bin/bash --login

#SBATCH --time=12:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=32G                  # memory required per node - amount of memory (in bytes)
#SBATCH --job-name F_costax_v2_fix_All_v4.2.2020_T         # you can give your job a name for easier identification (same as -J)
#SBATCH --output=%x-%j.SLURMout
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/
export PATH=$PATH:$HOME/anaconda3/bin
source activate CONSTAXv2
/mnt/home/belldere/Programs/CONSTAX/CONSTAXv2/constax.sh --input CONSTAX_out/REP_set_combined_merged_ITS1_2_GLBRC_otus.fasta --db /mnt/ufs18/home-087/belldere/Databases/UNITe8.2/sh_general_release_all_04.02.2020/sh_general_release_dynamic_all_04.02.2020_fix.fasta --trainfile /mnt/home/belldere/Databases/UNITe8.2/sh_general_release_all_04.02.2020/CONSTAX_v2_Nav_training_files/ --tax CONSTAX_out/constax_V2_classification_all_v4.2.2020_taxa/ -o CONSTAX_out/constax_V2_classification_all_v4.2.2020_taxa/ --conf 0.8 -b --pathfile /mnt/home/belldere/Programs/CONSTAX/CONSTAXv2/pathfile.txt 

scontrol show job $SLURM_JOB_ID     ### write job information to output file
####

sbatch CONSTAX_v2_GLBRC018_OTU_full.sbatch

#Submitted batch job 20674014
```
```
RunTime=00:38:31
```

### 9b) Classifying the taxa with the native install and no Bonito lab paths and the most recent CONSTAX [CONSTAX](https://constax.readthedocs.io/en/latest/index.html) 
```

cd /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/

nano CONSTAX_v2.1_GLBRC018_OTU_full.sbatch

#!/bin/bash --login

#SBATCH --time=3:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=128G                  # memory required per node - amount of memory (in bytes)
#SBATCH --job-name F_costax_v2.1_fix_All_v4.2.2020_T         # you can give your job a name for easier identification (same as -J)
#SBATCH --output=%x-%j.SLURMout
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/

source activate CONSTAXv2
constax --input CONSTAX_out/REP_set_combined_merged_ITS1_2_GLBRC_otus.fasta --db /mnt/research/EvansLab/Databases/UNITE_8.2_fix/sh_general_release_dynamic_all_04.02.2020_fix.fasta --train --trainfile /mnt/research/EvansLab/Databases/UNITE_8.2_fix/CONSTAX_v2_Nav_training_files/ --tax constax_V2.1_classification_all_v4.2.2020_taxa/  -o constax_V2.1_classification_all_v4.2.2020_taxa/ --conf 0.8 -b 

scontrol show job $SLURM_JOB_ID     ### write job information to output file
####

sbatch CONSTAX_v2.1_GLBRC018_OTU_full.sbatch

#Submitted batch job 23948967
```
```
Run time 00:31:14
```


### TEST TEST CONSTAX 
```
cd /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/

nano TEST_CONSTAX_v2.1_GLBRC018_OTU_full.sbatch

#!/bin/bash --login

#SBATCH --time=3:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=128G                  # memory required per node - amount of memory (in bytes)
#SBATCH --job-name TEST_F_costax_v2.1_fix_All_v4.2.2020_T         # you can give your job a name for easier identification (same as -J)
#SBATCH --output=%x-%j.SLURMout
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /mnt/research/EvansLab/Lukas/GLBRC018_Root_Soil/Fungal_libs/

source activate CONSTAXv2
constax --input CONSTAX_out/REP_set_combined_merged_ITS1_2_GLBRC_otus.fasta --db /mnt/ufs18/home-087/belldere/Databases/UNITe8.2/NICO_Vers/sh_general_release_eukaryotes_91074_RepS_04.02.2020.fasta --train --trainfile /mnt/ufs18/home-087/belldere/Databases/UNITe8.2/NICO_Vers/CONSTAX_v2_Nav_training_files/ --tax TEST_constax_V2.1_classification_all_v4.2.2020_taxa/  -o TEST_constax_V2.1_classification_all_v4.2.2020_taxa/ --conf 0.8 -b 

scontrol show job $SLURM_JOB_ID     ### write job information to output file
####

sbatch TEST_CONSTAX_v2.1_GLBRC018_OTU_full.sbatch

#Submitted batch job 24610279
```