##Let  extract the reads from 2018 (20180119_A_16S-V4_PE) because that has the merds reads are in it. This is considered our Run C for 2018.

```
cd /MMPRNT_2018/20190617_16S-V4_PE250/
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
!#/bin/bash
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
grep "^File" * > 20190617_16S-V4_PE250_fastq_lengths.txt
grep "^EE" * > 20190617_16S-V4_PE250_fastq_EE.txt
```

look for any forward and reverse reads that look especially bad in terms of quality (high E is bad quality). This info will also be really helpful for troubleshooting later on (e.g. why some samples have extremely low read numbers)

#### we need to remove the underscores in the MERDS and ZymoMock file names
cd /MMPRNT_2018/20190617_16S-V4_PE250/
rename "MERDS_SG_" "MERDSSG" MERDS_SG*
rename "ZymoMock_C" "ZymoMockC" ZymoMock_C*


##  2) Merge the forward and reverse sequences and trim adapters (for each run individually)
### 2a) Merge Pairs
#### Make sure you are in the folder with the extracted forward and reverse reads from Run 1
https://www.drive5.com/usearch/manual/merge_options.html
#### -alnout gives you a human readable text file of the alignment and misalignments for each pair merged.
#### -tabbedout give you extensive information on the quality of the merge
#### This step takes approximately 10 minutes
```

cd /Lukas/
mkdir MMPRNT018
cd MMPRNT018
mkdir 20190617_16S-V4_PE250
cd /MMPRNT_2018/20190617_16S-V4_PE250/

usearch10.0.240 -fastq_mergepairs *R1*.fastq -relabel @ -fastqout /Lukas/MMPRNT018/20190617_16S-V4_PE250/combined_merged_20190617_16S-V4_PE250.fastq  -tabbedout /Lukas/MMPRNT018/20190617_16S-V4_PE250/combined_merged_20190617_16S-V4_PE250_pair_report.txt -alnout /Lukas/MMPRNT018/20190617_16S-V4_PE250/combined_merged_20190617_16S-V4_PE250_pair_aln.txt

```
09:09 722Mb   100.0% 61.8% merged

Totals:
  10157906  Pairs (10.2M)
   6280731  Merged (6.3M, 61.83%)
   2029432  Alignments with zero diffs (19.98%)
   3825172  Too many diffs (> 5) (37.66%)
         0  Fwd tails Q <= 2 trimmed (0.00%)
        16  Rev tails Q <= 2 trimmed (0.00%)
     52003  No alignment found (0.51%)
         0  Alignment too short (< 16) (0.00%)
    103138  Staggered pairs (1.02%) merged & trimmed
    246.56  Mean alignment length
    252.87  Mean merged length
      0.34  Mean fwd expected errors
      1.46  Mean rev expected errors
      0.08  Mean merged expected errors

```


###Reverse sequences look bad quality
Lets trim the reverse by 80bp tail of the reverse read by 80bp because the forward read alignment to the reverse read tail was poor and the 505F to 815R has a higher overlap (i.e. ~210bp overlap)

```

mkdir trim_reverse_20190617_16S-V4_PE250

nano rev_trunc.sh

for fq in *R2_001.fastq
do
  usearch10.0.240 -fastx_truncate $fq -stripright 80 -fastqout trim_reverse_20190617_16S-V4_PE250/$fq

done
```
end of executable file

make the file executable and run the code.

```
chmod +x rev_trunc.sh
./rev_trunc.sh
```

Now try to merge the sequences again with the truncated reverse reads.

```
usearch10.0.240 -fastq_mergepairs *R1*.fastq -reverse trim_reverse_20190617_16S-V4_PE250/*R2* -relabel @ -fastqout /Lukas/MMPRNT018/20190617_16S-V4_PE250/combined_merged_20190617_16S-V4_PE250_trunc.fastq  -tabbedout /Lukas/MMPRNT018/20190617_16S-V4_PE250/combined_merged_20190617_16S-V4_PE250_trunc_pair_report.txt -alnout /Lukas/MMPRNT018/20190617_16S-V4_PE250/combined_merged_20190617_16S-V4_PE250_trunc_pair_aln.txt
```
Good merging 
```
07:43 722Mb   100.0% 78.0% merged

Totals:
  10157906  Pairs (10.2M)
   7925163  Merged (7.9M, 78.02%)
   4074894  Alignments with zero diffs (40.12%)
   2167067  Too many diffs (> 5) (21.33%)
     65676  No alignment found (0.65%)
         0  Alignment too short (< 16) (0.00%)
       820  Staggered pairs (0.01%) merged & trimmed
    166.87  Mean alignment length
    252.87  Mean merged length
      0.45  Mean fwd expected errors
      0.74  Mean rev expected errors
      0.14  Mean merged expected errors

```


### 2b) let's check sequence quality of the merged seqs using USEARCH's [fastq_eestats2](https://www.drive5.com/usearch/manual/cmd_fastq_eestats2.html).

```
usearch11.0.667 -fastq_eestats2 /Lukas/MMPRNT018/20190617_16S-V4_PE250/combined_merged_20190617_16S-V4_PE250_trunc.fastq -output /Lukas/MMPRNT018/20190617_16S-V4_PE250/combined_merged_20190617_16S-V4_PE250_trunc_eestats2.txt
```
Quality and length look good
```
7925163 reads, max len 404, avg 252.9

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    7913439( 99.9%)    7924778(100.0%)    7924987(100.0%)
   100    7794047( 98.3%)    7915961( 99.9%)    7924830(100.0%)
   150    7685716( 97.0%)    7894640( 99.6%)    7924418(100.0%)
   200    7579552( 95.6%)    7871484( 99.3%)    7922681(100.0%)
   250    7316608( 92.3%)    7736598( 97.6%)    7841017( 98.9%)
   300       4813(  0.1%)       6015(  0.1%)       6527(  0.1%)
   350        478(  0.0%)        851(  0.0%)       1147(  0.0%)
   400          9(  0.0%)         23(  0.0%)         40(  0.0%)
```



### 2c) Now remove any residual bases from adapter seqs using [cut adapt](http://cutadapt.readthedocs.io/en/stable/index.html)
This code removes the forward adapter of 515F and the reverse complement of 806R. If you’re using the hpcc, load cutadapt.

```
#updated module for cutadapt (necessary with the update to the HPCC)
module load icc/2018.1.163-GCC-6.4.0-2.28
module load impi/2018.1.163
module load cutadapt/1.16-Python-3.6.4
cutadapt -a ATTAGAWACCCBDGTAGTCC -a GTGCCAGCMGCCGCGGTAA -o /Lukas/MMPRNT018/20190617_16S-V4_PE250/cut_combined_merged_20190617_16S-V4_PE250_trunc.fastq /Lukas/MMPRNT018/20190617_16S-V4_PE250/combined_merged_20190617_16S-V4_PE250_trunc.fastq > /Lukas/MMPRNT018/20190617_16S-V4_PE250/cut_adpt_results_combined_merged_20190617_16S-V4_PE250_trunc.txt
```
```
=== Summary ===

Total reads processed:               7,925,163
Reads with adapters:                     5,370 (0.1%)
Reads written (passing filters):     7,925,163 (100.0%)

Total basepairs processed: 2,004,054,748 bp
Total written (filtered):  2,004,027,165 bp (100.0%)
```

###Rename the merged file to get rid of hyphens in sample names (e.g. "MMPRNT-#" should be "MMPRNT#")

```
sed -e 's/MMPRNT-/MMPRNT/g' /Lukas/MMPRNT018/20190617_16S-V4_PE250/cut_combined_merged_20190617_16S-V4_PE250_trunc.fastq > /Lukas/MMPRNT018/20190617_16S-V4_PE250/renamed_cut_combined_merged_20190617_16S-V4_PE250_trunc.fastq
sed -e 's/Freeze-1123/Freeze1123/g' /Lukas/MMPRNT018/20190617_16S-V4_PE250/renamed_cut_combined_merged_20190617_16S-V4_PE250_trunc.fastq > /Lukas/MMPRNT018/20190617_16S-V4_PE250/renamed1_cut_combined_merged_20190617_16S-V4_PE250_trunc.fastq
sed -e 's/LC2017G5/LC2017G5C/g' /Lukas/MMPRNT018/20190617_16S-V4_PE250/renamed1_cut_combined_merged_20190617_16S-V4_PE250_trunc.fastq > /Lukas/MMPRNT018/20190617_16S-V4_PE250/renamed2_cut_combined_merged_20190617_16S-V4_PE250_trunc.fastq
sed -e 's/LUX2017G5/LUX2017G5C/g' /Lukas/MMPRNT018/20190617_16S-V4_PE250/renamed2_cut_combined_merged_20190617_16S-V4_PE250_trunc.fastq > /Lukas/MMPRNT018/20190617_16S-V4_PE250/renamed3_cut_combined_merged_20190617_16S-V4_PE250_trunc.fastq
sed -e 's/Marshall-/Marshall/g' /Lukas/MMPRNT018/20190617_16S-V4_PE250/renamed3_cut_combined_merged_20190617_16S-V4_PE250_trunc.fastq > /Lukas/MMPRNT018/20190617_16S-V4_PE250/renamed_cut_combined_merged_20190617_16S-V4_PE250_trunc.fastq
```


```
###Check sample names to make sure hyphen is gone:

```
less /Lukas/MMPRNT018/20190617_16S-V4_PE250/renamed_cut_combined_merged_20190617_16S-V4_PE250_trunc.fastq
```
Get sample names to check that all of the samples are present
```
usearch10.0.240 -fastx_get_sample_names /Lukas/MMPRNT018/20190617_16S-V4_PE250/renamed_cut_combined_merged_20190617_16S-V4_PE250_trunc.fastq -output /Lukas/MMPRNT018/20190617_16S-V4_PE250/renamed_cut_combined_merged_20190617_16S-V4_PE250_trunc_samples.txt
```
```
00:30 37Mb    100.0% 317 samples found
```
Looks good


##Let  extract the reads from 2018 run 20190614_16S-V4_PE250 this is considered our Run A

```
cd /MMPRNT_2018/20190614_16S-V4_PE250/
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
!#/bin/bash
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
grep "^File" * > 20190614_16S-V4_PE250_fastq_lengths.txt
grep "^EE" * > 20190614_16S-V4_PE250_fastq_EE.txt
```

look for any forward and reverse reads that look especially bad in terms of quality (high E is bad quality). This info will also be really helpful for troubleshooting later on (e.g. why some samples have extremely low read numbers)

#### we need to remove the underscores in the LC2017G5,LUX2017G5, and ZymoMock file names
cd /MMPRNT_2018/20190614_16S-V4_PE250/
rename "LC2017G5_A" "LC2017G5A" LC2017G5_A*
rename "LUX2017G5_A" "LUX2017G5A" LUX2017G5_A*
rename "ZymoMock_A" "ZymoMockA" ZymoMock_A*


##  2) Merge the forward and reverse sequences and trim adapters (for each run individually)
### 2a) Merge Pairs
#### Make sure you are in the folder with the extracted forward and reverse reads from Run 1
https://www.drive5.com/usearch/manual/merge_options.html
#### -alnout gives you a human readable text file of the alignment and misalignments for each pair merged.
#### -tabbedout give you extensive information on the quality of the merge
#### This step takes approximately 10 minutes
```

cd /Lukas/MMPRNT018

mkdir 20190614_16S-V4_PE250
cd /MMPRNT_2018/20190614_16S-V4_PE250/

usearch10.0.240 -fastq_mergepairs *R1*.fastq -relabel @ -fastqout /Lukas/MMPRNT018/20190614_16S-V4_PE250/combined_20190614_16S-V4_PE250.fastq  -tabbedout /Lukas/MMPRNT018/20190614_16S-V4_PE250/combined_merged_20190614_16S-V4_PE250_pair_report.txt -alnout /Lukas/MMPRNT018/20190614_16S-V4_PE250/combined_merged_20190614_16S-V4_PE250_pair_aln.txt

```
07:42 722Mb   100.0% 72.3% merged

Totals:
   9834162  Pairs (9.8M)
   7113261  Merged (7.1M, 72.33%)
   3170291  Alignments with zero diffs (32.24%)
   2687234  Too many diffs (> 5) (27.33%)
     33667  No alignment found (0.34%)
         0  Alignment too short (< 16) (0.00%)
    110918  Staggered pairs (1.13%) merged & trimmed
    246.51  Mean alignment length
    252.88  Mean merged length
      0.36  Mean fwd expected errors
      1.24  Mean rev expected errors
      0.08  Mean merged expected errors


```


Looks like an okay merge



### 2b) let's check sequence quality of the merged seqs using USEARCH's [fastq_eestats2](https://www.drive5.com/usearch/manual/cmd_fastq_eestats2.html).

```
usearch11.0.667 -fastq_eestats2 /Lukas/MMPRNT018/20190614_16S-V4_PE250/combined_20190614_16S-V4_PE250.fastq -output /Lukas/MMPRNT018/20190614_16S-V4_PE250/combined_20190614_16S-V4_PE250_eestats2.txt
```
Quality and length look good
```
00:58 661Mb   100.0% Reading reads

7113261 reads, max len 484, avg 252.9

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    7087230( 99.6%)    7112440(100.0%)    7113071(100.0%)
   100    7044444( 99.0%)    7109672( 99.9%)    7112964(100.0%)
   150    7010254( 98.6%)    7104742( 99.9%)    7112690(100.0%)
   200    6977737( 98.1%)    7099657( 99.8%)    7111113(100.0%)
   250    6839665( 96.2%)    7006470( 98.5%)    7027783( 98.8%)
   300       6854(  0.1%)       8348(  0.1%)       8903(  0.1%)
   350       1948(  0.0%)       3348(  0.0%)       4234(  0.1%)
   400        713(  0.0%)       1779(  0.0%)       2828(  0.0%)
   450         10(  0.0%)         35(  0.0%)         80(  0.0%)

```



### 2c) Now remove any residual bases from adapter seqs using [cut adapt](http://cutadapt.readthedocs.io/en/stable/index.html)
This code removes the forward adapter of 515F and the reverse complement of 806R. If you’re using the hpcc, load cutadapt.
Probably not super necessary
```
#updated module for cutadapt (necessary with the 2018 update to the HPCC)
module load icc/2018.1.163-GCC-6.4.0-2.28
module load impi/2018.1.163
module load cutadapt/1.16-Python-3.6.4
cutadapt -a ATTAGAWACCCBDGTAGTCC -a GTGCCAGCMGCCGCGGTAA -o /Lukas/MMPRNT018/20190614_16S-V4_PE250/cut_combined_20190614_16S-V4_PE250.fastq /Lukas/MMPRNT018/20190614_16S-V4_PE250/combined_20190614_16S-V4_PE250.fastq > /Lukas/MMPRNT018/20190614_16S-V4_PE250/cut_adpt_results_combined_20190614_16S-V4_PE250.txt
```
```
=== Summary ===

Total reads processed:               7,113,261
Reads with adapters:                     2,678 (0.0%)
Reads written (passing filters):     7,113,261 (100.0%)

Total basepairs processed: 1,798,809,266 bp
Total written (filtered):  1,798,788,269 bp (100.0%)

```

###Rename the merged file to get rid of hyphens in sample names (e.g. "MMPRNT-#" should be "MMPRNT#")

```
sed -e 's/MMPRNT-/MMPRNT/g' /Lukas/MMPRNT018/20190614_16S-V4_PE250/cut_combined_20190614_16S-V4_PE250.fastq > /Lukas/MMPRNT018/20190614_16S-V4_PE250/renamed_cut_combined_20190614_16S-V4_PE250.fastq

```


```
###Check sample names to make sure hyphen is gone:

```
less /Lukas/MMPRNT018/20190614_16S-V4_PE250/renamed_cut_combined_20190614_16S-V4_PE250.fastq
```
Get sample names to check that all of the samples are present
```
usearch10.0.240 -fastx_get_sample_names /Lukas/MMPRNT018/20190614_16S-V4_PE250/renamed_cut_combined_20190614_16S-V4_PE250.fastq -output /Lukas/MMPRNT018/20190614_16S-V4_PE250/renamed_cut_combined_20190614_16S-V4_PE250_samples.txt
```
```
00:25 37Mb    100.0% 315 samples found
```
Looks good 


##Let  extract the reads from 2018 run 20190621_16S-V4_PE250 this is considered our Run B

```
cd /MMPRNT_2018/20190621_16S-V4_PE250/
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
grep "^File" * > 20190614_16S-V4_PE250_fastq_lengths.txt
grep "^EE" * > 20190614_16S-V4_PE250_fastq_EE.txt
```

look for any forward and reverse reads that look especially bad in terms of quality (high E is bad quality). This info will also be really helpful for troubleshooting later on (e.g. why some samples have extremely low read numbers)

#### we need to remove the underscores in the LC2017G5,LUX2017G5, and ZymoMock file names
cd /MMPRNT_2018/20190621_16S-V4_PE250/
rename "LC2017G5_B" "LC2017G5B" LC2017G5_B*
rename "LUX2017G5_B" "LUX2017G5B" LUX2017G5_B*
rename "ZymoMock_B" "ZymoMockB" ZymoMock_B*
rename "CTRL-" "CTRL" CTRL-*
rename "error-1584-" "error1584" error-1584-*
rename "Freeze-" "Freeze" Freeze-*

##  2) Merge the forward and reverse sequences and trim adapters (for each run individually)
### 2a) Merge Pairs
#### Make sure you are in the folder with the extracted forward and reverse reads from Run 1
https://www.drive5.com/usearch/manual/merge_options.html
#### -alnout gives you a human readable text file of the alignment and misalignments for each pair merged.
#### -tabbedout give you extensive information on the quality of the merge
#### This step takes approximately 10 minutes
```

cd /Lukas/MMPRNT018

mkdir 20190621_16S-V4_PE250
cd /MMPRNT_2018/20190621_16S-V4_PE250/

usearch10.0.240 -fastq_mergepairs *R1*.fastq -relabel @ -fastqout /Lukas/MMPRNT018/20190621_16S-V4_PE250/combined_20190621_16S-V4_PE250.fastq  -tabbedout /Lukas/MMPRNT018/20190621_16S-V4_PE250/combined_merged_20190621_16S-V4_PE250_pair_report.txt -alnout /Lukas/MMPRNT018/20190621_16S-V4_PE250/combined_merged_20190621_16S-V4_PE250_pair_aln.txt

```
07:16 307Mb   100.0% 48.7% merged

Totals:
  12547934  Pairs (12.5M)
   6110497  Merged (6.1M, 48.70%)
   1665922  Alignments with zero diffs (13.28%)
   6352682  Too many diffs (> 5) (50.63%)
         0  Fwd tails Q <= 2 trimmed (0.00%)
       616  Rev tails Q <= 2 trimmed (0.00%)
     84755  No alignment found (0.68%)
         0  Alignment too short (< 16) (0.00%)
    169194  Staggered pairs (1.35%) merged & trimmed
    246.29  Mean alignment length
    252.64  Mean merged length
      0.79  Mean fwd expected errors
      1.45  Mean rev expected errors
      0.08  Mean merged expected errors


```


Looks like bad merge

###Reverse sequences look bad quality
Lets trim the reverse by 80bp tail of the reverse read by 80bp because the forward read alignment to the reverse read tail was poor and the 505F to 815R has a higher overlap (i.e. ~210bp overlap)

```

mkdir trim_reverse_20190621_16S-V4_PE250

nano rev_trunc.sh

for fq in *R2_001.fastq
do
  usearch10.0.240 -fastx_truncate $fq -stripright 80 -fastqout trim_reverse_20190621_16S-V4_PE250/$fq

done
```
end of executable file

make the file executable and run the code.

```
chmod +x rev_trunc.sh
./rev_trunc.sh
```

Now try to merge the sequences again with the truncated reverse reads.

```
usearch10.0.240 -fastq_mergepairs *R1*.fastq -reverse trim_reverse_20190621_16S-V4_PE250/*R2* -relabel @ -fastqout /Lukas/MMPRNT018/20190621_16S-V4_PE250/combined_merged_20190621_16S-V4_PE250_trunc.fastq  -tabbedout /Lukas/MMPRNT018/20190621_16S-V4_PE250/combined_merged_20190621_16S-V4_PE250_trunc_pair_report.txt -alnout /Lukas/MMPRNT018/20190621_16S-V4_PE250/combined_merged_20190621_16S-V4_PE250_trunc_pair_aln.txt
```
Okay merging 
```
05:42 307Mb   100.0% 61.5% merged

Totals:
  12547934  Pairs (12.5M)
   7712456  Merged (7.7M, 61.46%)
   3148981  Alignments with zero diffs (25.10%)
   4716629  Too many diffs (> 5) (37.59%)
    118849  No alignment found (0.95%)
         0  Alignment too short (< 16) (0.00%)
      1158  Staggered pairs (0.01%) merged & trimmed
    166.84  Mean alignment length
    252.69  Mean merged length
      0.92  Mean fwd expected errors
      0.74  Mean rev expected errors
      0.14  Mean merged expected errors


```



### 2b) let's check sequence quality of the merged seqs using USEARCH's [fastq_eestats2](https://www.drive5.com/usearch/manual/cmd_fastq_eestats2.html).

```
usearch11.0.667 -fastq_eestats2 /Lukas/MMPRNT018/20190621_16S-V4_PE250/combined_merged_20190621_16S-V4_PE250_trunc.fastq -output /Lukas/MMPRNT018/20190621_16S-V4_PE250/combined_merged_20190621_16S-V4_PE250_trunc_eestats2.txt
```
Quality and length look good
```
7712456 reads, max len 404, avg 252.7

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    7710498(100.0%)    7712054(100.0%)    7712093(100.0%)
   100    7629849( 98.9%)    7707968( 99.9%)    7712004(100.0%)
   150    7511884( 97.4%)    7689663( 99.7%)    7711554(100.0%)
   200    7382833( 95.7%)    7669129( 99.4%)    7709805(100.0%)
   250    7032008( 91.2%)    7484008( 97.0%)    7579310( 98.3%)
   300       6040(  0.1%)       7791(  0.1%)       8703(  0.1%)
   350        433(  0.0%)        976(  0.0%)       1494(  0.0%)
   400         20(  0.0%)         69(  0.0%)        125(  0.0%)


```



### 2c) Now remove any residual bases from adapter seqs using [cut adapt](http://cutadapt.readthedocs.io/en/stable/index.html)
This code removes the forward adapter of 515F and the reverse complement of 806R. If you’re using the hpcc, load cutadapt.
Probably not super necessary
```
#updated module for cutadapt (necessary with the 2018 update to the HPCC)
module load icc/2018.1.163-GCC-6.4.0-2.28
module load impi/2018.1.163
module load cutadapt/1.16-Python-3.6.4
cutadapt -a ATTAGAWACCCBDGTAGTCC -a GTGCCAGCMGCCGCGGTAA -o /Lukas/MMPRNT018/20190621_16S-V4_PE250/cut_combined_merged_20190621_16S-V4_PE250_trunc.fastq /Lukas/MMPRNT018/20190621_16S-V4_PE250/combined_merged_20190621_16S-V4_PE250_trunc.fastq > /Lukas/MMPRNT018/20190621_16S-V4_PE250/cut_adpt_results_combined_merged_20190621_16S-V4_PE250_trunc.txt
```
```
=== Summary ===

Total reads processed:               7,712,456
Reads with adapters:                     4,432 (0.1%)
Reads written (passing filters):     7,712,456 (100.0%)

Total basepairs processed: 1,948,838,472 bp
Total written (filtered):  1,948,802,785 bp (100.0%)


```

###Rename the merged file to get rid of hyphens in sample names (e.g. "MMPRNT-#" should be "MMPRNT#")

```
sed -e 's/MMPRNT-/MMPRNT/g' /Lukas/MMPRNT018/20190621_16S-V4_PE250/cut_combined_merged_20190621_16S-V4_PE250_trunc.fastq > /Lukas/MMPRNT018/20190621_16S-V4_PE250/rename_cut_combined_merged_20190621_16S-V4_PE250_trunc.fastq

```


```
###Check sample names to make sure hyphen is gone:

```
less /Lukas/MMPRNT018/20190621_16S-V4_PE250/rename_cut_combined_merged_20190621_16S-V4_PE250_trunc.fastq
```
Get sample names to check that all of the samples are present
```
usearch10.0.240 -fastx_get_sample_names /Lukas/MMPRNT018/20190621_16S-V4_PE250/rename_cut_combined_merged_20190621_16S-V4_PE250_trunc.fastq -output /Lukas/MMPRNT018/20190621_16S-V4_PE250/rename_cut_combined_merged_20190621_16S-V4_PE250_trunc_samples.txt
```
```
00:25 37Mb    100.0% 315 samples found
```
Looks good 


#Bonito Lab Root Field Libraries

## Let  demultiplex the reads from Lib_1

```
#29377783 seqs
#cd /Lukas
#mkdir GLBRC018_Root_Soil
```


+ /GLBRC_SG_Roots/Field2018_Bac/Lib_1_backup_Bac/

Pedro created a bar code file in a text editor

### Let's use the QIIME2 demux
https://docs.qiime2.org/2020.2/plugins/available/demux/emp-paired/?highlight=demux

###### Import the multiplexed sequences
We need to rename and move the zipped sequence files
```
cd /GLBRC_SG_Roots/Field2018_Bac/Lib_1_backup_Bac/
mkdir QIIME2_demux
cd QIIME2_demux/
mkdir compressed_seqs
cd ..

cp GLBRC_Lib1_16s_root_2018_r_S1_L001_R1_001.fastq.gz  QIIME2_demux/compressed_seqs/forward.fastq.gz
cp GLBRC_Lib1_16s_root_2018_r_S1_L001_R2_001.fastq.gz  QIIME2_demux/compressed_seqs/reverse.fastq.gz
cp GLBRC_Lib1_16s_root_2018_r_S1_L001_I1_001.fastq.gz  QIIME2_demux/compressed_seqs/barcodes.fastq.gz

#Now we need to convert the files to Q2 format and demux 

#qiime tools import --show-importable-types


cd QIIME2_demux

nano Q2_demux_Lib1_GLBRC018_root_bact.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=10:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name Q2_demux_Lib1_GLBRC018_root      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd ${SLURM_SUBMIT_DIR}
export PATH="/Software/anaconda3/bin:$PATH"
source activate qiime2-2019.4

#convert the files to the qza format 
qiime tools import --type EMPPairedEndSequences --input-path compressed_seqs --output-path multiplexed-seqs_paired_end.qza


#demux the sequences
qiime demux emp-paired --m-barcodes-file Lib_1_mapping.txt --m-barcodes-column BarcodeSequence --p-rev-comp-mapping-barcodes --p-no-golay-error-correction --i-seqs multiplexed-seqs_paired_end.qza --o-per-sample-sequences Lib1_demultiplexed-seqs.qza --o-error-correction-details Lib1_demux-details.qza

#Let's also export the demux files

qiime tools export  --input-path Lib1_demultiplexed-seqs.qza  --output-path output/demultiplexed-seqs

#Let's try to visualize the demux sequences too

qiime demux summarize  --i-data Lib1_demultiplexed-seqs.qza  --o-visualization demux-seq_summary_Lib1.qzv

###End

sbatch Q2_demux_Lib1_GLBRC018_root_bact.sbatch
#Submitted batch job 20549085
#Run time 
```

Demux stats

```


```

Unzip the demultiplexed sequences

```

cd /GLBRC_SG_Roots/Field2018_Bac/Lib_1_backup_Bac/QIIME2_demux/output/demultiplexed-seqs
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
grep "^File" * > Lib1_Amplicon_PE300_fastq_lengths.txt
grep "^EE" * > Lib1_Amplicon_PE300_fastq_EE.txt
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

cd /Lukas/GLBRC018_Root_Soil
#mkdir Bacterial_libs
mkdir Bacterial_libs/Lib1_bact
cd Bacterial_libs/Lib1_bact

#I am going to submit a job

nano merg_seq_Lib1_GLBRC018_soil_bact.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=1:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name merge_Lib1_GLBRC018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /GLBRC_SG_Roots/Field2018_Bac/Lib_1_backup_Bac/QIIME2_demux/output/demultiplexed-seqs

usearch11.0.667 -fastq_mergepairs *_R1*.fastq -relabel @ -fastqout /Lukas/GLBRC018_Root_Soil/Bacterial_libs/Lib1_bact/merged_Lib1_16S.fastq  -tabbedout /Lukas/GLBRC018_Root_Soil/Bacterial_libs/Lib1_bact/merged_Lib1_16S_pair_report.txt -fastqout_notmerged_fwd /Lukas/GLBRC018_Root_Soil/Bacterial_libs/Lib1_bact/NO_Match_16S_Amplicon_PE300_ITS1_for.fastq -fastqout_notmerged_rev /Lukas/GLBRC018_Root_Soil/Bacterial_libs/Lib1_bact/NO_Match_16S_Amplicon_PE300_ITS1_rev.fastq

###End

sbatch merg_seq_Lib1_GLBRC018_soil_bact.sbatch

#Submitted batch job 20587238

```

Merge stats

```

09:36 1.4Gb   100.0% 89.8% merged

Totals:
  27744744  Pairs (27.7M)
  24905143  Merged (24.9M, 89.77%)
  15535015  Alignments with zero diffs (55.99%)
   2726787  Too many diffs (> 5) (9.83%)
         0  Fwd tails Q <= 2 trimmed (0.00%)
     57975  Rev tails Q <= 2 trimmed (0.21%)
    112814  No alignment found (0.41%)
         0  Alignment too short (< 16) (0.00%)
   4127629  Staggered pairs (14.88%) merged & trimmed
    255.57  Mean alignment length
    278.91  Mean merged length
      1.60  Mean fwd expected errors
      2.76  Mean rev expected errors
      0.04  Mean merged expected errors

```

###Check sample names to make sure all of the samples are present

```
usearch11.0.667 -fastx_get_sample_names /Lukas/GLBRC018_Root_Soil/Bacterial_libs/Lib1_bact/merged_Lib1_16S.fastq -output /Lukas/GLBRC018_Root_Soil/Bacterial_libs/Lib1_bact/merged_merged_Lib1_16S_samples.txt
```
```
#Should be 368 samples
#100.0% 368 samples found 

```




## Let  demultiplex the reads from Lib_2

```
#22209927 reads
#cd /Lukas
#mkdir GLBRC018_Root_Soil
```


+ /GLBRC_SG_Roots/Field2018_Bac/Lib_2_backup_Bac/

Pedro created a bar code file in a text editor

### Let's use the QIIME2 demux
https://docs.qiime2.org/2020.2/plugins/available/demux/emp-paired/?highlight=demux

###### Import the multiplexed sequences
We need to rename and move the zipped sequence files
```
cd /GLBRC_SG_Roots/Field2018_Bac/Lib_2_backup_Bac/
mkdir QIIME2_demux
cd QIIME2_demux/
mkdir compressed_seqs
cd ..

cp GLBRC_Lib2_16s_root_2018_S1_L001_R1_001.fastq.gz  QIIME2_demux/compressed_seqs/forward.fastq.gz
cp GLBRC_Lib2_16s_root_2018_S1_L001_R2_001.fastq.gz  QIIME2_demux/compressed_seqs/reverse.fastq.gz
cp GLBRC_Lib2_16s_root_2018_S1_L001_I1_001.fastq.gz  QIIME2_demux/compressed_seqs/barcodes.fastq.gz

#Now we need to convert the files to Q2 format and demux 

#qiime tools import --show-importable-types


cd QIIME2_demux

nano Q2_demux_Lib2_GLBRC018_root_bact.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=10:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name Q2_demux_Lib2_GLBRC018_root      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd ${SLURM_SUBMIT_DIR}
export PATH="/Software/anaconda3/bin:$PATH"
source activate qiime2-2019.4

#convert the files to the qza format 
qiime tools import --type EMPPairedEndSequences --input-path compressed_seqs --output-path multiplexed-seqs_paired_end.qza


#demux the sequences
qiime demux emp-paired --m-barcodes-file Lib_2_mapping.txt --m-barcodes-column BarcodeSequence --p-rev-comp-mapping-barcodes --p-no-golay-error-correction --i-seqs multiplexed-seqs_paired_end.qza --o-per-sample-sequences Lib2_demultiplexed-seqs.qza --o-error-correction-details Lib2_demux-details.qza

#Let's also export the demux files

qiime tools export  --input-path Lib2_demultiplexed-seqs.qza  --output-path output/demultiplexed-seqs

#Let's try to visualize the demux sequences too

qiime demux summarize  --i-data Lib2_demultiplexed-seqs.qza  --o-visualization demux-seq_summary_Lib2.qzv

###End

sbatch Q2_demux_Lib2_GLBRC018_root_bact.sbatch
#Submitted batch job 20549092
#Run time 03:47:32
```

Demux stats

---
Demultiplexed sequence counts summary
Minimum:	1
Median:	53909.0
Mean:	55093.02827763496
Maximum:	145773
Total:	21431188

---

Unzip the demultiplexed sequences

```

cd /GLBRC_SG_Roots/Field2018_Bac/Lib_2_backup_Bac/QIIME2_demux/output/demultiplexed-seqs
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
grep "^File" * > Lib2_Amplicon_PE300_fastq_lengths.txt
grep "^EE" * > Lib2_Amplicon_PE300_fastq_EE.txt
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

cd /Lukas/GLBRC018_Root_Soil
mkdir Bacterial_libs
mkdir Bacterial_libs/Lib2_bact
cd Bacterial_libs/Lib2_bact

#I am going to submit a job

nano merg_seq_Lib2_GLBRC018_soil_bact.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=1:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name merge_Lib2_GLBRC018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /GLBRC_SG_Roots/Field2018_Bac/Lib_2_backup_Bac/QIIME2_demux/output/demultiplexed-seqs

usearch11.0.667 -fastq_mergepairs *_R1*.fastq -relabel @ -fastqout /Lukas/GLBRC018_Root_Soil/Bacterial_libs/Lib2_bact/merged_Lib2_16S.fastq  -tabbedout /Lukas/GLBRC018_Root_Soil/Bacterial_libs/Lib2_bact/merged_Lib2_16S_pair_report.txt -fastqout_notmerged_fwd /Lukas/GLBRC018_Root_Soil/Bacterial_libs/Lib2_bact/NO_Match_16S_Amplicon_PE300_ITS1_for.fastq -fastqout_notmerged_rev /Lukas/GLBRC018_Root_Soil/Bacterial_libs/Lib2_bact/NO_Match_16S_Amplicon_PE300_ITS1_rev.fastq

###End

sbatch merg_seq_Lib2_GLBRC018_soil_bact.sbatch

#Submitted batch job 20577266


```

Merge stats

```
07:00 1.6Gb   100.0% 83.4% merged

Totals:
  21431188  Pairs (21.4M)
  17868153  Merged (17.9M, 83.37%)
   8541349  Alignments with zero diffs (39.85%)
   3484677  Too many diffs (> 5) (16.26%)
     78358  No alignment found (0.37%)
         0  Alignment too short (< 16) (0.00%)
   1846877  Staggered pairs (8.62%) merged & trimmed
    263.62  Mean alignment length
    288.70  Mean merged length
      1.48  Mean fwd expected errors
      2.88  Mean rev expected errors
      0.05  Mean merged expected errors
```

###Check sample names to make sure all of the samples are present

```
usearch11.0.667 -fastx_get_sample_names /Lukas/GLBRC018_Root_Soil/Bacterial_libs/Lib2_bact/merged_Lib2_16S.fastq -output /Lukas/GLBRC018_Root_Soil/Bacterial_libs/Lib2_bact/merged_merged_Lib2_16S_samples.txt
```
```
#Should be 408 samples 
#01:03 38Mb    100.0% 388 samples found
#Amp408 Amp420 Amp432 Amp444 Amp456 Amp468 Amp480 Amp723 Amp724 Amp725 Amp726 Amp727 Amp728 Amp785 Amp786 Amp787 Amp788 Amp789 Amp791 Amp792


```



###Combined the Bonito lab generated libs

```
cd  /Lukas/GLBRC018_Root_Soil/Bacterial_libs/ 
cat Lib2_bact/merged_Lib2_16S.fastq Lib1_bact/merged_Lib1_16S.fastq > merged_roots_comb_16s.fastq 

```

#### Let's remove the barcode codons based off of Pedros codes

```

nano cutA_comb_GLBRC018_root_bact.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=3:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name merge_comb_GLBRC018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

module load icc/2018.1.163-GCC-6.4.0-2.28
module load impi/2018.1.163
module load cutadapt/1.16-Python-3.6.4
cd  /Lukas/GLBRC018_Root_Soil/Bacterial_libs/ 

cutadapt -g GTGCCAGCMGCCGCGGTAA -a ATTAGAWACCCBDGTAGTCC -e 0.01 -f fastq -n 2 --discard-untrimmed --match-read-wildcards -o strip_merged_roots_comb_16s.fastq merged_roots_comb_16s.fastq 

#additional primer filtering 

cutadapt -g "CMGCCGCGGTAA;min_overlap=10" -a "ATTAGAWACCCBDGTAGTCC;min_overlap=10" -e 0.1 -f fastq -n 2 --match-read-wildcards -o stripped_merged_roots_comb_16s.fastq strip_merged_roots_comb_16s.fastq

####End

sbatch cutA_comb_GLBRC018_root_bact.sbatch
#Submitted batch job 20588604
```
Output from first primer removal
```
=== Summary ===

Total reads processed:              42,773,296
Reads with adapters:                42,443,688 (99.2%)
Reads written (passing filters):    42,443,688 (99.2%)

Total basepairs processed: 12,104,755,498 bp
Total written (filtered):  9,691,315,718 bp (80.1%)
```

I got an error with the second line of primer removal I think I need to a newer version of [cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html)

```



mv merged_strip_roots_comb_16s.fastq strip_merged_roots_comb_16s.fastq 

#install conda version of cutadapt
conda create -n cutadaptenv cutadapt

nano cutA2_comb_GLBRC018_root_bact.sbatch



#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=3:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name cut2_comb_GLBRC018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu


cd  /Lukas/GLBRC018_Root_Soil/Bacterial_libs/ 

#cutadapt -g GTGCCAGCMGCCGCGGTAA -a ATTAGAWACCCBDGTAGTCC -e 0.01 -f fastq -n 2 --discard-untrimmed --match-read-wildcards -o strip_merged_roots_comb_16s.fastq merged_roots_comb_16s.fastq 

#mv strip_merged_roots_comb_16s.fastq merged_strip_roots_comb_16s.fastq

source activate cutadaptenv
#additional primer filtering 

cutadapt -g "CMGCCGCGGTAA;min_overlap=10" -a "ATTAGAWACCCBDGTAGTCC;min_overlap=10" -e 0.1 -n 2 --match-read-wildcards -o stripped_merged_roots_comb_16s.fastq strip_merged_roots_comb_16s.fastq

conda deactivate
####End

sbatch cutA2_comb_GLBRC018_root_bact.sbatch
#Submitted batch job 20606598

```

Output for the second filtering 

```
=== Summary ===

Total reads processed:              42,443,688
Reads with adapters:                 1,909,486 (4.5%)
Reads written (passing filters):    42,443,688 (100.0%)

Total basepairs processed: 9,691,315,718 bp
Total written (filtered):  9,634,441,175 bp (99.4%)
```

## 3) Combine the seven merged sequence files 
### Make sure you do not use replicate sample names between the runs
### This step is only necessary if more than one Illumina run are being analyzed together.

```
cd  /Lukas/GLBRC018_Root_Soil/Bacterial_libs/ 


#merge the 2018 data 
cat stripped_merged_roots_comb_16s.fastq /Lukas/MMPRNT018/20190621_16S-V4_PE250/rename_cut_combined_merged_20190621_16S-V4_PE250_trunc.fastq /Lukas/MMPRNT018/20190614_16S-V4_PE250/renamed_cut_combined_20190614_16S-V4_PE250.fastq /Lukas/MMPRNT018/20190617_16S-V4_PE250/renamed_cut_combined_merged_20190617_16S-V4_PE250_trunc.fastq > combined_soil_roots_merged_16S_GLBRC018.fastq


```

Before we continue, you may want to check if the sample names are formatted correctly. USEARCH does some funny cutting during the merging step. Any hyphens or underscores can be problematic and you need to remove these (use sed command and merged_cut files)

Additionally, this is a good opportunity to double check that all of your samples merged and have unique IDs using [fastx_get_sample_names](https://www.drive5.com/usearch/manual/cmd_fastx_get_sample_names.html)

```
usearch11.0.667 -fastx_get_sample_names combined_soil_roots_merged_16S_GLBRC018.fastq -output combined_soil_roots_merged_16S_GLBRC018_samples.txt
```
```
#Shoud be 317 + 315 + 315 + 368 + 388 = 1703
#03:03 38Mb    100.0% 1703 samples found


```
Looks good

## 4) Filtering phix, quality, and truncate the merged seqs  to MaxEE and set length using [fastq_filter](https://www.drive5.com/usearch/manual/cmd_fastq_filter.html)

#### 250 bp is the expected overlaps with 515F and 806R
```

nano filter_comb_GLBRC018_root_bact.sbatch



#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=5:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name B_fil_comb_GLBRC018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd  /Lukas/GLBRC018_Root_Soil/Bacterial_libs/ 

usearch11.0.667 -fastq_filter combined_soil_roots_merged_16S_GLBRC018.fastq -fastq_maxee 0.1 -fastq_trunclen 250 -fastaout combined_soil_roots_merged_16S_GLBRC018_fil.fasta

usearch11.0.667 -filter_phix combined_soil_roots_merged_16S_GLBRC018_fil.fasta -output combined_soil_roots_merged_16S_GLBRC018_phix_fil.fasta  -alnout combined_soil_roots_merged_16S_GLBRC018_phix_hits.txt




sbatch filter_comb_GLBRC018_root_bact.sbatch

#Submitted batch job 20616131
#Submitted batch job 20640396

```
```
#Quality and length filtering
04:57 282Mb   100.0% Filtering, 81.6% passed
  65194568  Reads (65.2M)                   
   5864144  Discarded reads length < 250
   6137154  Discarded reads with expected errs > 0.10
  53193270  Filtered reads (53.2M, 81.6%)

#Phix filtering did not work correctly first time

04:57 523Mb   100.0% Filtering for phix, 1 hits (0.0%)

```

## 5) Filter so we only have unique sequences with [fastx_uniques](https://www.drive5.com/usearch/manual/cmd_fastx_uniques.html)



## 6) Cluster into OTUS and filter out singletons
There are two options here. **(A)** uses the traditional approach and clusters sequences into 0.97 identity cutoff OTUs. **(B)** uses unoise3 to identify ZOTUs.

### 6A) Cluster into 0.97 OTUs using UPARSE and [cluster_otus](https://www.drive5.com/usearch/manual/cmd_cluster_otus.html)
This step will also denovo chimera check and filter out singletons.You can remove single sequences prior to clustering but singletons are also removed at the OTU clustering step (cluster_otus filters out OTUs <2 and unoise3 filters ZOTUs <8)

This step takes approximately 30 minutes with ~330,000 unique sequences


```
cd  /Lukas/GLBRC018_Root_Soil/Bacterial_libs/ 
nano derep_OTU_B_clustering_GLBRC018.sbatch

#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=165:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=3-4                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=3                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   B_derep_OTU_GLBRC018  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu
cd  /Lukas/GLBRC018_Root_Soil/Bacterial_libs/ 

usearch11.0.667 -fastx_uniques combined_soil_roots_merged_16S_GLBRC018_phix_fil.fasta  -fastaout uniques_combined_soil_roots_merged_16S_GLBRC018_phix_fil.fa -sizeout

usearch11.0.667 -cluster_otus uniques_combined_soil_roots_merged_16S_GLBRC018_phix_fil.fa -otus rep_set_combined_soil_roots_merged_16S_GLBRC018_otus.fasta -uparseout combined_soil_roots_merged_16S_GLBRC018_otus_uparse.txt -relabel OTU



###end of .sbatch

sbatch derep_OTU_B_clustering_GLBRC018.sbatch
#Submitted batch job 20673830
```

```
#Uniques 

02:41 24.3Gb 53193269 seqs, 8867367 uniques, 6529709 singletons (73.6%)

#OTUS
01:01:04 628Mb   100.0% 48842 OTUs, 142059 chimeras


```




## 7) Map reads back to OTUs at a 97% similarity score using [otutab](https://www.drive5.com/usearch/manual/cmd_otutab.html)
**-id 0.97 -strand plus are defaults**
### 7A) Mapping reads to traditional 0.97 OTUS

This step takes approximately 1 hour 30 minutes with ~30,000 OTUs



```

cd  /Lukas/GLBRC018_Root_Soil/Bacterial_libs/ 
nano OTU_mapping_GLBRC018_16S.sbatch

#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=165:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   B_mapping_OTU_GLBRC018  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd  /Lukas/GLBRC018_Root_Soil/Bacterial_libs/ 

usearch11.0.667 -otutab combined_soil_roots_merged_16S_GLBRC018.fastq -otus rep_set_combined_soil_roots_merged_16S_GLBRC018_otus.fasta -otutabout combined_soil_roots_merged_16S_GLBRC018_OTU_table.txt -notmatchedfq combined_soil_roots_merged_16S_GLBRC018_otu_unmapped.fq





###end of .sbatch
sbatch OTU_mapping_GLBRC018_16S.sbatch
#Submitted batch job 20679882
```
```
#Output from OTUs 

06:44:17 1.1Gb   100.0% Searching, 85.3% matched
55627577 / 65194568 mapped to OTUs (85.3%)      
06:44:17 1.1Gb  Writing combined_soil_roots_merged_16S_GLBRC018_OTU_table.txt
06:44:17 1.1Gb  Writing combined_soil_roots_merged_16S_GLBRC018_OTU_table.txt ...done.

```



## 8) Classifying taxa against the reference database using [sintax](https://www.drive5.com/usearch/manual/cmd_sintax.html)
Currently used database is [silva version 123](https://www.drive5.com/usearch/manual/sintax_downloads.html).
I am also going to use RDP database. (https://drive5.com/sintax/rdp_16s_v16.fa.gz)
### 8A) Classifying the traditional 0.97 OTUs
This step takes approximately  18 hours and 20 minutes with ~30,000 OTUs against the Silva reference database, so you may want to submit a job for it.

```
cd  /Lukas/GLBRC018_Root_Soil/Bacterial_libs/ 
nano taxa_class_OTU_B_GLBRC018.sbatch

#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=65:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   B_taxa_class_OTU_GLBRC018  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd  /Lukas/GLBRC018_Root_Soil/Bacterial_libs/ 

usearch11.0.667 -sintax rep_set_combined_soil_roots_merged_16S_GLBRC018_otus.fasta -db /Databases/Silva132_release/silva_16s_v123.fa -tabbedout combined_soil_roots_merged_16S_GLBRC018_otus_silva123_taxonomy.sintax -strand both

usearch11.0.667 -sintax rep_set_combined_soil_roots_merged_16S_GLBRC018_otus.fasta -db /Databases/RDP_16S_v16_usearch/rdp_16s_v16.fa -tabbedout combined_soil_roots_merged_16S_GLBRC018_otus_RDP_v16_taxonomy.sintax -strand both

###end of .sbatch

sbatch taxa_class_OTU_B_GLBRC018.sbatch
#Submitted batch job 20679946
```
#For Silva v123 based classification
06:22:57 12.0Gb  100.0% Processing
#For RDP based classification
00:36 762Mb   100.0% Processing
```

### 9) Build tree using PASTA


````
cd /Lukas/GLBRC018_Root_Soil/Bacterial_libs/ 
nano OTU_pasta_tree_build_GLBRC018.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=96:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name tree_OTU_pasta_build_GLBRC018      # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=belldere@msu.edu

cd /Lukas/GLBRC018_Root_Soil/Bacterial_libs/ 

module load GNU/7.3.0-2.30  OpenMPI/3.1.1
module load Python/2.7.15
python /Software/pasta-code/pasta/run_pasta.py -i rep_set_combined_soil_roots_merged_16S_GLBRC018_otus.fasta -j tree_OTU_16S_GLBRC018 -o tree_pasta_OTU_16S_GLBRC018/
####end sbatch file


sbatch OTU_pasta_tree_build_GLBRC018.sbatch
#Submitted batch job 20680082

```
```
#Run time was 
Total time spent: 4304.00007606s
```



