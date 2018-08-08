# Customized protein database construction

## Introduction

Customprodbj is a Java-based tool for customized protein database construction:

  - Build a customized database based on single VCF file.
  - Build a customized database based on multiple VCF files from a sample.
  - Build a customized database based on multiple VCF files from multiple samples.

## Usage

Please download Customprodbj program from the release page: https://github.com/wenbostar/customprodbj/releases

java -jar customprodbj.jar 
```
 -d      mRNA fasta database
 -f      A file which includes multiple samples. This parameter is used to
         build a customized database for
 -h      Help
 -i      Annovar annotation result file. Multiple files are separated by
         ','.
 -o      Output folder
 -p1     The prefix of variant protein ID, default is VAR_
 -p2     The prefix of final output files, default is merge
 -r      Gene annotation data
 -ref    Output reference protein database file
 -t      Whether or not to add reference protein sequences to the output
         database file
 -v      Verbose
```


## Example

#### Build a customized database based on single VCF file

##### Step 1: Variant annotation using [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/):

```
perl table_annovar.pl test.vcf annovar_database/humandb/ -buildver hg19 -out out/test -protocol refGene -operation g -nastring . -vcfinput --thread 30 --maxgenethread 30 -polish

```

##### Step 2: Build customized protein database using Customprodbj:

```
java -jar customprodbj.jar -i test.hg19_multianno.txt -d annovar_database/humandb/hg19_refGeneMrna.fa -r annovar_database/humandb/hg19_refGene.txt -t -o out/
```

The input file "test.hg19_multianno.txt" is from ANNOVAR annotation result. The input files "annovar_database/humandb/hg19_refGeneMrna.fa" and "annovar_database/humandb/hg19_refGene.txt" are two files used by ANNOVAR. Before uses do variant annotation using ANNOVAR, users need to download these files using ANNOVAR. Please follow the instruction described here: [http://annovar.openbioinformatics.org/en/latest/user-guide/download/](http://annovar.openbioinformatics.org/en/latest/user-guide/download/).


#### Build a customized database based on multiple VCF files from a sample

##### Step 1: Variant annotation using [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/):

Perform variant annotation for each VCF file using the same method described in above section.

##### Step 2: Build customized protein database using Customprodbj:

```
java -jar customprodbj.jar -f input_variant_file_list.txt -d annovar_database/humandb/hg19_refGeneMrna.fa -r annovar_database/humandb/hg19_refGene.txt -t -o out/
```

The format of input_variant_file_list.txt looks like below:

```
sample	somatic	germline	rna	msi
sample1	s1a.hg19_multianno.txt	s1b.hg19_multianno.txt	s1c.hg19_multianno.txt	s1d.hg19_multianno.txt
```
Please note the columns are separated by "\t".


#### Build a customized database based on multiple VCF files from multiple samples

##### Step 1: Variant annotation using [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/):

Perform variant annotation for each VCF file using the same method described in above section.

##### Step 2: Build customized protein database using Customprodbj:

```
java -jar customprodbj.jar -f input_variant_file_list.txt -d annovar_database/humandb/hg19_refGeneMrna.fa -r annovar_database/humandb/hg19_refGene.txt -t -o out/
```

The format of input_variant_file_list.txt looks like below:

```
sample	somatic	germline	rna	msi
sample1	s1a.hg19_multianno.txt	s1b.hg19_multianno.txt	s1c.hg19_multianno.txt	s1d.hg19_multianno.txt
sample2	s2a.hg19_multianno.txt	s2b.hg19_multianno.txt	s2c.hg19_multianno.txt	s2d.hg19_multianno.txt
sample3	s3a.hg19_multianno.txt	s3b.hg19_multianno.txt	s3c.hg19_multianno.txt	s3d.hg19_multianno.txt
```


## Ouput

The final outputs consist of three files:

*-varInfo.txt: A table contains the detailed amino acid change information. Each row is a variant.

*-var.fasta: Customized protein database.

*-varStat.txt: Summary data.



