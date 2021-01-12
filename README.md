### Description

Pipeline for assembling large-insert clones (e.g. BACs, Fosmids) from Illumina data.

Software and documentation written by Daniel W. Bellott. For further discussion see:

Bellott, Cho, et al. "[Cost-effective, high-throughput, single-haplotype iterative mapping and sequencing for complex genomic structures](http://dx.doi.org/10.1038/nprot.2018.019)", *Nature Protocols* **13** 787-809 (2018).

### Prerequisites

Install the following tools used by the pipeline:

- [cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html)	(version 1.18)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2) (version 2.3.4.1)
- [SAMtools](http://www.htslib.org/download/) (version 1.9)
- [flash](https://ccb.jhu.edu/software/FLASH/) (version 1.2.11)
- [SPAdes](http://cab.spbu.ru/files/release3.10.1/manual.html#sec2) (version 3.14.1)
- [blasr](https://github.com/PacificBiosciences/blasr/wiki/Blasr-Installation-Qs-&-As) (version 5.3.2)
- [BESST](https://github.com/ksahlin/BESST/blob/master/docs/INSTALL.md) (version 2.2.8)
- [Gap2seq](https://www.cs.helsinki.fi/u/lmsalmel/Gap2Seq/) (version 2.0)
- [blat](http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads) (version 34)
- [Consed](http://www.phrap.org/consed/consed.html#howToGet) (version 29.0)
- [GATK](https://github.com/broadinstitute/gatk/releases) (version 3.7-0-gcfedb67)
- [Picard](https://github.com/broadinstitute/picard/releases/tag/2.24.0) (version 2.2.8)

Your computer should already have perl installed

All the required perl modules are packaged in the /vendor/cache directory for your convenience.

### Installation

Navigate to the directory where you want to install the pipeline and type:

```
git clone
https://github.com/dwbellott/shims2_assembly_pipeline.git
```

You will be prompted for your username and password.

After the download completes, type:

```
cd shims2_assembly_pipeline
vendor/bin/carton install --cached --deployment
```

To install the cached perl modules.

### Usage

```
$ ./shims2.pl -h
SHIMS Pipeline
(Version: 1.3.2)


USAGE: /home/bellott/perlscripts/shims2.pl -1 <upstream mates> -2 <downstream mates> -o <output directory> [optional arguments]

-1 and -2 may be comma spearated lists of files containing reads in fastq or fastq.gz

Optional  arguments:

	About this program:
		-h --help     prints this message
		-v --version  print version information
		--executables print default path to executables

	Screening Contamination:
		--vector  <bowtie index of vector sequences>
		--host    <bowtie index of host cell genome>
		--adapter <Illumina adapter sequences>

	Increasing Contiguity:
		--pacbioroi <fastq sequences of pacbio reads of insert>
		--pacbiofsr <fastq sequences of pacbio filtered subreads>
		--finished  <fasta sequences of finished neighbors>
		--draft     <fasta sequences of draft neighbors>
		--left_end  <fasta sequences of clone ends>
		--right_end <fasta sequences of clone ends>
		--peptides  <fasta sequences of peptides>

	Changing Executables:
		--spades        <path to SPAdes: /home/bellott/bin/spades.py>
		--samtools      <path to samtools: /usr/local/bin/samtools>
		--bowtie2build  <path to bowtie2-build: /usr/bin/bowtie2-build>
		--bowtie2       <path to bowtie2: /usr/local/bin/bowtie2>
		--blat          <path to blat: /usr/bin/blat>
		--cutadapt      <path to cutadapt: /usr/local/bin/python3.6/cutadapt>
		--blasr         <path to blasr: /home/bellott/miniconda3/bin/blasr>
		--consed        <path to consed: /usr/local/genome/bin/consed>
		--makeregions   <path to consed's makeRegionsFile.perl: /usr/local/genome/bin/makeRegionsFile.perl>
		--flash         <path to flash: /home/bellott/bin/flash>
		--besst         <path to besst: /home/bellott/bin/BESST/runBESST>
		--gap2seq       <path to gap2seq: /usr/local/bin/Gap2Seq.sh>
		--picardtools   <path to picardtools: /usr/local/share/picard-tools/picard.jar>
		--gatk          <path to gatk: /usr/local/bin/gatk>
		--java          <path to java: /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java>
  ```
