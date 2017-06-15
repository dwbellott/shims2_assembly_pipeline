### Prerequisites:

Install the following tools used by the pipeline:

- [cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2)
- [SAMtools](http://www.htslib.org/download/)
- [flash](https://ccb.jhu.edu/software/FLASH/)
- [SPAdes](http://cab.spbu.ru/files/release3.10.1/manual.html#sec2)
- [blasr](https://github.com/PacificBiosciences/blasr/wiki/Blasr-Installation-Qs-&-As)
- [BESST](https://github.com/ksahlin/BESST/blob/master/docs/INSTALL.md)
- [Gap2seq](https://www.cs.helsinki.fi/u/lmsalmel/Gap2Seq/)
- [blat](http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads)
- [Consed](http://www.phrap.org/consed/consed.html#howToGet)

If you use [Docker](https://www.docker.com/), you may be pleased to learn that several of these tools are available as [BioContainers](https://biocontainers.pro/#documentation):

- [cutadapt](https://quay.io/repository/biocontainers/cutadapt)
- [bowtie2](https://quay.io/repository/biocontainers/bowtie2)
- [SAMtools](https://quay.io/repository/biocontainers/samtools)
- [SPAdes](https://quay.io/repository/biocontainers/spades)
- [blasr](https://quay.io/repository/biocontainers/blasr)
- [BESST](https://quay.io/repository/biocontainers/BESST)
- [Gap2seq](https://quay.io/repository/biocontainers/gap2seq)
- [blat](https://quay.io/repository/biocontainers/blat)

Your computer should already have perl installed

All the required perl modules are packaged in the /vendor/cache directory for your convenience.

### Installation

Navigate to the directory where you want to install the pipeline and type:

```
git clone https://git.wi.mit.edu/winstons_projects/shims2_assembly_pipeline
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
(Version: 1.1)


USAGE: ./shims2.pl -1 <upstream mates> -2 <downstream mates> -o <output directory> [optional arguments]

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
		--pacbioroi <fastq sequences of pacbio reads of insert>>
		--pacbiofsr <fastq sequences of pacbio filtered subreads>
		--finished  <fasta sequences of finished neighbors>
		--draft     <fasta sequences of draft neighbors>
		--left_end  <fasta sequences of clone ends>
		--right_end <fasta sequences of clone ends>
		--peptides  <fasta sequences of peptides>

	Changing Executables:
		--spades        <path to SPAdes: >
		--samtools      <path to samtools: >
		--bowtie2build  <path to bowtie2-build: >
		--bowtie2       <path to bowtie2: >
		--blat          <path to blat: >
		--cutadapt      <path to cutadapt: >
		--blasr         <path to blasr: >
		--consed        <path to consed: >
		--makeregions   <path to consed's makeRegionsFile.perl: >
		--flash         <path to flash: >
		--besst         <path to besst: >
		--gap2seq       <path to gap2seq: >
  ```
