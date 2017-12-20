#!/usr/bin/env perl

############### Libraries
use FindBin qw($Bin);
use lib "$Bin/local/";

use strict;
use Cwd;
use Cwd 'abs_path';
use File::Which;
use File::Path 'rmtree';
use DateTime;
use Getopt::Long;
use Bio::Seq;
use Bio::Seq::SeqWithQuality;
use Bio::SeqIO;
use Statistics::Descriptive;
use IO::Compress::Gzip qw($GzipError);


use vars qw/$VERSION/;
use vars qw/$spades_default/;
use vars qw/$samtools_default/;
use vars qw/$bowtie2build_default/;
use vars qw/$bowtie2_default/;
use vars qw/$blat_default/;
use vars qw/$cutadapt_default/;
use vars qw/$blasr_default/;
use vars qw/$consed_default/;
use vars qw/$makeregions_default/;
use vars qw/$flash_default/;
use vars qw/$besst_default/;
use vars qw/$gap2seq_default/;

use vars qw/$spades_exec/;
use vars qw/$samtools_exec/;
use vars qw/$bowtie2build_exec/;
use vars qw/$bowtie2_exec/;
use vars qw/$blat_exec/;
use vars qw/$cutadapt_exec/;
use vars qw/$blasr_exec/;
use vars qw/$consed_exec/;
use vars qw/$makeregions_exec/;
use vars qw/$flash_exec/;
use vars qw/$besst_exec/;
use vars qw/$gap2seq_exec/;

use vars qw/$gap2seq_timeout/;

use vars qw/$bowtie2_screen_params/;
use vars qw/$cutadapt_params/;
use vars qw/$adapter_default/;
use vars qw/$phix_def/;
use vars qw/$keeptemp/;
use vars qw/$cleantemp/;
use vars qw/@temporary/;



BEGIN {
	$VERSION = '1.1.27';
	$spades_default = $ENV{'SHIMS_SPADES_EXEC'} || which('spades.py');
	$samtools_default = $ENV{'SHIMS_SAMTOOLS_EXEC'} || which('samtools');
	$bowtie2build_default = $ENV{'SHIMS_BOWTIE2BUILD_EXEC'} || which('bowtie2-build');
	$bowtie2_default = $ENV{'SHIMS_BOWTIE2_EXEC'} || which('bowtie2');
	$blat_default = $ENV{'SHIMS_BLAT_EXEC'} || which('blat');
	$cutadapt_default = $ENV{'SHIMS_CUTADAPT_EXEC'} || which('cutadapt');
	$blasr_default = $ENV{'SHIMS_BLASR_EXEC'} || which('blasr');
	$consed_default = $ENV{'SHIMS_CONSED_EXEC'} || which('consed');
	$makeregions_default = $ENV{'SHIMS_MAKEREGIONS_EXEC'} || which('makeRegionsFile.perl');
	$flash_default = $ENV{'SHIMS_FLASH_EXEC'} || which('flash');
	$besst_default = $ENV{'SHIMS_BESST_EXEC'} || which('runBESST');
	$gap2seq_default = $ENV{'SHIMS_GAP2SEQ_EXEC'} || which('Gap2Seq.sh');

	$gap2seq_timeout = 3600;

	$bowtie2_screen_params = "--very-sensitive-local --n-ceil L,0,1 -I 0 -X 2501";
	$cutadapt_params = "--mask-adapter --quiet --match-read-wildcards -q 10 --minimum-length 22";
	$adapter_default = "AGATCGGAAGAGC";
	$phix_def = "gi|216019|gb|J02482.1|PX1CG";
}

############### Declare Subroutines and Arguments
sub version ();
sub usage ();
sub executables ();
sub check_executable ($$$);
sub screen_contamination ($$$$$);
sub mask_adapters ($$$);
sub make_illumina_arguments ($$$);
sub merge_fastas($$$);
sub screen_pacbio ($$$$);
sub imagine ($$$$$$);
sub write_fasta_sequences ($$$);
sub quality_trim_contigs ($$);
sub load_fasta_contigs ($);
sub overlap_mates ($$$$$$);
sub estimate_coverage_cutoff_from_vector ($$$);
sub autofinish ($$$$$);
sub autofinish_status ($);
sub make_result_contig ($$$$$);
sub fasta_output_autofinished ($$);
sub main ();

############### Run the program

main();

############### Define Subroutines

sub main() {
	my $help = 1;
	my $date = DateTime->now->ymd;
	#declare the options
	my (@upstream_mates,
		@downstream_mates,
		@singles,
		$output_dir,$help,
		$version,
		$executables,
		$vector,
		$host,
		@adapter,
		@pacbioroi,
		@pacbiofsr,
		@finished,
		@draft,
		$left_end,
		$right_end,
		@peptides,
		$spades,
		$samtools,
		$bowtie2build,
		$bowtie2,
		$blat,
		$cutadapt,
		$blasr,
		$consed,
		$makeregions,
		$flash,
		$besst,
		$gap2seq,
		$keeptemp,
		$cleantemp);

        #Check for options
	if (!GetOptions(
		'1=s' => \@upstream_mates,
		'2=s' => \@downstream_mates,
		'o=s' => \$output_dir,
		'h|help' => \$help,
		'v|version' => \$version,
		'executables' => \$executables,
		'vector=s' => \$vector,
		'host=s' => \$host,
		'adapter=s' => \@adapter,
		'pacbioroi=s' => \@pacbioroi,
		'pacbiofsr=s' => \@pacbiofsr,
		'finished=s' => \@finished,
		'draft=s' => \@draft,
		'left_end=s' => \$left_end,
		'right_end=s' => \$right_end,
		'peptides=s' => \@peptides,
		'spades=s' => \$spades,
		'samtools=s' => \$samtools,
		'bowtie2build=s' => \$bowtie2build,
		'bowtie2=s' => \$bowtie2,
		'blat=s' => \$blat,
		'cutadapt=s' => \$cutadapt,
		'blasr=s' => \$blasr,
		'consed=s' => \$consed,
		'makeregions=s' => \$makeregions,
		'flash=s' => \$flash,
		'besst=s' => \$besst,
		'gap2seq=s' => \$gap2seq,
		'keeptemp' => \$keeptemp,
		'clean'	=> \$cleantemp)){
		usage;
		exit -1;
	}

	#make sure the executables are present and sane, keep track of how many are missing
	if ($cleantemp){
		$keeptemp = 0;
	}else{
		$keeptemp = 1;
	}
	($spades_exec, $executables) = check_executable($spades, $spades_default, $executables);
	($samtools_exec, $executables) = check_executable($samtools, $samtools_default, $executables);
	($bowtie2build_exec, $executables) = check_executable($bowtie2build, $bowtie2build_default, $executables);
	($bowtie2_exec, $executables) = check_executable($bowtie2, $bowtie2_default, $executables);
	($blat_exec, $executables) = check_executable($blat, $blat_default, $executables);
	($cutadapt_exec, $executables) = check_executable($cutadapt, $cutadapt_default, $executables);
	($blasr_exec, $executables) = check_executable($blasr, $blasr_default, $executables);
	($consed_exec, $executables) = check_executable($consed, $consed_default, $executables);
	($makeregions_exec, $executables) = check_executable($makeregions, $makeregions_default, $executables);
	($flash_exec, $executables) = check_executable($flash, $flash_default, $executables);
	($besst_exec, $executables) = check_executable($besst, $besst_default, $executables);
	($gap2seq_exec, $executables) = check_executable($gap2seq, $gap2seq_default, $executables);

	@upstream_mates = split(/,/,join(',',@upstream_mates));
	@downstream_mates = split(/,/,join(',',@downstream_mates));
	my $mates_warning = "";
	my $library_count = @upstream_mates;
	if (@upstream_mates == @downstream_mates){
		my @ups = @upstream_mates;
		my @downs = @downstream_mates;
		while (@ups || @downs){
			my $up = shift(@ups);
			my $down = shift(@downs);
			if (-e $up && -e $down){

			}elsif (-e $up){
				$mates_warning .= "\t\t$down\n";
			}elsif (-e $down){
				$mates_warning .= "\t$up\t\n";
			}else{
				$mates_warning .= "\t$up\t$down\n";
			}
		}
		if ($mates_warning){
			$mates_warning = "Some read files were not found\n".$mates_warning;
		}
	}else{
		$help = 1;
		$mates_warning .= "Some read files are not paired:\n";
		while (@upstream_mates || @downstream_mates){
			my $up = shift(@upstream_mates) || '';
			my $down = shift(@downstream_mates) || '';
			$mates_warning .= "\t$up\t$down\n";
		}
	}
	if ($mates_warning){
		$help = 1;
		warn $mates_warning;
	}
	my $output_warning = "";
	if (defined($output_dir)){
		if (-e $output_dir){
			if (-d $output_dir){
				if (-w $output_dir){
					$mates_warning .= "$output_dir already exists, files may be overwritten\n";
				}else{
					$help = 1;
					$mates_warning .= "$output_dir is not writable\n";
				}
			}else{
				$help = 1;
				$mates_warning .= "$output_dir is not a directory\n";
			}
		}else{
			$mates_warning .= "creating directory $output_dir\n";
			mkdir abs_path($output_dir);
		}
	}else{
		$help = 1;
		$mates_warning .= "No output directory specified!\n"
	}
	if ($mates_warning){
		warn $mates_warning;
	}

	if ($executables){
		executables();
		exit -1;
	}

	if ($help){
		usage();
		exit -1;
	}

	if ($version){
		version();
		exit -1;
	}

	#ready to go!

	my @uorig = @upstream_mates;
	my @dorig = @downstream_mates;

	push(@adapter, $adapter_default);


	#perform adatper and quality trimming

	if (@adapter){
		@adapter = split(/,/,join(',',@adapter));
		foreach my $adapter (@adapter){
			$cutadapt_params .= " -b $adapter -B $adapter";
		}
		my ($u, $d) = mask_adapters($output_dir, \@upstream_mates, \@downstream_mates);
		@upstream_mates = @{$u};
		@downstream_mates = @{$d};
		push(@temporary, @upstream_mates);
		push(@temporary, @downstream_mates);
	}

	my @utrimmed = @upstream_mates;
	my @dtrimmed = @downstream_mates;

	#screen out reads matching host strain; calculate average library fragment size and standard deviation, average read length

	my ($avg_frag, $std_frag, $avg_read, $cutoff) = ();

	if (defined($host)){
		my ($u, $d, $f, $s, $r) = screen_contamination("host", $output_dir, $host, \@upstream_mates, \@downstream_mates);
		@upstream_mates = @{$u};
		@downstream_mates = @{$d};
		push(@temporary, @upstream_mates);
		push(@temporary, @downstream_mates);
		($avg_frag, $std_frag, $avg_read) = ($f, $s, $r);
		print "Average fragment size: $avg_frag\nStandard Deviation: $std_frag\nAverage Read: $avg_read\n";
	}

	#screen out reads matching vector; caluculate average fold coverage

	my $cov_cutoff;

	if (defined($vector)){
		$cov_cutoff = estimate_coverage_cutoff_from_vector($vector, \@upstream_mates, \@downstream_mates);
		my ($u, $d, $f, $s, $r) = screen_contamination("vector", $output_dir, $vector, \@upstream_mates, \@downstream_mates);
		@upstream_mates = @{$u};
		@downstream_mates = @{$d};
		push(@temporary, @upstream_mates);
		push(@temporary, @downstream_mates);
		unless ($avg_frag){
			($avg_frag, $std_frag, $avg_read) = ($f, $s, $r);
			print "Average fragment size: $avg_frag\nStandard Deviation: $std_frag\nAverage Read: $avg_read\n";
		}
	}
	print "Average fragment size: $avg_frag\nStandard Deviation: $std_frag\nAverage Read: $avg_read\n";

	if ($cov_cutoff){
		print "predicted coverage: $cov_cutoff\n";
	}

	#decide whether to run flash to merge paired reads -- do many pairs overlap?

	if ($avg_frag > 0 && $avg_frag - $std_frag - 2*$avg_read < 0){
		my ($u, $d, $s) = overlap_mates($output_dir, \@upstream_mates, \@downstream_mates, $avg_frag, $std_frag, $avg_read);
		@upstream_mates = @{$u};
		@downstream_mates = @{$d};
		@singles = @{$s};
	}

	#build up arguments for spades

	#start with the reads

	my $spades_arguments = make_illumina_arguments(\@upstream_mates, \@downstream_mates, \@singles);

	#are there any finished sequences we should add?

	my $trusted = "";
	if (@finished){
		@finished = split(/,/,join(',',@finished));
		push(@finished, $left_end, $right_end);
		$trusted = merge_fastas("trusted_contigs", $output_dir, \@finished);
		$spades_arguments .= " --trusted-contigs $trusted";
		push(@temporary, $trusted);
	}

	#are there any less reliable draft sequences?

	my $untrusted = "";
	if (@draft){
		@draft = split(/,/,join(',',@draft));
		$untrusted = merge_fastas("untrusted_contigs", $output_dir, \@draft);
		$spades_arguments .= " --untrusted-contigs $untrusted";
		push(@temporary, $untrusted);
	}

	#we can use peptides for ordering and orienting contigs later

	my $peptide = "";
	if (@peptides){
		@peptides = split(/,/,join(',',@peptides));
		$peptide = merge_fastas("peptides", $output_dir, \@peptides);
		push(@temporary, $peptide);
	}

	#where should the spades output go?

	my $spades_dir = "$output_dir/$date"."_spades";
	push(@temporary, $spades_dir);

	#because we've done the quality trimming, run --only-assembler, error correction won't be helpful
	#running with --careful enables repeat resolution

	my $assembly_command = "$spades_exec $spades_arguments  --only-assembler --careful -o $spades_dir";

	#if we know the expected coverage, that can help get rid of low coverage junk contigs

	if ($cov_cutoff){
		$assembly_command .= " --cov-cutoff $cov_cutoff";
	}

	#let spades do it's thing

	system($assembly_command);

	#check to see if spades was successful; it should make scaffolds.fasta
	#if it fails, try again with only paired reads and one k-mer size
	#trusted-contigs, untrusted-contigs, and single ended reads sometimes cause crashes

  my $scaffolds = "$spades_dir/scaffolds.fasta";

  if (-e $scaffolds){
  	print "Successful assembly\n";
  }else{
    print "Warning: Assembly failed using command:\n\n$assembly_command\n\nFalling back to a simpler command\n";
		my @null_list = ();
    my $simple_arguments = make_illumina_arguments(\@upstream_mates, \@downstream_mates, \@null_list);
    $assembly_command = "$spades_exec $simple_arguments  --only-assembler -k 127 --careful -o $spades_dir";
		if ($cov_cutoff){
			$assembly_command .= " --cov-cutoff $cov_cutoff";
		}
		system($assembly_command);
		if (-e $scaffolds){
			print "Success with command:\n\n$assembly_command\n\nProceeding to next step\n";
		}else{
			die "Failed to generate an assembly with spades!\n";
		}
  }

	#if we have pacbio reads, we'll add them in if they match existing contigs
	#this lets you use pacbio reads from a pooled clones

	#check for "reads of insert" or high-accuracy CCS reads

	my $pacbioroifqgz = "";
	if (@pacbioroi){
		@pacbioroi = split(/,/,join(',',@pacbioroi));
		$pacbioroifqgz = screen_pacbio("scaffold_match_roi", $output_dir, "$spades_dir/scaffolds.fasta", \@pacbioroi);
		my $n = @upstream_mates + 1;
		$spades_arguments .= " --s$n $pacbioroifqgz";
		push(@temporary, $pacbioroifqgz);
	}

	#check for "filtered subreads" or low-accuracy single-pass reads

	my $pacbiofsrfqgz = "";
	if (@pacbiofsr){
		@pacbiofsr = split(/,/,join(',',@pacbiofsr));
		$pacbiofsrfqgz = screen_pacbio("scaffold_match_fsr", $output_dir, "$spades_dir/scaffolds.fasta", \@pacbiofsr);
		$spades_arguments .= " --pacbio $pacbiofsrfqgz";
		push(@temporary, $pacbiofsrfqgz);
	}


	$assembly_command = "$spades_exec $spades_arguments --careful -o $spades_dir";
	if ($cov_cutoff){
		$assembly_command .= " --cov-cutoff $cov_cutoff";
	}

	#run spades again with pacbio data, if you have it

	if (-e $pacbioroifqgz || -e $pacbiofsrfqgz) {
		system($assembly_command);
	}


	#collect up trimmed illumina reads

	my $utrim = join(',',@utrimmed);
	my $dtrim = join(',',@dtrimmed);

	my $besst_bowtie_index = "$output_dir/scaffolding";
	my $besst_bowtie_output = "$output_dir/scaffolding.srt.bam";

	#align reads back to scaffolds, then use that alignment as input to BESST
	#then fill gaps in BESST scaffolds with Gap2Seq

	if (-e $scaffolds){
		print "Scaffolding with illumina reads\n";
		system("$bowtie2build_exec -q $scaffolds $besst_bowtie_index");
		system("$bowtie2_exec -x $besst_bowtie_index -1 $utrim -2 $dtrim | $samtools_exec view -b -S - | $samtools_exec sort - >$besst_bowtie_output");
		system("$samtools_exec index $besst_bowtie_output");
		system("$besst_exec -c $scaffolds -f $besst_bowtie_output -o $output_dir --orientation fr");
		my $besst_scaffolds = "$output_dir/BESST_output/pass1/Scaffolds_pass1.fa";
		my $filled_gaps = "$output_dir/filled_gaps.fa";
		if (-e $besst_scaffolds){
			print "Scaffolding successful\n";
			print "Filling gaps with illumina reads\n";
			$scaffolds = $besst_scaffolds;
			eval {
				local $SIG{ALRM} = sub { die "alarm\n" };
				alarm $gap2seq_timeout;
				system("$gap2seq_exec -scaffolds $besst_scaffolds -filled $filled_gaps -reads $utrim,$dtrim");
				alarm 0;
		  };
			if ($@){
				if ($@ eq"alarm\n"){
					print "$gap2seq_exec took more than $gap2seq_timeout seconds; that's too long\n"
				}else{
					print "Unexpected error from $gap2seq_exec: $@\n";
				}
			}
			if (-e $filled_gaps){
				print "Gap Filling Successful\n";
				$scaffolds = $filled_gaps;
			}else{
				print "Gap Filling Failed\n";
			}
		}else{
			print "Scaffolding Failed\n";
		}
	}

	#collect up filtered illumina reads
	my $ufilt = join(',',@upstream_mates);
	my $dfilt = join(',',@downstream_mates);
	my $sfilt = join(',',@singles);
	my $final = "$output_dir/final.fasta";

	#load up the output of BESST and Gap2Seq and write in the output folder

	my $contigs = load_fasta_contigs($scaffolds);
	write_fasta_sequences($contigs, $final, $avg_read);
	my $input_contigs = load_fasta_contigs($final);

	#if there are end sequences, try to order and orient remaining contigs
	#overwrite the final.fasta file with ordered and oriented contigs

	if (-e $left_end && -e $right_end){
		my $autofinished_contigs = autofinish($input_contigs, $final, $left_end, $right_end, $peptide);
		my $auto_out = "$output_dir/Autofinished.fasta";
		fasta_output_autofinished($autofinished_contigs, $auto_out);
		fasta_output_autofinished($autofinished_contigs, $final);
	}

	#align reads to final.fasta to produce output for consed

	my $scaffoldssrt = "$output_dir/final.srt.bam";
	my $consed_dir = "$output_dir/consed";
	system("$bowtie2build_exec -q $final $output_dir/final");
	system("$bowtie2_exec -I 0 -X 2501 --rdg 502,502 --rfg 502,502 -x $output_dir/final -1 $ufilt -2 $dfilt | $samtools_exec view -b -S - | $samtools_exec sort -o $scaffoldssrt -");
	system("$samtools_exec index $scaffoldssrt");
	system("$makeregions_exec $final");

	#consed won't write to an exisitng directory, so wipe those out if they exist.

	if (-e $consed_dir) {
		rmtree([ "$consed_dir" ]) || print "$! : for $consed_dir\n";
	}
	system("$consed_exec -bam2ace -bamFile $scaffoldssrt -regionsFile $output_dir/finalRegions.txt -dir $consed_dir");

	my @bowties = glob("$output_dir"."*"."bt2");
	push(@temporary,@bowties);

	#unless we're keeping temporary files, remove them to save space

	unless ($keeptemp){
		rmtree([ @temporary ]);
	}

	#we're all done!

}


#use bowtie and samtools to get the average depth across the vector sequence

sub estimate_coverage_cutoff_from_vector ($$$) {
	my ($vector, $ulist, $dlist) = @_;
	my $ups = join(',',@{$ulist});
	my $downs = join(',',@{$dlist});
	my (%tot,%max,%ave) = ();
	open (DEPTH, "$bowtie2_exec $bowtie2_screen_params -x $vector -1 $ups -2 $downs | samtools view -S -b - | samtools sort -o - | samtools depth -aa - |");
	while (<DEPTH>){
		chomp;
		my ($s, $p, $d) = split(/\t/, $_);
		unless ($phix_def =~ m/$s/){
			$tot{$s} += $d;
			$max{$s} = $p;
			$ave{$s} = $tot{$s} / $p;
		}
	}
	close DEPTH;
	my @sorted = sort {$ave{$b} <=> $ave{$a}} (keys(%ave));
	return int($ave{$sorted[0]});
}


#uses bioperl to write sequences in fasta format

sub write_fasta_sequences ($$$){
	my ($contigs, $file, $avg_read_len) = @_;
	my $out = Bio::SeqIO->new(-file => ">$file", '-format' => 'Fasta');
	my $n = 0;
	foreach my $id (keys(%$contigs)){
		#get rid of 128bp junk contigs
		unless ($contigs->{$id}->length() < $avg_read_len){
			$out->write_seq($contigs->{$id});
			$n++
		}
	}
	return $n;
}


#loads a fasta file

sub load_fasta_contigs ($) {
	my $file = shift(@_);
	my %contigs;
	my $in = Bio::SeqIO->new(       -format	=> 'fasta',
					-file	=> $file);
	while (my $c = $in->next_seq()){
		 $contigs{$c->id} = $c;
	}
	return \%contigs;
}

#concatenates multiple fasta files into one for spades

sub merge_fastas ($$$){
	my ($text, $dir, $fr) = @_;
	my $file = "$dir/$text.fa";
	my $out = Bio::SeqIO->new(-file => ">$file", '-format' => 'Fasta');
	my @fastas = @{$fr};
	foreach my $f (@fastas){
			my $in = Bio::SeqIO->new(	-format	=> 'fasta',
																-file	=> $f);
			while (my $c = $in->next_seq()){
				$out->write_seq($c);
			}
	}
	return $file;
}

#uses blasr and samtools to extract pacbio reads matching contigs

sub screen_pacbio ($$$$){
	my ($text, $dir, $scaffolds, $pbr) = @_;
	my @pacbio = @{$pbr};
	my $output = "$dir/$text.fq.gz";
	my $z = new IO::Compress::Gzip $output, Level => 9 or die "couldn't write outfile: $output\nIO::Compress::Gzip failed: $GzipError\n";
	foreach my $lib (1 .. @pacbio){
		open (SAM, "$blasr_exec $pacbio[$lib-1] $scaffolds -bestn 1 -sam | $samtools_exec view -S -F 4 -F 256 - |") || die "failed aligning $pacbio[$lib-1] to $scaffolds\n";
		while (<SAM>){
			my @line = split(/\s+/, $_);
			print $z "@"."$line[0]\n$line[9]\n+\n$line[10]\n";
		}
		close SAM;

	}
	close $z;
	return $output;
}

#turn arrays of read files into string of arguments for spades

sub make_illumina_arguments ($$$){
	my ($upr, $downr, $singler) = @_;
	my @ups = @{$upr};
	my @downs = @{$downr};
	my @singles = @{$singler};
	my $count = @ups;
	my $string = "";
	foreach my $lib (1 .. $count){
		my ($u, $d) = ($ups[$lib-1], $downs[$lib-1]);
		if ($string =~ m/\S/){
			$string .= " ";
		}
		$string .= "--pe".$lib."-1 $u --pe".$lib."-2 $d";
	}
	if (@singles){
		foreach my $lib (1 .. $count){
			my $s = $singles[$lib-1];
			my $libn = $count + $lib;
			if ($string =~ m/\S/){
				$string .= " ";
			}
			$string .= "--s".$libn." $s";
		}
	}
	return $string;
}


#runs cutadapt to remove adapter sequences from reads

sub mask_adapters ($$$) {
	my ($dir, $upr, $downr) = @_;
	my @ups = @{$upr};
	my @downs = @{$downr};
	my $count = @ups;
	my (@firsts, @seconds);

	foreach my $lib (1 .. $count){
		my ($u, $d) = ($ups[$lib-1], $downs[$lib-1]);
		my ($utmp, $dtmp) = ("$dir/$lib.cutadapt.1.tmp.fq", "$dir/$lib.cutadapt.2.tmp.fq");
		my ($uout, $dout) = ("$dir/$lib.cutadapt.1.out.fq", "$dir/$lib.cutadapt.2.out.fq");
		my ($ufqgz, $dfqgz) = ("$dir/$lib.cutadapt.1.fq.gz", "$dir/$lib.cutadapt.2.fq.gz");

		print "$ufqgz, $dfqgz\n";
		system("$cutadapt_exec $cutadapt_params -o $ufqgz -p $dfqgz $u $d");
		push(@firsts, $ufqgz);
		push(@seconds, $dfqgz);
	}
	return (\@firsts, \@seconds);
}

#remove contaminating reads from host or vector
#uses bowtie to align, but perl to parse because
#I want to keep the reads that don't align
#but count the ones that do


sub screen_contamination ($$$$$) {
	my ($text, $dir, $index, $upr, $downr) = @_;
	my @ups = @{$upr};
	my @downs = @{$downr};
	my $count = @ups;
	my (@firsts, @seconds);

	my @line;
	my @frags;

	my $proper = 0x0002;
	my $unmaps = 0x0004;
	my $unmapp = 0x0008;
	my $firstp = 0x0040;
	my $secndp = 0x0080;


	my @frags;
	my @reads;

	foreach my $lib (1 .. $count){
		my ($u, $d) = ($ups[$lib-1], $downs[$lib-1]);
		my ($ufqgz, $dfqgz) = ("$dir/$lib.$text.1.fq.gz", "$dir/$lib.$text.2.fq.gz");
		print "$ufqgz, $dfqgz\n";

		my $uz = new IO::Compress::Gzip $ufqgz, Level => 9 or die "can't write $text screened reads to $ufqgz\nIO::Compress::Gzip failed: $GzipError\n";
		my $dz = new IO::Compress::Gzip $dfqgz, Level => 9 or die "can't write $text screened reads to $dfqgz\nIO::Compress::Gzip failed: $GzipError\n";

		open (SAM, "$bowtie2_exec $bowtie2_screen_params -x $index -1 $u -2 $d |") || die "can't run bowtie to screen $text\n";
		while (<SAM>){
			if (m/^\@/){
			}else{
				@line = split(/\t/, $_);
				if ($line[1] & $unmaps && $line[1] & $unmapp){
					push(@reads, length($line[9]));
					if ($line[1] & $firstp){
						print $uz "@"."$line[0]\n$line[9]\n+$line[0]\n$line[10]\n";
					}elsif ($line[1] & $secndp){
						print $dz "@"."$line[0]\n$line[9]\n+$line[0]\n$line[10]\n";
					}
				}elsif ($line[8]){
					push(@frags, abs($line[8]));
				}
			}
		}
		close SAM;
		close $uz;
		close $dz;
		push(@firsts, $ufqgz);
		push(@seconds, $dfqgz);
	}
  my $fstat = Statistics::Descriptive::Full->new();
	$fstat->add_data(@frags);
	my $rstat = Statistics::Descriptive::Full->new();
	$rstat->add_data(@reads);

  my $avg_frag = int(1 + $fstat->mean());
  my $std_frag = int(1 + $fstat->standard_deviation());
  my $avg_read = int(1 + $rstat->mean());

	return (\@firsts, \@seconds, $avg_frag, $std_frag, $avg_read);
}

#runs flash to overlap paired reads
#put a function prototype up top!

sub overlap_mates ($$$$$$) {
	my ($dir, $upr, $downr, $f, $s, $r) = @_;
	my @ups = @{$upr};
	my @downs = @{$downr};
	my $count = @ups;
	my (@firsts, @seconds, @singles);
	foreach my $lib (1 .. $count){
		my ($u, $d) = ($ups[$lib-1], $downs[$lib-1]);
		my ($unc, $dnc, $ext) = ("$dir/$lib.notCombined_1.fastq", "$dir/$lib.notCombined_2.fastq", "$dir/$lib.extendedFrags.fastq");
		print "flash command:\n$flash_exec $u $d -o $lib -d $dir -f $f -s $s -r $r\n";
		system("$flash_exec $u $d -o $lib -d $dir -f $f -s $s -r $r");
		push(@firsts, $unc);
		push(@seconds, $dnc);
		push(@singles, $ext);
	}
	return (\@firsts, \@seconds, \@singles);
}

#checks to make sure all executables are reachable and executable!

sub check_executable ($$$) {
	my ($option, $default, $execs) = @_;
	my $value = '';

	if (defined($option) && -x $option){
		$value = $option;
	}else{
		if (defined($option)){
			warn "$option is not executable\n";
			$execs++;
		}
		if (defined($default) && -x $default){
			$value = $default;
		}else{
			if (defined($default)){
				warn "$default is not executable\n";
				$execs++;
			}
		}
	}
	return ($value, $execs);
}

#prints version information

sub version() {
	print qq/SHIMS Pipeline
(Version: $VERSION)
/;
	return 1;
}

#prints paths to executables, and shows how to set environment variables

sub executables() {
	print qq/
Changing Executables:
	--spades        <path to SPAdes: $spades_exec>
	--samtools      <path to samtools: $samtools_exec>
	--bowtie2build  <path to bowtie2-build: $bowtie2build_exec>
	--bowtie2       <path to bowtie2: $bowtie2_exec>
	--blat          <path to blat: $blat_exec>
	--cutadapt      <path to cutadapt: $cutadapt_exec>
	--blasr         <path to blasr: $blasr_exec>
	--consed        <path to consed: $consed_exec>
	--makeregions   <path to consed's makeRegionsFile.perl: $makeregions_exec>
	--flash         <path to flash: $flash_exec>
	--besst         <path to besst: $besst_exec>
	--gap2seq       <path to gap2seq: $gap2seq_exec>

To make the current executables the default:

export SHIMS_SPADES_EXEC=$spades_exec
export SHIMS_SAMTOOLS_EXEC=$samtools_exec
export SHIMS_BOWTIE2BUILD_EXEC=$bowtie2build_exec
export SHIMS_BOWTIE2_EXEC=$bowtie2_exec
export SHIMS_BLAT_EXEC=$blat_exec
export SHIMS_CUTADAPT_EXEC=$cutadapt_exec
export SHIMS_BLASR_EXEC=$blasr_exec
export SHIMS_CONSED_EXEC=$consed_exec
export SHIMS_MAKEREGIONS_EXEC=$makeregions_exec
export SHIMS_FLASH_EXEC=$flash_exec
export SHIMS_BESST_EXEC=$besst_exec
export SHIMS_GAP2SEQ_EXEC=$gap2seq_exec
	/;

	return 1;
}

#prints usage information.

sub usage() {
	version();
	print qq/

USAGE: $0 -1 <upstream mates> -2 <downstream mates> -o <output directory> [optional arguments]

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
		--spades        <path to SPAdes: $spades_exec>
		--samtools      <path to samtools: $samtools_exec>
		--bowtie2build  <path to bowtie2-build: $bowtie2build_exec>
		--bowtie2       <path to bowtie2: $bowtie2_exec>
		--blat          <path to blat: $blat_exec>
		--cutadapt      <path to cutadapt: $cutadapt_exec>
		--blasr         <path to blasr: $blasr_exec>
		--consed        <path to consed: $consed_exec>
		--makeregions   <path to consed's makeRegionsFile.perl: $makeregions_exec>
		--flash         <path to flash: $flash_exec>
		--besst         <path to besst: $besst_exec>
		--gap2seq       <path to gap2seq: $gap2seq_exec>
/;

	return 1;

}

#autofinisher -- complex code here below, need to double check it.

sub autofinish ($$$$$) {
	my ($c, $contigs, $le, $re, $peptides) = @_;
	my ($bestleft, $bestright);
	my $r;
	print "#autofinisher: Left: $le Right: $re\n";
	print "#autofinisher: running blat\n";
	open (LE, "$blat_exec $contigs $le stdout | ") || die "can't read blat output for left end: $le";
	while (<LE>){
		chomp;
		my @psl = split(/\s+/, $_);
		my $score = ($psl[0] + ($psl[2] >> 1) - $psl[1] - $psl[4] - $psl[6]);
		if ($score > 50 && $score > $bestleft->{'score'}){
			print "#autofinisher: L: $_\n";
			$bestleft = {
				'score'		=> $score,
				'strand'	=> $psl[8],
				'tname'		=> $psl[13],
				'tstart'	=> $psl[15],
				'tend'		=> $psl[16]
			};
		}
	}
	close LE;
	open (RE, "$blat_exec $contigs $re stdout | ") || die "can't read blat output for right end: $re";
	while (<RE>){
		chomp;
		my @psl = split(/\s+/, $_);
		my $score = ($psl[0] + ($psl[2] >> 1) - $psl[1] - $psl[4] - $psl[6]);
		if ($score > 50 && $score > $bestright->{'score'}){
			print "#autofinisher: R: $_\n";
			$bestright = {
				'score'		=> $score,
				'strand'	=> $psl[8],
				'tname'		=> $psl[13],
				'tstart'	=> $psl[15],
				'tend'		=> $psl[16]
			};
		}
	}
	close RE;

	print "#autofinisher: after blat L: $bestleft->{'tname'} R: $bestright->{'tname'}\n";

	### Decide if this is simple or not

	if ($bestleft->{'tname'} && $bestright->{'tname'} && $bestleft->{'tname'} eq $bestright->{'tname'}){
		### this is the simple case -- one contig matching both ends
		print "#autofinisher: simple case of single contig\n";
		my $id = $bestright->{'tname'};
		### okay -- is it forward or reverse?
		if ($bestleft->{'strand'} eq "+" && $bestright->{'strand'} eq "-" && $bestleft->{'tstart'} < $bestright->{'tend'}){
			### forward case
			print "#autofinisher: contig is in forward orientation\n";
			$r->{$id} = make_result_contig($c->{$id}, 1, 1, 1, 1);
			return $r;
		}elsif ($bestleft->{'strand'} eq "-" && $bestright->{'strand'} eq "+" && $bestright->{'tstart'} < $bestleft->{'tend'}){
			### reverse case
			print "#autofinisher: contig is in reverse orientation\n";
			$r->{$id} = make_result_contig($c->{$id}, 1, 1, 1, -1);
			return $r;
		}else{
			### bad case -- lets toss out the lower scoring end
			print "#autofinisher: contig ends match in an unusual way\n";
			if ($bestleft->{'score'} > $bestright->{'score'}){
				$bestright = {};
			}else{
				$bestleft = {};
			}
		}
	}

	### Still here? Oh, then it's complicated.

	print "#autofinisher: complex case\n";

	### First lets dump all the contigs into the result hash, but mark them as unordered and unoriented

	foreach my $id (keys(%{$c})){
		$r->{$id} = make_result_contig($c->{$id}, 0, 0, 0, 1);
	}

	### Now lets see if the end matches let us find the first and last contig

	my ($first_id, $last_id);
	if ($bestleft->{'tname'}){
		$first_id = $bestleft->{'tname'};
		$r->{$first_id}->{'number'} = 1;
		$r->{$first_id}->{'ordered'} = 1;
		$r->{$first_id}->{'oriented'} = 1;
		if ($bestleft->{'strand'} eq '+'){
		}else{
			$r->{$first_id}->{'strand'} = -1
		}
		print "#autofinisher: $first_id contig is on the left\n";
	}

	### Remember that a plus stranded hit with the right end means that the last contig should be reversed

	if ($bestright->{'tname'}){
		$last_id = $bestright->{'tname'};
		$r->{$last_id}->{'number'} = -1;
		$r->{$last_id}->{'ordered'} = 1;
		$r->{$last_id}->{'oriented'} = 1;
		if ($bestright->{'strand'} eq '+'){
			$r->{$last_id}->{'strand'} = -1
		}else{
		}
		print "#autofinisher: $last_id contig is on the right\n";
	}


	### We can align to peptides to try to orient the other contigs

	my ($peptide_matches, $matrix, $peptide_info);;

	if ($peptides){
		print "#autofinisher: conducting alignments to peptides in $peptides\n";

		open (PEP, "$blat_exec $contigs $peptides -t=dnax -q=prot stdout |") || die "can't blat peptides in $peptides versus contigs in $contigs\n";
		while (<PEP>){
			chomp;
			my @psl = split(/\s+/, $_);
			my $score = ($psl[0] + ($psl[2] >> 1) - $psl[1] - $psl[4] - $psl[6]);
                        #median exon is ~90bp
			if ($score > 50){
				print "#autofinisher:\t$_\n";
				foreach my $aa ($psl[11] .. $psl[12]){
					if ($peptide_matches->{$psl[9]}->{$aa}->{'score'} < $score){
						$peptide_matches->{$psl[9]}->{$aa}->{'score'} =  $score;
						$peptide_matches->{$psl[9]}->{$aa}->{'best'} = $psl[13];
						$peptide_matches->{$psl[9]}->{$aa}->{'strand'} = $psl[8];
					}
				}
			}
		}
		close PEP;
		foreach my $pep (keys(%{$peptide_matches})){
			foreach my $aa (sort {$a <=> $b} (keys(%{$peptide_matches->{$pep}}))){
				my $contig = $peptide_matches->{$pep}->{$aa}->{'best'};
				print "#autofinisher: adding information about $contig and $pep\n";
				$matrix->{$contig}->{$pep}->{'aa'} = $aa;
				$matrix->{$contig}->{$pep}->{'score'} = $peptide_matches->{$pep}->{$aa}->{'score'};
				if ($peptide_matches->{$pep}->{$aa}->{'strand'} eq '++'){
					$matrix->{$contig}->{$pep}->{'strand'} = 1;
				}else{
					$matrix->{$contig}->{$pep}->{'strand'} = -1;
				}
			}
		}
	}

	### We'll be done when we assign a number to every contig
	### When we check that, it's a good time to give a status update

	my $innerleft_id = $first_id;
	my $innerright_id = $last_id;
	my $middleright_id;
	my $middleleft_id;

	my $m;
	my $middle_count;

	while (autofinish_status($r)){
		if ($innerleft_id) {
			### start on the left side, going right
			print "#autofinisher: trying to extend inner left contig $innerleft_id to the right\n";
			my $nextleft_id;
			my @peps = sort {$matrix->{$innerleft_id}->{$b}->{'score'} <=> $matrix->{$innerleft_id}->{$a}->{'score'}} keys(%{$matrix->{$innerleft_id}});
			foreach my $pep (@peps){
				print "#autofinisher: checking peptide: $pep\n";
				if ($r->{$innerleft_id}->{'strand'} == 1){
					### case 1: inner left contig is forward orientation
					if ($matrix->{$innerleft_id}->{$pep}->{'strand'} == 1){
						### case 1a: inner left contig is forward orientation; peptide is in forward orientation
						foreach my $contig (sort {$matrix->{$a}->{$pep}->{'aa'} <=> $matrix->{$b}->{$pep}->{'aa'}} (keys(%{$r}))){
							if ($matrix->{$contig}->{$pep}->{'aa'} > $matrix->{$innerleft_id}->{$pep}->{'aa'}){
								### look for the next contig on the peptide
								if ($r->{$contig}->{'ordered'} == 1){
								}else{
									print "#autofinisher: peptide $pep bridges to contig $contig\n";
									$r->{$contig}->{'number'} = $r->{$innerleft_id}->{'number'} + 1;
									$r->{$contig}->{'ordered'} = 1;
									if ($matrix->{$contig}->{$pep}->{'strand'} == 1){
										$r->{$contig}->{'oriented'} = 1;
										$r->{$contig}->{'strand'} = 1;
									}else{
										$r->{$contig}->{'oriented'} = 1;
										$r->{$contig}->{'strand'} = -1;
									}
									$nextleft_id = $contig;
								}
								last;
							}
						}
					}else{
						### case 1b: inner left contig is forward orientation; peptide is in reverse orientation
						foreach my $contig (sort {$matrix->{$b}->{$pep}->{'aa'} <=> $matrix->{$a}->{$pep}->{'aa'}} (keys(%{$r}))){
							if ($matrix->{$contig}->{$pep}->{'aa'} > 0 && $matrix->{$contig}->{$pep}->{'aa'} < $matrix->{$innerleft_id}->{$pep}->{'aa'}){
								### look for the previous contig on the peptide
								if ($r->{$contig}->{'ordered'} == 1){
								}else{
									print "#autofinisher: peptide $pep bridges to contig $contig\n";
									$r->{$contig}->{'number'} = $r->{$innerleft_id}->{'number'} + 1;
									$r->{$contig}->{'ordered'} = 1;
									if ($matrix->{$contig}->{$pep}->{'strand'} == 1){
										$r->{$contig}->{'oriented'} = 1;
										$r->{$contig}->{'strand'} = -1;
									}else{
										$r->{$contig}->{'oriented'} = 1;
										$r->{$contig}->{'strand'} = 1;
									}
									$nextleft_id = $contig;
								}
								last;
							}
						}
					}
				}else{
					### case 2: inner left contig is reverse orientation
					if ($matrix->{$innerleft_id}->{$pep}->{'strand'} == 1){
						### case 2a: inner left contig is reverse orientation; peptide is in forward orientation
						foreach my $contig (sort {$matrix->{$b}->{$pep}->{'aa'} <=> $matrix->{$a}->{$pep}->{'aa'}} (keys(%{$r}))){
            	if ($matrix->{$contig}->{$pep}->{'aa'} > 0 && $matrix->{$contig}->{$pep}->{'aa'} < $matrix->{$innerleft_id}->{$pep}->{'aa'}){
								### look for the previous contig on the peptide
              	if ($r->{$contig}->{'ordered'} == 1){
                }else{
									print "#autofinisher: peptide $pep bridges to contig $contig\n";
	                $r->{$contig}->{'number'} = $r->{$innerleft_id}->{'number'} + 1;
                  $r->{$contig}->{'ordered'} = 1;
                  if ($matrix->{$contig}->{$pep}->{'strand'} == 1){
                  	$r->{$contig}->{'oriented'} = 1;
                    $r->{$contig}->{'strand'} = -1;
                  }else{
                    $r->{$contig}->{'oriented'} = 1;
                    $r->{$contig}->{'strand'} = 1;
                  }
                  $nextleft_id = $contig;
                }
                last;
              }
            }
					}else{
						### case 2b: inner left contig is reverse orientation; peptide is in reverse orientation
						foreach my $contig (sort {$matrix->{$a}->{$pep}->{'aa'} <=> $matrix->{$b}->{$pep}->{'aa'}} (keys(%{$r}))){
            	if ($matrix->{$contig}->{$pep}->{'aa'} > $matrix->{$innerleft_id}->{$pep}->{'aa'}){
								### look for the next contig on the peptide
                if ($r->{$contig}->{'ordered'} == 1){
                }else{
									print "#autofinisher: peptide $pep bridges to contig $contig\n";
                  $r->{$contig}->{'number'} = $r->{$innerleft_id}->{'number'} + 1;
                  $r->{$contig}->{'ordered'} = 1;
                  if ($matrix->{$contig}->{$pep}->{'strand'} == 1){
                  	$r->{$contig}->{'oriented'} = 1;
                    $r->{$contig}->{'strand'} = 1;
                  }else{
                    $r->{$contig}->{'oriented'} = 1;
                    $r->{$contig}->{'strand'} = -1;
                  }
                  $nextleft_id = $contig;
                }
              	last;
              }
            }
					}
				}
				if ($nextleft_id){
					print "#autofinisher: contig $nextleft_id chosen\n";
					last;
				}
			}
			print "#autofinisher: moving on from contig $innerleft_id to $nextleft_id\n";
			$innerleft_id = $nextleft_id;

		}elsif ($innerright_id){

			### start on the right side, going left
			print "#autofinisher: trying to extend inner left contig $innerright_id to the left\n";
			my $nextright_id;
			my @peps = sort {$matrix->{$innerright_id}->{$b}->{'score'} <=> $matrix->{$innerright_id}->{$a}->{'score'}} keys(%{$matrix->{$innerright_id}});
			foreach my $pep (@peps){
				print "#autofinisher: checking peptide: $pep\n";
				if ($r->{$innerright_id}->{'strand'} == 1){
					### case 1: inner right contig is forward orientation
					if ($matrix->{$innerright_id}->{$pep}->{'strand'} == 1){
						### case 1a: inner right contig is forward orientation; peptide is in forward orientation
						foreach my $contig (sort {$matrix->{$b}->{$pep}->{'aa'} <=> $matrix->{$a}->{$pep}->{'aa'}} (keys(%{$r}))){
							if ($matrix->{$contig}->{$pep}->{'aa'} > 0 && $matrix->{$contig}->{$pep}->{'aa'} < $matrix->{$innerright_id}->{$pep}->{'aa'}){
								### look for the previous contig on the peptide
								if ($r->{$contig}->{'ordered'} == 1){
								}else{
									print "#autofinisher: peptide $pep bridges to contig $contig\n";
									$r->{$contig}->{'number'} = $r->{$innerright_id}->{'number'} - 1;
									$r->{$contig}->{'ordered'} = 1;
									if ($matrix->{$contig}->{$pep}->{'strand'} == 1){
										$r->{$contig}->{'oriented'} = 1;
										$r->{$contig}->{'strand'} = 1;
									}else{
										$r->{$contig}->{'oriented'} = 1;
										$r->{$contig}->{'strand'} = -1;
									}
									$nextright_id = $contig;
								}
								last;
							}
						}
					}else{
						### case 1b: inner right contig is forward orientation; peptide is in reverse orientation
						foreach my $contig (sort {$matrix->{$a}->{$pep}->{'aa'} <=> $matrix->{$b}->{$pep}->{'aa'}} (keys(%{$r}))){
							if ($matrix->{$contig}->{$pep}->{'aa'} > $matrix->{$innerright_id}->{$pep}->{'aa'}){
								### look for the next contig on the peptide
								if ($r->{$contig}->{'ordered'} == 1){
								}else{
									print "#autofinisher: peptide $pep bridges to contig $contig\n";
									$r->{$contig}->{'number'} = $r->{$innerright_id}->{'number'} - 1;
									$r->{$contig}->{'ordered'} = 1;
									if ($matrix->{$contig}->{$pep}->{'strand'} == 1){
										$r->{$contig}->{'oriented'} = 1;
										$r->{$contig}->{'strand'} = -1;
									}else{
										$r->{$contig}->{'oriented'} = 1;
										$r->{$contig}->{'strand'} = 1;
									}
									$nextright_id = $contig;
								}
								last;
							}
						}
					}

				}else{
					### case 2: inner right contig is reverse orientation
					if ($matrix->{$innerright_id}->{$pep}->{'strand'} == 1){
						### case 2a: inner right contig is reverse orientation; peptide is in forward orientation
						foreach my $contig (sort {$matrix->{$a}->{$pep}->{'aa'} <=> $matrix->{$b}->{$pep}->{'aa'}} (keys(%{$r}))){
                if ($matrix->{$contig}->{$pep}->{'aa'} > $matrix->{$innerright_id}->{$pep}->{'aa'}){
								### look for the next contig on the peptide
                	if ($r->{$contig}->{'ordered'} == 1){
	                }else{
										print "#autofinisher: peptide $pep bridges to contig $contig\n";
                    $r->{$contig}->{'number'} = $r->{$innerright_id}->{'number'} - 1;
                    $r->{$contig}->{'ordered'} = 1;
                    if ($matrix->{$contig}->{$pep}->{'strand'} == 1){
                    	$r->{$contig}->{'oriented'} = 1;
                      $r->{$contig}->{'strand'} = -1;
                    }else{
                      $r->{$contig}->{'oriented'} = 1;
                      $r->{$contig}->{'strand'} = 1;
                    }
                    $nextright_id = $contig;
                  }
                last;
              }
            }
					}else{
						### case 2b: inner right contig is reverse orientation; peptide is in reverse orientation
						foreach my $contig (sort {$matrix->{$b}->{$pep}->{'aa'} <=> $matrix->{$a}->{$pep}->{'aa'}} (keys(%{$r}))){
                if ($matrix->{$contig}->{$pep}->{'aa'} > 0 && $matrix->{$contig}->{$pep}->{'aa'} < $matrix->{$innerright_id}->{$pep}->{'aa'}){
								### look for the previous contig on the peptide
                if ($r->{$contig}->{'ordered'} == 1){
                }else{
									print "#autofinisher: peptide $pep bridges to contig $contig\n";
                  $r->{$contig}->{'number'} = $r->{$innerright_id}->{'number'} - 1;
                  $r->{$contig}->{'ordered'} = 1;
	                if ($matrix->{$contig}->{$pep}->{'strand'} == 1){
                  	$r->{$contig}->{'oriented'} = 1;
                    $r->{$contig}->{'strand'} = 1;
                  }else{
                    $r->{$contig}->{'oriented'} = 1;
                    $r->{$contig}->{'strand'} = -1;
                  }
                  $nextright_id = $contig;
                }
                last;
              }
            }
					}
				}
				if ($nextright_id){
					print "#autofinisher: contig $nextright_id chosen\n";
					last;
				}
			}
			print "#autofinisher: moving on from contig $innerright_id to $nextright_id\n";
			$innerright_id = $nextright_id;

		}elsif ($middleright_id){

			### start from middle on the right side, going left
			print "#autofinisher: trying to extend middle right contig $middleright_id to the left\n";
			my $nextright_id;
			my @peps = sort {$matrix->{$middleright_id}->{$b}->{'score'} <=> $matrix->{$middleright_id}->{$a}->{'score'}} keys(%{$matrix->{$middleright_id}});
			foreach my $pep (@peps){
				print "#autofinisher: checking peptide: $pep\n";
				if ($r->{$middleright_id}->{'strand'} == 1){
					### case 1: inner right contig is forward orientation
					if ($matrix->{$middleright_id}->{$pep}->{'strand'} == 1){
						### case 1a: inner right contig is forward orientation; peptide is in forward orientation
						foreach my $contig (sort {$matrix->{$b}->{$pep}->{'aa'} <=> $matrix->{$a}->{$pep}->{'aa'}} (keys(%{$r}))){
							if ($matrix->{$contig}->{$pep}->{'aa'} > 0 && $matrix->{$contig}->{$pep}->{'aa'} < $matrix->{$middleright_id}->{$pep}->{'aa'}){
								### look for the previous contig on the peptide
								if ($r->{$contig}->{'ordered'} == 1){
								}else{
									print "#autofinisher: peptide $pep bridges to contig $contig\n";
									$m->{$contig}->{'number'} = $m->{$middleright_id}->{'number'} - 1;
									$r->{$contig}->{'ordered'} = 1;
									if ($matrix->{$contig}->{$pep}->{'strand'} == 1){
										$r->{$contig}->{'oriented'} = 1;
										$r->{$contig}->{'strand'} = 1;
									}else{
										$r->{$contig}->{'oriented'} = 1;
										$r->{$contig}->{'strand'} = -1;
									}
									$nextright_id = $contig;
								}
								last;
							}
						}
					}else{
						### case 1b: inner right contig is forward orientation; peptide is in reverse orientation
						foreach my $contig (sort {$matrix->{$a}->{$pep}->{'aa'} <=> $matrix->{$b}->{$pep}->{'aa'}} (keys(%{$r}))){
							if ($matrix->{$contig}->{$pep}->{'aa'} > $matrix->{$middleright_id}->{$pep}->{'aa'}){
								### look for the next contig on the peptide
								if ($r->{$contig}->{'ordered'} == 1){
								}else{
									print "#autofinisher: peptide $pep bridges to contig $contig\n";
									$m->{$contig}->{'number'} = $m->{$middleright_id}->{'number'} - 1;
									$r->{$contig}->{'ordered'} = 1;
									if ($matrix->{$contig}->{$pep}->{'strand'} == 1){
										$r->{$contig}->{'oriented'} = 1;
										$r->{$contig}->{'strand'} = -1;
									}else{
										$r->{$contig}->{'oriented'} = 1;
										$r->{$contig}->{'strand'} = 1;
									}
									$nextright_id = $contig;
								}
								last;
							}
						}
					}

				}else{
					### case 2: inner right contig is reverse orientation
					if ($matrix->{$middleright_id}->{$pep}->{'strand'} == 1){
						### case 2a: inner right contig is reverse orientation; peptide is in forward orientation
						foreach my $contig (sort {$matrix->{$a}->{$pep}->{'aa'} <=> $matrix->{$b}->{$pep}->{'aa'}} (keys(%{$r}))){
              if ($matrix->{$contig}->{$pep}->{'aa'} > $matrix->{$middleright_id}->{$pep}->{'aa'}){
								### look for the next contig on the peptide
	              if ($r->{$contig}->{'ordered'} == 1){
                }else{
									print "#autofinisher: peptide $pep bridges to contig $contig\n";
                  $m->{$contig}->{'number'} = $m->{$middleright_id}->{'number'} - 1;
                  $r->{$contig}->{'ordered'} = 1;
                  if ($matrix->{$contig}->{$pep}->{'strand'} == 1){
                    $r->{$contig}->{'oriented'} = 1;
                    $r->{$contig}->{'strand'} = -1;
                  }else{
                    $r->{$contig}->{'oriented'} = 1;
                    $r->{$contig}->{'strand'} = 1;
                  }
                  $nextright_id = $contig;
                }
                last;
              }
            }
					}else{
						### case 2b: inner right contig is reverse orientation; peptide is in reverse orientation
						foreach my $contig (sort {$matrix->{$b}->{$pep}->{'aa'} <=> $matrix->{$a}->{$pep}->{'aa'}} (keys(%{$r}))){
              if ($matrix->{$contig}->{$pep}->{'aa'} > 0 && $matrix->{$contig}->{$pep}->{'aa'} < $matrix->{$middleright_id}->{$pep}->{'aa'}){
								### look for the previous contig on the peptide
                if ($r->{$contig}->{'ordered'} == 1){
                }else{
									print "#autofinisher: peptide $pep bridges to contig $contig\n";
                  $m->{$contig}->{'number'} = $m->{$middleright_id}->{'number'} - 1;
                  $r->{$contig}->{'ordered'} = 1;
                  if ($matrix->{$contig}->{$pep}->{'strand'} == 1){
                    $r->{$contig}->{'oriented'} = 1;
                    $r->{$contig}->{'strand'} = 1;
                  }else{
                    $r->{$contig}->{'oriented'} = 1;
                    $r->{$contig}->{'strand'} = -1;
                  }
                  $nextright_id = $contig;
                }
                last;
              }
            }
					}
				}
				if ($nextright_id){
					print "#autofinisher: contig $nextright_id chosen\n";
					last;
				}
			}
			print "#autofinisher: moving on from contig $middleright_id to $nextright_id\n";
			$middleright_id = $nextright_id;

		}elsif ($middleleft_id){
			### start on the left side, going right
			print "#autofinisher: trying to extend middle left contig $middleleft_id to the right\n";
			my $nextleft_id;
			my @peps = sort {$matrix->{$middleleft_id}->{$b}->{'score'} <=> $matrix->{$middleleft_id}->{$a}->{'score'}} keys(%{$matrix->{$middleleft_id}});
			foreach my $pep (@peps){
				print "#autofinisher: checking peptide: $pep\n";
				if ($r->{$middleleft_id}->{'strand'} == 1){
					### case 1: inner left contig is forward orientation
					if ($matrix->{$middleleft_id}->{$pep}->{'strand'} == 1){
						### case 1a: inner left contig is forward orientation; peptide is in forward orientation
						foreach my $contig (sort {$matrix->{$a}->{$pep}->{'aa'} <=> $matrix->{$b}->{$pep}->{'aa'}} (keys(%{$r}))){
							if ($matrix->{$contig}->{$pep}->{'aa'} > $matrix->{$middleleft_id}->{$pep}->{'aa'}){
								### look for the next contig on the peptide
								if ($r->{$contig}->{'ordered'} == 1){
								}else{
									print "#autofinisher: peptide $pep bridges to contig $contig\n";
									$m->{$contig}->{'number'} = $m->{$middleleft_id}->{'number'} + 1;
									$r->{$contig}->{'ordered'} = 1;
									if ($matrix->{$contig}->{$pep}->{'strand'} == 1){
										$r->{$contig}->{'oriented'} = 1;
										$r->{$contig}->{'strand'} = 1;
									}else{
										$r->{$contig}->{'oriented'} = 1;
										$r->{$contig}->{'strand'} = -1;
									}
									$nextleft_id = $contig;
								}
								last;
							}
						}
					}else{
						### case 1b: inner left contig is forward orientation; peptide is in reverse orientation
						foreach my $contig (sort {$matrix->{$b}->{$pep}->{'aa'} <=> $matrix->{$a}->{$pep}->{'aa'}} (keys(%{$r}))){
							if ($matrix->{$contig}->{$pep}->{'aa'} > 0 && $matrix->{$contig}->{$pep}->{'aa'} < $matrix->{$middleleft_id}->{$pep}->{'aa'}){
								### look for the previous contig on the peptide
								if ($r->{$contig}->{'ordered'} == 1){
								}else{
									print "#autofinisher: peptide $pep bridges to contig $contig\n";
									$m->{$contig}->{'number'} = $m->{$middleleft_id}->{'number'} + 1;
									$r->{$contig}->{'ordered'} = 1;
									if ($matrix->{$contig}->{$pep}->{'strand'} == 1){
										$r->{$contig}->{'oriented'} = 1;
										$r->{$contig}->{'strand'} = -1;
									}else{
										$r->{$contig}->{'oriented'} = 1;
										$r->{$contig}->{'strand'} = 1;
									}
									$nextleft_id = $contig;
								}
								last;
							}
						}
					}
				}else{
					### case 2: inner left contig is reverse orientation
					if ($matrix->{$middleleft_id}->{$pep}->{'strand'} == 1){
						### case 2a: inner left contig is reverse orientation; peptide is in forward orientation
						foreach my $contig (sort {$matrix->{$b}->{$pep}->{'aa'} <=> $matrix->{$a}->{$pep}->{'aa'}} (keys(%{$r}))){
              if ($matrix->{$contig}->{$pep}->{'aa'} > 0 && $matrix->{$contig}->{$pep}->{'aa'} < $matrix->{$middleleft_id}->{$pep}->{'aa'}){
								### look for the previous contig on the peptide
                if ($r->{$contig}->{'ordered'} == 1){
                }else{
									print "#autofinisher: peptide $pep bridges to contig $contig\n";
                  $m->{$contig}->{'number'} = $m->{$middleleft_id}->{'number'} + 1;
                  $r->{$contig}->{'ordered'} = 1;
                  if ($matrix->{$contig}->{$pep}->{'strand'} == 1){
                    $r->{$contig}->{'oriented'} = 1;
                    $r->{$contig}->{'strand'} = -1;
                  }else{
                    $r->{$contig}->{'oriented'} = 1;
                    $r->{$contig}->{'strand'} = 1;
                  }
                  $nextleft_id = $contig;
                }
                last;
              }
            }
					}else{
						### case 2b: inner left contig is reverse orientation; peptide is in reverse orientation
						foreach my $contig (sort {$matrix->{$a}->{$pep}->{'aa'} <=> $matrix->{$b}->{$pep}->{'aa'}} (keys(%{$r}))){
              if ($matrix->{$contig}->{$pep}->{'aa'} > $matrix->{$middleleft_id}->{$pep}->{'aa'}){
								### look for the next contig on the peptide
                if ($r->{$contig}->{'ordered'} == 1){
                }else{
									print "#autofinisher: peptide $pep bridges to contig $contig\n";
                  $m->{$contig}->{'number'} = $m->{$middleleft_id}->{'number'} + 1;
                  $r->{$contig}->{'ordered'} = 1;
                  if ($matrix->{$contig}->{$pep}->{'strand'} == 1){
                    $r->{$contig}->{'oriented'} = 1;
                    $r->{$contig}->{'strand'} = 1;
                  }else{
                    $r->{$contig}->{'oriented'} = 1;
                    $r->{$contig}->{'strand'} = -1;
                  }
                  $nextleft_id = $contig;
                }
                last;
              }
            }
					}
				}
				if ($nextleft_id){
					print "#autofinisher: contig $nextleft_id chosen\n";
					last;
				}
			}
			print "#autofinisher: moving on from contig $middleleft_id to $nextleft_id\n";
			$middleleft_id = $nextleft_id;

		}else{
			my @mkeys = (sort {$m->{$a}->{'number'} <=> $m->{$b}->{'number'}} (keys(%{$m})));
			my $mcount = @mkeys;
			if ($mcount){
				print "#autofinisher: clean up former middle contigs\n";
				my @ids = (sort {$r->{$b}->{'number'} <=> $r->{$a}->{'number'}} (keys(%{$r})));
				my $max = $r->{$ids[0]}->{'number'};
				print "#autofinisher: innermost left contig was number: $max\n";
				foreach my $c (@mkeys){
					$max++;
					$r->{$c}->{'number'} = $max;
					$r->{$c}->{'middle'} = $middle_count;
					if ($mcount == 1){
						$r->{$c}->{'ordered'} = "0";
						$r->{$c}->{'oriented'} = "0";
					}
					print "#autofinisher: innermost left contig is number: $max\n";
					delete $m->{$c};
				}
			}
			print "#autofinisher: look for a new middle contig\n";
			$middle_count++;
			foreach my $contig (sort {$r->{$b}->{'bioseq'}->length <=> $r->{$a}->{'bioseq'}->length} (keys(%{$r}))){
				if ($r->{$contig}->{'number'} == 0){
					print "#autofinisher: new middle contig is $contig\n";
					$middleright_id = $contig;
					$middleleft_id = $contig;
					$m->{$contig}->{'number'} = "0";
					$r->{$contig}->{'ordered'} = 1;
					$r->{$contig}->{'oriented'} = 1;
					$r->{$contig}->{'strand'} = 1;
					$r->{$contig}->{'number'} = 0;
					last;
				}
			}
		}
	}


	### We'll be done when we assign a number to every contig
	### When we check that, it's a good time to give a status update

	while (autofinish_status($r)){



		### some magic happens in here

	}

	### Okay, we just need to give positive numbers to the contigs oriented from the right side

	my @keys = keys(%{$r});
	my $nall = @keys;
	foreach my $id (keys(%{$r})){
		if ($r->{$id}->{'number'} < 0){
			$r->{$id}->{'number'} = $nall + $r->{$id}->{'number'} + 1;
		}
	}

	return $r;
}

#gives a progress report on ordering and orienting contigs.

sub autofinish_status ($) {
	my ($r) = @_;
	my ($nori, $nord, $nall, $nnum);
	foreach my $id (keys(%{$r})){
		$nall++;
		if ($r->{$id}->{'number'} ne 0){
			$nnum++;
		}
		if ($r->{$id}->{'ordered'}){
			$nord++;
		}
		if ($r->{$id}->{'oriented'}){
			$nori++;
		}
	}
	print "#autofinisher: processed $nnum of $nall contigs: $nord ordered and $nori oriented\n";
	if ($nnum == $nall){
		return 0;
	}else{
		return 1;
	}
}

#assign order and orientation to contig for autofinisher

sub make_result_contig ($$$$$) {
	my ($bioseq, $num, $ord, $ori, $str) = @_;
	my $x = {
		'bioseq'	=>	$bioseq,
		'number'	=>	$num,
		'ordered'	=>	$ord,
		'oriented'	=>	$ori,
		'strand'	=>	$str
		};
	return $x;
}

#writes output of autofinisher to fasta file.

sub fasta_output_autofinished ($$) {
	my ($c, $o) = @_;
	my $out = Bio::SeqIO->new(-file => ">$o", '-format' => 'Fasta');

	print "#autofinisher: writing contigs to file: $o\n";

	my @contigs = sort {$c->{$a}->{'number'} <=> $c->{$b}->{'number'}} (keys(%{$c}));
	foreach my $contig (@contigs){
		my $seq;
		if ($c->{$contig}->{'strand'} > 0){
			$seq = $c->{$contig}->{'bioseq'};
		}else{
			$seq = $c->{$contig}->{'bioseq'}->revcom();
		}
		unless ($c->{$contig}->{'ordered'}){
			$c->{$contig}->{'ordered'} = "0";
		}
		unless ($c->{$contig}->{'oriented'}){
			$c->{$contig}->{'oriented'} = "0";
		}
		unless ($c->{$contig}->{'middle'}){
			$c->{$contig}->{'middle'} = "0";
		}

		$seq->primary_id($c->{$contig}->{'number'});
		$seq->display_id($c->{$contig}->{'number'});
		my $description = "Length:".$seq->length()." Ordered:".$c->{$contig}->{'ordered'}." Oriented:".$c->{$contig}->{'oriented'}." Middle:".$c->{$contig}->{'middle'};
		print "#autofinisher: printing $c->{$contig}->{'number'}|$description\n";
		$seq->desc($description);
		$out->write_seq($seq)
	}
}
