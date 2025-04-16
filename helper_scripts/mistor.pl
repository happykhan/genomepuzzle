#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/basename/;
use File::Temp qw/tempdir/;
use File::Copy qw/cp mv/;
use File::Slurp qw/read_file/;
use List::Util qw/shuffle/;

# Use bioperl goodness to make the gbk files
use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

use threads;

$0=basename $0;

sub logmsg{print STDERR "$0: @_\n"}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s mlstdir=s numcpus=i)) or die $!;
  $$settings{tempdir}||=tempdir("$0.XXXXXX",TMPDIR=>1,CLEANUP=>1);
  $$settings{mlstdir} || die "ERROR: need --mlstdir\n".usage();
  $$settings{numcpus}||=0;

  my($asm) =@ARGV;
  die usage() if(!$asm);

  # Get all MLST loci but also shuffle them to help with
  # thread balance.
  my @mlstLocus = shuffle glob("$$settings{mlstdir}/*.fna");
  my $numLoci = scalar(@mlstLocus);
  if(!$numLoci){
    die "ERROR: no loci were found in $$settings{mlstdir}";
  } elsif($numLoci < 7){
    logmsg "WARNING: there are only $numLoci loci in this scheme.";
  }
  
  # Make the blast database
  my $blastdb="$$settings{tempdir}/assembly.fna";
  cp($asm,$blastdb) or die "ERROR: could not copy $asm to $blastdb: $!";
  system("makeblastdb -dbtype nucl -in $blastdb 1>&2"); die if $?;
  
  my $numLociPerThread = int($numLoci / $$settings{numcpus}) + 1;
  my @thr;
  for(0..$$settings{numcpus}-1){
    my @threadLocus = splice(@mlstLocus, 0, $numLociPerThread);
    $thr[$_] = threads->new(\&annotationWorker, $blastdb, \@threadLocus, $settings);
    logmsg scalar(@threadLocus)." loci being compared in thread ".$thr[$_]->tid;
  }

  # Go ahead and read in the assembly
  my %seq;
  my $seqin=Bio::SeqIO->new(-file=>$asm);
  while(my $seq=$seqin->next_seq){
    # Add a source tag
    my $feature=Bio::SeqFeature::Generic->new(
      -start       =>1,
      -end         =>$seq->length,
      -primary     =>"source",
    );
    $seq->add_SeqFeature($feature);

    $seq{$seq->id}=$seq;
  }
  $seqin->close;

  # Join threads and add the sequence features as we go.
  my @feature;
  for(my $t=@thr-1;$t>=0;$t--){
    my $tmp = $thr[$t]->join;

    my $percentDone = sprintf("%0.2f",100-($t+1)/$$settings{numcpus} * 100);
    logmsg "Adding features for ".scalar(@$tmp)." loci ($percentDone% done)";
    push(@feature, @$tmp);
  }

  for my $feature(sort {$a->start <=> $b->start} @feature){
    my $seqid = $feature->seq_id;
    $seq{$seqid}->add_SeqFeature($feature);
  }
  logmsg "Done adding all ".scalar(@feature)." features.";

  # Print the genbank file
  logmsg "Printing to stdout";
  my $outseq=Bio::SeqIO->new(-format=>"genbank");
  # I didn't realize it but the features have to be added in order
  for my $seq( sort{ $b->length <=> $a->length } values(%seq)){
    $outseq->write_seq($seq);
  }

  return 0;
}

sub annotationWorker{
  my($blastdb, $loci, $settings)=@_;
  my $numLoci = scalar(@$loci);

  my @feature;

  my $blastResult = "$$settings{tempdir}/blast".threads->tid.".tsv";
  # truncate and test for r/w
  open(my $blastResultFh, ">", $blastResult) or die "ERROR: could not write to $blastResult: $!";

  my $resultCounter=0;
  for(my $i=0;$i<$numLoci;$i++){
    my $locusname=basename($$loci[$i], qw(.fna));
    my $blastResult="$$settings{tempdir}/$locusname.tsv";

    # Run the blast and sort by highest score
    my $bestHit = `blastn -query $$loci[$i] -db $blastdb -num_threads 1 -outfmt 6 | sort -k12,12n | head -n 1`;
    if($?){
      die "ERROR blasting $$loci[$i] vs $blastdb: $!";
    }
    print $blastResultFh $bestHit;
  }
  close $blastResultFh;

  # Create seq features
  open(my $blsFh, $blastResult) or die "ERROR: could not read $blastResult: $!";
  while(<$blsFh>){
    chomp;
    my($allele,$contig,$identity,$alnLength,$mismatch, $gap, $qstart, $qend, $sstart, $send, $evalue, $score)=split /\t/;
    my $locusname=$allele;

    # check strand
    my $strand=1;
    my($feat_start,$feat_end)=($sstart,$send);
    if($sstart > $send){
      $strand=-1;
      ($feat_start,$feat_end)=($send,$sstart)
    }
    
    my $feature=Bio::SeqFeature::Generic->new(
      -start         =>$feat_start,
      -end           =>$feat_end,
      -strand        =>$strand,
      -primary       =>"gene",
      -source_tag    => $0,
      -display_name  =>$locusname,
      -score         =>$score,
      -seq_id        =>$contig,
      -tag           =>{ 
                         locus=>$locusname,
                         gene =>$allele,
                         note =>"Evidence:$0",
                       },
    );

    push(@feature, $feature);

  }

  @feature = sort{$a->start <=> $b->start} @feature;

  return \@feature;
}

sub usage{
  "$0: annotate a genome quickly using blast and an MLST scheme.
  Usage: $0 assembly.fasta --mlstdir MLST > annotation.gbk
  --numcpus         1      number of cpus to use
  --mlstdir         ''     Location of MLST files. All 
                           files must have .fna extension.
  --tempdir         ''     Location of temporary files (optional)
  "
}