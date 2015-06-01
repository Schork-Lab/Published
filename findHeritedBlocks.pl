#!/usr/bin/perl
# findHeritedBlocks.pl
#  made to find blocks of variants inherited by two cousins from
#  their mothers (who are sisters)
# Original Code by Tristan M. Carland, PhD (TSRI-STSI, March 2012)
use warnings; use strict;

# command line parameter processing
unless (scalar @ARGV == 2)
{
  print STDERR "$0 - used to find inherited haplotype blocks\n";
  print STDERR "Usage: $0 cg_phased_data outFile\n";
  exit(0);
}

# input and output files
my $dataFile  = $ARGV[0];
my $cataFile  = $ARGV[1];
open(INFILE,   "$dataFile") or die "Could not open $dataFile, $!\n";
open(OUTFILE, ">$cataFile") or die "Could not open $cataFile, $!\n";

# these were determined to yield an appropriate amount of shared haplotype blocks
my $minLength   = 600;  # minimum length of set
my $minPercent  = .99;  # threshold percentage to be a set

# other variables for later use
my @win       = ();   # stores the current window (1,0)
my @winData   = ();   # stores the current window data
my $inSet     = 0;    # are we in a set (fake boolean)
my $winTally  = 0;    # running sum of the window
my @set       = ();   # a running set that meets criterion
my $numRegChr = 0;    # tally of the number of regions per chr
my $curChr    = "";   # the current chromosome
my $percent   = 0;

print "#Chr\tStart\tStop\tLength\tHits\tNot\t%\tDensity\n";

# read through the entire file, line by line
while (defined(my $line = <INFILE>))
{
  chomp($line);                 # remove spaces and such
  next unless($line =~ /\S/);   # skip ahead unless not empty
  next if($line =~ /^\#/);      # skip comment lines
  next if($line =~ /MISSING/);  # skip missing data

  # first split by tabs, location then each cousin/child
  my @lineData = split /\t/, $line;

  # 0 vartype 1 chrom 2 start 3 end 4 ref-allele 5 variant-allele
  my @loc = split /\s/, $lineData[0];

  # 0 child-maternal-haplotype 1 child-paternal-haplotype 2 child-genotypes
  # 3 mother-genotypes 4 father-genotypes 5 TRIO-status 
  my @c1  = split /\s/, $lineData[1];
  my @c2  = split /\s/, $lineData[2];

  # skip missing data
  next if( ($c1[0] eq "N") || ($c2[0] eq "N") );
  
############################################################
############ Filters complete, begin processing ############
############################################################
  ### In the event of a new chromosome, reset window/set stats
  if($loc[1] ne $curChr)
  {
    print STDERR $curChr."\t".$numRegChr."\n";
    printSet(\@set);
    @win        = ();   # stores the current window (1,0)
    @winData    = ();   # stores the current window data, all of it
    $winTally   = 0;    # running sum of the window
    $inSet      = 0;    # are we in a set (fake boolean)
    @set        = ();   # a growable window that meets criterion
    $curChr     = $loc[1];
    $numRegChr  = 0;    # tally of the number of regions per chr
  }

  ## when the window has reached it's full size
  if( scalar(@win) == $minLength )
  {
    # remove the bottom of the stacks
    $winTally = $winTally - $win[0];
    shift(@win);
    shift(@winData);
  }

  # update the stacks
  if($c1[0] eq $c2[0])
  {
    push(@win, 1);
    $winTally++;
  }
  else
  {
    push(@win, 0);
  }
  push(@winData, $line);
  #################################
  ### Haplotype Threshold Check ###
  $percent = $winTally/scalar(@win);
  if(($percent >= $minPercent) && (scalar(@win) >= $minLength))
  {
    if($inSet)  # if still in a set, just add to it
    {
      push(@set, $line);
    }
    elsif( getStart(\@winData) <= getStop(\@set) )
    {
      # this set overlaps with the last one, add to it
      $inSet = 1;
      @set = mergeSets(\@set, \@winData);
#     push(@set, @winData);
    } 
    else # must be a new set, add the entire window, print last set
    { 
      $numRegChr++;
      printSet(\@set);
      @set = ();
      push(@set, @winData);
      $inSet = 1;
    } 
  }   
  else
  {   
    if($inSet)  # if we've fallen below thresh while in a window
    {
      $inSet = 0;
    }
  }   
}   
if($inSet)  # print any remaining sets
{     
  $numRegChr++;
  printSet(\@set);
} 
print STDERR $curChr."\t".$numRegChr."\n";
close INFILE;
close OUTFILE;

#############################################################
############ Subroutines for printing of Regions ############
#############################################################
sub printSet
{
  my @set   = @{$_[0]};
# my $file  = $_[1];

  if(scalar(@set) > 2)
  {
    my $chr   = getChr(\@set);
    my $start = getStart(\@set);
    my $stop  = getStop(\@set);
    my $len   = $stop-$start;
    my $hits  = getHits(\@set);
    my $miss  = scalar(@set)-$hits;
    my $score = $hits/scalar(@set);
    my $dens  = $hits/$len;
    
    print   "$chr\t$start\t$stop\t$len\t$hits\t$miss\t";
    printf  "%.2f\t", $score;
    printf  "%.5f\n", $dens;
    
    printCatalog(\@set);
  }
}   
#########################################################
# Finds the number of hits/matches in a set
sub printCatalog
{
  my @set       = @{$_[0]};
# my $cataFile  = $_[1];
# open(OUTFILE, ">$cataFile") or die "Could not open $cataFile, $!\n";

  foreach (@set)
  {
    # first split by tabs, location then each cousin/child
    my @lineData = split /\t/, $_;

    # 0 vartype 1 chrom 2 start 3 end 4 ref-allele 5 variant-allele
    my @loc = split /\s/, $lineData[0];

    # 0 child-maternal-haplotype 1 child-paternal-haplotype 2 child-genotypes
    # 3 mother-genotypes 4 father-genotypes 5 TRIO-status 
    my @c1  = split /\s/, $lineData[1];
    my @c2  = split /\s/, $lineData[2];

    if($c1[0] eq $c2[0])
    {
      # print the shared variant locations to the catalog file
      print OUTFILE "$loc[1]\t$loc[2]\t$loc[3]\n";
    }
  }
}
#########################################################
# Designed to merge two overlapping sets into one
sub mergeSets
{
  my @set1 = @{$_[0]};
  my @set2 = @{$_[1]};

  while($set1[-1] ne $set2[0])
  {
    pop(@set1);
  
    if( scalar(@set1) == 0 )  # cheap fix to an infrequent error
    {
      return (@{$_[0]}, @{$_[1]});
    }
  } 
  pop(@set1);
  return (@set1,@set2);
} 
#########################################################
# Finds the number of hits/matches in a set
sub getHits
{
  my @set   = @{$_[0]};
  my $count = 0;

  foreach (@set)
  {
    # first split by tabs, location then each cousin/child
    my @lineData = split /\t/, $_;

    # 0 vartype 1 chrom 2 start 3 end 4 ref-allele 5 variant-allele
    my @loc = split /\s/, $lineData[0];

    # 0 child-maternal-haplotype 1 child-paternal-haplotype 2 child-genotypes
    # 3 mother-genotypes 4 father-genotypes 5 TRIO-status 
    my @c1  = split /\s/, $lineData[1];
    my @c2  = split /\s/, $lineData[2];

    if($c1[0] eq $c2[0])
    {
      $count++; # keep count of our matches/hits
    }
  }
  return $count;
}
#########################################################
# Returns the start location of the first match/hit
sub getStart
{
  my @set   = @{$_[0]};

  foreach (@set)
  {
    # first split by tabs, location then each cousin/child
    my @lineData = split /\t/, $_;
    
    # 0 vartype 1 chrom 2 start 3 end 4 ref-allele 5 variant-allele
    my @loc = split /\s/, $lineData[0];
    
    # 0 child-maternal-haplotype 1 child-paternal-haplotype 2 child-genotypes
    # 3 mother-genotypes 4 father-genotypes 5 TRIO-status 
    my @c1  = split /\s/, $lineData[1];
    my @c2  = split /\s/, $lineData[2];
    
    if($c1[0] eq $c2[0])
    {
      return $loc[2];
    }
  }   
}   
#########################################################
# Finds the stop location of the last match/hit
sub getStop
{
  unless(@{$_[0]})
  {
    return 0;
  }
  my @set   = @{$_[0]};

  # iterate through the list backwards until you hit a match
  for (my $i = scalar(@set)-1; $i >= 0; $i--)
  {
    # first split by tabs, location then each cousin/child
    my @lineData = split /\t/, $set[$i];

    # 0 vartype 1 chrom 2 start 3 end 4 ref-allele 5 variant-allele
    my @loc = split /\s/, $lineData[0];
    
    # 0 child-maternal-haplotype 1 child-paternal-haplotype 2 child-genotypes
    # 3 mother-genotypes 4 father-genotypes 5 TRIO-status 
    my @c1  = split /\s/, $lineData[1];
    my @c2  = split /\s/, $lineData[2];
    
    if($c1[0] eq $c2[0])
    {
      return $loc[3];
    }
  }   
  return 0;
} 
#########################################################
# Returns the chrom of the first entry
sub getChr
{
  my @set = @{$_[0]};

  # first split by tabs, location then each cousin/child
  my @lineData = split /\t/, $set[0];

  # 0 vartype 1 chrom 2 start 3 end 4 ref-allele 5 variant-allele
  my @loc = split /\s/, $lineData[0];

  return $loc[1];
}
