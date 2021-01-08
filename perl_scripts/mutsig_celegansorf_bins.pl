#!/usr/bin/perl

use strict;
use warnings;
use FindBin;
use lib $FindBin::Bin;

use CelegansGeneCoord;

# ask for probe filename to analyze
print STDERR "Enter filename of sorted mutation bedfile\n";
my $mutbedfile = <STDIN>;
chomp($mutbedfile);

my $outputfile = $mutbedfile;
$outputfile =~ s/\.bed/_TrxStrandBins\.txt/ || die "Wrong file type!\n";

open (MUT, "$mutbedfile" ) || die "Couldn't open file: $mutbedfile\n";

print STDERR "Loading Gene coordinates\n";
my $genes = CelegansGeneCoord->new();

my %tsbin = ();
my %ntsbin = ();
my %intergenicbin = ();

my %trinucref = ();
my %subref = ();

my %chromosomes = $genes->get_chromosomes();
my %trxstart = $genes->get_tss();
my %trxend = $genes->get_tts();
my %trxstrand = $genes->get_strand();

my $chr = "";
my %genelookup;
my $mitomutcount = 0;
while ( my $line = <MUT> )
{
	chomp $line;
	my @fields = split /\t/, $line;
	if ( $fields[0] =~ /(chr[IXV]+)/ )
	{
		my $temp = $1;
		if ( $temp ne $chr )
		{
                	print STDERR "Starting to process $temp\n";
			$chr = $temp;
			%genelookup = ();
			my $skipped = 0;
			foreach my $acc ( @{$chromosomes{$chr}} )
			{
				my $strand = $trxstrand{$acc};
                		if ( exists $trxend{$acc} )
                		{
					;
                		}
                		else
                		{
					$skipped++;
                        		next;
                		}

                                my $start;
                                my $end;
				if ( $strand eq "+" )
				{
					$start = $trxstart{$acc};
					$end = $trxend{$acc};
				}
				elsif ( $strand eq "-" )
				{
					$start = $trxend{$acc};
					$end = $trxstart{$acc};
				}
				else
				{
					die "No strand info!\n";
				}
		
				if ( $start >= $end )
				{
					die "Error with gene coords!\n";
				}
	
				for ( my $i = $start; $i <= $end; $i++ )
				{
					if ( exists $genelookup{$i} )
					{
						if ( $genelookup{$i} eq "+" || $genelookup{$i} eq "-" || $genelookup{$i} eq "AMBIG" )
						{
							if ( $genelookup{$i} ne $strand )
							{
								$genelookup{$i} = "AMBIG";
							}
						}
						else
						{
							die "Weird gene lookup: $genelookup{$i}\n";
						}
					}
					else
					{	
						$genelookup{$i} = $strand;
					}
				} 			
			}

			print STDERR "$skipped skipped genes\n";
		}
		my $mutpos = $fields[2];  #use 1-based end coord;
		if ( $mutpos != $fields[1] + 1 )
		{
			die "Error in mut coordinates for line: $line\n";
		}
		my $trinuc = $fields[3];
		my $mutbase = $fields[4];
		my $mutstrand = $fields[5];
		if ( $mutstrand ne "+" && $mutstrand ne "-" )
		{
			die "Error in mutation strand in line: $line\n";
		}
		my $refbase = substr $trinuc, 1, 1; 
		my $substitution = $refbase . ">" . $mutbase;

		$trinucref{$trinuc} = 1;
		$subref{$substitution} = 1;
		if ( exists $genelookup{$mutpos} )
		{			
			if ( $genelookup{$mutpos} ne "+" && $genelookup{$mutpos} ne "-" )
			{
				if ( $genelookup{$mutpos} eq "AMBIG" )	
				{
					$intergenicbin{$trinuc}{$substitution}++;
				}
				else
				{
					die "Error with gene strand\n";
				}
			}
			elsif ( $genelookup{$mutpos} eq $mutstrand )
			{
				$ntsbin{$trinuc}{$substitution}++; # Shouldn't this be the opposite?
			}	
			else
			{
				$tsbin{$trinuc}{$substitution}++;					
			}
		}
		else
		{
			$intergenicbin{$trinuc}{$substitution}++;
		}
	}
	elsif ( $fields[0] =~ /chrM/ )
	{
		$mitomutcount++;
	}
	else
	{
		die "Error in bed file line: $line\n";
	}
}

open (OUT, ">$outputfile") || die "couldn't open file\n";
#print header
print OUT "Data from file: $mutbedfile\n";

print OUT "Trinucleotide context\tMutant base\tMutation substitution\tTS count\tNTS count\tIntergenic count\n";
foreach my $tri ( sort keys %trinucref )
{
	my $trirefbase = substr $tri, 1, 1;
	foreach my $sub ( sort keys %subref ) 
	{
		my $refbase = substr $sub, 0, 1;		
		my $mutantbase = substr $sub, 2, 1;

		# skip if middle base of trinuc is not same as current sub reference base
		if ( $refbase ne $trirefbase )
		{
			next;
		}

		my $ts = 0;
		my $nts = 0;
		my $int = 0;
		if ( exists $tsbin{$tri}{$sub} )
		{
			$ts = $tsbin{$tri}{$sub};
		}
		if ( exists $ntsbin{$tri}{$sub} )
		{
			$nts = $ntsbin{$tri}{$sub};
		}
                if ( exists $intergenicbin{$tri}{$sub} )
                {
                        $int = $intergenicbin{$tri}{$sub};
                }
                print OUT "$tri\t$mutantbase\t$sub\t$ts\t$nts\t$int\n";
	}
}

