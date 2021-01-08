#!/usr/bin/perl

use strict;
use warnings;
use FindBin;
use File::Spec;

package CelegansGeneCoord;

sub new
{
	my ($class) = @_;
	
	my $self = bless {}, $class;

	# open file with gene positions
	open( GENE, File::Spec->catfile($FindBin::Bin,"..","data","test_gene_designations.bed")) || die "Couldn't open file\n";
	my $header = <GENE>;
	
	my %chromosome;
	my %tss;
	my %tts;
	my %strand;

	while( my $line = <GENE> )
	{
		chomp($line);
		my @fields = split /\t/, $line;
		my $acc = "$fields[0]:$fields[1]-$fields[2]";
		my $str = $fields[5];
		if ( $str eq "+" )
		{
			$tss{$acc} = $fields[1];
			$tts{$acc} = $fields[2];
		}
		elsif ( $str eq "-" )
		{
			$tss{$acc} = $fields[2];
			$tts{$acc} = $fields[1];
		} 
		else
		{
			die "No strand information for gene $acc\n";
		}
		
		push @{$chromosome{$fields[0]}}, $acc;

		$strand{$acc} = $str;
	}

	close ( GENE );

	$self->{'chromosome'} = \%chromosome;
	$self->{'tss'} = \%tss;
	$self->{'tts'} = \%tts;
	$self->{'strand'} = \%strand;

	return $self;	
}

sub get_chromosomes
{
	my ($self) = @_;

	return %{$self->{'chromosome'}};

}

sub get_tss
{
        my ($self) = @_;

        return %{$self->{'tss'}};
}

sub get_tts
{
        my ($self) = @_;

        return %{$self->{'tts'}};
}

sub get_strand
{
        my ($self) = @_;

        return %{$self->{'strand'}};
}
	
	

1;	
