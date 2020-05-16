# tRNAscanSE/SprinzlAlign.pm
# This class contains functions to align tRNAs for Sprinzl positioning.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2017 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::SprinzlAlign;

use strict;
use tRNAscanSE::tRNA;
use tRNAscanSE::SprinzlPos;
use tRNAscanSE::Configuration;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(do_sprinzl_pos);

sub do_sprinzl_pos
{
	my ($global_vars, $tRNA, $clade) = @_;
	
	my $sprinzl = $global_vars->{sprinzl};
	my @ar_sprinzl_pos = $sprinzl->sprinzl_pos();
	my $cmalign_out = &align_tRNA($global_vars, $tRNA, $clade);
	my ($alignment, $ss) = &read_cmalign_out($cmalign_out);
	($alignment, $ss) = &fix_alignment($alignment, $ss);
	$tRNA->sprinzl_align($alignment);
	$tRNA->sprinzl_ss($ss);
	$tRNA->map_sprinzl_pos(\@ar_sprinzl_pos);
}	

sub align_tRNA
{
	my ($global_vars, $tRNA, $clade) = @_;
	my $global_constants = $global_vars->{global_constants};

    my $temp_dir = $global_constants->get("temp_dir");
    my $infernal_dir = $global_constants->get("infernal_dir");
	my $cmalign_bin = $infernal_dir."cmalign";	
	my $cm = $global_constants->get_subvalue("sprinzl_cm", lc($clade));
	
	my $temp_fa = "$temp_dir/talign$$".".fa";
	my $temp_out = "$temp_dir/talign$$".".out";

	open(FA_FILE, ">$temp_fa") or die "Fail to open $temp_fa\n";
	print FA_FILE ">".$tRNA->tRNAscan_id()."\n";
	print FA_FILE $tRNA->mat_seq()."\n";
	close(FA_FILE);
	
	my $cmd = $cmalign_bin." -g --notrunc --dnaout ".$cm." ".$temp_fa." > ".$temp_out;
	system($cmd);
	
	return $temp_out;
}

sub read_cmalign_out
{
	my ($cmalign_out) = @_;
	my $line = "";
	my $alignment = "";
	my $ss = "";
	
	open (FILE_IN, "$cmalign_out") or die "Fail to open $cmalign_out\n";
	while ($line = <FILE_IN>)
	{
		chomp($line);
		if ($line !~ /^#/)
		{
			if ($line =~ /^\S+\s+(\S+)$/)
			{			
				$alignment .= $1;
			}
		}
		elsif ($line =~ /^#=GC SS_cons\s+(\S+)$/)
		{
			$ss .= $1;
		}
	}
	close (FILE_IN);
	
	return ($alignment, $ss);
}

sub fix_alignment
{
	my ($alignment, $ss) = @_;
	
	my @bases = split(//, $alignment);
	my $alignment_mod = $alignment;
	my $ss_mod = $ss;
	
	my $last_ins = -999;
	my $ins_start = -1;
	my $ins_count = 0;
	my %ins = ();
	for (my $i = 0; $i < scalar(@bases); $i++)
	{
		if ($bases[$i] =~ /^[acgt]$/)
		{
			if ($last_ins != $i - 1)
			{
				if ($ins_start != -1)
				{
					$ins{$ins_start} = $ins_count;
				}
				$ins_start = $i;
				$ins_count = 0;
			}
			$ins_count++;
			$last_ins = $i;
		}
	}
	if ($ins_start != -1)
	{
		$ins{$ins_start} = $ins_count;
	}
	
	my $del_before_count = 0;
	my $del_after_count = 0;
	foreach $ins_start (sort keys %ins)
	{
		for (my $i = $ins_start - 1; $i >= 0; $i--)
		{
			if ($bases[$i] ne "-")
			{
				last;
			}
			else
			{
				$del_before_count++;
			}
		}
		for (my $i = $ins_start + $ins{$ins_start}; $i < scalar(@bases); $i++)
		{
			if ($bases[$i] ne "-")
			{
				last;
			}
			else
			{
				$del_after_count++;
			}
		}
		if ($del_before_count == $ins{$ins_start})
		{
			for (my $i = $ins_start; $i < $ins_start + $ins{$ins_start}; $i++)
			{
				$bases[$i] = uc $bases[$i];
			}
			splice(@bases, $ins_start - $ins{$ins_start}, $ins{$ins_start});
			$alignment_mod = join("", @bases);
			$ss_mod = substr($ss, 0, $ins_start).substr($ss, $ins_start + $ins{$ins_start});
		}
		elsif ($del_after_count == $ins{$ins_start})
		{
			for (my $i = $ins_start; $i < $ins_start + $ins{$ins_start}; $i++)
			{
				$bases[$i] = uc $bases[$i];
			}
			splice(@bases, $ins_start + $ins{$ins_start}, $ins{$ins_start});
			$alignment_mod = join("", @bases);
			$ss_mod = substr($ss, 0, $ins_start).substr($ss, $ins_start + $ins{$ins_start});
		}
		last;
	}
	return ($alignment_mod, $ss_mod);
}


1;
