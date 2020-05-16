# tRNAscanSE/SprinzlAlignResults.pm
# This class contains functions to write Sprinzl position alignment results.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2017 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::SprinzlAlignResults;

use strict;
use tRNAscanSE::tRNA;
use tRNAscanSE::ArraytRNA;
use tRNAscanSE::LogFile;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(write_sprinzl_pos write_sprinzl_pos_db);

sub write_sprinzl_pos
{
	my ($global_vars, $file, $tRNA) = @_;
	my $log = $global_vars->{log_file};
	
	$log->status("Writing tRNA Sprinzl position $file");
	open(FILE_OUT, ">$file") or die "Error: Fail to open $file\n";
	
	print FILE_OUT "# ".$tRNA->tRNAscan_id()."\n";
	print FILE_OUT "#\n";
	
	my $sprinzl_align = $tRNA->sprinzl_align();
	my @ar = $tRNA->ar_pos_sprinzl_map();
	for (my $i = 0; $i < scalar(@ar); $i++)
	{
		print FILE_OUT $ar[$i]."\t".substr($sprinzl_align, $i, 1)."\n";
	}
	
	close(FILE_OUT);
}

sub write_sprinzl_pos_db
{
	my ($global_vars, $file) = @_;
	my $log = $global_vars->{log_file};
	my $tRNAs = $global_vars->{tRNAs};
	
	$tRNAs->sort_array("tRNAscan_id");
	
	$log->status("Writing tRNA Sprinzl position db $file");
	open(FILE_OUT, ">$file") or die "Error: Fail to open $file\n";

	for (my $i = 0; $i < $tRNAs->get_count(); $i++)
	{
		my $tRNA = $tRNAs->get($i);
		print FILE_OUT "# ".$tRNA->tRNAscan_id()."\n";
		print FILE_OUT "#\n";
		
		my $sprinzl_align = $tRNA->sprinzl_align();
		my @ar = $tRNA->ar_pos_sprinzl_map();
		for (my $i = 0; $i < scalar(@ar); $i++)
		{
			print FILE_OUT $ar[$i]."\t".substr($sprinzl_align, $i, 1)."\n";
		}
		print FILE_OUT "//\n\n";
	}
	
	close(FILE_OUT);
}

1;
