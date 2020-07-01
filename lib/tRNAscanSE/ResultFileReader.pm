# tRNAscanSE/ResultFileReader.pm
# This class contains functions to read result files produced by tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2017 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::ResultFileReader;

use strict;
use tRNAscanSE::tRNA;
use tRNAscanSE::ArraytRNA;
use tRNAscanSE::LogFile;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(read_out_file read_ss_file read_sprinzl_pos_map read_noncanomical read_name_map);

sub read_out_file
{
	my ($global_vars, $file) = @_;
	my $log = $global_vars->{log_file};
	my $tRNAs = $global_vars->{tRNAs};
	my $tRNA = undef;
	my %header = ();
	my $intron = {};
	my $line = "";
	my ($startpos, $endpos, $note);
	my @columns = ();
	my $domain_model = "";
	
	$log->status("Reading tRNAs from $file");
	open(OUTFILE, "$file") || die "Error: Fail to open $file\n";	
	while ($line = <OUTFILE>)
	{
		chomp($line);
		if ($line =~ /^Name/)
		{
			$line =~ s/tRNA #/tRNA#/;
		}
		if ($line =~ /^Sequence/)
		{
			$line =~ s/Intron Bounds/Intron\tBound/;
		}

		@columns = split(/\t/, $line, -1);
		for (my $i = 0; $i < scalar(@columns); $i++)
		{
			$columns[$i] = &trim($columns[$i]);
		}
		if ($columns[0] =~ /^Sequence/ || $columns[0] =~ /^Name/ || $columns[0] =~ /^-----/)
		{
			if ($columns[0] =~ /^Sequence/)
			{
				for (my $i = 0; $i < scalar(@columns); $i++)
				{
					if ($columns[$i] eq "Sequence")
					{
						$header{seqname} = $i;
					}
					elsif ($columns[$i] eq "Anti")
					{
						$header{anticodon} = $i;
					}
					elsif ($columns[$i] eq "Intron")
					{
						$header{intron_start} = $i;
						$header{intron_end} = $i+1;
					}
					elsif ($columns[$i] eq "Inf")
					{
						$header{score} = $i;
						$domain_model = "infernal";
					}
					elsif ($columns[$i] eq "Cove")
					{
						$header{score} = $i;
						$domain_model = "cove";
					}
					elsif ($columns[$i] eq "HMM")
					{
						$header{hmm_score} = $i;
					}
					elsif ($columns[$i] eq "2'Str")
					{
						$header{ss_score} = $i;
					}
					elsif ($columns[$i] eq "Hit")
					{
						$header{hit_origin} = $i;
					}
					elsif ($columns[$i] eq "HMM")
					{
						$header{hmm_score} = $i;
					}
				}
			}
			elsif ($columns[0] =~ /^Name/)
			{
				for (my $i = 0; $i < scalar(@columns); $i++)
				{
					if ($columns[$i] eq "tRNA#")
					{
						$header{trna_id} = $i;
					}
					elsif ($columns[$i] eq "Begin" and !defined $header{start})
					{
						$header{start} = $i;
						$header{end} = $i+1;
					}
					elsif ($columns[$i] eq "Type")
					{
						$header{isotype} = $i;
					}
					elsif ($columns[$i] eq "CM")
					{
						$header{isotype_cm} = $i;
						$header{isotype_score} = $i+1;
					}
					elsif ($columns[$i] eq "Note")
					{
						$header{note} = $i;
					}
					elsif ($columns[$i] eq "Count")
					{
						$header{intron_count} = $i;
					}
				}
			}
		}
		else
		{
			$tRNA = tRNAscanSE::tRNA->new;
			$tRNA->seqname($columns[$header{seqname}]);
			$tRNA->tRNAscan_id($columns[$header{seqname}].".trna".$columns[$header{trna_id}]);
			$startpos = $columns[$header{start}];
			$endpos = $columns[$header{end}];
			if ($startpos < $endpos)
			{
				$tRNA->start($startpos);
				$tRNA->end($endpos);
				$tRNA->strand("+");
			}
			else
			{
				$tRNA->end($startpos);
				$tRNA->start($endpos);
				$tRNA->strand("-");
			}
			if ($tRNA->seqname() eq "chrM" or $tRNA->seqname() eq "M" or $tRNA->seqname() eq "chrMT" or $tRNA->seqname() eq "MT")
			{
				$tRNA->category("mt");
			}
			else
			{
				$tRNA->category("cyto");
			}
			$tRNA->isotype($columns[$header{isotype}]);
			$tRNA->anticodon($columns[$header{anticodon}]);
			$tRNA->score($columns[$header{score}]);
			$tRNA->set_domain_model($domain_model, $tRNA->score());
			$tRNA->tRNAscan_id($tRNA->tRNAscan_id()."-".$tRNA->isotype().$tRNA->anticodon());
			if ($tRNA->isotype() eq "Undet")
			{
				$tRNA->category("und");
			}			
			$tRNA->hmm_score($columns[$header{hmm_score}]);
			$tRNA->ss_score($columns[$header{ss_score}]);
			$tRNA->update_domain_model($domain_model, $tRNA->score(), $tRNA->score(), $tRNA->hmm_score(), $tRNA->ss_score());
			
			$note = $columns[$header{Note}];
			if (index($note, "pseudo") > -1)
			{
				$tRNA->pseudo(1);
			}
			
			my $trunc_label = "";
			if ($note =~ /(trunc_start:\d+)/)
			{
				$trunc_label = $1;
			}			
			if ($note =~ /(trunc_end:\d+)/)
			{
				if ($trunc_label ne "")
				{
					$trunc_label .=" ";
				}
				$trunc_label .= $1;
			}
			$tRNA->trunc($trunc_label);
			
			if ($columns[$header{intron_start}] != 0)
			{
				my @intron_starts = split(/\,/, $columns[$header{intron_start}]);
				my @intron_ends = split(/\,/, $columns[$header{intron_end}]);
				for (my $j = 0; $j < scalar(@intron_starts); $j++)
				{
					my ($rel_start, $rel_end) = (0,0); 
					if ($tRNA->strand() eq "+")
					{
						$rel_start = $intron_starts[$j] - $startpos + 1;
						$rel_end = $intron_ends[$j] - $startpos + 1;
					}
					else
					{
						$rel_start = $startpos - $intron_starts[$j] + 1;
						$rel_end = $startpos - $intron_ends[$j] + 1;
					}
					$tRNA->add_intron($rel_start, $rel_end, $intron_starts[$j], $intron_ends[$j], "", "");
				}
			}
			
			$tRNAs->put($tRNA);
		}
	}
	close(OUTFILE);
}

sub read_ss_file
{
	my ($global_vars, $file) = @_;
	my $log = $global_vars->{log_file};
	my $tRNAs = $global_vars->{tRNAs};
	my $tRNA = undef;
	my $intron = {};
	my $line = "";
	my ($startpos, $endpos, $coords);
	
	$log->status("Reading tRNAs from $file");
	open(SSFILE, "$file") || die "Error: Fail to open $file\n";	
	while ($line = <SSFILE>)
	{
		if ($line =~ /^(\S+)\s+\((\S+)\)\s+Length:\s(\d+)\sbp/)
		{
			if (defined $tRNA and $tRNA->tRNAscan_id() ne "")
			{
				$tRNA->set_mature_tRNA();
				$tRNAs->put($tRNA);
			}
			my $id = $1;
			$coords = $2;
			$tRNA = tRNAscanSE::tRNA->new;
			$tRNA->tRNAscan_id($1);
			$tRNA->seqname(substr($tRNA->tRNAscan_id(), 0, rindex($tRNA->tRNAscan_id(), ".")));
			if (index($coords, ",") == -1)
			{
				($startpos, $endpos) = split(/-/, $coords);
				if ($startpos < $endpos)
				{
					$tRNA->start($startpos);
					$tRNA->end($endpos);
					$tRNA->strand("+");
				}
				else
				{
					$tRNA->end($startpos);
					$tRNA->start($endpos);
					$tRNA->strand("-");
				}
			}
			else
			{
				my @loci = split(/\,/, $coords);
				for (my $exon_count = 0; $exon_count < scalar(@loci); $exon_count++)
				{
					($startpos, $endpos) = split(/-/, $loci[$exon_count]);
					if ($startpos < $endpos)
					{
						$tRNA->exon_start($exon_count+1, $startpos);
						$tRNA->exon_end($exon_count+1, $endpos);
						$tRNA->exon_strand($exon_count+1, "+");
					}
					else
					{
						$tRNA->exon_end($exon_count+1, $startpos);
						$tRNA->exon_start($exon_count+1, $endpos);
						$tRNA->exon_strand($exon_count+1, "-");
					}
				}
			}
			if ($tRNA->seqname() eq "chrM" or $tRNA->seqname() eq "M" or $tRNA->seqname() eq "chrMT" or $tRNA->seqname() eq "MT")
			{
				$tRNA->category("mito");
			}
			else
			{
				$tRNA->category("cyto");
			}
		}
		elsif ($line =~ /^Type:\s(\S+)\s+Anticodon:\s(\S+)\sat\s(.+)\s(\(.+\))\s+Score:\s(\S+)/)
		{
			$tRNA->isotype($1);
			$tRNA->anticodon($2);
			my $pos = $3;
			my @ac_pos = split(/\,/, $pos);
			for (my $i = 0; $i < scalar(@ac_pos); $i++)
			{
				if ($ac_pos[$i] =~ /^(\d+)-(\d+)$/)
				{
					$tRNA->add_ac_pos($1, $2);
				}
			}
			$tRNA->score($5);
			$tRNA->set_domain_model("infernal", $tRNA->score());
			my $frag = index($tRNA->isotype(), "-exon");
			if ($frag > -1)
			{
				$tRNA->isotype(substr($tRNA->isotype(), 0, $frag));
			}
			$tRNA->tRNAscan_id($tRNA->tRNAscan_id()."-".$tRNA->isotype().$tRNA->anticodon());
			$tRNA->gtrnadb_id($tRNA->tRNAscan_id());
			if ($tRNA->isotype() eq "Undet")
			{
				$tRNA->category("undetermined_ac");
			}			
		}
		elsif ($line =~ /^HMM Sc=(\S+)\s+Sec struct Sc=(\S+)/)
		{
			$tRNA->hmm_score($1);
			$tRNA->ss_score($2);
			$tRNA->update_domain_model("infernal", $tRNA->score(), $tRNA->score(), $tRNA->hmm_score(), $tRNA->ss_score());
		}
		elsif ($line =~ /pseudogene:\s+HMM Sc=(\S+)\s+Sec struct Sc=(\S+)/)
		{
			$tRNA->hmm_score($1);
			$tRNA->ss_score($2);
			$tRNA->update_domain_model("infernal", $tRNA->score(), $tRNA->score(), $tRNA->hmm_score(), $tRNA->ss_score());
			$tRNA->pseudo(1);
		}
		elsif ($line =~ /^Possible intron: (\d+)-(\d+) \((\d+)-(\d+)\)/)
		{
			$tRNA->add_intron($1, $2, $3, $4, "", "");
		}
		elsif ($line =~ /^Seq:\s(\S+)$/)
		{
			$tRNA->seq($1);
		}
		elsif ($line =~ /^Str:\s(\S+)$/)
		{
			$tRNA->ss($1);
		}
		elsif ($line =~ /^Pre:\s(\S+)$/ or $line =~ /^PRE:\s(\S+)$/ or $line =~ /^BHB:\s(\S+)$/)
		{
			my $temp = $1;
			$temp =~ s/[\[\]]//g;
			$tRNA->mat_seq($tRNA->seq());
			$tRNA->seq($temp);
		}
	}
	if (defined $tRNA and $tRNA->tRNAscan_id() ne "")
	{
		$tRNA->set_mature_tRNA();
		$tRNAs->put($tRNA);
	}
	
	close(SSFILE);
}

sub read_sprinzl_pos_map
{
	my ($global_vars, $file) = @_;
	my $log = $global_vars->{log_file};
	my $tRNAs = $global_vars->{tRNAs};
	my $sprinzl = $global_vars->{sprinzl};
	my @ar_sprinzl_pos = $sprinzl->sprinzl_pos();
	my $line = "";
	my $seq_name = "";
	my $index = -1;
	
	$tRNAs->sort_array("tRNAscan_id");
	
	$log->status("Reading tRNA Sprinzl position map from $file");
	open(FILE_IN, "$file") or die "Error: Fail to open $file\n";
	while ($line = <FILE_IN>)
	{
		chomp($line);
		if ($line !~ /^#/ && $line ne "")
		{
			if ($line =~ /^>(\S+)$/)
			{
				$seq_name = $1;
				$index = $tRNAs->bsearch_id($seq_name, "tRNAscan_id");
				
				$line = <FILE_IN>;
				chomp($line);
				if ($index > -1)
				{
					$tRNAs->get($index)->sprinzl_align($line);
				}

				$line = <FILE_IN>;
				chomp($line);
				if ($index > -1)
				{
					$tRNAs->get($index)->sprinzl_ss($line);
					$tRNAs->get($index)->map_sprinzl_pos(\@ar_sprinzl_pos);
				}

				$line = <FILE_IN>;
				chomp($line);
				if ($index > -1)
				{
					if ($line ne "")
					{
						my @nc = split(//, $line, -1);
						for (my $j = 0; $j < scalar(@nc); $j++)
						{
							if ($nc[$j] eq " ")
							{
								$nc[$j] = "";
							}						
							$tRNAs->get($index)->add_non_canonical($j, $nc[$j]);
						}
					}
					else
					{
						for (my $j = 0; $j < length($tRNAs->get($index)->sprinzl_align()); $j++)
						{
							$tRNAs->get($index)->add_non_canonical($j, "");
						}						
					}
				}
			}
		}
	}
	
	close(FILE_IN);
}

sub read_noncanomical
{
	my ($global_vars, $file) = @_;
	my $log = $global_vars->{log_file};
	my $tRNAs = $global_vars->{tRNAs};
	my $line = "";
	my @columns = ();
	my $index = -1;
	
	$tRNAs->sort_array("tRNAscan_id");
	
	$log->status("Reading tRNA noncanonical features from $file");
	open(FILE_IN, "$file") or die "Fail to open $file\n";

	$line = <FILE_IN>;
	chomp($line);
	my @header = split(/\t/, $line);
	
	while ($line = <FILE_IN>)
	{
		chomp($line);
		@columns = split(/\t/, $line, -1);
		$index = $tRNAs->bsearch_id($columns[0], "tRNAscan_id");
		if ($index > -1)
		{		
			for (my $i = 8; $i < scalar(@columns); $i++)
			{
				$tRNAs->get($index)->add_non_canonical($i - 8, $columns[$i]);
			}
		}
	}

	close(FILE_IN);
}

sub read_name_map
{
	my ($global_vars, $file) = @_;
	my $log = $global_vars->{log_file};
	my $tRNAs = $global_vars->{tRNAs};
	my $line = "";
	my @columns = ();
	my $index = -1;
	
	$tRNAs->sort_array("tRNAscanid");
	
	$log->status("Reading GtRNAdb name map from $file");
	open(FILE_IN, "$file") or die "Fail to open $file\n";

	$line = <FILE_IN>;
	chomp($line);
	my @header = split(/\t/, $line);
	
	while ($line = <FILE_IN>)
	{
		chomp($line);
		@columns = split(/\t/, $line, -1);
		$index = $tRNAs->bsearch_id($columns[0], "tRNAscanid");
		if ($index > -1)
		{		
			$tRNAs->get($index)->gtrnadb_id($columns[1]);
		}
	}

	close(FILE_IN);
}

1;
