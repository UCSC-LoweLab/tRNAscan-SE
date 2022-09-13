# tRNAscanSE/ScanResult.pm
# This class describes the outputs of scan results used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2017 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::ScanResult;

use strict;
use tRNAscanSE::Utils;
use tRNAscanSE::Sequence;
use tRNAscanSE::tRNA;
use tRNAscanSE::ArraytRNA;
use tRNAscanSE::IntResultFile;
use tRNAscanSE::MultiResultFile;
use tRNAscanSE::Options;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(save_Acedb_from_firstpass write_tRNA write_tRNAs output_tRNA write_bed write_gff output_split_fragments print_results_header);

our $printed_header = 0;            # keeps track of whether or
                                    # or not results column header
                                    # has been printed yet
our ($max_seq_name_width, $max_seq_len_width);


sub save_Acedb_from_firstpass
{
    my ($output_codon, $r_one_let_trans_map, $fp_tRNAs, $out_file) = @_;
    my $triplet = "";

    &open_for_append(\*FILE_H, $out_file);

	for (my $i = 0; $i < $fp_tRNAs->get_count(); $i++)
	{
		my $trna = $fp_tRNAs->get($i);
		
		printf FILE_H "Sequence\t%s\nSubsequence\t%s.t%d %d %d\n\n",
			$trna->seqname(), $trna->seqname(), $i + 1, $trna->start(), $trna->end();		
		printf FILE_H "Sequence\t%s.t%d\nSource\t\t%s\n", $trna->seqname(), $i + 1, $trna->seqname();

		if ($trna->get_intron_count() > 0)
		{
			my @ar_introns = $trna->ar_introns();
			
			if ($ar_introns[0]->{rel_start} < $ar_introns[0]->{rel_end})
			{
				printf FILE_H "Source_Exons\t1 %d\n", ($ar_introns[0]->{rel_start} - $trna->start());
				printf FILE_H "Source_Exons\t%d %d\n",
					$ar_introns[0]->{rel_end} - $trna->start() + 2,
					$trna->end() - $trna->start() + 1;
			}
			else
			{
				printf FILE_H "Source_Exons\t1 %d\n", ($trna->start() - $ar_introns[0]->{rel_start} + 1);
				printf FILE_H "Source_Exons\t%d %d\n",
					$trna->start() - $ar_introns[0]->{rel_end} + 2,
					$trna->start() - $trna->end() + 1;
			}
		}	 
		printf FILE_H "Brief_identification tRNA-%s\n", $trna->isotype();
		
		# either output Codon or Anticodon for tRNA
		$triplet = uc($trna->anticodon());
		if ($output_codon)
		{
			$triplet = &rev_comp_seq($triplet);
		}
	
		printf FILE_H "Transcript tRNA \"%s %s %s\"\n\n",
			$triplet, $trna->isotype(), $r_one_let_trans_map->{$trna->isotype()};
		
    }
    close(FILE_H);
}

sub print_results_header
{    
    my ($TABOUT, $opts, $get_hmm_score, $seq_name_width, $seq_len_width, $fp) = @_;
    my ($label, $codon_label) = "";
    
    if ($opts->cove_mode())
	{
		$label = "\tCove";
    }
    elsif ($opts->infernal_mode())
	{
		$label = "\tInf";
    }
    elsif ($opts->eufind_mode() && !$opts->tscan_mode())
	{
		$label = "\tEufind";
    }

    if ($opts->output_codon())
	{
		$codon_label = "   "; 
    }
    else
	{
		$codon_label = "Anti";
    }
    
    if (!($opts->ace_output()))
	{
		printf {$$TABOUT} "%-".$seq_name_width."s\t\t","Sequence";
		printf {$$TABOUT} "%-".$seq_len_width."s\t","tRNA";
		printf {$$TABOUT} "%-".$seq_len_width."s\t","Bounds";
		print  {$$TABOUT} "tRNA\t$codon_label\tIntron Bounds",$label;
	
		if ($get_hmm_score)
		{ 
			print {$$TABOUT} "\tHMM\t2'Str";
		}
		if ($opts->infernal_score())
		{
			print {$$TABOUT} "\tInf";
		}
		if ($opts->save_source())
		{
			print {$$TABOUT} "\tHit";
		}		
		if (!$fp and ((($opts->euk_mode() or $opts->bact_mode() or $opts->arch_mode()) and !$opts->no_isotype() and $opts->detail()) or $opts->metagenome_mode()))
		{
			print {$$TABOUT} "\tIsotype\tIsotype";
			if ($opts->euk_mode() and $opts->mito_model() ne "")
			{
				print {$$TABOUT} "\tType";
			}			
		}
		print {$$TABOUT} "\t      ";
		print {$$TABOUT} "\n";

		printf {$$TABOUT} "%-".$seq_name_width."s\t","Name";
		print  {$$TABOUT} "tRNA \#\t";
		printf {$$TABOUT} "%-".$seq_len_width."s\t","Begin";
		printf {$$TABOUT} "%-".$seq_len_width."s\t","End";
	
		print {$$TABOUT} "Type\tCodon\tBegin\tEnd\tScore";
	
		if  ($get_hmm_score)
		{ 
			print {$$TABOUT} "\tScore\tScore";
		}
		if ($opts->infernal_score())
		{
			print {$$TABOUT} "\tScore";
		}
		if ($opts->save_source())
		{
			print {$$TABOUT} "\tOrigin";
		}
		if (!$fp and ((($opts->euk_mode() or $opts->bact_mode() or $opts->arch_mode()) and !$opts->no_isotype() and $opts->detail()) or $opts->metagenome_mode()))
		{
			print {$$TABOUT} "\tCM\tScore";
			if ($opts->euk_mode() and $opts->mito_model() ne "")
			{
				print {$$TABOUT} "\t         ";
			}			
		}
		print {$$TABOUT} "\tNote";
		print {$$TABOUT} "\n";
	
		printf {$$TABOUT} "%-".$seq_name_width."s\t","--------";
		print  {$$TABOUT} "------\t";
		printf {$$TABOUT} "%-".$seq_len_width."s\t","-----";
		printf {$$TABOUT} "%-".$seq_len_width."s\t","------";
		print  {$$TABOUT} "----\t-----\t-----\t----\t------";
	
		if  ($get_hmm_score)
		{ 
			print {$$TABOUT} "\t-----\t-----";
		}
		if ($opts->infernal_score())
		{
			print {$$TABOUT} "\t-----";
		}
		if ($opts->save_source())
		{
			print {$$TABOUT} "\t------";
		}
		if (!$fp and ((($opts->euk_mode() or $opts->bact_mode() or $opts->arch_mode()) and !$opts->no_isotype() and $opts->detail()) or $opts->metagenome_mode()))
		{
			print {$$TABOUT} "\t-------\t-------";
			if ($opts->euk_mode() and $opts->mito_model() ne "")
			{
				print {$$TABOUT} "\t---------";
			}			
		}
		print {$$TABOUT} "\t------";
		print {$$TABOUT} "\n";
    }
}

sub write_tRNA
{    
    my ($file_name, $seq_name, $seq_desc, $seq, $overwrite) = @_;
    
    my $trna_file = tRNAscanSE::Sequence->new;
    my $write_mode = "append";
    if ($overwrite) {
		$write_mode = "write";
    }
    $trna_file->set_seq_info($seq_name, $seq_desc, length($seq), $seq);
    $trna_file->open_file($file_name, $write_mode);
    $trna_file->write_fasta();
    $trna_file->close_file();
}

sub write_tRNAs
{    
    my ($file_name, $sp_int_results, $mat_seq, $overwrite, $model) = @_;
    
	my $count = 0;
    my $trna_file = tRNAscanSE::Sequence->new;
    my $write_mode = "append";
    if ($overwrite)
	{
		$write_mode = "write";
    }
	my $seq = "";
	my $tRNA = tRNAscanSE::tRNA->new;
    $trna_file->open_file($file_name, $write_mode);
	if ($sp_int_results->open_file("read"))
	{
        my @record_indexes = $sp_int_results->get_indexes();
        
        for (my $i = 0; $i < scalar(@record_indexes); $i++)
        {
			$sp_int_results->get_tRNA($record_indexes[$i]->[0], $tRNA);
			if ($model eq "" or $tRNA->model() eq $model)
			{
				if ($mat_seq)
				{
					$seq = $tRNA->mat_seq();
				}
				else
				{
					$seq = $tRNA->seq();
				}
				$trna_file->set_seq_info($tRNA->seqname().".t".&pad_num($tRNA->id(), 6), "", length($seq), $seq);
				$trna_file->write_fasta();
				$count++;
			}
		}
	}
    $trna_file->close_file();
	
	return $count;
}

# Write final tRNA prediction to various selected output sources/files
# Sets globals $MaxSeqNameWidth and $MaxSeqLenWidth and $printed_header

sub output_tRNA
{
    my ($global_vars, $cm, $r_tab_results, $get_hmm_score, $program_id) = @_;         
	my $opts = $global_vars->{options};
	my $log = $global_vars->{log_file};
	my $gc = $global_vars->{gc};
    my $sp_int_results = $global_vars->{sp_int_results};
    my $iso_int_results = $global_vars->{iso_int_results};
    my $tRNA = tRNAscanSE::tRNA->new; 
	
    my $results_line = "";
	my $isotype_line = "";
	my ($iso_models, $mito_models);

    my @sp_indexes = $sp_int_results->get_indexes();
    if ($sp_int_results->open_file("read"))
    {
		if (!$opts->no_isotype() and $iso_int_results->open_file("read"))
		{
			($iso_models, $mito_models) = &get_models($opts, $cm);
			$iso_int_results->read_models();
		}
		if ($opts->isotype_specific_file() ne "")
		{
			&open_for_append(\*ISOTYPE, $opts->isotype_specific_file());
		}
		if ($opts->save_all_struct())
		{
			&open_for_append(\*SECSTRUCT, $opts->all_struct_file());
		}
		if ($opts->ace_output())
		{			
		    &open_for_append(\*ACEOUT, $opts->out_file());
		}
		else
		{
			&open_for_append(\*TABOUT, $opts->out_file());
		}
		if ($opts->output_fasta_file() ne "")
		{
			&open_for_append(\*FASTA, $opts->output_fasta_file());
		}
		
		for (my $i = 0; $i < scalar(@sp_indexes); $i++)
		{
			$sp_int_results->get_tRNA($sp_indexes[$i]->[0], $tRNA);
			if (!$opts->no_isotype())
			{
				$iso_int_results->get_next_tRNA($tRNA);
			}

			my ($type, $model, $score, $ss) = $tRNA->get_highest_score_model();
			if ($tRNA->isotype() eq "Met" and $type eq "cyto" and ($model eq "iMet" or $model eq "fMet" or $model eq "Ile2"))
			{
				$tRNA->isotype($model);
			}
			elsif ($tRNA->isotype() eq "Met" and $type eq "cyto" and $model ne "Met" and $model ne "iMet" and $model ne "fMet")
			{
				$tRNA->sort_multi_models("model");
				my ($met_iso_model, $met_iso_score, $met_iso_ss) = $tRNA->get_model_hit("cyto", $tRNA->isotype());
				my ($ile2_iso_model, $ile2_iso_score, $ile2_iso_ss) = $tRNA->get_model_hit("cyto", "Ile2");
				if ($ile2_iso_score > 0 and $met_iso_score > 0)
				{
					if (($score - $ile2_iso_score) <= 5 and ($ile2_iso_score - $met_iso_score) >= 5 and $tRNA->score() > 50)
					{
						$tRNA->isotype("Ile2");
					}
				}
			}
			
			if (!$opts->results_to_stdout())
			{
				$log->broadcast($tRNA->tRNAscan_id().":  ".$opts->second_pass_label()." type= ".$tRNA->isotype()."\t ".
				"Score= ".$tRNA->score());
			}
			if ($opts->save_all_struct())
			{
				&save_allStruct_output(\*SECSTRUCT, $opts, $gc, $get_hmm_score, $tRNA);
			}
			if ($opts->output_fasta_file() ne "")
			{
				&write_tRNA_sequence(\*FASTA, $tRNA);
			}
			    
			# Create tabular results line, ready for output    
			if (!$printed_header)
			{
				$max_seq_name_width = &max(length($tRNA->src_seqid()) + 1, 8);
				$max_seq_len_width  = length($tRNA->src_seqlen());
			}    
		    $results_line = &construct_tab_output($opts, $get_hmm_score, $tRNA, $max_seq_name_width, $max_seq_len_width);
			
			# Internal copy of results saved for later uses
			push(@$r_tab_results, $results_line);
			
			if ($opts->ace_output())
			{       
				&save_Acedb_from_secpass(\*ACEOUT, $opts, $gc, $tRNA, $program_id);
			}
			else 
			{    
				if (!($opts->brief_output() || $printed_header))
				{
					&print_results_header(\*TABOUT, $opts, $get_hmm_score, $max_seq_name_width, $max_seq_len_width, 0);
					if ($opts->isotype_specific_file() ne "")
					{
						&print_isotype_specific_header(\*ISOTYPE, $opts, $iso_models, $mito_models);
					}
					
					$printed_header = 1;
				}	    
				print TABOUT $results_line;
				
				if ($opts->isotype_specific_file() ne "")
				{
					$isotype_line = &construct_isotype_specific_output($opts, $iso_models, $mito_models, $tRNA);
					print ISOTYPE $isotype_line;
				}
			}
		}
		
		if ($opts->ace_output())
		{
			close ACEOUT;
		}
		else
		{
			close TABOUT;
		}
		if ($opts->save_all_struct())
		{
			close SECSTRUCT;
		}
		if ($opts->isotype_specific_file() ne "")
		{
			close ISOTYPE;
		}
		if (!$opts->no_isotype())
		{
			$iso_int_results->close_file();
		}
		if (!$opts->output_fasta_file() ne "")
		{
			close FASTA;
		}
        $sp_int_results->close_file();		
	}
}

sub save_allStruct_output
{    
    my ($SECSTRUCT, $opts, $gc, $get_hmm_score, $tRNA) = @_;

    my $ruler = '    *    |' x 20;     # ruler printed out with
                                       #  secondary structure output

    my $seqlen = length($tRNA->seq());
	my ($type, $model, $score, $ss) = $tRNA->get_highest_score_model();
    
	if ($tRNA->strand() eq "+")
	{
		print {$$SECSTRUCT} $tRNA->seqname().".trna".$tRNA->id()." (".$tRNA->start()."-".$tRNA->end().")\t",
			"Length: $seqlen bp\nType: ".$tRNA->isotype()."\t";
	}
	else
	{
		print {$$SECSTRUCT} $tRNA->seqname().".trna".$tRNA->id()." (".$tRNA->end()."-".$tRNA->start().")\t",
			"Length: $seqlen bp\nType: ".$tRNA->isotype()."\t";
	}

    if ($opts->output_codon())
	{
		print {$$SECSTRUCT} "Codon: ", &rev_comp_seq($tRNA->anticodon()), " at ";
    }
    else
	{
		print {$$SECSTRUCT} "Anticodon: ".$tRNA->anticodon()." at ";
    }

    if ($tRNA->anticodon() eq $gc->undef_anticodon())
	{
		print {$$SECSTRUCT} "0-0 (0-0)\t";
    }
    else
	{
		my @ar_ac_pos = $tRNA->ar_ac_pos();
		for (my $i = 0; $i < scalar(@ar_ac_pos); $i++)
		{
			if ($i > 0)
			{
				print {$$SECSTRUCT} ",";
			}
			print {$$SECSTRUCT} $ar_ac_pos[$i]->{rel_start}."-".$ar_ac_pos[$i]->{rel_end};
		}
		for (my $i = 0; $i < scalar(@ar_ac_pos); $i++)
		{
			if ($i == 0)
			{
				print {$$SECSTRUCT} " (";
			}
			elsif ($i > 0)
			{
				print {$$SECSTRUCT} ",";
			}
			
			if ($tRNA->strand() eq "+")
			{
				print {$$SECSTRUCT} $ar_ac_pos[$i]->{rel_start} + $tRNA->start() - 1, "-",
					$ar_ac_pos[$i]->{rel_start} + $tRNA->start() + 1;
			}
			else
			{
				print {$$SECSTRUCT} $tRNA->end() - $ar_ac_pos[$i]->{rel_start} + 1, "-",
					$tRNA->end() - $ar_ac_pos[$i]->{rel_start} - 1;
			}
			if ($i == (scalar(@ar_ac_pos) - 1))
			{
				print {$$SECSTRUCT} ")\t";
			}
			
		}
	}

    print {$$SECSTRUCT} "Score: ".$tRNA->score()."\n";
	my @ar_introns = ();
	my $nci_count = 0;
    if ($tRNA->get_intron_count() > 0)
	{
		@ar_introns = $tRNA->ar_introns();
		foreach my $intron (@ar_introns)
		{
			if (defined $intron)
			{
				print {$$SECSTRUCT} "Possible intron: $intron->{rel_start}-$intron->{rel_end} ";
				if ($tRNA->strand() eq "+")
				{
					print {$$SECSTRUCT} "(".$intron->{start}."-".$intron->{end}.")\n";
				}
				else
				{
					print {$$SECSTRUCT} "(".$intron->{end}."-".$intron->{start}.")\n";
				}
			}
		}
	}
	
	my $line = "";
	my $note = "";
    if ($tRNA->is_pseudo() and $tRNA->is_trunc())
	{
		$note = "Possible truncation, pseudogene";
    }
	elsif ($tRNA->is_pseudo())
	{
		$note = "Possible pseudogene";
	}
	elsif ($tRNA->is_trunc())
	{
		$note = "Possible truncation";
	}
	
	if ($note ne "")
	{
		$line = sprintf("%s", $note);
	}
	
	if ($get_hmm_score)
	{
		if ($note ne "")
		{
			$line .= ": ";
		}
		
		$line .= sprintf("HMM Sc=%.2f\tSec struct Sc=%.2f", $tRNA->hmm_score(), $tRNA->ss_score());
    }
	
	if ($opts->infernal_score())
	{
		my $inf = $tRNA->get_domain_model("infernal");
		if (defined $inf)
		{
			$line .= "\tInfernal Sc=".$inf->{score};
		}				
	}
	if ($opts->mito_mode())
	{
		$note = "";
		if ($tRNA->category() ne "")
		{
			$note = $tRNA->category();
			$note =~ s/_/ /g;
			$note =~ s/mito //g;
			$note =~ s/ac/anticodon/g;
			$note = uc(substr($note, 0, 1)).substr($note, 1);
			if ($tRNA->note() ne "")
			{
				if ($tRNA->note() =~ /^\(/)
				{
					$note = $note . " ". $tRNA->note();
				}
				else
				{
					$note = $note . ";". $tRNA->note();
				}
			}
			$line .= "Note: ".$note;
		}
		elsif ($tRNA->note() ne "")
		{
			$note = $tRNA->note();
			$line .= "Note: ".$note;
		}
	}

	if ($line ne "")
	{
		print {$$SECSTRUCT} $line."\n";
	}

	if (!$opts->arch_mode())
	{
	    print {$$SECSTRUCT} "     ",substr($ruler, 0, $seqlen - 1),"\n";
	    print {$$SECSTRUCT} "Seq: ".$tRNA->seq()."\nStr: ".$tRNA->ss()."\n\n";
	}
	else
	{
	    print {$$SECSTRUCT} "     ",substr($ruler, 0, length($tRNA->mat_seq()) - 1),"\n";
	    print {$$SECSTRUCT} "Seq: ".$tRNA->mat_seq()."\nStr: ".$tRNA->mat_ss()."\n";
		if (uc($tRNA->seq()) ne uc($tRNA->mat_seq()))
		{
			my $precursor_seq = uc($tRNA->seq());
			foreach my $intron (@ar_introns)
			{
				if (defined $intron)
				{
					my $intron_seq = uc(substr($tRNA->seq(), $intron->{rel_start} - 1, $intron->{rel_end} - $intron->{rel_start} + 1));
					if ($intron_seq ne "")
					{
						$precursor_seq =~ s/$intron_seq/\[$intron_seq\]/;
					}
				}
			}
			print {$$SECSTRUCT} "Pre: ". $precursor_seq ."\n\n";		
		}
		else
		{
			print {$$SECSTRUCT} "\n";
		}
	}
}

sub write_tRNA_sequence
{    
    my ($FASTA, $tRNA) = @_;
	
	my $tRNAscan_id = $tRNA->seqname().".trna".$tRNA->id();
	print {$FASTA} ">".$tRNAscan_id." ".
		$tRNA->seqname().":".$tRNA->start()."-".$tRNA->end()." (".$tRNA->strand().") ".
		$tRNA->isotype()." (".$tRNA->anticodon().") ".length($tRNA->seq())." bp Sc: ".$tRNA->score();
	if ($tRNA->is_pseudo())
	{
		print {$FASTA} " Possible pseudogene\n";
	}
	else
	{
		print {$FASTA} "\n";
	}
	my $parts = int(length($tRNA->seq()) / 60);
	my $remain = length($tRNA->seq()) % 60;
	for (my $j = 0; $j < $parts; $j++)
	{
		print {$FASTA} uc(substr($tRNA->seq(), $j * 60, 60))."\n";
	}
	if ($remain > 0)
	{
		print {$FASTA} uc(substr($tRNA->seq(), $parts * 60))."\n";
	}			
}

# Save tRNA hits in Tabular output
sub construct_tab_output
{
    my ($opts, $get_hmm_score, $tRNA, $max_seq_name_width, $max_seq_len_width) = @_;
    
    my ($result_line, $tRNA_type);
	my ($type, $model, $score, $ss) = $tRNA->get_highest_score_model();
	$tRNA->sort_multi_models("model");
	my ($iso_model, $iso_score, $iso_ss) = $tRNA->get_model_hit("cyto", $tRNA->isotype());
    
    $result_line = sprintf "%-".$max_seq_name_width."s\t", $tRNA->seqname();
    $result_line .= $tRNA->id()."\t";
    
	if ($tRNA->strand() eq "+")
	{
		$result_line .= sprintf "%-".$max_seq_len_width."d\t", $tRNA->start();
		$result_line .= sprintf "%-".$max_seq_len_width."d\t", $tRNA->end();
	}
	else
	{
		$result_line .= sprintf "%-".$max_seq_len_width."d\t", $tRNA->end();
		$result_line .= sprintf "%-".$max_seq_len_width."d\t", $tRNA->start();		
	}
		
    $result_line .= $tRNA->isotype()."\t";

    if ($opts->output_codon())
	{
		$result_line .= (&rev_comp_seq($tRNA->anticodon()))."\t";
    }
    else
	{
		$result_line .= $tRNA->anticodon()."\t";
    }

	my @ar_introns = ();
    if ($tRNA->get_intron_count() == 0)
	{
		$result_line .= "0\t0"; 
    }
    else
	{
		my $intron_ct = 0;
		@ar_introns = $tRNA->ar_introns();
		for (my $i = 0; $i < scalar(@ar_introns); $i++)
		{
			if (defined $ar_introns[$i])
			{
				if ($intron_ct > 0)
				{
					$result_line .= ",";
				}
				if ($tRNA->strand() eq "+")
				{
					$result_line .= $ar_introns[$i]->{start};
				}
				else
				{
					$result_line .= $ar_introns[$i]->{end};
				}
				$intron_ct++;
			}
		}
		$result_line .= "\t";
		$intron_ct = 0;
		for (my $i = 0; $i < scalar(@ar_introns); $i++)
		{
			if (defined $ar_introns[$i])
			{
				if ($intron_ct > 0)
				{
					$result_line .= ",";
				}
				if ($tRNA->strand() eq "+")
				{
					$result_line .= $ar_introns[$i]->{end};
				}
				else
				{
					$result_line .= $ar_introns[$i]->{start};
				}
				$intron_ct++;
			}
		}
	}			
    $result_line .= "\t".$tRNA->score();
 
    if ($get_hmm_score)
	{
		$result_line .= sprintf "\t%.2f\t%.2f", $tRNA->hmm_score(), $tRNA->ss_score();
    }
	if ($opts->infernal_score())
	{
		my $inf = $tRNA->get_domain_model("infernal");
		if (defined $inf)
		{
			$result_line .= "\t".$inf->{score};
		}		
	}
    if ($opts->save_source())
	{
		$result_line .= "\t".$tRNA->hit_source();
    }
	if ((($opts->euk_mode() or $opts->bact_mode() or $opts->arch_mode()) and !$opts->no_isotype() and $opts->detail()) or $opts->metagenome_mode())
	{
		$result_line .= "\t".$model;
		$result_line .= "\t".$score;
		if ($opts->euk_mode() and $opts->mito_model() ne "")
		{
			$tRNA->category($type);
			if ($type eq "cyto")
			{
				$result_line .= "\tcytosolic";
			}
			else
			{
				$result_line .= "\t".$type;
			}
		}		
	}
	if (!$opts->mito_mode() and !$opts->numt_mode())
	{
		my $note = "\t";
		if ($tRNA->is_pseudo())
		{
			$note .= "pseudo";
		}

		if ($opts->detail())
		{
			if ((($opts->euk_mode() or $opts->bact_mode() or $opts->arch_mode()) and !$opts->no_isotype()) or $opts->metagenome_mode())
			{
				if ($model ne "" and $tRNA->isotype() ne "Undet")
				{
					if ($model ne $tRNA->isotype())
					{
						if (($model eq "iMet" or $model eq "fMet" or $model eq "Ile2") and $tRNA->isotype() eq "Met")
						{}
						elsif (($model eq "LeuTAA" or $model eq "LeuTAG") and ($tRNA->isotype() eq "Leu" or $tRNA->isotype() eq "Undet") and $type eq "mito")
						{}
						elsif (($model eq "SerGCT" or $model eq "SerTGA") and ($tRNA->isotype() eq "Ser" or $tRNA->isotype() eq "Undet") and $type eq "mito")
						{}
						else
						{
							if ($note ne "\t")
							{
								$note .= ",";
							}							
							$note .= sprintf("IPD:%0.2f", ($iso_score - $score));
						}
					}		
				}
			}
		
			if ($tRNA->is_trunc())
			{
				if ($note ne "\t")
				{
					$note .= ",";
				}
				$note .= $tRNA->trunc();
			}
			
			if ($opts->search_mode() eq "archaea")
			{
				if ($tRNA->get_intron_count() > 0)
				{
					my $ci_count = 0;
					my $nci_count = 0;
					for (my $i = 0; $i < scalar(@ar_introns); $i++)
					{
						if ($ar_introns[$i]->{type} eq "CI")
						{
							$ci_count++;
						}
						elsif ($ar_introns[$i]->{type} eq "NCI")
						{
							$nci_count++;
						}
					}
					if ($ci_count > 0)
					{
						if ($note ne "\t")
						{
							$note .= ",";
						}
						$note .= "CI";
					}
					if ($nci_count > 0)
					{
						if ($note ne "\t")
						{
							$note .= ",";
						}
						$note .= "NCI:".$nci_count;
					}
				}
			}
		}

		$result_line .= $note;
	}
	if ($opts->mito_mode())
	{
		my $note = "";
		if ($tRNA->category() ne "" and $opts->detail())
		{
			$note = $tRNA->category();
			$note =~ s/mito_//g;
			if ($tRNA->note() ne "")
			{
				if ($tRNA->note() =~ /^\(/)
				{
					$note = $note . " " . $tRNA->note();
				}
				else
				{
					$note = $note . ";" . $tRNA->note();
				}
			}
		}
		elsif ($tRNA->note() ne "" and $opts->detail())
		{
			$note = $tRNA->note();
		}
		$result_line .= "\t".$note;
	}
    $result_line .= "\n";
    
    return $result_line;
}

sub get_models
{
	my ($opts, $cm) = @_;
	my %iso_models = ();
	my %mito_models = ();
	
	my $models = $cm->get_models_from_db($cm->isotype_cm_db_file_path());
    foreach my $cur_iso_cm (sort @$models)
    {
		my $model = $cur_iso_cm;
		if ($cur_iso_cm =~ /^arch-(\S+)/ || $cur_iso_cm =~ /^euk-(\S+)/ || $cur_iso_cm =~ /^bact-(\S+)/)
		{
			$model = $1;
		}
		$iso_models{$model} = 1;
	}

	if ($opts->euk_mode() and $opts->mito_model() ne "")
	{
		$models = $cm->get_models_from_db($cm->mito_isotype_cm_db_file_path());
		foreach my $cur_iso_cm (sort @$models)
		{
			$mito_models{$cur_iso_cm} = 1;
		}
	}	
	
	return (\%iso_models, \%mito_models);
}

sub print_isotype_specific_header
{    
    my ($ISOTYPE, $opts, $iso_models, $mito_models) = @_;

	print {$$ISOTYPE} "tRNAscanID\tAnticodon_predicted_isotype";
    foreach my $cur_iso_cm (sort keys %$iso_models)
    {
		print {$$ISOTYPE} "\t".$cur_iso_cm;
	}

	if ($opts->euk_mode() and $opts->mito_model() ne "")
	{
		foreach my $cur_iso_cm (sort keys %$mito_models)
		{
			print {$$ISOTYPE} "\tmito_".$cur_iso_cm;
		}
	}	

	print {$$ISOTYPE} "\n";		
}

sub construct_isotype_specific_output
{    
    my ($opts, $iso_models, $mito_models, $tRNA) = @_;
	
	my $result_line = $tRNA->seqname().".trna".$tRNA->id();
	$result_line .= "\t".$tRNA->isotype();
	$tRNA->sort_multi_models("model");
	
    foreach my $cur_iso_cm (sort keys %$iso_models)
    {
		my ($model, $score, $ss) = $tRNA->get_model_hit("cyto", $cur_iso_cm);
		if ($model ne "")
		{
			$result_line .= "\t".$score;
		}
		else
		{
			$result_line .= "\t-999";
		}
	}
	
	if ($opts->euk_mode() and $opts->mito_model() ne "")
	{
		foreach my $cur_iso_cm (sort keys %$mito_models)
		{
			my ($model, $score, $ss) = $tRNA->get_model_hit("mito", $cur_iso_cm);
			if ($model ne "")
			{
				$result_line .= "\t".$score;
			}
			else
			{
				$result_line .= "\t-999";
			}
		}
	}	
	
	$result_line .= "\n";
	
	return $result_line;
}

sub save_Acedb_from_secpass
{
    my ($ACEOUT, $opts, $gc, $tRNA, $program_id) = @_;                     

    print {$$ACEOUT} "Sequence\t".$tRNA->seqname()."\nSubsequence\t".$tRNA->tRNAscan_id()." ".$tRNA->start()." ".$tRNA->end()."\n\n";
    print {$$ACEOUT} "Sequence\t".$tRNA->tRNAscan_id()."\nSource\t\t".$tRNA->seqname()."\n";
    if ($tRNA->get_intron_count() > 0)
	{
		my @ar_introns = $tRNA->ar_introns();
		print {$$ACEOUT} "Source_Exons\t1 ", $ar_introns[0]->{rel_start} - 1,"\n";
		print {$$ACEOUT} "Source_Exons\t", $ar_introns[0]->{rel_end} + 1," ", abs($tRNA->end() - $tRNA->start()) + 1,"\n";
    }	   
    print {$$ACEOUT} "Brief_identification tRNA-".$tRNA->isotype()."\n",
        "Transcript tRNA \"";

    if ($opts->output_codon())
	{
		print {$$ACEOUT} &rev_comp_seq($tRNA->anticodon());
    }
    else
	{
		print {$$ACEOUT} $tRNA->anticodon();
    }
    
    print {$$ACEOUT} " ".$tRNA->isotype()." ", $gc->one_let_trans_map()->{$tRNA->isotype()},
        "\"\nScore $program_id ".$tRNA->score()."\n";

    if ($tRNA->is_pseudo())
	{
		printf {$$ACEOUT} "Remark \"Likely pseudogene (HMM Sc=%.2f / Sec struct Sc=%.2f)\"\n",
            $tRNA->hmm_score(), $tRNA->ss_score();
    }
    print {$$ACEOUT} "\n";
}

sub write_bed
{
	my ($global_vars) = @_;
	my $opts = $global_vars->{options};
    my $sp_int_results = $global_vars->{sp_int_results};
    my $iso_int_results = $global_vars->{iso_int_results};
	
	$sp_int_results->sort_records("bed_output");
	if (!$opts->no_isotype())
	{
		$iso_int_results->sort_records("tRNAscan_id");
	}
	
	my $tRNA = tRNAscanSE::tRNA->new;
    &open_for_append(\*FILE_OUT, $opts->bed_file());
    my @sp_indexes = $sp_int_results->get_indexes();
    my @ iso_indexes = $iso_int_results->get_indexes();
    if ($sp_int_results->open_file("read"))
    {		
		if (!$opts->no_isotype())
		{
			$iso_int_results->open_file("read");
		}
		
		for (my $i = 0; $i < scalar(@sp_indexes); $i++)
		{
			$sp_int_results->get_tRNA($sp_indexes[$i]->[0], $tRNA);
			
			if (!$opts->no_isotype())
			{
				my $id = $tRNA->seqname().".t".&pad_num($tRNA->id(), 6);
				my $index = $iso_int_results->bsearch_tRNAscan_id($id);
				if ($index > -1)
				{
					$iso_int_results->get_tRNA($iso_indexes[$index]->[0], $tRNA);
					my ($type, $model, $score, $ss) = $tRNA->get_highest_score_model();
					if ($tRNA->isotype() eq "Met" and $type eq "cyto" and ($model eq "iMet" or $model eq "fMet" or $model eq "Ile2"))
					{
						$tRNA->isotype($model);
						$tRNA->tRNAscan_id($tRNA->seqname().".tRNA".$tRNA->id()."-".$tRNA->isotype().$tRNA->anticodon());
					}
					elsif ($tRNA->isotype() eq "Met" and $type eq "cyto" and $model ne "Met" and $model ne "iMet" and $model ne "fMet")
					{
						$tRNA->sort_multi_models("model");
						my ($met_iso_model, $met_iso_score, $met_iso_ss) = $tRNA->get_model_hit("cyto", $tRNA->isotype());
						my ($ile2_iso_model, $ile2_iso_score, $ile2_iso_ss) = $tRNA->get_model_hit("cyto", "Ile2");
						if ($ile2_iso_score > 0 and $met_iso_score > 0)
						{
							if (($score - $ile2_iso_score) <= 5 and ($ile2_iso_score - $met_iso_score) >= 5 and $tRNA->score() > 50)
							{
								$tRNA->isotype("Ile2");
								$tRNA->tRNAscan_id($tRNA->seqname().".tRNA".$tRNA->id()."-".$tRNA->isotype().$tRNA->anticodon());
							}
						}
					}
				}				
			}			
		
			print FILE_OUT $tRNA->seqname()."\t".($tRNA->start() - 1)."\t".$tRNA->end()."\t".$tRNA->tRNAscan_id()."\t".&convert_bed_score($tRNA->score())."\t".
				$tRNA->strand()."\t".($tRNA->start() - 1)."\t".$tRNA->end()."\t0\t".($tRNA->get_intron_count() + 1)."\t";
			if ($tRNA->get_intron_count() == 0)
			{
				print FILE_OUT ($tRNA->end() - $tRNA->start() + 1).",\t0,\n";
			}
			else
			{
				my @ar_introns = $tRNA->ar_introns();
				my $block_sizes = "";
				my $block_starts = "0,";
				my $prev_start = 1;
				if ($tRNA->strand() eq "+")
				{
					for (my $i = 0; $i < scalar(@ar_introns); $i++)
					{
						$block_sizes .= ($ar_introns[$i]->{rel_start} - $prev_start).",";
						$block_starts .= $ar_introns[$i]->{rel_end}.",";
						$prev_start = $ar_introns[$i]->{rel_end} + 1;
					}
					$block_sizes .= ($tRNA->end() - $ar_introns[(scalar(@ar_introns)-1)]->{end}).",";
				}
				else
				{
					$prev_start = length($tRNA->seq());
					for (my $i = (scalar(@ar_introns)-1); $i >= 0; $i--)
					{
						$block_sizes .= ($prev_start - $ar_introns[$i]->{rel_end}).",";
						$block_starts .= ($prev_start - $ar_introns[$i]->{rel_start} + 1).",";
						$prev_start = $ar_introns[$i]->{rel_start};
					}
					$block_sizes .= ($ar_introns[0]->{rel_start} - 1).",";
				}
				print FILE_OUT $block_sizes."\t".$block_starts."\n";
			}			
		}
		
		if (!$opts->no_isotype())
		{
			$iso_int_results->close_file();
		}
		$sp_int_results->close_file();
	}
	close(FILE_OUT);
}

sub convert_bed_score
{
	my ($cm_score) = @_;
	
	my $bed_score = $cm_score * 10;
	if ($bed_score > 1000)
	{
		$bed_score = 1000;
	}
	elsif ($bed_score < 0)
	{
		$bed_score = 0;
	}
	
	return $bed_score;
}

sub write_gff
{
	my ($global_vars) = @_;
	my $opts = $global_vars->{options};
    my $sp_int_results = $global_vars->{sp_int_results};
    my $iso_int_results = $global_vars->{iso_int_results};
	
	$sp_int_results->sort_records("bed_output");
	if (!$opts->no_isotype())
	{
		$iso_int_results->sort_records("tRNAscan_id");
	}
	
	my $tRNA = tRNAscanSE::tRNA->new;
    &open_for_append(\*FILE_OUT, $opts->gff_file());
	print FILE_OUT "##gff-version 3\n";
    my @sp_indexes = $sp_int_results->get_indexes();
    my @ iso_indexes = $iso_int_results->get_indexes();
    if ($sp_int_results->open_file("read"))
    {		
		if (!$opts->no_isotype())
		{
			$iso_int_results->open_file("read");
		}
		
		for (my $i = 0; $i < scalar(@sp_indexes); $i++)
		{
			$sp_int_results->get_tRNA($sp_indexes[$i]->[0], $tRNA);
			
			if (!$opts->no_isotype())
			{
				my $id = $tRNA->seqname().".t".&pad_num($tRNA->id(), 6);
				my $index = $iso_int_results->bsearch_tRNAscan_id($id);
				if ($index > -1)
				{
					$iso_int_results->get_tRNA($iso_indexes[$index]->[0], $tRNA);
					my ($type, $model, $score, $ss) = $tRNA->get_highest_score_model();
					if ($tRNA->isotype() eq "Met" and $type eq "cyto" and ($model eq "iMet" or $model eq "fMet" or $model eq "Ile2"))
					{
						$tRNA->isotype($model);
						$tRNA->tRNAscan_id($tRNA->seqname().".tRNA".$tRNA->id()."-".$tRNA->isotype().$tRNA->anticodon());
					}
					elsif ($tRNA->isotype() eq "Met" and $type eq "cyto" and $model ne "Met" and $model ne "iMet" and $model ne "fMet")
					{
						$tRNA->sort_multi_models("model");
						my ($met_iso_model, $met_iso_score, $met_iso_ss) = $tRNA->get_model_hit("cyto", $tRNA->isotype());
						my ($ile2_iso_model, $ile2_iso_score, $ile2_iso_ss) = $tRNA->get_model_hit("cyto", "Ile2");
						if ($ile2_iso_score > 0 and $met_iso_score > 0)
						{
							if (($score - $ile2_iso_score) <= 5 and ($ile2_iso_score - $met_iso_score) >= 5 and $tRNA->score() > 50)
							{
								$tRNA->isotype("Ile2");
								$tRNA->tRNAscan_id($tRNA->seqname().".tRNA".$tRNA->id()."-".$tRNA->isotype().$tRNA->anticodon());
							}
						}
					}
				}				
			}
			my $biotype = "tRNA";			
			if ($tRNA->is_pseudo())
			{
				$biotype = "pseudogene";
			}
			print FILE_OUT $tRNA->seqname()."\ttRNAscan-SE\t".$biotype."\t".$tRNA->start()."\t".$tRNA->end()."\t".$tRNA->score()."\t".$tRNA->strand()."\t.\t".
				"ID=".$tRNA->seqname().".trna".$tRNA->id().";Name=".$tRNA->tRNAscan_id().";isotype=".$tRNA->isotype().";anticodon=".$tRNA->anticodon().
				";gene_biotype=".$biotype.";\n";

			if ($tRNA->get_intron_count() == 0)
			{
				print FILE_OUT $tRNA->seqname()."\ttRNAscan-SE\texon\t".$tRNA->start()."\t".$tRNA->end()."\t.\t".$tRNA->strand()."\t.\t".
					"ID=".$tRNA->seqname().".trna".$tRNA->id().".exon1;Parent=".$tRNA->seqname().".trna".$tRNA->id().";\n";
			}
			else
			{
				my @ar_introns = $tRNA->ar_introns();
				if ($tRNA->strand() eq "+")
				{
					print FILE_OUT $tRNA->seqname()."\ttRNAscan-SE\texon\t".$tRNA->start()."\t";
					for (my $i = 0; $i < scalar(@ar_introns); $i++)
					{
						print FILE_OUT ($ar_introns[$i]->{start}-1)."\t.\t".$tRNA->strand()."\t.\t".
							"ID=".$tRNA->seqname().".trna".$tRNA->id().".exon".($i+1).";Parent=".$tRNA->seqname().".trna".$tRNA->id().";\n";
						print FILE_OUT $tRNA->seqname()."\ttRNAscan-SE\texon\t".($ar_introns[$i]->{end}+1)."\t";
					}
					print FILE_OUT $tRNA->end()."\t.\t".$tRNA->strand()."\t.\t".
						"ID=".$tRNA->seqname().".trna".$tRNA->id().".exon".(scalar(@ar_introns)+1).";Parent=".$tRNA->seqname().".trna".$tRNA->id().";\n";
				}
				else
				{
					my $end = $tRNA->end();
					for (my $i = 0; $i < scalar(@ar_introns); $i++)
					{
						print FILE_OUT $tRNA->seqname()."\ttRNAscan-SE\texon\t".($ar_introns[$i]->{end}+1)."\t".$end."\t.\t".$tRNA->strand()."\t.\t".
							"ID=".$tRNA->seqname().".trna".$tRNA->id().".exon".($i+1).";Parent=".$tRNA->seqname().".trna".$tRNA->id().";\n";
						$end = $ar_introns[$i]->{start} - 1;
					}
					print FILE_OUT $tRNA->seqname()."\ttRNAscan-SE\texon\t".$tRNA->start()."\t".$end."\t.\t".$tRNA->strand()."\t.\t".
						"ID=".$tRNA->seqname().".trna".$tRNA->id().".exon".(scalar(@ar_introns)+1).";Parent=".$tRNA->seqname().".trna".$tRNA->id().";\n";
				}
			}			
		}
		
		if (!$opts->no_isotype())
		{
			$iso_int_results->close_file();
		}
		$sp_int_results->close_file();
	}
	close(FILE_OUT);
}

sub output_split_fragments
{
    my ($opts, $r_pairs, $r_5half_hits, $r_3half_hits) = @_;                     
	
	my ($r_5half, $r_3half);
	
	&open_for_append(\*SPLITFILE, $opts->split_fragment_file());
	printf SPLITFILE "Fragment1\tFragment2\tSeqName1\tStartPos1\tEndPos1\tSeqName2\tStartPos2\tEndPos2\tScore1\tScore2\n";
	
	foreach my $r_pair (@$r_pairs)
	{
		if (defined $r_pair->{"5h"} && defined $r_pair->{"3h"})
		{
			$r_5half = $r_5half_hits->[$r_pair->{"5h"}];
			$r_3half = $r_3half_hits->[$r_pair->{"3h"}];
			print SPLITFILE $r_5half->{seq}."\t".$r_3half->{seq}."\t",
				$r_5half->{hit_seqname}."\t".$r_5half->{tRNA_start}."\t".$r_5half->{tRNA_end}."\t",
				$r_3half->{hit_seqname}."\t".$r_3half->{tRNA_start}."\t".$r_3half->{tRNA_end}."\t",
				$r_5half->{score}."\t".$r_3half->{score}."\n";
			print SPLITFILE $r_5half->{ss}."\t".$r_3half->{ss}."\t\t\t\t\t\t\t\t\n";
		}
		elsif (defined $r_pair->{"5h"} && !defined $r_pair->{"3h"})
		{
			$r_5half = $r_5half_hits->[$r_pair->{"5h"}];
			print SPLITFILE $r_5half->{seq}."\t\t",
				$r_5half->{hit_seqname}."\t".$r_5half->{tRNA_start}."\t".$r_5half->{tRNA_end}."\t",
				"\t\t\t",
				$r_5half->{score}."\t\n";
			print SPLITFILE $r_5half->{ss}."\t\t\t\t\t\t\t\t\t\n";
		}
		elsif (!defined $r_pair->{"5h"} && defined $r_pair->{"3h"})
		{
			$r_3half = $r_3half_hits->[$r_pair->{"3h"}];
			print SPLITFILE "\t".$r_3half->{seq}."\t",
				"\t\t\t",
				$r_3half->{hit_seqname}."\t".$r_3half->{tRNA_start}."\t".$r_3half->{tRNA_end}."\t",
				"\t".$r_3half->{score}."\n";
			print SPLITFILE "\t".$r_3half->{ss}."\t\t\t\t\t\t\t\t\n";
		}
	}
	close(SPLITFILE);
}

1;
