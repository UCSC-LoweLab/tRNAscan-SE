# tRNAscanSE/Options.pm
# This class defines options used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2017 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::Options;

use strict;
use tRNAscanSE::Utils;
use tRNAscanSE::CM;

sub new {
    my $class = shift;
    my $self = {};

    initialize($self);

    bless ($self, $class);
    return $self;
}

sub DESTROY
{
    my $self = shift;
}

sub initialize
{
    my $self = shift;
    $self->{fafile} = "";
    $self->{fasta_file} = "";           # input sequence file
    $self->{multiple_files} = 0;        # multiple input sequence files
    $self->{out_file} = "-";            # output result file -- send to 
                                        #  stdout ("-") by default 

    $self->{results_to_stdout} = 1;     # send results to stdout by default

    $self->{ace_output} = 0;            # output in ACeDB format if non-zero
    $self->{brief_output} = 0;          # don't print tabular output column headers
                                        #  if non-zero
    $self->{quiet_mode} = 0;            # don't print credits & selected run options
                                        #  if non-zero
    $self->{display_progress} = 0;      # print program progress info if non-zero
    $self->{save_progress} = 0;         # save progress to log file if non-zero
    $self->{log_file} = "";             # name of log file

    $self->{seq_key} = "";              # require seq names to match this key
    $self->{raw_seq_key} = "";          # unmodified user-input key
    $self->{start_at_key} = 0;          # read all seqs after finding seqname=KEY?

    $self->{tscan_mode} = 1;            # run tRNAscan if non-zero
    $self->{eufind_mode} = 1;           # run eufindtRNA (pavesi) if non-zero
    $self->{strict_params} = 1;         # use original strict tRNAscan params if non-zero
    $self->{infernal_fp} = 0;           # run infernal hmm as first pass
                                    
    $self->{CM_mode} = "cove";          # run covariance model search in second pass
                                        #  cove - run Cove if non-zero
                                        #  infernal - run Infernal if non-zero
                                        #             run Infernal by default
    
    $self->{hmm_filter} = 0;            # run infernal with hmm filter in second pass
    
    $self->{second_pass_label} = "Cove";    # Second scan pass label: Cove by default

    $self->{search_mode} = "";          # tRNA search mode
                                        # euk - run eukaryotic cov model if set
                                        # bacteria - run covariance model for bacteria if set
                                        # archaea - run archaea cov model if set
                                        # organelle - run general cov model 
                                        # general - run general cov model (combines tRNAs from all 3 domains)
                                        # mito mammal - run mammalian mito cov models
                                        # mito vert - run vertebrate mito cov models
                                        # metagenome
                                        # numt
                                        
    $self->{no_isotype} = 0;            # flag for disabling isotype specific cm scan
    
    $self->{mito_model} = "";           # mito CM for isotype specific scan
    
    $self->{org_mode} = 0;              # run in organellar mode
                                        # run eukaryotic model by default

    $self->{alt_gcode} = 0;             # use alternate genetic translation table
                                        #  file if non-zero
    $self->{gc_file} = "";              # alternate transl table file

    $self->{save_stats} = 0;            # save statistics for search
    $self->{stats_file} = "";

    $self->{detail} = 0;                 # display detailed prediction info
    
    $self->{save_odd_struct} = 0;        # save structures for which Cove
                                         #  was unable to determine anticodon
    $self->{odd_struct_file} = "";

    $self->{save_all_struct} = 0;        # save secondary structures if nonzero
    
    $self->{bed_file} = "";              # bed file
    $self->{gff_file} = "";              # gff file
    $self->{output_fasta_file} = "";     # predicted tRNA fasta file
    $self->{all_struct_file} = "";       # sec struct file, set with -f option
    $self->{isotype_specific_file} = ""; # isotype specific result file

    $self->{split_fragment_file} = "";   # split fragment file, set with --frag option

    $self->{save_verbose} = 0;           # save verbose output from tRNAscan
    $self->{verb_file} = "";

    $self->{save_firstpass_res} = 0;   # save tabular tRNAscan results
    $self->{firstpass_result_file} = "";
    $self->{firstpass_flanking_file} = "";
    $self->{secondpass_int_result_file} = "";
    $self->{isotype_int_result_file} = "";
    $self->{truncated_int_result_file} = "";

    $self->{use_prev_ts_run} = 0;   # specify result file from previous
                                    # tRNA search for Cove-confirmation

    $self->{default_Padding} = 10;
    $self->{padding} = $self->{default_Padding}; # pad both ends of first-pass hits with this
                                                 # many extra bases before passing to Cove

    $self->{save_falsepos} = 0;     # save false positive tRNAs in 
                                    # fasta file
    $self->{falsepos_file} = "";

    $self->{save_missed} = 0;       # save seqs without a hit
    $self->{missed_seq_file} = "";

    $self->{save_source} = 0;       # save source of first-pass hit

    $self->{output_codon} = 0;      # output tRNA codon instead of anticodon
                                    # (off by default)
    
    $self->{infernal_score} = 0;    # flag to output infernal score
 
#    $self->{use_orig_cm} = 0;       # use original covariance model that
#                                    # contains tRNAS from all three domains

    $self->{def_max_int_len} = 200;     # default MAX intron+variable loop region size
                                        # used in EufindtRNA

    $self->{max_int_len} = $self->{def_max_int_len};

    $self->{prompt_for_overwrite} = 1;  # prompt user before overwriting a pre-existing 
                                        # output file, disabled with -Q option
}

sub fafile
{
    my $self = shift;
    if (@_) { $self->{fafile} = shift; }
    return $self->{fafile};
}

sub fasta_file
{
    my $self = shift;
    if (@_) { $self->{fasta_file} = shift; }
    return $self->{fasta_file};
}

sub multiple_files
{
    my $self = shift;
    if (@_) { $self->{multiple_files} = shift; }
    return $self->{multiple_files};
}

sub out_file
{
    my $self = shift;
    if (@_) { $self->{out_file} = shift; }
    return $self->{out_file};
}

sub results_to_stdout
{
    my $self = shift;
    if (@_) { $self->{results_to_stdout} = shift; }
    return $self->{results_to_stdout};
}

sub ace_output
{
    my $self = shift;
    if (@_) { $self->{ace_output} = shift; }
    return $self->{ace_output};
}

sub brief_output
{
    my $self = shift;
    if (@_) { $self->{brief_output} = shift; }
    return $self->{brief_output};
}

sub quiet_mode
{
    my $self = shift;
    if (@_) { $self->{quiet_mode} = shift; }
    return $self->{quiet_mode};
}

sub display_progress
{
    my $self = shift;
    if (@_) { $self->{display_progress} = shift; }
    return $self->{display_progress};
}

sub save_progress
{
    my $self = shift;
    if (@_) { $self->{save_progress} = shift; }
    return $self->{save_progress};
}

sub log_file
{
    my $self = shift;
    if (@_) { $self->{log_file} = shift; }
    return $self->{log_file};
}

sub seq_key
{
    my $self = shift;
    if (@_) { $self->{seq_key} = shift; }
    return $self->{seq_key};
}

sub raw_seq_key
{
    my $self = shift;
    if (@_) { $self->{raw_seq_key} = shift; }
    return $self->{raw_seq_key};
}

sub start_at_key
{
    my $self = shift;
    if (@_) { $self->{start_at_key} = shift; }
    return $self->{start_at_key};
}

sub tscan_mode
{
    my $self = shift;
    if (@_) { $self->{tscan_mode} = shift; }
    return $self->{tscan_mode};
}

sub eufind_mode
{
    my $self = shift;
    if (@_) { $self->{eufind_mode} = shift; }
    return $self->{eufind_mode};
}

sub strict_params
{
    my $self = shift;
    if (@_) { $self->{strict_params} = shift; }
    return $self->{strict_params};
}

sub infernal_fp
{
    my $self = shift;
    if (@_) { $self->{infernal_fp} = shift; }
    return $self->{infernal_fp};
}

sub CM_mode
{
    my $self = shift;
    if (@_) { $self->{CM_mode} = shift; }
    return $self->{CM_mode};
}

sub cove_mode
{
    my $self = shift;
    return ($self->{CM_mode} eq 'cove');
}

sub infernal_mode
{
    my $self = shift;
    return ($self->{CM_mode} eq 'infernal');
}

sub hmm_filter
{
    my $self = shift;
    if (@_) { $self->{hmm_filter} = shift; }
    return $self->{hmm_filter};    
}

sub second_pass_label
{
    my $self = shift;
    if (@_) { $self->{second_pass_label} = shift; }
    return $self->{second_pass_label};
}

sub search_mode
{
    my $self = shift;
    if (@_) { $self->{search_mode} = shift; }
    return $self->{search_mode};
}

sub euk_mode
{
    my $self = shift;
    return ($self->{search_mode} eq 'euk');    
}

sub bact_mode
{
    my $self = shift;
    return ($self->{search_mode} eq 'bacteria');
}

sub arch_mode
{
    my $self = shift;
    return ($self->{search_mode} eq 'archaea');
}

sub mito_mammal_mode
{
    my $self = shift;
    return ($self->mito_mode() and $self->{mito_model} eq 'mammal');
}

sub mito_vert_mode
{
    my $self = shift;
    return ($self->mito_mode() and $self->{mito_model} eq 'vert');
}

sub mito_mode
{
    my $self = shift;
    return ($self->{search_mode} eq 'mito');
}

sub metagenome_mode
{
    my $self = shift;
    return ($self->{search_mode} eq 'metagenome');
}

sub numt_mode
{
    my $self = shift;
    return ($self->{search_mode} eq 'numt');
}

sub alternate_mode
{
    my $self = shift;
    return ($self->{search_mode} eq 'alt');
}

sub no_isotype
{
    my $self = shift;
    if (@_) { $self->{no_isotype} = shift; }
    return $self->{no_isotype};
}

sub mito_model
{
    my $self = shift;
    if (@_) { $self->{mito_model} = shift; }
    return $self->{mito_model};
}

sub organelle_mode
{
    my $self = shift;
    return ($self->{search_mode} eq 'organelle');
}

sub general_mode
{
    my $self = shift;
    return ($self->{search_mode} eq 'general');
}

sub org_mode
{
    my $self = shift;
    if (@_) { $self->{org_mode} = shift; }
    return $self->{org_mode};
}

sub alt_gcode
{
    my $self = shift;
    if (@_) { $self->{alt_gcode} = shift; }
    return $self->{alt_gcode};
}

sub gc_file
{
    my $self = shift;
    if (@_) { $self->{gc_file} = shift; }
    return $self->{gc_file};
}

sub save_stats
{
    my $self = shift;
    if (@_) { $self->{save_stats} = shift; }
    return $self->{save_stats};
}

sub stats_file
{
    my $self = shift;
    if (@_) { $self->{stats_file} = shift; }
    return $self->{stats_file};
}

sub save_odd_struct
{
    my $self = shift;
    if (@_) { $self->{save_odd_struct} = shift; }
    return $self->{save_odd_struct};
}

sub odd_struct_file
{
    my $self = shift;
    if (@_) { $self->{odd_struct_file} = shift; }
    return $self->{odd_struct_file};
}

sub save_all_struct
{
    my $self = shift;
    if (@_) { $self->{save_all_struct} = shift; }
    return $self->{save_all_struct};
}

sub all_struct_file
{
    my $self = shift;
    if (@_) { $self->{all_struct_file} = shift; }
    return $self->{all_struct_file};
}

sub bed_file
{
    my $self = shift;
    if (@_) { $self->{bed_file} = shift; }
    return $self->{bed_file};
}

sub gff_file
{
    my $self = shift;
    if (@_) { $self->{gff_file} = shift; }
    return $self->{gff_file};
}

sub output_fasta_file
{
    my $self = shift;
    if (@_) { $self->{output_fasta_file} = shift; }
    return $self->{output_fasta_file};
}

sub isotype_specific_file
{
    my $self = shift;
    if (@_) { $self->{isotype_specific_file} = shift; }
    return $self->{isotype_specific_file};
}

sub split_fragment_file
{
    my $self = shift;
    if (@_) { $self->{split_fragment_file} = shift; }
    return $self->{split_fragment_file};
}

sub save_verbose
{
    my $self = shift;
    if (@_) { $self->{save_verbose} = shift; }
    return $self->{save_verbose};
}

sub verb_file
{
    my $self = shift;
    if (@_) { $self->{verb_file} = shift; }
    return $self->{verb_file};
}

sub save_firstpass_res
{
    my $self = shift;
    if (@_) { $self->{save_firstpass_res} = shift; }
    return $self->{save_firstpass_res};
}

sub firstpass_result_file
{
    my $self = shift;
    if (@_) { $self->{firstpass_result_file} = shift; }
    return $self->{firstpass_result_file};
}

sub firstpass_flanking_file
{
    my $self = shift;
    if (@_) { $self->{firstpass_flanking_file} = shift; }
    return $self->{firstpass_flanking_file};
}

sub secondpass_int_result_file
{
    my $self = shift;
    if (@_) { $self->{secondpass_int_result_file} = shift; }
    return $self->{secondpass_int_result_file};
}

sub isotype_int_result_file
{
    my $self = shift;
    if (@_) { $self->{isotype_int_result_file} = shift; }
    return $self->{isotype_int_result_file};
}

sub truncated_int_result_file
{
    my $self = shift;
    if (@_) { $self->{truncated_int_result_file} = shift; }
    return $self->{truncated_int_result_file};
}

sub use_prev_ts_run
{
    my $self = shift;
    if (@_) { $self->{use_prev_ts_run} = shift; }
    return $self->{use_prev_ts_run};
}

sub default_Padding
{
    my $self = shift;
    if (@_) { $self->{default_Padding} = shift; }
    return $self->{default_Padding};
}

sub padding
{
    my $self = shift;
    if (@_) { $self->{padding} = shift; }
    return $self->{padding};
}

sub save_falsepos
{
    my $self = shift;
    if (@_) { $self->{save_falsepos} = shift; }
    return $self->{save_falsepos};
}

sub falsepos_file
{
    my $self = shift;
    if (@_) { $self->{falsepos_file} = shift; }
    return $self->{falsepos_file};
}

sub save_missed
{
    my $self = shift;
    if (@_) { $self->{save_missed} = shift; }
    return $self->{save_missed};
}

sub missed_seq_file
{
    my $self = shift;
    if (@_) { $self->{missed_seq_file} = shift; }
    return $self->{missed_seq_file};
}

sub save_source
{
    my $self = shift;
    if (@_) { $self->{save_source} = shift; }
    return $self->{save_source};
}

sub output_codon
{
    my $self = shift;
    if (@_) { $self->{output_codon} = shift; }
    return $self->{output_codon};
}

sub detail
{
    my $self = shift;
    if (@_) { $self->{detail} = shift; }
    return $self->{detail};
}

sub infernal_score
{
    my $self = shift;
    if (@_) { $self->{infernal_score} = shift; }
    return $self->{infernal_score};
}

sub def_max_int_len
{
    my $self = shift;
    if (@_) { $self->{def_max_int_len} = shift; }
    return $self->{def_max_int_len};
}

sub max_int_len
{
    my $self = shift;
    if (@_) { $self->{max_int_len} = shift; }
    return $self->{max_int_len};
}

sub prompt_for_overwrite
{
    my $self = shift;
    if (@_) { $self->{prompt_for_overwrite} = shift; }
    return $self->{prompt_for_overwrite};
}

sub display_run_options
{
    my $self = shift;
    my ($cm, $tscan, $eufind, $global_constants, $FHAND) = @_;

    print $FHAND ('-' x 60,"\n", "Sequence file(s) to search:        ",join(', ',@ARGV),"\n");
    if ($self->{seq_key} ne '\S*')
    {
        if ($self->{start_at_key})
        {
            print $FHAND "Starting at sequence name:         $self->{raw_seq_key}\n";
        }
        else
        {
            print $FHAND "Search only names matching:        $self->{raw_seq_key}\n";
        }
    }

    print $FHAND "Search Mode:                       ";
    if ($self->euk_mode())
    {
        print $FHAND "Eukaryotic";
    }
    elsif ($self->bact_mode())
    {
        print $FHAND "Bacterial";
    }
    elsif ($self->arch_mode())
    {
        print $FHAND "Archaeal";
    }	
    elsif ($self->mito_mammal_mode())
    {
        print $FHAND "Mitochondrial in mammals";
    }	
    elsif ($self->mito_vert_mode())
    {
        print $FHAND "Mitochondrial in vertebrates";
    }	
    elsif ($self->organelle_mode())
    {
        print $FHAND "Organellar";
    }	
    elsif ($self->general_mode())
    {
        print $FHAND "General";
    }	
    elsif ($self->metagenome_mode())
    {
        print $FHAND "Metagenome";
    }	
    elsif ($self->numt_mode())
    {
        print $FHAND "Numt";
    }	
    else
    {
        print $FHAND "Eukaryotic";
    }	
    print $FHAND "\n";

    print $FHAND "Results written to:                ", &print_filename($self->{out_file}),"\n";

    print $FHAND "Output format:                     ";
    if ($self->{ace_output})
    {
        print $FHAND "ACeDB\n";
    }
    else
    {
        print $FHAND "Tabular\n";
    }

    print $FHAND "Searching with:                    ";
    if ($self->{eufind_mode})
    {
        if ($self->{tscan_mode})
        {
            if ($self->{CM_mode} =~ /infernal|cove/)
            {
                print $FHAND "tRNAscan + EufindtRNA -> $self->{second_pass_label}\n";
            }
            else
            {
                print $FHAND "tRNAscan + EufindtRNA\n";
            }
        }
        elsif ($self->{CM_mode} =~ /infernal|cove/)
        {
            print $FHAND "EufindtRNA->$self->{second_pass_label}\n";
        }
        else
        {
            print $FHAND "EufindtRNA only\n";
        }
    }
    elsif ($self->{tscan_mode})
    {
        if ($self->{CM_mode} =~ /infernal|cove/)
        {
            print $FHAND "tRNAscan->$self->{second_pass_label}\n";
        }
        else
        {
            print $FHAND "tRNAscan only\n";
        }
    }
    elsif ($self->{infernal_fp})
    {
        if ($self->{CM_mode} =~ /infernal|cove/)
        {
            print $FHAND "Infernal First Pass->$self->{second_pass_label}\n";
        }
        else
        {
            print $FHAND "Infernal First Pass only\n";
        }        
    }
    else
    {
        print $FHAND "$self->{second_pass_label} single-pass scan\n";
    }

    if ($self->{CM_mode} eq "infernal" and !$self->{infernal_fp})
    {
        if ($self->{hmm_filter})
        {
            print $FHAND "                                   Fast mode\n";
        }
        else
        {
            print $FHAND "                                   Maximum sensitivity mode\n";
        }
    }
    
    if (!$self->{no_isotype} and ($self->euk_mode() or $self->bact_mode() or $self->arch_mode()))
    {
        print $FHAND "Isotype-specific model scan:       Yes\n";        
    }
    elsif ($self->{no_isotype} and ($self->euk_mode() or $self->bact_mode() or $self->arch_mode()))
    {
        print $FHAND "Isotype-specific model scan:       No\n";      
    }
    
    if ($self->euk_mode() and $self->mito_model() ne "")
    {
        print $FHAND "Mito isotype model scan:           ".$self->mito_model()."\n";        
    }
    
    if ($cm->CM_check_for_introns())
    {
        print $FHAND "Scan for noncanonical introns\n";        
    }
    if ($cm->CM_check_for_split_halves())
    {
        print $FHAND "Scan for fragments of split tRNAs\n";        
    }

    my $ct = 0;
    my %cms = $cm->main_cm_file_path();
    foreach my $cur_cm (sort keys %cms)
    {
        $ct++;
        if ($ct == 1)
        {
            print $FHAND "Covariance model:                  ".$cms{$cur_cm}."\n";
        }
        else
        {
            print $FHAND "                                   ".$cms{$cur_cm}."\n";
        }
    }

#    if ($self->mito_mode())
#    {
#        if ($cm->organelle_cm_cutoff() != $global_constants->get("organelle_cm_cutoff"))
#        {
#            print $FHAND "tRNA covariance model search       \n"."    cutoff score:                  ".$cm->organelle_cm_cutoff()."\n";
#        }        
#    }
#    else
#    {
        if (((!$self->mito_mode() and !$self->organelle_mode()) and $cm->cm_cutoff() != $global_constants->get("cm_cutoff")) or 
            (($self->mito_mode() or $self->organelle_mode()) and $cm->cm_cutoff() != $global_constants->get("organelle_cm_cutoff")))
        {
            print $FHAND "tRNA covariance model search       \n"."    cutoff score:                  ".$cm->cm_cutoff()."\n";
        }
#    }

    if ($self->{use_prev_ts_run})
    {
        print $FHAND "Using previous\n",
            "tabular output file:               $self->{firstpass_result_file}\n";
    }

    if ($self->{tscan_mode})
    {
        print $FHAND "tRNAscan parameters:               ";
        if ($self->{strict_params})
        {
            print $FHAND "Strict\n";
        }
        else
        {
            print $FHAND "Relaxed\n";
        }
    }

    if ($self->{eufind_mode})
    {
        print $FHAND "EufindtRNA parameters:             ";
        if ($eufind->eufind_params() eq $global_constants->get_subvalue("eufind", "relaxed_param"))
        {
            print $FHAND "Relaxed (Int Cutoff= ".$eufind->eufind_intscore().")\n";
        }
        elsif ($eufind->eufind_params() eq "")
        {
            print $FHAND "Normal\n";
        }
        elsif  ($eufind->eufind_params() eq $global_constants->get_subvalue("eufind", "strict_param"))
        {
            print $FHAND "Strict\n";
        }
        else
        { 
            print $FHAND "?\n";
        }  
    }

    if ($self->{infernal_fp})
    {
        print $FHAND "Infernal first pass cutoff score:  ".$cm->infernal_fp_cutoff()."\n";
    }
    
    if ($self->{padding} != $self->{default_Padding})
    {
        print $FHAND "First-pass tRNA hit padding:       $self->{padding} bp\n";
    }

    if ($self->{alt_gcode})
    {
        print $FHAND "Alternate transl code used:        ", "from file $self->{gc_file}\n";  
    }
    print $FHAND "\nTemporary directory:               ".$global_constants->get("temp_dir")."\n";
    if ($self->{save_all_struct})
    {
        print $FHAND "tRNA secondary structure\n", "    predictions saved to:          ";
        if ($self->{all_struct_file} eq "-")
        {
            print $FHAND "Standard output\n";
        }
        else
        {
            print $FHAND "$self->{all_struct_file}\n";
        }
    }
    if ($self->{bed_file} ne "")
    {
        print $FHAND "tRNA predictions saved to:         ";
        if ($self->{bed_file} eq "-")
        {
            print $FHAND "Standard output\n";
        }
        else
        {
            print $FHAND "$self->{bed_file}\n";
        }
    }
    if ($self->{gff_file} ne "")
    {
        print $FHAND "tRNA predictions saved to:         ";
        if ($self->{gff_file} eq "-")
        {
            print $FHAND "Standard output\n";
        }
        else
        {
            print $FHAND "$self->{gff_file}\n";
        }
    }
    if ($self->{output_fasta_file} ne "")
    {
        print $FHAND "Predicted tRNA sequences\n", "    saved to:                      ";
        if ($self->{output_fasta_file} eq "-")
        {
            print $FHAND "Standard output\n";
        }
        else
        {
            print $FHAND "$self->{output_fasta_file}\n";
        }
    }
    if ($self->{isotype_specific_file} ne "")
    {
        print $FHAND "Isotype specific\n", "    predictions saved to:          ";
        if ($self->{isotype_specific_file} eq "-")
        {
            print $FHAND "Standard output\n";
        }
        else
        {
            print $FHAND "$self->{isotype_specific_file}\n";
        }
    }
    if ($self->{split_fragment_file} ne "")
    {
        print $FHAND "Split tRNA fragment\n", "    predictions saved to:          ";
        if ($self->{split_fragment_file} eq "-")
        {
            print $FHAND "Standard output\n";
        }
        else
        {
            print $FHAND "$self->{split_fragment_file}\n";
        }
    }
    if ($self->{save_odd_struct})
    {
        print $FHAND "Sec structures for tRNAs\n", " with no anticodon predictn: $self->{odd_struct_file}\n";
    }
    if ($self->{save_firstpass_res})
    {
        print $FHAND "First-pass results saved in:      ", "$self->{firstpass_result_file}\n";
    }
    if ($self->{log_file} ne "")
    {
        print $FHAND "Search log saved in:               $self->{log_file}\n";
    }
    if ($self->{save_stats})
    {
        print $FHAND "Search statistics saved in:        $self->{stats_file}\n";
    }
    if ($self->{save_falsepos})
    {
        print $FHAND "False positives saved in:          $self->{falsepos_file}\n";
    }
    if ($self->{save_missed})
    {
        print $FHAND "Seqs with 0 hits saved in:         $self->{missed_seq_file}\n";
    }
    if ($cm->skip_pseudo_filter() | $cm->get_hmm_score() | $tscan->keep_tscan_repeats())
    {
        print $FHAND "\n";
    }
    if ($self->{max_int_len} != $self->{def_max_int_len})
    {
        print $FHAND "Max intron + var. length:          $self->{max_int_len}\n";
    }
    if ($cm->skip_pseudo_filter())
    {
        print $FHAND "Pseudogene checking disabled\n";
    }
    if ($cm->get_hmm_score())
    {
        print $FHAND "Reporting HMM/2' structure score breakdown\n";
    }
    if ($tscan->keep_tscan_repeats())
    {
        print $FHAND "Redundant tRNAscan hits not merged\n";
    } 

    print $FHAND ('-' x 60,"\n\n");
}

1;

