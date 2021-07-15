# tRNAscanSE/Tscan.pm
# This class contains parameters and functions for running tRNAscan used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2017 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::Tscan;

use strict;
use tRNAscanSE::Configuration;
use tRNAscanSE::Utils;
use tRNAscanSE::ArraytRNA;
use tRNAscanSE::tRNA;
use tRNAscanSE::GeneticCode;
use tRNAscanSE::Stats;

sub new
{
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
    
    # set to non-zero if you do NOT want redundant, overlapping hits
    #  found by tRNAscan merged into one hit
    $self->{keep_tscan_repeats} = 0;
        
    $self->{tscan_version} = 1.4;   # version of tRNAscan used by tRNAscan-SE
    $self->{tscan_bin} = "trnascan-1.4";
    $self->{tscan_mask} = 1;        # Bit-wise masks for source of tRNA hits
}

sub set_defaults
{
    my $self = shift;
    my $global_vars = shift;
    my $global_constants = $global_vars->{global_constants};
    
    $self->{tscan_params} = $global_constants->get_subvalue("tscan", "strict_param");
}

sub keep_tscan_repeats
{
    my $self = shift;
    if (@_) { $self->{keep_tscan_repeats} = shift; }
    return $self->{keep_tscan_repeats};
}

sub tscan_params
{
    my $self = shift;
    if (@_) { $self->{tscan_params} = shift; }
    return $self->{tscan_params};
}

sub tscan_version
{
    my $self = shift;
    if (@_) { $self->{tscan_version} = shift; }
    return $self->{tscan_version};
}

sub tscan_bin
{
    my $self = shift;
    if (@_) { $self->{tscan_bin} = shift; }
    return $self->{tscan_bin};
}

sub tscan_mask
{
    my $self = shift;
    return $self->{tscan_mask};
}

sub set_bin
{    
    my $self = shift;
    my $bindir = shift;
    
    # choose correct name for version being run
    # only version 1.4 is provided with distribution

    if ($self->{tscan_version} == 1.4)
    {
        $self->{tscan_bin} = "trnascan-1.4";
    }
    elsif ($self->{tscan_version} == 1.39)
    {
        $self->{tscan_bin} = "trnascan-1.39";
    }
    elsif ($self->{tscan_version} == 2)
    {
        $self->{tscan_bin} = "TRNAscan";
    }
    elsif ($self->{tscan_version} == 1.3)
    {             
        $self->{tscan_bin} = "trnascan-1.3";
    }
    else
    {
        die "FATAL:  Illegal tRNAscan version.\n\n";
    }

    if ($^O =~ /^MSWin/)
    {
        $self->{tscan_bin} .= ".exe";
    }

    if (!(-x $self->{tscan_bin}))
    {
        $self->{tscan_bin} = $bindir."/".$self->{tscan_bin};
        if (!(-x $self->{tscan_bin}))
        {
            die "FATAL: Unable to find ".$self->{tscan_bin}." executable\n\n";
        }
    }
}

sub run_tRNAscan
{    
    my $self = shift;
    my ($tmp_fa, $tmp_raw, $start_index, $lib_dir, $seq_name) = @_;
    my $tscan_version = $self->{tscan_version};
    my $tscan_bin = $self->{tscan_bin};
    my $tscan_params = $self->{tscan_params};

    # version provided with distribution
    if ($tscan_version == 1.4)
    {
        # run default tRNAscan 1.4 using selected param set
        system ("$tscan_bin -i $start_index -c $tscan_params $tmp_fa > $tmp_raw");
        if (&error_exit_status("tRNAscan", $seq_name))
        {
            return -1;
        }
    }
    
    # run tRNAscan without conservative ambiguous base pairing rules
    # not available in distribution version
    elsif ($tscan_version == 1.39)
    {
        system ("$tscan_bin -c $tscan_params $tmp_fa > $tmp_raw"); 
    }

    # run tRNAscan v2.0, not available in distribution version
    elsif ($tscan_version == 2)
    {
        system ("$tscan_bin -SEQ $tmp_fa -TEMPLATE SEtemplate -OUTPUT $tmp_raw > /dev/null");
    }

    # run original tRNAscan 1.3, not available in distribution version
    elsif ($tscan_version == 1.3)
    {             
        if (!(-r "./TPCsignal"))
        {
            system ("ln -s ".$lib_dir."TPCsignal TPCsignal");
        }
        if (!(-r "./Dsignal"))
        {
            system ("ln -s ".$lib_dir."Dsignal Dsignal");
        }
        system ("reformat -ld genbank $tmp_fa > tmp.gb");
        system ("$tscan_bin tmp.gb $tmp_raw > /dev/null");
        system ("rm tmp.gb");
    }
    else
    {
        die "FATAL:  Illegal tRNAscan version.\n\n";
    }
}

# Append tRNAscan verbose output to 
#   result file with header tag
sub append_verbfile
{    
    my $self = shift;
    my ($verb_file, $tmp_fa, $seq_name) = @_;

    &open_for_append(\*TSCANVERB, $verb_file);    
    print TSCANVERB "\n>>>> tRNA-Scan verbose output for <$seq_name>\n\n";
    close TSCANVERB;
    system ("cat tscan.verb.out >> $verb_file");
}

# extract trna hits from raw result file while weeding out repeated hits
# save non-redundant hits in "hit_list" array
sub process_tRNAscan_hits
{    
    my $self = shift;
    my ($global_vars, $seq_name) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $gc = $global_vars->{gc};
    my $fp_tRNAs = $global_vars->{fp_tRNAs};
    my $stats = $global_vars->{stats};
    
    my $tmp_raw = $global_constants->get("tmp_raw");
    my $trna = tRNAscanSE::tRNA->new;
    my $i = 0;

    # open trnascan raw output file for current seq    
    open (TSCANRAW, "$tmp_raw")  || die ("FATAL: Unable to open temp raw output file $tmp_raw\n\n");
    
    # parse one complete hit per call 
    while ($self->parse_tscan_hit($global_constants, $gc, \*TSCANRAW, $trna))
    {
        # if NOT a repeat hit, put it on the hit list 
        if ($self->{keep_tscan_repeats} ||
            (!$self->merge_repeat_hit($global_vars->{stats}, $global_vars->{fp_tRNAs}, $trna)))
        {            
            # check to see if tscan 1.3 has incorrectly reported
            #  start/end index (happens occassionally)             
            if ((abs($trna->end() - $trna->start()) + 1) != length($trna->seq()))
            {
                if ($trna->strand() eq "+")
                {
                    $trna->end($trna->start() + length($trna->seq()) - 1);
                }
                else
                {
                    $trna->end($trna->start() - length($trna->seq()) + 1);
                }
            }
            
            $i = 0;
            while (($i < $fp_tRNAs->get_count()) && ($fp_tRNAs->get($i)->position() < $trna->position()))
            {
                $i++;
            }
            
            $trna->seqname($seq_name);
            $trna->hit_source($self->{tscan_mask});
            $fp_tRNAs->insert($trna, $i);
            $stats->increment_trnatotal();
        }         
        $trna = tRNAscanSE::tRNA->new;
    }        # while (&Parse_tscan_hit), more hits to process for cur seq    
}

sub parse_tscan_hit
{
    my $self = shift;
    my ($global_constants, $gc, $TSCANRAW, $trna) = @_;

    if ($self->{tscan_version} <= 1.4)
    {
        while (<$TSCANRAW>)
        {
            if (/^start position=\s*(\d+)\s*end position=\s*(\d+)/o)
            {
                $trna->start($1);
                $trna->end($2);
                if ($trna->start() < $trna->end())
                {
                    $trna->strand("+");
                    $trna->position($trna->start());
                }
                else
                {
                    $trna->strand("-");
                    $trna->position($global_constants->get("really_big_number") - $trna->start() + 1);
                }
            }
            elsif (/^potential tRNA sequence=\s(.+)\n/o)
            {
                $trna->seq($1);
            }                        
            elsif (/^tRNA predict as a tRNA-\s*(\S+)\s*: anticodon (\S+)/o)
            {
                $trna->isotype($1);
                $trna->anticodon($2);
            }
            elsif (/^anticodon includes unknown bases/o)
            {
                $trna->isotype($gc->undef_isotype());
                $trna->anticodon($gc->undef_anticodon());
            }
            elsif (/^potential intron between positions\s*(\d+)\s*(\d+)/o)
            {
                $trna->add_intron($1, $2, $1, $2, "CI");
            }
            # flag for end of current tRNA hit info
            elsif (/^number of base pairing in the anticodon/o)
            {
                return 1;
            } 
            elsif (/^number of predicted tRNA=(\d+)/o)
            {
                return 0;        # end of hits for this seq 
            }
        }
        return 0;                # reached end of raw hits file
    }                               
    else
    {
        die "FATAL: Illegal tRNAscan version selected.\n\n";
    }
}        

# check current hit for redundancy against all previous hits in hitlist
# if it IS a repeat, merge it with overlapping hit and return 1
# if it doesn't overlap with any hits, return 0
sub merge_repeat_hit
{
    my $self = shift;
    my ($stats, $fp_tRNAs, $trna) = @_;

    for (my $i = 0; $i < $fp_tRNAs->get_count(); $i++)
    {
        if ($trna->strand() eq "+")
        {
            if (($fp_tRNAs->get($i)->strand() eq "+") &&
                (&seg_overlap($trna->start(), $trna->end(), $fp_tRNAs->get($i)->start(), $fp_tRNAs->get($i)->end()))) 
            {
                $fp_tRNAs->get($i)->start(&min($trna->start(), $fp_tRNAs->get($i)->start()));
                $fp_tRNAs->get($i)->end(&max($trna->end(), $fp_tRNAs->get($i)->end()));
                $fp_tRNAs->get($i)->hit_source($fp_tRNAs->get($i)->hit_source() | $self->{tscan_mask});
                $fp_tRNAs->get($i)->isotype($trna->isotype());
                $fp_tRNAs->get($i)->score($trna->score());
    
                # check to see if extended endpoint overlaps i+1 hit's start boundary
                # if so, combine hit[i] and hit[i+1] into one hit and delete hit[i+1]
                if (($i != ($fp_tRNAs->get_count() - 1)) and ($fp_tRNAs->get($i+1)->strand() eq "+")
                    and ($fp_tRNAs->get($i)->end() >= $fp_tRNAs->get($i+1)->start())) 
                {
                    $fp_tRNAs->get($i)->end(&max($fp_tRNAs->get($i)->end(), $fp_tRNAs->get($i+1)->end()));
                    $fp_tRNAs->get($i)->hit_source($fp_tRNAs->get($i)->hit_source() | $fp_tRNAs->get($i+1)->hit_source());

#                    splice(@$r_hit_list,$i+1,1);          # toss out overlapping hit
                    $fp_tRNAs->remove($i+1);
#                    $$r_trnact--;
                    $stats->decrement_trnatotal();
                }   
                return 1;                                 # exit loop immediately
            }
        }
        else         # else (antisense) strand 
        {                
            if (($fp_tRNAs->get($i)->strand() eq "-") &&
                (&seg_overlap($trna->end(), $trna->start(), $fp_tRNAs->get($i)->end(), $fp_tRNAs->get($i)->start()))) 
            {
                $fp_tRNAs->get($i)->start(&max($trna->start(), $fp_tRNAs->get($i)->start()));
                $fp_tRNAs->get($i)->end(&min($trna->end(), $fp_tRNAs->get($i)->end()));
                $fp_tRNAs->get($i)->hit_source($fp_tRNAs->get($i)->hit_source() | $self->{tscan_mask});
                $fp_tRNAs->get($i)->isotype($trna->isotype());
                $fp_tRNAs->get($i)->score($trna->score());

                if (($i != ($fp_tRNAs->get_count() - 1)) and
                    ($fp_tRNAs->get($i)->end() <= $fp_tRNAs->get($i+1)->start()))
                {
                    $fp_tRNAs->get($i)->end(&min($fp_tRNAs->get($i)->end(), $fp_tRNAs->get($i+1)->end()));
                    $fp_tRNAs->get($i)->hit_source($fp_tRNAs->get($i)->hit_source() | $fp_tRNAs->get($i+1)->hit_source());

#                    splice(@$r_hit_list,$i+1,1);          # toss out overlapping hit
                    $fp_tRNAs->remove($i+1);
#                    $$r_trnact--;
                    $stats->decrement_trnatotal();
                }
                return 1;                                 # exit loop immediately
            }
        } # else (antisense) strand
        
    }  # for each (hit)                        

    return 0;                                             # current hit is not a repeat
}

1;