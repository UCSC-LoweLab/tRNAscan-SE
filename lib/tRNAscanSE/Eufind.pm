# tRNAscanSE/Eufind.pm
# This class contains parameters and functions for running eufindtRNA used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2017 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::Eufind;

use strict;
use tRNAscanSE::Configuration;
use tRNAscanSE::Utils;
use tRNAscanSE::ArraytRNA;
use tRNAscanSE::tRNA;
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
        
    $self->{eufind_bin} = "eufindtRNA";
    $self->{eufind_mask} = 2;           # Bit-wise masks for source of tRNA hits
}

sub eufind_params
{
    my $self = shift;
    if (@_) { $self->{eufind_params} = shift; }
    return $self->{eufind_params};
}

sub eufind_intscore
{
    my $self = shift;
    if (@_) { $self->{eufind_intscore} = shift; }
    return $self->{eufind_intscore};
}

sub eufind_bin
{
    my $self = shift;
    if (@_) { $self->{eufind_bin} = shift; }
    return $self->{eufind_bin};
}

sub eufind_mask
{
    my $self = shift;
    return $self->{eufind_mask};
}

sub set_defaults
{
    my $self = shift;
    my $global_vars = shift;
    my $global_constants = $global_vars->{global_constants};
    
    $self->{eufind_params} = $global_constants->get_subvalue("eufind", "relaxed_param");
    $self->{eufind_intscore} = $global_constants->get_subvalue("eufind", "intscore");
}

sub set_bin
{    
    my $self = shift;
    my $bindir = shift;
    
    if ($^O =~ /^MSWin/)
    {
        $self->{eufind_bin} .= ".exe";
    }

    if (!(-x $self->{eufind_bin}))
    {
        $self->{eufind_bin} = $bindir."/".$self->{eufind_bin};
        if (!(-x $self->{eufind_bin}))
        {
            die "FATAL: Unable to find ".$self->{eufind_bin}." executable\n\n";
        }
    }
}

sub run_eufind
{    
    my $self = shift;
    my ($tmp_fa, $start_index, $max_int_len, $seq_name) = @_;
    my $eufind_bin = $self->{eufind_bin};
    my $eufind_intscore = $self->{eufind_intscore};
    my $eufind_params = $self->{eufind_params};

    # run default Eufind using selected param set
    my $eufind_output = `$eufind_bin -i $start_index -F -I $eufind_intscore -l $max_int_len $eufind_params $tmp_fa`;
    if (&error_exit_status("EufindtRNA",$seq_name))
    {
        $eufind_output = "";
    }
    return $eufind_output;
}

sub process_Eufind_hits
{
    my $self = shift;
    my ($global_vars, $eufind_output) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $fp_tRNAs = $global_vars->{fp_tRNAs};
    my $stats = $global_vars->{stats};

    my $trna = undef;
    my $i = 0;
    
    my (@eufind_lines) = split(/\n/, $eufind_output);
    foreach (@eufind_lines)
    {
        if (/^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)/o)
        {
            $trna = tRNAscanSE::tRNA->new;
            $trna->seqname($1);
#            $trnact = $2; 
            $trna->start($3);
            $trna->end($4);
            $trna->isotype($5);
            $trna->anticodon($6);
            $trna->score($9);
            
            if ($trna->start() < $trna->end())
            {
                $trna->position($trna->start());                
                $trna->strand("+");
            }
            else
            { 
                $trna->position($global_constants->get("really_big_number") - $trna->start() + 1);
                $trna->strand("-");
            }
            
            if ($trna->start() == $trna->end())
            {
                print STDERR "Error reading EufindtRNA results: tRNA of length 0"; 
            }
            
            if (!$self->merge_repeat_hit($global_vars->{stats}, $global_vars->{fp_tRNAs}, $trna))
            {
                # insert non-redundant hit in order
                # 'Merge_repeat_hits' depends on list being in order
                $i = 0;
                while (($i < $fp_tRNAs->get_count()) && ($fp_tRNAs->get($i)->position() < $trna->position()))
                {
                    $i++;
                }
                       
                $trna->hit_source($self->{eufind_mask});
                $fp_tRNAs->insert($trna, $i);                
                $stats->increment_trnatotal();
            }
        }
    }
}

# check current hit for redundancy against all previous hits in hitlist
#
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
                $fp_tRNAs->get($i)->hit_source($fp_tRNAs->get($i)->hit_source() | $self->{eufind_mask});
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
                $fp_tRNAs->get($i)->hit_source($fp_tRNAs->get($i)->hit_source() | $self->{eufind_mask});
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
        
    } # for each (hit)                        

    return 0;                                             # current hit is not a repeat
}
1;
