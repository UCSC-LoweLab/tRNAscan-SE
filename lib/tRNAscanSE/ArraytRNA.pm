# tRNAscanSE/ArraytRNA.pm
# This class contains parameters and functions describing an array of tRNA genes used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2017 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::ArraytRNA;

use strict;
use tRNAscanSE::tRNA;
use tRNAscanSE::LogFile;
use tRNAscanSE::Options;
use tRNAscanSE::Configuration;

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
    
    $self->{ar_tRNAs} = [];
    $self->{h_ordered_seqname} = {};
}

sub clear
{
    my $self = shift;
    @{ $self->{ar_tRNAs} } = ();
    %{$self->{h_ordered_seqname}} = ();
}

sub array
{
    my $self = shift;
    if (@_) { @{ $self->{ar_tRNAs} } = @_; }
    return $self->{ar_tRNAs};
}

sub put
{
    my $self = shift;
	my ($obj) = @_;
	push(@{$self->{ar_tRNAs}}, $obj);
}

sub get
{
    my $self = shift;
	my ($index) = @_;
    if ($index > -1 && $index < $self->get_count())
    {
        return $self->{ar_tRNAs}->[$index];
    }
    else
    {
        return undef;
    }
}

sub insert
{
    my $self = shift;
	my ($obj, $index) = @_;
    
    if ($index > -1 and $index <= $self->get_count())
    {
        splice(@{ $self->{ar_tRNAs} }, $index, 0, $obj);
    }    
}

sub remove
{
    my $self = shift;
    my $index = shift;
    
    if ($index > -1 and $index < $self->get_count())
    {
        splice(@{ $self->{ar_tRNAs} }, $index, 1);
    }
}

sub get_count
{
	my $self = shift;
	return (scalar(@{$self->{ar_tRNAs}}));
}

sub get_seq_order_index
{
    my $self = shift;
    my $seqname = shift;
    
    my $index = -1;
    if (defined $self->{h_ordered_seqname}->{$seqname})
    {
        $index = $self->{h_ordered_seqname}->{$seqname};
    }
    return $index;
}

sub reorder_all_tRNA_id
{
    my $self = shift;
    
    my $ct = 0;    
    my $prev_seqname = "";
    $self->sort_array("tRNAscanSE_output");
    for (my $i = 0; $i < $self->get_count(); $i++)
    {
        if ($self->{ar_tRNAs}->[$i]->seqname() ne $prev_seqname)
        {
            $ct = 0;
        }

        $ct++;
        $self->{ar_tRNAs}->[$i]->id($ct);
        $self->{ar_tRNAs}->[$i]->tRNAscan_id($self->{ar_tRNAs}->[$i]->seqname().".tRNA".$ct."-".$self->{ar_tRNAs}->[$i]->isotype().$self->{ar_tRNAs}->[$i]->anticodon());
        $prev_seqname = $self->{ar_tRNAs}->[$i]->seqname();
    }
}

sub reorder_tRNA_id
{
    my $self = shift;
    my $seqname = shift;
    
    my $ct = 0;    
    $self->sort_array("tRNAscanSE_output");
    for (my $i = 0; $i < $self->get_count(); $i++)
    {
        if ($self->{ar_tRNAs}->[$i]->{seqname} lt $seqname)
        {}
        elsif ($self->{ar_tRNAs}->[$i]->{seqname} gt $seqname)
        {
            last;
        }
        else
        {
            $ct++;
            $self->{ar_tRNAs}->[$i]->{id} = $ct;
            $self->{ar_tRNAs}->[$i]->{tRNAscan_id} = $seqname.".tRNA".$ct."-".$self->{ar_tRNAs}->[$i]->{isotype}.$self->{ar_tRNAs}->[$i]->{anticodon};
        }
    }
}

sub bsearch_num_field
{
    my $self = shift;
    my $x = shift;
    my $key = shift;
    my ($l, $u) = (0, @{$self->{ar_tRNAs}} - 1);  
    my $i;                       
    while ($l <= $u)
	{
		$i = int(($l + $u)/2);
        if ($key eq "id")
        {
            if ($self->{ar_tRNAs}->[$i]->id() < $x)
            {
                $l = $i+1;
            }
            elsif ($self->{ar_tRNAs}->[$i]->id() > $x)
            {
                $u = $i-1;
            } 
            else
            {
                return $i; 
            }
        }
        elsif ($key eq "position")
        {
            if ($self->{ar_tRNAs}->[$i]->position() < $x)
            {
                $l = $i+1;
            }
            elsif ($self->{ar_tRNAs}->[$i]->position() > $x)
            {
                $u = $i-1;
            } 
            else
            {
                return $i; 
            }
        }
    }
    return -1;         	
}

sub bsearch_id
{
    my $self = shift;
    my $x = shift;
    my $key = shift;
    my ($l, $u) = (0, @{$self->{ar_tRNAs}} - 1);  
    my $i;

    my @x_parts = ();
    if ($key eq "gtrnadb_id_parts")
    {
        @x_parts = split(/\-/, $x);
    }

    while ($l <= $u)
	{
		$i = int(($l + $u)/2);
        if ($key eq "tRNAscan_id")
        {
            if ($self->{ar_tRNAs}->[$i]->tRNAscan_id() lt $x)
            {
                $l = $i+1;
            }
            elsif ($self->{ar_tRNAs}->[$i]->tRNAscan_id() gt $x)
            {
                $u = $i-1;
            } 
            else
            {
                return $i; 
            }
        }
        elsif ($key eq "tRNAscanid")
        {
            if (substr($self->{ar_tRNAs}->[$i]->tRNAscan_id(), 0, rindex($self->{ar_tRNAs}->[$i]->tRNAscan_id(), "-")) lt $x)
            {
                $l = $i+1;
            }
            elsif (substr($self->{ar_tRNAs}->[$i]->tRNAscan_id(), 0, rindex($self->{ar_tRNAs}->[$i]->tRNAscan_id(), "-")) gt $x)
            {
                $u = $i-1;
            } 
            else
            {
                return $i; 
            }
        }
        elsif ($key eq "gtrnadb_id")
        {
            if ($self->{ar_tRNAs}->[$i]->gtrnadb_id() lt $x)
            {
               $l = $i+1;
            }
            elsif ($self->{ar_tRNAs}->[$i]->gtrnadb_id() gt $x)
            {
                $u = $i-1;
            } 
            else
            {
                return $i; 
            }
        }
        elsif ($key eq "gtrnadb_id_parts")
        {
            my @ar_parts = split(/\-/, $self->{ar_tRNAs}->[$i]->gtrnadb_id());
            if ($ar_parts[0] lt $x_parts[0])
            {
                $l = $i+1;
            }
            elsif ($ar_parts[0] gt $x_parts[0])
            {
                $u = $i-1;
            }
            else
            {
                if ($ar_parts[1] lt $x_parts[1])
                {
                    $l = $i+1;
                }
                elsif ($ar_parts[1] gt $x_parts[1])
                {
                    $u = $i-1;
                }
                else
                {
                    if ($ar_parts[2] lt $x_parts[2])
                    {
                        $l = $i+1;
                    }
                    elsif ($ar_parts[2] gt $x_parts[2])
                    {
                        $u = $i-1;
                    }
                    else
                    {
                        if ($ar_parts[3] < $x_parts[3])
                        {
                            $l = $i+1;
                        }
                        elsif ($ar_parts[3] > $x_parts[3])
                        {
                            $u = $i-1;
                        }
                        else
                        {
                            if ($ar_parts[4] < $x_parts[4])
                            {
                                $l = $i+1;
                            }
                            elsif ($ar_parts[4] > $x_parts[4])
                            {
                                $u = $i-1;
                            }
                            else
                            {
                                return $i;
                            }
                        }
                    }
                }
            }
        }
        elsif ($key eq "extdb_id")
        {
            if ($self->{ar_tRNAs}->[$i]->extdb_id() lt $x)
            {
                $l = $i+1;
            }
            elsif ($self->{ar_tRNAs}->[$i]->extdb_id() gt $x)
            {
                $u = $i-1;
            } 
            else
            {
                return $i; 
            }
        }
    }
    return -1;         	
}

sub sort_array
{
    my $self = shift;
    my $key = shift;
    
    if ($key eq "tRNAscan_id")
    {
        @{$self->{ar_tRNAs}} = sort sort_by_tRNAscan_id @{$self->{ar_tRNAs}};
    }
    elsif ($key eq "tRNAscanid")
    {
        @{$self->{ar_tRNAs}} = sort sort_by_tRNAscanid @{$self->{ar_tRNAs}};
    }
    elsif ($key eq "gtrnadb_id")
    {
        @{$self->{ar_tRNAs}} = sort sort_by_gtrnadb_id @{$self->{ar_tRNAs}};
    }
    elsif ($key eq "gtrnadb_id_parts")
    {
        @{$self->{ar_tRNAs}} = sort sort_by_gtrnadb_id_parts @{$self->{ar_tRNAs}};
    }
    elsif ($key eq "extdb_id")
    {
        @{$self->{ar_tRNAs}} = sort sort_by_extdb_id @{$self->{ar_tRNAs}};
    }
    elsif ($key eq "seqname_start")
    {
        @{$self->{ar_tRNAs}} = sort sort_by_seqname_start @{$self->{ar_tRNAs}};
    }
    elsif ($key eq "coord")
    {
        @{$self->{ar_tRNAs}} = sort sort_by_coord @{$self->{ar_tRNAs}};
    }
    elsif ($key eq "strand")
    {
        @{$self->{ar_tRNAs}} = sort sort_by_strand_coord @{$self->{ar_tRNAs}};
    }
    elsif ($key eq "score")
    {
        @{$self->{ar_tRNAs}} = sort sort_by_score_coord @{$self->{ar_tRNAs}};
    }
    elsif ($key eq "isotype")
    {
        @{$self->{ar_tRNAs}} = sort sort_by_isotype @{$self->{ar_tRNAs}};
    }
    elsif ($key eq "ac")
    {
        @{$self->{ar_tRNAs}} = sort sort_by_ac_matscore @{$self->{ar_tRNAs}};
    }
    elsif ($key eq "tRNAscanSE_output")
    {
        @{$self->{ar_tRNAs}} = sort sort_by_tRNAscanSE_output @{$self->{ar_tRNAs}};
    }
}

sub sort_by_tRNAscan_id
{
    return ($a->tRNAscan_id() cmp $b->tRNAscan_id()); 
}

sub sort_by_tRNAscanid
{
    return (substr($a->tRNAscan_id(), 0, rindex($a->tRNAscan_id(), "-")) cmp substr($b->tRNAscan_id(), 0, rindex($b->tRNAscan_id(), "-"))); 
}

sub sort_by_gtrnadb_id
{
    return ($a->gtrnadb_id() cmp $b->gtrnadb_id()); 
}

sub sort_by_gtrnadb_id_parts
{
    my @a_parts = split(/\-/, $a->gtrnadb_id());
    my @b_parts = split(/\-/, $b->gtrnadb_id());

    return ($a_parts[0] cmp $b_parts[0] ||
            $a_parts[1] cmp $b_parts[1] ||
            $a_parts[2] cmp $b_parts[2] ||
            int($a_parts[3]) <=> int($b_parts[3]) ||
            int($a_parts[4]) <=> int($b_parts[4]));
}

sub sort_by_extdb_id
{
    return ($a->extdb_id() cmp $b->extdb_id()); 
}

sub sort_by_seqname_start
{
    return ($a->seqname() cmp $b->seqname() ||
            $a->start() <=> $b->start());
}

sub sort_by_coord
{
    return ($a->ordered_seqname() <=> $b->ordered_seqname() ||
        $a->start() <=> $b->start());
}

sub sort_by_strand_coord
{
    return ($a->strand() cmp $b->strand() ||
            $a->ordered_seqname() <=> $b->ordered_seqname() ||
            $a->start() <=> $b->start());
}

sub sort_by_score_coord
{
    return ($a->is_mito() <=> $b->is_mito() ||
            $a->score() <=> $b->score() ||
            $a->ordered_seqname() <=> $b->ordered_seqname() ||
            $a->start() <=> $b->start());
}

sub sort_by_ac_matscore
{
    return ($a->is_mito() <=> $b->is_mito() ||
            $a->is_numt() <=> $b->is_numt() ||
            $a->isotype() cmp $b->isotype() ||
            $a->anticodon() cmp $b->anticodon() ||
            $a->mat_score() <=> $b->mat_score() ||
            $a->seqname() cmp $b->seqname() ||
            $a->start() <=> $b->start());
}

sub sort_by_isotype
{
    return ($a->isotype() cmp $b->isotype() ||
            $a->anticodon() cmp $b->anticodon() ||
            $a->mat_score() <=> $b->mat_score() ||
            $a->seqname() cmp $b->seqname() ||
            $a->start() <=> $b->start());
}

sub sort_by_tRNAscanSE_output
{
    if ((($a->strand() eq $b->strand()) && ($a->strand() eq "+")) ||
        ($a->strand() ne $b->strand()))
    {
        return ($a->ordered_seqname() <=> $b->ordered_seqname() ||
                $a->strand() cmp $b->strand() ||
                $a->start() <=> $b->start());
    }
    if (($a->strand() eq $b->strand()) && ($a->strand() eq "-"))
    {
        return ($a->ordered_seqname() <=> $b->ordered_seqname() ||
                $b->strand() cmp $a->strand() ||
                $b->end() <=> $a->end());        
    }    
}

sub set_ext_seqname_order
{
    my $self = shift;
    my $ar = shift;
    
    if (scalar(keys %{$self->{h_ordered_seqname}}) > 0)
    {
        for (my $i = 0; $i < scalar(@$ar); $i++)
        {
            my $seqname = $ar->[$i]->[0];
            unshift(@{$ar->[$i]}, $self->{h_ordered_seqname}->{$seqname});
        }
    }
}

sub set_seqname_order
{
    my $self = shift;
    my $file = shift;
    my $log = shift;

    $self->read_seqname_order_list($file, $log);
    
    if (scalar(keys %{$self->{h_ordered_seqname}}) > 0)
    {
        for (my $i = 0; $i < scalar(@{$self->{ar_tRNAs}}); $i++)
        {
            my $seqname = $self->{ar_tRNAs}->[$i]->seqname();
            $self->{ar_tRNAs}->[$i]->ordered_seqname($self->{h_ordered_seqname}->{$seqname});
        }
    }
}

sub read_seqname_order_list
{
    my $self = shift;
    my $file = shift;
    my $log = shift;
    
    my $ct = 0;
    my $line = "";
    
	$log->status("Reading sequnece name order list from $file");
	open(FILE_IN, "$file") or die "Error: Fail to open $file\n";
	while ($line = <FILE_IN>)
	{
		chomp($line);
		if ($line !~ /^#/ && $line ne "")
		{
            $ct++;
            $self->{h_ordered_seqname}->{$line} = $ct;
		}
	}
	
	close(FILE_IN);    
}

sub intersect
{
    my $self = shift;
    my $ar = shift;
    my $mode = shift;
    
    my $pair = [];
    my @pairs = ();
    
    # assume input array is pre-sorted
    # sort tRNAs by coord
    $self->sort_array("coord");
    
    my $tRNA_ct = 0;
    for (my $i = 0; $i < scalar(@$ar); $i++)
    {
        if ($tRNA_ct < $self->get_count())
        {
            if ($ar->[$i]->[0] == $self->{ar_tRNAs}->[$tRNA_ct]->ordered_seqname())
            {
                if ($mode eq "full")
                {                
                    if ($ar->[$i]->[2] >= $self->{ar_tRNAs}->[$tRNA_ct]->start() && $ar->[$i]->[3] <= $self->{ar_tRNAs}->[$tRNA_ct]->end())
                    {
                        $pair = [];
                        push(@$pair, $i);
                        push(@$pair, $tRNA_ct);
                        push(@pairs, $pair);
                        $tRNA_ct++;
                    }
                }
                elsif($mode eq "any")
                {
                    if (($ar->[$i]->[2] >= $self->{ar_tRNAs}->[$tRNA_ct]->start() && $ar->[$i]->[2] <= $self->{ar_tRNAs}->[$tRNA_ct]->end()) ||
                        ($ar->[$i]->[3] >= $self->{ar_tRNAs}->[$tRNA_ct]->start() && $ar->[$i]->[3] <= $self->{ar_tRNAs}->[$tRNA_ct]->end()))
                    {
                        $pair = [];
                        push(@$pair, $i);
                        push(@$pair, $tRNA_ct);
                        push(@pairs, $pair);
                        $tRNA_ct++;
                    }                    
                }
            }
            elsif ($ar->[$i]->[0] > $self->{ar_tRNAs}->[$tRNA_ct]->ordered_seqname())
            {
                $tRNA_ct++;
            }
        }
    }
    
    return \@pairs;
}

sub strand_specific_intersect
{
    my $self = shift;
    my $ar = shift;
    my $mode = shift;
    
    my $pair = [];
    my @pairs = ();
    
    # assume input array is pre-sorted
    # sort tRNAs by strand and coord
    $self->sort_array("strand");
    
    my $tRNA_ct = 0;
    for (my $i = 0; $i < scalar(@$ar); $i++)
    {
        if ($tRNA_ct < $self->get_count())
        {
            if ($ar->[$i]->[4] eq $self->{ar_tRNAs}->[$tRNA_ct]->strand())
            {
                if ($ar->[$i]->[0] == $self->{ar_tRNAs}->[$tRNA_ct]->ordered_seqname())
                {
                    if ($mode eq "full")
                    {                
                        if ($ar->[$i]->[2] >= $self->{ar_tRNAs}->[$tRNA_ct]->start() && $ar->[$i]->[3] <= $self->{ar_tRNAs}->[$tRNA_ct]->end())
                        {
                            $pair = [];
                            push(@$pair, $i);
                            push(@$pair, $tRNA_ct);
                            push(@pairs, $pair);
                            $tRNA_ct++;
                        }
                    }
                    elsif($mode eq "any")
                    {
                        if (($ar->[$i]->[2] >= $self->{ar_tRNAs}->[$tRNA_ct]->start() && $ar->[$i]->[2] <= $self->{ar_tRNAs}->[$tRNA_ct]->end()) ||
                            ($ar->[$i]->[3] >= $self->{ar_tRNAs}->[$tRNA_ct]->start() && $ar->[$i]->[3] <= $self->{ar_tRNAs}->[$tRNA_ct]->end()))
                        {
                            $pair = [];
                            push(@$pair, $i);
                            push(@$pair, $tRNA_ct);
                            push(@pairs, $pair);
                            $tRNA_ct++;
                        }                    
                    }
                }
                elsif ($ar->[$i]->[0] > $self->{ar_tRNAs}->[$tRNA_ct]->ordered_seqname())
                {
                    $tRNA_ct++;
                }
            }
            elsif ($ar->[$i]->[4] gt $self->{ar_tRNAs}->[$tRNA_ct]->strand())
            {
                $tRNA_ct++;
            }
        }
    }
    
    return \@pairs;
}


1;
