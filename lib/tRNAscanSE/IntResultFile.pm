# tRNAscanSE/IntResultFile.pm
# This class defines the intermediate result file used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2017 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::IntResultFile;

use strict;
use tRNAscanSE::Utils;
use tRNAscanSE::tRNA;

sub new
{
    my $class = shift;
    my $file_name = shift;
    my $self = {
        file_name => $file_name,
        FILE_H => undef,
        log_file => 0
    };

    $self->{process_id} = $$;

    if ($file_name eq "")
    {
        $self->{file_name} = "tscan$$"."_sp.out";
    }
    
    $self->{indexes} = [];
    
    bless ($self, $class);
    return $self;
}

sub DESTROY
{
    my $self = shift;
}

sub clear
{
    my $self = shift;
    $self->{file_name} = "";
    @{$self->{indexes}} = ();
}

sub clear_index
{
    my $self = shift;
    @{$self->{indexes}} = ();
}

sub file_name
{
    my $self = shift;
    if (@_) { $self->{file_name} = shift; }
    return $self->{file_name};
}

sub get_indexes
{
    my $self = shift;
    return @{$self->{indexes}};    
}

sub get_count
{
    my $self = shift;
    return scalar(@{$self->{indexes}});        
}

sub get_pos
{
    my $self = shift;
    my $index = shift;
    my $pos = undef;
    
    if ($index > -1 and $index < scalar(@{$self->{indexes}}))
    {
        $pos = $self->{indexes}->[$index]->[0];
    }
    return $pos;
}

sub open_file
{
    my $self = shift;
    my $mode = shift;
    
    my $success = 0;
    
    if ($self->{file_name} ne "")
    {
        if ($mode eq "read")
        {
            &open_for_read(\$self->{FILE_H}, $self->{file_name});
        }
        elsif ($mode eq "write")
        {
            &open_for_write(\$self->{FILE_H}, $self->{file_name});        
            my $fh = $self->{FILE_H};
            print $fh "seqname\tordered_seqname\tid\ttRNAscan_id\tstart\tend\tstrand\tstart2\tend2\tstrand2\tstart3\tend3\tstrand3\tposition\tisotype\tanticodon\t".
                "ac_pos\tintron\tmodel\tscore\tmat_score\thmm_score\tss_score\tcove_score\tcove_mat_score\tcove_hmm_score\tcove_ss_score\t".
                "inf_score\tinf_mat_score\tinf_hmm_score\tinf_ss_score\tpseudo\ttrunc\tcategory\t".
                "hit_source\tsrc_seqid\tsrc_seq_len\tseq\tmat_seq\tss\tmat_ss\n";
        }
        elsif ($mode eq "append")
        {
            &open_for_append(\$self->{FILE_H}, $self->{file_name});        
        }
        $success = 1;
    }
    else
    {
        die "Intermediate result file name is not set.\n"
    }

    return $success;
}

sub close_file
{
    my $self = shift;
    
    if (defined $self->{FILE_H})
    {
        close($self->{FILE_H});
    }
}

sub write_line
{
    my $self = shift;
    my $line = shift;
    
    my $fh = $self->{FILE_H};
    
    print $fh $line . "\n";
}

sub write_tRNA
{
    my $self = shift;
    my $tRNA = shift;
    
    my $fh = $self->{FILE_H};

    my $filepos = tell($fh);
    my $record = [];
    push(@$record, $filepos);
    push(@$record, $tRNA->strand()); # strand
    push(@$record, $tRNA->ordered_seqname()); # ordered_seqname
    push(@$record, $tRNA->start()); # start
    push(@$record, $tRNA->end()); # end
    push(@{$self->{indexes}}, $record);
    
    my @ar_ac_pos = $tRNA->ar_ac_pos();
    my @ar_introns = $tRNA->ar_introns();
    my $cove = $tRNA->get_domain_model("cove");
    my $inf = $tRNA->get_domain_model("infernal");
    
    print $fh $tRNA->seqname()."\t".$tRNA->ordered_seqname()."\t".$tRNA->id()."\t".$tRNA->tRNAscan_id()."\t".$tRNA->start()."\t".$tRNA->end()."\t".$tRNA->strand()."\t".
        $tRNA->exon_start(2)."\t".$tRNA->exon_end(2)."\t".$tRNA->exon_strand(2)."\t".$tRNA->exon_start(3)."\t".$tRNA->exon_end(3)."\t".$tRNA->exon_strand(3)."\t".
        $tRNA->position()."\t".$tRNA->isotype()."\t".$tRNA->anticodon()."\t";
    
    for (my $i = 0; $i < scalar(@ar_ac_pos); $i++)
    {
        print $fh $ar_ac_pos[$i]->{rel_start}."-".$ar_ac_pos[$i]->{rel_end}.";";
    }
    print $fh "\t";
    
    for (my $i = 0; $i < scalar(@ar_introns); $i++)
    {
        print $fh $ar_introns[$i]->{rel_start}."-".$ar_introns[$i]->{rel_end}.",".$ar_introns[$i]->{start}."-".$ar_introns[$i]->{end}.",".$ar_introns[$i]->{type}.",".$ar_introns[$i]->{seq}.";";
    }
    print $fh "\t";

    print $fh $tRNA->model()."\t".$tRNA->score()."\t".$tRNA->mat_score()."\t".$tRNA->hmm_score()."\t".$tRNA->ss_score()."\t";
    if (defined $cove)
    {
        print $fh $cove->{score}."\t".$cove->{mat_score}."\t".$cove->{hmm_score}."\t".$cove->{ss_score}."\t";
    }
    else
    {
        print $fh "\t\t\t\t";
    }
    if (defined $inf)
    {
        print $fh $inf->{score}."\t".$inf->{mat_score}."\t".$inf->{hmm_score}."\t".$inf->{ss_score}."\t";
    }
    else
    {
        print $fh "\t\t\t\t";
    }
    print $fh $tRNA->pseudo()."\t".$tRNA->trunc()."\t".$tRNA->category().
    "\t".$tRNA->hit_source()."\t".$tRNA->src_seqid()."\t".$tRNA->src_seqlen()."\t".$tRNA->seq()."\t".$tRNA->mat_seq()."\t".$tRNA->ss()."\t".$tRNA->mat_ss()."\t".$tRNA->note()."\n";
}

sub get_tRNA
{
    my $self = shift;
    my $filepos = shift;
    my $tRNA = shift;
    
    my $fh = $self->{FILE_H};
    
    $tRNA->clear();
    seek($fh, $filepos, 0);
    my $line = <$fh>;
    chomp($line);
    my @columns = split(/\t/, $line, -1);
    
    $tRNA->seqname($columns[0]);
    $tRNA->ordered_seqname($columns[1]);
    $tRNA->id($columns[2]);
    $tRNA->tRNAscan_id($columns[3]);
    $tRNA->start($columns[4]);
    $tRNA->end($columns[5]);
    $tRNA->strand($columns[6]);
    $tRNA->exon_start(2, $columns[7]);
    $tRNA->exon_end(2, $columns[8]);
    $tRNA->exon_strand(2, $columns[9]);
    $tRNA->exon_start(3, $columns[10]);
    $tRNA->exon_end(3, $columns[11]);
    $tRNA->exon_strand(3, $columns[12]);
    $tRNA->position($columns[13]);
    $tRNA->isotype($columns[14]);
    $tRNA->anticodon($columns[15]);
    my @ar_ac_pos = split(/\;/, $columns[16]);
    for (my $i = 0; $i < scalar(@ar_ac_pos); $i++)
    {
        if ($ar_ac_pos[$i] =~ /^(\d+)-(\d+)$/)
        {
            $tRNA->add_ac_pos($1, $2);
        }        
    }
    if ($columns[17] ne "")
    {
        my @ar_introns = split(/\;/, $columns[17]);
        for (my $i = 0; $i < scalar(@ar_introns); $i++)
        {
            my @parts = split(/\,/, $ar_introns[$i], -1);
            my ($rel_start, $rel_end, $start, $end) = (0, 0, 0, 0);
            if ($parts[0] =~ /^(\d+)-(\d+)$/)
            {
                $rel_start = $1;
                $rel_end = $2;
            }
            if ($parts[1] =~ /^(\d+)-(\d+)$/)
            {
                $start = $1;
                $end = $2;
            }
            $tRNA->add_intron($rel_start, $rel_end, $start, $end, $parts[2], $parts[3]);
        }
    }
    $tRNA->model($columns[18]);
    $tRNA->score($columns[19]);
    $tRNA->mat_score($columns[20]);
    $tRNA->hmm_score($columns[21]);
    $tRNA->ss_score($columns[22]);
    if ($columns[23] ne "")
    {
        $tRNA->set_domain_model("cove", $columns[23]);
        $tRNA->update_domain_model("cove", $columns[23], $columns[24], $columns[25], $columns[26])
    }
    if ($columns[27] ne "")
    {
        $tRNA->set_domain_model("infernal", $columns[27]);
        $tRNA->update_domain_model("infernal", $columns[27], $columns[28], $columns[29], $columns[30])
    }
    $tRNA->pseudo($columns[31]);
    $tRNA->trunc($columns[32]);
    $tRNA->category($columns[33]);
    $tRNA->hit_source($columns[34]);
    $tRNA->src_seqid($columns[35]);
    $tRNA->src_seqlen($columns[36]);
    $tRNA->seq($columns[37]);
    $tRNA->mat_seq($columns[38]);
    $tRNA->ss($columns[39]);
    $tRNA->mat_ss($columns[40]);
    $tRNA->note($columns[41]);}

sub sort_records
{
    my $self = shift;
    my $key = shift;
    
    if ($key eq "tRNAscan_id")
    {
        $self->build_tRNAscan_id_index();
    }
    elsif ($key eq "tRNAscanSE_output")
    {
        $self->build_tRNAscanSE_output_index();
    }
    elsif ($key eq "bed_output")
    {
        $self->sort_bed_output_index();
    }
}

sub build_tRNAscan_id_index
{
    my $self = shift;
    my $line = "";
    my @columns = ();
    my $filepos = undef;
    my $record = [];
    my @indexes = ();
    
    if ($self->open_file("read"))
    {
        my $fh = $self->{FILE_H};
        $line = <$fh>;
        $filepos = tell($fh);
        while ($line = <$fh>)
        {
            $record = [];
            @columns = split(/\t/, $line);
            push(@$record, $filepos);
            push(@$record, $columns[3]);
            push(@indexes, $record);
            $filepos = tell($fh);
        }        
        $self->close_file();            
    }
    
    @{$self->{indexes}} = sort {$a->[1] cmp $b->[1]} @indexes;
}

sub build_tRNAscanSE_output_index
{
    my $self = shift;
    my $line = "";
    my @columns = ();
    my $filepos = undef;
    my $record = [];
    my @indexes = ();
    
    if ($self->open_file("read"))
    {
        my $fh = $self->{FILE_H};
        $line = <$fh>;
        $filepos = tell($fh);
        while ($line = <$fh>)
        {
            $record = [];
            @columns = split(/\t/, $line);
            push(@$record, $filepos);
            push(@$record, $columns[6]); # strand
            push(@$record, $columns[1]); # ordered_seqname
            push(@$record, $columns[4]); # start
            push(@$record, $columns[5]); # end
            push(@indexes, $record);
            $filepos = tell($fh);
        }        
        $self->close_file();            
    }
    
    @{$self->{indexes}} = sort sort_by_tRNAscanSE_output @indexes;    
}

sub sort_bed_output_index
{
    my $self = shift;    
    @{$self->{indexes}} = sort sort_by_bed_output @{$self->{indexes}};    
}

sub sort_by_tRNAscanSE_output
{
    if ((($a->[1] eq $b->[1]) && ($a->[1] eq "+")) ||
        ($a->[1] ne $b->[1]))
    {
        return ($a->[2] <=> $b->[2] ||
                $a->[1] cmp $b->[1] ||
                $a->[3] <=> $b->[3]);
    }
    if (($a->[1] eq $b->[1]) && ($a->[1] eq "-"))
    {
        return ($a->[2] <=> $b->[2] ||
                $b->[1] cmp $a->[1] ||
                $b->[4] <=> $a->[4]);        
    }    
}

sub sort_by_bed_output
{
    return ($a->[2] <=> $b->[2] ||
            $a->[3] <=> $b->[3]);
}

sub bsearch_tRNAscan_id
{
    my $self = shift;
    my $x = shift;
    my ($l, $u) = (0, @{$self->{indexes}} - 1);  
    my $i;                       
    while ($l <= $u)
	{
		$i = int(($l + $u)/2);
		if ($self->{indexes}->[$i]->[1] lt $x)
		{
		    $l = $i+1;
		}
		elsif ($self->{indexes}->[$i]->[1] gt $x)
		{
		    $u = $i-1;
		} 
		else
		{
			return $i; 
		}
    }
    return -1;         	    
}

1;
