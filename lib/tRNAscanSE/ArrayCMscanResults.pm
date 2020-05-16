# tRNAscanSE/ArrayCMscanResultFile.pm
# This class defines the array of covariance model scanning result files used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2017 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::ArrayCMscanResults;

use strict;
use tRNAscanSE::CMscanResultFile;
use tRNAscanSE::Utils;
use tRNAscanSE::LogFile;

sub new
{
    my $class = shift;
    my $log = shift;
    my $self = {};
    
    $self->{files} = [];    
    $self->{indexes} = [];
    $self->{FILE_H} = undef;
    $self->{file_name} = "";
    $self->{log} = $log;
    
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
    @{$self->{files}} = ();
    @{$self->{indexes}} = ();
}

sub file_name
{
    my $self = shift;
    if (@_) { $self->{file_name} = shift; }
    return $self->{file_name};
}

sub add_file
{
    my $self = shift;
    my ($file_name, $type) = @_;
    my $file = tRNAscanSE::CMscanResultFile->new($file_name, $type);
    push(@{$self->{files}}, $file);
    return (scalar(@{$self->{files}}) - 1);
}

sub get_files
{
    my $self = shift;
    return @{$self->{files}};    
}

sub get_indexes
{
    my $self = shift;
    return @{$self->{indexes}};    
}

sub append_index
{
    my $self = shift;
    my ($file_index, $record_indexes) = @_;
    my $record = [];
    
    for (my $i = 0; $i < scalar(@$record_indexes); $i++)
    {
        $record = [];
        push(@$record, $file_index);
        push(@$record, @{$record_indexes->[$i]});
        push(@{$self->{indexes}}, $record);
    }
}

sub get_result_count
{
    my $self = shift;
    return scalar(@{$self->{indexes}});
}

sub merge_results
{
    my $self = shift;
    my $overlap_range = shift;
    $self->sort_result_files();
    $self->set_initial_indexes();
    $self->merge_result_files($overlap_range);
}

sub sort_result_file
{
    my $self = shift;
    my $file_idx = shift;
    
    if ($file_idx > -1 and $file_idx < scalar(@{$self->{files}}))
    {
        $self->{files}->[$file_idx]->sort_cmsearch_records();
    }
}

sub sort_result_files
{
    my $self = shift;
    
    for (my $i = 0; $i < scalar(@{$self->{files}}); $i++)
    {
        $self->{files}->[$i]->sort_cmsearch_records();
    }
}

sub set_initial_indexes
{
    my $self = shift;
    my $record = [];
    
    if (scalar(@{$self->{files}}) > 0)
    {
        my @record_indexes = $self->{files}->[0]->get_indexes();
        $self->append_index(0, \@record_indexes);
        $self->{files}->[0]->clear_indexes();
    }    
}

sub merge_result_file
{
    my $self = shift;
    my $file_idx = shift;
    my $overlap_range = shift;
    my @record_indexes = ();
    
    if ($file_idx > -1 and $file_idx < scalar(@{$self->{files}}))
    {
        $self->sort_result_file($file_idx);
        if ($file_idx == 0)
        {
            $self->set_initial_indexes();
        }
        else
        {
            @record_indexes = $self->{files}->[$file_idx]->get_indexes();
            $self->append_index($file_idx, \@record_indexes);
            @{$self->{indexes}} = sort sort_by_tRNAscanSE_output @{$self->{indexes}};
            $self->{files}->[$file_idx]->clear_indexes();
        }
        $self->merge_indexes($overlap_range);
    }    
}

sub merge_result_files
{
    my $self = shift;
    my $overlap_range = shift;
    my @record_indexes = ();
    
    for (my $i = 1; $i < scalar(@{$self->{files}}); $i++)
    {
        @record_indexes = $self->{files}->[$i]->get_indexes();
        $self->append_index($i, \@record_indexes);
        @{$self->{indexes}} = sort sort_by_tRNAscanSE_output @{$self->{indexes}};
        $self->merge_indexes($overlap_range);
        $self->{files}->[$i]->clear_indexes();
    }    
}

sub merge_indexes
{
    my $self = shift;
    my $overlap_range = shift;
    
    my $hit_overlap = 0;
    for (my $i = scalar(@{$self->{indexes}}) - 2; $i >= 0; $i--)
    {
        $hit_overlap = eval(($self->{indexes}->[$i]->[2] eq $self->{indexes}->[$i+1]->[2]) && 
                 (&seg_overlap($self->{indexes}->[$i]->[3], $self->{indexes}->[$i]->[4], $self->{indexes}->[$i+1]->[3], $self->{indexes}->[$i+1]->[4], $overlap_range)));
        if ($hit_overlap)
        {
            if ($self->{indexes}->[$i]->[6] >= $self->{indexes}->[$i+1]->[6])
            {
                splice(@{$self->{indexes}}, $i + 1, 1);
            }
            else
            {
                splice(@{$self->{indexes}}, $i, 1);                
            }
        }
    }
}

sub sort_by_tRNAscanSE_output
{
    if ((($a->[5] eq $b->[5]) && ($a->[5] eq "+")) ||
        ($a->[5] ne $b->[5]))
    {
        return ($a->[2] cmp $b->[2] ||
                $a->[5] cmp $b->[5] ||
                $a->[3] <=> $b->[3] ||
                $b->[6] <=> $a->[6]);
    }
    if (($a->[5] eq $b->[5]) && ($a->[5] eq "-"))
    {
        return ($a->[2] cmp $b->[2] ||
                $b->[5] cmp $a->[5] ||
                $b->[4] <=> $a->[4] ||
                $b->[6] <=> $a->[6]);        
    }    
}

sub write_merge_file
{
    my $self = shift;
    my $merge_file = shift;
    my $format = shift;
    
    $self->{file_name} = $merge_file;
    &open_for_write(\$self->{FILE_H}, $self->{file_name});
    my $fh = $self->{FILE_H};
    $self->open_result_files();
    $self->parse_result_files($format);
    $self->close_result_files();
    close($fh);
}

sub open_result_files
{
    my $self = shift;

    for (my $i = 0; $i < scalar(@{$self->{files}}); $i++)
    {
        if (!$self->{files}->[$i]->open_file())
        {
            $self->{log}->error("Fail to open cmsearch output file ".$self->{files}->[$i]->file_name());
        }
    }
}

sub close_result_files
{
    my $self = shift;

    for (my $i = 0; $i < scalar(@{$self->{files}}); $i++)
    {
        $self->{files}->[$i]->close_file();
    }    
}

sub parse_result_files
{
    my $self = shift;
    my $format = shift;
    my $fh = $self->{FILE_H};
    my $file = undef;
    my ($ss, $seq, $model, $nc);
    
    for (my $i = 0; $i < scalar(@{$self->{indexes}}); $i++)
    {
        $file = $self->{files}->[$self->{indexes}->[$i]->[0]];
        ($ss, $seq, $model, $nc) = $file->get_cmsearch_record($self->{indexes}->[$i]->[1], $format);
        
        for (my $j = 2; $j < scalar(@{$self->{indexes}->[$i]}); $j++)
        {
            print $fh $self->{indexes}->[$i]->[$j]."\t";
        }
        print $fh $ss."\t".$seq."\t".$model."\t".$nc."\t".$file->type()."\n";
    }    
}


sub open_file
{
    my $self = shift;
    my $file_name = shift;
    my $success = 0;
    
    if ($file_name ne "")
    {
        $self->{file_name} = $file_name;
        &open_for_read(\$self->{FILE_H}, $self->{file_name});
        $success = 1;
    }
    else
    {
        die "Merge result file name is not set.\n"
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

sub get_next_cmsearch_hit
{
    my $self = shift;
    my ($trna) = @_;
    my $fh = $self->{FILE_H};
    
    my $line = "";
    my $cm_model = "";
    my @columns = ();
    $trna->clear();
    
    if (defined $fh and $line = <$fh>)
    {
        chomp($line);
        @columns = split(/\t/, $line);

        $trna->seqname($columns[0]);
        $trna->score($columns[4]);
        if ($columns[3] eq "+")
        {
            $trna->start($columns[1]);
            $trna->end($columns[2]);
        }
        else
        {
            $trna->start($columns[2]);
            $trna->end($columns[1]);
        }
        $trna->strand($columns[3]);
        $trna->trunc($columns[5]);
        $trna->ss($columns[6]);
        $trna->seq($columns[7]);
        $trna->model($columns[10]);
        $trna->hit_source("Inf");
        $cm_model = $columns[8];
    }
    
    return $cm_model;
}

1;
