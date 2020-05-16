# tRNAscanSE/CMscanResultFile.pm
# This class defines the covariance model scanning result file used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2017 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::CMscanResultFile;

use strict;
use tRNAscanSE::Utils;

sub new
{
    my $class = shift;
    my $file_name = shift;
    my $type = shift;
    my $self = {
        file_name => $file_name,
        type => $type,
        FILE_H => undef
    };
    
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
    $self->{type} = "";
    @{$self->{indexes}} = ();
}

sub file_name
{
    my $self = shift;
    if (@_) { $self->{file_name} = shift; }
    return $self->{file_name};
}

sub type
{
    my $self = shift;
    if (@_) { $self->{type} = shift; }
    return $self->{type};
}

sub get_indexes
{
    my $self = shift;
    return @{$self->{indexes}};    
}

sub clear_indexes
{
    my $self = shift;
    @{$self->{indexes}} = ();
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
    
    my $success = 0;
    
    if ($self->{file_name} ne "")
    {
        &open_for_read(\$self->{FILE_H}, $self->{file_name});
        $success = 1;
    }
    else
    {
        die "$self->{type} result file name is not set.\n"
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

sub sort_cmsearch_records
{
    my $self = shift;
    my $line = "";
    my $filepos = undef;
    my $record = [];
    my @indexes = ();
    my $seq_name = "";
    my $score = 0;
    my $start = 0;
    my $end = 0;
    my $strand = "";
    my $trunc = "";
    my $pos = 0;
    
    if ($self->open_file("read"))
    {
        my $fh = $self->{FILE_H};
        $filepos = tell($fh);
        while ($line = <$fh>)
        {
            if ($line =~ /Hit alignments:/)
            {
                $start = 1;
            }
            
            if ($start)
            {
                if ($line =~ /^>>\s*(\S+)/)
                {
                    $seq_name = $1;
                    $pos = $filepos;
                }
                elsif ($line =~ /^\s+\(\d+\)\s+\S+\s+([e0-9.\-]+)\s+([0-9.\-]+)\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)\s+([+-])\s+\S+\s+\S+\s+(\S+)\s+([0-9.]+)/)
                {
                    $score = $2;
                    $start = $5;
                    $end = $6;
                    $strand = $7;
                    $trunc = $8;
                    
                    if ($strand eq "-")
                    {
                        my $temp = $start;
                        $start = $end;
                        $end = $temp;
                    }
                    if ($trunc eq "no")
                    {
                        $trunc = "";
                    }
                    
                    $record = [];
                    push(@$record, $pos);
                    push(@$record, $seq_name);
                    push(@$record, $start);
                    push(@$record, $end);
                    push(@$record, $strand);
                    push(@$record, $score);                    
                    push(@$record, $trunc);
                    push(@indexes, $record);
                }
            }
            $filepos = tell($fh);
        }        
        $self->close_file();            
    }
    
    @{$self->{indexes}} = sort sort_by_tRNAscanSE_output @indexes;
}

sub sort_by_tRNAscanSE_output
{
    if ((($a->[4] eq $b->[4]) && ($a->[4] eq "+")) ||
        ($a->[4] ne $b->[4]))
    {
        return ($a->[1] cmp $b->[1] ||
                $a->[4] cmp $b->[4] ||
                $a->[2] <=> $b->[2] ||
                $b->[5] <=> $a->[5]);
    }
    if (($a->[4] eq $b->[4]) && ($a->[4] eq "-"))
    {
        return ($a->[1] cmp $b->[1] ||
                $b->[4] cmp $a->[4] ||
                $b->[3] <=> $a->[3] ||
                $b->[5] <=> $a->[5]);        
    }    
}

sub sort_by_model
{
    if ((($a->[4] eq $b->[4]) && ($a->[4] eq "+")) ||
        ($a->[4] ne $b->[4]))
    {
        return ($a->[1] cmp $b->[1] ||
                $a->[0] cmp $b->[0] ||
                $b->[5] <=> $a->[5] ||
                $a->[4] cmp $b->[4] ||
                $a->[2] <=> $b->[2]);
    }
    if (($a->[4] eq $b->[4]) && ($a->[4] eq "-"))
    {
        return ($a->[1] cmp $b->[1] ||
                $a->[0] cmp $b->[0] ||
                $b->[5] <=> $a->[5] ||
                $b->[4] cmp $a->[4] ||
                $b->[3] <=> $a->[3]);        
    }    
}

sub get_cmsearch_record
{
    my $self = shift;
    my $filepos = shift;
    my $format = shift;
    my $fh = $self->{FILE_H};
    my $line = "";
    my $seq_name = "";
    my $score = 0;
    my $start = 0;
    my $end = 0;
    my $strand = "";
    my $seq_line = "";
    my $nc_line = "";
    my $ss = "";
    my $seq = "";
    my $model = "";
    my $nc = "";
    
    seek($fh, $filepos, 0);
    while ($line = <$fh>)
    {
        chomp($line);
        if (index($line, "Internal CM pipeline") > -1)
        {
            last;
        }        
        elsif ($line =~ /^>>\s*(\S+)/)
        {
            if ($seq_name ne "")
            {
                last;
            }            
            $seq_name = $1;
        }
        elsif ($line =~ /^\s+\(\d+\)\s+\S+\s+([e0-9.\-]+)\s+([0-9.\-]+)\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)\s+([+-])\s+\S+\s+\S+\s+(\S+)\s+([0-9.]+)/)
        {
            $score = $2;
            $start = $5;
            $end = $6;
            $strand = $7;
            
            if ($strand eq "-")
            {
                my $temp = $start;
                $start = $end;
                $end = $temp;
            }            
        }
        # Parse model structure line         
        elsif ($line =~ /^(.+) NC$/)
        {
            $nc_line = $1;
        }
        elsif ($line =~ /^\s{5,}([(),<>._\-,\[\]\{\}\:\~]{1,250}) CS$/)
        {
            $ss .= $1;

            # Parse model sequence line
            $line = <$fh>;
            if ($line =~ /^\s+\S+\s+\d+\s+([a-zA-Z\.0-9\>\<\[\]\*]{1,250})\s+\d+/)
            {  
                $model .= $1;  
            } 
            # Parse target sequence line
            $line = <$fh>;
            $line = <$fh>;
            $line =~ s/\*\[0\]\*/-----/g;
            if (($line =~ /^\s+\S+\s+\d+\s+([a-zA-Z\-]{1,250})\s+\d+/) ||
                ($line =~ /^\s+\S+\s+\-\s+([a-zA-Z\-]{1,250})\s+\-/))
            {
                $seq_line = $1;
                $seq .= $seq_line;

                # NC line
                $nc .= substr($nc_line, length($nc_line) - length($seq_line));
            }

            # Advance to next set of alignment info
            $line = <$fh>;
        }
        elsif ($line =~ /^(.+) PP$/)
        {
        }
        elsif ($line =~ /^(.+) RF$/)
        {
        }
    }

    if ($format)
    {
        ($ss, $seq) = $self->format_cmsearch_output($ss, $seq, $nc);
    }
    
    return ($ss, $seq, $model, $nc);
}

sub format_cmsearch_output
{    
    my $self = shift;
    my $cmsearch_ss = shift;
    my $cmsearch_seq = shift;
    my $cmsearch_nc = shift;
    $cmsearch_seq =~ s/U/T/g; 
    $cmsearch_seq =~ s/u/t/g;
    
    $cmsearch_ss = $self->fix_mismatch_ss($cmsearch_ss, $cmsearch_seq, $cmsearch_nc);
    
    for (my $index = 0; $index < length($cmsearch_seq); $index++)
    {
        if (substr($cmsearch_seq, $index, 1) eq '-')
        {
            substr($cmsearch_seq, $index, 1) = '*';
            if (length($cmsearch_ss) > $index)
            {
                substr($cmsearch_ss, $index, 1) = '*';
            }
        }
    }
    $cmsearch_seq =~ s/\*//g;
    $cmsearch_ss =~ s/\*//g;

    $cmsearch_ss =~ s/[,_\-:]/./g;
    $cmsearch_ss =~ s/[>)]/@/g;
    $cmsearch_ss =~ s/[(<]/>/g;
    $cmsearch_ss =~ s/@/</g;
    
    my $diff = length($cmsearch_seq) - length($cmsearch_ss);
    for (my $ct = 0; $ct < $diff; $ct++)
    {
        $cmsearch_ss .= ".";
    }
    
    return ($cmsearch_ss, $cmsearch_seq);
}

sub fix_mismatch_ss
{
    my $self = shift;
    my ($cmsearch_ss, $cmsearch_seq, $cmsearch_nc) = @_;
    
    my @left = ();
    my @right = ();
    my %pairs = ();
    my $left_index = -1;
    
    for (my $index = 0; $index < length($cmsearch_nc); $index++)
    {
        if (substr($cmsearch_nc, $index, 1) eq "v")
        {
            substr($cmsearch_ss, $index, 1) = ".";
        }
    }

    for (my $pos = 0; $pos < length($cmsearch_ss); $pos++)
	{
		if (substr($cmsearch_ss, $pos, 1) eq "<" or substr($cmsearch_ss, $pos, 1) eq "(")
		{
			push(@left, $pos);
			$pairs{$#left} = -1;
		}
		elsif (substr($cmsearch_ss, $pos, 1) eq ">" or substr($cmsearch_ss, $pos, 1) eq ")")
		{
			push(@right, $pos);
			$left_index = scalar(@left) - 1;
			while (($pairs{$left_index} > -1) && ($left_index > -1))
			{
				$left_index--;
			}
			if (($left_index > -1) && ($pairs{$left_index} == -1))
			{
				$pairs{$left_index} = scalar(@right) - 1;
			}
		}
	}

    foreach $left_index (sort keys %pairs)
    {
        my $left_base = uc(substr($cmsearch_seq, $left[$left_index], 1));
        my $right_base = uc(substr($cmsearch_seq, $right[$pairs{$left_index}], 1));
        
        if (($left_base eq "A" and $right_base ne "U" and $right_base ne "T") ||
            ($left_base eq "T" and $right_base ne "A" and $right_base ne "G") ||
            ($left_base eq "U" and $right_base ne "A" and $right_base ne "G") ||
            ($left_base eq "G" and $right_base ne "C" and $right_base ne "U" and $right_base ne "T") ||
            ($left_base eq "C" and $right_base ne "G") ||
            ($left_base eq "-") ||
            ($right_base eq "-"))
        {
            substr($cmsearch_ss, $left[$left_index], 1) = ".";
            substr($cmsearch_ss, $right[$pairs{$left_index}], 1) = ".";
        }
    }
    
    return $cmsearch_ss;
}

sub get_next_tab_seq_hits
{
    my $self = shift;
    my $line = "";
    my $filepos = undef;
    my $record = [];
    my @hits = ();
    my $seq_name = "";
    my $score = 0;
    my $start = 0;
    my $end = 0;
    my $strand = "";
    my @columns = ();
    my $last_seq = "";
    
    my $fh = $self->{FILE_H};
    $filepos = tell($fh);
    while ($line = <$fh>)
    {
        chomp($line);
        if ($line !~ /^#/)
        {
            @columns = split(/\s+/, $line);

            if ($columns[3] ne $last_seq)
            {
                if (scalar(@hits) > 0)
                {
                    seek($fh, $filepos, 0);
                    last;
                }
                $last_seq = $columns[3];
            }
            
            $start = $columns[9];
            $end = $columns[10];
            $strand = $columns[11];
            
            if ($strand eq "-")
            {
                my $temp = $start;
                $start = $end;
                $end = $temp;
            }
            
            $record = [];
            push(@$record, $columns[1]);
            push(@$record, $columns[3]);
            push(@$record, $start);
            push(@$record, $end);
            push(@$record, $strand);
            push(@$record, $columns[16]);
            push(@hits, $record);
        }
        $filepos = tell($fh);
    }
    
    my @sorted_hits = sort sort_by_model @hits;
    return \@sorted_hits;
}

1;
