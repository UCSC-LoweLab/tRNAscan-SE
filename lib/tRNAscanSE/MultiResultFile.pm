# tRNAscanSE/MultiResultFile.pm
# This class defines the intermediate multi-model result file used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2017 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::MultiResultFile;

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
        $self->{file_name} = "tscan$$"."_multisp.out";
    }
    
    $self->{models} = [];
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
        }
        elsif ($mode eq "append")
        {
            &open_for_append(\$self->{FILE_H}, $self->{file_name});        
        }
        $success = 1;
    }
    else
    {
        die "Multi-model intermediate result file name is not set.\n"
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

sub get_line
{
    my $self = shift;
    my $fh = $self->{FILE_H};
    if (my $line = <$fh>)
    {
        chomp($line);
        return ($line);
    }
    else
    {
        return undef;
    }
}

sub read_models
{
    my $self = shift;
    my $line = $self->get_line();
    if (defined $line)
    {
        my @columns = split(/\t/, $line);
        for (my $i = 1; $i < scalar(@columns); $i++)
        {
            push(@{$self->{models}}, $columns[$i]);
        }
    }    
}

sub get_next_record
{
    my $self = shift;
    my $fh = $self->{FILE_H};
    if (defined $fh and my $line = <$fh>)
    {
        chomp($line);
        my @columns = split(/\t/, $line);
        return ($columns[0], $line);
    }
    else
    {
        return ("", "");
    }
}

sub get_next_tRNA
{
    my $self = shift;
    my $tRNA = shift;
    
    $self->retrieve_record($tRNA);
}

sub get_tRNA
{
    my $self = shift;    
    my $filepos = shift;
    my $tRNA = shift;

    my $fh = $self->{FILE_H};
    
    seek($fh, $filepos, 0);
    $self->retrieve_record($tRNA);
}

sub retrieve_record
{
    my $self = shift;    
    my $tRNA = shift;

    my $fh = $self->{FILE_H};
    
    if (defined $fh and my $line = <$fh>)
    {
        chomp($line);
        my @columns = split(/\t/, $line, -1);
        if (defined $columns[0] and $columns[0] eq $tRNA->seqname().".t".&pad_num($tRNA->id(), 6))
        {
            for (my $i = 1; $i < scalar(@columns); $i++)
            {
                if ($columns[$i] ne "")
                {
                    if ($self->{models}->[$i-1] =~ /^mito_(\S+)$/)
                    {
                        $tRNA->add_model_hit("mito", $1, $columns[$i], "");
                    }
                    else
                    {
                        $tRNA->add_model_hit("cyto", $self->{models}->[$i-1], $columns[$i], "");
                    }
                }
            }
        }
    }    
}

sub sort_records
{
    my $self = shift;
    my $key = shift;
    
    if ($key eq "tRNAscan_id")
    {
        $self->build_tRNAscan_id_index();
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
            push(@$record, $columns[0]);
            push(@indexes, $record);
            $filepos = tell($fh);
        }        
        $self->close_file();            
    }
    
    @{$self->{indexes}} = sort {$a->[1] cmp $b->[1]} @indexes;
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
