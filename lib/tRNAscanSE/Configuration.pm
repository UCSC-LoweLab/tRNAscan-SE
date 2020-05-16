# tRNAscanSE/Configuration.pm
# This class contains functions for parsing tRNAscan-SE configuration file.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2017 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::Configuration;

use strict;

sub new {
    my $class = shift;
    my $self = { @_ };

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
	
	$self->{conf} = "";
	$self->{home_dir} = $ENV{'HOME'};
}

sub config_file
{
    my $self = shift;
    if (@_) { $self->{conf} = shift; }
    return $self->{conf};
}

sub get
{
    my $self = shift;
	my $key = shift;
	
	if (defined $self->{$key})
	{
		return $self->{$key};
	}
	else
	{
		return undef;
	}
}

sub get_subvalue
{
    my $self = shift;
	my $key = shift;
	my $subkey = shift;
	
	if (defined $self->{$key}->{$subkey})
	{
		return $self->{$key}->{$subkey};
	}
	else
	{
		return undef;
	}
}

# Read configuration file
sub read_configuration_file
{
	my $self = shift;
	my $line = "";
    my ($key, $value, $temp, $key2, $value2, $subkey1, $subkey2, $index);

	open (FILE_IN, "$self->{conf}") || die "Fail to open configuration file $self->{conf}. Please use -c to specify a local configuration file.";
    while ($line = <FILE_IN>)
    {
        chomp($line);
        if ($line !~ /^#/)
        {
            if ($line =~ /^(\S+):\s+(.+)\s*$/)
            {
                $key = $1;
                $value = $2;
                $temp = $value;
                $index = index($temp, "}");
                while ($index > -1)
                {
                    $key2 = substr($temp, index($temp,"{") + 1, $index - index($temp,"{") - 1);
                    $value2 = $self->{$key2};
                    $value =~ s/{$key2}/$value2/;
                    $temp = substr($temp, $index + 1);
                    $index = index($temp, "}");
                }
                
				if (index($value, "\$\$") > -1)
				{
					$temp = "$$";
					$value =~ s/\$\$/$temp/;
				}
				
                if ($key =~ /^(\S+)\.(\S+)/)
                {
                    $subkey1 = $1;
                    $subkey2 = $2;
                    $self->{$subkey1}->{$subkey2} = $value;
                }
                else
                {
					if ($key eq "temp_dir" and defined $self->{$key}) {}
					else
					{
						$self->{$key} = $value;
					}
                }
            }
        }
    }
    close (FILE_IN);
}

sub set_temp_dir
{
	my $self = shift;
    if (@_) { $self->{temp_dir} = shift; }	
}

1;
