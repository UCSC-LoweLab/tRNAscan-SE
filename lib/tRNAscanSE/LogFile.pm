# tRNAscanSE/LogFile.pm
# This class defines the log file used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2017 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::LogFile;

use strict;
use POSIX 'strftime';

sub new {
    my $class = shift;
    my $log_name = shift;
    my $self = {
        file_name => $log_name,
        FILE_H => undef,
        log_file => 0
    };
    if ($log_name ne "")
    {
        $self->{log_file} = 1;
    }

    $self->{process_id} = $$;
    $self->{quiet_mode} = 0;
    
    bless ($self, $class);
    return $self;
}

sub DESTROY
{
    my $self = shift;
}

sub quiet_mode
{
    my $self = shift;
    if (@_) { $self->{quiet_mode} = shift; }
    return $self->{quiet_mode};
}

sub file_name
{
    my $self = shift;
    if (@_) { $self->{file_name} = shift; }
    return $self->{file_name};
}

sub open_file
{
    my $self = shift;
    my $file = shift;
    
    my $success = 0;

    if (($file eq "-") || ($file eq "/dev/null"))
    {
        $success = open($self->{FILE_H}, ">$file");
        $self->{file_name} = $file;
    }
    else
    {
        open($self->{FILE_H}, ">$file") || die "Unable to open $file for writing. Aborting program.\n";
        $|=1;
        $self->{file_name} = $file;
        $success = 1;
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

sub initialize_log
{
    my $self = shift;
    my $app = shift;
    my $log_dir = shift;
    my $cmd = shift;

    if ($self->{log_file})
    {
        my $time_str = strftime('%Y%m%d-%H:%M:%S', localtime);
        my $file = $log_dir."/$app-$time_str-".$self->{process_id}.".log";
        if ($self->{file_name} ne "default")
        {
            $file = $self->{file_name};
        }
        if (open_file($self, $file))
        {
            write_line($self, "Application: $app");
            write_line($self, "User: ".getlogin);
            write_line($self, "Host: ".`hostname`);
            write_line($self, "Start Time: ".localtime);
            write_line($self, "");
            if ($cmd ne "")
            {
                write_line($self, "Command: $cmd");
                write_line($self, "");
            }
        }
    }
}

sub finish_process()
{
    my $self = shift;

    if ($self->{log_file})
    {
        write_line($self, "");
        write_line($self, "End Time: ".localtime);
        close_file($self);
    }
}

sub command
{
    my $self = shift;
    my $line = shift;    

    if ($self->{log_file})
    {
        write_line($self, "Command: ".$line);
    }
}

sub broadcast
{
    my $self = shift;
    my $line = shift;    

    if ($self->{log_file})
    {
        write_line($self, $line);
    }
}

sub status
{
    my $self = shift;
    my $line = shift;    
    if ($self->{log_file})
    {
        write_line($self, "Status: ".$line);
    }
    if (!$self->{quiet_mode})
    {
        print "Status: ".$line."\n";
    }
}

sub error
{
    my $self = shift;
    my $line = shift;    

    if ($self->{log_file})
    {
        write_line($self, "Error: ".$line);
    }
    if (!$self->{quiet_mode})
    {
        print "Error: ".$line."\n";
    }
}

sub warning
{
    my $self = shift;
    my $line = shift;    

    if ($self->{log_file})
    {
        write_line($self, "Warning: ".$line);
    }
    if (!$self->{quiet_mode})
    {
        print "Warning: ".$line."\n";
    }
}

sub debug
{
    my $self = shift;
    my $line = shift;    

    if ($self->{log_file})
    {
        write_line($self, "Debug: ".$line);
    }
}

1;
