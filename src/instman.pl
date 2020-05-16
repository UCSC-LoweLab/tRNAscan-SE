#! /usr/bin/perl

# Sean Eddy, Wed Jul 29 15:24:43 1992

# instman - "cp" replacement for formatting and installing man pages
#           Requires that GNU groff is installed: geqn, gtbl, groff.
#
# Usage:  instman <man-formatted file> <destination for formatted file>
#


$usage = "Usage: instman <infile> <outfile>\n  where <infile> is a man-formatted input file,\n  and <outfile> is the name to store the formatted file under\n  (such as /usr/catman/local/cat1/foo)\n";

if ($#ARGV != 1) { die "Incorrect argument number.\n$usage"; } 
$man  = shift(@ARGV);	
$dest = shift(@ARGV);
print "Reading $man, installing as $dest\n";

die "Error: can't read $man.\n$usage" unless -r $man;
die "Error: can't write $dest.\n$usage" unless -w $man;
die "Error: $dest is a directory, not a file name\n$usage" if -d $dest;

system "geqn $man | gtbl | groff -man -Tascii > $dest";
