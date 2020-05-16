# tRNAscanSE/FpScanResultFile.pm
# This class defines the first pass scanning result file used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2017 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::FpScanResultFile;

use strict;
use tRNAscanSE::Utils;
use tRNAscanSE::Sequence;
use tRNAscanSE::tRNA;
use tRNAscanSE::ArraytRNA;
use tRNAscanSE::Options;

sub new
{
    my $class = shift;
    my $file_name = shift;
    my $self = {
        file_name => $file_name,
        FILE_H => undef,
        FLANKING_H => undef
    };
    
    $self->{indexes} = [];
    $self->{flanking_file} = "";
    
    bless ($self, $class);
    $self->reset_current_seq();
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
    $self->reset_current_seq();
}

sub file_name
{
    my $self = shift;
    if (@_) { $self->{file_name} = shift; }
    return $self->{file_name};
}

sub flanking_file
{
    my $self = shift;
    if (@_) { $self->{flanking_file} = shift; }
    return $self->{flanking_file};
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

sub get_hit_count
{
    my $self = shift;
    my $count = 0;
    
    for (my $i = 0; $i < scalar(@{$self->{indexes}}); $i++)
    {
        $count += $self->{indexes}->[$i]->[2];
    }
    return $count;
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

sub reset_current_seq
{
    my $self = shift;
    $self->{current_seq} = -1;
    $self->{current_record} = 0;    
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
        die "First pass result file name is not set.\n"
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

sub init_fp_result_file
{
    my $self = shift;
    my ($file) = @_;   
    $self->{file_name} = $file;
    
    &open_for_append(\$self->{FILE_H}, $self->{file_name});
    my $fh = $self->{FILE_H};
    
    print $fh "Sequence\t\ttRNA Bounds\ttRNA\tAnti\t\n";
	print $fh "Name     \ttRNA #\tBegin\tEnd\tType\tCodon\tSeqID\tSeqLen\tScore\n";
	print $fh "--------\t------\t-----\t---\t----\t-----\t-----\t------\t-----\n";
    
    $self->close_file();
}

sub save_firstpass_output
{
    my $self = shift;
    my ($opts, $fp_tRNAs, $r_fpass_trna_base_ct, $seq_len, $seq_id) = @_;
    my $triplet = "";

    my @source_tab = ('Inf', 'Ts', 'Eu', 'Bo');
    
    if (!$opts->CM_mode())
	{
		&open_for_append(\$self->{FILE_H}, $opts->out_file());
    }
    else
	{		       
		&open_for_append(\$self->{FILE_H}, $self->{file_name});	
    }
    my $fh = $self->{FILE_H};
    
	for (my $i = 0; $i < $fp_tRNAs->get_count(); $i++)
	{
		my $trna = $fp_tRNAs->get($i);

		$triplet = uc($trna->anticodon());
		if ($opts->output_codon())
		{
			$triplet = &rev_comp_seq($triplet);
		}
		
		printf $fh "%-10s\t%d\t%d\t%d\t%s\t%s\t",
			$trna->seqname(), $i + 1, $trna->start(), $trna->end(), $trna->isotype(), $triplet;
		
		# save intron bounds if not doing covariance model analysis		
		if (!$opts->CM_mode())
		{
			if ($trna->get_intron_count() > 0)
			{
				my @ar_introns = $trna->ar_introns();
				printf $fh "%d\t%d\t", $ar_introns[0]->{rel_start}, $ar_introns[0]->{rel_end};				
			}
			else
			{
				printf $fh "0\t0\t";
			}
			printf $fh "%.2f", $trna->score();
		}	
		# save seq id number and source seq length if needed for covariance model analysis 	
		else
		{
			printf $fh "%d\t%d\t%.2f", $seq_id, $seq_len, $trna->score();
			if ($trna->model() ne "")
			{
				printf $fh "\t%s", $trna->model();
			}
			else
			{
				printf $fh "\tnone";
			}
		}
		
		if ($opts->save_source())
		{
			print $fh "\t", $source_tab[$trna->hit_source()];
		}
		print $fh "\n";
		
		$$r_fpass_trna_base_ct += abs($trna->end() - $trna->start()) + 1;
    }
    $self->close_file();
}				

# Create dummy first-pass result file with all sequences
sub prep_for_secpass_only
{
    my $self = shift;
    my ($opts, $stats, $seq_file) = @_;

    $seq_file->open_file($opts->fasta_file(), "read");
    &open_for_append(\$self->{FILE_H}, $self->{file_name});
    my $fh = $self->{FILE_H};
    
    # Don't look for a specific Seq number 
    while ($seq_file->read_fasta($opts, 0))
	{
		print $fh $seq_file->seq_name()."\t1\t1\t".$seq_file->seq_length()."\t???\t???\t".$seq_file->seq_id()."\t".$seq_file->seq_length()." C none\n";
		$stats->increment_numscanned();
    }
    
    $self->close_file();
    $seq_file->close_file();
}

sub index_results
{
    my $self = shift;
    my ($r_seqinfo_flag) = @_;
    my $line = "";
    my $filepos = undef;
    my $line_ct = 0;
    my $prev_seqname = "";
    my $seqname = "";
    my $record = [];
    
    if ($self->open_file())
    {
        my $fh = $self->{FILE_H};
        $filepos = tell($fh);
        while ($line = <$fh>) 
        {	
            if ($line =~ /Type\tCodon\tSeqID\tSeqLen/)
            {
                # Column header present if we record seqID's and lengths
                $$r_seqinfo_flag = 1;
            }
            elsif ($line =~ /^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)/o)  
            {
                $seqname = $1;
                if ($seqname ne $prev_seqname)
                {
                    $record = [];
                    push(@$record, $filepos);
                    push(@$record, $seqname);
                    push(@$record, 0);
                    push(@{$self->{indexes}}, $record);
                    
                    if ($prev_seqname ne "")
                    {
                        $self->{indexes}->[scalar(@{$self->{indexes}})-2]->[2] = $line_ct;
                    }
                    
                    $prev_seqname = $seqname;
                    $line_ct = 0;
                }
                $line_ct++;
            }
            $filepos = tell($fh);
        }
        if ($prev_seqname ne "")
        {
            $self->{indexes}->[scalar(@{$self->{indexes}})-1]->[2] = $line_ct;
        }
        $self->close_file();
    }
}

sub get_next_tRNA_candidate
{
    my $self = shift;
    my ($opts, $seqinfo_flag, $seq_ct, $trna) = @_;
    my $fh = $self->{FILE_H};
    my $line = "";
    my $padding = $opts->padding();
    my ($trnact, $hit_source);

    $trna->clear();
    my $record = $self->{indexes}->[$self->{current_seq}];
    if ($self->{current_seq} != $seq_ct)
    {
        $self->{current_seq} = $seq_ct;
        $self->{current_record} = 0;
        $record = $self->{indexes}->[$seq_ct];
        seek($fh, $record->[0], 0);
    }
    
    if (defined $fh and $line = <$fh> and ($self->{current_record} < $record->[2])) 
    {
		if ($line =~ /^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)/o)  
		{
			$trna->seqname($1);
			$trnact = $2;
            $trna->id($trnact);
			$trna->tRNAscan_id($trna->seqname().".t".$trnact);
			$trna->start($3);	        # trna subseq absolute start index
			$trna->end($4);				# trna subseq absolute end index
			$trna->isotype($5);
			$trna->anticodon($6);
			$trna->src_seqid($7);
			$trna->src_seqlen($8);
			$trna->score($9);
			$trna->model($10);
			$hit_source = $';
			$hit_source =~ s/[\s\t\n]//g;
			$trna->hit_source($hit_source);
			
			# if seqinfo_flag not set, file does not have SeqID info in
			#  7th column of output, don't mistake number read for SeqID			
			if (!$seqinfo_flag)
			{
				$trna->src_seqid(0);
			}
			
			if ($trna->end() > $trna->start())  
			{
				$trna->strand("+");   
			
				# pad ends of sequence only if EufindtRNA or Infernal first-pass is being used
				#  and $seqinfo_flag is set (we know the seq lengths)	
				if (($opts->eufind_mode() || $opts->infernal_fp()) && $seqinfo_flag) 
				{
					$trna->start(&max(1, $trna->start() - $padding));
					$trna->end(&min($trna->src_seqlen(), $trna->end() + $padding));
				}
			}
			else
			{ 
				$trna->strand("-");
				
				if (($opts->eufind_mode() || $opts->infernal_fp()) && $seqinfo_flag)
				{
					$trna->start(&min($trna->src_seqlen(), $trna->start() + $padding));
					$trna->end(&max(1, $trna->end() - $padding));
				}
		    }
	    
			if ($trna->end() == $trna->start())
			{
				print STDERR "Error reading $self->{file_name}: tRNA of length 0\n"; 
			}
			
		}
        $self->{current_record}++;
    }
}

sub open_flanking
{
    my $self = shift;
    my $mode = shift;
    
    my $success = 0;
    
    if ($self->{flanking_file} ne "")
    {
        if ($mode eq "write")
        {
            &open_for_write(\$self->{FLANKING_H}, $self->{flanking_file});
        }
        else
        {
            &open_for_read(\$self->{FLANKING_H}, $self->{flanking_file});
        }
        $success = 1;
    }
    else
    {
        die "First pass flanking file name is not set.\n"
    }

    return $success;
}

sub close_flanking
{
    my $self = shift;
    
    if (defined $self->{FLANKING_H})
    {
        close($self->{FLANKING_H});
    }
}

sub write_tRNA_flanking
{
    my $self = shift;
    my $trna = shift;
    
    if (defined $self->{FLANKING_H})
    {
        my $fh = $self->{FLANKING_H};
        print $fh $trna->seqname()."\t".$trna->id()."\t".$trna->seq()."\t".$trna->upstream()."\t".$trna->downstream()."\n";
    }
}

sub read_tRNA_flanking
{
    my $self = shift;
    my $trna = shift;
    
    my $fh = $self->{FLANKING_H};
    my $line = "";
    
    if (defined $fh and $line = <$fh>) 
    {
        chomp($line);
        my @columns = split(/\t/, $line);
        if ($trna->seqname() eq $columns[0] and $trna->id() == $columns[1])
        {
            $trna->seq($columns[2]);
            $trna->upstream($columns[3]);
            $trna->downstream($columns[4]);
        }
    }
}

1;
