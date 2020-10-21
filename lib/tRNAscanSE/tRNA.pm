# tRNAscanSE/tRNA.pm
# This class contains parameters and functions describing a tRNA gene used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2017 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::tRNA;

use strict;
use tRNAscanSE::SprinzlPos;

use constant TYPE2_TRNA =>
{
	Eukaryota => "Leu;Ser",
	Archaea => "Leu;Ser",
	Bacteria => "Leu;Ser;Tyr"
};

use constant
{
	_seqname => 0,
	_ordered_seqname => 1,
	_start => 2,
	_end => 3,
	_strand => 4,
	_start2 => 5,
	_end2 => 6,
	_strand2 => 7,
	_start3 => 8,
	_end3 => 9,
	_strand3 => 10,
	_position => 11,
	_id => 12,
	_tRNAscan_id => 13,
	_gtrnadb_id => 14,
	_extdb_id => 15,
	_isotype => 16,
	_anticodon => 17,
	_ar_ac_pos => 18,
	_ar_introns => 19,
	_model => 20,
	_clade => 21,
	_score => 22,
	_mat_score => 23,
	_hmm_score => 24,
	_ss_score => 25,
	_h_domain_models => 26,
	_seq => 27,
	_mat_seq => 28,
	_mat_ss => 29,
	_ss => 30,
	_sprinzl_align => 31,
	_sprinzl_ss => 32,
	_pseudo => 33,
	_trunc => 34,
	_category => 35,
	_upstream => 36,
	_downstream => 37,
	_hit_source => 38,
	_src_seqid => 39,
	_src_seqlen => 40,
	_ar_multi_models => 41,
	_ar_pos_sprinzl_map => 42,
	_h_sprinzl_pos_map => 43,
	_ar_non_canonical => 44,
	_note => 45
};

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

	$self->{data} = [];

	$self->{data}->[_seqname] = "";
	$self->{data}->[_ordered_seqname] = 0;
	$self->{data}->[_start] = 0;
	$self->{data}->[_end] = 0;
	$self->{data}->[_strand] = "";
	$self->{data}->[_start2] = 0;
	$self->{data}->[_end2] = 0;
	$self->{data}->[_strand2] = "";
	$self->{data}->[_start3] = 0;
	$self->{data}->[_end3] = 0;
	$self->{data}->[_strand3] = "";
	$self->{data}->[_position] = 0;
	$self->{data}->[_id] = 0;
	$self->{data}->[_tRNAscan_id] = "";
	$self->{data}->[_gtrnadb_id] = "";
	$self->{data}->[_extdb_id] = "";
	$self->{data}->[_isotype] = "";
	$self->{data}->[_anticodon] = "";
	$self->{data}->[_ar_ac_pos] = [];
	$self->{data}->[_ar_introns] = [];
	$self->{data}->[_model] = "";
	$self->{data}->[_clade] = "";
	$self->{data}->[_score] = 0;
	$self->{data}->[_mat_score] = 0;
	$self->{data}->[_hmm_score] = 0;
	$self->{data}->[_ss_score] = 0;
	$self->{data}->[_h_domain_models] = {};
	$self->{data}->[_seq] = "";
	$self->{data}->[_mat_seq] = "";
	$self->{data}->[_mat_ss] = "";
	$self->{data}->[_ss] = "";
	$self->{data}->[_sprinzl_align] = "";
	$self->{data}->[_sprinzl_ss] = "";
	$self->{data}->[_pseudo] = 0;
	$self->{data}->[_trunc] = "";
	$self->{data}->[_category] = "";
	$self->{data}->[_upstream] = "";
	$self->{data}->[_downstream] = "";
	$self->{data}->[_hit_source] = "";
	$self->{data}->[_src_seqid] = 0;
	$self->{data}->[_src_seqlen] = 0;
	$self->{data}->[_ar_multi_models] = [];
	$self->{data}->[_ar_pos_sprinzl_map] = [];
	$self->{data}->[_h_sprinzl_pos_map] = {};
	$self->{data}->[_ar_non_canonical] = [];
	$self->{data}->[_note] = "";
}

sub clear
{
    my $self = shift;

	@{$self->{data}} = ();
	
	$self->{data}->[_seqname] = "";
	$self->{data}->[_ordered_seqname] = 0;
	$self->{data}->[_start] = 0;
	$self->{data}->[_end] = 0;
	$self->{data}->[_strand] = "";
	$self->{data}->[_start2] = 0;
	$self->{data}->[_end2] = 0;
	$self->{data}->[_strand2] = "";
	$self->{data}->[_start3] = 0;
	$self->{data}->[_end3] = 0;
	$self->{data}->[_strand3] = "";
	$self->{data}->[_position] = 0;
	$self->{data}->[_id] = 0;
	$self->{data}->[_tRNAscan_id] = "";
	$self->{data}->[_gtrnadb_id] = "";
	$self->{data}->[_extdb_id] = "";
	$self->{data}->[_isotype] = "";
	$self->{data}->[_anticodon] = "";
	@{$self->{data}->[_ar_ac_pos]} = ();
	@{$self->{data}->[_ar_introns]} = ();
	$self->{data}->[_model] = "";
	$self->{data}->[_clade] = "";
	$self->{data}->[_score] = 0;
	$self->{data}->[_mat_score] = 0;
	$self->{data}->[_hmm_score] = 0;
	$self->{data}->[_ss_score] = 0;
	%{$self->{data}->[_h_domain_models]} = ();
	$self->{data}->[_seq] = "";
	$self->{data}->[_mat_seq] = "";
	$self->{data}->[_mat_ss] = "";
	$self->{data}->[_ss] = "";
	$self->{data}->[_sprinzl_align] = "";
	$self->{data}->[_sprinzl_ss] = "";
	$self->{data}->[_pseudo] = 0;
	$self->{data}->[_trunc] = "";
	$self->{data}->[_category] = "";
	$self->{data}->[_upstream] = "";
	$self->{data}->[_downstream] = "";
	$self->{data}->[_hit_source] = "";
	$self->{data}->[_src_seqid] = 0;
	$self->{data}->[_src_seqlen] = 0;
	@{$self->{data}->[_ar_multi_models]} = ();
	@{$self->{data}->[_ar_pos_sprinzl_map]} = ();
	%{$self->{data}->[_h_sprinzl_pos_map]} = ();
	@{$self->{data}->[_ar_non_canonical]} = ();
	$self->{data}->[_note] = "";
}

sub copy
{
	my $self = shift;
	my $obj = shift;
	
	$self->{data}->[_seqname] = $obj->{data}->[_seqname];
	$self->{data}->[_ordered_seqname] = $obj->{data}->[_ordered_seqname];
	$self->{data}->[_start] = $obj->{data}->[_start];
	$self->{data}->[_end] = $obj->{data}->[_end];
	$self->{data}->[_strand] = $obj->{data}->[_strand];
	$self->{data}->[_start2] = $obj->{data}->[_start2];
	$self->{data}->[_end2] = $obj->{data}->[_end2];
	$self->{data}->[_strand2] = $obj->{data}->[_strand2];
	$self->{data}->[_start3] = $obj->{data}->[_start3];
	$self->{data}->[_end3] = $obj->{data}->[_end3];
	$self->{data}->[_strand3] = $obj->{data}->[_strand3];
	$self->{data}->[_position] = $obj->{data}->[_position];
	$self->{data}->[_id] = $obj->{data}->[_id];
	$self->{data}->[_tRNAscan_id] = $obj->{data}->[_tRNAscan_id];
	$self->{data}->[_gtrnadb_id] = $obj->{data}->[_gtrnadb_id];
	$self->{data}->[_extdb_id] = $obj->{data}->[_extdb_id];
	$self->{data}->[_isotype] = $obj->{data}->[_isotype];
	$self->{data}->[_anticodon] = $obj->{data}->[_anticodon];
	@{$self->{data}->[_ar_ac_pos]} = @{$obj->{data}->[_ar_ac_pos]};
	@{$self->{data}->[_ar_introns]} = @{$obj->{data}->[_ar_introns]};
	$self->{data}->[_model] = $obj->{data}->[_model];
	$self->{data}->[_clade] = $obj->{data}->[_clade];
	$self->{data}->[_score] = $obj->{data}->[_score];
	$self->{data}->[_mat_score] = $obj->{data}->[_mat_score];
	$self->{data}->[_hmm_score] = $obj->{data}->[_hmm_score];
	$self->{data}->[_ss_score] = $obj->{data}->[_ss_score];
	%{$self->{data}->[_h_domain_models]} = %{$obj->{data}->[_h_domain_models]};
	$self->{data}->[_seq] = $obj->{data}->[_seq];
	$self->{data}->[_mat_seq] = $obj->{data}->[_mat_seq];
	$self->{data}->[_mat_ss] = $obj->{data}->[_mat_ss];
	$self->{data}->[_ss] = $obj->{data}->[_ss];
	$self->{data}->[_sprinzl_align] = $obj->{data}->[_sprinzl_align];
	$self->{data}->[_sprinzl_ss] = $obj->{data}->[_sprinzl_ss];
	$self->{data}->[_pseudo] = $obj->{data}->[_pseudo];
	$self->{data}->[_trunc] = $obj->{data}->[_trunc];
	$self->{data}->[_category] = $obj->{data}->[_category];
	$self->{data}->[_upstream] = $obj->{data}->[_upstream];
	$self->{data}->[_downstream] = $obj->{data}->[_downstream];
	$self->{data}->[_hit_source] = $obj->{data}->[_hit_source];
	$self->{data}->[_src_seqid] = $obj->{data}->[_src_seqid];
	$self->{data}->[_src_seqlen] = $obj->{data}->[_src_seqlen];
	@{$self->{data}->[_ar_multi_models]} = @{$obj->{data}->[_ar_multi_models]};
	@{$self->{data}->[_ar_pos_sprinzl_map]} = @{$obj->{data}->[_ar_pos_sprinzl_map]};
	%{$self->{data}->[_h_sprinzl_pos_map]} = %{$obj->{data}->[_h_sprinzl_pos_map]};
	@{$self->{data}->[_ar_non_canonical]} = @{$obj->{data}->[_ar_non_canonical]};
	$self->{data}->[_note] = $obj->{data}->[_note];
}

sub set_covels_hit
{
	my $self = shift;
	my ($covels_hit) = @_;
	$self->{data}->[_seqname] = $covels_hit->{hit_seqname};
	$self->set_domain_model("cove", $covels_hit->{score});
	if ($covels_hit->{tRNA_start} < $covels_hit->{tRNA_end})
	{
		$self->{data}->[_strand] = "+";
		$self->{data}->[_start] = $covels_hit->{tRNA_start};
		$self->{data}->[_end] = $covels_hit->{tRNA_end};
	}
	else
	{
		$self->{data}->[_strand] = "-";
		$self->{data}->[_start] = $covels_hit->{tRNA_end};
		$self->{data}->[_end] = $covels_hit->{tRNA_start};		
	}
}

sub set_cmsearch_hit
{
	my $self = shift;
	my ($cmsearch_hit) = @_;
	$self->{data}->[_seqname] = $cmsearch_hit->{hit_seqname};
	$self->set_domain_model("infernal", $cmsearch_hit->{score});
	$self->{data}->[_strand] = $cmsearch_hit->{strand};
	if ($cmsearch_hit->{strand} eq "+")
	{
		$self->{data}->[_start] = $cmsearch_hit->{tRNA_start};
		$self->{data}->[_end] = $cmsearch_hit->{tRNA_end};
	}
	else
	{
		$self->{data}->[_start] = $cmsearch_hit->{tRNA_end};
		$self->{data}->[_end] = $cmsearch_hit->{tRNA_start};		
	}
	$self->{data}->[_seq] = $cmsearch_hit->{seq};
	$self->{data}->[_ss] = $cmsearch_hit->{ss};
}

sub seqname
{
    my $self = shift;
    if (@_) { $self->{data}->[_seqname] = shift; }
    return $self->{data}->[_seqname];
}

sub ordered_seqname
{
    my $self = shift;
    if (@_) { $self->{data}->[_ordered_seqname] = shift; }
    return $self->{data}->[_ordered_seqname];
}

sub start
{
    my $self = shift;
    if (@_) { $self->{data}->[_start] = shift; }
    return $self->{data}->[_start];
}

sub end
{
    my $self = shift;
    if (@_) { $self->{data}->[_end] = shift; }
    return $self->{data}->[_end];
}

sub strand
{
    my $self = shift;
    if (@_) { $self->{data}->[_strand] = shift; }
    return $self->{data}->[_strand];
}

sub exon_start
{
    my $self = shift;
	my $exon = shift;
	my $start = 0;
    if (@_)
	{
		$start = shift;
		$self->{data}->[_start2] = $start if ($exon == 2);
		$self->{data}->[_start3] = $start if ($exon == 3);
	}
	if ($exon == 2)
	{
		return $self->{data}->[_start2];
	}
	if ($exon == 3)
	{
		return $self->{data}->[_start3];
	}	
}

sub exon_end
{
    my $self = shift;
	my $exon = shift;
	my $end = 0;
    if (@_)
	{
		$end = shift;
		$self->{data}->[_end2] = $end if ($exon == 2);
		$self->{data}->[_end3] = $end if ($exon == 3);
	}
	if ($exon == 2)
	{
		return $self->{data}->[_end2];
	}
	if ($exon == 3)
	{
		return $self->{data}->[_end3];
	}	
}

sub exon_strand
{
    my $self = shift;
	my $exon = shift;
	my $strand = "";
    if (@_)
	{
		$strand = shift;
		$self->{data}->[_strand2] = $strand if ($exon == 2);
		$self->{data}->[_strand3] = $strand if ($exon == 3);
	}
	if ($exon == 2)
	{
		return $self->{data}->[_strand2];
	}
	if ($exon == 3)
	{
		return $self->{data}->[_strand3];
	}	
}

sub position
{
    my $self = shift;
    if (@_) { $self->{data}->[_position] = shift; }
    return $self->{data}->[_position];
}

sub len
{
    my $self = shift;
	my $len = 0;
	if ($self->{data}->[_start] < $self->{data}->[_end])
	{
		$len = $self->{data}->[_end] - $self->{data}->[_start] + 1;
	}
	else
	{
		$len = $self->{data}->[_start] - $self->{data}->[_end] + 1;
	}
    return $len;
}

sub id
{
    my $self = shift;
    if (@_) { $self->{data}->[_id] = shift; }
    return $self->{data}->[_id];
}

sub tRNAscan_id
{
    my $self = shift;
    if (@_) { $self->{data}->[_tRNAscan_id] = shift; }
    return $self->{data}->[_tRNAscan_id];
}

sub gtrnadb_id
{
    my $self = shift;
    if (@_) { $self->{data}->[_gtrnadb_id] = shift; }
    return $self->{data}->[_gtrnadb_id];
}

sub extdb_id
{
    my $self = shift;
    if (@_) { $self->{data}->[_extdb_id] = shift; }
    return $self->{data}->[_extdb_id];
}

sub clade
{
    my $self = shift;
    if (@_) { $self->{data}->[_clade] = shift; }
    return $self->{data}->[_clade];
}

sub isotype
{
    my $self = shift;
    if (@_) { $self->{data}->[_isotype] = shift; }
    return $self->{data}->[_isotype];
}

sub is_type2
{
    my $self = shift;
	my $type2 = 0;
	
	if (defined TYPE2_TRNA->{$self->{data}->[_clade]})
	{
		if (TYPE2_TRNA->{$self->{data}->[_clade]} =~ /$self->{data}->[_isotype]/)
		{
			$type2 = 1;
		}
	}
	
	return $type2;
}

sub anticodon
{
    my $self = shift;
    if (@_) { $self->{data}->[_anticodon] = shift; }
    return $self->{data}->[_anticodon];
}

sub ar_ac_pos
{
    my $self = shift;
    if (@_) { @{ $self->{data}->[_ar_ac_pos] } = @_; }
    return @{ $self->{data}->[_ar_ac_pos] };
}

sub add_ac_pos
{
    my $self = shift;
	my ($rel_start, $rel_end) = @_;
	my $ac = {};
	$ac->{rel_start} = $rel_start;
	$ac->{rel_end} = $rel_end;
	push(@{$self->{data}->[_ar_ac_pos]}, $ac);
}

sub get_ac_pos_count
{
	my $self = shift;
	return (scalar(@{$self->{data}->[_ar_ac_pos]}));
}

sub remove_ac_pos
{
	my $self = shift;
	my $index = shift;
	
	if ($index >= 0 and $index < $self->get_ac_pos_count())
	{
		splice(@{$self->{data}->[_ar_ac_pos]}, $index, 1);
	}	
}

sub ar_introns
{
    my $self = shift;
    if (@_) { @{ $self->{data}->[_ar_introns] } = @_; }
    return @{ $self->{data}->[_ar_introns] };
}

sub add_rel_intron
{
    my $self = shift;
	my ($rel_start, $rel_end, $type, $seq) = @_;
	my ($start, $end);
	if ($self->{data}->[_strand] eq "+")
	{
		$start = $self->{data}->[_start] + $rel_start - 1;
		$end = $self->{data}->[_start] + $rel_end - 1;
	}
	else
	{
		$end = $self->{data}->[_end] - $rel_start + 1;
		$start = $self->{data}->[_end] - $rel_end + 1;		
	}
	$self->add_intron($rel_start, $rel_end, $start, $end, $type, $seq);
}

sub add_intron
{
    my $self = shift;
	my ($rel_start, $rel_end, $abs_start, $abs_end, $type, $seq) = @_;
	my $intron = {};
	$intron->{rel_start} = $rel_start;
	$intron->{rel_end} = $rel_end;
	$intron->{start} = $abs_start;
	$intron->{end} = $abs_end;
	$intron->{type} = $type;
	$intron->{seq} = uc($seq);
	push(@{$self->{data}->[_ar_introns]}, $intron);
}

sub set_intron
{
    my $self = shift;
	my ($index, $rel_start, $rel_end, $type, $seq) = @_;

	$self->{data}->[_ar_introns]->[$index]->{rel_start} = $rel_start;
	$self->{data}->[_ar_introns]->[$index]->{rel_end} = $rel_end;
	$self->{data}->[_ar_introns]->[$index]->{type} = $type;
	$self->{data}->[_ar_introns]->[$index]->{seq} = $seq;
	if ($self->{data}->[_strand] eq "+")
	{
		$self->{data}->[_ar_introns]->[$index]->{start} = $self->{data}->[_start] + $rel_start - 1;
		$self->{data}->[_ar_introns]->[$index]->{end} = $self->{data}->[_start] + $rel_end - 1;
	}
	else
	{
		$self->{data}->[_ar_introns]->[$index]->{end} = $self->{data}->[_end] - $rel_start + 1;
		$self->{data}->[_ar_introns]->[$index]->{start} = $self->{data}->[_end] - $rel_end + 1;		
	}
}

sub get_intron
{
	my $self = shift;
	my $index = shift;
	
	if ($index >= 0 and $index < $self->get_intron_count())
	{
		return $self->{data}->[_ar_introns]->[$index];
	}
	else
	{
		return undef;		
	}
}

sub remove_intron
{
	my $self = shift;
	my $index = shift;
	
	if ($index >= 0 and $index < $self->get_intron_count())
	{
		splice(@{$self->{data}->[_ar_introns]}, $index, 1);
	}	
}

sub sort_introns
{
	my $self = shift;
	
	@{$self->{data}->[_ar_introns]} = sort {$a->{rel_start} <=> $b->{rel_start}} @{$self->{data}->[_ar_introns]};
}

sub merge_introns
{
	my $self = shift;
	my $last_end = -1;

	$self->sort_introns();
	for (my $ct = 0; $ct < $self->get_intron_count(); $ct++)
	{
		if ($last_end == $self->{data}->[_ar_introns]->[$ct]->{rel_start} - 1)
		{
			$last_end = $self->{data}->[_ar_introns]->[$ct]->{rel_end};
			$self->{data}->[_ar_introns]->[$ct]->{rel_start} = $self->{data}->[_ar_introns]->[$ct - 1]->{rel_start};
			if ($self->{data}->[_strand] eq "+")
			{
				$self->{data}->[_ar_introns]->[$ct]->{start} = $self->{data}->[_ar_introns]->[$ct - 1]->{start};
			}
			else
			{
				$self->{data}->[_ar_introns]->[$ct]->{end} = $self->{data}->[_ar_introns]->[$ct - 1]->{end};
			}
			if ($self->{data}->[_ar_introns]->[$ct - 1]->{type} eq "CI")
			{
				$self->{data}->[_ar_introns]->[$ct]->{type} = "CI";
			}
			
			$self->{data}->[_ar_introns]->[$ct - 1] = undef;
		}
		$last_end = $self->{data}->[_ar_introns]->[$ct]->{rel_end};
	}
	for (my $ct = ($self->get_intron_count() - 1); $ct >= 0; $ct--)
	{
		if (!defined $self->{data}->[_ar_introns]->[$ct])
		{
			splice(@{$self->{data}->[_ar_introns]}, $ct, 1);
		}		
	}	
}

sub get_intron_count
{
	my $self = shift;
	return (scalar(@{$self->{data}->[_ar_introns]}));
}

sub get_intron_seq
{
	my $self = shift;
	my $index = shift;
	my $seq = "";
	
	if ($index >= 0 and $index < $self->get_intron_count())
	{
		my $intron = $self->{data}->[_ar_introns]->[$index];
		$seq = uc(substr($self->{data}->[_seq], $intron->{rel_start} - 1, ($intron->{rel_end} - $intron->{rel_start} + 1)));
	}
	
	return $seq;
}

sub model
{
    my $self = shift;
    if (@_) { $self->{data}->[_model] = shift; }
    return $self->{data}->[_model];
}

sub score
{
    my $self = shift;
    if (@_) { $self->{data}->[_score] = shift; }
    return $self->{data}->[_score];
}

sub mat_score
{
    my $self = shift;
    if (@_) { $self->{data}->[_mat_score] = shift; }
    return $self->{data}->[_mat_score];
}

sub hmm_score
{
    my $self = shift;
    if (@_) { $self->{data}->[_hmm_score] = shift; }
    return $self->{data}->[_hmm_score];
}

sub ss_score
{
    my $self = shift;
    if (@_) { $self->{data}->[_ss_score] = shift; }
    return $self->{data}->[_ss_score];
}

sub h_domain_models
{
    my $self = shift;
    if (@_) { %{ $self->{data}->[_h_domain_models] } = @_; }
    return %{ $self->{data}->[_h_domain_models] };
}

sub set_domain_model
{
    my $self = shift;
	my ($alg, $score) = @_;
	my $model_hit = {};
	$model_hit->{score} = $score;
	$model_hit->{mat_score} = 0;
	$model_hit->{hmm_score} = 0;
	$model_hit->{ss_score} = 0;	
	$self->{data}->[_h_domain_models]->{$alg} = $model_hit;
}

sub get_domain_model
{
	my $self = shift;
	my ($alg) = @_;
	if (defined $self->{data}->[_h_domain_models]->{$alg})
	{
		return $self->{data}->[_h_domain_models]->{$alg};
	}
	else
	{
		return undef;
	}
}

sub update_domain_model
{
	my $self = shift;
	my ($alg, $score, $mat_score, $hmm_score, $ss_score) = @_;
	if (defined $self->{data}->[_h_domain_models]->{$alg})
	{
		$self->{data}->[_h_domain_models]->{$alg}->{score} = $score;
		$self->{data}->[_h_domain_models]->{$alg}->{mat_score} = $mat_score;
		$self->{data}->[_h_domain_models]->{$alg}->{hmm_score} = $hmm_score;
		$self->{data}->[_h_domain_models]->{$alg}->{ss_score} = $ss_score;
	}
}

sub set_default_scores
{
	my $self = shift;
	
	if (defined $self->{data}->[_h_domain_models]->{cove})
	{
		$self->{data}->[_score] = $self->{data}->[_h_domain_models]->{cove}->{score};
		$self->{data}->[_mat_score] = $self->{data}->[_h_domain_models]->{cove}->{mat_score};
	}
	else
	{
		$self->{data}->[_score] = $self->{data}->[_h_domain_models]->{infernal}->{score};
		$self->{data}->[_mat_score] = $self->{data}->[_h_domain_models]->{infernal}->{mat_score};		
	}
	if (defined $self->{data}->[_h_domain_models]->{infernal})
	{
		$self->{data}->[_hmm_score] = $self->{data}->[_h_domain_models]->{infernal}->{hmm_score};
		$self->{data}->[_ss_score] = $self->{data}->[_h_domain_models]->{infernal}->{ss_score};		
	}
	else
	{
		$self->{data}->[_hmm_score] = $self->{data}->[_h_domain_models]->{cove}->{hmm_score};
		$self->{data}->[_ss_score] = $self->{data}->[_h_domain_models]->{cove}->{ss_score};		
	}
}

sub get_default_score_type
{
	my $self = shift;
	my $type = "";
	if (defined $self->{data}->[_h_domain_models]->{cove})
	{
		$type = "cove";
	}
	elsif (defined $self->{data}->[_h_domain_models]->{infernal})
	{
		$type = "infernal";
	}
	return $type;
}

sub seq
{
    my $self = shift;
    if (@_) { $self->{data}->[_seq] = shift; }
    return $self->{data}->[_seq];
}

sub mat_seq
{
    my $self = shift;
    if (@_) { $self->{data}->[_mat_seq] = shift; }
    return $self->{data}->[_mat_seq];
}

sub set_mature_tRNA
{
    my $self = shift;

	my $set_ss = 0;
	if ($self->{data}->[_mat_ss] eq "")
	{
		$set_ss = 1;
	}
	if ($self->{data}->[_mat_seq] eq "")
	{
		if ($self->get_intron_count() == 0)
		{
			$self->{data}->[_mat_seq] = $self->{data}->[_seq];
			$self->{data}->[_mat_ss] = $self->{data}->[_ss];
		}
		else
		{
			my @introns = @{$self->{data}->[_ar_introns]};
			for (my $i = 0; $i < $self->get_intron_count(); $i++)
			{
				if ($i == 0)
				{
					$self->{data}->[_mat_seq] = substr($self->{data}->[_seq], 0, $introns[$i]->{rel_start} - 1);
					if ($set_ss)
					{
						$self->{data}->[_mat_ss] = substr($self->{data}->[_ss], 0, $introns[$i]->{rel_start} - 1);
					}
				}
				else
				{
					$self->{data}->[_mat_seq] .= substr($self->{data}->[_seq], $introns[$i-1]->{rel_end}, $introns[$i]->{rel_start} - $introns[$i-1]->{rel_end} - 1);
					if ($set_ss)
					{
						$self->{data}->[_mat_ss] .= substr($self->{data}->[_ss], $introns[$i-1]->{rel_end}, $introns[$i]->{rel_start} - $introns[$i-1]->{rel_end} - 1);
					}
				}
			}
			$self->{data}->[_mat_seq] .= substr($self->{data}->[_seq], $introns[$self->get_intron_count()-1]->{rel_end});
			if ($set_ss)
			{
				$self->{data}->[_mat_ss] .= substr($self->{data}->[_ss], $introns[$self->get_intron_count()-1]->{rel_end});
			}
		}
	}
}

sub mat_ss
{
    my $self = shift;
    if (@_) { $self->{data}->[_mat_ss] = shift; }
    return $self->{data}->[_mat_ss];
}

sub ss
{
    my $self = shift;
    if (@_) { $self->{data}->[_ss] = shift; }
    return $self->{data}->[_ss];
}

sub sprinzl_align
{
    my $self = shift;
    if (@_) { $self->{data}->[_sprinzl_align] = shift; }
    return $self->{data}->[_sprinzl_align];
}

sub sprinzl_ss
{
    my $self = shift;
    if (@_) { $self->{data}->[_sprinzl_ss] = shift; }
    return $self->{data}->[_sprinzl_ss];
}

sub pseudo
{
    my $self = shift;
    if (@_) { $self->{data}->[_pseudo] = shift; }
    return $self->{data}->[_pseudo];
}

sub is_pseudo
{
    my $self = shift;
	if ($self->{data}->[_pseudo])
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub trunc
{
    my $self = shift;
    if (@_) { $self->{data}->[_trunc] = shift; }
    return $self->{data}->[_trunc];
}

sub is_trunc
{
    my $self = shift;
	if ($self->{data}->[_trunc] ne "")
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub category
{
    my $self = shift;
    if (@_) { $self->{data}->[_category] = shift; }
    return $self->{data}->[_category];
}

sub is_numt
{
    my $self = shift;
	if ($self->{data}->[_category] eq "numt")
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub is_mito
{
    my $self = shift;
	if ($self->{data}->[_category] eq "mito")
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub is_undetermined
{
    my $self = shift;
	if ($self->{data}->[_category] eq "undetermined_ac")
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub is_cytosolic
{
    my $self = shift;
	if ($self->{data}->[_category] eq "cyto")
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub is_organelle
{
    my $self = shift;
	if ($self->{data}->[_category] eq "org")
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub is_mito_iso_conflict
{
    my $self = shift;
	if ($self->{data}->[_category] eq "mito_iso_conflict")
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub is_mito_noncanonical_ac
{
    my $self = shift;
	if ($self->{data}->[_category] eq "mito_noncanonical_ac")
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub is_mito_ac_mislocation
{
    my $self = shift;
	if ($self->{data}->[_category] eq "mito_ac_mislocation")
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub is_mito_mismatch_ac
{
    my $self = shift;
	if ($self->{data}->[_category] eq "mito_mismatch_ac")
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub note
{
    my $self = shift;
    if (@_) { $self->{data}->[_note] = shift; }
    return $self->{data}->[_note];
}

sub upstream
{
    my $self = shift;
    if (@_) { $self->{data}->[_upstream] = shift; }
    return $self->{data}->[_upstream];
}

sub downstream
{
    my $self = shift;
    if (@_) { $self->{data}->[_downstream] = shift; }
    return $self->{data}->[_downstream];
}

sub hit_source
{
    my $self = shift;
    if (@_) { $self->{data}->[_hit_source] = shift; }
    return $self->{data}->[_hit_source];
}

sub src_seqid
{
    my $self = shift;
    if (@_) { $self->{data}->[_src_seqid] = shift; }
    return $self->{data}->[_src_seqid];
}

sub src_seqlen
{
    my $self = shift;
    if (@_) { $self->{data}->[_src_seqlen] = shift; }
    return $self->{data}->[_src_seqlen];
}

sub ar_multi_models
{
    my $self = shift;
    if (@_) { @{ $self->{data}->[_ar_multi_models] } = @_; }
    return @{ $self->{data}->[_ar_multi_models] };
}

sub add_model_hit
{
    my $self = shift;
	my ($type, $model, $score, $ss) = @_;
	my $model_hit = {};
	$model_hit->{type} = $type;	
	$model_hit->{model} = $model;
	$model_hit->{score} = $score;
	$model_hit->{ss} = $ss;
	push(@{$self->{data}->[_ar_multi_models]}, $model_hit);
}

sub sort_multi_models
{
	my $self = shift;
	my $key = shift;
	if ($key eq "model")
	{
		@{$self->{data}->[_ar_multi_models]} = sort {$a->{type} cmp $b->{type} || $a->{model} cmp $b->{model}} @{$self->{data}->[_ar_multi_models]};
	}
	elsif ($key eq "score")
	{
		@{$self->{data}->[_ar_multi_models]} = sort {$b->{score} <=> $a->{score}} @{$self->{data}->[_ar_multi_models]};		
	}
}

sub get_model_hit
{
	# assume pre-sort
	my $self = shift;
	my ($type, $search_name) = @_;
	my $model = "";
	my $score = 0;
	my $ss = "";
	
	my $index = $self->bsearch_model_hit($type, $search_name);
	if ($index > -1)
	{
		$model = $self->{data}->[_ar_multi_models]->[$index]->{model};
		$score = $self->{data}->[_ar_multi_models]->[$index]->{score};
		$ss = $self->{data}->[_ar_multi_models]->[$index]->{ss};		
	}
	return ($model, $score, $ss);
}

sub get_highest_score_model
{
	my $self = shift;
	my $type = "None";
	my $model = "None";
	my $score = 0;
	my $ss = "";

	if (scalar(@{$self->{data}->[_ar_multi_models]}) > 0)
	{
		$self->sort_multi_models("score");
		$type = $self->{data}->[_ar_multi_models]->[0]->{type};
		$model = $self->{data}->[_ar_multi_models]->[0]->{model};
		$score = $self->{data}->[_ar_multi_models]->[0]->{score};
		$ss = $self->{data}->[_ar_multi_models]->[0]->{ss};		
	}
	return ($type, $model, $score, $ss);	
}

sub bsearch_model_hit
{
	my $self = shift;
    my ($x, $y) = @_;            
    my ($l, $u) = (0, @{$self->{data}->[_ar_multi_models]} - 1);  
    my $i;                       
    while ($l <= $u)
	{
		$i = int(($l + $u)/2);
		if ($self->{data}->[_ar_multi_models]->[$i]->{type} lt $x)
		{
		    $l = $i+1;
		}
		elsif ($self->{data}->[_ar_multi_models]->[$i]->{type} gt $x)
		{
		    $u = $i-1;
		} 
		else
		{
			if ($self->{data}->[_ar_multi_models]->[$i]->{model} lt $y)
			{
				$l = $i+1;
			}
			elsif ($self->{data}->[_ar_multi_models]->[$i]->{model} gt $y)
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

sub ar_pos_sprinzl_map
{
    my $self = shift;
    return @{ $self->{data}->[_ar_pos_sprinzl_map] };
}

sub h_sprinzl_pos_map
{
    my $self = shift;
    return %{ $self->{data}->[_h_sprinzl_pos_map] };
}

sub map_sprinzl_pos
{
	my $self = shift;
	my $ar_sprinzl_pos = shift;
	
	$self->{data}->[_ar_pos_sprinzl_map] = [];
	$self->{data}->[_h_sprinzl_pos_map] = {};
	
	my @bases = split(//, $self->{data}->[_sprinzl_align]);
	
	my $pos_count = 0;
	my $ins_count = 0;
	for (my $i = 0; $i < scalar(@bases); $i++)
	{
		if ($bases[$i] =~ /^[ACGT-]$/)
		{
			push(@{$self->{data}->[_ar_pos_sprinzl_map]}, $ar_sprinzl_pos->[$pos_count]);
			$self->{data}->[_h_sprinzl_pos_map]->{$ar_sprinzl_pos->[$pos_count]} = $i;
			$pos_count++;
			$ins_count = 0;
		}
		else
		{
			$ins_count++;
			push(@{$self->{data}->[_ar_pos_sprinzl_map]}, $ar_sprinzl_pos->[$pos_count-1].":i".$ins_count);
			$self->{data}->[_h_sprinzl_pos_map]->{$ar_sprinzl_pos->[$pos_count-1].":i".$ins_count} = $i;
		}
	}
}

sub get_sprinzl_pos
{
	my $self = shift;
	my $rel_pos = shift;
	
	my $pos = "";
	
	if ($rel_pos >= 0 and $rel_pos < scalar(@{$self->{data}->[_ar_pos_sprinzl_map]}))
	{
		$pos = $self->{data}->[_ar_pos_sprinzl_map]->[$rel_pos];
	}
	return $pos;
}

sub get_rel_align_pos
{
	my $self = shift;
	my $sprinzl_pos = shift;
	
	my $rel_pos = -1;
	
	if (defined $self->{data}->[_h_sprinzl_pos_map]->{$sprinzl_pos})
	{
		$rel_pos = $self->{data}->[_h_sprinzl_pos_map]->{$sprinzl_pos};
	}
	return $rel_pos;
}

sub get_base_at_sprinzl
{
	my $self = shift;
	my $sprinzl_pos = shift;

	my $rel_pos = $self->get_rel_align_pos($sprinzl_pos);
	return ($self->get_base_at_rel_align_pos($rel_pos));
}

sub get_base_at_rel_align_pos
{
	my $self = shift;
	my $rel_pos = shift;

	my $base = "";
	if ($rel_pos != -1)
	{
		$base = substr($self->{data}->[_sprinzl_align], $rel_pos, 1);
	}
	return $base;
}

sub get_seq_at_rel_align_pos
{
	my $self = shift;
	my $rel_pos_start = shift;
	my $rel_pos_end = shift;

	my $len = $rel_pos_end - $rel_pos_start + 1;

	my $seq = "";
	if ($rel_pos_start != -1 and $len > 0)
	{
		$seq = substr($self->{data}->[_sprinzl_align], $rel_pos_start, $len);
	}
	return $seq;
}

sub get_rel_pos
{
	my $self = shift;
	my $coord = shift;
	
	my $rel_pos = -1;
	
	if ($coord >= $self->{data}->[_start] && $coord <= $self->{data}->[_end])
	{
		if ($self->{data}->[_strand] eq "+")
		{
			$rel_pos = $coord - $self->{data}->[_start] + 1;
		}
		else
		{
			$rel_pos = $self->{data}->[_end] - $coord + 1;
		}
	}
	
	return $rel_pos;
}

sub get_rel_mature_pos
{
	my $self = shift;
	my $coord = shift;

	my $rel_pos = $self->get_rel_pos($coord);
	my $rel_mature_pos = $rel_pos;
	
	for (my $i = ($self->get_intron_count() - 1); $i >= 0; $i--)
	{
		my $intron = $self->{data}->[_ar_introns]->[$i];
		if ($rel_pos > $intron->{rel_end})
		{
			$rel_mature_pos = $rel_mature_pos - ($intron->{rel_end} - $intron->{rel_start} + 1);
		}
	}
	
	return $rel_mature_pos;
}

sub is_intronic
{
	my $self = shift;
	my $coord = shift;
	
	my $intronic = 0;
	my $rel_pos = $self->get_rel_pos($coord);
	if ($rel_pos != -1)
	{
		for (my $i = 0; $i < $self->get_intron_count(); $i++)
		{
			my $intron = $self->{data}->[_ar_introns]->[$i];
			if ($rel_pos >= $intron->{rel_start} && $rel_pos <= $intron->{rel_end})
			{
				$intronic = 1;
			}
		}
	}
	return $intronic;
}

sub convert_sprinzl_pos
{
	my $self = shift;
	my $coord = shift;
	
	my $sprinzl_pos = "";
	
	if (!$self->is_intronic($coord))
	{
		my $rel_mature_pos = $self->get_rel_mature_pos($coord);
		my @bases = split(//, $self->{data}->[_sprinzl_align]);
		my $ct = 0;
		for (my $i = 0; $i < scalar(@bases); $i++)
		{
			if ($bases[$i] ne "-")
			{
				$ct++;
				if ($ct == $rel_mature_pos)
				{
					$sprinzl_pos = $self->get_sprinzl_pos($i);
					$i = scalar(@bases);
				}
			}
		}
	}
	return $sprinzl_pos;
}

sub get_seq_fragment
{
	my $self = shift;
	my ($rel_start, $rel_end) = @_;
	
	my $seq = "";
	if ($rel_start > 0 and $rel_end <= &len($self))
	{
		$seq = substr($self->{data}->[_seq], $rel_start - 1, $rel_end - $rel_start + 1);
	}
	return $seq;
}

sub ar_non_canonical
{
    my $self = shift;
    return @{ $self->{data}->[_ar_non_canonical] };
}

sub add_non_canonical
{
	my $self = shift;
	my $rel_pos = shift;
	my $nc = shift;
	
	$self->{data}->[_ar_non_canonical]->[$rel_pos] = $nc;
}

sub get_non_canonical
{
	my $self = shift;
	my $sprinzl_pos = shift;
	
	my $nc = "";
	my $rel_pos = $self->get_rel_align_pos($sprinzl_pos);
	if ($rel_pos != -1)
	{
		$nc = $self->{data}->[_ar_non_canonical]->[$rel_pos];
	}
	return $nc;
}

sub get_mismatch_count
{
	my $self = shift;
	
	my $count = 0;
	for (my $i = 0; $i < scalar(@{$self->{data}->[_ar_non_canonical]}); $i++)
	{
		if ($self->{data}->[_ar_non_canonical]->[$i] eq "M")
		{
			$count++;
		}		
	}
	return (int($count / 2));
}

sub get_noncanonical_count
{
	my $self = shift;
	
	my $count = 0;
	for (my $i = 0; $i < scalar(@{$self->{data}->[_ar_non_canonical]}); $i++)
	{
		if ($self->{data}->[_ar_non_canonical]->[$i] eq "N")
		{
			$count++;
		}		
	}
	return (int($count / 2));
}

sub is_mismatch
{
	my $self = shift;
	my $rel_pos = shift;

	my $mismatch = 0;
	if ($self->{data}->[_ar_non_canonical]->[$rel_pos] eq "M")
	{
		$mismatch = 1;
	}
	return $mismatch;
}

sub is_insertion
{
	my $self = shift;
	my $rel_pos = shift;

	my $insert = 0;
	if ($self->{data}->[_ar_non_canonical]->[$rel_pos] eq "I")
	{
		$insert = 1;
	}
	return $insert;
}

sub is_deletion
{
	my $self = shift;
	my $rel_pos = shift;

	my $delete = 0;
	if ($self->{data}->[_ar_non_canonical]->[$rel_pos] eq "D")
	{
		$delete = 1;
	}
	return $delete;
}

sub is_non_canonical
{
	my $self = shift;
	my $rel_pos = shift;

	my $nc = 0;
	if ($self->{data}->[_ar_non_canonical]->[$rel_pos] eq "N")
	{
		$nc = 1;
	}
	return $nc;
}

sub get_annotated_mismatches
{
	my $self = shift;
	my $h_sprinzl_pairs = shift;
	my $value = "";
	
	for (my $i = 0; $i < scalar(@{$self->{data}->[_ar_non_canonical]}); $i++)
	{
		if ($self->{data}->[_ar_non_canonical]->[$i] eq "M")
		{
			if (defined $h_sprinzl_pairs->{$self->get_sprinzl_pos($i)})
			{
				my $pos2 = $h_sprinzl_pairs->{$self->get_sprinzl_pos($i)};
				my $base2 = $self->get_base_at_sprinzl($pos2);
				if ($value ne "")
				{
					$value .= " ";
				}				
				$value .= $self->get_base_at_rel_align_pos($i).$self->get_sprinzl_pos($i).":".$base2.$pos2;
			}
		}		
	}
	return $value;
}

sub get_annotated_noncanonical
{
	my $self = shift;
	my $value = "";
	
	for (my $i = 0; $i < scalar(@{$self->{data}->[_ar_non_canonical]}); $i++)
	{
		if ($self->{data}->[_ar_non_canonical]->[$i] eq "N")
		{
			if ($value ne "")
			{
				$value .= " ";
			}				
			$value .= $self->get_base_at_rel_align_pos($i).$self->get_sprinzl_pos($i);
		}		
	}
	return $value;
}

sub check_stem_mismatches
{
	my $self = shift;
	my $h_sprinzl_pairs = shift;
	
	my @bases = split(//, $self->{data}->[_sprinzl_align]);
	
	for (my $i = 0; $i < scalar(@bases); $i++)
	{
		if (index($self->{data}->[_ar_pos_sprinzl_map]->[$i], ":i") == -1)
		{
			if (defined $h_sprinzl_pairs->{$self->{data}->[_ar_pos_sprinzl_map]->[$i]} && $self->{data}->[_ar_pos_sprinzl_map]->[$i] ne '13' && index($self->{data}->[_ar_pos_sprinzl_map]->[$i], "e") == -1)
			{
				my $b1 = $bases[$i];
				my $pos2 = $self->get_rel_align_pos($h_sprinzl_pairs->{$self->{data}->[_ar_pos_sprinzl_map]->[$i]});
				my $b2 = $bases[$pos2];
				if ($b1 eq "A" && $b2 eq "T") {}
				elsif ($b1 eq "T" && $b2 eq "A") {}
				elsif ($b1 eq "T" && $b2 eq "G") {}
				elsif ($b1 eq "G" && $b2 eq "C") {}
				elsif ($b1 eq "G" && $b2 eq "T") {}
				elsif ($b1 eq "C" && $b2 eq "G") {}
				elsif ($b1 eq "-" || $b2 eq "-") {}
				else
				{
					$self->{data}->[_ar_non_canonical]->[$i] = "M";
					$self->{data}->[_ar_non_canonical]->[$pos2] = "M";
				}
			}
		}
	}
}

sub is_WC_base_pair
{
	my $self = shift;
	my $b1 = shift;
	my $b2 = shift;
	
	my $wc = 1;	
	if ($b1 eq "A" && $b2 eq "T") {}
	elsif ($b1 eq "A" && $b2 eq "U") {}
	elsif ($b1 eq "T" && $b2 eq "A") {}
	elsif ($b1 eq "U" && $b2 eq "A") {}
	elsif ($b1 eq "G" && $b2 eq "C") {}
	elsif ($b1 eq "C" && $b2 eq "G") {}
	elsif ($b1 eq "-" || $b2 eq "-") {}
	else
	{
		$wc = 0;
	}
	return $wc;
}

sub is_GU_base_pair
{
	my $self = shift;
	my $b1 = shift;
	my $b2 = shift;
	
	my $gu = 1;	
	if ($b1 eq "G" && $b2 eq "T") {}
	elsif ($b1 eq "G" && $b2 eq "U") {}
	elsif ($b1 eq "T" && $b2 eq "G") {}
	elsif ($b1 eq "U" && $b2 eq "G") {}
	else
	{
		$gu = 0;
	}
	return $gu;
}

1;
