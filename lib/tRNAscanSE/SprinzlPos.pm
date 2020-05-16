# tRNAscanSE/SprinzlPos.pm
# This class contains parameters and functions describing Sprinzl positions in tRNAs.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2017 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::SprinzlPos;

use strict;

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
	
	$self->{sprinzl_pos} = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16',
					'17','17a','18','19','20','20a','20b','21','22','23','24','25',
					'26','27','28','29','30','31','32','33','34','35','36','37','38',
					'39','40','41','42','43','44','45',
					'e11','e12','e13','e14','e15','e16','e17',
					'e1','e2','e3','e4','e5',
					'e27','e26','e25','e24','e23','e22','e21',
					'46','47','48','49','50','51','52','53','54','55','56','57','58',
					'59','60','61','62','63','64','65','66','67','68','69','70','71',
					'72','73','74','75','76'];
	
	$self->{sprinzl_pairs} = {'1'=>'72', '2'=>'71', '3'=>'70', '4'=>'69', '5'=>'68', '6'=>'67', '7'=>'66',
					  '10'=>'25', '11'=>'24', '12'=>'23', '13'=>'22',
					  '27'=>'43', '28'=>'42', '29'=>'41', '30'=>'40', '31'=>'39',
					  'e11'=>'e21', 'e12'=>'e22', 'e13'=>'e23', 'e14'=>'e24', 'e15'=>'e25', 'e16'=>'e26', 'e17'=>'e27',
					  '49'=>'65', '50'=>'64', '51'=>'63', '52'=>'62', '53'=>'61'};

	$self->{rev_sprinzl_pairs} = {'72'=>'1', '71'=>'2', '70'=>'3', '69'=>'4', '68'=>'5', '67'=>'6', '66'=>'7',
					  '25'=>'10', '24'=>'11', '23'=>'12', '22'=>'13',
					  '43'=>'27', '42'=>'28', '41'=>'29', '40'=>'30', '39'=>'31',
					  'e21'=>'e11', 'e22'=>'e12', 'e23'=>'e13', 'e24'=>'e14', 'e25'=>'e15', 'e26'=>'e16', 'e27'=>'e17',
					  '65'=>'49', '64'=>'50', '63'=>'51', '62'=>'52', '61'=>'53'};
	
	$self->{stem_pos} = {'1'=>'A1', '2'=>'A2', '3'=>'A3', '4'=>'A4', '5'=>'A5', '6'=>'A6', '7'=>'A7',
					  '10'=>'D1', '11'=>'D2', '12'=>'D3', '13'=>'D4',
					  '27'=>'C1', '28'=>'C2', '29'=>'C3', '30'=>'C4', '31'=>'C5',
					  'e11'=>'V1', 'e12'=>'V2', 'e13'=>'V3', 'e14'=>'V4', 'e15'=>'V5', 'e16'=>'V6', 'e17'=>'V7',
					  '49'=>'T1', '50'=>'T2', '51'=>'T3', '52'=>'T4', '53'=>'T5'};
	
	$self->{regions} = {"5P1"=>"5p Acceptor Stem", "3P1"=>"3p Acceptor Stem", "L1"=>"Acceptor-D-arm-linker", "5P2"=>"5p D-arm", "3P2"=>"3p D-arm",
						"L2"=>"D-loop", "L3"=>"D-arm-Anticodon-linker", "5P3"=>"5p Anticodon Stem", "3P3"=>"3p Anticodon Stem", "L4"=>"Anticodon Loop",
						"L5"=>"Variable Loop", "P4"=>"Variable Stem", "5P5"=>"5p T-arm", "3P5"=>"3p T-arm", "L6"=>"T-Psi-C Loop", "L7"=>"3p end"};
	
	$self->{ss_pos} = {'1'=>'5P1', '2'=>'5P1', '3'=>'5P1', '4'=>'5P1', '5'=>'5P1', '6'=>'5P1', '7'=>'5P1',
					   '72'=>'3P1', '71'=>'3P1', '70'=>'3P1', '69'=>'3P1', '68'=>'3P1', '67'=>'3P1', '66'=>'3P1',
					   '8'=>'L1', '9'=>'L1', 
					  '10'=>'5P2', '11'=>'5P2', '12'=>'5P2', '13'=>'5P2',
					  '25'=>'3P2', '24'=>'3P2', '23'=>'3P2', '22'=>'3P2',
					  '14'=>'L2', '15'=>'L2', '16'=>'L2', '17'=>'L2', '17a'=>'L2', '18'=>'L2', '19'=>'L2', '20'=>'L2', '20a'=>'L2', '20b'=>'L2', '21'=>'L2',
					  '26'=>'L3',
					  '27'=>'5P3', '28'=>'5P3', '29'=>'5P3', '30'=>'5P3', '31'=>'5P3',
					  '43'=>'3P3', '42'=>'3P3', '41'=>'3P3', '40'=>'3P3', '39'=>'3P3',
					  '32'=>'L4', '33'=>'L4', '34'=>'L4', '35'=>'L4', '36'=>'L4', '37'=>'L4', '38'=>'L4',
					  '44'=>'L5', '45'=>'L5', '46'=>'L5', '47'=>'L5', '48'=>'L5', 
					  'e11'=>'P4', 'e12'=>'P4', 'e13'=>'P4', 'e14'=>'P4', 'e15'=>'P4', 'e16'=>'P4', 'e17'=>'P4',
					  'e21'=>'P4', 'e22'=>'P4', 'e23'=>'P4', 'e24'=>'P4', 'e25'=>'P4', 'e26'=>'P4', 'e27'=>'P4',
					  'e1'=>'P4', 'e2'=>'P4', 'e3'=>'P4', 'e4'=>'P4', 'e5'=>'P4',
					  '49'=>'5P5', '50'=>'5P5', '51'=>'5P5', '52'=>'5P5', '53'=>'5P5',
					  '65'=>'3P5', '64'=>'3P5', '63'=>'3P5', '62'=>'3P5', '61'=>'3P5',
					  '54'=>'L6', '55'=>'L6', '56'=>'L6', '57'=>'L6', '58'=>'L6', '59'=>'L6', '60'=>'L6',
					  '73'=>'L7', '74'=>'L7', '75'=>'L7', '76'=>'L7'};

	$self->{universal} = {'8'=>'T', '14'=>'A', '18'=>'G', '19'=>'G', '21'=>'A', '33'=>'T', '53'=>'G',
						  '54'=>'T', '55'=>'T', '56'=>'C', '58'=>'A', '74'=>'C', '75'=>'C', '76'=>'A'};
}

sub sprinzl_pos
{
    my $self = shift;
    return @{ $self->{sprinzl_pos} };
}

sub sprinzl_pairs
{
    my $self = shift;
    return %{ $self->{sprinzl_pairs} };
}

sub sprinzl_pair
{
    my $self = shift;
	my $pos = shift;
	
	my $value = "";
	if (defined $self->{sprinzl_pairs}->{$pos})
	{
		$value = $self->{sprinzl_pairs}->{$pos};
	}	
    return $value;
}

sub rev_sprinzl_pairs
{
    my $self = shift;
    return %{ $self->{rev_sprinzl_pairs} };
}

sub rev_sprinzl_pair
{
    my $self = shift;
	my $pos = shift;
	
	my $value = "";
	if (defined $self->{rev_sprinzl_pairs}->{$pos})
	{
		$value = $self->{rev_sprinzl_pairs}->{$pos};
	}	
    return $value;
}

sub stem_pos
{
    my $self = shift;
    return %{ $self->{stem_pos} };
}

sub stem_pos_code
{
    my $self = shift;
	my $pos = shift;
	
	my $value = "";
	if (defined $self->{stem_pos}->{$pos})
	{
		$value = $self->{stem_pos}->{$pos};
	}	
    return $value;
}

sub regions
{
    my $self = shift;
    return %{ $self->{regions} };
}

sub region_desc
{
    my $self = shift;
	my $region = shift;
	
	my $value = "";
    if (defined $self->{regions}->{$region})
	{
		$value = $self->{regions}->{$region};
	}
	return $value;
}

sub ss_pos
{
    my $self = shift;
    return %{ $self->{ss_pos} };
}

sub pos_region
{
    my $self = shift;
	my $pos = shift;
	
	my $value = "";
	if (defined $self->{ss_pos}->{$pos})
	{
		$value = $self->{ss_pos}->{$pos};
	}
    return $value;
}

sub pos_region_desc
{
    my $self = shift;
	my $pos = shift;
	
	my $value = "";
	my $search = $pos;
	if (index($pos, ":i") > -1)
	{
		$search = substr($pos, 0, index($pos, ":i"));
	}	
	if (defined $self->{ss_pos}->{$search})
	{
		$value = $self->{regions}->{$self->{ss_pos}->{$search}};
	}
    return $value;
}

sub universal_base
{
	my $self = shift;
	my $pos = shift;
	
	my $value = "";
	if (defined $self->{universal}->{$pos})
	{
		$value = $self->{universal}->{$pos};
	}
	return $value;
}

1;

