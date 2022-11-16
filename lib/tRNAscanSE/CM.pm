# tRNAscanSE/CM.pm
# This class contains parameters and functions for running CM tRNA search used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2020 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::CM;

use strict;
use File::Copy;
use tRNAscanSE::Configuration;
use tRNAscanSE::Options;
use tRNAscanSE::Utils;
use tRNAscanSE::LogFile;
use tRNAscanSE::ScanResult;
use tRNAscanSE::SS;
use tRNAscanSE::Sequence;
use tRNAscanSE::tRNA;
use tRNAscanSE::ArraytRNA;
use tRNAscanSE::IntResultFile;
use tRNAscanSE::MultiResultFile;
use tRNAscanSE::CMscanResultFile;
use tRNAscanSE::ArrayCMscanResults;
use tRNAscanSE::FpScanResultFile;

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
    
    $self->{CM_check_for_introns} = 0;          # check for non-canonical introns
    $self->{CM_check_for_split_halves} = 0;     # check for split tRNA fragments                                                
    $self->{skip_pseudo_filter} = 0;            # enable filter for psuedogenes (Cove score <40,
                                                # primary struct score <10 bits, secondary 
                                                # structure score < 5 bits)
    
    $self->{get_hmm_score} = 0;                 # also score tRNA with covariance model
                                                #  without sec structure info, similar
                                                #  to getting hmm score for match of 
                                                #  seq to tRNA hmm  (-H option)
                                                
    
    # Convariance model file path
    $self->{main_cm_file_path} = {};            
    $self->{mainNS_cm_file_path} = {};
    $self->{cove_cm_file_path} = '';            
    $self->{intron_cm_file_path} = {};    
    $self->{arch_five_half_cm_file_path} = '';
    $self->{arch_three_half_cm_file_path} = '';
    $self->{Pselc_cm_file_path} = '';
    $self->{Eselc_cm_file_path} = '';
    $self->{isotype_cm_db_file_path} = '';
    $self->{mito_isotype_cm_db_file_path} = '';    
    $self->{isotype_cm_file_paths} = {};
    $self->{mito_isotype_cm_file_paths} = {};    
    
    $self->{isotype_cm_cutoff} = {};

    $self->{covels_bin} = "covels-SE";          # Application executable name
    $self->{coves_bin} = "coves-SE";
    $self->{cmsearch_bin} = "cmsearch";
    $self->{cmscan_bin} = "cmscan";

    $self->{infernal_thread} = -1;
    
    $self->{tab_results} = +[];
}

sub CM_mode
{
    my $self = shift;
    if (@_) { $self->{CM_mode} = shift; }
    return $self->{CM_mode};
}

sub cove_mode
{
    my $self = shift;
    return ($self->{CM_mode} eq 'cove');
}

sub infernal_mode
{
    my $self = shift;
    return ($self->{CM_mode} eq 'infernal');
}

sub cm_cutoff
{
    my $self = shift;
    if (@_) { $self->{cm_cutoff} = shift; }
    return $self->{cm_cutoff};
}

sub organelle_cm_cutoff
{
    my $self = shift;
    if (@_) { $self->{organelle_cm_cutoff} = shift; }
    return $self->{organelle_cm_cutoff};
}

sub infernal_fp_cutoff
{
    my $self = shift;
    if (@_) { $self->{infernal_fp_cutoff} = shift; }
    return $self->{infernal_fp_cutoff};
}

sub BHB_cm_cutoff
{
    my $self = shift;
    if (@_) { $self->{BHB_cm_cutoff} = shift; }
    return $self->{BHB_cm_cutoff};
}

sub max_tRNA_length
{
    my $self = shift;
    if (@_) { $self->{max_tRNA_length} = shift; }
    return $self->{max_tRNA_length};
}

sub max_cove_tRNA_length
{
    my $self = shift;
    if (@_) { $self->{max_cove_tRNA_length} = shift; }
    return $self->{max_cove_tRNA_length};
}

sub max_cmsearch_tRNA_length
{
    my $self = shift;
    if (@_) { $self->{max_cmsearch_tRNA_length} = shift; }
    return $self->{max_cmsearch_tRNA_length};
}

sub CM_check_for_introns
{
    my $self = shift;
    if (@_) { $self->{CM_check_for_introns} = shift; }
    return $self->{CM_check_for_introns};
}

sub CM_check_for_split_halves
{
    my $self = shift;
    if (@_) { $self->{CM_check_for_split_halves} = shift; }
    return $self->{CM_check_for_split_halves};
}

sub min_tRNA_no_intron
{
    my $self = shift;
    if (@_) { $self->{min_tRNA_no_intron} = shift; }
    return $self->{min_tRNA_no_intron};
}

sub min_intron_length
{
    my $self = shift;
    if (@_) { $self->{min_intron_length} = shift; }
    return $self->{min_intron_length};
}

sub skip_pseudo_filter
{
    my $self = shift;
    if (@_) { $self->{skip_pseudo_filter} = shift; }
    return $self->{skip_pseudo_filter};
}

sub min_pseudo_filter_score
{
    my $self = shift;
    if (@_) { $self->{min_pseudo_filter_score} = shift; }
    return $self->{min_pseudo_filter_score};
}

sub min_ss_score
{
    my $self = shift;
    if (@_) { $self->{min_ss_score} = shift; }
    return $self->{min_ss_score};
}

sub min_hmm_score
{
    my $self = shift;
    if (@_) { $self->{min_hmm_score} = shift; }
    return $self->{min_hmm_score};
}

sub get_hmm_score
{
    my $self = shift;
    if (@_) { $self->{get_hmm_score} = shift; }
    return $self->{get_hmm_score};
}

sub main_cm_file_path
{
    my $self = shift;
    if (@_) { %{$self->{main_cm_file_path}} = shift; }
    return %{$self->{main_cm_file_path}};
}

sub add_main_cm_file_path
{
    my $self = shift;
    my ($key, $file) = @_;
    $self->{main_cm_file_path}->{$key} = $file;
}

sub mainNS_cm_file_path
{
    my $self = shift;
    if (@_) { %{$self->{mainNS_cm_file_path}} = shift; }
    return %{$self->{mainNS_cm_file_path}};
}

sub add_mainNS_cm_file_path
{
    my $self = shift;
    my ($key, $file) = @_;
    $self->{mainNS_cm_file_path}->{$key} = $file;
}

sub cove_cm_file_path
{
    my $self = shift;
    if (@_) { $self->{cove_cm_file_path} = shift; }
    return $self->{cove_cm_file_path};
}

sub isotype_cm_db_file_path
{
    my $self = shift;
    if (@_) { $self->{isotype_cm_db_file_path} = shift; }
    return $self->{isotype_cm_db_file_path};
}

sub mito_isotype_cm_db_file_path
{
    my $self = shift;
    if (@_) { $self->{mito_isotype_cm_db_file_path} = shift; }
    return $self->{mito_isotype_cm_db_file_path};
}

sub isotype_cm_file_paths
{
    my $self = shift;
    if (@_) { %{$self->{isotype_cm_file_paths}} = shift; }
    return %{$self->{isotype_cm_file_paths}};
}

sub add_isotype_cm_file_path
{
    my $self = shift;
    my ($isotype, $file) = @_;
    $self->{isotype_cm_file_paths}->{$isotype} = $file;
}

sub mito_isotype_cm_file_paths
{
    my $self = shift;
    if (@_) { %{$self->{mito_isotype_cm_file_paths}} = shift; }
    return %{$self->{mito_isotype_cm_file_paths}};
}

sub add_mito_isotype_cm_file_path
{
    my $self = shift;
    my ($isotype, $file) = @_;
    $self->{mito_isotype_cm_file_paths}->{$isotype} = $file;
}

sub intron_cm_file_path
{
    my $self = shift;
    if (@_) { %{$self->{intron_cm_file_path}} = shift; }
    return %{$self->{intron_cm_file_path}};
}

sub add_intron_cm_file_path
{
    my $self = shift;
    my ($key, $file) = @_;
    $self->{intron_cm_file_path}->{$key} = $file;
}

sub Pselc_cm_file_path
{
    my $self = shift;
    if (@_) { $self->{Pselc_cm_file_path} = shift; }
    return $self->{Pselc_cm_file_path};
}

sub Eselc_cm_file_path
{
    my $self = shift;
    if (@_) { $self->{Eselc_cm_file_path} = shift; }
    return $self->{Eselc_cm_file_path};
}

sub add_isotype_cm_cutoff
{
    my $self = shift;
    my ($type, $cutoff) = @_;
    $self->{isotype_cm_cutoff}->{$type} = $cutoff;
}

sub covels_bin
{
    my $self = shift;
    if (@_) { $self->{covels_bin} = shift; }
    return $self->{covels_bin};
}

sub coves_bin
{
    my $self = shift;
    if (@_) { $self->{coves_bin} = shift; }
    return $self->{coves_bin};
}

sub cmsearch_bin
{
    my $self = shift;
    if (@_) { $self->{cmsearch_bin} = shift; }
    return $self->{cmsearch_bin};
}

sub infernal_thread
{
    my $self = shift;
    if (@_) { $self->{infernal_thread} = shift; }
    return $self->{infernal_thread};    
}

sub tab_results
{
    my $self = shift;
    if (@_) { $self->{tab_results} = shift; }
    return $self->{tab_results};
}

sub set_defaults
{
    my $self = shift;
    my $global_vars = shift;
    my $global_constants = $global_vars->{global_constants};
    
    $self->{cm_mode} = $global_constants->get("cm_mode");
    $self->{cm_cutoff} = $global_constants->get("cm_cutoff");
    $self->{organelle_cm_cutoff} = $global_constants->get("organelle_cm_cutoff");
    $self->{infernal_fp_cutoff} = $global_constants->get("infernal_fp_cutoff");
    $self->{max_tRNA_length} = $global_constants->get("max_tRNA_length");
    $self->{max_cove_tRNA_length} = $global_constants->get("max_cove_tRNA_length");
    $self->{max_cmsearch_tRNA_length} = $global_constants->get("max_cmsearch_tRNA_length");
    $self->{min_tRNA_no_intron} = $global_constants->get("min_tRNA_no_intron");
    $self->{min_intron_length} = $global_constants->get("min_intron_length");
    $self->{min_cove_pseudo_filter_score} = $global_constants->get("min_cove_pseudo_filter_score");
    $self->{min_cmsearch_pseudo_filter_score} = $global_constants->get("min_cmsearch_pseudo_filter_score");
    $self->{min_ss_score} = $global_constants->get("min_ss_score");
    $self->{min_hmm_score} = $global_constants->get("min_hmm_score");
    $self->{nci_scan_cutoff} = $global_constants->get("nci_scan_cutoff");
    $self->{BHB_cm_cutoff} = $global_constants->get("BHB_cm_cutoff");
    $self->{split_tRNA_scan_cutoff} = $global_constants->get("split_tRNA_scan_cutoff");
    $self->{half_tRNA_cutoff} = $global_constants->get("half_tRNA_cutoff");
    $self->{left_splicing_len} = $global_constants->get("left_splicing_len");
    $self->{right_splicing_len} = $global_constants->get("right_splicing_len");
    my $search_types = $global_constants->get("isotype_cm_cutoff");
    foreach my $type (sort keys %$search_types)
    {
        $self->add_isotype_cm_cutoff($type, $search_types->{$type});
    }
}

sub set_file_paths
{    
    my $self = shift;
    my $global_vars = shift;
    my $opts = $global_vars->{options};
    my $global_constants = $global_vars->{global_constants};
    my $isotype_cms = "";
    
    if ($opts->euk_mode())
    {
        if ($self->infernal_mode())
        {
            $self->{main_cm_file_path}->{Domain} = $global_constants->get_subvalue("cm", "eukaryota");     # default to eukar model 
            $self->{main_cm_file_path}->{SeC} = $global_constants->get_subvalue("euk_cm", "SeC");          # Euk SeC model
            $self->{mainNS_cm_file_path}->{Domain} = $global_constants->get_subvalue("cm", "eukaryota-ns");          # no secondary struct
            $self->{cove_cm_file_path} = $global_constants->get_subvalue("cove_cm", "eukaryota");          # default to eukar cove model 
        }
        elsif ($self->cove_mode())
        {
            $self->{main_cm_file_path}->{Domain} = $global_constants->get_subvalue("cove_cm", "eukaryota");          # default to eukar cove model 
            $self->{mainNS_cm_file_path}->{Domain} = $global_constants->get_subvalue("cove_cm", "eukaryota-ns");     # no secondary struct
        }
        $self->{isotype_cm_db_file_path} = $global_constants->get_subvalue("isotype_cm", "eukaryota");
        if ($opts->mito_model() ne "")
        {
            $self->{mito_isotype_cm_db_file_path} = $global_constants->get_subvalue("mito_cm", $opts->mito_model());
        }
    }                           
    elsif ($opts->bact_mode())
    {
        if ($self->infernal_mode())
        {
            $self->{main_cm_file_path}->{Domain} = $global_constants->get_subvalue("cm", "bacteria");               # use bacterial covariance model 
            $self->{main_cm_file_path}->{SeC} = $global_constants->get_subvalue("bact_cm", "SeC");        # Bacterial SeC model
            $self->{mainNS_cm_file_path}->{Domain} = $global_constants->get_subvalue("cm", "bacteria-ns");          # no sec struct
            $self->{cove_cm_file_path} = $global_constants->get_subvalue("cove_cm", "bacteria");          # use bacterial covariance model 
        }
        elsif ($self->cove_mode())
        {
            $self->{main_cm_file_path}->{Domain} = $global_constants->get_subvalue("cove_cm", "bacteria");          # use bacterial covariance model 
            $self->{mainNS_cm_file_path}->{Domain} = $global_constants->get_subvalue("cove_cm", "bacteria-ns");     # no sec struct
        }
        $self->{isotype_cm_db_file_path} = $global_constants->get_subvalue("isotype_cm", "bacteria");
    }
    elsif ($opts->arch_mode())
    {
        my $nci_cms = $global_constants->get("nci_cm");
        foreach my $type (sort keys %$nci_cms)
        {
            $self->add_intron_cm_file_path($type, $nci_cms->{$type});
        }
        $self->{arch_five_half_cm_file_path} = $global_constants->get_subvalue("cm", "arch_5h");          # model for finding 5'half
        $self->{arch_three_half_cm_file_path} = $global_constants->get_subvalue("cm", "arch_3h");         # model for finding 3'half
        if ($self->infernal_mode())
        {
            $self->{main_cm_file_path}->{Domain} = $global_constants->get_subvalue("cm", "archaea");                # use archaea covariance model 
            $self->{main_cm_file_path}->{SeC} = $global_constants->get_subvalue("arch_cm", "SeC");        # Archaeal SeC model
            $self->{mainNS_cm_file_path}->{Domain} = $global_constants->get_subvalue("cm", "archaea-ns");           # no sec struct
            $self->{cove_cm_file_path} = $global_constants->get_subvalue("cove_cm", "archaea");           # use archaea covariance model 
        }
        elsif ($opts->cove_mode())
        {
            $self->{main_cm_file_path}->{Domain} = $global_constants->get_subvalue("cove_cm", "archaea");           # use archaea covariance model 
            $self->{mainNS_cm_file_path}->{Domain} = $global_constants->get_subvalue("cove_cm", "archaea-ns");      # no sec struct
        }
        $self->{isotype_cm_db_file_path} = $global_constants->get_subvalue("isotype_cm", "archaea");
    }
    elsif ($opts->mito_mode())
    {
        if ($self->infernal_mode())
        {
            $isotype_cms = $global_constants->get("mito_cm_".$opts->mito_model());
            foreach my $isotype (sort keys %$isotype_cms)
            {
                $self->{main_cm_file_path}->{$isotype} = $isotype_cms->{$isotype};
            }
        }
    }
    elsif ($opts->metagenome_mode())
    {
        if ($self->infernal_mode())
        {
            $self->{main_cm_file_path}->{general} = $global_constants->get_subvalue("cm", "general");               # use general covariance model 
            $self->{mainNS_cm_file_path}->{general} = $global_constants->get_subvalue("cm", "general-ns");          # no sec struct
            $self->{main_cm_file_path}->{euk} = $global_constants->get_subvalue("cm", "eukayota");                # use eukayota covariance model 
            $self->{mainNS_cm_file_path}->{euk} = $global_constants->get_subvalue("cm", "eukayota-ns");           # no sec struct
            $self->{main_cm_file_path}->{arch} = $global_constants->get_subvalue("cm", "archaea");                  # use archaea covariance model 
            $self->{mainNS_cm_file_path}->{arch} = $global_constants->get_subvalue("cm", "archaea-ns");             # no sec struct
            $self->{main_cm_file_path}->{bact} = $global_constants->get_subvalue("cm", "bacteria");                # use bacteria covariance model 
            $self->{mainNS_cm_file_path}->{bact} = $global_constants->get_subvalue("cm", "bacteria-ns");           # no sec struct
        }
    }
    elsif ($opts->numt_mode())
    {
        $isotype_cms = $global_constants->get("mito_cm_mammal");
        foreach my $isotype (sort keys %$isotype_cms)
        {
            $self->add_isotype_cm_file_path($isotype, $isotype_cms->{$isotype});
        }
    }
    elsif ($opts->general_mode())
    {
        if ($self->infernal_mode())
        {
            $self->{main_cm_file_path}->{general} = $global_constants->get_subvalue("cm", "general");                # use original covariance model 
            $self->{mainNS_cm_file_path}->{general} = $global_constants->get_subvalue("cm", "general-ns");           # no sec struct
            $self->{cove_cm_file_path} = $global_constants->get_subvalue("cove_cm", "general");           # use original covariance model 
        }
        elsif ($self->cove_mode())
        {
            $self->{main_cm_file_path}->{general} = $global_constants->get_subvalue("cove_cm", "general");           # use original covariance model 
            $self->{mainNS_cm_file_path}->{general} = $global_constants->get_subvalue("cove_cm", "general-ns");      # no sec struct
        }
    }
    elsif ($opts->organelle_mode())
    {
        if ($self->infernal_mode())
        {
            $self->{main_cm_file_path}->{general} = $global_constants->get_subvalue("cm", "general1415");                # use original covariance model 
            $self->{mainNS_cm_file_path}->{general} = $global_constants->get_subvalue("cm", "general1415-ns");           # no sec struct
            $self->{cove_cm_file_path} = $global_constants->get_subvalue("cove_cm", "general");           # use original covariance model 
        }
        elsif ($self->cove_mode())
        {
            $self->{main_cm_file_path}->{general} = $global_constants->get_subvalue("cove_cm", "general");           # use original covariance model 
            $self->{mainNS_cm_file_path}->{general} = $global_constants->get_subvalue("cove_cm", "general-ns");      # no sec struct
        }
    }
    elsif ($opts->alternate_mode())
    {
        my $cms = $global_constants->get("alt_cm");
        foreach my $key (sort keys %$cms)
        {
            $self->add_main_cm_file_path($key, $cms->{$key});
        }
    }
        
    if ($self->infernal_mode())
    {
        $self->{Pselc_cm_file_path} = $self->{main_cm_file_path}->{SeC};
        $self->{Eselc_cm_file_path} = $self->{main_cm_file_path}->{SeC};
    }
    elsif ($self->cove_mode())
    {
        $self->{Pselc_cm_file_path} = $global_constants->get_subvalue("cove_cm", "PSELC");
        $self->{Eselc_cm_file_path} = $global_constants->get_subvalue("cove_cm", "ESELC");
    }
}

sub check_lib_files
{    
    my $self = shift;
    my ($opts) = @_;
    
    foreach my $cur_cm (sort keys %{$self->{main_cm_file_path}})
    {
        if ($self->{main_cm_file_path}->{$cur_cm} ne "" and !-r $self->{main_cm_file_path}->{$cur_cm})
        {
            die "FATAL: Unable to open ".$self->{main_cm_file_path}->{$cur_cm}." covariance model file\n\n";
        }
    }

    foreach my $cur_cm (sort keys %{$self->{mainNS_cm_file_path}})
    {
        if ($self->{mainNS_cm_file_path}->{$cur_cm} ne "" and !-r $self->{mainNS_cm_file_path}->{$cur_cm})
        {
            die "FATAL: Unable to open ".$self->{mainNS_cm_file_path}->{$cur_cm}." covariance model file\n\n";
        }
    }
    if ($self->{cove_cm_file_path} ne "" and !-r $self->{cove_cm_file_path})
    {
        die "FATAL: Unable to open ".$self->{cove_cm_file_path}." covariance model file\n\n";
    }

    if ($self->{Pselc_cm_file_path} ne "" and !-r $self->{Pselc_cm_file_path})
    {
        die "FATAL: Unable to open ".$self->{Pselc_cm_file_path}." covariance model file\n\n";
    }

    if ($self->{Eselc_cm_file_path} ne "" and !-r $self->{Eselc_cm_file_path})
    {
        die "FATAL: Unable to open ".$self->{Eselc_cm_file_path}." covariance model file\n\n";
    }

    if ($self->{isotype_cm_db_file_path} ne "" and (!-r $self->{isotype_cm_db_file_path} || !-e "$self->{isotype_cm_db_file_path}".".i1f"))
    {
        die "FATAL: Unable to open ".$self->{isotype_cm_db_file_path}." covariance model file\n\n";
    }

    if ($self->{mito_isotype_cm_db_file_path} ne "" and (!-r $self->{mito_isotype_cm_db_file_path} || !-e "$self->{mito_isotype_cm_db_file_path}".".i1f"))
    {
        die "FATAL: Unable to open ".$self->{mito_isotype_cm_db_file_path}." covariance model file\n\n";
    }

    foreach my $isotype (sort keys %{$self->{isotype_cm_file_paths}})
    {
        if ($self->{isotype_cm_file_paths}->{$isotype} ne "" and !-r $self->{isotype_cm_file_paths}->{$isotype})
        {
            die "FATAL: Unable to open ".$self->{isotype_cm_file_paths}->{$isotype}." covariance model file\n\n";
        }
    }
    
    if ($opts->arch_mode() && ($self->infernal_mode() ||  $self->cove_mode()))
    {
        foreach my $cur_cm (sort keys %{$self->{intron_cm_file_path}})
        {
            if ($self->{intron_cm_file_path}->{$cur_cm} ne "" and !-r $self->{intron_cm_file_path}->{$cur_cm})
            {
                die "FATAL: Unable to open ".$self->{intron_cm_file_path}->{$cur_cm}." covariance model file\n\n";
            }
        }
        if ($self->{arch_five_half_cm_file_path} ne "" and !-r $self->{arch_five_half_cm_file_path})
        {
            die "FATAL: Unable to open ".$self->{arch_five_half_cm_file_path}." covariance model file\n\n";
        }
        if ($self->{arch_three_half_cm_file_path} ne "" and !-r $self->{arch_three_half_cm_file_path})
        {
            die "FATAL: Unable to open ".$self->{arch_three_half_cm_file_path}." covariance model file\n\n";
        }
    }
}

sub set_bin
{    
    my $self = shift;
    my $bindir = shift;
    
    if ($^O =~ /^MSWin/)
    {
        $self->{cmsearch_bin} .= ".exe";
        $self->{covels_bin} .= ".exe";
        $self->{coves_bin} .= ".exe";
    }
    if (!(-x $self->{covels_bin}))
    {
        $self->{covels_bin} = $bindir."/".$self->{covels_bin};
        if (!(-x $self->{covels_bin}))
        {
            die "FATAL: Unable to find ".$self->{covels_bin}." executable\n\n";
        }
    }
    if (!(-x $self->{coves_bin}))
    {
        $self->{coves_bin} = $bindir."/".$self->{coves_bin};
        if (!(-x $self->{coves_bin}))
        {
            die "FATAL: Unable to find ".$self->{coves_bin}." executable\n\n";
        }
    }
}

sub set_infernal_bin
{    
    my $self = shift;
    my $bindir = shift;
    
    $self->{cmsearch_bin} = $bindir."/".$self->{cmsearch_bin};
    if (!(-x $self->{cmsearch_bin}))
    {
        die "FATAL: Unable to find ".$self->{cmsearch_bin}." executable\n\n";
    }
    else
    {
        $self->check_infernal_version($self->{cmsearch_bin});
    }

    $self->{cmscan_bin} = $bindir."/".$self->{cmscan_bin};
    if (!(-x $self->{cmscan_bin}))
    {
        die "FATAL: Unable to find ".$self->{cmscan_bin}." executable\n\n";
    }
    else
    {
        $self->check_infernal_version($self->{cmscan_bin});
    }
}

sub check_infernal_version
{
    my $self = shift;
    my $exec_bin = shift;

    my $cmd = $exec_bin." -h | grep INFERNAL";
    my $line = `$cmd`;
    my $version = 0;
    if ($line =~ /INFERNAL 1.1.(\d+) /)
    {
        $version = $1;
    }
    if ($version < 2)
    {
        die "FATAL: ".$exec_bin." is incompatible with tRNAscan-SE. Please install Infernal version 1.1.2 or above.\n\n";
    } 
}

sub set_search_params
{
    my $self = shift;
    my ($opts, $r_scan_len, $r_cur_cm_file, $max_search_tRNA_length, $trna) = @_;

    # don't set '-W' param over 200 bp if a pre-scanner is being used,
    #  use max window of 150 bp if cmsearch only (too slow otherwise)    
    if ($opts->eufind_mode() || $opts->tscan_mode() || $opts->use_prev_ts_run())
    {
        $$r_scan_len = &min($trna->len(), $self->{max_tRNA_length});
    }
    else
    {
        $$r_scan_len = $max_search_tRNA_length;
    }        

    # set correct CM file for current tRNA
    $$r_cur_cm_file = $self->{main_cm_file_path}->{Domain};
    if ($opts->general_mode() or $opts->organelle_mode())
    {
        $$r_cur_cm_file = $self->{main_cm_file_path}->{general};
    }
    
    if ($opts->eufind_mode())
    {
        # use prok selcys model
        if ($trna->isotype() eq "SeCp")
        {       
            $$r_cur_cm_file = $self->{Pselc_cm_file_path};
        }
        # use euk selcys model
        elsif  ($trna->isotype() eq "SeCe")
        {    
            $$r_cur_cm_file = $self->{Eselc_cm_file_path};
        }            
    }
}

# find anticodon loop & a-codon from coves or cmsearch output
sub find_anticodon
{                
    my $self = shift;
    my ($opts, $trna, $gc) = @_;
    my ($antiloop, $antiloop_end, $ac_index, $anticodon, $verify_ac);

    my $seq = $trna->seq();
    my $ss = $trna->ss();
    my $undef_anticodon = $gc->undef_anticodon();
    
    my $antiloop_index = 0;
    my $antiloop_len = 0;

    # Match pattern in secondary structure output, 
    # looking for second stem-loop structure ">>>>...<<<<"
    # that should be the anitocodon stem-loop 

    if ($ss =~ /^([>.]+<[<.]+>[>.]*)>([.]{4,})<+.+[>.]+<[<.]+/o)
    {
        # set to index position of first base in anticodon loop
        $antiloop_index = length($1) + 1;
        $antiloop_len = length($2);   # anticodon loop length
    }
    
    if ($antiloop_index != 0 and $antiloop_len != 0)
    {
        # index of end of anticodon loop
        $antiloop_end = $antiloop_index + $antiloop_len - 1;
    
        $antiloop = substr($seq, $antiloop_index, $antiloop_len);
    
        # remove '-' gaps from loop
        $antiloop =~ s/[\-]//g;      

        # Don't guess if even number of bp in 
        # anticodon loop

        # remove introns & non-canonical bases
        $antiloop =~ s/[a-z]//g;      

        if ((length($antiloop) < 5) || ((length($antiloop) % 2) == 0))
        {
            return ($undef_anticodon, -1, -1, -1);
        }
        # get anticodon 
        $ac_index = (length($antiloop) - 3) / 2;
        $anticodon = substr($antiloop, $ac_index, 3);
        $verify_ac = substr($seq, $ac_index + $antiloop_index, 3);
            
        # check to see if anticodon extracted from the entire
        #  trna sequence (coveseq) is same as that extracted from
        #  just the anticodon loop sequence (antiloop)
    
        if ($verify_ac ne $anticodon)
        {
            $trna->category("undetermined_ac");
            return ($undef_anticodon, -1, -1, -1);            
        }
        return ($anticodon, $antiloop_index, $antiloop_end, $ac_index + $antiloop_index + 1);
    }
    else
    {
        $trna->category("undetermined_ac");
        return ($undef_anticodon, -1, -1, -1);
    }
}

# find anticodon loop & a-codon from coves or cmsearch output
sub find_mito_anticodon
{                
    my $self = shift;
    my ($opts, $trna, $gc) = @_;
    my ($antiloop, $antiloop_end, $ac_index, $anticodon, $verify_ac);

    my $seq = $trna->seq();
    my $ss = $trna->ss();
    my $model = $trna->model();
    my $undef_anticodon = $gc->undef_anticodon();
    
    my $antiloop_index = 0;
    my $antiloop_len = 0;

    # Match pattern in secondary structure output, 
    # looking for second stem-loop structure ">>>>...<<<<"
    # that should be the anitocodon stem-loop 

    if (($model eq "SerGCT" or $model eq "Cys_NoDarm") and $ss =~ /^([>.]+)>([.]{4,})<+.+[>.]+<[<.]+/o)
    {
        # set to index position of first base in anticodon loop
        $antiloop_index = length($1) + 1;
        $antiloop_len = length($2);   # anticodon loop length
        $trna->note("No D-arm");
    }
    elsif ($ss =~ /^([>.]+<[<.]+[>.]*)>([.]{4,})<[<.]+[.]{4,}<[<.]+$/o)
    {
        # set to index position of first base in anticodon loop
        $antiloop_index = length($1) + 1;
        $antiloop_len = length($2);   # anticodon loop length          
        $trna->note("No T-arm");
    }
    elsif ($ss =~ /^([>.]+[.]{4,}[>.]+)>([.]{4,})<[<.]+\.+[>.]+<[<.]+$/o)
    {
        # set to index position of first base in anticodon loop
        $antiloop_index = length($1) + 1;
        $antiloop_len = length($2);   # anticodon loop length
        $trna->note("No D-arm");
    }
    elsif ($ss =~ /^([>.]+<[<.]+[>.]*)>([.]{4,})<+.+[>.]+<[<.]+/o)
    {
        # set to index position of first base in anticodon loop
        $antiloop_index = length($1) + 1;
        $antiloop_len = length($2);   # anticodon loop length          
    }
    
    if ($antiloop_index != 0 and $antiloop_len != 0)
    {
        # index of end of anticodon loop
        $antiloop_end = $antiloop_index + $antiloop_len - 1;
    
        $antiloop = substr($seq, $antiloop_index, $antiloop_len);
    
        # remove '-' gaps from loop
        $antiloop =~ s/[\-]//g;      

        # Don't guess if even number of bp in 
        # anticodon loop
        $antiloop = uc($antiloop);
        if (length($antiloop) < 5)
        {
            $trna->category("undetermined_ac");
            return ($undef_anticodon, -1, -1, -1);
        }
        elsif ((length($antiloop) % 2) == 0)
        {
            my $found = 0;
            my $j = 0;
            for (my $i = int((length($antiloop) - 3) / 2); $i <= (length($antiloop) - 3) and $i >= 0; $i = $i + $j)
            {
                $ac_index = $i;
                $anticodon = substr($antiloop, $ac_index, 3);
                my $isotype = $gc->get_tRNA_type($self, $anticodon, $self->{main_cm_file_path}->{$model}, $model, $self->cove_mode());
                if ($model eq $isotype or ($model eq "SerGCT" and $isotype eq "Ser") or ($model eq "SerTGA" and $isotype eq "Ser")
                    or ($model eq "LeuTAG" and $isotype eq "Leu") or ($model eq "LeuTAA" and $isotype eq "Leu") 
                    or ($model eq "Cys_NoDarm" and $isotype eq "Cys"))
                {
                    $trna->category("mito_ac_mislocation");
                    $found = 1;
                    last;
                }
                $j = abs($j);
                $j++;
                if ($j % 2 == 0)
                {
                    $j = $j * -1;
                }
            }

            if ($found)
            {
                $verify_ac = uc(substr($seq, $ac_index + $antiloop_index, 3));
            }
            else
            {
                $trna->category("undetermined_ac");
                return ($undef_anticodon, -1, -1, -1);
            }
        }
        else
        {
            # get anticodon 
            $ac_index = (length($antiloop) - 3) / 2;
            $anticodon = substr($antiloop, $ac_index, 3);
            $verify_ac = uc(substr($seq, $ac_index + $antiloop_index, 3));
            my $isotype = $gc->get_tRNA_type($self, $anticodon, $self->{main_cm_file_path}->{$model}, $model, $self->cove_mode());
            my $model_iso = $model;
            if (length($model) > 3)
            {
                $model_iso = substr($model, 0, 3);
            }
            if ($isotype ne $model_iso)
            {
                $verify_ac = "";
                my $found = 0;
                my $j = 1;
                for (my $i = (length($antiloop) - 3) / 2 - 1; $i <= (length($antiloop) - 3) and $i >= 0; $i = $i + $j)
                {
                    $ac_index = $i;
                    $anticodon = substr($antiloop, $ac_index, 3);
                    my $isotype = $gc->get_tRNA_type($self, $anticodon, $self->{main_cm_file_path}->{$model}, $model, $self->cove_mode());
                    if ($model eq $isotype or ($model eq "SerGCT" and $isotype eq "Ser") or ($model eq "SerTGA" and $isotype eq "Ser")
                        or ($model eq "LeuTAG" and $isotype eq "Leu") or ($model eq "LeuTAA" and $isotype eq "Leu") 
                        or ($model eq "Cys_NoDarm" and $isotype eq "Cys"))
                    {
                        $trna->category("mito_ac_mislocation");
                        $found = 1;
                        last;
                    }
                    $j = abs($j);
                    $j++;
                    if ($j % 2 == 1)
                    {
                        $j = $j * -1;
                    }
                }
                if ($found)
                {
                    $verify_ac = uc(substr($seq, $ac_index + $antiloop_index, 3));
                }
            }
        }
            
        # check to see if anticodon extracted from the entire
        #  trna sequence (coveseq) is same as that extracted from
        #  just the anticodon loop sequence (antiloop)
    
        if ($verify_ac ne $anticodon)
        {
            $trna->category("undetermined_ac");
            return ($undef_anticodon, -1, -1, -1);            
        }
        return ($anticodon, $antiloop_index, $antiloop_end, $ac_index + $antiloop_index + 1);
    }
    else
    {
        $trna->category("undetermined_ac");
        return ($undef_anticodon, -1, -1, -1);
    }
}

sub find_intron
{
    my $self = shift;
    my ($trna_seq, $antiloop_index, $antiloop_end) = @_;
    my ($intron, $istart, $iend, $tmpstr, $antiloop_seq);
    my $min_intron_length = $self->{min_intron_length};

    # check to see if it was unable 
    # to determine the anticodon loop
    if ($antiloop_index == -1)
    {
        return(0, 0, 0);
    }
    # get subsequence from start of anticodon loop
    # to end of anticodon loop -- look for intron in it
    $antiloop_seq = substr($trna_seq, $antiloop_index, $antiloop_end - $antiloop_index + 1);
    
    if ($antiloop_seq =~ /^(.*[^a-z]+)([a-z]{$min_intron_length,})[^a-z]+/o)
    {
        $intron = $2;
    
        # make sure to get the base index for the last (not nec. only) occurrence
        # of the intron sequence string up to end of anticodon loop
        $tmpstr = substr($trna_seq, 0, $antiloop_end+1);
        $istart = index($tmpstr, $intron) + 1; 
        $iend = length($intron) + $istart - 1;
    }
    else
    {
        $intron = 0; 
    }
    return ($intron, $istart, $iend);
}                        

# is_pseudo_gene
#
# Runs a covariance model without secondary structure information on predicted tRNA, puts this value
# in "hmm_score". Contribution to total score from secondary structure derived by subtracting hmm_score from total score
# Returns non-zero if tRNA scores fall below minima for either primary or secondary structure components of score
sub is_pseudo_gene
{    
    my $self = shift;
    my $global_constants = shift;
    my $opts = shift;
    my $log = shift;
    my ($trna, $prescan_tRNA_id, $tmp_trnaseq_file) = @_;
    my ($dummy1, $dummy2, $hmm_score, $hit_start, $hit_end, $hit_ct, $min_pseudo_filter_score);

    $dummy1 = $dummy2 = "";                 # return values not used
    
    # skip check for pseudo gene if score is above minimum
    # -D (disable pseudogene checking) is specified 
    # AND -H option (get hmm scores) is NOT specified
    if ($self->cove_mode())
    {
        $min_pseudo_filter_score = $self->{min_cove_pseudo_filter_score};
    }
    elsif ($self->infernal_mode())
    {
        $min_pseudo_filter_score = $self->{min_cmsearch_pseudo_filter_score};
    }
    
    if ((($trna->score() >= $min_pseudo_filter_score) || $self->{skip_pseudo_filter}) && !$self->{get_hmm_score})
    {
        return 0;
    }

    $hmm_score = 0;
    my $ss_score = 0;
    my $score = 0;
    my $ns_cm = "";
    if ($self->cove_mode())
    {
        $ns_cm = $self->{mainNS_cm_file_path}->{Domain};
        if ($opts->general_mode() or $opts->organelle_mode())
        {
            $ns_cm = $self->{mainNS_cm_file_path}->{general}
        }
        
        ($dummy1, $dummy2, $hmm_score) = 
            $self->run_coves($log, $tmp_trnaseq_file, $prescan_tRNA_id, $ns_cm);
        $score = $trna->get_domain_model("cove")->{score};
        $ss_score = $score - $hmm_score;
        $trna->update_domain_model("cove", $score, $score, $hmm_score, $ss_score);
    }
    elsif ($self->infernal_mode()) 
    {
        $ns_cm = $self->{mainNS_cm_file_path}->{$trna->model()};
        my $cms_output_file = $global_constants->get("temp_dir")."/tscan$$"."_cm_hmm.out";

       ($hmm_score, $hit_start, $hit_end, $hit_ct) = 
            $self->cmsearch_scoring($opts, $log, $tmp_trnaseq_file, $prescan_tRNA_id, $ns_cm, $cms_output_file);
        $score = $trna->get_domain_model("infernal")->{score};
        $ss_score = $score - $hmm_score;
        $trna->update_domain_model("infernal", $score, $score, $hmm_score, $ss_score);
    }
    else
    {
        return -1;                              # Error - no second pass scanner selected
    }

    if ((($ss_score < $self->{min_ss_score}) || ($hmm_score < $self->{min_hmm_score})) &&
        ($score < $min_pseudo_filter_score))
    {
        return 1;
    }
}    

# Get anticodon, isotype, intron location, and pseudogene status for tRNA with noncanonical introns
sub decode_nci_tRNA_properties
{
    my $self = shift;
    my ($global_vars, $trna) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $opts = $global_vars->{options};
    my $gc = $global_vars->{gc};
    my $log = $global_vars->{log_file};
    
    my ($anticodon, $acodon_index, $acodon_end1, $acodon_index2, $trna_type, $intron, $istart, $iend,
          $antiloop_index, $antiloop_end, $trna_len, $scan_len);

    $anticodon = "ERR";

    my $mat_seq = "";
    my $precursor_seq = $trna->seq();
    my @introns = $trna->ar_introns();
    for (my $i = 0; $i < scalar(@introns); $i++)
    {
        if ($i == 0)
        {
            $mat_seq = substr($precursor_seq, 0, $introns[$i]->{rel_start} - 1);
        }
        else
        {
            $mat_seq .= substr($precursor_seq, $introns[$i-1]->{rel_end}, $introns[$i]->{rel_start} - $introns[$i-1]->{rel_end} - 1);
        }
    }
    $mat_seq .= substr($precursor_seq, $introns[scalar(@introns)-1]->{rel_end});

    if (uc($mat_seq) eq uc($trna->mat_seq()))
    {
		$trna->seq($trna->mat_seq());
	}
	else
    {
        $log->warning("Mismatch mature sequence for tRNA for noncanonical intron\t". $trna->tRNAscan_id());
    }
    
    if ($self->cove_mode() || $self->infernal_mode())
    {
        ($anticodon, $antiloop_index, $antiloop_end, $acodon_index) = $self->find_anticodon($opts, $trna, $gc);
        $acodon_index2 = 0;
        $trna->anticodon($anticodon);
        
        for (my $i = 0; $i < scalar(@introns); $i++)
        {
            if ($acodon_index >= $introns[$i]->{rel_start})
            {
				$acodon_index += ($introns[$i]->{rel_end} - $introns[$i]->{rel_start} + 1);
			}
			elsif ($acodon_index < $introns[$i]->{rel_start} and ($acodon_index + 2) >= $introns[$i]->{rel_start})
            {
                $acodon_end1 = $introns[$i]->{rel_start} - 1;
                $acodon_index2 = $introns[$i]->{rel_end} + 1; 
            }
            elsif ($acodon_index + 2 < $introns[$i]->{rel_start})
            {
                last;
            }
        }
        if ($trna->get_ac_pos_count() > 0)
        {
            $trna->remove_ac_pos(0);
        }
        if ($acodon_index2 == 0)
        {
            $trna->add_ac_pos($acodon_index, $acodon_index + 2);
        }
        else
        {
            $trna->add_ac_pos($acodon_index, $acodon_end1);
            $trna->add_ac_pos($acodon_index2, (3 - ($acodon_end1 - $acodon_index + 1)) + $acodon_index2 - 1);
        }
    }
    else
    {
        die "Second pass mode not selected -- can't decode tRNA type\n\n";
    }
    
    # check for problem parsing anticodon loop 
    if (($anticodon eq $gc->undef_anticodon()) || ($trna->seq()  eq 'Error'))    
    {
        $trna->anticodon($gc->undef_anticodon());
        $trna->isotype($gc->undef_isotype());
        $intron = 0;        
        
        if ($opts->save_odd_struct())
        {     
            open(ODDTRNA, ">>".$opts->odd_struct_file()) ||
                die "FATAL: Can't open ".$opts->odd_struct_file()." to save seconary structures\n\n"; 
            print ODDTRNA $trna->tRNAscan_id()." (".$trna->start()."-".$trna->end()."):\n".$trna->seq()."\n".$trna->ss()."\n\n"; 
            close(ODDTRNA);
        }
    }
    # continue tRNA struct parsing
    else
    {                               
        ($intron, $istart, $iend) = $self->find_intron($trna->seq(), $antiloop_index, $antiloop_end);

        if ($intron)
        {
            $mat_seq = substr($trna->seq(), 0, $istart - 1).substr($trna->seq(), $iend);
            $trna->mat_ss(substr($trna->ss(), 0, $istart - 1).substr($trna->ss(), $iend));
            $trna->ss($trna->mat_ss());
            
            for (my $i = 0; $i < scalar(@introns); $i++)
            {
                if ($istart > $introns[$i]->{rel_start})
                {
                    $istart += ($introns[$i]->{rel_end} - $introns[$i]->{rel_start} + 1);
                    $iend += ($introns[$i]->{rel_end} - $introns[$i]->{rel_start} + 1);
                }
                elsif ($iend < $introns[$i]->{rel_start})
                {
                    last;
                }
            }
            my $intron_start = 0;
            my $intron_end = 0;
            if ($trna->strand() eq "+")
            {
                $intron_start = $istart + $trna->start() - 1;
                $intron_end = $iend + $trna->start() - 1
            }
            else
            {
                $intron_start = $trna->end() - $iend + 1;
                $intron_end = $trna->end() - $istart + 1;
            }
            $trna->add_intron($istart, $iend, $intron_start, $intron_end, "CI", $intron);
            $trna->merge_introns();
        }
                
        $trna->isotype($gc->get_tRNA_type($self, $trna->anticodon(), $self->{main_cm_file_path}->{$trna->model()}, $trna->model(), $self->cove_mode()));
    }

    # Reset CI type if marked as NCI
    @introns = $trna->ar_introns();
    my $has_CI = 0;
    for (my $i = 0; $i < scalar(@introns); $i++)
    {
        if ($introns[$i]->{type} eq "CI")
        {
			$has_CI = 1;
            last;
		}		
    }
    
    for (my $i = 0; !$has_CI and $i < scalar(@introns); $i++)
    {
        if ($acodon_index2 == 0)
        {
            if ($acodon_index + 2 + 2 == $introns[$i]->{rel_start} and $introns[$i]->{type} eq "NCI")
            {
                $trna->set_intron($i, $introns[$i]->{rel_start}, $introns[$i]->{rel_end}, "CI", $introns[$i]->{seq});
                last;
            }
        }
        else
        {
            if (((3 - ($acodon_end1 - $acodon_index + 1)) + $acodon_index2 - 1 + 2) == $introns[$i]->{rel_start} and $introns[$i]->{type} eq "NCI")
            {
                $trna->set_intron($i, $introns[$i]->{rel_start}, $introns[$i]->{rel_end}, "CI", $introns[$i]->{seq});
                last;
            }
        }        
    }
    
    # Write current tRNA to temp file for re-analysis with other models
    $trna_len = length($trna->seq());
    &write_tRNA($global_constants->get("tmp_trnaseq_file"), $trna->tRNAscan_id(), "", $trna->seq(), 1);
    
    if (($trna->model() ne "SeC") &&
        ($self->is_pseudo_gene($global_constants, $opts, $log, $trna, $trna->tRNAscan_id(), $global_constants->get("tmp_trnaseq_file"))) &&
        (!$self->{skip_pseudo_filter}))
    {
        $trna->pseudo(1);
    }
    else
    {
        $trna->pseudo(0);
    }
    
    $trna->mat_seq($mat_seq);
    $trna->seq($precursor_seq);
}

# Get tRNA anticodon, isotype, intron location, and pseudogene status
sub decode_tRNA_properties
{
    my $self = shift;
    my ($global_vars, $trna, $prescan_tRNA, $cur_cm_file) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $opts = $global_vars->{options};
    my $gc = $global_vars->{gc};
    my $log = $global_vars->{log_file};
    
    my ($anticodon, $acodon_index, $trna_type, $intron, $istart, $iend,
          $antiloop_index, $antiloop_end, $trna_len, $scan_len);

    $anticodon = "ERR";

    if ($self->cove_mode() || $self->infernal_mode())
    {
        ($anticodon, $antiloop_index, $antiloop_end, $acodon_index) = $self->find_anticodon($opts, $trna, $gc);
        $trna->anticodon($anticodon);
        $trna->add_ac_pos($acodon_index, $acodon_index + 2);
    }
    else
    {
        die "Second pass mode not selected -- can't decode tRNA type\n\n";
    }
    
    # check for problem parsing anticodon loop 
    if (($anticodon eq $gc->undef_anticodon()) || ($trna->seq()  eq 'Error'))    
    {
        $trna->anticodon($gc->undef_anticodon());
        $trna->isotype($gc->undef_isotype());
        $intron = 0;        
        
        if ($opts->save_odd_struct())
        {     
            open(ODDTRNA, ">>".$opts->odd_struct_file()) ||
                die "FATAL: Can't open ".$opts->odd_struct_file()." to save seconary structures\n\n"; 
            print ODDTRNA $prescan_tRNA->tRNAscan_id()." (".$trna->start()."-".$trna->end()."):\n".$trna->seq()."\n".$trna->ss()."\n\n"; 
            close(ODDTRNA);
        }
    }
    # continue tRNA struct parsing
    else
    {                               
        ($intron, $istart, $iend) = $self->find_intron($trna->seq(), $antiloop_index, $antiloop_end);
        
        if ($intron)
        {
            my $intron_start = 0;
            my $intron_end = 0;
            if ($trna->strand() eq "+")
            {
                $intron_start = $istart + $trna->start() - 1;
                $intron_end = $iend + $trna->start() - 1
            }
            else
            {
                $intron_start = $trna->end() - $iend + 1;
                $intron_end = $trna->end() - $istart + 1;
            }
            $trna->add_intron($istart, $iend, $intron_start, $intron_end, "CI", $intron);
        }
        
        if (defined $prescan_tRNA->anticodon() and ($prescan_tRNA->anticodon() ne "???"))
        {
            if (($anticodon ne (uc($prescan_tRNA->anticodon()))) && 
                ($opts->tscan_mode() || $opts->eufind_mode()) && ($opts->strict_params()))
            {
                $log->broadcast($prescan_tRNA->tRNAscan_id()." - anticondon conflict\t".$opts->second_pass_label().": $anticodon\tfirstpass (".$prescan_tRNA->hit_source().")".
                    ": ".$prescan_tRNA->anticodon()."\n".$trna->seq()."\n".$trna->ss()."\n"); 
            }
        }
        
        $trna->isotype($gc->get_tRNA_type($self, $trna->anticodon(), $cur_cm_file, $trna->model(), $self->cove_mode()));
        
        if ($anticodon ne "TCA" and $trna->isotype() eq "SeC")
        {
            $trna->anticodon($gc->undef_anticodon());
            $trna->isotype($gc->undef_isotype());
        }
    }
    
    # Write current tRNA to temp file for re-analysis with other models
    $trna_len = length($trna->seq());
    &write_tRNA($global_constants->get("tmp_trnaseq_file"), $prescan_tRNA->tRNAscan_id(), "", $trna->seq(), 1);
    
    if (($trna->model() ne "SeC") &&
        ($self->is_pseudo_gene($global_constants, $opts, $log, $trna, $prescan_tRNA->tRNAscan_id(), $global_constants->get("tmp_trnaseq_file"))) &&
        (!$self->{skip_pseudo_filter}))
    {
        $trna->pseudo(1);
    }   
}

# Fix fMet with missing first base
sub fix_fMet
{
    my $self = shift;
    my ($global_vars, $trna) = @_;
    my $opts = $global_vars->{options};
    my $rescore = 0;
    
    my $inf = $trna->get_domain_model("infernal");
    
    if ($inf->{score} > 40 and $trna->isotype() eq "Met" and $opts->bact_mode())
    {
		if ((substr($trna->seq(), length($trna->seq())-3) eq "CCA" and substr($trna->ss(), length($trna->ss())-5) eq ".....") or
            (substr($trna->ss(), length($trna->ss())-5) eq "<<<.."))
        {
			if (substr($trna->ss(), 0, 1) ne ".")
            {
                if (($trna->strand() eq "+" and $trna->start() > 1) or
                    ($trna->strand() eq "-" and $trna->end() < $trna->src_seqlen()))
                {
                    $trna->seq(substr($trna->upstream(), length($trna->upstream())-1) . $trna->seq());
                    $trna->upstream(substr($trna->upstream(), 0, length($trna->upstream())-1));
                    $trna->ss(".".$trna->ss());
                    if ($trna->strand() eq "+")
                    {
                        $trna->start($trna->start() - 1);
                    }
                    else
                    {
                        $trna->end($trna->end() + 1);
                    }
                    my @ar_ac_pos = $trna->ar_ac_pos();
                    if (scalar(@ar_ac_pos) > 0)
                    {
                        $ar_ac_pos[0]->{rel_start} += 1;
                        $ar_ac_pos[0]->{rel_end} += 1;
                        $trna->ar_ac_pos(@ar_ac_pos);
                    }
                    
                    $rescore = 1;
                }
            }
            elsif (substr($trna->ss(), 0, 4) eq ".>.>" and substr($trna->seq(), 0, 2) eq "CG")
            {
                my $pos71 = "";
                if (substr($trna->seq(), length($trna->seq())-3) eq "CCA" and substr($trna->ss(), length($trna->ss())-5) eq ".....")
                {
					$pos71 = substr($trna->seq(), length($trna->ss())-6, 1);
				}
				elsif (substr($trna->ss(), length($trna->ss())-5) eq "<<<..")
                {
                    $pos71 = substr($trna->seq(), length($trna->ss())-3, 1);
                }
                if (&rev_comp_seq(uc($pos71)) eq uc(substr($trna->seq(), 2, 1)))
                {
					$trna->seq(substr($trna->seq(), 1, 1).uc(substr($trna->seq(), 2, 1)).substr($trna->seq(), 3));
                    $trna->ss(substr($trna->ss(), 0, 2).substr($trna->ss(), 3));
                    if ($trna->strand() eq "+")
                    {
                        $trna->start($trna->start() + 1);
                    }
                    else
                    {
                        $trna->end($trna->end() - 1);
                    }
                    my @ar_ac_pos = $trna->ar_ac_pos();
                    if (scalar(@ar_ac_pos) > 0)
                    {
                        $ar_ac_pos[0]->{rel_start} -= 1;
                        $ar_ac_pos[0]->{rel_end} -= 1;
                        $trna->ar_ac_pos(@ar_ac_pos);
                    }
                    $rescore = 1;
				}				
            }
		}		
	}
    
    if ($rescore)
    {
		$self->rescore_tRNA($global_vars, $trna, $trna);
	}	
}

# Fix archaeal His for incorrect bulge
sub fix_His
{
    my $self = shift;
    my ($global_vars, $trna) = @_;
    my $opts = $global_vars->{options};
    my $rescore = 0;
    
    my $inf = $trna->get_domain_model("infernal");
    
    if ($inf->{score} > 35 and $trna->isotype() eq "His" and $opts->arch_mode())
    {
        my $pos5 = "";
        my $pos68 = "";
        my $mid = "";
		if (substr($trna->ss(), 0, 9) eq ">>>>.>>>." and substr($trna->ss(), length($trna->ss())-9) eq "<<<.<<<<.")
        {
            $pos5 = uc(substr($trna->seq(), 4, 1));
            $pos68 = uc(substr($trna->seq(), length($trna->seq())-6, 1));
            $mid = substr($trna->seq(), 5, length($trna->seq())-11);
		}
		if (($pos5 eq "A" and $pos68 eq "T") or ($pos5 eq "T" and $pos68 eq "A") or
            ($pos5 eq "G" and $pos68 eq "C") or ($pos5 eq "C" and $pos68 eq "G") or
            ($pos5 eq "G" and $pos68 eq "T") or ($pos5 eq "T" and $pos68 eq "G"))
        {
            $trna->upstream($trna->upstream().substr($trna->seq(), 0, 1));
            $trna->downstream(substr($trna->seq(), length($trna->seq()) - 1).$trna->downstream());
			$trna->seq(substr($trna->seq(), 1, 3).$pos5.$mid.$pos68.substr($trna->seq(), length($trna->seq())-5, 4));
			$trna->ss(substr($trna->ss(), 1, 3).">".substr($trna->ss(), 5, length($trna->ss())-11)."<".substr($trna->ss(), length($trna->ss())-5, 3).".");
            $trna->start($trna->start() + 1);
            $trna->end($trna->end() - 1);
            my @ar_ac_pos = $trna->ar_ac_pos();
            if (scalar(@ar_ac_pos) > 0)
            {
                for (my $i = 0; $i < scalar(@ar_ac_pos); $i++)
                {
                    $ar_ac_pos[$i]->{rel_start} -= 1;
                    $ar_ac_pos[$i]->{rel_end} -= 1;
                }
                $trna->ar_ac_pos(@ar_ac_pos);
            }
            $self->rescore_tRNA($global_vars, $trna, $trna);
		}
	}
}

# Get tRNA anticodon, isotype, and intron location
sub decode_mito_tRNA_properties
{
    my $self = shift;
    my ($global_vars, $trna, $prescan_tRNA, $cur_cm_file) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $opts = $global_vars->{options};
    my $gc = $global_vars->{gc};
    my $log = $global_vars->{log_file};
    
    my ($anticodon, $acodon_index, $trna_type, $intron, $istart, $iend,
          $antiloop_index, $antiloop_end, $trna_len, $scan_len);

    $anticodon = "ERR";

    ($anticodon, $antiloop_index, $antiloop_end, $acodon_index) = $self->find_mito_anticodon($opts, $trna, $gc);
    $trna->anticodon($anticodon);
    $trna->add_ac_pos($acodon_index, $acodon_index + 2);
    
    # check for problem parsing anticodon loop 
    if (($anticodon eq $gc->undef_anticodon()) || ($trna->seq() eq 'Error'))    
    {
        $trna->anticodon($gc->undef_anticodon());
        $intron = 0;        
        
        if ($opts->save_odd_struct())
        {     
            open(ODDTRNA, ">>".$opts->odd_struct_file()) ||
                die "FATAL: Can't open ".$opts->odd_struct_file()." to save seconary structures\n\n"; 
            print ODDTRNA $prescan_tRNA->tRNAscan_id()." (".$trna->start()."-".$trna->end()."):\n".$trna->seq()."\n".$trna->ss()."\n\n"; 
            close(ODDTRNA);
        }
    }
    
    ($intron, $istart, $iend) = $self->find_intron($trna->seq(), $antiloop_index, $antiloop_end);
    
    if ($intron)
    {
        my $intron_start = 0;
        my $intron_end = 0;
        if ($trna->strand() eq "+")
        {
            $intron_start = $istart + $trna->start() - 1;
            $intron_end = $iend + $trna->start() - 1
        }
        else
        {
            $intron_start = $trna->end() - $iend + 1;
            $intron_end = $trna->end() - $istart + 1;
        }
    }
    
    my $isotype = $gc->get_tRNA_type($self, $trna->anticodon(), $cur_cm_file, $trna->model(), 0);
    my $vert_mito_isotype = $gc->get_vert_mito_type($trna->anticodon());
    
    my $model_iso = $trna->model();
    my $model_ac = "";
    
    if (length($trna->model()) > 3)
    {
        $model_iso = substr($trna->model(), 0, 3);
        if ($trna->model() =~ /^(\S+)_/)
        {
            my $temp = $1;
            if (length($temp) > 3)
            {
                $model_ac = substr($temp, 3);
            }
        }
        else
        {
            $model_ac = substr($trna->model(), 3);
        }
    }
    if ($model_iso ne $isotype)
    {
        if ($trna->category() eq "")
        {
            $trna->category("mito_inconsistent_isotype");
            if ($trna->note() ne "")
            {
                $trna->note("(".$isotype."); ".$trna->note());
            }
            else
            {
                $trna->note("(".$isotype.")");
            }
        }
        $log->broadcast($prescan_tRNA->tRNAscan_id().".trna".$trna->id()." - isotype/anticondon conflict\t Model: ".$trna->model()." Detected: ".$isotype.$anticodon."\n".
                        $trna->seq()."\n".$trna->ss()."\n");
    }
    if ($model_ac ne "" and $model_ac ne $anticodon)
    {
        if ($trna->category() eq "")
        {
            $trna->category("mito_inconsistent_ac");
            if ($trna->note() ne "")
            {
                $trna->note("(".$model_ac."); ".$trna->note());
            }
            else
            {
                $trna->note("(".$model_ac.")");
            }
        }
    }
    if ($vert_mito_isotype eq "")
    {
        if ($trna->category() eq "")
        {
            $trna->category("mito_unexpected_ac");
        }
    }
    
    $trna->isotype($model_iso);
}

sub scan_noncanonical_introns
{
    my $self = shift;
    my ($global_vars, $seq_name) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $opts = $global_vars->{options};
    my $log = $global_vars->{log_file};
    my $sp_int_results = $global_vars->{sp_int_results};
    my $sp_tRNAs = $global_vars->{sp_tRNAs};
    my $tmp_trnaseq_file = $global_constants->get("tmp_trnaseq_file");
    
    my $tRNA = tRNAscanSE::tRNA->new;
    my $tRNA_copy = tRNAscanSE::tRNA->new;
    my $cm_intron = tRNAscanSE::tRNA->new;
    my $cm_intron2 = tRNAscanSE::tRNA->new;
    my $previous_intron_len = 0;
    my $rnd2 = 0;
    my $add_ci = 0;
        
    my $arrayCMscanResults = tRNAscanSE::ArrayCMscanResults->new;
    my $cms_merge_file_rnd1 = $global_constants->get("temp_dir")."/tscan$$"."_intron_merge.out";
    my $count = $self->run_cmsearch_intron($global_vars, $seq_name, $arrayCMscanResults, 0, $cms_merge_file_rnd1, 1);
    if ($cms_merge_file_rnd1 eq "")
    {
        $log->error("Fail to run infernal for $seq_name");
        return 0;
    }

    $sp_int_results->clear_index();
    $sp_int_results->open_file("write");

    $arrayCMscanResults->open_file($cms_merge_file_rnd1);
    $arrayCMscanResults->get_next_cmsearch_hit($cm_intron);
    
    for (my $i = 0; $i < $sp_tRNAs->get_count(); $i++)
    {
        $rnd2 = 0;
        $add_ci = 0;
        $previous_intron_len = 0;
        $tRNA = $sp_tRNAs->get($i);
        $tRNA_copy->copy($tRNA);
        my $padded_seq = $tRNA->upstream().$tRNA->seq().$tRNA->downstream();
        if ($cm_intron->seqname() ne $tRNA->seqname().".trna".&pad_num($tRNA->id(), 6))
        {
			if ($tRNA->get_intron_count() == 0)
            {
                $sp_int_results->write_tRNA($tRNA);
                next;
			}
			else
            {
                $padded_seq = $tRNA->upstream().$tRNA->mat_seq().$tRNA->downstream();
                $rnd2 = 1;
                $add_ci = 1;
            }
		}
        else
        {        
            while (($cm_intron->seqname() ne "") and ($cm_intron->seqname() eq $tRNA->seqname().".trna".&pad_num($tRNA->id(), 6)))
            {
                my ($ret_value, $duplicate, $clip_seq, $intron_len) = $self->check_intron_validity($global_vars, $tRNA, $cm_intron, $padded_seq, $previous_intron_len);
                if ($ret_value)
                {
					$padded_seq = $clip_seq;
                    $previous_intron_len += $intron_len;
                    $rnd2 = 1;
                    if ($duplicate)
                    {
                        $add_ci = 1;
                    }
				}
                $arrayCMscanResults->get_next_cmsearch_hit($cm_intron);
            }
        }
        
        if ($rnd2)
        {
            my $trna_file = tRNAscanSE::Sequence->new;
            $trna_file->open_file($global_constants->get("tmp_trnaseq_file"), "write");
            $trna_file->set_seq_info($tRNA->seqname().".trna".&pad_num($tRNA->id(), 6), $tRNA->tRNAscan_id(), length($padded_seq), $padded_seq);
            $trna_file->write_fasta();
            $trna_file->close_file();
            
            my $rnd2IntronScanResults = tRNAscanSE::ArrayCMscanResults->new;
            my $cms_merge_file_rnd2 = $global_constants->get("temp_dir")."/tscan$$"."_intron_rnd2_merge.out";
            my $count = $self->run_cmsearch_intron($global_vars, $tRNA->tRNAscan_id(), $rnd2IntronScanResults, 0, $cms_merge_file_rnd2, 2);
            if ($cms_merge_file_rnd2 eq "")
            {
                $log->error("Fail to run infernal for ".$tRNA->tRNAscan_id());
                return 0;
            }
            
            $previous_intron_len = 0;
            $rnd2IntronScanResults->open_file($cms_merge_file_rnd2);
            $rnd2IntronScanResults->get_next_cmsearch_hit($cm_intron2);
            while ($cm_intron2->seqname() ne "")
            {
                my ($ret_value, $duplicate, $clip_seq, $intron_len) = $self->check_intron_validity($global_vars, $tRNA, $cm_intron2, $padded_seq, $previous_intron_len);
                if ($ret_value)
                {
                    $padded_seq = $clip_seq;
                    $previous_intron_len += $intron_len;
                    if ($duplicate)
                    {
                        $add_ci = 1;
                    }
                }
                
                $rnd2IntronScanResults->get_next_cmsearch_hit($cm_intron2);
            }
            $rnd2IntronScanResults->close_file();
        }
        
        my @ar_introns = $tRNA->ar_introns();
        my $nci_count = 0;
        my $ci_index = -1;
        for (my $i = 0; $i < scalar(@ar_introns); $i++)
        {
            if ($ar_introns[$i]->{type} eq "NCI")
            {
                $nci_count++;
            }
            elsif ($ar_introns[$i]->{type} eq "CI")
            {
                $ci_index = $i;
            }
        }
        
        if ($nci_count > 0)
        {
            my $ci_seq = "";
			if ($ci_index > -1 and $tRNA->model() ne "SeC")
            {
                if ($add_ci)
                {
					$ci_seq = $ar_introns[$ci_index]->{seq};
				}				
				$tRNA->remove_intron($ci_index);
			}
			$tRNA->sort_introns();
            
            if ($add_ci)
            {
                $self->add_canonical_intron($global_vars, $tRNA, $ci_seq);
            }
            $self->decode_nci_tRNA_properties($global_vars, $tRNA);

            $tRNA->tRNAscan_id($tRNA->seqname().".tRNA".$tRNA->id()."-".$tRNA->isotype().$tRNA->anticodon());
            $tRNA->set_default_scores();
            $sp_int_results->write_tRNA($tRNA);
		}
		else
        {
            $sp_int_results->write_tRNA($tRNA_copy);
        }
    }
    
    $sp_int_results->close_file();
    $arrayCMscanResults->close_file();
}

sub add_canonical_intron
{
    my $self = shift;
    my ($global_vars, $tRNA, $ci_seq) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $log = $global_vars->{log_file};
    my $ret_value = 1;
    
    my $mat_seq = "";
    my $precursor_seq = $tRNA->seq();
    my @introns = $tRNA->ar_introns();

    my $index = index(uc($precursor_seq), uc($ci_seq));
    if ($index > -1)
    {
        $precursor_seq = substr($precursor_seq, 0, $index).lc($ci_seq).substr($precursor_seq, $index + length($ci_seq));
    }
    for (my $i = 0; $i < scalar(@introns); $i++)
    {
        if ($i == 0)
        {
            $mat_seq = substr($precursor_seq, 0, $introns[$i]->{rel_start} - 1);
        }
        else
        {
            $mat_seq .= substr($precursor_seq, $introns[$i-1]->{rel_end}, $introns[$i]->{rel_start} - $introns[$i-1]->{rel_end} - 1);
        }
    }
    $mat_seq .= substr($precursor_seq, $introns[scalar(@introns)-1]->{rel_end});

    if (uc($mat_seq) ne uc($tRNA->mat_seq()))
    {
        my $trna_file = tRNAscanSE::Sequence->new;
        $trna_file->open_file($global_constants->get("tmp_trnaseq_file"), "write");
        $trna_file->set_seq_info($tRNA->seqname().".trna".&pad_num($tRNA->id(), 6), $tRNA->tRNAscan_id(), length($mat_seq), $mat_seq);
        $trna_file->write_fasta();
        $trna_file->close_file();
        
        my $scan_flag = 0;
        my $cms_output_file = "";
        my $file_idx = -1;
        my $arrayCMscanResults = tRNAscanSE::ArrayCMscanResults->new;
        my $cms_merge_file = $global_constants->get("temp_dir")."/tscan$$"."_add_intron_merge.out";
        my $cm_tRNA = tRNAscanSE::tRNA->new;
        
        foreach my $cur_cm (sort keys %{$self->{main_cm_file_path}})
        {
            $cms_output_file = $global_constants->get("temp_dir")."/tscan$$"."_add_intron_".$cur_cm.".out";
            if ($self->exec_cmsearch($self->{main_cm_file_path}->{$cur_cm}, $scan_flag, $global_constants->get("tmp_trnaseq_file"), $tRNA->tRNAscan_id(), $cms_output_file, $log, $self->{cm_cutoff}))
            {
                $file_idx = $arrayCMscanResults->add_file($cms_output_file, $cur_cm);
                $arrayCMscanResults->merge_result_file($file_idx, 0);
            }
            else
            {
                $ret_value = 0;
            }
        }
        $arrayCMscanResults->write_merge_file($cms_merge_file, 1);

        if ($arrayCMscanResults->get_result_count() > 0)
        {	
            $arrayCMscanResults->open_file($cms_merge_file);
            $arrayCMscanResults->get_next_cmsearch_hit($cm_tRNA);
            
            if ($cm_tRNA->seqname() ne "")
            {
                if ((uc(substr($cm_tRNA->seq(), length($cm_tRNA->seq()) - 3)) ne "CCA") and
                    (substr($cm_tRNA->ss(), length($cm_tRNA->ss()) - 4) eq "...."))
                {
                    if ($cm_tRNA->strand() eq "+")
                    {
                        $cm_tRNA->end($cm_tRNA->end() - 3);
                    }
                    else
                    {
                        $cm_tRNA->start($cm_tRNA->start() + 3);
                    }
                    $cm_tRNA->seq(substr($cm_tRNA->seq(), 0, length($cm_tRNA->seq()) - 3));
                    $cm_tRNA->ss(substr($cm_tRNA->ss(), 0, length($cm_tRNA->ss()) - 3));
                }
                
                $tRNA->ss($cm_tRNA->ss());
                $tRNA->mat_seq($cm_tRNA->seq());
                $tRNA->mat_ss($tRNA->ss());
                $tRNA->update_domain_model("infernal", $cm_tRNA->score(), $cm_tRNA->score(), 0, 0);
                $tRNA->model($cm_tRNA->model());
            }
            $arrayCMscanResults->close_file();
       }
    }
    return $ret_value;
}

# Run Infernal cmsearch for noncanonical introns
sub run_cmsearch_intron
{
    my $self = shift;
    my ($global_vars, $seq_name, $arrayCMscanResults, $merge_range, $cms_merge_file, $round) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $log = $global_vars->{log_file};
    
    my $scan_flag = 5;
    my $cms_output_file = "";
    my $file_idx = -1;
    
    # run cmsearch    
    foreach my $cur_cm (sort keys %{$self->{intron_cm_file_path}})
    {
        $cms_output_file = $global_constants->get("temp_dir")."/tscan$$"."_intron_".$cur_cm."_rnd".$round.".out";
        if ($self->exec_cmsearch($self->{intron_cm_file_path}->{$cur_cm}, $scan_flag, $global_constants->get("tmp_trnaseq_file"), $seq_name, $cms_output_file, $log, $self->{BHB_cm_cutoff}))
        {
            $file_idx = $arrayCMscanResults->add_file($cms_output_file, $cur_cm);
            $arrayCMscanResults->merge_result_file($file_idx, $merge_range);
        }
        else
        {
            return ("", 0);
        }
    }
    $arrayCMscanResults->write_merge_file($cms_merge_file, 0);
    my $count = $arrayCMscanResults->get_result_count();

    return $count;
}

sub check_intron_validity
{
    my $self = shift;
    my ($global_vars, $tRNA, $cm_intron, $padded_seq, $previous_intron_len) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $log = $global_vars->{log_file};
    my $ret_value = 1;
    my $duplicate = 0;

    my ($pre_intron, $intron, $post_intron, $pre_intron_seq, $intron_seq, $post_intron_seq);
    if ($cm_intron->ss() =~ /^([\<\-\.]{11,})(\-\<[<.]+[_.]{4,}[>.]{9,}\-[.]*\-)([-.>]+)$/)
    {
        $pre_intron  = $1;
        $intron      = $2;
        $post_intron = $3;
        my $full_intron_seq = $cm_intron->seq();
        $full_intron_seq =~ s/U/T/g; 
        $full_intron_seq =~ s/u/t/g;
        $cm_intron->seq($full_intron_seq);
        $intron_seq = substr($cm_intron->seq(), length($pre_intron), length($intron));
        $pre_intron_seq = substr($cm_intron->seq(), 0, length($pre_intron));
        $post_intron_seq = substr($cm_intron->seq(), length($pre_intron) + length($intron));
        $intron_seq =~ s/-//g;
        $pre_intron_seq =~ s/-//g;
        $post_intron_seq =~ s/-//g;
        $log->debug("Found intron ".$intron_seq." for ".$tRNA->tRNAscan_id());
    }
    else
    {
        $log->warning("Fail to parse intron sequence for ".$cm_intron->seqname()." from ".$cm_intron->seq()."\t".$cm_intron->ss());
        return 0;
    }

    my $padded_full_seq = $tRNA->upstream().$tRNA->seq().$tRNA->downstream();
    my ($upstream_start, $upstream_end, $downstream_start, $downstream_end);
    if ($tRNA->strand() eq "+")
    {
        $upstream_start = $tRNA->start() - length($tRNA->upstream());
        $downstream_end = $tRNA->end() + length($tRNA->downstream());
    }
    else
    {
        $downstream_start = $tRNA->start() - length($tRNA->downstream());
        $upstream_end = $tRNA->end() + length($tRNA->upstream());
    }
    
    my $clip_seq = substr($padded_seq, 0, $cm_intron->start() - $previous_intron_len + length($pre_intron_seq) - 1) . substr($padded_seq, $cm_intron->end() - $previous_intron_len - length($post_intron_seq));

    my $trna_file = tRNAscanSE::Sequence->new;
    $trna_file->open_file($global_constants->get("tmp_trnaseq_file"), "write");
    $trna_file->set_seq_info($tRNA->seqname().".trna".&pad_num($tRNA->id(), 6), $tRNA->tRNAscan_id(), length($clip_seq), $clip_seq);
    $trna_file->write_fasta();
    $trna_file->close_file();

    my $scan_flag = 0;
    my $cms_output_file = "";
    my $file_idx = -1;
    my $arrayCMscanResults = tRNAscanSE::ArrayCMscanResults->new;
    my $cms_merge_file = $global_constants->get("temp_dir")."/tscan$$"."_clip_intron_merge.out";
    my $cm_tRNA = tRNAscanSE::tRNA->new;
    
    foreach my $cur_cm (sort keys %{$self->{main_cm_file_path}})
    {
        $cms_output_file = $global_constants->get("temp_dir")."/tscan$$"."_clip_intron_".$cur_cm.".out";
        if ($self->exec_cmsearch($self->{main_cm_file_path}->{$cur_cm}, $scan_flag, $global_constants->get("tmp_trnaseq_file"), $tRNA->tRNAscan_id(), $cms_output_file, $log, $self->{cm_cutoff}))
        {
            $file_idx = $arrayCMscanResults->add_file($cms_output_file, $cur_cm);
            $arrayCMscanResults->merge_result_file($file_idx, 0);
        }
        else
        {
            $ret_value = 0;
        }
    }
    $arrayCMscanResults->write_merge_file($cms_merge_file, 1);
    
    if ($arrayCMscanResults->get_result_count() > 0 and $ret_value)
    {	
        $arrayCMscanResults->open_file($cms_merge_file);
        $arrayCMscanResults->get_next_cmsearch_hit($cm_tRNA);
        
        while ($cm_tRNA->seqname() ne "")
        {
            $duplicate = 0;
            $ret_value = 1;

            if ((uc(substr($cm_tRNA->seq(), length($cm_tRNA->seq()) - 3)) ne "CCA") and
                (substr($cm_tRNA->ss(), length($cm_tRNA->ss()) - 4) eq "...."))
            {
                if ($cm_tRNA->strand() eq "+")
                {
                    $cm_tRNA->end($cm_tRNA->end() - 3);
                }
                else
                {
                    $cm_tRNA->start($cm_tRNA->start() + 3);
                }
                $cm_tRNA->seq(substr($cm_tRNA->seq(), 0, length($cm_tRNA->seq()) - 3));
                $cm_tRNA->ss(substr($cm_tRNA->ss(), 0, length($cm_tRNA->ss()) - 3));
            }

            my ($intron_start, $intron_end);
            my $upstream_len = $cm_tRNA->start() - 1;
            my $downstream_len = length($clip_seq) - $cm_tRNA->end();
            my $seq = substr($padded_full_seq, $cm_tRNA->start() - 1, length($padded_full_seq) - $upstream_len - $downstream_len);
            my $pos = index(uc($seq), uc($intron_seq));
            if ($pos == -1)
            {
                $ret_value = 0;
                $log->debug("Intron out of boundary for ".$tRNA->tRNAscan_id()." ".$intron_seq);
            }
            else
            {
                $intron_start = $pos + 1;
                $intron_end = length($intron_seq) + $pos;
                my @ar_introns_src = $tRNA->ar_introns();
                my @ar_introns = ();
                my $intron_rec = {};
                # Adjust intron relative start and end when round1 seq is different from round2 seq
                for (my $i = 0; $i < scalar(@ar_introns_src); $i++)
                {
                    $intron_rec = {};
                    $intron_rec->{seq} = $ar_introns_src[$i]->{seq};
                    $intron_rec->{start} = $ar_introns_src[$i]->{start};
                    $intron_rec->{end} = $ar_introns_src[$i]->{end};
                    $intron_rec->{type} = $ar_introns_src[$i]->{type};
                    if (uc(substr($seq, $ar_introns_src[$i]->{rel_start} - 1, $ar_introns_src[$i]->{rel_end} - $ar_introns_src[$i]->{rel_start} + 1)) ne uc($ar_introns_src[$i]->{seq}))
                    {
                        my $pos_i = index(uc($seq), uc($ar_introns_src[$i]->{seq}));
                        $intron_rec->{rel_start} = $pos_i + 1;
                        $intron_rec->{rel_end} = length($ar_introns_src[$i]->{seq}) + $pos_i;
                    }
                    else
                    {
                        $intron_rec->{rel_start} = $ar_introns_src[$i]->{rel_start};
                        $intron_rec->{rel_end} = $ar_introns_src[$i]->{rel_end};
                    }
                    $ar_introns[$i] = $intron_rec;
                }
                # Check intron duplication
                for (my $i = 0; $i < scalar(@ar_introns); $i++)
                {
                    if ($ar_introns[$i]->{rel_start} == $intron_start and $ar_introns[$i]->{rel_end} == $intron_end)
                    {
                        $duplicate = 1;
                        $log->debug("Duplicate intron detected for ".$tRNA->tRNAscan_id()." ".$intron_seq);
                        last;
					}
                    elsif ($ar_introns[$i]->{type} eq "CI" and
                           length($ar_introns[$i]->{seq}) > 40 and
                           $ar_introns[$i]->{rel_start} < $intron_start and $ar_introns[$i]->{rel_end} > $intron_start and
                           $ar_introns[$i]->{rel_start} < $intron_end and $ar_introns[$i]->{rel_end} > $intron_end)
                    {
                        $ret_value = 0;
                        $log->debug("Inclusion of intron ".$tRNA->tRNAscan_id()." ".$intron_seq);
                        last;
                    }
                    elsif ($ar_introns[$i]->{type} eq "CI" and
#                           length($ar_introns[$i]->{seq}) > 40 and
                            $ar_introns[$i]->{rel_start} == $intron_start and
                           &seg_overlap($ar_introns[$i]->{rel_start}, $ar_introns[$i]->{rel_end}, $intron_start, $intron_end, 0))
                    {
                        $ret_value = 0;
                        $log->debug("Overlap with intron ".$tRNA->tRNAscan_id()." ".$intron_seq);
                        last;
                    }
                    elsif ($ar_introns[$i]->{type} eq "NCI" and
                           &seg_overlap($ar_introns[$i]->{rel_start}, $ar_introns[$i]->{rel_end}, $intron_start, $intron_end, 0))
                    {
                        $ret_value = 0;
                        $log->debug("Overlap with noncanonical intron ".$tRNA->tRNAscan_id()." ".$intron_seq);
                        last;
                    }
                }
            }
            
            my ($start, $end);
            if ($tRNA->strand() eq "+")
            {
                $upstream_end = $upstream_start + $upstream_len - 1;
                $downstream_start = $downstream_end - $downstream_len + 1;
                $start = $upstream_end + 1;
                $end = $downstream_start - 1;
            }
            else
            {
                $downstream_end = $downstream_start + $downstream_len - 1;
                $upstream_start = $upstream_end - $upstream_len + 1;
                $start = $downstream_end + 1;
                $end = $upstream_start - 1;
            }
            
            my $hit_overlap = &seg_overlap($tRNA->start(), $tRNA->end(), $start, $end, 40);
            if ($ret_value and $hit_overlap and $cm_tRNA->score() > $tRNA->score() and length($cm_tRNA->seq()) >= $self->{min_tRNA_no_intron})
            {
                # Adjust intron relative start and end when round1 seq is different from round2 seq
                my @ar_introns = $tRNA->ar_introns();
                for (my $i = 0; $i < scalar(@ar_introns); $i++)
                {
                    if (uc(substr($seq, $ar_introns[$i]->{rel_start} - 1, $ar_introns[$i]->{rel_end} - $ar_introns[$i]->{rel_start} + 1)) ne uc($ar_introns[$i]->{seq}))
                    {
                        my $pos_i = index(uc($seq), uc($ar_introns[$i]->{seq}));
                        $tRNA->set_intron($i, $pos_i + 1, length($ar_introns[$i]->{seq}) + $pos_i, $ar_introns[$i]->{type}, $ar_introns[$i]->{seq});
                    }
                }
                $tRNA->upstream(substr($clip_seq, 0, $upstream_len));
                $tRNA->downstream(substr($clip_seq, $cm_tRNA->end()));
                $tRNA->seq($seq);
                $tRNA->ss($cm_tRNA->ss());
                $tRNA->mat_seq($cm_tRNA->seq());
                $tRNA->mat_ss($tRNA->ss());
                $tRNA->start($start);
                $tRNA->end($end);
                $tRNA->update_domain_model("infernal", $cm_tRNA->score(), $cm_tRNA->score(), 0, 0);
                $tRNA->set_default_scores();
                $tRNA->model($cm_tRNA->model());
                if (!$duplicate)
                {
                    $tRNA->add_rel_intron($intron_start, $intron_end, "NCI", $intron_seq);
                }
				last;
			}
            else
            {
                $ret_value = 0;
            }
            $arrayCMscanResults->get_next_cmsearch_hit($cm_tRNA);
        }
        $arrayCMscanResults->close_file();
    }

    return ($ret_value, $duplicate, $clip_seq, length($intron_seq));
}

sub parse_covels_hit
{
    my $self = shift;
    my ($opts, $covels_hit, $r_covels_hit_elements, $prescan_tRNA) = @_;

    my $covels_hit_found = 0;

    if ($covels_hit =~ /^\s*(\S+)\s+(\d+)\s+(\d+).+: (\S+)\s*/o)
    {
        $r_covels_hit_elements->{score} = $1;
        $r_covels_hit_elements->{subseq_start} = $2;
        $r_covels_hit_elements->{subseq_end} = $3;
        $r_covels_hit_elements->{hit_seqname} = $4;
        $covels_hit_found = 1;        
    }

    if ($covels_hit_found)
    {
        if ($prescan_tRNA->strand() eq "+")
        {
            $r_covels_hit_elements->{tRNA_start} = $prescan_tRNA->start() + $r_covels_hit_elements->{subseq_start} - 1;        
            $r_covels_hit_elements->{tRNA_end} = $prescan_tRNA->start() + $r_covels_hit_elements->{subseq_end} - 1;
        }
        else
        {
            $r_covels_hit_elements->{tRNA_start} = $prescan_tRNA->start() - $r_covels_hit_elements->{subseq_start} + 1;        
            $r_covels_hit_elements->{tRNA_end} = $prescan_tRNA->start() - $r_covels_hit_elements->{subseq_end} + 1;
        }                
        
        return 1;                
    }
    else
    {
        return 0;
    }
}                                

# Run covels, return hits in $covels_hit_list array
sub run_covels
{
    my $self = shift;
    my ($global_vars, $r_covels_hit_list, $r_cur_cm_file, $prescan_tRNA) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $opts = $global_vars->{options};
    my $stats = $global_vars->{stats};
    my $log = $global_vars->{log_file};

    my $scan_len = 0;
    my %covels_hit_elements = ();

    $self->set_search_params($opts, \$scan_len, $r_cur_cm_file, $self->{max_cove_tRNA_length}, $prescan_tRNA);

    # set covels reporting threshold below 0 (default) if -X param is
    # set below 0 by user
    my $report_cutoff = &min(0, $self->{cm_cutoff});
    
    my $complement = "-c";
    if ($opts->tscan_mode() || $opts->eufind_mode() || $opts->infernal_fp())
    {
        $complement = "";
    }
    
    # run Covels
    my $covels_cmd = $self->covels_bin()." $complement -w$scan_len -t$report_cutoff $$r_cur_cm_file ".$global_constants->get("tmp_trnaseq_file");
    $log->broadcast($covels_cmd);
    my $covels_output = `$covels_cmd`;

    if ($? != 0)
    {
        $log->error("Covels-SE cannot be completed successfully for ".$prescan_tRNA->seqname().". (Exit Code: $?)");
        print "Exit first loop at 1\n";
        return 0;
    }
    
    my ($junk, $allhits) = split(/----------\n\n/, $covels_output);
    @$r_covels_hit_list = split(/\n/, $allhits);

    # count no. of hits over cutoff
    my $total_hits = 0;   
    foreach my $covels_hit (@$r_covels_hit_list)
    {
        %covels_hit_elements = ();
        if (($self->parse_covels_hit($opts, $covels_hit, \%covels_hit_elements, $prescan_tRNA)) &&
            ($covels_hit_elements{score} >= $self->{cm_cutoff}))
        {
            $total_hits++;
        }        
    }
    
    # if no tRNAs detected when using a selenocysteine cove model,
    #  try main model and run again before giving up
    if (($total_hits == 0) && 
        (($$r_cur_cm_file eq $self->{Pselc_cm_file_path}) || ($$r_cur_cm_file eq $self->{Eselc_cm_file_path})))
    {
        $$r_cur_cm_file = $self->{main_cm_file_path}->{Domain};
        
        # re-run Covels with main model    
        $covels_cmd = $self->covels_bin()." -w$scan_len -t$report_cutoff $$r_cur_cm_file ".$global_constants->get("tmp_trnaseq_file");
        $covels_output = `$covels_cmd`;
        if ($? != 0)
        {
            $log->error("Covels-SE cannot be completed successfully for ".$prescan_tRNA->seqname().". (Exit Code: $?)");
            print "Exit first loop at 2\n";
            return 0;
        }
        
        ($junk,$allhits) = split(/----------\n\n/, $covels_output);
        @$r_covels_hit_list = split(/\n/, $allhits);
    }

    # Go thru hit list, save info for tRNA hits with sub-cutoff scores
    my $ct = 0;
    my $over_cutoff = 0;
    my $trna_desc = "";

    foreach my $covels_hit (@$r_covels_hit_list)
    {
        %covels_hit_elements = ();
        if ($self->parse_covels_hit($opts, $covels_hit, \%covels_hit_elements, $prescan_tRNA))
        {
            $ct++;
            if ($covels_hit_elements{score} >= $self->{cm_cutoff})
            {
                $over_cutoff++;
            }
            else
            {
                $log->broadcast("Low covels score for ".$prescan_tRNA->tRNAscan_id().".$ct: $covels_hit_elements{score}");
                $trna_desc .= "(Cove Hit#$ct: $covels_hit_elements{tRNA_start}-$covels_hit_elements{tRNA_end},".
                    " Sc: $covels_hit_elements{score},  Len: ".(abs($covels_hit_elements{tRNA_start} - $covels_hit_elements{tRNA_end}) + 1).") ";
            }
        }
    }        
    
    # report if no scores over 0 bit reporting threshold
    if ($over_cutoff == 0)
    {
        if ((!$opts->results_to_stdout()) &&
            ($opts->eufind_mode() || $opts->tscan_mode() || $opts->use_prev_ts_run()))
        {
            $log->broadcast("Covels score(s) below cutoff for ".$prescan_tRNA->tRNAscan_id().". Skipping...");
        }
        if ($opts->save_falsepos())
        {
            my $fulltrnaDesc = "(Fp Hit: ".$prescan_tRNA->start()."-".$prescan_tRNA->end().", ".
                (abs($prescan_tRNA->start() - $prescan_tRNA->end()) + 1)." bp, Src: ".$prescan_tRNA->hit_source().") ".$trna_desc;
        
            $stats->increment_fpos_base_ct(length($prescan_tRNA->seq()));          
            &write_tRNA($opts->falsepos_file(), $prescan_tRNA->tRNAscan_id(), $fulltrnaDesc, $prescan_tRNA->seq(), 0);
        }           
    }

    return 1;
}

sub run_coves
{
    my $self = shift;
    my ($log, $tmp_trnaseq_file, $seq_name, $cm_file) = @_;

    my ($covseq, $covss, $coves_output, $junk, @coves_lines, $sec_struct, $coves_score);
    
    my $coves_cmd = $self->{coves_bin}." -s $cm_file $tmp_trnaseq_file";
    $log->broadcast($coves_cmd);
    $coves_output = `$coves_cmd`;
    if ($? != 0)
    {
        $log->error("Coves-SE cannot be completed successfully for $seq_name. (Exit Code: $?)");
        return ("Error", "", -1);
    }

    ($junk, $sec_struct) = split(/----------\n\n/, $coves_output);
    @coves_lines = split(/\n/,$sec_struct);
    $covseq = '';
    $covss = '';
    $coves_score = -1000;
    $seq_name =~ s/(\W)/\\$1/g;

    foreach (@coves_lines)
    {
        if (/^\s+$seq_name\s([a-zA-Z\-]{1,60})\s*/)
        {
            $covseq .= $1;
        } 
        if (/^\s+$seq_name\s([\.\<\>\ ]{1,60})/)
        {
            $covss .= $1;
        }
        if (/^\s*(\S+)\sbits\s:\s$seq_name/)
        {
            $coves_score = $1; 
        }
    }
    $covss =~ s/\s//g;     #  take spaces out of alignment        
    $covseq =~ s/-//g;     #  take '-' gaps out of seq

    if (($covseq eq '') || ($covss eq ''))
    {
        print STDERR "Could not complete coves successfully for $seq_name ",
        "because unable to parse coves secondary structure string. ",
        "Skipping tRNA anticodon & type prediction\n";
        return ("Error", "", -1);
    }

    return ($covseq, $covss, $coves_score);
}

sub analyze_with_cove
{
    my $self = shift;
    my ($global_vars, $prescan_tRNA, $r_curseq_trnact) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $opts = $global_vars->{options};
    my $stats = $global_vars->{stats};
    my $log = $global_vars->{log_file};
   
    my (@covels_hit_list, $cur_cm_file, $cove_confirmed_ct);
    my ($covseq, $covss, $coves_score);
    my %covels_hit_elements = ();
    my $covels_tRNA = undef;
    my $cov_hit = {};    
    my $rescore = 0;
    
    $cove_confirmed_ct = 0;
    if (!$self->run_covels($global_vars, \@covels_hit_list, \$cur_cm_file, $prescan_tRNA))
    {
        return 0;
    }
        
    # Loop to parse covels tRNA hit(s) and run Coves on each tRNA    
    foreach my $covels_hit (@covels_hit_list)
    {
        $rescore = 0;
        %covels_hit_elements = ();
        if ((!$self->parse_covels_hit($opts, $covels_hit, \%covels_hit_elements, $prescan_tRNA)) ||
            ($covels_hit_elements{score} < $self->{cm_cutoff}))
        {
            next; 
        }                       
        $covels_tRNA = tRNAscanSE::tRNA->new;
        $$r_curseq_trnact++;

        $covels_tRNA->id($$r_curseq_trnact);
        $covels_tRNA->set_covels_hit(\%covels_hit_elements);
        $covels_tRNA->upstream($prescan_tRNA->upstream());
        $covels_tRNA->downstream($prescan_tRNA->downstream());
        if (($covels_hit_elements{subseq_start} == 1) && ($covels_hit_elements{subseq_end} == $prescan_tRNA->len()))
        {
            $covels_hit_elements{tRNA_len} = $prescan_tRNA->len();
        }
        else
        {
            # get correct subseq for coves & save to file            
            if ($covels_hit_elements{subseq_start} < $covels_hit_elements{subseq_end})
            {
                $covels_hit_elements{tRNA_len} = $covels_hit_elements{subseq_end} - $covels_hit_elements{subseq_start} + 1;
                &write_tRNA($global_constants->get("tmp_trnaseq_file"), $covels_tRNA->seqname(), " ",
                        substr($prescan_tRNA->seq(), $covels_hit_elements{subseq_start} - 1, $covels_hit_elements{tRNA_len}), 1);
                if ($covels_hit_elements{subseq_start} > 1)
                {
                    $covels_tRNA->upstream($covels_tRNA->upstream() . substr($prescan_tRNA->seq(), 0, $covels_hit_elements{subseq_start} - 1));
                }
                if ($covels_hit_elements{subseq_end} < $prescan_tRNA->len())
                {
                    $covels_tRNA->downstream(substr($prescan_tRNA->seq(), $covels_hit_elements{subseq_end}) . $covels_tRNA->downstream());
                }
            }
            else
            {
                $covels_hit_elements{tRNA_len} = $covels_hit_elements{subseq_start} - $covels_hit_elements{subseq_end} + 1;
                &write_tRNA($global_constants->get("tmp_trnaseq_file"), $covels_tRNA->seqname(), " ",
                        &rev_comp_seq(substr($prescan_tRNA->seq(), $covels_hit_elements{subseq_end} - 1, $covels_hit_elements{tRNA_len})), 1);
                if ($covels_hit_elements{subseq_end} > 1)
                {
                    $covels_tRNA->upstream($covels_tRNA->upstream() . &rev_comp_seq(substr($prescan_tRNA->seq(), $covels_hit_elements{subseq_start} - 1, $global_constants->get("upstream_len"))));
                }
                if ($covels_hit_elements{subseq_start} < $prescan_tRNA->len())
                {
                    $covels_tRNA->downstream(&rev_comp_seq(substr($prescan_tRNA->seq(), $covels_hit_elements{subseq_end} - $global_constants->get("downstream_len"), $global_constants->get("downstream_len")) . $covels_tRNA->downstream()));
                }
            }
        }                       
        $stats->increment_coves_base_ct($covels_hit_elements{tRNA_len});
    
        ($covseq, $covss, $coves_score) = $self->run_coves($log, $global_constants->get("tmp_trnaseq_file"), $prescan_tRNA->seqname(), $cur_cm_file);
        $covels_tRNA->seq($covseq);
        $covels_tRNA->ss($covss);
        $covels_tRNA->update_domain_model("cove", $coves_score, $coves_score, 0, 0);        
        
        if ((uc(substr($covels_tRNA->seq(), length($covels_tRNA->seq()) - 3)) ne "CCA") &&
            (substr($covels_tRNA->ss(), length($covels_tRNA->ss()) - 4) eq "....") &&
            ($self->{main_cm_file_path}->{Domain} ne $global_constants->get_subvalue("cove_cm", "eukaryota")) &&
            ($self->{main_cm_file_path}->{Domain} ne $global_constants->get_subvalue("cove_cm", "general")))
        {
            $covels_hit_elements{subseq_end} -= 3;
            if ($covels_tRNA->strand() eq "+")
            {
                $covels_tRNA->end($covels_tRNA->end() - 3);
            }
            else
            {
                $covels_tRNA->start($covels_tRNA->start() + 3);
            }
            $covels_tRNA->seq(substr($covels_tRNA->seq(), 0, length($covels_tRNA->seq()) - 3));
            $covels_tRNA->ss(substr($covels_tRNA->ss(), 0, length($covels_tRNA->ss()) - 3));
            $rescore = 1;
        }
        elsif ((uc(substr($covels_tRNA->seq(), length($covels_tRNA->seq()) - 3)) eq "CCA") &&
            (((substr($covels_tRNA->ss(), length($covels_tRNA->ss()) - 6) eq "<.....") &&
            (substr($covels_tRNA->seq(), length($covels_tRNA->seq()) - 5, 2) =~ /[acgtn][ACGTN]/)) ||
            ((substr($covels_tRNA->ss(), length($covels_tRNA->ss()) - 7) eq "<......") &&
            (substr($covels_tRNA->seq(), length($covels_tRNA->seq()) - 6, 3) =~ /[ACGTN][acgtn][ACGTN]/)) ||
            ((substr($covels_tRNA->ss(), length($covels_tRNA->ss()) - 7) eq "<......") &&
            (substr($covels_tRNA->seq(), length($covels_tRNA->seq()) - 6, 3) =~ /[acgtn]{2}[ACGTN]/))) &&
            ($self->{main_cm_file_path}->{Domain} ne $global_constants->get_subvalue("cove_cm", "eukaryota")) &&
            ($self->{main_cm_file_path}->{Domain} ne $global_constants->get_subvalue("cove_cm", "general")))
        {
            my $trim_len = 4;
            if (substr($covels_tRNA->ss(), length($covels_tRNA->ss()) - 7) eq "<......" && substr($covels_tRNA->seq(), length($covels_tRNA->seq()) - 6, 3) =~ /[acgtn]{2}[ACGTN]/)
            {
                $trim_len = 5;
            }
            
            if ($covels_tRNA->strand() eq "+")
            {
                $covels_tRNA->end($covels_tRNA->end() - $trim_len);
            }
            else
            {
                $covels_tRNA->start($covels_tRNA->start() + $trim_len);
            }
            $covels_tRNA->seq(substr($covels_tRNA->seq(), 0, length($covels_tRNA->seq()) - $trim_len));
            $covels_tRNA->seq(substr($covels_tRNA->seq(), 0, length($covels_tRNA->seq()) - 1).uc(substr($covels_tRNA->seq(), length($covels_tRNA->seq()) - 1, 1)));
            $covels_tRNA->ss(substr($covels_tRNA->ss(), 0, length($covels_tRNA->ss()) - $trim_len));
            $rescore = 1;
        }

        if ($rescore)
        {
            my $trna_file = tRNAscanSE::Sequence->new;
            $trna_file->open_file($global_constants->get("tmp_trnaseq_file"), "write");
            $trna_file->set_seq_info($prescan_tRNA->seqname(), $prescan_tRNA->seqname(), length($covels_tRNA->seq()), $covels_tRNA->seq());
            $trna_file->write_fasta();
            $trna_file->close_file();
            ($covseq, $covss, $coves_score) = $self->run_coves($log, $global_constants->get("tmp_trnaseq_file"), $prescan_tRNA->seqname(), $cur_cm_file);
            $covels_tRNA->update_domain_model("cove", $coves_score, $coves_score, 0, 0);  
		}
		
        # look for intron
        $self->decode_tRNA_properties($global_vars, $covels_tRNA, $prescan_tRNA, $cur_cm_file);
        $covels_tRNA->set_mature_tRNA();
        
        $covels_tRNA->tRNAscan_id($covels_tRNA->seqname().".tRNA".$$r_curseq_trnact."-".$covels_tRNA->isotype().$covels_tRNA->anticodon());
        $covels_tRNA->src_seqlen($prescan_tRNA->src_seqlen());
        $covels_tRNA->src_seqid($covels_hit_elements{hit_seqname});
        $covels_tRNA->hit_source($prescan_tRNA->hit_source());
        $covels_tRNA->model($cur_cm_file);
        $covels_tRNA->ordered_seqname($prescan_tRNA->src_seqid());
        $covels_tRNA->set_default_scores();
        
        if (!$self->{CM_check_for_introns})
        {
            $cove_confirmed_ct++;
        }
        $global_vars->{sp_int_results}->write_tRNA($covels_tRNA);
    }            # while more covels_hits
   
    return $cove_confirmed_ct;
}

# Format command and run Infernal cmscan
sub exec_cmscan
{
    my $self = shift;
    my ($cm_db, $scan_flag, $tmp_trnaseq_file, $seq_name, $cms_output_file, $cms_tab_file, $log) = @_;
    
    my $cm_options = "-g --nohmm --toponly --notrunc";
    if ($scan_flag == 1)
    {
        $cm_options = "-g --mid --notrunc";
    }
    elsif ($scan_flag == 2)
    {
        $cm_options = "-g --mid --toponly --notrunc";
    }
    elsif ($scan_flag == 3)
    {
        $cm_options = "-g --mid --notrunc -T ".$self->{infernal_fp_cutoff};
    }
    elsif ($scan_flag == 4)
    {
        $cm_options = "-g --nohmm --notrunc";
    }
    if ($self->{infernal_thread} != -1)
    {
		$cm_options .= " --cpu ".$self->{infernal_thread};
	}
    
    my $cms_cmd = "$self->{cmscan_bin} $cm_options --fmt 2 --tblout $cms_tab_file -o $cms_output_file $cm_db $tmp_trnaseq_file";
    $log->broadcast($cms_cmd);
    system($cms_cmd);

    if ($? != 0)
    {
        $log->error("Infernal cmscan cannot be completed successfully for ".$seq_name.". $cms_cmd (Exit Code: $? $!)");
        return 0;
    }
    return 1;
}

# Format command and run Infernal cmsearch
sub exec_cmsearch
{
    my $self = shift;
    my ($cm_file, $scan_flag, $tmp_trnaseq_file, $seq_name, $cms_output_file, $log, $score_cutoff) = @_;
    
    my $cm_options = "-g --nohmm --toponly --notrunc";
    if ($scan_flag == 1)
    {
        $cm_options = "-g --mid --notrunc";
    }
    elsif ($scan_flag == 2)
    {
        $cm_options = "-g --mid --toponly --notrunc";
    }
    elsif ($scan_flag == 3)
    {
        $cm_options = "-g --mid --notrunc -T ".$self->{infernal_fp_cutoff};
    }
    elsif ($scan_flag == 4)
    {
        $cm_options = "-g --nohmm --notrunc";
    }
    elsif ($scan_flag == 5)
    {
        $cm_options = "-g --max --toponly --notrunc --notextw -T ".$self->{BHB_cm_cutoff};
    }
    elsif ($scan_flag == 6)
    {
        $cm_options = "-g --toponly --notextw";
    }
    elsif ($scan_flag == 7)
    {
        $cm_options = "-g --max --toponly --notrunc --notextw -T 0";
    }
    
    if ($scan_flag != 3 and $scan_flag != 5 and $scan_flag != 7)
    {
		if ($score_cutoff <= 10)
        {
			$cm_options .= " -T ".$score_cutoff;
		}
	}
    if ($self->{infernal_thread} != -1)
    {
		$cm_options .= " --cpu ".$self->{infernal_thread};
	}
    
    my $cms_cmd = "$self->{cmsearch_bin} $cm_options $cm_file $tmp_trnaseq_file > $cms_output_file";
    $log->broadcast($cms_cmd);
    system($cms_cmd);

    if ($? != 0)
    {
        $log->error("Infernal cmsearch cannot be completed successfully for ".$seq_name.". $cms_cmd (Exit Code: $? $!)");
        return 0;
    }
    
    return 1;
}

# Run Infernal cmsearch, return results in $r_cms_hit_list array reference
sub run_cmsearch
{
    my $self = shift;
    my ($global_vars, $seq_name, $arrayCMscanResults, $merge_range) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $opts = $global_vars->{options};
    my $stats = $global_vars->{stats};
    my $log = $global_vars->{log_file};
    
    my $score_cutoff = -1;
    my $scan_flag = -1;
    my $cms_output_file = "";
    my $file_idx = -1;
    my $cms_merge_file = $global_constants->get("temp_dir")."/tscan$$"."_cm_merge.out";
    
    # run cmsearch
    if ($opts->tscan_mode() || $opts->eufind_mode() || $opts->infernal_fp())
    {
        $scan_flag = 0;
        if ($opts->hmm_filter())
        {
            $scan_flag = 2;
        }
    }
    else
    {
        $scan_flag = 4;
        if ($opts->hmm_filter())
        {
            $scan_flag = 1;
        }
    }
    
#    if ($opts->mito_mode())
#    {
#		$score_cutoff = $self->{organelle_cm_cutoff};
#	}
#	else
#    {
        $score_cutoff = $self->{cm_cutoff};
#    }
    
    foreach my $cur_cm (sort keys %{$self->{main_cm_file_path}})
    {
        $cms_output_file = $global_constants->get("temp_dir")."/tscan$$"."_cm_".$cur_cm.".out";
        if ($self->exec_cmsearch($self->{main_cm_file_path}->{$cur_cm}, $scan_flag, $global_constants->get("tmp_trnaseq_file"), $seq_name, $cms_output_file, $log, $score_cutoff))
        {
            $file_idx = $arrayCMscanResults->add_file($cms_output_file, $cur_cm);
            $arrayCMscanResults->merge_result_file($file_idx, $merge_range);
        }
        else
        {
            return "";
        }
    }
    $arrayCMscanResults->write_merge_file($cms_merge_file, 1);
    my $count = $arrayCMscanResults->get_result_count();

    return ($cms_merge_file, $count);
}

# Run Infernal for scanning truncation in predicited tRNA 
sub truncated_tRNA_search
{
    my $self = shift;
    my ($global_vars, $seq_name) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $opts = $global_vars->{options};
    my $log = $global_vars->{log_file};
    my $arrayCMscanResults = tRNAscanSE::ArrayCMscanResults->new;
    my $sp_int_results = $global_vars->{sp_int_results};
    my $tmp_trnaseq_file = $global_constants->get("tmp_trnaseq_file");
    my $tRNA = tRNAscanSE::tRNA->new;
    my $cm_tRNA = tRNAscanSE::tRNA->new;

    my %cms_output_file = ();
    my $scan_flag = 6;
    my $file_idx = -1;
    my $cms_merge_file = $global_constants->get("temp_dir")."/tscan$$"."_trunc_cm_merge.out";

    my $sp_trunc_int_results = tRNAscanSE::IntResultFile->new;
    $sp_trunc_int_results->file_name($opts->truncated_int_result_file());

    foreach my $cur_cm (sort keys %{$self->{main_cm_file_path}})
    {
        my $count = &write_tRNAs($tmp_trnaseq_file, $sp_int_results, 1, 1, $cur_cm);

        if ($count > 0)
        {
            $cms_output_file{$cur_cm} = $global_constants->get("temp_dir")."/tscan$$"."_trunc_cm_".$cur_cm.".out";
            if ($self->exec_cmsearch($self->{main_cm_file_path}->{$cur_cm}, $scan_flag, $global_constants->get("tmp_trnaseq_file"), $seq_name, $cms_output_file{$cur_cm}, $log, $self->{cm_cutoff}))
            {
                $file_idx = $arrayCMscanResults->add_file($cms_output_file{$cur_cm}, $cur_cm);
                $arrayCMscanResults->merge_result_file($file_idx, 0);
            }
            else
            {
                return "";
            }
        }
    }        
    $arrayCMscanResults->write_merge_file($cms_merge_file, 0);
    
    my @sp_indexes = $sp_int_results->get_indexes();
    if ($sp_int_results->open_file("read") and $sp_trunc_int_results->open_file("write"))
    {
        $arrayCMscanResults->open_file($cms_merge_file);
        my $cm_model = $arrayCMscanResults->get_next_cmsearch_hit($cm_tRNA);

        for (my $i = 0; $i < scalar(@sp_indexes); $i++)
        {
            $sp_int_results->get_tRNA($sp_indexes[$i]->[0], $tRNA);
            if (($cm_tRNA->seqname() ne "") and ($cm_tRNA->seqname() eq $tRNA->seqname().".t".&pad_num($tRNA->id(), 6)))
            {
                my $trunc_label = $self->check_truncation($tRNA, $cm_tRNA, $cm_model, $global_constants);
                $tRNA->trunc($trunc_label);
                $cm_model = $arrayCMscanResults->get_next_cmsearch_hit($cm_tRNA);
            }
            $sp_trunc_int_results->write_tRNA($tRNA);
        }
        $sp_int_results->close_file();
        $sp_trunc_int_results->close_file();
        $sp_int_results->clear_index();
        $global_vars->{sp_int_results} = $sp_trunc_int_results;
    }
}

sub check_truncation
{
    my $self = shift;
    my ($tRNA, $trunc_tRNA, $cm_model, $global_constants) = @_;
    
    my $label = "";
    if ($trunc_tRNA->trunc() eq "" or $trunc_tRNA->trunc() eq "no")
    {}
    else
    {
        if (index($trunc_tRNA->trunc(),"5'") > -1)
        {
            if ($cm_model =~/^\<\[(\d+)\]\*/)
            {
                $label = "trunc_start:".$1;
            }		
        }
        if (index($trunc_tRNA->trunc(), "3'") > -1)
        {
            if ($cm_model =~/\*\[(\d+)\]\>$/)
            {
                my $diff = $1;
                if ($diff <= 3 and ((uc(substr($tRNA->mat_seq(), length($tRNA->mat_seq()) - 3)) ne "CCA" and index($trunc_tRNA->trunc(),"5'") == -1) or
                        index($trunc_tRNA->trunc(),"5'") > -1) and
                    ($self->{main_cm_file_path}->{Domain} ne $global_constants->get_subvalue("cm", "eukaryota")) and
                    ($self->{main_cm_file_path}->{Domain} ne $global_constants->get_subvalue("cm", "general")))
                {}
                else
                {
                    if ($label ne "")
                    {
						$label .= ",";
					}					
                    $label .= "trunc_end:".$diff;
                }
            }
        }
    }

    return $label;
}

# Run Infernal with isotype specific CMs for predicited tRNA
sub isotype_cmsearch
{
    my $self = shift;
    my ($global_vars) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $opts = $global_vars->{options};
    my $sp_int_results = $global_vars->{sp_int_results};
    my $iso_int_results = $global_vars->{iso_int_results};
    my $tmp_trnaseq_file = $global_constants->get("tmp_trnaseq_file");
    my $tRNA = tRNAscanSE::tRNA->new;
    my $seqname = "";

    &write_tRNAs($tmp_trnaseq_file, $sp_int_results, 1, 1, "");

    my @sp_indexes = $sp_int_results->get_indexes();
    if ($sp_int_results->open_file("read") && $iso_int_results->open_file("write"))
    {
        $iso_int_results->write_line("");
        for (my $i = 0; $i < scalar(@sp_indexes); $i++)
        {
            $sp_int_results->get_tRNA($sp_indexes[$i]->[0], $tRNA);
            $seqname = $tRNA->seqname();
            $iso_int_results->write_line($tRNA->seqname().".t".&pad_num($tRNA->id(), 6));
        }
        $iso_int_results->close_file();
        $sp_int_results->close_file();
    }
    
    # run cmscan
    $self->scan_isotype_cm($global_vars, $self->{isotype_cm_db_file_path}, "cyto", $seqname);
    if ($opts->euk_mode() and $opts->mito_model() ne "")
    {
        $self->scan_isotype_cm($global_vars, $self->{mito_isotype_cm_db_file_path}, "mito", $seqname);        
    }
}

sub scan_isotype_cm
{
    my $self = shift;
    my ($global_vars, $cm_db, $type, $seq_name) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $opts = $global_vars->{options};
    my $log = $global_vars->{log_file};
    my $iso_int_results = $global_vars->{iso_int_results};
    my $tmp_trnaseq_file = $global_constants->get("tmp_trnaseq_file");
    my $cmscanResults = undef;
    my $tRNA = tRNAscanSE::tRNA->new; 

    my $scan_flag = 2;

    my $iso_result_file = $global_constants->get("temp_dir")."/tscan$$"."_iso_temp.out";
    my $iso_results_out = tRNAscanSE::MultiResultFile->new;
    my ($cur_iso_cm, $cur_iso_cm_file);
    my $cms_output_file = "";
    my $cms_tab_file = "";
    my $line = "";
    my $tRNAid = "";
    my $iso_index_ct = -1;
    my %found_models = ();
    
    my $models = $self->get_models_from_db($cm_db);
    
    $iso_results_out->file_name($iso_result_file);
    
    $cms_output_file = $global_constants->get("temp_dir")."/tscan$$"."_iso_cm.out";
    $cms_tab_file = $global_constants->get("temp_dir")."/tscan$$"."_iso_cm.tab";
    if ($self->exec_cmscan($cm_db, $scan_flag, $global_constants->get("tmp_trnaseq_file"), $seq_name, $cms_output_file, $cms_tab_file, $log))
    {
        $cmscanResults = tRNAscanSE::CMscanResultFile->new($cms_tab_file, "multi");
    }
    else
    {
        return 0;
    }
    
    if ($iso_int_results->open_file("read") && $iso_results_out->open_file("write"))
    {
        # header line
        $line = $iso_int_results->get_line();
        foreach $cur_iso_cm (sort @$models)
        {
			my $model = $cur_iso_cm;
        	if ($cur_iso_cm =~ /^arch-(\S+)/ || $cur_iso_cm =~ /^euk-(\S+)/ || $cur_iso_cm =~ /^bact-(\S+)/)
            {
				$model = $1;
			}
			
            if ($type eq "mito")
            {
                $line = $line."\tmito_".$model;
            }
            else
            {
                $line = $line."\t".$model;
            }
        }
        $iso_results_out->write_line($line);

        #content
        if ($cmscanResults->open_file())
        {
            my $hits = $cmscanResults->get_next_tab_seq_hits();
            ($tRNAid, $line) = $iso_int_results->get_next_record();
            while ($tRNAid ne "")
            {
                %found_models = ();
                if (scalar(@$hits) > 0)
                {
                    $iso_index_ct = 0;
                    if ($tRNAid eq $hits->[0]->[1])
                    {
                        foreach $cur_iso_cm (sort @$models)
                        {
                            if ($iso_index_ct < scalar(@$hits))
                            {
                                if ($cur_iso_cm eq $hits->[$iso_index_ct]->[0])
                                {
                                    if (!defined $found_models{$cur_iso_cm})
                                    {
                                        $line .= "\t".$hits->[$iso_index_ct]->[5];
                                        $found_models{$cur_iso_cm} = 1;
                                    }
                                    $iso_index_ct++;
                                }
                                elsif ($cur_iso_cm lt $hits->[$iso_index_ct]->[0])
                                {
                                    $line .= "\t";
                                }
                                elsif ($cur_iso_cm gt $hits->[$iso_index_ct]->[0])
                                {
                                    $iso_index_ct++;
                                }
                            }
                            else
                            {
                                $line .= "\t";
                            }
                        }
                        $iso_results_out->write_line($line);

                        $hits = $cmscanResults->get_next_tab_seq_hits();
                        ($tRNAid, $line) = $iso_int_results->get_next_record();
                    }
                    elsif ($tRNAid lt $hits->[0]->[1])
                    {
                        for (my $i = 0; $i < scalar(@$models); $i++)
                        {
                            $line .= "\t";
                        }
                        $iso_results_out->write_line($line);
                        ($tRNAid, $line) = $iso_int_results->get_next_record();
                    }
                    elsif ($tRNAid gt $hits->[0]->[1])
                    {
                        $hits = $cmscanResults->get_next_tab_seq_hits();
                    }
                }
                else
                {
                    for (my $i = 0; $i < scalar(@$models); $i++)
                    {
                        $line .= "\t";
                    }
                    $iso_results_out->write_line($line);
                    ($tRNAid, $line) = $iso_int_results->get_next_record();
                }
            }
            $cmscanResults->close_file();
        }
        $iso_results_out->close_file();
        $iso_int_results->close_file();
        copy($iso_results_out->file_name(), $iso_int_results->file_name());
    }
}

sub get_models_from_db
{
    my $self = shift;
    my ($cm_db) = @_;
    my $line = "";
    my @models = ();
    my $start = 0;
    
    open(FILE_IN, "$cm_db") or die "Fail to open $cm_db\n";
    while ($line = <FILE_IN>)
    {
		chomp($line);
        if ($line =~ /^INFERNAL/)
        {
			$start = 1;
		}		
        elsif ($start and $line =~ /^NAME\s+(\S+)/)
        {
			push(@models, $1);
            $start = 0;
		}
	}
	close(FILE_IN);
    
    my @sorted_models = sort @models;
    return \@sorted_models;
}

sub rescore_tRNA_cove
{
    my $self = shift;
    my ($global_vars, $sp_int_results) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $opts = $global_vars->{options};
    my $log = $global_vars->{log_file};
    my $tmp_trnaseq_file = $global_constants->get("tmp_trnaseq_file");

#    $sp_int_results->sort_records("tRNAscan_id");

    &write_tRNAs($tmp_trnaseq_file, $sp_int_results, 1, 1, "");

    # run cove
    my $cove_output_file = $global_constants->get("temp_dir")."/tscan$$"."_cove.out";
    my $coves_cmd = $self->{coves_bin}." -s ".$self->{cove_cm_file_path}." ".$tmp_trnaseq_file." > ".$cove_output_file;
    system($coves_cmd);
    if ($? != 0)
    {
        $log->error("Coves-SE for rescoring tRNAs cannot be completed successfully. (Exit Code: $?)");
        return ("Error", "", -1);
    }
   
    my $line = "";
    my ($score, $trna_id);
    my $index = -1;
    my $cove_model = undef;
    my $new_int_result_file = $global_constants->get("temp_dir")."/tscan$$"."_sp_rescore.out";
    my $int_results_out = tRNAscanSE::IntResultFile->new;
    my $tRNA = tRNAscanSE::tRNA->new; 
    
    $int_results_out->file_name($new_int_result_file);
    &open_for_read(\*FILE_H, $cove_output_file);
    if ($sp_int_results->open_file("read") && $int_results_out->open_file("write"))
    {
        while ($line = <FILE_H>)
        {
            chomp($line);
            if ($line =~ /^\s*([0-9\.]+) bits : (.+)$/)
            {
                $score = $1;
                $trna_id = $2;
                $index = $sp_int_results->bsearch_tRNAscan_id($trna_id);
                if ($index > -1)
                {
                    $sp_int_results->get_tRNA($sp_int_results->get_pos($index), $tRNA);
                    $tRNA->tRNAscan_id($tRNA->seqname().".tRNA".$tRNA->id()."-".$tRNA->isotype().$tRNA->anticodon());
                    if ($tRNA->model() ne $self->{Pselc_cm_file_path} and $tRNA->model() ne $self->{Eselc_cm_file_path} and
                        $tRNA->model() ne $self->{isotype_cm_file_paths}->{SeC})
                    {
                        $cove_model = $tRNA->get_domain_model("cove");
                        if (!defined $cove_model)
                        {
                            $tRNA->set_domain_model("cove", $score);
                        }
                        elsif ($score > $cove_model->{score})
                        {
                            $tRNA->update_domain_model("cove", $score, $score, 0, 0);
                        }
                    }
                    $tRNA->set_default_scores();
                    $int_results_out->write_tRNA($tRNA);
                }            
            }        
        }
        $sp_int_results->close_file();
        $int_results_out->close_file();
        $opts->secondpass_int_result_file($new_int_result_file);
        $sp_int_results->clear();
        $sp_int_results->file_name($new_int_result_file);
    }
    close(FILE_H);
}

sub rescore_tRNA
{
    my $self = shift;
    my ($global_vars, $cm_tRNA, $prescan_tRNA) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $opts = $global_vars->{options};
    my $log = $global_vars->{log_file};

    my $trna_file = tRNAscanSE::Sequence->new;
    $trna_file->open_file($global_constants->get("tmp_trnaseq_file"), "write");
    $trna_file->set_seq_info($cm_tRNA->seqname().".trna".&pad_num($cm_tRNA->id(), 6), $cm_tRNA->seqname(), length($cm_tRNA->seq()), $cm_tRNA->seq());
    $trna_file->write_fasta();
    $trna_file->close_file();

    my $cm_file = $self->{main_cm_file_path}->{$cm_tRNA->model()};
    my $cms_output_file = $global_constants->get("temp_dir")."/tscan$$"."_cm_rescore.out";

    my ($hit_score, $hit_start, $hit_end, $hit_ct) = 
        $self->cmsearch_scoring($opts, $log, $global_constants->get("tmp_trnaseq_file"), $prescan_tRNA->tRNAscan_id(), $cm_file, $cms_output_file);
    $cm_tRNA->update_domain_model("infernal", $hit_score, $hit_score, 0, 0);
}

# Runs Infernal cmsearch, and returns the score
sub cmsearch_scoring
{
    my $self = shift;
    my ($opts, $log, $tmp_trnaseq_file, $seq_name, $cur_cm_file, $cms_output_file) = @_;

    my (@cms_lines, $subseq_start, $subseq_end, $evalue, $score, $besthit_score, $line,
        $besthit_start, $besthit_end, $hit_ct);

    my $scan_flag = 0;
    if ($opts->arch_mode() or $opts->bact_mode())
    {
        $scan_flag = 7;
    }
    if (!$self->exec_cmsearch($cur_cm_file, $scan_flag, $tmp_trnaseq_file, $seq_name, $cms_output_file, $log, 0))
    {
        return 0;
    }

    @cms_lines = split(/\n/, `cat $cms_output_file`);
    $hit_ct = 0;
    $score = 0;
    $besthit_score = 0;

    foreach $line (@cms_lines)
    {
        if ($line =~ /^\s+\(\d+\)\s+\S+\s+([e0-9.\-]+)\s+([0-9.\-]+)\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)\s+([+-])\s+\S+\s+\S+\s+(\S+)\s+([0-9.]+)/)
        {
            $evalue = $1;
            $score = $2;
            $subseq_start  = $5;
            $subseq_end    = $6;
            $hit_ct++;
            
            if ($score > $besthit_score)
            {
                $besthit_score  = $score;   
                $besthit_start  = $subseq_start;
                $besthit_end    = $subseq_end;
            }
        }
    }
    return ($besthit_score, $besthit_start, $besthit_end, $hit_ct); 
}

sub run_first_pass_cmsearch
{
    my $self = shift;
    my ($global_vars, $seq_name) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $log = $global_vars->{log_file};
    
    my $scan_flag = 3;
    my $cms_output_file = "";
    my $cms_merge_file = $global_constants->get("temp_dir")."/tscan$$"."_fp_cm_merge.out";
    my $arrayCMscanResults = tRNAscanSE::ArrayCMscanResults->new;
    my $file_idx = -1;

    foreach my $cur_cm (sort keys %{$self->{main_cm_file_path}})
    {
        $cms_output_file = $global_constants->get("temp_dir")."/tscan$$"."_fp_cm_".$cur_cm.".out";
        if ($self->exec_cmsearch($self->{main_cm_file_path}->{$cur_cm}, $scan_flag, $global_constants->get("tmp_fa"), $seq_name, $cms_output_file, $log, $self->{infernal_fp_cutoff}))
        {
            $file_idx = $arrayCMscanResults->add_file($cms_output_file, $cur_cm);
            $arrayCMscanResults->merge_result_file($file_idx, 0);
        }
        else
        {
            return "";
        }
    }
    $arrayCMscanResults->write_merge_file($cms_merge_file, 1);
    my $count = $arrayCMscanResults->get_result_count();
    
    return ($cms_merge_file, $count);
}

sub process_fp_cmsearch_hits
{
    my $self = shift;
    my ($global_vars, $start_index, $cms_merge_file, $result_count) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $fp_tRNAs = $global_vars->{fp_tRNAs};
    my $stats = $global_vars->{stats};
    
    my $line = "";
    my @columns = ();
    my $trna = undef;
    my $i = 0;
    my $ct = 0;
    
    &open_for_read(\*FILE_H, $cms_merge_file);
    while ($line = <FILE_H>)
    {
        $ct++;
        chomp($line);
        @columns = split(/\t/, $line);

        $trna = tRNAscanSE::tRNA->new;
        $trna->seqname($columns[0]);
        $trna->score($columns[4]);
        if ($columns[3] eq "+")
        {
            $trna->start($columns[1] + $start_index);
            $trna->end($columns[2] + $start_index);
        }
        else
        {
            $trna->start($columns[2] + $start_index);
            $trna->end($columns[1] + $start_index);
        }
        $trna->strand($columns[3]);
        $trna->isotype("???");
        $trna->anticodon("???");
        $trna->ss($columns[5]);
        $trna->seq($columns[6]);
        $trna->model($columns[10]);
        
        if ($trna->start() < $trna->end())
        {
            $trna->position($trna->start());                
        }
        else
        { 
            $trna->position($global_constants->get("really_big_number") - $trna->start() + 1);
        }

        $trna->hit_source(0);
        
        my $merge = 0;
        if ($ct < 3 or $ct > ($result_count - 3))
        {
            $merge = $self->merge_repeat_hit($global_vars->{stats}, $global_vars->{fp_tRNAs}, $trna);
        }
        if (!$merge)
        {
            # insert non-redundant hit in order
            # 'Merge_repeat_hits' depends on list being in order
            $i = &max(0, $i - 10);
            while (($i < $fp_tRNAs->get_count()) && ($fp_tRNAs->get($i)->position() < $trna->position()))
            {
                $i++;
            }
                   
            $fp_tRNAs->insert($trna, $i);
            $stats->increment_trnatotal();
        }
    }
    close(FILE_H);
}

# check current hit for redundancy against all previous hits in hitlist
#
# if it IS a repeat, merge it with overlapping hit and return 1
# if it doesn't overlap with any hits, return 0
sub merge_repeat_hit
{
    my $self = shift;
    my ($stats, $fp_tRNAs, $trna) = @_;

    for (my $i = 0; $i < $fp_tRNAs->get_count(); $i++)
    {    
        if ($trna->strand() eq "+")
        {
            if (($fp_tRNAs->get($i)->strand() eq "+") &&
                (&seg_overlap($trna->start(), $trna->end(), $fp_tRNAs->get($i)->start(), $fp_tRNAs->get($i)->end(), 0))) 
            {
                $fp_tRNAs->get($i)->start(&min($trna->start(), $fp_tRNAs->get($i)->start()));
                $fp_tRNAs->get($i)->end(&max($trna->end(), $fp_tRNAs->get($i)->end()));
                $fp_tRNAs->get($i)->hit_source($fp_tRNAs->get($i)->hit_source());
                $fp_tRNAs->get($i)->isotype($trna->isotype());
                $fp_tRNAs->get($i)->score($trna->score());
    
                # check to see if extended endpoint overlaps i+1 hit's start boundary
                # if so, combine hit[i] and hit[i+1] into one hit and delete hit[i+1]
                if (($i != ($fp_tRNAs->get_count() - 1)) and ($fp_tRNAs->get($i+1)->strand() eq "+")
                    and ($fp_tRNAs->get($i)->end() >= $fp_tRNAs->get($i+1)->start())) 
                {
                    $fp_tRNAs->get($i)->end(&max($fp_tRNAs->get($i)->end(), $fp_tRNAs->get($i+1)->end()));
                    $fp_tRNAs->get($i)->hit_source($fp_tRNAs->get($i)->hit_source() | $fp_tRNAs->get($i+1)->hit_source());

                    $fp_tRNAs->remove($i+1);
                    $stats->decrement_trnatotal();
                }   
                return 1;                                 # exit loop immediately
            }
        }
        else         # else (antisense) strand 
        {                
            if (($fp_tRNAs->get($i)->strand() eq "-") &&
                (&seg_overlap($trna->end(), $trna->start(), $fp_tRNAs->get($i)->end(), $fp_tRNAs->get($i)->start(), 0))) 
            {
                $fp_tRNAs->get($i)->start(&max($trna->start(), $fp_tRNAs->get($i)->start()));
                $fp_tRNAs->get($i)->end(&min($trna->end(), $fp_tRNAs->get($i)->end()));
                $fp_tRNAs->get($i)->hit_source($fp_tRNAs->get($i)->hit_source() | $self->{eufind_mask});
                $fp_tRNAs->get($i)->isotype($trna->isotype());
                $fp_tRNAs->get($i)->score($trna->score());

                if (($i != ($fp_tRNAs->get_count() - 1)) and
                    ($fp_tRNAs->get($i)->end() <= $fp_tRNAs->get($i+1)->start()))
                {
                    $fp_tRNAs->get($i)->end(&min($fp_tRNAs->get($i)->end(), $fp_tRNAs->get($i+1)->end()));
                    $fp_tRNAs->get($i)->hit_source($fp_tRNAs->get($i)->hit_source() | $fp_tRNAs->get($i+1)->hit_source());

                    $fp_tRNAs->remove($i+1);
                    $stats->decrement_trnatotal();
                }
                return 1;                                 # exit loop immediately
            }
        } # else (antisense) strand
        
    } # for each (hit)                        

    return 0;                                             # current hit is not a repeat
}

sub first_pass_scan
{
    my $self = shift;
    my ($global_vars, $start_index, $seq_name) = @_;
    my $opts = $global_vars->{options};
    my $log = $global_vars->{log_file};
    
    if (!$opts->numt_mode())
    {
        my ($cms_merge_file, $count) = $self->run_first_pass_cmsearch($global_vars, $seq_name);
        if ($cms_merge_file eq "")
        {
            $log->error("Fail to run infernal for first pass scan");
            return 0;
        }
        else
        {
            $self->process_fp_cmsearch_hits($global_vars, $start_index, $cms_merge_file, $count);
        }
    }
    
    return 1;
}

sub analyze_alternate
{
    my $self = shift;
    my ($global_vars, $seqinfo_flag, $seq_ct, $seq_name, $r_curseq_trnact) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $opts = $global_vars->{options};
    my $stats = $global_vars->{stats};
    my $log = $global_vars->{log_file};
    my $fp_result_file = $global_vars->{fp_result_file};
    my $arrayCMscanResults = tRNAscanSE::ArrayCMscanResults->new;
    my $prescan_tRNA = tRNAscanSE::tRNA->new;
    my $cm_tRNA = tRNAscanSE::tRNA->new;
    
    my $over_cutoff = 0;
    my $cms_hit_count = 0;
    my ($trnaDesc, $fulltrnaDesc, $subseq_start, $subseq_end);
    my ($cms_confirmed_ct, $cur_cm_file);    
    my $cm_hit = {};
    
    $cms_confirmed_ct = 0;

    if (scalar(keys %{$self->{main_cm_file_path}}) == 0)
    {
        $log->error("No covariance models are found for scanning $seq_name. Please check the configuration file.");
        return 0;
    }
    
    my ($cms_merge_file, $count) = $self->run_cmsearch($global_vars, $seq_name, $arrayCMscanResults, 0);
    if ($cms_merge_file eq "")
    {
        $log->error("Fail to run infernal for $seq_name");
        return 0;
    }
    
    $arrayCMscanResults->open_file($cms_merge_file);
    $fp_result_file->reset_current_seq();
    $fp_result_file->get_next_tRNA_candidate($opts, $seqinfo_flag, $seq_ct, $prescan_tRNA);
    $arrayCMscanResults->get_next_cmsearch_hit($cm_tRNA);
    while ($cm_tRNA->seqname() ne "")
    {
        $over_cutoff = 0;
        while (($cm_tRNA->seqname() ne "") and ($cm_tRNA->seqname() eq $prescan_tRNA->seqname().".t".&pad_num($prescan_tRNA->id(), 6)))
        {
            $cms_hit_count++;
            
            $subseq_start = $cm_tRNA->start();
            $subseq_end = $cm_tRNA->end();
            $cm_tRNA->start($cm_tRNA->start() + $prescan_tRNA->start() - 1);
            $cm_tRNA->end($cm_tRNA->end() + $prescan_tRNA->start() - 1);

            if ($cm_tRNA->score() < $self->{cm_cutoff})
            {
                $log->broadcast("Low cmsearch score for ".$prescan_tRNA->tRNAscan_id().".$cms_hit_count: ".$cm_tRNA->score());
                $trnaDesc .= "(CMSearch Hit#$cms_hit_count: ".$cm_tRNA->start()."-".$cm_tRNA->end().",".
                    " Sc: ".$cm_tRNA->score().",  Len: ".(abs($cm_tRNA->start() - $cm_tRNA->end()) + 1).") ";
            }
            else
            {
                $over_cutoff++;
                $$r_curseq_trnact++;
                
                $cm_tRNA->id($$r_curseq_trnact);
                $cm_tRNA->seqname($prescan_tRNA->seqname());
                $cm_tRNA->set_domain_model("infernal", $cm_tRNA->score());
                if ($cm_tRNA->strand() eq "-")
                {
                    my $temp = $cm_tRNA->start();
                    $cm_tRNA->start($cm_tRNA->end());
                    $cm_tRNA->end($temp);
                }
                
                if ((length($prescan_tRNA->seq()) >= $subseq_end + 2) &&
                    (uc(substr($prescan_tRNA->seq(), $subseq_end, 3)) eq "CCA") &&
                    (substr($cm_tRNA->ss(), length($cm_tRNA->ss()) - 4) ne "...."))
                {
                    $subseq_end += 3;
                    if ($cm_tRNA->strand() eq "+")
                    {
                        $cm_tRNA->end($cm_tRNA->end() + 3);
                    }
                    else
                    {
                        $cm_tRNA->start($cm_tRNA->start() - 3);
                    }
                    $cm_tRNA->seq($cm_tRNA->seq()."CCA");
                    $cm_tRNA->ss($cm_tRNA->ss()."...");
                }
                
                $cm_tRNA->upstream($prescan_tRNA->upstream());
                $cm_tRNA->downstream($prescan_tRNA->downstream());
                if ($subseq_start > 1)
                {
                    $cm_tRNA->upstream($cm_tRNA->upstream() . substr($prescan_tRNA->seq(), 0, $subseq_start - 1));
                }
                if ($subseq_end < $prescan_tRNA->len())
                {
                    $cm_tRNA->downstream(substr($prescan_tRNA->seq(), $subseq_end) . $cm_tRNA->downstream());
                }
                my $upstream_len   = $global_constants->get("upstream_len");
                my $downstream_len = $global_constants->get("downstream_len");
                
                $cm_tRNA->upstream(substr($cm_tRNA->upstream(), &max((length($cm_tRNA->upstream()) - $upstream_len), 0)));
                $cm_tRNA->downstream(substr($cm_tRNA->downstream(), 0, &min(length($cm_tRNA->downstream()), $downstream_len)));
                $self->decode_tRNA_properties($global_vars, $cm_tRNA, $prescan_tRNA, $self->{main_cm_file_path}->{$cm_tRNA->model()});
                
                $cm_tRNA->tRNAscan_id($cm_tRNA->seqname().".tRNA".$$r_curseq_trnact."-".$cm_tRNA->isotype().$cm_tRNA->anticodon());
                $cm_tRNA->set_mature_tRNA();
                $cm_tRNA->src_seqlen($prescan_tRNA->src_seqlen());
                $cm_tRNA->src_seqid($cm_tRNA->seqname());
                $cm_tRNA->ordered_seqname($prescan_tRNA->src_seqid());
                $cm_tRNA->set_default_scores();
                
                $cms_confirmed_ct++;
                $global_vars->{sp_int_results}->write_tRNA($cm_tRNA);
            }
                        
            $arrayCMscanResults->get_next_cmsearch_hit($cm_tRNA);
        }
        
        if ($over_cutoff == 0)
        {
            if ((!$opts->results_to_stdout()) && ($opts->eufind_mode() || $opts->tscan_mode() || $opts->use_prev_ts_run()))
            {
                $log->broadcast("CMSearch score(s) below cutoff for ".$prescan_tRNA->tRNAscan_id().". Skipping...");
            }
            if ($opts->save_falsepos())
            {
                $fulltrnaDesc = "(Fp Hit: ".$prescan_tRNA->start()."-".$prescan_tRNA->end().", ".
                    (abs($prescan_tRNA->start() - $prescan_tRNA->end()) + 1)." bp, Src: ".$prescan_tRNA->hit_source().") ".$trnaDesc;
            
                $stats->increment_fpos_base_ct(length($prescan_tRNA->seq()));          
                &write_tRNA($opts->falsepos_file(), $prescan_tRNA->tRNAscan_id(), $fulltrnaDesc, $prescan_tRNA->seq(), 0);
            }           
        }
        
        $fp_result_file->get_next_tRNA_candidate($opts, $seqinfo_flag, $seq_ct, $prescan_tRNA);
    }
    $arrayCMscanResults->close_file();

    return $cms_confirmed_ct;
}

sub analyze_mito
{
    my $self = shift;
    my ($global_vars, $seqinfo_flag, $seq_ct, $seq_name, $r_curseq_trnact) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $opts = $global_vars->{options};
    my $stats = $global_vars->{stats};
    my $log = $global_vars->{log_file};
    my $fp_result_file = $global_vars->{fp_result_file};
    my $arrayCMscanResults = tRNAscanSE::ArrayCMscanResults->new;
    my $prescan_tRNA = tRNAscanSE::tRNA->new;
    my $cm_tRNA = tRNAscanSE::tRNA->new;
    
    my $over_cutoff = 0;
    my $cms_hit_count = 0;
    my ($trnaDesc, $fulltrnaDesc, $subseq_start, $subseq_end);
    my ($cms_confirmed_ct, $cur_cm_file);    
    my $cm_hit = {};
    
    $cms_confirmed_ct = 0;

    if (scalar(keys %{$self->{main_cm_file_path}}) == 0)
    {
        $log->error("No covariance models are found for scanning $seq_name. Please check the configuration file.");
        return 0;
    }
    
    my ($cms_merge_file, $count) = $self->run_cmsearch($global_vars, $seq_name, $arrayCMscanResults, 10);
    if ($cms_merge_file eq "")
    {
        $log->error("Fail to run infernal for $seq_name");
        return 0;
    }
    
    $arrayCMscanResults->open_file($cms_merge_file);
    $fp_result_file->reset_current_seq();
    $fp_result_file->get_next_tRNA_candidate($opts, $seqinfo_flag, $seq_ct, $prescan_tRNA);
    $arrayCMscanResults->get_next_cmsearch_hit($cm_tRNA);
    while ($cm_tRNA->seqname() ne "")
    {
        $over_cutoff = 0;
        while (($cm_tRNA->seqname() ne "") and ($cm_tRNA->seqname() eq $prescan_tRNA->seqname().".t".&pad_num($prescan_tRNA->id(), 6)))
        {
            $cms_hit_count++;
            
            $subseq_start = $cm_tRNA->start();
            $subseq_end = $cm_tRNA->end();
            $cm_tRNA->start($cm_tRNA->start() + $prescan_tRNA->start() - 1);
            $cm_tRNA->end($cm_tRNA->end() + $prescan_tRNA->start() - 1);

            if ($cm_tRNA->score() < $self->{cm_cutoff})
            {
                $log->broadcast("Low cmsearch score for ".$prescan_tRNA->tRNAscan_id().".$cms_hit_count: ".$cm_tRNA->score());
                $trnaDesc .= "(CMSearch Hit#$cms_hit_count: ".$cm_tRNA->start()."-".$cm_tRNA->end().",".
                    " Sc: ".$cm_tRNA->score().",  Len: ".(abs($cm_tRNA->start() - $cm_tRNA->end()) + 1).") ";
            }
            else
            {
                $over_cutoff++;
                $$r_curseq_trnact++;
                
                $cm_tRNA->id($$r_curseq_trnact);
                $cm_tRNA->seqname($prescan_tRNA->seqname());
                $cm_tRNA->set_domain_model("infernal", $cm_tRNA->score());
                if ($cm_tRNA->strand() eq "-")
                {
                    my $temp = $cm_tRNA->start();
                    $cm_tRNA->start($cm_tRNA->end());
                    $cm_tRNA->end($temp);
                }
                
                if ((length($prescan_tRNA->seq()) >= $subseq_end + 2) &&
                    (uc(substr($prescan_tRNA->seq(), $subseq_end, 3)) eq "CCA") &&
                    (substr($cm_tRNA->ss(), length($cm_tRNA->ss()) - 4) ne "...."))
                {
                    $subseq_end += 3;
                    if ($cm_tRNA->strand() eq "+")
                    {
                        $cm_tRNA->end($cm_tRNA->end() + 3);
                    }
                    else
                    {
                        $cm_tRNA->start($cm_tRNA->start() - 3);
                    }
                    $cm_tRNA->seq($cm_tRNA->seq()."CCA");
                    $cm_tRNA->ss($cm_tRNA->ss()."...");
                }
                
                $cm_tRNA->upstream($prescan_tRNA->upstream());
                $cm_tRNA->downstream($prescan_tRNA->downstream());
                if ($subseq_start > 1)
                {
                    $cm_tRNA->upstream($cm_tRNA->upstream() . substr($prescan_tRNA->seq(), 0, $subseq_start - 1));
                }
                if ($subseq_end < $prescan_tRNA->len())
                {
                    $cm_tRNA->downstream(substr($prescan_tRNA->seq(), $subseq_end) . $cm_tRNA->downstream());
                }
                my $upstream_len   = $global_constants->get("upstream_len");
                my $downstream_len = $global_constants->get("downstream_len");
                
                $cm_tRNA->upstream(substr($cm_tRNA->upstream(), &max((length($cm_tRNA->upstream()) - $upstream_len), 0)));
                $cm_tRNA->downstream(substr($cm_tRNA->downstream(), 0, &min(length($cm_tRNA->downstream()), $downstream_len)));
                $self->decode_mito_tRNA_properties($global_vars, $cm_tRNA, $prescan_tRNA, $self->{main_cm_file_path}->{$cm_tRNA->model()});
                
                $cm_tRNA->tRNAscan_id($cm_tRNA->seqname().".tRNA".$$r_curseq_trnact."-".$cm_tRNA->isotype().$cm_tRNA->anticodon());
                $cm_tRNA->set_mature_tRNA();
                $cm_tRNA->src_seqlen($prescan_tRNA->src_seqlen());
                $cm_tRNA->src_seqid($cm_tRNA->seqname());
                $cm_tRNA->ordered_seqname($prescan_tRNA->src_seqid());
                $cm_tRNA->set_default_scores();
                
                $cms_confirmed_ct++;
                $global_vars->{sp_int_results}->write_tRNA($cm_tRNA);
            }
                        
            $arrayCMscanResults->get_next_cmsearch_hit($cm_tRNA);
        }
        
        if ($over_cutoff == 0)
        {
            if ((!$opts->results_to_stdout()) && ($opts->eufind_mode() || $opts->tscan_mode() || $opts->use_prev_ts_run()))
            {
                $log->broadcast("CMSearch score(s) below cutoff for ".$prescan_tRNA->tRNAscan_id().". Skipping...");
            }
            if ($opts->save_falsepos())
            {
                $fulltrnaDesc = "(Fp Hit: ".$prescan_tRNA->start()."-".$prescan_tRNA->end().", ".
                    (abs($prescan_tRNA->start() - $prescan_tRNA->end()) + 1)." bp, Src: ".$prescan_tRNA->hit_source().") ".$trnaDesc;
            
                $stats->increment_fpos_base_ct(length($prescan_tRNA->seq()));          
                &write_tRNA($opts->falsepos_file(), $prescan_tRNA->tRNAscan_id(), $fulltrnaDesc, $prescan_tRNA->seq(), 0);
            }           
        }
        
        $fp_result_file->get_next_tRNA_candidate($opts, $seqinfo_flag, $seq_ct, $prescan_tRNA);
    }
    $arrayCMscanResults->close_file();

    return $cms_confirmed_ct;
}

sub analyze_with_cmsearch
{ 
    my $self = shift;
    my ($global_vars, $seqinfo_flag, $seq_ct, $seq_name, $r_curseq_trnact) = @_;
    my $global_constants = $global_vars->{global_constants};
    my $opts = $global_vars->{options};
    my $stats = $global_vars->{stats};
    my $log = $global_vars->{log_file};
    my $fp_result_file = $global_vars->{fp_result_file};
    my $arrayCMscanResults = tRNAscanSE::ArrayCMscanResults->new;
    my $prescan_tRNA = tRNAscanSE::tRNA->new;
    my $cm_tRNA = tRNAscanSE::tRNA->new;
    my $flanking_exist = 0;
    
    my $over_cutoff = 0;
    my $cms_hit_count = 0;
    my $rescore = 0;
    my ($trnaDesc, $fulltrnaDesc, $subseq_start, $subseq_end);
    my ($cms_confirmed_ct, $cur_cm_file);    
    my $cm_hit = {};
    
    $cms_confirmed_ct = 0;

    my ($cms_merge_file, $count) = $self->run_cmsearch($global_vars, $seq_name, $arrayCMscanResults, 0);
    if ($cms_merge_file eq "")
    {
        $log->error("Fail to run infernal for $seq_name");
        return 0;
    }

    $arrayCMscanResults->open_file($cms_merge_file);
    $fp_result_file->reset_current_seq();
    $fp_result_file->get_next_tRNA_candidate($opts, $seqinfo_flag, $seq_ct, $prescan_tRNA);
    if ($fp_result_file->open_flanking("read"))
    {
        $fp_result_file->read_tRNA_flanking($prescan_tRNA);
        $flanking_exist = 1;
    }
    $arrayCMscanResults->get_next_cmsearch_hit($cm_tRNA);
    while ($cm_tRNA->seqname() ne "")
    {
        $over_cutoff = 0;
        while (($cm_tRNA->seqname() ne "") and ($cm_tRNA->seqname() eq $prescan_tRNA->seqname().".t".&pad_num($prescan_tRNA->id(), 6)))
        {
            $rescore = 0;
            $cms_hit_count++;
            if ($prescan_tRNA->hit_source() ne "")
            {
                $cm_tRNA->hit_source($prescan_tRNA->hit_source());
            }
            $subseq_start = $cm_tRNA->start();
            $subseq_end = $cm_tRNA->end();
            if ($opts->infernal_fp())
            {
                $cm_tRNA->strand($prescan_tRNA->strand());
                if ($prescan_tRNA->strand() eq "+")
                {
                    $cm_tRNA->start($cm_tRNA->start() + $prescan_tRNA->start() - 1);
                    $cm_tRNA->end($cm_tRNA->end() + $prescan_tRNA->start() - 1);
                }
                else
                {
                    if ($prescan_tRNA->start() > $prescan_tRNA->end())
                    {
                        $cm_tRNA->start($prescan_tRNA->start() - $cm_tRNA->start() + 1);	
                        $cm_tRNA->end($prescan_tRNA->start() - $cm_tRNA->end() + 1);
                    }
                    else
                    {
                        $cm_tRNA->start($prescan_tRNA->end() - $cm_tRNA->start() + 1);	
                        $cm_tRNA->end($prescan_tRNA->end() - $cm_tRNA->end() + 1);
                    }
                }
            }
            else
            {
                $cm_tRNA->start($cm_tRNA->start() + $prescan_tRNA->start() - 1);
                $cm_tRNA->end($cm_tRNA->end() + $prescan_tRNA->start() - 1);
            }

            if ($cm_tRNA->score() < $self->{cm_cutoff})
            {
                $log->broadcast("Low cmsearch score for ".$prescan_tRNA->tRNAscan_id().".$cms_hit_count: ".$cm_tRNA->score());
                $trnaDesc .= "(CMSearch Hit#$cms_hit_count: ".$cm_tRNA->start()."-".$cm_tRNA->end().",".
                    " Sc: ".$cm_tRNA->score().",  Len: ".(abs($cm_tRNA->start() - $cm_tRNA->end()) + 1).") ";
            }
            else
            {
                $over_cutoff++;
                $$r_curseq_trnact++;
                
                $cm_tRNA->id($$r_curseq_trnact);
                $cm_tRNA->seqname($prescan_tRNA->seqname());
                $cm_tRNA->set_domain_model("infernal", $cm_tRNA->score());
                if ($cm_tRNA->strand() eq "-")
                {
                    my $temp = $cm_tRNA->start();
                    $cm_tRNA->start($cm_tRNA->end());
                    $cm_tRNA->end($temp);
                }
                
                if ((length($prescan_tRNA->seq()) >= $subseq_end + 2) &&
                    (uc(substr($prescan_tRNA->seq(), $subseq_end, 3)) eq "CCA") &&
                    (substr($cm_tRNA->ss(), length($cm_tRNA->ss()) - 4) ne "....") &&
                    ($self->{main_cm_file_path}->{Domain} ne $global_constants->get_subvalue("cm", "eukaryota")) &&
                    ($self->{main_cm_file_path}->{Domain} ne $global_constants->get_subvalue("cm", "general")))
                {
                    $subseq_end += 3;
                    if ($cm_tRNA->strand() eq "+")
                    {
                        $cm_tRNA->end($cm_tRNA->end() + 3);
                    }
                    else
                    {
                        $cm_tRNA->start($cm_tRNA->start() - 3);
                    }
                    $cm_tRNA->seq($cm_tRNA->seq()."CCA");
                    $cm_tRNA->ss($cm_tRNA->ss()."...");
                    $rescore = 1;
                }

                if ((uc(substr($cm_tRNA->seq(), length($cm_tRNA->seq()) - 3)) ne "CCA") &&
                    (substr($cm_tRNA->ss(), length($cm_tRNA->ss()) - 4) eq "....") &&
                    ($self->{main_cm_file_path}->{Domain} ne $global_constants->get_subvalue("cm", "eukaryota")) &&
                    ($self->{main_cm_file_path}->{Domain} ne $global_constants->get_subvalue("cm", "general")))
                {
                    $subseq_end -= 3;
                    if ($cm_tRNA->strand() eq "+")
                    {
                        $cm_tRNA->end($cm_tRNA->end() - 3);
                    }
                    else
                    {
                        $cm_tRNA->start($cm_tRNA->start() + 3);
                    }
                    $cm_tRNA->seq(substr($cm_tRNA->seq(), 0, length($cm_tRNA->seq()) - 3));
                    $cm_tRNA->ss(substr($cm_tRNA->ss(), 0, length($cm_tRNA->ss()) - 3));
                    $rescore = 1;
                }
                elsif ((uc(substr($cm_tRNA->seq(), length($cm_tRNA->seq()) - 3)) eq "CCA") &&
                    (((substr($cm_tRNA->ss(), length($cm_tRNA->ss()) - 6) eq "<.....") &&
                    (substr($cm_tRNA->seq(), length($cm_tRNA->seq()) - 5, 2) =~ /[acgtn][ACGTN]/)) ||
                    ((substr($cm_tRNA->ss(), length($cm_tRNA->ss()) - 7) eq "<......") &&
                    (substr($cm_tRNA->seq(), length($cm_tRNA->seq()) - 6, 3) =~ /[ACGTN][acgtn][ACGTN]/)) ||
                    ((substr($cm_tRNA->ss(), length($cm_tRNA->ss()) - 7) eq "<......") &&
                    (substr($cm_tRNA->seq(), length($cm_tRNA->seq()) - 6, 3) =~ /[acgtn]{2}[ACGTN]/))) &&
                    ($self->{main_cm_file_path}->{Domain} ne $global_constants->get_subvalue("cm", "eukaryota")) &&
                    ($self->{main_cm_file_path}->{Domain} ne $global_constants->get_subvalue("cm", "general")))
                {
                    my $trim_len = 4;
                    if (substr($cm_tRNA->ss(), length($cm_tRNA->ss()) - 7) eq "<......" && substr($cm_tRNA->seq(), length($cm_tRNA->seq()) - 6, 3) =~ /[acgtn]{2}[ACGTN]/)
                    {
						$trim_len = 5;
					}
					
                    $subseq_end -= $trim_len;
                    if ($cm_tRNA->strand() eq "+")
                    {
                        $cm_tRNA->end($cm_tRNA->end() - $trim_len);
                    }
                    else
                    {
                        $cm_tRNA->start($cm_tRNA->start() + $trim_len);
                    }
                    $cm_tRNA->seq(substr($cm_tRNA->seq(), 0, length($cm_tRNA->seq()) - $trim_len));
                    $cm_tRNA->seq(substr($cm_tRNA->seq(), 0, length($cm_tRNA->seq()) - 1).uc(substr($cm_tRNA->seq(), length($cm_tRNA->seq()) - 1, 1)));
                    $cm_tRNA->ss(substr($cm_tRNA->ss(), 0, length($cm_tRNA->ss()) - $trim_len));
                    $rescore = 1;
                }
                
                if ($rescore)
                {
                    $self->rescore_tRNA($global_vars, $cm_tRNA, $prescan_tRNA);
                }
                
                $cm_tRNA->upstream($prescan_tRNA->upstream());
                $cm_tRNA->downstream($prescan_tRNA->downstream());
                if ($subseq_start > 1)
                {
                    $cm_tRNA->upstream($cm_tRNA->upstream() . substr($prescan_tRNA->seq(), 0, $subseq_start - 1));
                }
                if ($subseq_end < $prescan_tRNA->len())
                {
                    $cm_tRNA->downstream(substr($prescan_tRNA->seq(), $subseq_end) . $cm_tRNA->downstream());
                }
                if (!$opts->tscan_mode() and !$opts->eufind_mode() and !$opts->infernal_fp())
                {
                    my $upstream_len   = $global_constants->get("upstream_len");
                    my $downstream_len = $global_constants->get("downstream_len");
                    
                    $cm_tRNA->upstream(substr($cm_tRNA->upstream(), &max((length($cm_tRNA->upstream()) - $upstream_len), 0)));
                    $cm_tRNA->downstream(substr($cm_tRNA->downstream(), 0, &min(length($cm_tRNA->downstream()), $downstream_len)));
                }
                $self->decode_tRNA_properties($global_vars, $cm_tRNA, $prescan_tRNA, $self->{main_cm_file_path}->{$cm_tRNA->model()});
                
                $cm_tRNA->tRNAscan_id($cm_tRNA->seqname().".tRNA".$$r_curseq_trnact."-".$cm_tRNA->isotype().$cm_tRNA->anticodon());
                $cm_tRNA->src_seqlen($prescan_tRNA->src_seqlen());
                $cm_tRNA->src_seqid($cm_tRNA->seqname());
                $cm_tRNA->ordered_seqname($prescan_tRNA->src_seqid());
                $cm_tRNA->set_default_scores();
                $self->fix_fMet($global_vars, $cm_tRNA);
                $self->fix_His($global_vars, $cm_tRNA);
                $cm_tRNA->set_mature_tRNA();

                $cms_confirmed_ct++;
                $global_vars->{sp_int_results}->write_tRNA($cm_tRNA);
            }
                        
            $arrayCMscanResults->get_next_cmsearch_hit($cm_tRNA);
        }
        
        if ($over_cutoff == 0)
        {
            if ((!$opts->results_to_stdout()) && ($opts->eufind_mode() || $opts->tscan_mode() || $opts->use_prev_ts_run()))
            {
                $log->broadcast("CMSearch score(s) below cutoff for ".$prescan_tRNA->tRNAscan_id().". Skipping...");
            }
            if ($opts->save_falsepos())
            {
                $fulltrnaDesc = "(Fp Hit: ".$prescan_tRNA->start()."-".$prescan_tRNA->end().", ".
                    (abs($prescan_tRNA->start() - $prescan_tRNA->end()) + 1)." bp, Src: ".$prescan_tRNA->hit_source().") ".$trnaDesc;
            
                $stats->increment_fpos_base_ct(length($prescan_tRNA->seq()));          
                &write_tRNA($opts->falsepos_file(), $prescan_tRNA->tRNAscan_id(), $fulltrnaDesc, $prescan_tRNA->seq(), 0);
            }           
        }
        
        $fp_result_file->get_next_tRNA_candidate($opts, $seqinfo_flag, $seq_ct, $prescan_tRNA);
        if ($flanking_exist)
        {
			$fp_result_file->read_tRNA_flanking($prescan_tRNA);
		}		
    }
    $arrayCMscanResults->close_file();

    return $cms_confirmed_ct;
}	    

sub sort_cm_hits_by_start
{    
    my $self = shift;
    my $cms_hits = shift;
    
    my $a_start = $a->{start};
    my $b_start = $b->{start};
    
    if ($a->{strand} == 0) {
        $a_start = $a->{end};
    }
    if ($b->{strand} == 0) {
        $b_start = $b->{end};
    }
    
    return ($a->{seqname} cmp $b->{seqname} ||
            $a_start <=> $b_start);    
}

1;
