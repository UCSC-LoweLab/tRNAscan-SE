# tRNAscanSE/Stats.pm
# This class describes the statistics of each run in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2017 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#
package tRNAscan::Stats;

use strict;
use tRNAscanSE::Utils;
use tRNAscanSE::Options;
use tRNAscanSE::GeneticCode;

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
    $self->{file_name} = "";             # name of log file
    $self->{FILE_H} = undef;             # file handle
    
    $self->{fp_start_time} = +[];        # save first pass starting time
    $self->{fp_end_time} = +[];          # save first pass ending time
    $self->{sp_end_time} = +[];          # save second pass ending time
    $self->{seqs_hit} = 0;               # num seqs with at least one trna hit
    $self->{numscanned} = 0;             # total sequences scanned
    $self->{trnatotal} = 0;              # total trnas found by tscan

    $self->{first_pass_base_ct} = 0;     # no bases in all seqs in first pass scans
    $self->{fpass_trna_base_ct} = 0;     # no bases in tRNAs in first pass scans
    $self->{fpos_base_ct} = 0;           # no bases in false positive tRNAs 
    $self->{secpass_base_ct} = 0;
    $self->{coves_base_ct} = 0;
    $self->{total_secpass_ct} = 0;
}

sub file_name
{
    my $self = shift;
    if (@_) { $self->{file_name} = shift; }
    return $self->{file_name};
}

sub FILE_H
{
    my $self = shift;
    return $self->{FILE_H};
}

sub start_fp_timer
{
    my $self = shift;
    @{$self->{fp_start_time}} = (times)[0,2,1,3];
    $self->{fp_end_time} = +[];
    $self->{sp_end_time} = +[];
}

sub end_fp_timer
{
    my $self = shift;
    @{$self->{fp_end_time}} = (times)[0,2,1,3];            
}

sub start_sp_timer
{
    my $self = shift;
    if (!defined $self->{fp_end_time}->[0]) {
        push (@{$self->{fp_end_time}}, @{$self->{fp_start_time}});
    }
}

sub end_sp_timer
{
    my $self = shift;
    @{$self->{sp_end_time}} = (times)[0,2,1,3];            
}

sub seqs_hit
{
    my $self = shift;
    if (@_) { $self->{seqs_hit} = shift; }
    return $self->{seqs_hit};
}

sub increment_seqs_hit
{
    my $self = shift;
    if (@_) { $self->{seqs_hit} += shift; }
    else { $self->{seqs_hit} += 1;}
}

sub numscanned
{
    my $self = shift;
    if (@_) { $self->{numscanned} = shift; }
    return $self->{numscanned};
}

sub increment_numscanned
{
    my $self = shift;
    if (@_) { $self->{numscanned} += shift; }
    else { $self->{numscanned} += 1;}
}

sub trnatotal
{
    my $self = shift;
    if (@_) { $self->{trnatotal} = shift; }
    return $self->{trnatotal};
}

sub increment_trnatotal
{
    my $self = shift;
    if (@_) { $self->{trnatotal} += shift; }
    else { $self->{trnatotal} += 1;}
}

sub decrement_trnatotal
{
    my $self = shift;
    if (@_) { $self->{trnatotal} -= shift; }
    else { $self->{trnatotal} -= 1;}
}

sub first_pass_base_ct
{
    my $self = shift;
    if (@_) { $self->{first_pass_base_ct} = shift; }
    return $self->{first_pass_base_ct};
}

sub increment_first_pass_base_ct
{
    my $self = shift;
    if (@_) { $self->{first_pass_base_ct} += shift; }
    else { $self->{first_pass_base_ct} += 1;}
}

sub fpass_trna_base_ct
{
    my $self = shift;
    if (@_) { $self->{fpass_trna_base_ct} = shift; }
    return $self->{fpass_trna_base_ct};
}

sub fpos_base_ct
{
    my $self = shift;
    if (@_) { $self->{fpos_base_ct} = shift; }
    return $self->{fpos_base_ct};
}

sub increment_fpos_base_ct
{
    my $self = shift;
    if (@_) { $self->{fpos_base_ct} += shift; }
    else { $self->{fpos_base_ct} += 1;}
}

sub secpass_base_ct
{
    my $self = shift;
    if (@_) { $self->{secpass_base_ct} = shift; }
    return $self->{secpass_base_ct};
}

sub increment_secpass_base_ct
{
    my $self = shift;
    if (@_) { $self->{secpass_base_ct} += shift; }
    else { $self->{secpass_base_ct} += 1;}
}

sub coves_base_ct
{
    my $self = shift;
    if (@_) { $self->{coves_base_ct} = shift; }
    return $self->{coves_base_ct};
}

sub increment_coves_base_ct
{
    my $self = shift;
    if (@_) { $self->{coves_base_ct} += shift; }
    else { $self->{coves_base_ct} += 1;}
}

sub total_secpass_ct
{
    my $self = shift;
    if (@_) { $self->{total_secpass_ct} = shift; }
    return $self->{total_secpass_ct};
}

sub increment_total_secpass_ct
{
    my $self = shift;
    if (@_) { $self->{total_secpass_ct} += shift; }
    else { $self->{total_secpass_ct} += 1;}
}

sub open_file
{
    my $self = shift;
    
    my $success = 0;
    
    if ($self->{file_name} ne "")
    {
        &open_for_append(\$self->{FILE_H}, $self->{file_name});
        $success = 1;
    }
    else
    {
        die "Statistics file name is not set.\n"
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

sub save_firstpass_stats
{    
    my $self = shift;
    my $fh = $self->{FILE_H};

    print $fh "First-pass Stats:\n",
        "---------------\n";
    print $fh  "Sequences read:         $self->{numscanned}\n";
    print $fh  "Seqs w/at least 1 hit:  $self->{seqs_hit}\n"; 
    print $fh  "Bases read:             $self->{first_pass_base_ct} (x2 for both strands)\n";
    print $fh  "Bases in tRNAs:         $self->{fpass_trna_base_ct}\n";
    print $fh  "tRNAs predicted:        $self->{trnatotal}\n";
    printf $fh "Av. tRNA length:        %d\n",
        int($self->{fpass_trna_base_ct} / &max(1, $self->{trnatotal}));
    printf $fh "Script CPU time:        %.2f s\n",
        $self->{fp_end_time}->[0] - $self->{fp_start_time}->[0];
    printf $fh "Scan CPU time:          %.2f s\n",
        $self->{fp_end_time}->[1] - $self->{fp_start_time}->[1];
    printf $fh "Scan speed:             %.1f Kbp/sec\n", $self->{first_pass_base_ct}*2/
         (&max(0.001, $self->{fp_end_time}->[1] - $self->{fp_start_time}->[1]))/1000;
    print $fh "\nFirst pass search(es) ended: ",`date`,"\n";
}

sub save_final_stats
{
    my $self = shift;
    my $opts = shift;
    my $gc = shift;
    my $prescan_count = shift;
    my $r_tab_results = shift;
    my $fh = $self->{FILE_H};
    my $second_pass_label = $opts->second_pass_label();

    if ($opts->CM_mode() ne "")
    {
        print $fh "$second_pass_label Stats:\n-----------\n";
        
        if ($opts->tscan_mode() || $opts->eufind_mode() || $opts->infernal_fp())
        {
            print $fh "Candidate tRNAs read:     ", $prescan_count,"\n"; 
        }
        else
        {
            print $fh "Sequences read:           $self->{numscanned}\n";
        } 
        print $fh  "$second_pass_label","-confirmed tRNAs:     $self->{total_secpass_ct}\n";
        print $fh  "Bases scanned by $second_pass_label:  $self->{secpass_base_ct}\n";    
        printf $fh "%% seq scanned by $second_pass_label:  %2.1f %%\n",
            &min(($self->{secpass_base_ct} / &max(1, $self->{first_pass_base_ct} * 2)) * 100,100);
        printf $fh "Script CPU time:          %2.2f s\n", $self->{sp_end_time}->[0] - $self->{fp_end_time}->[0];
        printf $fh "$second_pass_label CPU time:            %2.2f s\n", $self->{sp_end_time}->[1] - $self->{fp_end_time}->[1];
        printf $fh "Scan speed:               %.1f bp/sec\n", $self->{secpass_base_ct}/
            &max(0.001, $self->{sp_end_time}->[1] - $self->{fp_end_time}->[1]);
        print $fh "\n$second_pass_label analysis of tRNAs ended: ",`date`,"\n";
        if ($opts->tscan_mode() || $opts->eufind_mode())
        {        
            print $fh "Summary\n--------\n";
        }
    }                                
    my $total_time = ($self->{sp_end_time}->[0] - $self->{fp_start_time}->[0]) + 
        ($self->{sp_end_time}->[1] - $self->{fp_start_time}->[1]);
    printf $fh "Overall scan speed: %.1f bp/sec\n",
        &max($self->{first_pass_base_ct} * 2, $self->{secpass_base_ct}) / &max(0.001, $total_time);

    $self->output_summary($opts, $gc, $r_tab_results);                
}

sub output_summary
{
    my $self = shift;
    my $opts = shift;
    my $gc = shift;
    my $r_tab_results = shift;
    my $fh = $self->{FILE_H};
    
    my ($trna_ct, $selcys_ct, $stop_sup_ct, $undet_ct, $pseudo_ct, $mismatch_ct,
           $total, $intron_ct, $line);
    my (%iso_AR, %ac_AR, %intron_ac_AR, %iso_cm_AR);
    my ($iso, $ac, $acset, $iso_count, $iso_cm_count, $istart, $aa, $iso_cm);
    my @columns = ();
    my $note_col = -1;
    my $iso_cm_col = -1;
           
    $trna_ct   = 0;
    $selcys_ct = 0;
    $pseudo_ct = 0;
    $mismatch_ct = 0;
    $undet_ct  = 0;
    $intron_ct = 0;
    $stop_sup_ct = 0;
    $total = 0;
    
    $line = shift(@$r_tab_results);

    while ($line ne '')
    {
        if ($line =~ /^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+([0-9\,]+)\s+([0-9\,]+)\s+(\S+)/)
        {
            @columns = split(/\t/, $line, -1);
            
            $iso     = $columns[4];
            $ac      = $columns[5];
            $istart  = &trim($columns[6]);
            
            if ($note_col == -1)
            {
                $note_col = scalar(@columns) - 1;
            }
            if ($iso_cm_col == -1)
            {
                if ((($opts->euk_mode() or $opts->bact_mode() or $opts->arch_mode()) and !$opts->no_isotype() and $opts->detail()) or $opts->metagenome_mode())
                {
					$iso_cm_col = $note_col - 2;
				}
                if ($opts->euk_mode() and $opts->mito_model() ne "")
                {
                    $iso_cm_col--;
                }
			}
			if ($iso_cm_col > -1)
            {
				$iso_cm = $columns[$iso_cm_col];
			}
			
            if ($note_col > -1 and index($columns[$note_col], "pseudo") > -1)
            {
                $pseudo_ct++;
                $iso_AR{"Pseudo"}++;
            }
            elsif ($note_col > -1 and index($columns[$note_col], "IPD") > -1)
            {
                $mismatch_ct++;
            }
            elsif ($iso eq $gc->undef_isotype())
            {
                $undet_ct++;
            }            
            elsif ($iso =~ /SeC/)
            {
                $selcys_ct++;
            }
            elsif ($iso eq "Sup")
            {
                $stop_sup_ct++;
            }           
            else
            {
                $trna_ct++;
            }
            
            if ($iso =~ /SeC/)
            {
                $iso_AR{"SelCys"}++;
                $ac_AR{"SelCys"}->{$ac}++;
            }
            elsif ($iso eq "Sup")
            {
                $iso_AR{"Supres"}++;
                $ac_AR{"Supres"}->{$ac}++;
            }
            else
            {
                $iso_AR{$iso}++;
                $ac_AR{$iso}->{$ac}++;                
            }
            if (defined $iso_cm)
            {
				$iso_cm_AR{$iso_cm}++;
			}
			
            if ($istart ne "0")
            {
                my @introns = split(/\,/, $istart);
                $intron_ct += scalar(@introns);
                my $tag = $iso;
                if ($iso =~ /SeC/)
                {
                    $tag = "SelCys";
                }
                elsif ($iso eq "Sup")
                {
                    $tag = "Supres";
                }
                $intron_ac_AR{$tag}->{$ac} += scalar(@introns);
            }            
        }
        $line = shift(@$r_tab_results);        
    }
    
    $total = $trna_ct + $selcys_ct + $pseudo_ct + $mismatch_ct + $undet_ct + $stop_sup_ct;
    
    
    print $fh "\n",
    "tRNAs decoding Standard 20 AA:              $trna_ct\n",
    "Selenocysteine tRNAs (TCA):                 $selcys_ct\n",
    "Possible suppressor tRNAs (CTA,TTA,TCA):    $stop_sup_ct\n",
    "tRNAs with undetermined/unknown isotypes:   $undet_ct\n";
    
    if ((($opts->euk_mode() or $opts->bact_mode() or $opts->arch_mode()) and !$opts->no_isotype() and $opts->detail()) or $opts->metagenome_mode())
    {
        print $fh "tRNAs with mismatch isotypes:               $mismatch_ct\n";
    }
    
    print $fh "Predicted pseudogenes:                      $pseudo_ct\n",
    "                                            -------\n",
    "Total tRNAs:                                $total\n\n",
    
    "tRNAs with introns:     \t$intron_ct\n\n";

    foreach $aa (@{$gc->isotypes()})
    {
        foreach $acset ($gc->ac_list()->{$aa})
        {
            foreach $ac (@$acset)
            {
                if (defined($intron_ac_AR{$aa}->{$ac}))
                {
                    print $fh "| $aa-$ac: $intron_ac_AR{$aa}->{$ac} "; 
                }
            }
        }
    }
    if (scalar(keys %intron_ac_AR) > 0)
    {	
        print $fh "|\n\n";
    }
    print $fh "Isotype / Anticodon Counts:\n\n";
    
    foreach $aa (@{$gc->isotypes()})
    {
        my $label = $aa;
        $iso_count = $iso_AR{$aa} + 0;
        $iso_cm_count = $iso_cm_AR{$aa} + 0;
        if ($aa eq "Met")
        {
            if (defined $iso_AR{iMet})
            {
				$label = "Met/iMet";
                $iso_count += $iso_AR{iMet};
			}
			elsif (defined $iso_AR{fMet})
            {
				$label = "Met/fMet";
                $iso_count += $iso_AR{fMet};
			}
            if (defined $iso_cm_AR{iMet})
            {
				$iso_cm_count += $iso_cm_AR{iMet};
			}
            elsif (defined $iso_cm_AR{fMet})
            {
				$iso_cm_count += $iso_cm_AR{fMet};
			}
		}
        elsif ($aa eq "Ile")
        {
            if (defined $iso_AR{Ile2})
            {
                $iso_count += $iso_AR{Ile2};
			}
            if (defined $iso_cm_AR{Ile2})
            {
				$iso_cm_count += $iso_cm_AR{Ile2};
			}
		}
		
        if ((($opts->euk_mode() or $opts->bact_mode() or $opts->arch_mode()) and !$opts->no_isotype() and $opts->detail()) or $opts->metagenome_mode())
        {
            printf $fh ("%-8s: %d (%d)\t", $label, $iso_count, $iso_cm_count);
        }
        else
        {
            printf $fh ("%-8s: %d\t", $label, $iso_count);
        }
        
        foreach $acset ($gc->ac_list()->{$aa})
        {
            foreach $ac (@$acset)
            {
                if ($ac eq "&nbsp")
                {
                    print $fh "             ";
                }
                else
                {
                    my $count = $ac_AR{$aa}->{$ac};
                    if ($aa eq "Met")
                    {
						if (defined $ac_AR{iMet}->{$ac})
                        {
                            $count += $ac_AR{iMet}->{$ac};
                        }
                        elsif (defined $ac_AR{fMet}->{$ac})
                        {
                            $count += $ac_AR{fMet}->{$ac};
                        }
					}
					elsif ($aa eq "Ile")
                    {
						if (defined $ac_AR{Ile2}->{$ac})
                        {
                            $count += $ac_AR{Ile2}->{$ac};
                        }
                    }
                    printf $fh ("%5s: %-6s",$ac,$count);
                }
            }
        }
        
        print $fh "\n";
        
    }
    print $fh "\n";
}

1;
