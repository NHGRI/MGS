#!/usr/bin/env perl
# $Id: prepare_data_stouffer.plx 6244 2018-04-22 14:00:37Z narisu $

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use GTB::File qw(Open);
our %Opt;
our $EMPTY = q{};

=head1 NAME

prepare_data_stouffer.plx - a script to prepare input for Stouffer test

=head1 SYNOPSIS

This script prepares data for Stouffer meta analysis for ATAC-seq allelic imbalance. 

  prepare_data_stouffer.plx [options] file1, file2 ...

For complete documentation, run C<prepare_data_stouffer.plx -man>

=head1 DESCRIPTION

This script currently takes bed files of ATAC-seq counts, mostly generated from WASP and some post processing. 
One file per sample. Expected columns are:
    1   #chr
    2   start
    3   stop
    4   variant_ID
    5   ref_base
    6   alt_base
    7   ref_count
    8   alt_count
    9   total_count
    10  binom_p_val
    11  binom_q_val ??
    12  ref.effectSize
    13  beta_binom_p_val T2D signal
    14  beta_binom_q_val
    15  beta.binom.prob
    16  beta.binom.rho

=cut

#------------
# Begin MAIN 
#------------

process_commandline();
my $ofh = Open($Opt{output}, "w");

my %snps;
while (<>) {
    next if /^#/;
    my @data=split;
    # ref, alt, total, beta binom_p, q
    my @val = @data[6..10];
    push @{$snps{join("-",@data[0..5])}}, \@val; 
}

foreach my $site(sort keys %snps) {
    my $values = join("\t", split(/-/, $site))."\t";
    my @ra = @{$snps{$site}};
    my $total_reads = 0;
    foreach my $ra_pvalue(@ra) {
        my @col = @$ra_pvalue;
        $values = $values.join(":", @col).";";
        $total_reads += $col[2]; 
    }
    $values =~ s/;$//;
    print $values."\n" if @ra >= $Opt{count} && $total_reads >= $Opt{read};
}

#------------
# End MAIN
#------------

sub process_commandline {
    %Opt = (output  => '-',
            );
    GetOptions(\%Opt, qw(output=s
                count=i
                read=i
                help+ manual version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) {
        die "prepare_data_stouffer.plx, ", q$Revision: 6244 $, "\n";
    }
    # This script requires NO list of non-option arguments
    if (@ARGV) {
#        pod2usage("No arguments expected, only named options");
    }
    $Opt{count} ||= 3;
    # Required options
    my @miss = map { "-$_" } grep { !$Opt{$_} } qw();
    if (@miss) {
        pod2usage("Please supply required @miss on command line");
    }
    # Add criteria to test other options below
    # Whenever possible, Globals should be initialized in main code
}

=head1 OPTIONS

=over 4

=item B<--count> 

Minimum number of samples with hets for Stouffer test. Default to 3.

=item B<--reads> 

Minimum number of total reads for B<--count> samples with >= this many total reads. 

=item B<--output> FILE

Destination file.  Otherwise, output is written to STDOUT.

=item B<--help|--manual>

Display documentation.  C<--help> gives a brief synopsis, C<-h -h> shows
all options, C<--manual> provides complete documentation.

=back

=head1 AUTHOR

 Narisu Narisu - narisu@mail.nih.gov

=head1 LEGAL

This software/database is "United States Government Work" under the terms of
the United States Copyright Act.  It was written as part of the authors'
official duties for the United States Government and thus cannot be
copyrighted.  This software/database is freely available to the public for
use without a copyright notice.  Restrictions cannot be placed on its present
or future use. 

Although all reasonable efforts have been taken to ensure the accuracy and
reliability of the software and data, the National Human Genome Research
Institute (NHGRI) and the U.S. Government does not and cannot warrant the
performance or results that may be obtained by using this software or data.
NHGRI and the U.S.  Government disclaims all warranties as to performance,
merchantability or fitness for any particular purpose. 

In any work or product derived from this material, proper attribution of the
authors as the source of the software or data should be made, using "NHGRI
FUSION Research Group" as the citation. 

=cut
