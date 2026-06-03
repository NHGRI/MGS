#!/usr/bin/env perl
# $Id: bam_stat.plx 6244 2018-04-22 14:00:37Z narisu $

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use GTB::File qw(Open);
our %Opt;
our $EMPTY = q{};

=head1 NAME

bam_stat.plx - this script calculates reads statistics of a bam file

=head1 SYNOPSIS

This script extracts read statistics of a bam file. Following fields will be 
generated:
    total number of reads                   
    total number of primary reads           
    % of reads mapped as primary 
    total number of primary reads in a region 
    % primary reads mapped in the region            
    average depth in the region

  bam_stat.plx [options] -bam bam_file -output stat.out

For complete documentation, run C<bam_stat.plx -man>

=head1 DESCRIPTION

This is a generic descriptive statistics of a bam file. 

=cut

#------------
# Begin MAIN 
#------------

process_commandline();
my $ofh = Open($Opt{output}, "w");

# total number of reads, exclude secondary alignment and
# supplementary alignment. This number matches records in
# fastq files.
my $cmd = "samtools view -c -F 2304 $Opt{bam}";
my $total = `$cmd` || die "Can not execute '$cmd'";
chomp($total);
warn "Sample: $Opt{sample}\n";
warn "Total reads: $total\n";

# get read length
$cmd = "samtools view $Opt{bam} | head -100";
my $read_length = 0;
open(READ, "$cmd |") || die "Can not execute '$cmd'";
while (<READ>) {
    my @d = split(/\t/);
    $read_length = length($d[9]) if $read_length < length($d[9]) ;
}
if ($read_length) {
    print "Read length is $read_length (bp)\n";
} else {
    die "No read length found for the bam file\n";
}

# total number of reads with primary alignments 
$cmd = "samtools view -c -F 2308 $Opt{bam}";
my $primary = `$cmd` || die "Can not execute '$cmd'";
chomp($primary);
warn "Primary aligned reads: $primary\n";
my $primary_per = sprintf("%.02f", $primary*100/$total);

my @header = (
        "Sample",
        "# of reads",
        "Length of reads",
        "# of reads mapped as primary alignment",
        "% of reads mapped as primary alignment"
        );
# total number of reads with primary alignments in a region
my ($in_region, $depth, $length_of_region);
if ($Opt{region}) {
    $cmd = "samtools view -c -F 2308 -L $Opt{region} $Opt{bam} ";
    $in_region = `$cmd` || die "Can not execute '$cmd'";
    chomp($in_region);
    warn "In region aligned reads: $in_region\n";
    my $cmd_length = "bed_union $Opt{region} |bed_stats - ";
    my @output = `$cmd_length 2>&1` || die "Can not execute '$cmd_length'";
    foreach my $o(@output) {
        if($o =~ /Total length of features: (\d+)/) {
            $length_of_region = $1;
        }
    }
    warn "Length of region: $length_of_region\n";
    $depth = sprintf("%.02f", $in_region*$read_length/$length_of_region);
    my $in_region_per = sprintf("%.02f", $in_region*100/$primary);

    # coverage of the primary
    $cmd = "samtools view -b -F 2308 $Opt{bam} | "
        ."bedtools intersect -abam - -b $Opt{region} -bed |"
        ."sort-bed  --max-mem 10G --tmpdir ./tmp - | bedops --merge - |bed_stats - ";
    @output = `$cmd 2>&1` || die "Can not execute '$cmd'";
    my $covered_region;
    foreach my $o(@output) {
        if($o =~ /Total length of features: (\d+)/) {
            $covered_region = $1;
        }
    }
    warn "Covered regions: $covered_region\n";
    my $coverage_per = sprintf("%.02f", $covered_region*100/$length_of_region);

    print $ofh join("\t", @header, 
        "# of reads mapped to consensus coding region",
        "% of reads mapped to consensus coding region",
        "Depth of sequence in the consensus coding region",
        "Size of consensus coding region (bp)",
        "Size of the consensus coding region covered at least by 1x",
        "% of consensus coding region covered by at least 1x")."\n";
    print $ofh join("\t", $Opt{sample}, $total, $read_length, $primary, $primary_per, 
        $in_region, $in_region_per, $depth,
        $length_of_region, $covered_region, $coverage_per)."\n";
}
else {
    print $ofh join("\t", @header)."\n";
    print $ofh join("\t", $Opt{sample}, $total, $read_length, $primary, $primary_per)."\n";
}
#------------
# End MAIN
#------------

sub process_commandline {
    %Opt = (output  => '-',
            );
    GetOptions(\%Opt, qw(output=s
                region=s
                bam=s
                sample=s
                help+ manual version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) {
        die "bam_stat.plx, ", q$Revision: 6244 $, "\n";
    }
    # This script requires NO list of non-option arguments
    if (@ARGV) {
        pod2usage("No arguments expected, only named options");
    }
    # Required options
    my @miss = map { "-$_" } grep { !$Opt{$_} } qw();
    if (@miss) {
        pod2usage("Please supply required @miss on command line");
    }
}

=head1 OPTIONS

=over 4

=item B<--bam> FILE

Input bam file.  

=item B<--region> FILE

A bed file of regions where reads overlap. 

=item B<--sample> STRING

Sample name.

=item B<--output> FILE

Destination read statistics summary file. Otherwise, output is written to STDOUT.

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
