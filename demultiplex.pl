# This script was written by Wei Wang

#! /bin/env perl

use strict;
use Text::LevenshteinXS;
use PerlIO::gzip;

open GZ, "<::gzip", "$ARGV[0]" or die $!;
open BC, "$ARGV[1]" or die $!;
open SEQ, ">$ARGV[2]" or die $!;
my %files;
my %tallyup;
my %samples;
my %duplicate;

$samples{'discard_nomat'}{'id'} = 'discard_nomat';
$samples{'discard_multi'}{'id'} = 'discard_multi';
$samples{'discard_nomat'}{'nm'} = 0;
$samples{'discard_multi'}{'nm'} = 0;

while (<BC>) {
    chomp;
    my ($id, $code) = split(/\t/);
    open ($files{$code}, ">::gzip", "$id.fastq.gz") or die $!;
    $samples{$code}{'id'} = $id;
    $samples{$code}{'nm'} = 0;
}
close(BC);
my $mis_file = $ARGV[0];
$mis_file =~ s/\.fastq/_mis\.fastq/;
open DIS, ">::gzip", "$mis_file";

while (my $line = <GZ>) {
    my $output = $line;
    my $ns = 0;
    $line = <GZ>;
    my $tmp_code = substr $line, 0, 4;
    my @numbers = $tmp_code =~ m/(N)/g;
    my $ns = scalar(@numbers);
    my $trim = substr $line, 4;
    $output .= $trim;
    chomp($trim);
    $line = <GZ>;
    $output .= $line;
    $line = <GZ>;
    $line = substr $line, 4;
    $output .= $line;
    my %mis_count={};
    my $mis0 = 0;
    my $mis1 = 0;
    foreach my $code (keys %files) {
        my $dist = distance($tmp_code, $code) - $ns;
        if ($dist == 0) {
            $mis0++;
            $mis_count{'0'} = $code;
        }
        elsif ($dist == 1) {
            $mis1++;
            $mis_count{'1'} = $code;
        }
    }
    if ($mis0 == 1 && $mis1 == 0) {
        $tallyup{'perfect_m'}{$tmp_code}++;
        $samples{$mis_count{'0'}}{'nm'}++;
        print {$files{$mis_count{'0'}}} $output;
        if (not exists $duplicate{$mis_count{'0'}}{$trim}) {
            $duplicate{$mis_count{'0'}}{$trim} = 1;
        }
        else {
            $duplicate{$mis_count{'0'}}{$trim}++;
        }
    }
    elsif ($mis0 == 0 && $mis1 == 1) {
        $tallyup{'1_mismat'}{$tmp_code}++;
        $samples{$mis_count{'1'}}{'nm'}++;
        print {$files{$mis_count{'1'}}} $output;
        if (not exists $duplicate{$mis_count{'1'}}{$trim}) {
            $duplicate{$mis_count{'1'}}{$trim} = 1;
        }
        else {
            $duplicate{$mis_count{'1'}}{$trim}++;
        }
    }
    elsif ($mis0 >=1 && $mis1 >= 1) {
        print DIS $output;
        $tallyup{'multi_mat'}{$tmp_code}++;
        $samples{'discard_multi'}{'nm'}++;
    }
    else {
        print DIS $output;
        $tallyup{'no_match'}{$tmp_code}++;
        $samples{'discard_nomat'}{'nm'}++;
    }
}

print "Barcode Frequency:\n";
foreach my $sit (keys %tallyup) {
    foreach my $code (keys %{$tallyup{$sit}}) {
        print $sit,"\t",$code,"\t",$tallyup{$sit}{$code},"\n";
    }
}

print "small RNA frequency:\n";
foreach my $code (keys %duplicate) {
    print SEQ $samples{$code}{'id'},"\n";
    foreach my $seq (keys %{$duplicate{$code}}) {
        print SEQ $seq,"\t",$duplicate{$code}{$seq},"\n";
    }
}

print "\nSample Frequency:\nSample_id\tReads_num\tUniq_Reads\n";
foreach my $tmp (keys %samples) {
    print $samples{$tmp}{'id'},"\t",$samples{$tmp}{'nm'},"\t";
    print scalar keys %{$duplicate{$tmp}};
    print "\n";
}
