#!/usr/bin/perl
# test and compare different versions of ijkMCtable_invert.
# created by R. Wenger March 18, 2024

use strict;

my $temp0_xit = "temp0.xit";
my $temp1_xit = "temp1.xit";

my @proglist = @ARGV;
my @flag;
my $fastflag = 0;
my $diff_options = "";

while (scalar(@proglist) > 0 &&
       $proglist[0] =~ /-.*/) {

  my $newflag = shift(@proglist);

  if ($newflag eq "-fast") {
    $fastflag = 1;
    next;
  }
  elsif ($newflag eq "-ignore_blanks") {
    $diff_options = "-b";
  }

  push(@flag, $newflag);
}

if ((scalar(@proglist) < 1) || (scalar(@proglist) > 2)) {
  print "testinvert.pl requires one or two input files.\n";
  usage_error();
}

my %testdata;

$testdata{c2}{poly} = "cube";
$testdata{c2}{dimension} = 2;

$testdata{c3}{poly} = "cube";
$testdata{c3}{dimension} = 3;

$testdata{s2}{poly} = "simplex";
$testdata{s2}{dimension} = 2;

$testdata{s3}{poly} = "simplex";
$testdata{s3}{dimension} = 3;

$testdata{s4}{poly} = "simplex";
$testdata{s4}{dimension} = 4;

$testdata{s5}{poly} = "simplex";
$testdata{s5}{dimension} = 5;

$testdata{pyr2}{poly} = "pyramid";
$testdata{pyr2}{dimension} = 2;

$testdata{pyr3}{poly} = "pyramid";
$testdata{pyr3}{dimension} = 3;

$testdata{pyr4}{poly} = "pyramid";
$testdata{pyr4}{dimension} = 4;


my $prog0 = shift(@proglist);

if (scalar(@proglist) < 1) {

  foreach my $tdata (keys %testdata) {

    my $poly = $testdata{$tdata}{poly};
    my $dimension = $testdata{$tdata}{dimension};

    my $isotablefile = "iso." . "$poly." . "$dimension" . "D.xit";
    test_prog($prog0, "", $poly, $dimension, $isotablefile);

# check interval volumes
    if (!$fastflag || ($tdata ne "c3" && $tdata ne "pyr4")) {
      my $ivoltablefile = "ivol." . "$poly." . "$dimension" . "D.xit";
      test_prog($prog0, "", $poly, $dimension, $ivoltablefile);
    }

# check nep
    $isotablefile = "iso.nep." . "$poly." . "$dimension" . "D.xit";
    test_prog($prog0, "", $poly, $dimension, $isotablefile);
  }
}
else {
  my $prog1 = shift(@proglist);

  foreach my $tdata (keys %testdata) {

    my $poly = $testdata{$tdata}{poly};
    my $dimension = $testdata{$tdata}{dimension};

    my $isotablefile = "iso." . "$poly." . "$dimension" . "D.xit";
    compare_prog($prog0, $prog1, "", $poly, $dimension, $isotablefile);

# check interval volumes
    if (!$fastflag || ($tdata ne "c3" && $tdata ne "pyr4")) {
      my $ivoltablefile = "ivol." . "$poly." . "$dimension" . "D.xit";
      compare_prog($prog0, $prog1, "", $poly, $dimension, $ivoltablefile);
    }

# check nep
    $isotablefile = "iso.nep." . "$poly." . "$dimension" . "D.xit";
    compare_prog($prog0, $prog1, "", $poly, $dimension, $isotablefile);
  }
}

if (-e "$temp0_xit") { system "rm $temp0_xit"; };
if (-e "$temp1_xit") { system "rm $temp1_xit"; };


# ***********************************

sub test_prog() {

  my $num_param = 5;
  scalar(@_) == $num_param ||
    die "Error in sub compare_prog. Requires $num_param parameters.\n";

  my $prog0 = $_[0];
  my $param_list = $_[1];
  my $poly = $_[2];
  my $dimension = $_[3];
  my $isotablefile = $_[4];

  run_ijkMCtable_invert
      ($prog0, $param_list, $poly, $dimension, $isotablefile, $temp0_xit);

  # Invert table again.
  run_ijkMCtable_invert
      ($prog0, $param_list, $poly, $dimension, $temp0_xit, $temp1_xit);

  system "diff $diff_options --brief $isotablefile $temp1_xit";

}


# ***********************************

sub compare_prog() {

  my $num_param = 6;
  scalar(@_) == $num_param ||
    die "Error in sub compare_prog. Requires $num_param parameters.\n";

  my $prog0 = $_[0];
  my $prog1 = $_[1];
  my $param_list = $_[2];
  my $poly = $_[3];
  my $dimension = $_[4];
  my $isotablefile = $_[5];

  run_ijkMCtable_invert
      ($prog0, $param_list, $poly, $dimension, $isotablefile, $temp0_xit);
  run_ijkMCtable_invert
      ($prog1, $param_list, $poly, $dimension, $isotablefile, $temp1_xit);
  system "diff $diff_options --brief $temp0_xit $temp1_xit";
}

# ***********************************

sub usage_error() {
 print "Usage: testinvert.pl {prog1} [{prog2}]\n";
 exit(10);
}

sub run_ijkMCtable_invert() {

  my $num_param = 6;
  scalar(@_) == $num_param ||
    die "Error in sub run_ijkMCtable_inverte. Requires $num_param parameters.\n";

  my $ijkMCtable_invert = $_[0];
  my $param_list = $_[1];
  my $poly = $_[2];
  my $dimension = $_[3];
  my $isotablefile = $_[4];
  my $outfile = $_[5];

  my $command_line =
      "./$ijkMCtable_invert $param_list $isotablefile $outfile";

  print "$command_line\n";
  system "$command_line";
}


