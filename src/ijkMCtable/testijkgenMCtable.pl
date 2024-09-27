#!/usr/bin/perl
# test and compare different versions of ijkgenMCtable
# created by R. Wenger Jan. 13, 2008

use strict;

my $temp_xit = "temp.xit";

my @proglist = @ARGV;
my @flag;
my $orientflag = 0;
my $fastflag = 0;
my $check_all_pairs = 0;
my $diff_options = "";

while (scalar(@proglist) > 0 &&
       $proglist[0] =~ /-.*/) {

  my $newflag = shift(@proglist);

  if ($newflag eq "-orient") {
    $orientflag = 1;
    next;
  }
  elsif ($newflag eq "-fast") {
    $fastflag = 1;
    next;
  }
  elsif ($newflag eq "-check_all_pairs") {
    $check_all_pairs = 1;
    next;
  }
  elsif ($newflag eq "-ignore_blanks") {
    $diff_options = "-b";
  }

  push(@flag, $newflag);
}

if (scalar(@proglist) != 2) {
  print "testijkgentable.pl requires two input files.\n";
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
my $prog1 = shift(@proglist);

foreach my $tdata (keys %testdata) {

  my $poly = $testdata{$tdata}{poly};
  my $dimension = $testdata{$tdata}{dimension};
  my $fname_prefix = "";
  
  $fname_prefix = "iso." . "$poly";
  run_on_all_orient_and_sep
      ($prog0, $prog1, "", $poly, $dimension, "$fname_prefix");

  # check nep
  $fname_prefix = "iso.nep." . "$poly";
  run_on_all_orient_and_sep
      ($prog0, $prog1, "-nep", $poly, $dimension, $fname_prefix);

# check interval volumes
  if (!$fastflag || ($tdata ne "c3" && $tdata ne "pyr4")) {
    $fname_prefix = "ivol." . "$poly";
    run_on_all_orient_and_sep
      ($prog0, $prog1, "-ivol", $poly, $dimension, $fname_prefix);
  }

# check -edge_groups, -chull and -no_sep_opp with dimension 3 cube.
  my $dimension = $testdata{$tdata}{dimension};

  if (($dimension == 3) && ("$poly" eq "cube")) {
    $fname_prefix = "iso.cube.edgeGroups";
    run_on_all_orient_and_sep
      ($prog0, $prog1, "-edge_groups", $poly, $dimension, $fname_prefix);
    $fname_prefix = "iso.cube.cHull";
    run_on_all_orient_and_sep
      ($prog0, $prog1, "-chull -no_sep_opp", $poly, $dimension, $fname_prefix);
    run_on_all_orient_and_sep
      ($prog0, $prog1, "-chull", $poly, $dimension, $fname_prefix);
  }
  
}
if (-e "$temp_xit") { system "rm $temp_xit"; };

# ***********************************

sub compare_isotables {

  my $num_param = 6;
  scalar(@_) == $num_param ||
    die "Error in sub compare_isotables. Requires $num_param parameters.\n";

  my $prog0 = $_[0];
  my $prog1 = $_[1];
  my $param_list = $_[2];
  my $poly = $_[3];
  my $dimension = $_[4];
  my $isotablefile = $_[5];

  run_ijkgentable($prog0, $param_list, $poly, $dimension, $isotablefile);
  run_ijkgentable($prog1, $param_list, $poly, $dimension, $temp_xit);
  system "diff $diff_options --brief $isotablefile $temp_xit";
}

# ***********************************

sub run_on_all_orient_and_sep_old {

  my $num_param = 6;
  scalar(@_) == $num_param ||
    die "Error in sub run_on_all_orient_and_sep. Requires $num_param parameters.\n";

  my $prog0 = $_[0];
  my $prog1 = $_[1];
  my $param_list = $_[2];
  my $poly = $_[3];
  my $dimension = $_[4];
  my $fname_prefix = $_[5];

  my $isotablefile;
  my $paramII;
  my $labels = "";
  
  if ("$param_list" eq "-edge_groups") {
    $labels = ".edgeGroups";
  }
  elsif ("$param_list" eq "-chull") {
    $labels = ".cHull";
  }

  my $fname_prefixII = $fname_prefix . ".$poly" . "$labels";
  
  $isotablefile = $fname_prefixII . ".sepNeg.negO." . "$dimension" . "D.xit";
  $paramII = $param_list . " -sep_neg -neg_orient";
  compare_isotables($prog0, $prog1, $paramII, $poly, $dimension, $isotablefile);

  $isotablefile = $fname_prefixII . ".sepNeg.posO." . "$dimension" . "D.xit";
  $paramII = $param_list . " -sep_neg -pos_orient";
  compare_isotables($prog0, $prog1, $paramII, $poly, $dimension, $isotablefile);

  $isotablefile = $fname_prefixII . ".sepPos.negO." . "$dimension" . "D.xit";
  $paramII = $param_list . " -sep_pos -neg_orient";
  compare_isotables($prog0, $prog1, $paramII, $poly, $dimension, $isotablefile);

  $isotablefile = $fname_prefixII . ".sepPos.posO." . "$dimension" . "D.xit";
  $paramII = $param_list . " -sep_pos -pos_orient";
  compare_isotables($prog0, $prog1, $paramII, $poly, $dimension, $isotablefile);
}

# ***********************************

sub run_on_all_orient_and_sep {

  my $num_param = 6;
  scalar(@_) == $num_param ||
    die "Error in sub run_on_all_orient_and_sep. Requires $num_param parameters.\n";

  my $prog0 = $_[0];
  my $prog1 = $_[1];
  my $param_list = $_[2];
  my $poly = $_[3];
  my $dimension = $_[4];
  my $fname_prefix = $_[5];

  my $isotablefile;
  my $paramII;
  my $labels = "";

  $isotablefile = $fname_prefix . ".sepNeg.negO." . "$dimension" . "D.xit";
  $paramII = $param_list . " -sep_neg -neg_orient";
  compare_isotables($prog0, $prog1, $paramII, $poly, $dimension, $isotablefile);

  $isotablefile = $fname_prefix . ".sepNeg.posO." . "$dimension" . "D.xit";
  $paramII = $param_list . " -sep_neg -pos_orient";
  compare_isotables($prog0, $prog1, $paramII, $poly, $dimension, $isotablefile);

  $isotablefile = $fname_prefix . ".sepPos.negO." . "$dimension" . "D.xit";
  $paramII = $param_list . " -sep_pos -neg_orient";
  compare_isotables($prog0, $prog1, $paramII, $poly, $dimension, $isotablefile);

  $isotablefile = $fname_prefix . ".sepPos.posO." . "$dimension" . "D.xit";
  $paramII = $param_list . " -sep_pos -pos_orient";
  compare_isotables($prog0, $prog1, $paramII, $poly, $dimension, $isotablefile);
}

# ***********************************

sub usage_error() {
 print "Usage: testijkgentable.pl {prog1} {prog2} [{prog3} ...]\n";
 exit(10);
}

sub run_ijkgentable {

  my $num_param = 5;
  scalar(@_) == $num_param ||
    die "Error in sub run_ijkgentable. Requires $num_param parameters.\n";

  my $ijkgentable = $_[0];
  my $param_list = $_[1];
  my $poly = $_[2];
  my $dimension = $_[3];
  my $isotablefile = $_[4];

  my $command_line;
  if ($check_all_pairs) {
    $command_line = 
        "./$ijkgentable $param_list -check_all_pairs -o $isotablefile -q -d $dimension -poly $poly";
  }
  else {
    $command_line = 
      "./$ijkgentable $param_list -o $isotablefile -q -d $dimension -poly $poly";
  }


  print "$command_line\n";
  system("$command_line") == 0 ||
    die "Program $ijkgentable abnormally terminated.\n";


  if ($orientflag) {
    $command_line = "ijkorient $isotablefile";
    system "$command_line";
  }

}
