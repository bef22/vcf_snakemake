#!/usr/bin/perl

# setPassFilter.pl

use warnings;
use strict;
use Getopt::Long;
use IO::Zlib;


my $version = '1.0.0';
my $medianDP = undef;                  # genome-wide median
my $percentDP = undef;                 
my $lowDP = undef;                     # lower DP threshold medianDP - percentDP
my $highDP = undef;                    # upper DP threshold medianDP + precentDP
my $minQ = undef;                      # min QUAL
my $minGQ = undef;                     # min mean GQ
my $maxMissing = undef;                # % max missing per site
my $vcf = undef;

die(qq/
setPassFilter.pl version $version

bcftools view <bcf file> | setPassFilter.pl [options] - | bcftools view -O b -o <bcf file>

--medianDP  	genome-wide median DP, format '<INT>'
--percentDP  	% of the medianDP to define lowDP and HighDP bounds of the medianDP to annotate FILTER field, format '<INT>'
--minQ		minimum QUAL value to annotate FILTER field, format '<INT>'
--minGQ   	minimum mean GQ value to annotate FILTER field, format '<INT>'
--maxMissing	maximum % missing to annotate FILTER field, format '<INT>'

Notes: 
* setting any of above parameters to 0 switches them off

* this script removes\/overwrites already existing filter settings

\n/) if (-t STDIN && !@ARGV);

GetOptions('medianDP=i' => \$medianDP, 'percentDP=i' => \$percentDP, 'minQ=i' => \$minQ, 'minGQ=i' => \$minGQ, 'maxMissing=i' => \$maxMissing);

# this is for printing the command to the VCF header
my %userargs = ('--medianDP' => $medianDP, '--percentDP' => $percentDP, '--minQ' => $minQ, '--minGQ' => $minGQ, '--maxMissing' => $maxMissing);



### check options are valid
if ($medianDP !~ m/\d+/) {
	die("ERROR: --medianDP requires an integer - use genome wide median\n");
}
if ($percentDP !~ m/\d+/) {
	die("ERROR: --percentDP requires an integer threshold\n");
}
if ($minQ !~ m/\d+/) {
	die("ERROR: --minQ requires an integer threshold\n");
}
if ($minGQ !~ m/\d+/) {
	die("ERROR: --minGQ requires an integer threshold\n");
}
if ($maxMissing !~ m/\d+/) {
	die("ERROR: --maxMissing requires an integer threshold\n");
}

# low and high DP thresholds
$lowDP = $medianDP - ($medianDP * ($percentDP / 100));
$highDP = $medianDP + ($medianDP * ($percentDP / 100));


### print pretty VCF header - took this from Tylers script

my $datestr = localtime();
my $command = "##insertFilterCommand=<ID=setPassFilter.pl,Version=${version},Date=\"$datestr\",CommandLineOptions=\"";
foreach my $arg (keys %userargs) {
	$command .= "$arg $userargs{$arg} " if (defined $userargs{$arg});
}
$command .= "$vcf" if ($vcf);
$command =~ s/\s$//;
$command .= "\">\n";

my $lowdp_string = "##FILTER=<ID=LowDP,Description=\"Site DP less than genome-wide mean site DP -$percentDP%\">\n";
my $highdp_string = "##FILTER=<ID=HighDP,Description=\"Site DP greater than genome-wide median site DP +$percentDP%\">\n";
my $minQ_string = "##FILTER=<ID=lowQ,Description=\"QUAL score below $minQ\">\n";
my $minGQ_string = "##FILTER=<ID=lowGQ,Description=\"Site mean GQ score below $minGQ\">\n";
my $miss_string = "##FILTER=<ID=highMiss,Description=\"Site with missing individuals >$maxMissing%\">\n";

my @headorder = ('fileformat', 'reference', 'contig', 'INFO', 'FILTER', 'FORMAT', 'ALT', 'other'); # header order

my %header = (fileformat => undef, ALT => undef, FILTER => undef, INFO => undef, FORMAT => undef, reference => undef, contig => undef, other => undef);

### read in VCF
$vcf = pop @ARGV;
my $vcffh;
my $nSamples = 0;
my $passCnt = 0;
my $failCnt = 0;


open($vcffh, "<", $vcf) or die("Unable to open VCF $vcf: $!\n");

my $hline;
my $filter;
while (($hline = <$vcffh>) =~ /^##/) {
	if ($hline =~ /^##([^=]+)/i) {
		my $annotation = $1;
		# remove existing filters
		if ($hline =~ /^##FILTER=<ID=([^,]+)/) { 
			next;
		} elsif ($annotation =~ /^fileformat|ALT|FILTER|INFO|FORMAT|reference|contig/) {
			push @{$header{$annotation}}, $hline;
		} else {
			push @{$header{other}}, $hline;
		}
	}
}

push @{$header{FILTER}}, $lowdp_string if ($medianDP && $percentDP);
push @{$header{FILTER}}, $highdp_string if ($medianDP && $percentDP);
push @{$header{FILTER}}, $minQ_string if ($minQ);
push @{$header{FILTER}}, $minGQ_string if ($minGQ);
push @{$header{FILTER}}, $miss_string if ($maxMissing);

push @{$header{other}}, $command;

foreach my $tag (@headorder) {
	foreach (@{$header{$tag}}) {
		print STDOUT $_;
	}
}

print STDOUT $hline;

#exit;


################################

while (<$vcffh>) {
	my @filter_annotations;
	my @tabs = split(/\t/, $_);
	$nSamples = scalar(@tabs) - 9;
	
	### check depth
	if($medianDP && $percentDP) {
		my $dp = $1 if ($tabs[7] =~ /DP=(\d+)/);
		if (defined $dp) {
			if ($dp < $lowDP) {
				push @filter_annotations, "LowDP";
			} elsif ($dp > $highDP) {
				push @filter_annotations, "HighDP";
			}
		} else {
			print STDERR "No DP info for $tabs[0] $tabs[1], INFO=<$tabs[7]>\n";
		}
	}
	
	### check missing
	if ($maxMissing) {
		my $ns = $1 if ($tabs[7] =~ m/NS=(\d+)/);
		if(defined $ns) {
			if ($ns < ($nSamples * ($maxMissing/100))) {
				push @filter_annotations, "highMiss";
			}
		}
	}

	### check qual
	if ($minQ) {
		if($tabs[5] ne "." && $tabs[5] < $minQ) {
			push @filter_annotations, "lowQ";
		}
	}
	
	### check mean GQ - only if GQ is present
	if ($minGQ) {
		my $gqFlag = 1 if ($tabs[8] =~ m/GQ/);
		if(defined $gqFlag) {
			my $formatString = "GQ";
			my $meanGQ = getMeanGQ($_, $formatString);
			if ($meanGQ < $minGQ) {
				push @filter_annotations, "lowGQ";
			}
		}
	}
	
	### update FILTER field with new annotations
	if (@filter_annotations) {
		if ($tabs[6] eq "." || $tabs[6] eq "PASS") {
			$tabs[6] = join(';',@filter_annotations);
		}
	} else {
		$tabs[6] = "PASS" if ($tabs[6] eq '.');
	}
	
	#print "$tabs[6]\t$tabs[5]\t$tabs[7]\n";
	
	
	### count pass|fail
	if($tabs[6] eq "PASS") {
		$passCnt++;
	}
	if($tabs[6] ne "PASS") {
		$failCnt++;
	}
	
	my $out = join("\t", @tabs);
	print STDOUT "$out";
}


print STDERR "PASS positions: $passCnt\n";
print STDERR "Fail positions: $failCnt\n";


##### subroutines ###

sub getFormatPos {
	my ($format, $formatString) = @_;
	my $cnt = -1;
	my $stringPos;
	my @tmp = split(":", $format);
	foreach my $tmp(@tmp) {
		$cnt++;
		if ($tmp eq $formatString) {
			$stringPos = $cnt;
		}
	}
	return $stringPos;
}

sub getMeanGQ {
	my ($line, $formatString) = @_;
	my $sumGQ = 0;
	my $cnt = 0;
	my $meanGQ;
	
	my @tabs = split("\t", $line);
	
	# get the position of the GQ tag
	my $GQpos = getFormatPos($tabs[8], $formatString);
	
	# get the sum of GQ over present calls
	my @samples = @tabs[ 9 .. $#tabs ];
	foreach my $tab (@samples) {
		my @fields = split(":", $tab);
		my $GQ = $fields[$GQpos];
		if($GQ > 0) {
			$sumGQ += $GQ;
			$cnt++;
		}
	}

	$meanGQ = $sumGQ/$cnt;
	#print "$sumGQ\t$cnt\t$meanGQ\n";
	return $meanGQ;
	
}


