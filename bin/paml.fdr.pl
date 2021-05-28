#!/usr/bin/perl -w
use strict;
#use lib '/ifs5/PC_PA_UN/ANIMAL/USER/GROUP2/zhangpei/perl/pm';
#use lib '/ifs5/PC_PA_UN/ANIMAL/USER/GROUP2/zhangpei/perl/pm/CDF/lib64/perl5/site_perl/5.8.5/x86_64-linux-thread-multi/';
use lib '/hwfssz1/ST_DIVERSITY/PUB/USER/panhailin/bin/perl_modul/Statistics-Multtest-0.14/lib/';
#use Pvalue_adj;
use Statistics::Multtest qw(:all);
#use MLN;
if (@ARGV < 3){
	print "perl $0 <pvalues> <method> <threshold>";
	exit();
}
my $p_file=shift;
my $Padjust=shift; #bonferroni
my $Threshold=shift;

open(IN,$p_file)||die"$!\n";
my %pvalue;
my %tag;
while(<IN>){
	chomp;
	my @a=split /[\s+\t+]/;
	$pvalue{$a[1]}=$a[-1];
	$tag{$a[1]} = $a[0];
}
close IN;

my $qvalue_p;
#if($Padjust eq "fdr"){%qvalue=&fdr(%pvalue)}
if($Padjust eq "BH"){$qvalue_p=&BH(\%pvalue)}
if($Padjust eq "hochberg"){$qvalue_p=&hochberg(\%pvalue)}
if($Padjust eq "BY"){$qvalue_p=&BY(%pvalue)}
if($Padjust eq "bonferroni"){$qvalue_p=&bonferroni(\%pvalue)}
#if($Padjust eq "none"){%qvalue=%pvalue}

foreach my $key(sort {$pvalue{$a}<=>$pvalue{$b}} keys %pvalue){
	if($qvalue_p->{$key}<$Threshold){
		print "$tag{$key}\t$key\tsignificant\t$pvalue{$key}\t$qvalue_p->{$key}\n"
	}else{
		print "$tag{$key}\t$key\tnon-significant\t$pvalue{$key}\t$qvalue_p->{$key}\n"
	}
}
