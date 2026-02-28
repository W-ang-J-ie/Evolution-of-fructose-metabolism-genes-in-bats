#!/usr/bin/perl
use warnings ; use strict;

my ($nw,$rst,$pre,@cls) = @ARGV;

#open OUT, ">$pre.conver.txt" or die $!;

my %ha;
my %cha;
my %cnm;


for ( my $i = 0; $i < @cls -1 ; $i ++ ){
	my $c1 = $cls[$i];
	
	for (my $j = $i+1; $j < @cls; $j ++){

		my $c2 = $cls[$j];

		warn "perl ~/Chuan/fruc/Codes/convergent.sub.pl $nw $c1 $c2  $rst\n";
		open IN , "perl ~/Chuan/fruc/Codes/convergent.sub.pl $nw $c1 $c2 $rst | " or die $!;
		while(<IN>){
			chomp;
			my @ar = split /\t/;
			next unless ($ar[1] =~ /convergent/);
			next unless ($ar[7] eq $ar[9]);
	
			my $sub = $ar[7];
			my $pos = $ar[0];
			
			my ($n1,$c1) = split /:/,$ar[2];
			my ($n2,$c2) = split /:/,$ar[3];

			$cnm{$c1} = $n1;
			$cnm{$c2} = $n2;
			
			my ($b1) = $ar[4] =~ (/^(\d+)\-/);
			my ($b2) = $ar[5] =~ (/^(\d+)\-/);
			
			$cha{$b1}{N} = "$n1";
			$cha{$b1}{C} = $c1;
			$cha{$b2}{N} = $n2;
			$cha{$b2}{C} = $c2;
			
			$ha{$pos}{$sub}{$c1} = $b1;
			$ha{$pos}{$sub}{$c2} = $b2;
		}
	}
}



my %cfa;
open RST, $rst or die $!;
while(<RST>){
	chomp;
	if(/^node #(\d+)\s+(.+)/){
		my $nid = $1;
		my $seq = $2;
		#		warn "seen $nid\n";
		if($cha{$nid}){
			my $c = $cha{$nid}{C};
			my $n = $cha{$nid}{N};
			$seq =~ s/ //g;
			open FA, ">07_proven/$pre.$c.fa" or die $!;
			warn "$nid $c sequence\n";
			print FA ">$pre.$n.$c.$nid\n$seq\n";
			close FA;
		}
	}
}
close RST;


my %vha;

foreach my $pos( sort {$a <=> $b} keys %ha){
	foreach my $sub ( keys %{$ha{$pos}} ){
		warn "in loop\n";
		my @cls = sort keys %{$ha{$pos}{$sub}};
		my @brs = map{$ha{$pos}{$sub}{$_}} @cls;
		my @nms = map{$cnm{$_}} @cls;
		my $cl = join "=", @cls;
		$sub =~ s/-/$pos/;
		
		foreach my $c (@cls){
			my $fh;
			if($vha{$c}){
				$fh = $vha{$c};
			}else{
				open my $f, ">07_proven/$pre.$c.var" or die $!;
				$fh = $f;
				$vha{$c} = $f;
			}
			warn "$pos $sub $c  $sub\n";
			print $fh "$sub\n";
		}
		print "$pre\t$pos\t$sub\t@cls\t@nms\t@brs\n";
	}
}


