#!/usr/bin/perl
use warnings; use strict;
use Bio::SeqIO;

my($par_tree, $fg , $bg , $seq_f  , $pre, $test ,$para ) = @ARGV;

######
srand(1); # 2022年5月8日
#srand($seq_f);
my $rand = rand();
######

die "$0\t'A,B  C,D' 'E,F +G -I' seq_file back_trans_f prefix\n" unless ($fg);

unless(-d $pre){
	warn "making $pre\n";
	mkdir $pre;
}
unless(-d "$pre/transver_dir"){
	warn "making $pre/transver_dir\n";
	mkdir "$pre/transver_dir";
}

#####

if ($para){
	$para = "--codeml_param $para";
}else{
	$para = "";
}

#######  backgroud  ######

my %labels;
my @og = split /\s+/, $bg;
my $fbg = shift @og;
if($fbg =~ /ROOT/g){  ### fbg = the first background
	$fbg = "-r ." ;
}else{
	$fbg =~ s/,/ /;
}

warn "nw_clade  $par_tree  $fbg | nw_labels  -I - |\n";
open NW, "nw_clade  $par_tree  $fbg | nw_labels  -I - | " or die $!;
while(<NW>){
	chomp;
	$labels{$_} = 1;
}
close NW;

foreach my $b (@og){
	my $sign = $b =~/^-/ ? 0 : 1;
	$b =~ s/^([\+-]*)//;
	$b =~ s/,/ /;
	#print "nw_clade  $par_tree  $b | nw_labels  -I - |\n";
	open NW, "nw_clade  $par_tree  $b | nw_labels  -I - | " or die $!;
	while(<NW>){
		chomp;
		$labels{$_} = $sign;
	}
	close NW;

}

######  backgroud end  ######

######   fasta   ######
warn  "$pre.aln.fa\n";
open FA, ">$pre.raw.fa" or die $!;
open CMD, ">$pre.cmds.txt" or die $!;
open LOS, ">$pre.loss.txt" or die $!;

my %allbgs;

my $seq_in = Bio::SeqIO -> new(-file => $seq_f, -format => "fasta");
while(my $seq_obj = $seq_in -> next_seq){
	my $id = $seq_obj -> id;
	my $desc = $seq_obj -> desc;
	if($labels{$id}){
		$allbgs{$id} = 1;
		my $alb = $seq_obj -> seq;
		print FA ">$id\n$alb\n";
	}
}
close FA;

if(-e "$pre/aln.success.$rand"){
	1;
}else{
	my $cmd = "/storage/Codes/prank/bin_new/prank -d=$pre.raw.fa -o=$pre.aln -codon "; 
	print CMD "$cmd\n";
	my $out = system ($cmd);
	print "$out";
	`touch $pre/aln.success.$rand`;
}

my $triml = "trimal  -in $pre.aln.best.fas  -phylip_paml -out $pre.aln.phy  -htmlout  $pre.aln.html  -resoverlap 0.80 -seqoverlap 90  -gt 0.9 -st 0.001 -colnumbering >$pre.cols 2>/dev/null ";

print CMD "$triml\n";
my $out = system ($triml);

#########################

##### find gene ids in phylip alignment   ######
my $chi_aln = "$pre.aln.phy" unless ($out);

my %kps;
open PHY, "$pre.aln.phy" or die $!;
<PHY>;
while(<PHY>){
	chomp;
	my ($s) =  $_ =~ /^(\w+)/;
	$kps{$s} = 1  if $s;

}
close PHY;

##### find lost species due to triming  ####
foreach my $s (keys %allbgs){
	unless($kps{$s}){
		print LOS "$s\tTRIM\n";
	}
}
close LOS;



#############################

####    prune trees   ####
my @disc;
my @allbgs_neg;
open PP, "nw_labels -I $par_tree | " or die $1;
while(<PP>){
	chomp;
	push @disc, $_ unless $kps{$_};
	push @allbgs_neg, $_  unless $allbgs{$_};
}

my $cmd  ;
if(@allbgs_neg){
	$cmd = "nw_prune  $par_tree @allbgs_neg > $pre.all.nw";
}else{
	$cmd = "cp $par_tree $pre.all.nw";
}
print CMD "$cmd\n";
system($cmd);

if(@disc){
	$cmd = "nw_prune  $par_tree @disc > $pre.aln.nw";
}else{
	$cmd = "cp $par_tree $pre.aln.nw";
}
print CMD "$cmd\n";
system($cmd);


my $chi_tree = "$pre.aln.nw" unless ($out);

##################################         
#########  MODELS sorter   #######

my(@models,@models_s,@models_b,@models_c);

@models = split /[,\s]/,$test;
@models  = do { my %seen; grep { !$seen{$_}++ } @models};
#my @models  = do { my %seen; grep { !$seen{$_}++ } (split /[,\s]/,$test) };

foreach my $m(@models){
	if($m =~ /^M|fb|XX/){
		push @models_s, $m;
	}else{
		if($m =~ /^b_|bsC|bsD/){
			push @models_c,$m;
		}
		push @models_b,$m unless ($m =~ /bsC|bsD/);
	}
}

warn "\tMODELS:
	B:@models_b
	S:@models_s
	C:@models_c\n";

#################################
########   TEST types   #########
my (@tests ,@tests_s, @tests_b, @tests_c);


@tests = split /\s/,$test;

foreach my $t (@tests){
	next unless ($t =~ /,/);
	if($t =~ /bs/ and $t !~ /bsC|bsD/ ){
		push @tests_b,$t;
	}
	if($t =~ /b_|bsC|bsD/){
		push @tests_c,$t;	
	}
	unless($t =~/bs|b_/){
		push @tests_s,$t;
	}
}

my $test_spar = @tests_s ? "--test @tests_s" : "" ;
my $test_bpar = @tests_b ? "@tests_b" : "" ;
my $test_cpar = @tests_c ? "@tests_c" : "" ;

warn "\tTEST:
	B:@tests_b
	S:@tests_s
	C:@tests_c\n";



##############

my @fors = split /\s+/, $fg ;
#warn "@fors\n";


my @for_brns;

for( my $i = 0 ; $i < @fors; $i++ ){
	my $for = $fors[$i];
	warn "fore $for\n";

	if($for eq "LEAVES"){
		push @for_brns , "leaves";
	}elsif($for eq "ROOT"){
		push @for_brns , "ROOT";	
	}else{
		my @ca = rep_gener($for,$par_tree,$chi_tree );
		print CMD "#$ca[0],$ca[1]\n";
		push @for_brns,\@ca;
	}
}





#### test models container ####
my @models_t;
####################





########### SITE  models  #####
my $site_join = join ".",@models_s;
print "ete3 evol -t $chi_tree --alg $chi_aln --models @models_s  $test_spar  --clear_tree --resume -C 0 -o $pre $para> $pre/st.$site_join.log  # SITE \n" if ( @models_s );


###### BRANCH MODEL  ###

if (@models_b){
	for(my$i = 0 ; $i < @for_brns ; $i++){
		my $for = $for_brns[$i];
		warn "ffffore @$for\n";
		
		my $bmark;
		if($for eq "leaves"){
			$bmark = "--leaves";
		}elsif($for eq "ROOT"){
			next;
		}else{
			$bmark = "--mark $$for[0],,$$for[1]";
		}
		print "ete3 evol -t $chi_tree --alg $chi_aln $bmark  --models  @models_b  --clear_tree --resume -C 0 -o $pre $para > /dev/null   # BRANCH_SITE\n";
	
		@models_t  = do { my %seen; grep { !$seen{$_}++ } (split /[,\s]/,$test_bpar) };	
		
		print "ete3 evol -t $chi_tree --alg $chi_aln --models @models_t  $bmark  --test $test_bpar  --clear_tree --resume -C 0 -o $pre $para >$pre/bs.$$for[0]--$$for[1].log # TEST_BRANCH_SITE \n" if ($test_bpar);
	}
}

### CLADE model  ###
my @fma;
my @clas;
if(@models_c){
	for(my$i = 0 ; $i < @for_brns ; $i++){
		my $for = $for_brns[$i];

		my $bmark;
		if($for eq "leaves"){
			next;
		}else{
			$bmark = "--mark $$for[0],,,$$for[1]";
			push @clas, "$$for[0],,,$$for[1]";
			push @fma,  "$$for[0]---$$for[1]";
		}
		print "ete3 evol -t $chi_tree --alg $chi_aln $bmark  --models  @models_c  --clear_tree --resume -C 0 -o $pre $para > /dev/null   # CLADE_EACH\n";
		@models_t  = do { my %seen; grep { !$seen{$_}++ } (split /[,\s]/,$test_cpar) };
		print "ete3 evol -t $chi_tree --alg $chi_aln --models @models_t  $bmark  --test $test_cpar  --clear_tree --resume -C 0 -o $pre $para >$pre/cl.$$for[0]---$$for[1].log # TEST_CLADE_EACG \n" if ($test_bpar);
	}
}


my $fmj = join "-",@fma;
my $cmj = join ",", @clas;


print "ete3 evol  -t $chi_tree --alg $chi_aln --mark $cmj  --models @models_c --clear_tree --resume -C 0 -o $pre  $para > $pre/cl.$fmj.log # CLADE_ALL \n " if @models_c ;

@models_t  = do { my %seen; grep { !$seen{$_}++ } (split /[,\s]/,$test_cpar) };

print "ete3 evol  -t $chi_tree --alg $chi_aln --mark $cmj  --models @models_t  --test  $test_cpar   --clear_tree --resume -C 0 -o $pre $para  > $pre/cl.$fmj.log # TEST_CLADE_ALL \n " if($test_cpar);

#################################
#
###### generate transversiong  branch site  script #####
my %trans;

foreach my $for (@for_brns){
	if($for eq "ROOT"){
		$trans{A} = 1;
		last;
	}elsif($for eq "leaves"){
		next;
	}else{
		open my $fh, "nw_clade $chi_tree  @{$for} | nw_labels - -I |" or die $!;
		while(<$fh>){
			chomp;
			$trans{$_} = 1;
		}
		close $fh;
	}
}

open SH, ">$pre.transversing_bs.sh" or die $!;

my @nodes;
if($trans{A}){
    chomp(my $b_c = `nw_labels $chi_tree | wc -l `);
    $b_c  ++;
    @nodes = 1..$b_c;
}else{
    open ETE, "ete3 evol -t $chi_tree --node_ids | " or die $!;
 	while(<ETE>){
	chomp;
    my ($nid, $snm) ;
        if(/\s+(\d+)\s+:\s+([\w,\s]+)/){
            $nid = $1;
            $snm = $2;
        }else{
            next;
		}
        my $tag = 1;
        while( $snm =~ /(\w+)/g ){
            #print "====='$1'\n";
            unless($trans{$1}){
                $tag = 0;
                 #print "00000\n";
                last;
            }
        }
        if($tag){
           push @nodes, $nid;
        }
    }
}


print SH "ete3 evol -t $chi_tree --alg $chi_aln --models  @models_s  --clear_tree --resume -C 0 -o $pre >/dev/null # SITE \n";
print SH "ete3 evol -t $chi_tree --alg $chi_aln --models  @models_s --test $test  --clear_tree --resume -C 0 -o $pre > $pre/transver_dir/bs.evol.log  # TEST \n" unless(@models_b);

foreach my $node ( @nodes ){
    print SH "ete3 evol -t $chi_tree --alg $chi_aln --models @models_b --mark $node  --clear_tree --resume -C 0 -o $pre >/dev/null  #BRANCH_SITE\n" if (@models_b);
    print SH "ete3 evol -t $chi_tree --alg $chi_aln --models @models   --mark $node   --test $test  --clear_tree --resume -C 0 -o $pre > $pre/transver_dir/$node.log  # TEST \n" if (@models_b);
}

close SH;
    warn "scrits in $pre.transversing_bs.sh. script ending\n";

    ##############


sub rep_gener{

	my($for,$par_t,$chi_t ) = @_;
	warn "begin rep_genr....\n";    
	unless($for =~ /,/){
		return ($for,$for);
	}
	$for =~ s/,/ /g;

    my %labs;

    open my $fh, "nw_clade $par_t  $for | nw_labels - -I |" or die $!;
    while(<$fh>){
        chomp;
        $labs{$_} = 1;
    }

    close $fh;

    open ETE, "ete3 evol -t $chi_t --node_ids | " or die $!;
    my %iha;

    while(<ETE>){
        chomp;
        my ($nid, $snm) ;

        if(/\s+(\d+)\s+:\s+([\w,\s]+)/){
            $nid = $1;
            $snm = $2;
        }else{
            next;
        }

        my ($ct,$pure);
        $pure = 1;
        while( $snm =~ /(\w+)/g ){
            if($pure){
                if($labs{$1}){
                    $ct ++;
                }else{
                    $pure = 0;
                }
            }
            push @{$iha{S}{$1}} ,$nid;
            push @{$iha{N}{$nid}} ,$1;
        }
        if($pure){
            $iha{C}{$ct} = $nid;
        }
    }


    my ($bc) = (sort {$b <=> $a} keys %{$iha{C}})[0];

    return  unless $bc;

    my $bn = $iha{C}{$bc};

    my @ss = @{$iha{N}{$bn}};
    if (grep /^ROOT$/,@ss){
        @ss = keys %{$iha{S}};
    }
    my @ca ;
    my $rcd = 1000000;
    foreach my $i(@ss){
        foreach my $j(@ss){
            my @a = @{$iha{S}{$i}};
            my @b = @{$iha{S}{$j}};
            my @inter = grep { defined } @{ { map { lc ,=> $_ } @a } }  { map {lc } @b };
            if(@inter < $rcd){
                @ca = $i eq $j ?  $i : ($i,$j);
                $rcd = @inter;
            }
        }
    }

    if(@ca < 2){
        warn  "small node $for \n" ;
		return;
    }
	return (@ca);

}
