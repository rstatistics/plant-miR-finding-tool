#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;

use constant{
	WNMIN => 1,
	WNMAX => 10000,
	MAXIT => 500,
	EPS => 3.0E-30,
	FPMIN => 1.0E-30
};

############################## Variable declaration ################################

my $version = "1.0";
my ($mature_file,$clean_reads,$combined_reads,$list_file,$out_file);
my %opts=();
getopts("c:l:r:m:o:hv",\%opts);

var_check();

my $rn=0;
my %hash_mature;
my %hash_known;
my %hash_total;
my %hash_count;
my @list;

###################################  Main  #######################################
open my $fh, ">", $out_file or die "Cannot create file $out_file: $!\n";
my @samples = split /,/,$clean_reads;
Read_Mature($mature_file) if (defined $mature_file);
Check_Known_MiRNA($combined_reads) if (defined $combined_reads);
Check_MiRNA_Sequence($list_file) if (defined $list_file);
for my $file (@samples){
	$rn ++;
	Find_Total_Reads_Count($file,$rn);
	Find_Known_Count($file,$rn);
}
if ($combined_reads){
	for my $seq (keys %hash_known){
		next unless exists $hash_count{$seq};
		my $line .= "$hash_known{$seq}\t$seq"; 
		for my $i (1..@samples){
			my $x_raw;
			if (exists $hash_count{$seq}{$i}){
				$x_raw = $hash_count{$seq}{$i};
			}else{
				$x_raw = 0;
			}
			my $x_std = sprintf ("%.2f", 1000000 * $x_raw / $hash_total{$i});
			$x_std = 0.01 if $x_std < 0.01;
			$line .= "\t$x_raw\t$x_std";
		}
		push @list, $line;
	}

	if (@samples == 2){
		print $fh "mature_id\tsequence\tctl_exp\ttrt_exp\tctl_std\ttrt_std\tLFC\tP-value\tFDR\tSig\n";
		my @tmp;
		for my $line (@list){
			chomp $line;
			my ($id,$seq,$x,$y,$xstd,$ystd);
			if ($line =~ /^(\S+)\s+(\w+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\S+)$/){
				$id = $1;
				$seq = $2;
				$x = $3;
				$xstd = $4;
				$y = $5;
				$ystd = $6;
				if ($x < 0 || $y < 0 || $x + $y == 0){
					&var_error();
				}
				my $LFC = log($ystd/$xstd) / log(2);
				my $pvalue = Calculate_P_value($x, $y, $hash_total{1}, $hash_total{2});
				#my $flag = Test_Significant($LFC, $pvalue);
				my $str = "$id\t$seq\t$x\t$y\t$xstd\t$ystd\t$LFC\t$pvalue";
				push @tmp, $str;
			}else{
				print STDERR "The line $. is not recognised\n";
				next;
			}
		}
		my $index = 1;
		my %pvalue;
		my %logfc;
		my %tmps;
		for my $line (@tmp){
			if ($line =~ /^\S+\s+\w+\s+\d+\s+\d+\s+\S+\s+\S+\s+(\S+)\s+(\S+)/){
				$logfc{$index} = $1;
				$pvalue{$index} = $2;
				$tmps{$index} = $line;
				$index ++;
			}
		}

		my @ids = sort { $pvalue{$a} <=> $pvalue{$b} } keys %pvalue;
		my $number = @ids;
		my $rank = 1;
		my $fdr;
		my %fdrs;
		for my $id (@ids){
			my $p = $pvalue{$id};
			$fdr = $p * $number / $rank;
			$fdr = 1 if $fdr > 1;
			$fdrs{$id} = $fdr;
			$rank ++;
		}

		for my $i (1 .. $number){
			my $flag = Test_Significant($logfc{$i}, $fdrs{$i});
			print $fh $tmps{$i}, "\t", $fdrs{$i}, "\t", $flag, "\n";
		}
	}else{
		my $rp = @samples;
		print $fh "mature_id\tsequence";
		for my $rp (1..@samples){
			print $fh "\t", "No${rp}_exp", "\t", "No${rp}_std";
		}
		print $fh "\n";
		for my $line (@list){
			print $fh $line, "\n";
		}
	}

}elsif ($list_file){
	open my $fh_tmp, "<", $list_file or die "Cannot open file $list_file: $!\n";
	while(my $line = <$fh_tmp>){
		my ($mature_id,$star_id,$mature_arm,$star_arm,$mature_seq,$star_seq,$precursor_info);
		my ($mature_raw,$mature_std,$star_raw,$star_std);
		chomp $line;
		if ($line =~ /^(\S+)\s+(\S+)\s+(\w+)\s+(\w+)\s+\d+\s+\d+\s+(\w+)\s+(\w+)\s+(.*)/){
			$mature_id = $1;
			$star_id = $2;
			$mature_arm = $3;
			$star_arm = $4;
			$mature_seq = $5;
			$star_seq = $6;
			$precursor_info = $7;
			next unless (exists $hash_count{$mature_seq} or exists $hash_count{$star_seq});
			if (defined $hash_mature{$mature_seq}){
				$mature_id = $hash_mature{$mature_seq};
			}
			if (defined $hash_mature{$star_seq}){
				$star_id = $hash_mature{$star_seq};
			}
			my $line_out .= "$mature_id\t$mature_arm\t$mature_seq\t$star_id\t$star_arm\t$star_seq";
			for my $i (1..@samples){
				if (exists $hash_count{$mature_seq}{$i}){
					$mature_raw = $hash_count{$mature_seq}{$i};
				}else{
					$mature_raw = 0;
				}
				if (exists $hash_count{$star_seq}{$i}){
					$star_raw = $hash_count{$star_seq}{$i};
				}else{
					$star_raw = 0;
				}
				$mature_std = sprintf ("%.2f", 1000000 * $mature_raw / $hash_total{$i});
				$star_std = sprintf ("%.2f", 1000000 * $star_raw / $hash_total{$i});
				$mature_std = 0.01 if $mature_std < 0.01;
				$star_std = 0.01 if $star_std < 0.01;
				$line_out .= "\t$mature_raw\t$star_raw\t$mature_std\t$star_std";	
			}
			$line_out .= "\t$precursor_info";
			push @list, $line_out;
		}else{
			print STDERR "The line $. is ignored\n";
			next;
		}
	}
	close $fh_tmp;

	if (@samples == 2){
		print $fh "Mature_id\tMature_arm\tMature_seq\tStar_id\tStar_arm\tStar_seq\tCtl_Mature_exp\tTrt_Mature_exp\tCtl_Mature_std\tTrt_Mature_std\tMature_LFC\tMature_Pvalue\tMature_Sig\t";
		print $fh "Ctl_Star_exp\tTrt_Star_exp\tCtl_Star_std\tTrt_Star_std\tStar_LFC\tStar_Pvalue\tStar_Sig\tPrecursor_site\tPrecursor_seq\n";
		for my $line (@list){
			chomp $line;
			my ($left,$right);
			my ($ctl_mr,$ctl_sr,$ctl_ms,$ctl_ss,$trt_mr,$trt_sr,$trt_ms,$trt_ss);
			if ($line =~/^(\S+\s+\w{2}\s+\w+\s+\S+\s+\w{2}\s+\w+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(.*)/){
				$left = $1;
				$ctl_mr = $2;
				$ctl_sr = $3;
				$ctl_ms = $4;
				$ctl_ss = $5;
				$trt_mr = $6;
				$trt_sr = $7;
				$trt_ms = $8;
				$trt_ss = $9;
				$right = $10;
				if ($ctl_mr < 0 || $trt_mr < 0 || $ctl_mr + $trt_mr == 0){
					&var_error();
				}
				my $mLFC = log($trt_ms/$ctl_ms) / log(2);
				my $mpvalue = Calculate_P_value($ctl_mr, $trt_mr, $hash_total{1}, $hash_total{2});
				my $mflag = Test_Significant($mLFC, $mpvalue);
				
				my $sLFC = log($trt_ss/$ctl_ss) / log(2);
				my $spvalue = Calculate_P_value($ctl_sr, $trt_sr, $hash_total{1}, $hash_total{2});
				my $sflag = Test_Significant($sLFC, $spvalue);
				print $fh "$left\t$ctl_mr\t$trt_mr\t$ctl_ms\t$trt_ms\t$mLFC\t$mpvalue\t$mflag\t";
				print $fh "$ctl_sr\t$trt_sr\t$ctl_ss\t$trt_ss\t$sLFC\t$spvalue\t$sflag\t$right\n"
			}else{
				print STDERR "The line $. is ignored\n";
				next;
			}
		}
	}else{
		my $rp = @samples;
		print $fh "Mature_id\tMature_arm\tMature_seq\tStar_id\tStar_arm\tStar_seq";
		for my $rp (1..@samples){
			print $fh "\t", "No${rp}_mature_exp", "\t", "No${rp}_star_exp", "\t", "No${rp}_mature_std", "\t", "No${rp}_star_std";
		}
		print $fh "\tPrecursor_site\tPrecursor_seq\n";
		for my $line (@list){
			print $fh $line, "\n";
		}
	}
}
close $fh;
exit 0;

################################## Subroutine ####################################

sub Test_Significant {
	my ($LFC, $pvalue) = @_;
	my $flag = " ";
	if ( abs($LFC) >= 1 && $pvalue <= 0.01 ){
		$flag = "**";
	}elsif ( abs($LFC) >= 1 && $pvalue > 0.01 && $pvalue <= 0.05 ){
		$flag = "*";
	}
	return $flag;
}

sub Calculate_P_value {
	my ($x, $y, $n1, $n2) = @_;

	my $temp;
	my $thisproba;
	my $thisproba2;
	my $thisy;
	my $ratio;
	my $t1;
	my $t2;
	my $t3;
	my $t4;
	my $sum;
	my $ymin;
	my $ymax;
	my $argcount = 0;
	my $noup;
	my $sig;
	my $show = 0;
	my $p;

	# Check arguments and invoke the right procedure 
	if ($x > -1 && $y > -1){
		# Both x and y are defined so we compute the significance window
		$ymin =  WNMIN ;
		$ymax =  WNMAX ;
  
		$sum = 0 ; 
		$noup = 1 ; 
		$ratio = $n1 / $n2 ; 

		$p = $n1 / ($n1 + $n2)  ; 

		$thisproba = betai( ($x + 1), ($y + 1), $p ) ; 
		$thisproba2 = betai( ($y + 1), ($x + 1), (1 - $p) ) ; 
		# print "P( y <= $y | x = $x ) = $thisproba \n"; 
		# print "P( y >= $y | x = $x ) = $thisproba2 \n";
		return ($thisproba < $thisproba2 ? 2*$thisproba : 2*$thisproba2);
    } elsif( $x > -1 && $y == -1 ){

		$y = 0 ; 
		$ymin =  WNMIN ;
		$ymax =  WNMAX ;
  
		$sum = 0 ; 
		$noup = 1 ; 
		$ratio = $n1 / $n2 ; 

		$p = $n1 / ($n1 + $n2); 

		$y = 0 ; 
		# print "x = \$x y= \$y sig = \$sig\n" 
		while( 1 ){
			$thisproba = betai( ($x + 1), ($y + 1), $p) ; 
			$thisproba2 = betai( ($y + 1), ($x + 1), (1 - $p) ) ; 

			if( $thisproba < $sig / 2.0 ){
				$ymin = $y ; 
			}

			if( $ymax == WNMAX && $thisproba2 < $sig / 2.0 ){
				$ymax = $y ; 
				last; 
			}
			$y++ ;
		}
		$thisproba = betai( ($x + 1), ($ymin + 1), $p ) ; 
		$thisproba2 = betai( ($ymax + 1), ($x + 1), (1 - $p) ) ; 
		# print "P( y <= $ymin | x = $x ) = $thisproba \n"; 
		# print "P( y >= $ymax | x = $x ) = $thisproba2 \n"; 
		return ($thisproba < $thisproba2 ? 2*$thisproba : 2*$thisproba2);

		if( $ymin == -1 ){
			print "$x *--$ymax\n";
		} else {
			print "$x $ymin--$ymax\n";
		}
	
	} elsif( $show == 1 ){
	
		print "$x $y C($y | $x ) = $thisproba    $thisproba2\n"; 
		$y++ ; 
		if( $thisproba > 0.9999999999 ){

		}
	}
}

sub betacf {
	my ($a, $b, $x) = @_;
	my $m;
	my $m2;
	my $aa;
	my $c;
	my $d;
	my $del;
	my $h;
	my $qab;
	my $qam;
	my $qap;
	
	$qab = $a + $b;
	$qap = $a + 1.0;
	$qam = $a - 1.0;
	$c = 1.0;
	$d = 1.0 - $qab * $x / $qap;
	if (abs($d) < FPMIN) { $d = FPMIN}
	$d = 1.0 / $d;
	$h = $d;
	for ($m=1;$m<=MAXIT;$m++) {
		$m2 = 2 * $m;
		$aa = $m * ($b - $m) * $x / (($qam + $m2) * ($a + $m2));
		$d = 1.0 + $aa * $d;
		if (abs($d) < FPMIN) { $d = FPMIN }
		$c = 1.0 + $aa / $c;
		if (abs($c) < FPMIN) { $c = FPMIN }
		$d = 1.0 / $d;
		$h *= $d * $c;
		$aa = -($a + $m) * ($qab + $m) * $x / (($a + $m2) * ($qap + $m2));
		$d = 1.0 + $aa * $d;
		if (abs($d) < FPMIN) { $d = FPMIN }
		$c = 1.0 + $aa / $c;
		if (abs($c) < FPMIN) { $c = FPMIN }
		$d = 1.0 / $d;
		$del = $d * $c;
		$h *= $del;
		if (abs($del-1.0) < EPS) { last; }
	}
	if ($m > MAXIT){
	  print "a or b too big, or MAXIT too small in betacf";
	  exit; 
	}
	return $h;
}

sub gammln {
	my ($xx) = @_;
	my $x;
	my $y;
	my $tmp;
	my $ser;
	my @cof = (76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5
			);
	my $j;

	$y = $x = $xx;
	$tmp = $x + 5.5;
	$tmp -= ($x+0.5) * log($tmp);
	$ser = 1.000000000190015;
	for ($j=0;$j<=5;$j++){
		$ser += $cof[$j] / ++$y;
	}
	return ( -$tmp + log(2.5066282746310005*$ser/$x) );
}

sub betai {
	my ($a, $b, $x) = @_;
	my $bt;

	if ($x < 0.0 || $x > 1.0) { 
	  print "Bad x in routine betai"; 
	  exit; 
	}
	if ($x == 0.0 || $x == 1.0){
		$bt = 0.0;
	}else{
		$bt = exp(gammln($a+$b)-gammln($a)-gammln($b)+$a*log($x)+$b*log(1.0-$x));
	}
	if ($x < ($a+1.0)/($a+$b+2.0)){
		return $bt * betacf($a, $b, $x) / $a;
	}else{
		return 1.0 - $bt * betacf($b, $a, 1.0-$x) / $b;
	}
}

sub Read_Mature{
	my ($file) = @_;
	my ($id, $desc, $seq) = ();
	open FILE, "<$file" or die "Can not open $file: $!\n";
	while (<FILE>){
		chomp;
		if (/^>(\S+)\s*(.*)/){
			$id = $1;
			$desc = $2;
			$seq = "";
			while (<FILE>){
				chomp;
				if (/^>(\S+)\s*(.*)/){
					if (exists $hash_mature{$seq}){
						$hash_mature{$seq} .= ",$id";
					}else{
						$hash_mature{$seq} = $id;
					}
					$id = $1;
					$desc = $2;
					$seq = "";
					next;
				}
				$seq .= $_;
			}
		}
	}
	if (exists $hash_mature{$seq}){
		$hash_mature{$seq} .= ",$id";
	}else{
		$hash_mature{$seq} = $id;
	}
	close FILE;
	return;
}

sub Check_MiRNA_Sequence{
	my ($file) = @_;
	open my $fh, "<", $file or die "Cannot open file $file: $!\n";
	while(<$fh>){
		chomp;
		my ($mature_id,$star_id,$mature_arm,$star_arm,$mature_seq,$star_seq,$precursor_info);
		my ($mature_raw,$mature_std,$star_raw,$star_std);
		if (/^(\S+)\s+(\S+)\s+(\w{2})\s+(\w{2})\s+\d+\s+\d+\s+(\w+)\s+(\w+)\s+(.*)/){
			$mature_id = $1;
			$star_id = $2;
			$mature_arm = $3;
			$star_arm = $4;
			$mature_seq = $5;
			$star_seq = $6;
			$precursor_info = $7;
			$hash_known{$mature_seq} = 1;
			$hash_known{$star_seq} = 1;
		}else{
			print STDERR "The line $. is ignored\n";
			next;
		}
	}
	close $fh;
	return;
}

sub Check_Known_MiRNA{
	my ($file) = @_;
	my $seq;
	open my $fh, "<$file" or die "Cannot open file $file: $!\n";
	while (<$fh>){
		chomp;
		if (/^>\S+.*/){
			$seq = "";
			while (<$fh>){
				chomp;
				if (/^>\S+.*/){
					if (defined $hash_mature{$seq}){
						$hash_known{$seq}=$hash_mature{$seq};
					}
					$seq = "";
					next;
				}
				$seq .= $_;
			}
		}
	}
	if (defined $hash_mature{$seq}){
	    $hash_known{$seq}=$hash_mature{$seq};
	}
	close $fh;
	return;
}

sub Find_Known_Count{
	my ($file,$rn) = @_;
	my ($cnt,$seq) = ();
	open FILE, "<$file" or die "Cannot open $file: $!\n";
	while (<FILE>){
		chomp;
		if (/^>\S+\_x(\d+).*/){
			$cnt = $1;
			$seq = "";
			while (<FILE>){
				chomp;
				if (/^>\S+\_x(\d+).*/){
					if (defined $hash_known{$seq}){
						$hash_count{$seq}{$rn} = $cnt;
					}
					$cnt = $1;
					$seq = "";
					next;
				}
				$seq .= $_;
			}
		}
	}
	if (defined $hash_known{$seq}){
	    $hash_count{$seq}{$rn} = $cnt;
	}
	close FILE;
	return;
}

sub Find_Total_Reads_Count{
	my ($file,$rn) = @_;
	my $count_reads = 0;
	open FILE, "<", $file or die "Could not open file $file: $!\n";
	while (<FILE>){
		if (/_x(\d+)/){
			$count_reads += $1;
		}
	}
	close FILE;
	$hash_total{$rn} = $count_reads;
	return;
}

sub var_check{

    if (exists $opts{v}){
        die "\ndiff_miRNA_expression_analysis.pl version = $version\n\n";
    }
    if (exists $opts{h}){
        &var_error();
    }
    if ($opts{c}){
        $combined_reads = $opts{c};
    }
    if ($opts{m}){
        $mature_file = $opts{m};
    }
    if ($opts{l}){
        $list_file = $opts{l};
    }
    unless ($opts{'r'} && $opts{'l'} && $opts{'o'} ||
            $opts{'m'} && $opts{'r'} && $opts{'c'} && $opts{'o'})
    {
        &var_error();
    }
    if ($opts{c} && $opts{l}){
        &var_error();
    }
    if ($opts{r}){
        $clean_reads = $opts{r};
    }else{
        &var_error();
    }
    if ($opts{o}){
        $out_file = $opts{o};
    }else{
        &var_error();
    }
    return;
}

sub var_error{

	print "\n";
	print " diff_miRNA_expression_analysis.pl aims at analyzing differential microRNA expression \n";
	print " data. It will output a list with standard gene count, LFC (log2 transformed foldhange),\n";
    print " pvalue, FDR and significance, etc.\n";
	print "\n";
	print " Usage:\n\n";
    print "     diff_miRNA_expression_analysis.pl -m <miRNA.fa> -c <combined_reads.fa> \\\n";
	print "                                       -r file_1.fa[,...,file_N.fa] -o <file_out>\n";
	print " or\n";
	print "     diff_miRNA_expression_analysis.pl -l <exp_list> -r file1.fa[,...,fileN.fa] \\\n";
	print "                                       -o <file_out> [-m <miRNA.fa>]\n";
    print "\n";
	print " Options:\n";
	print " -c <reads.fa>                   a combined multi-fasta file of clean reads generated\n";
    print "                                 by format_clean_reads.pl\n";
	print " -l <exp.list>                   a tab-delimited file contains microRNAs sequence info\n";
	print " -m <miRNA.fa>                   a multi-fasta format file contains known microRNAs\n";
	print " -r <file_1.fa[,...,file_N.fa>   a comma-separated list of files with small RNA reads\n";
    print "                                 in multi-fasta format\n";
	print " -o <exp.out>                    a tab-delimited output file\n";
	print "\n";
	print " Note:\n";
	print " I.    option 'c' and option 'l' are mutually exclusive;\n";
	print " II.   The input reads should be formated by format_clean_reads.pl. The reads\n";
	print "       file should have each entry with unique sequence, the entry must be as\n";
	print "       '>CR_xN', while C represents one or more letters, R represents a non-redundant\n";
    print "       running number, and N represents current read count, _x is the separator.\n";
	print " III.  Pvalues, FDR and significance will be calculated only in the case of two samples.\n";
	print "\n";
	exit 1;

}

############################################  End  ##################################################



