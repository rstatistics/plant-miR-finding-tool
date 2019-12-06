#!/usr/bin/perl

use warnings;
use strict;
use Bio::DB::Fasta;
use Bio::DB::Sam;
use Getopt::Std;

my %opts = ();
my $version = "1.0";
getopts("b:l:m:o:r:hv", \%opts);
&var_check();

my ($bam_file, $min_freq, $reference, $out_file, $max_len);
my (%hash_mature, %hash_island, %hash_pos, %hash_len);
my $count_excisions = 0;

my $database = Bio::DB::Fasta->new($reference);
# get all subject ids
my @ids = $database->get_all_ids;
my $sam = Bio::DB::Sam->new(-bam => $bam_file, 
					        -fasta_file => $reference);

open PRECURSOR, ">", $out_file or die "Cannot create $out_file: $!\n";

get_subject_length(\@ids,\%hash_len);
parse_bam_file($bam_file,$reference);

close PRECURSOR;

exit 0;

############################################## SUBROUTINE ##########################################

sub parse_bam_file{
	my ($bam,$ref) = @_;
	open (BAM, "samtools view $bam |");
	my (%hash_plus,%hash_minus);
    my $subject_prev = "";
	while (<BAM>){		
		if (/^\S+\_x(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+\d+\s+(\d+)M\s+\*\s+0\s+0\s+(\w+)\s+/){
			my $freq = $1;
			my $strand;
			if ($2 == 0){ $strand = '+';}elsif ($2 == 16){ $strand = '-';}
			my $subject = $3;
			my $beg = $4;
			my $len = $5;
			my $end = $beg + $len -1;
			my $beg_min;
			my $end_max;
			my $freq_max;
            if ($subject_prev eq ""){
                $subject_prev = $subject;
            }elsif ($subject_prev eq $subject){
            }else{
                if (%hash_plus){
	            	my $subject_pre = (keys %hash_plus)[0];
            		my $beg_min = (sort {$a<=>$b} keys %{$hash_plus{$subject_pre}})[0];
            		my $end_max = (sort {$b<=>$a} keys %{$hash_plus{$subject_pre}{$beg_min}})[0];
            		my $freq_max = $hash_plus{$subject_pre}{$beg_min}{$end_max};
            		$hash_island{$subject_pre}{"+"}{$beg_min}{$end_max} = $freq_max;
            		%hash_plus = ();
            	}
            	if (%hash_minus){
            		my $subject_pre = (keys %hash_minus)[0];
            		my $beg_min = (sort {$a<=>$b} keys %{$hash_minus{$subject_pre}})[0];
            		my $end_max = (sort {$b<=>$a} keys %{$hash_minus{$subject_pre}{$beg_min}})[0];
            		my $freq_max = $hash_minus{$subject_pre}{$beg_min}{$end_max};
            		$hash_island{$subject_pre}{"-"}{$beg_min}{$end_max} = $freq_max;
            		%hash_minus = ();
            	}
                my $subject_prev_lng = $hash_len{$subject_prev};
                get_precursors($subject_prev,$subject_prev_lng,\%hash_island);    
                %hash_island = ();
                $subject_prev = $subject;
            }
			if ($strand eq '+'){
				if (exists $hash_plus{$subject}){
					$beg_min = (sort {$a<=>$b} keys %{$hash_plus{$subject}})[0];
					$end_max = (sort {$b<=>$a} keys %{$hash_plus{$subject}{$beg_min}})[0];
					$freq_max = $hash_plus{$subject}{$beg_min}{$end_max};
					if (overlap($beg_min,$end_max,$beg,$end)){
						$beg_min = $beg_min <= $beg ? $beg_min : $beg;
						$end_max = $end_max >= $end ? $end_max : $end;
						$freq_max = $freq_max >= $freq ? $freq_max : $freq;
						$hash_plus{$subject}{$beg_min}{$end_max} = $freq_max;
					}elsif ($end_max-$beg_min>100){
						%hash_plus = ();
						$hash_plus{$subject}{$beg}{$end} = $freq;
					}else{
						$hash_island{$subject}{$strand}{$beg_min}{$end_max} = $freq_max;
						%hash_plus = ();
						$hash_plus{$subject}{$beg}{$end} = $freq;
					}
				}else{
					if (%hash_plus){
						my $subject_pre = (keys %hash_plus)[0];
						$beg_min = (sort {$a<=>$b} keys %{$hash_plus{$subject_pre}})[0];
						$end_max = (sort {$b<=>$a} keys %{$hash_plus{$subject_pre}{$beg_min}})[0];
						$freq_max = $hash_plus{$subject_pre}{$beg_min}{$end_max};
						$hash_island{$subject_pre}{"+"}{$beg_min}{$end_max} = $freq_max;
						%hash_plus = ();
					}
					$hash_plus{$subject}{$beg}{$end} = $freq;
				}
			
			}elsif($strand eq '-'){
				if (exists $hash_minus{$subject}){
					$beg_min = (sort {$a<=>$b} keys %{$hash_minus{$subject}})[0];
					$end_max = (sort {$b<=>$a} keys %{$hash_minus{$subject}{$beg_min}})[0];
					$freq_max = $hash_minus{$subject}{$beg_min}{$end_max};
					if (overlap($beg_min,$end_max,$beg,$end)){
						$beg_min = $beg_min <= $beg ? $beg_min : $beg;
						$end_max = $end_max >= $end ? $end_max : $end;
						$freq_max = $freq_max >= $freq ? $freq_max : $freq;
						$hash_minus{$subject}{$beg_min}{$end_max} = $freq_max;
					}elsif ($end_max-$beg_min>100){
						%hash_minus = ();
						$hash_minus{$subject}{$beg}{$end} = $freq;
					}else{
						$hash_island{$subject}{$strand}{$beg_min}{$end_max} = $freq_max;
						%hash_minus = ();
						$hash_minus{$subject}{$beg}{$end} = $freq;
					}
				}else{
					if (%hash_minus){
						my $subject_pre = (keys %hash_minus)[0];
						$beg_min = (sort {$a<=>$b} keys %{$hash_minus{$subject_pre}})[0];
						$end_max = (sort {$b<=>$a} keys %{$hash_minus{$subject_pre}{$beg_min}})[0];
						$freq_max = $hash_minus{$subject_pre}{$beg_min}{$end_max};
						$hash_island{$subject_pre}{"-"}{$beg_min}{$end_max} = $freq_max;
						%hash_minus = ();
					}
					$hash_minus{$subject}{$beg}{$end} = $freq;
				}
			
			}
		}
	}
	if (%hash_plus){
		my $subject_pre = (keys %hash_plus)[0];
		my $beg_min = (sort {$a<=>$b} keys %{$hash_plus{$subject_pre}})[0];
		my $end_max = (sort {$b<=>$a} keys %{$hash_plus{$subject_pre}{$beg_min}})[0];
		my $freq_max = $hash_plus{$subject_pre}{$beg_min}{$end_max};
		$hash_island{$subject_pre}{"+"}{$beg_min}{$end_max} = $freq_max;
		%hash_plus = ();
	}
	if (%hash_minus){
		my $subject_pre = (keys %hash_minus)[0];
		my $beg_min = (sort {$a<=>$b} keys %{$hash_minus{$subject_pre}})[0];
		my $end_max = (sort {$b<=>$a} keys %{$hash_minus{$subject_pre}{$beg_min}})[0];
		my $freq_max = $hash_minus{$subject_pre}{$beg_min}{$end_max};
		$hash_island{$subject_pre}{"-"}{$beg_min}{$end_max} = $freq_max;
		%hash_minus = ();
	}
    my $subject_prev_lng = $hash_len{$subject_prev};
    get_precursors($subject_prev,$subject_prev_lng,\%hash_island);    
	close BAM;
	return;
}

sub get_precursors{
	my ($db,$db_lng,$island) = @_;
	my $flag;
		my @strands = sort keys %{$$island{$db}};
		for my $strand (@strands){
			if ($strand eq '+'){$flag = 1;}elsif($strand eq '-'){$flag = -1;}
			my ($db_beg,$db_end,$sub_beg,$sub_end,$beg_min,$end_max);
			my @db_begs = sort {$a<=>$b} keys %{$$island{$db}{$strand}};
			for my $i (0..$#db_begs){
				$db_beg = $db_begs[$i];
				$db_end = (keys %{$$island{$db}{$strand}{$db_beg}})[0];
				my $primary_freq = $$island{$db}{$strand}{$db_beg}{$db_end};
				unless (defined $beg_min && defined $end_max){
					$beg_min = $db_beg;
				}
				
				for my $j ($i+1..$#db_begs){
					$sub_beg = $db_begs[$j];
					next if $sub_beg <= $beg_min;
					$sub_end = (keys %{$$island{$db}{$strand}{$sub_beg}})[0];
					my $present_freq = $$island{$db}{$strand}{$sub_beg}{$sub_end};
					my $max_freq = $primary_freq >= $present_freq ? $primary_freq : $present_freq;
					if ($sub_end - $beg_min <= $max_len){
						next if $max_freq < $min_freq;
						$end_max = $sub_end;
						my $excise_beg = ($beg_min-5) >= 1 ? ($beg_min-5) : 1;
						my $excise_end = $db_lng <= ($end_max+5) ? $db_lng : ($end_max+5);
						my $strand_ratio = calc_strand_ratio($db,$strand,$excise_beg,$excise_end);
                        next if ($strand_ratio < 0.7);
                        my $excise_seq = $database->seq("$db:$excise_beg,$excise_end/$flag");
						if ($excise_seq !~ /N/){
							print PRECURSOR ">$db\_$count_excisions strand:$strand excise_beg:$excise_beg excise_end:$excise_end\n$excise_seq\n";
							$count_excisions ++;
						}

					}else{
						#if stack length beyond 500nt, then terminate current loop
						($beg_min, $end_max)=();
						last;
					}
			
				}
                #if outer for loop changes, initialize the variables
                ($beg_min, $end_max)=();
		}
		
	}
    return;
}

sub calc_strand_ratio{
	my ($primary_id,$subject_strand,$subject_beg,$subject_end) = @_;
	my $forward_cnt = 0;
	my $reversed_cnt = 0;
	my $total_cnt = 0;
	my @alignments = $sam->get_features_by_location(-seq_id => $primary_id,
													-start => $subject_beg,
													-end => $subject_end);
	for my $line (@alignments){
		my $query    		= $line->qname;
		my $query_seq 		= $line->query->dna;
		my $db_start		= $line->start;
		my $db_stop			= $line->end;
		my $query_strand	= $line->strand;
		next if $db_stop > $subject_end || $db_start < $subject_beg;
		if ($query_strand==1){$query_strand='+'}elsif ($query_strand==-1){$query_strand='-'}
        my $query_lng = length($query_seq);
        my $query_freq = 0;
        if ($query =~ /_x(\d+)$/){
            $query_freq = $1;
        }else{
            die "Error occured when calculating the freq of query $query: $!\n";
            exit;
        }
		unless ($query_strand eq $subject_strand){
			$reversed_cnt += $query_freq;
			next;
		}
		$forward_cnt += $query_freq;
	}
	$total_cnt = $forward_cnt + $reversed_cnt;
	if ($total_cnt == 0){
		return (0);
	}else{
		my $strandness = sprintf ("%.2f", $forward_cnt / $total_cnt);
		return $strandness;
	}

}

sub get_subject_length{
	
	my ($subjects) = @_;
	for my $subject (@{$subjects}){
		my $length = $database->length($subject);
		$hash_len{$subject} = $length;
	}
	return;
}

sub overlap{

	my ($beg1,$end1,$beg2,$end2) = @_;
	
	unless($beg1<$end1 and $beg2<$end2){

		print STDERR "beg can not be larger than end\n";

		exit;
	
	}
	
	if( ($beg1<=$beg2 and $beg2<$end1) or ($beg1<$end2 and $end2<=$end1)){

		return 1;
	
	}else{

		return 0;

	}

}

sub var_check{
	if (exists $opts{v}){
        die "\nexcise_potential_precursors.pl version = $version\n\n";
    }
    if (exists $opts{h}){
        &var_error();
    }
    if($opts{b}){
		$bam_file = $opts{b};
	}else{
		&var_error();
	}
    if($opts{l}){
        $max_len = $opts{l};
    }else{
        $max_len = 277;
    }
	if($opts{m}){
		$min_freq = $opts{m}; 
	}else{
		$min_freq = 15;
	}
	if($opts{r}){
		$reference = $opts{r}; 
	}else{
		&var_error();
	}
	if($opts{o}){
		$out_file = $opts{o};
	}else{
		&var_error();
	}
    return;
}

sub var_error {
    print "\n";
    print " excise_potential_precursors.pl excises potential miRNA precursors with high\n";
    print " throughput sequencing.\n\n";
    print " Version = $version\n\n";
    print " WARNING: You did not provide enough information!\n\n" unless exists $opts{h};
    print " Usage: excise_potential_precursors.pl -b <reads_vs_genome_sorted.bam> -r <genome.fa> \\\n";
    print "                                       -o <precursors.fa> [-m <min_freq>]\n";
    print " REQUIRED:\n";
    print " -b <sorted.bam>    profile of reads mapped to the related genome in a formatted, sorted and\n";
    print "                    indexed BAM file\n";
    print " -r <genome.fa>     a related genome in multi-fasta format\n";
    print " -o <precursors.fa> an output file in multi-fasta format\n";
    print "\n";
    print " OPTIONAL:\n";
    print " -l <int>           maxmium length of the potential precursor to be excised (default: 277)\n";
	print " -m <int>           minimum frequency of reads to trigger a precursor excising\n";
    print "\n";
    exit 1;
}

#######################################################################################################
