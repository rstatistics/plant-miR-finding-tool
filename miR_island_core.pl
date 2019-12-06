#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;
use Bio::DB::Fasta;
use Bio::DB::Sam;


################################## MIR-ISLAND #################################################

################################### INPUT ##################################################

my ($sorted_bam_file,$file_struct,$genome,$cotyledon,$file_mature,$result_out,$table_out,$col1_width,$coords_file);
my $version = "1.0";
my %options = ();
getopts("a:b:c:g:hl:m:r:s:t:v", \%options);

&var_check();

#############################################################################################

############################# GLOBAL VARIABLES ##############################################


#parameters
my $seed_lng=11; #for plant criteria (2..12)

#for plant criteria (20..24)
my $mature_lng_min=20;
my $mature_lng_max=24;
my $star_lng_max=24;

my $score_star=3.9;
my $score_star_not=-1.3;
my $score_seed=7.635;
my $score_seed_not=-1.17;
my @scores_stem=(-3.1,-2.3,-2.2,-1.6,-1.5,0.1,0.6,0.8,0.9,0.9,0);
my $score_min=1;

my $e=2.718281828;

#hashes
my %hash_desc;
my %hash_seq;
my %hash_struct;
my %hash_mfe;
my %hash_seeds;
my %hash_mirs;
my %hash_query;
my %hash_comp;
my %hash_bp;
my %hash_stars;
my %hash_stored;
my %hash_list;
my %hash_site;
my %hash_dead;
my %hash_value;
my %hash_figure;
#other variables
my $param_a;
my $param_b;
my $param_c;
my $flag = 0;
my $db_old;
my $lines;
my $out_of_bound;

##############################################################################################

################################  MAIN  ######################################################

#different cotyledon type has different parametres
if($cotyledon eq 'monocot'){
    $param_a = 1.339e-12;
    $param_b = 2.77826e-13;
    $param_c = 45.8426;
}elsif($cotyledon eq 'dicot'){
    $param_a = 4.46e-4;
    $param_b = 9.125e-5;
    $param_c = 26.9293;
}else{
    die "Option 'c' $cotyledon is not recognised\nOnly 'monocot' or 'dicot' is allowed.\n";
}

#if conservation is scored, the fasta file of known miRNA sequences is parsed
if($file_mature){create_hash_seeds($file_mature)};

#parse signature file in bam format and resolve each potential precursor
parse_signature_and_structure($file_struct,$sorted_bam_file,$genome);

#remove false result
error_checking();

#output detected miRNAs
print_results();

exit 0;

##############################################################################################

############################## SUBROUTINES ###################################################

sub parse_signature_and_structure{
    #parses the output from RNAfold and reads it into hashes
	
    my($struct_file,$bam_file,$reference) = @_;
	my $database = Bio::DB::Fasta->new($reference);
	my $sam = Bio::DB::Sam->new(-bam => $bam_file, 
					            -fasta_file => $reference);
	my @targets = $sam->seq_ids;
	
	##>SL2.40ch00_23 strand:+ excise_beg:12897 excise_end:12979
	##UGUUAAUGCAUUUUAUAUGUAUUGAACGGGACCAAGUAGAGAUAGUUGUUUUAGACUGACUAACAAAAACAUAUUCUCUAAAC
	##...((((((((.....))))))))((((((((....((((((....((((((.............))))))...))))))... (-27.54)
    ## $subject=SL2.40ch00_23; $db=SL2.40ch00; $strand='+'; $beg=12897; $end=12979; $mfe=-27.54;
	my ($subject_id, $primary_id, $subject_strand, $subject_beg);
	my ($subject_end, $subject_seq, $subject_struct, $subject_mfe);
    open (FILE_STRUCT, "<$struct_file") or die "can not open file_struct\n";
    while (<FILE_STRUCT>){
        chomp;
        if (/^>((\S+)_\S+)\s+strand:(\S)\s+excise_beg:(\d+)\s+excise_end:(\d+)(.*)/){
			$subject_id    		= $1;
			$primary_id	 		= $2;
            $subject_strand     = $3;
            $subject_beg        = $4;
            $subject_end        = $5;
			$subject_seq        = "";
			$subject_struct     = "";
			$subject_mfe        = 0;
			while (<FILE_STRUCT>){
                chomp;
                if (/^>((\S+)_\S+)\s+strand:(\S)\s+excise_beg:(\d+)\s+excise_end:(\d+)(.*)/){
					$db_old = $subject_id;
					my @alignments = $sam->get_features_by_location(-seq_id => $primary_id,
																	-start => $subject_beg,
																	-end => $subject_end);
					for my $line (@alignments){
						my $query    		= $line->qname;
						my $query_seq 		= $line->query->dna;
						my $db_start		= $line->start;
						my $db_stop			= $line->end;
						my $query_strand	= $line->strand;
                        # samtools' output only stands for reads' overlaps the range [start,end]
                        next if $db_stop > $subject_end || $db_start < $subject_beg;
						if ($query_strand==1){$query_strand='+'}elsif ($query_strand==-1){$query_strand='-'}
						unless ($query_strand eq $subject_strand){
                            my $inverse_strand_query = $query;
                            $inverse_strand_query =~ s/.*_x//g;
                            my $inverse_strand_count += $inverse_strand_query;
                            next;
                        }
                        # if query on minus strand, revcom it
						if ($query_strand eq '-'){ $query_seq = revcom($query_seq);}
						my $query_tmp = $query;
						$query_tmp =~ s/.*_x//g;
						my $query_freq = $query_tmp;

                        ## Not allow same read map to precursor with multiple sites
                        if (exists $hash_query{$query} && $flag==0){
							$flag = 1;
                            #print STDERR "$db_old\tSame read maps to multiple sites!\n";
                        }
						#read information of the query (deep sequence) into hash
                        $hash_query{$query}{"db_start"}=$db_start;
						$hash_query{$query}{"db_stop"}=$db_stop;
						if ($subject_strand eq '+'){
							$hash_query{$query}{"db_beg"}=$db_start-$subject_beg+1;
							$hash_query{$query}{"db_end"}=$db_stop-$subject_beg+1;
						}elsif($subject_strand eq '-'){
							$hash_query{$query}{"db_beg"}=$subject_end-$db_stop+1;
							$hash_query{$query}{"db_end"}=$subject_end-$db_start+1;
						}
						$hash_query{$query}{"strand"}='+';
						$hash_query{$query}{"freq"}=$query_freq;
						# translate Upper letter to lower
						$query_seq =~ tr/ATUCG/attcg/;
						$hash_query{$query}{"seq"}=$query_seq;
						#save the signature information
                        my $tmp_line = "$query\t$hash_query{$query}{'db_beg'}\t$hash_query{$query}{'db_end'}\n";
                        $lines .= $tmp_line;
					}

					$hash_seq{$subject_id}		= $subject_seq;
					$hash_struct{$subject_id}	= $subject_struct;
					$hash_mfe{$subject_id}		= $subject_mfe;
					$hash_comp{'subject_strand'}= $subject_strand;
					$hash_comp{'primary_id'}	= $primary_id;
					$hash_comp{'subject_beg'}	= $subject_beg;
					$hash_comp{'subject_end'}	= $subject_end;
					resolve_potential_precursor();

					$subject_id    		= $1;
					$primary_id    		= $2;
					$subject_strand     = $3;
					$subject_beg        = $4;
					$subject_end        = $5;
					$subject_seq        = "";
					$subject_struct     = "";
					$subject_mfe        = 0;

					next;
                }
				if(/^\w/){
					tr/aucgU/ATCGT/;
					$subject_seq .= $_;
				}
				if(/((\.|\(|\))+)/){
					$subject_struct .=$1;
				}
				if(/\((\s*-\d+\.\d+)\)/){
					$subject_mfe = $1;
				}
			}
        }
	}
	$db_old = $subject_id;
	my @alignments = $sam->get_features_by_location(-seq_id => $primary_id,
													-start => $subject_beg,
													-end => $subject_end);
	for my $line (@alignments){
		my $query    		= $line->qname;
		my $query_seq 		= $line->query->dna;
		my $db_start		= $line->start;
		my $db_stop			= $line->end;
		my $query_strand	= $line->strand;
        # samtools' output only stands for reads' overlaps the range [start,end]
        next if $db_stop > $subject_end || $db_start < $subject_beg;
		if ($query_strand==1){$query_strand='+'}elsif ($query_strand==-1){$query_strand='-'}
		unless ($query_strand eq $subject_strand){
            my $inverse_strand_query = $query;
            $inverse_strand_query =~ s/.*_x//g;
            my $inverse_strand_count += $inverse_strand_query;
            next;
        }
		# if query on minus strand, revcom it
		if ($query_strand eq '-'){ $query_seq = revcom($query_seq);}
		my $query_tmp = $query;
		$query_tmp =~ s/.*_x//g;
		my $query_freq = $query_tmp;

        ## Not allow same read map to precursor with multiple sites
        if (exists $hash_query{$query} && $flag==0){
			$flag = 1;
            #print STDERR "$db_old\tSame read maps to multiple sites!\n";

        }
		#read information of the query (deep sequence) into hash
		$hash_query{$query}{"db_start"}=$db_start;
		$hash_query{$query}{"db_stop"}=$db_stop;
		if ($subject_strand eq '+'){
			$hash_query{$query}{"db_beg"}=$db_start-$subject_beg+1;
			$hash_query{$query}{"db_end"}=$db_stop-$subject_beg+1;
		}elsif($subject_strand eq '-'){
			$hash_query{$query}{"db_beg"}=$subject_end-$db_stop+1;
			$hash_query{$query}{"db_end"}=$subject_end-$db_start+1;
		}
		$hash_query{$query}{"strand"}='+';
		$hash_query{$query}{"freq"}=$query_freq;
		# translate Upper letter to lower
		$query_seq =~ tr/ATUCG/attcg/;
		$hash_query{$query}{"seq"}=$query_seq;
		#save the signature information
        my $tmp_line = "$query\t$hash_query{$query}{'db_beg'}\t$hash_query{$query}{'db_end'}\n";
        $lines .= $tmp_line;
	}
	$hash_seq{$subject_id}		= $subject_seq;
	$hash_struct{$subject_id}	= $subject_struct;
	$hash_mfe{$subject_id}		= $subject_mfe;
	$hash_comp{'subject_strand'}= $subject_strand;
	$hash_comp{'primary_id'}	= $primary_id;
	$hash_comp{'subject_beg'}	= $subject_beg;
	$hash_comp{'subject_end'}	= $subject_end;
	resolve_potential_precursor();

    close FILE_STRUCT;
    return;
}

sub parse_file_arf{

#    read through the signature blastparsed file, fills up a hash with information on queries
#    (deep sequences) mapping to the current db (potential precursor) and resolve each
#    potential precursor in turn

    my($file)=@_;

    open(FILENAME, $file) or die "Could not open file $file";

    while (my $line = <FILENAME>){
	if($line=~/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){

	    my $query=$1;
	    my $query_lng=$2;
	    my $query_beg=$3;
	    my $query_end=$4;
	    my $query_seq=$5;
	    my $db=$6;
	    my $db_lng=$7;
	    my $db_beg=$8;
	    my $db_end=$9;
	    my $db_seq=$10;
	    my $strand=$11;
	    my $edits=$12;
	    my $edit_string=$13;

	    #only reads that map sense to the potential precursor are considered
	    if($strand eq "-"){next;}

	    #if the new line concerns a new db (potential precursor) then the old db must be resolved
	    if($db_old and $db_old ne $db){
		resolve_potential_precursor();
	    }

	    #resolve the number of reads that the deep sequence represents
	    my $freq=find_freq($query);

	    #read information of the query (deep sequence) into hash
	    $hash_query{$query}{"db_beg"}=$db_beg;
	    $hash_query{$query}{"db_end"}=$db_end;
	    $hash_query{$query}{"strand"}=$strand;
	    $hash_query{$query}{"freq"}=$freq;

	    #save the signature information
	    $lines.=$line;

	    $db_old=$db;
	}
    }
    resolve_potential_precursor();
}

sub resolve_potential_precursor{

#   dissects the potential precursor in parts by filling hashes, and tests if it passes the
#   initial filter and the scoring filter
#   binary variable whether the potential precursor is still viable
    
    # return if one read maps to potential precursor with multiple sites
    if ($flag){ reset_variables(); return; }

    my $ret=1;

    fill_structure();

    fill_pri();

    fill_mature();

    fill_star();

    fill_loop();

    fill_lower_flanks();

	#this is the actual classification
    unless(pass_filtering_initial() and pass_threshold_score()){$ret=0;}

    store_results($ret);

    reset_variables();

    return;

}

sub store_results{
    my ($ret)=@_;

    if($ret){
	    store_hash_comp();
	    return;
	}
}

sub pass_threshold_score{

#    this is the scoring

    #minimum free energy of the potential precursor
    my $score_mfe=score_mfe($hash_comp{"pri_mfe"},$hash_comp{"pri_end"});

    #count of reads that map in accordance with Dicer processing
    my $score_freq=score_freq($hash_comp{"freq_total"});

    #basic score
    my $score=$score_mfe+$score_freq;

    #scoring of conserved seed/seed (optional)
    if($file_mature){

	#if the seed is conserved
	if(test_seed_conservation()){

	    #seed from position 2-14
	    my $seed=substr($hash_comp{"mature_seq"},1,$seed_lng);

	    #resolve DNA/RNA ambiguities
	    $seed=~tr/[T]/[U]/;

	    #print score contribution
	    $hash_comp{"score_seed"}=$score_seed;

	    #add to score
	    $score+=$score_seed;

	#if the seed is not conserved
	}else{
	    #print (negative) score contribution
	    $hash_comp{"score_seed"}=$score_seed_not;

	    #add (negative) score contribution
	    $score+=$score_seed_not;
	}
    }

    #if the majority of potential star reads fall as expected from Dicer processing
    if($hash_comp{"star_read"}){
	$hash_comp{"score_star"}=$score_star;
	$score+=$score_star;
    }else{
	$hash_comp{"score_star"}=$score_star_not;
	$score+=$score_star_not;
    }

    #round off values to one decimal
    my $round_mfe=round($score_mfe*10)/10;
    my $round_freq=round($score_freq*10)/10;
    my $round=round($score*10)/10;

    #print scores
    $hash_comp{"score_mfe"}=$round_mfe;
    $hash_comp{"score_freq"}=$round_freq;
    $hash_comp{"score"}=$round;

    #return 1 if the potential precursor is accepted, return 0 if discarded
    unless($score>=$score_min){return 0;}
    return 1;
}

sub print_file{

    #print string to file

    my($file,$string)=@_;

    open(FILE, ">$file");
    print FILE "$string";
    close FILE;
}

sub test_seed_conservation{

    #test if seed is identical to seed from known miRNA, return 1 or 0

    my $seed=substr($hash_comp{"mature_seq"},1,$seed_lng);
    $seed=~tr/[T]/[U]/;
    if($hash_seeds{$seed}){

	my $seed_family=$hash_mirs{$seed};
	if($seed_family=~/^(\S+)/){
	    #example of a miRNA with the same seed
	    $hash_comp{"seed_family"}=$1;
	}

	return 1;
    }

    return 0;
}

sub pass_filtering_initial{

    #test if the structure forms a plausible hairpin
    unless(pass_filtering_structure()){$hash_comp{"problem_structure"}="The candidate has been discarded because the structure is inconsistent with Dicer processing:\n"; return 0;}

    #test if >90% of reads map to the hairpin in consistence with Dicer processing
    unless(pass_filtering_signature()){$hash_comp{"problem_signature"}="The candidate has been discarded because the signature is inconsistent with Dicer processing:\n";return 0;}

    return 1;

}

sub pass_filtering_signature{

    #if the putative mature sequence is longer than the designated number of nts, discard
    my $mature_lng=$hash_comp{"mature_end"}-$hash_comp{"mature_beg"}+1;
    my $star_lng=$hash_comp{"star_end_exp"}-$hash_comp{"star_beg_exp"}+1;
    if($mature_lng>$mature_lng_max){$hash_comp{"problem_mature_lng"}="the candidate mature seq is > $mature_lng_max nts\n"; return 0;}
    if($mature_lng<$mature_lng_min){$hash_comp{"problem_mature_lng"}="the candidate mature seq is < $mature_lng_min nts\n"; return 0;}
    if($star_lng>$star_lng_max){$hash_comp{"problem_star_lng"}="the candidate star seq is > $star_lng_max nts\n"; return 0;}
    #number of reads that map in consistence with Dicer processing
    my $consistent=0;

    #number of reads that map inconsistent with Dicer processing
    my $inconsistent=0;

#   number of potential star reads map in good consistence with Drosha/Dicer processing
#   (3' overhangs relative to mature product)
    my $star_perfect=0;

#   number of potential star reads that do not map in good consistence with 3' overhang
    my $star_fuzzy=0;

    my $freq_mature=0;
    my $freq_loop=0;
    my $freq_star=0;
    my $freq_span=0;


    #sort queries (deep sequences) by their position on the hairpin
    my @queries=sort {$hash_query{$a}{"db_beg"} <=> $hash_query{$b}{"db_beg"}} keys %hash_query;

    foreach my $query(@queries){

	#number of reads that the deep sequence represents
	unless(defined($hash_query{$query}{"freq"})){next;}
	my $query_freq=$hash_query{$query}{"freq"};

	#test which Dicer product (if any) the deep sequence corresponds to
	my $product=test_query($query);

	#and add the appropriate read counts
	if($product eq "mature"){$freq_mature+=$query_freq;}
	if($product eq "loop"){$freq_loop+=$query_freq;}
	if($product eq "star"){$freq_star+=$query_freq;}
    if($product eq "span"){$freq_span+=$query_freq;}
	#if the deep sequence corresponds to a Dicer product, add to the 'consistent' variable
	if(($product ne "others") && ($product ne "span")){$consistent+=$query_freq;}

	#if the deep sequence do not correspond to a Dicer product, add to the 'inconsistent' variable
	else{$inconsistent+=$query_freq;}

	#test a potential star sequence has good 3' overhang
	if($product eq "star"){
	    if(test_star($query)){$star_perfect+=$query_freq;}
	    else{$star_fuzzy+=$query_freq;}
	}
    }

#   if the majority of potential star sequences map in good accordance with 3' overhang
#    score for the presence of star evidence
    if($star_perfect>$star_fuzzy){$hash_comp{"star_read"}=1;}
    #find the most frequent star positions (this is opposed to star_expm which are the
    #positions that would be expected from the positions of the mature sequence and
    #the model of Dicer processing
    if(0<$star_perfect+$star_fuzzy){observed_star()};

    #number of reads that correspond to mature, loop, star and the sum of these
    $hash_comp{"freq_mature"}=$freq_mature;
    $hash_comp{"freq_loop"}=$freq_loop;
    $hash_comp{"freq_star"}=$freq_star;
    $hash_comp{"freq_total"}=$consistent;

    unless($consistent+$inconsistent>0){$hash_comp{"problem_read_cnt"}="no reads map to the candidate\n"; return 0;}
    #if precursor has no star supports with random read falls over star_span locus
    #this is not a miRNA
    if ($freq_star==0 && $freq_span>0){ return 0; }
    #unless >85% of the reads map in consistence with Dicer processing, the hairpin is discarded
    my $inconsistent_fraction=$inconsistent/($inconsistent+$consistent);
    unless($inconsistent_fraction<=0.15){$hash_comp{"problem_Dicer_signature"}="more than 15% of the reads map inconsistent with Dicer processing\n"; return 0;}
	#test if >70% of reads map to the hairpin in consistence with Dicer processing
	my $dicer_ratio = ($freq_mature+$freq_star)/($inconsistent+$consistent);
	unless($dicer_ratio>=0.70){$hash_comp{"problem_dicer_processing"}="more than 30% of the reads map inconsistent with dicer processing\n"; return 0;}
    #the hairpin is retained
    return 1;
}

sub test_star{

    #test if a deep sequence maps in good consistence with 3' overhang

    my $query=shift;

    #5' begin and 3' end positions
    my $beg=$hash_query{$query}{"db_beg"};
    my $end=$hash_query{$query}{"db_end"};
    my $freq=$hash_query{$query}{"freq"};

    $hash_stars{$beg}{$end}+=$freq;

    #the difference between observed and expected begin positions must be 0 or 1
    my $offset=$beg-$hash_comp{"star_beg_exp"};
    if($offset==0 or $offset==1 or $offset==-1){return 1;}

    return 0;
}

sub observed_star{

    my $beg_max;
    my $end_max;
    my $freq_max=0;

    my @begs=sort {$a<=>$b} keys %hash_stars;
    foreach my $beg(@begs){

	my @ends=sort {$b<=>$a} keys %{$hash_stars{$beg}};
	foreach my $end(@ends){

	    my $freq=$hash_stars{$beg}{$end};
	    if($freq_max<$freq){

		$beg_max=$beg;
		$end_max=$end;
		$freq_max=$freq;
	    }
	}
    }

    $hash_comp{"star_beg_obs"}=$beg_max;
    $hash_comp{"star_end_obs"}=$end_max;
    $hash_comp{"star_freq_consensus"}=$freq_max;

    return;
}

sub test_query{

    #test if deep sequence maps in consistence with Dicer processing

    my $query=shift;

    #begin, end, strand and read count
    my $beg=$hash_query{$query}{"db_beg"};
    my $end=$hash_query{$query}{"db_end"};
    my $strand=$hash_query{$query}{"strand"};
    my $freq=$hash_query{$query}{"freq"};
	my $mature_lng=$hash_comp{"mature_end"}-$hash_comp{"mature_beg"}+1;
    #should not be on the minus strand (although this has in fact anecdotally been observed for known miRNAs)
    if($strand eq '-'){return 0;}

    #the deep sequence is allowed to stretch 1 nt beyond the expected 5' end
    my $fuzz_beg=1;
    #the deep sequence is allowed to stretch 2 nt beyond the expected 3' end
    my $fuzz_end=2;
	
	#if mature sequence length equals 24, bring in strict criteria
	if ($mature_lng==24){
		$fuzz_beg=0;
		$fuzz_end=0;
	}
	
    #if in accordance with Dicer processing, return the type of Dicer product
    if(contained($beg,$end,$hash_comp{"mature_beg"}-$fuzz_beg,$hash_comp{"mature_end"}+$fuzz_end)){return "mature";}
    if(contained($beg,$end,$hash_comp{"star_beg_exp"}-$fuzz_beg,$hash_comp{"star_end_exp"}+$fuzz_end)){return "star";}
    if(contained($beg,$end,$hash_comp{"loop_beg"}-$fuzz_beg,$hash_comp{"loop_end"}+$fuzz_end)){return "loop";}
    if(between($hash_comp{"star_beg_exp"},$beg,$end) || between($hash_comp{"star_end_exp"},$beg,$end)){
        return "span";
    }
    #if not in accordance, return 0
    return "others";
}

sub between{
    my ($coord,$beg,$end) = @_;
    unless ($beg<$end){
        print STDERR "beg greater than end\n";
        return;
    }
    if (($beg<$coord-3) && ($end>$coord+3)){
        return 1;
    }else{
        return 0;
    }
}

sub pass_filtering_structure{

    #The potential precursor must form a hairpin with miRNA precursor-like characteristics

    #return value
    my $ret=1;

    #potential mature, star, loop and lower flank parts must be identifiable
    unless(test_components()){return 0;}

	#fill precurso into hash
	unless(fill_precursor()){return 0;}

    #minimum 60% base pairings in duplex
    #unless(bp_duplex()>=0.6){$ret=0; $hash_comp{"problem_duplex_bp"}="less than 60% of the nts in the candidate mature seq are base paired in the duplex\n";}
	
	## test if mature_seq has a low complex
    #unless (test_read_complex($hash_comp{"mature_seq"},6)){$ret=0; $hash_comp{"problem_mature_complex"}="mature sequence with low complex!\n"; }
	
    #maximum 5 mismatches in mature and star duplex
    unless(test_mismatch()){$ret=0; $hash_comp{"problem_duplex_mismatch"}="too many mismatches in mature and star duplex\n";}

    #not more than 3 nt difference between mature and star length
    unless(-3<=diff_lng() and diff_lng()<=3){$ret=0; $hash_comp{"problem_duplex_lng"}="the difference between the candidate mature and star sequence is > 3 nts\n";}
    
    #not less than 6 nt of the loop length
    unless(loop_lng()>=6){$ret=0; $hash_comp{"problem_loop_lng"}="loop length is less than 6 nts\n";}

    return $ret;
}

sub test_read_complex{
    my ($seq,$cnt) = @_;
    if ($seq =~ /A{$cnt}|T{$cnt}|C{$cnt}|G{$cnt}|U{$cnt}/i){
        return 0;
    }
    return 1;
}

sub test_mismatch{
    my $mature_struct=$hash_comp{'mature_struct'};
    my $star_struct=$hash_comp{'star_struct'};
    my %hash_duplex;
    my $struct = "";
    my $lng = 0;
    if ($mature_struct =~ /^(\(|\.)+$/ && $star_struct =~ /^(\)|\.)+$/){
        chop($mature_struct);
        chop($mature_struct);
        chop($star_struct);
        chop($star_struct);
        my $mature_brack_cnt = $mature_struct =~ tr/\(/\(/;
        my $star_brack_cnt = $star_struct =~ tr/\)/\)/;
        return 0 unless $mature_brack_cnt == $star_brack_cnt;
        $lng = length $mature_struct;
        $struct = $mature_struct . $star_struct;
    }elsif ($mature_struct =~ /^(\)|\.)+$/ && $star_struct =~ /^(\(|\.)+$/){
        chop($mature_struct);
        chop($mature_struct);
        chop($star_struct);
        chop($star_struct);
        my $mature_brack_cnt = $mature_struct =~ tr/\)/\)/;
        my $star_brack_cnt = $star_struct =~ tr/\(/\(/;
        return 0 unless $mature_brack_cnt == $star_brack_cnt;
        $lng = length $star_struct;
        $struct = $star_struct . $mature_struct;
    }
    my @bps;
    my $struct_lng = length $struct;
    for (my $pos=1; $pos<=$struct_lng; $pos++){
        my $struct_pos = excise_struct($struct,$pos,$pos,"+");
        if ($struct_pos eq "("){
            push @bps, $pos;
        }
        if ($struct_pos eq ")"){
            my $pos_prev = pop @bps;
            print STDERR "struct:\t",$struct unless $pos_prev;
            $hash_duplex{$pos_prev} = $pos;
        }
    }
    my $left_offset = 0;
    my $right_offset = 0;
    my $left_pos = 0;
    my $right_pos = $struct_lng + 1;
    my $mismatch = 0;
    my $mismatch_cnt = 0;
    my $bulge = 0;
    my $bulge_cnt = 0;
    for my $i (1 .. $lng){
        unless (defined $hash_duplex{$i}){
        }else{
            $right_offset = $right_pos - $hash_duplex{$i} - 1;
            $left_offset = $i - $left_pos - 1;
            $right_pos = $hash_duplex{$i};
            $left_pos = $i;
            if ($right_offset + $left_offset == 0){
            }elsif ($right_offset * $left_offset == 0){
                $bulge += ($right_offset + $left_offset);
                $bulge_cnt ++;
            }elsif ($right_offset + $left_offset == 2){
                $mismatch += 1;
                $mismatch_cnt ++;
            }elsif ($right_offset + $left_offset == 3){
                $mismatch += 1;
                $mismatch_cnt ++;
                $bulge += 1;
                $bulge_cnt ++;
            }elsif ($right_offset + $left_offset == 4){
                if ($right_offset * $left_offset < 4){
                    $bulge += 2;
                    $bulge_cnt ++;
                    $mismatch += 1;
                    $mismatch_cnt ++;
                }else{
                    $mismatch += 2;
                    $mismatch_cnt ++;
                }
            }elsif ($right_offset + $left_offset == 5){
                if ($right_offset * $left_offset < 5){
                    $bulge += 3;
                    $bulge_cnt ++;
                    $mismatch += 1;
                    $mismatch_cnt ++;
                }else{
                    $bulge += 1;
                    $bulge_cnt ++;
                    $mismatch += 2;
                    $mismatch_cnt ++;
                }
            }elsif ($right_offset + $left_offset == 6){
                if ($right_offset * $left_offset < 6){
                    $bulge += 4;
                    $bulge_cnt ++;
                    $mismatch += 1;
                    $mismatch_cnt ++;
                }elsif ($right_offset * $left_offset > 8){
                    $mismatch += 3;
                    $mismatch_cnt ++;
                }else{
                    $bulge += 2;
                    $bulge_cnt ++;
                    $mismatch += 2;
                    $mismatch_cnt ++;
                }
            }else{
                return 0;
            }
            
        }
    
    }
    #for the last loop if last is mismatched
    unless ($left_pos == $lng){
    $left_offset = $lng - $left_pos;
    $right_offset = $right_pos - ($lng + 1);
    if ($right_offset + $left_offset == 0){
    }elsif ($right_offset * $left_offset == 0){
        $bulge += ($right_offset + $left_offset);
        $bulge_cnt ++;
    }elsif ($right_offset + $left_offset == 2){
        $mismatch += 1;
        $mismatch_cnt ++;
    }elsif ($right_offset + $left_offset == 3){
        $mismatch += 1;
        $mismatch_cnt ++;
        $bulge += 1;
        $bulge_cnt ++;
    }elsif ($right_offset + $left_offset == 4){
        if ($right_offset * $left_offset < 4){
            $bulge += 2;
            $bulge_cnt ++;
            $mismatch += 1;
            $mismatch_cnt ++;
        }else{
            $mismatch += 2;
            $mismatch_cnt ++;
        }
    }elsif ($right_offset + $left_offset == 5){
        if ($right_offset * $left_offset < 5){
            $bulge += 3;
            $bulge_cnt ++;
            $mismatch += 1;
            $mismatch_cnt ++;
        }else{
            $bulge += 1;
            $bulge_cnt ++;
            $mismatch += 2;
            $mismatch_cnt ++;
        }
    }elsif ($right_offset + $left_offset == 6){
        if ($right_offset * $left_offset < 6){
            $bulge += 4;
            $bulge_cnt ++;
            $mismatch += 1;
            $mismatch_cnt ++;
        }elsif ($right_offset * $left_offset > 8){
            $mismatch += 3;
            $mismatch_cnt ++;
        }else{
            $bulge += 2;
            $bulge_cnt ++;
            $mismatch += 2;
            $mismatch_cnt ++;
        }
    }else{
        return 0;
    }
    }
    if ($bulge_cnt == 0 && $mismatch <= 5){
        return 1;
    }elsif ($bulge_cnt == 1 && $bulge <=2 && $bulge + $mismatch <=4){
        return 1;
    }elsif ($bulge_cnt == 2 && $bulge <= 2 && $mismatch <= 1){
        return 1;
    }else{
        return 0;
    }
}

sub loop_lng{
    my $loop_beg=$hash_comp{'loop_beg'};
    my $loop_end=$hash_comp{'loop_end'};
    unless($loop_beg <= $loop_end){
        return 0;
    }
    my $loop_lng = $loop_end - $loop_beg + 1;

    return $loop_lng;
}

sub test_components{

    #tests whether potential mature, star, loop and lower flank parts are identifiable

    unless($hash_comp{"mature_struct"}){
	$hash_comp{"problem_struct_mature"}="no mature sequence could be identified\n";
	return 0;
    }

    unless($hash_comp{"star_struct"}){
	$hash_comp{"problem_struct_star"}="the candidate mature sequence does not form part of a likely miRNA duplex\n";
	return 0;
    }

    unless($hash_comp{"loop_struct"}){
	$hash_comp{"problem_struct_loop"}="no loop sequence could be identified\n";
   	return 0;
    }

    unless(defined $hash_comp{"flank_first_struct"}){
    	$hash_comp{"problem_struct_us_flanks"}="no upstream flanking sequence could be identified\n";
        return 0;
    }

    unless(defined $hash_comp{"flank_second_struct"}){
        $hash_comp{"problem_struct_ds_flanks"}="no downstream flanking sequence could be identified\n";
        return 0;
    }
    return 1;
}

sub bp_duplex{

    #fraction of nts of the mature sequence that are base paired in the duplex

    my $lng=$hash_comp{"mature_end"}-$hash_comp{"mature_beg"}+1;

    my $duplex_bps=0;
    my $mature_struct=$hash_comp{"mature_struct"};

    #simple pattern matching
    while($mature_struct=~/(\(|\))/g){
	$duplex_bps++;
    }
    return ($duplex_bps/$lng);
}

sub diff_lng{

    #find difference between mature and star lengths

    my $mature_lng=length $hash_comp{"mature_struct"};
    my $star_lng=length $hash_comp{"star_struct"};
    my $diff_lng=$mature_lng-$star_lng;
    return $diff_lng;
}

sub fill_precursor{
    #assembles the potential precursor sequence and structure from the expected Dicer products
    #this is the expected biological precursor, in contrast with 'subject_seq' that includes
    #some genomic flanks on both sides

    my $pre_struct;
    my $pre_seq;
	my $key_start;
	my $key_stop;
	my $offset;
	#key_start means the left five prime of "M-L-S" or "S-L-M"
	#key_stop means the right three prime of "M-L-S" or "S-L-M"
	my $key_lng = length($hash_comp{"mature_seq"})+length($hash_comp{"loop_seq"})+length($hash_comp{"star_seq"});
    if(defined $hash_comp{"flank_first_struct"} and $hash_comp{"mature_struct"} and $hash_comp{"loop_struct"} and $hash_comp{"star_struct"} and defined $hash_comp{"flank_second_struct"}){
	if($hash_comp{"mature_arm"} eq "first"){
		$pre_struct.=$hash_comp{"flank_first_struct"}.$hash_comp{"mature_struct"}.$hash_comp{"loop_struct"}.$hash_comp{"star_struct"}.$hash_comp{"flank_second_struct"};
		$pre_seq.=$hash_comp{"flank_first_seq"}.$hash_comp{"mature_seq"}.$hash_comp{"loop_seq"}.$hash_comp{"star_seq"}.$hash_comp{"flank_second_seq"};
		$offset = length($hash_comp{"mature_seq"})+length($hash_comp{"loop_seq"})+length($hash_comp{"star_seq"});
		if($hash_comp{'subject_strand'} eq '+'){
			$key_start = $hash_comp{"mature_start"};
			$key_stop = $hash_comp{"mature_start"} + $offset - 1;
		}elsif($hash_comp{'subject_strand'} eq '-'){
            $key_start = $hash_comp{"mature_stop"} - $offset + 1;
            $key_stop = $hash_comp{"mature_stop"};
		}
	}else{
		$pre_struct.=$hash_comp{"flank_first_struct"}.$hash_comp{"star_struct"}.$hash_comp{"loop_struct"}.$hash_comp{"mature_struct"}.$hash_comp{"flank_second_struct"};
		$pre_seq.=$hash_comp{"flank_first_seq"}.$hash_comp{"star_seq"}.$hash_comp{"loop_seq"}.$hash_comp{"mature_seq"}.$hash_comp{"flank_second_seq"};
		$offset = length($hash_comp{"loop_seq"})+length($hash_comp{"star_seq"})+length($hash_comp{"mature_seq"});
		if ($hash_comp{'subject_strand'} eq '+'){
			$key_start = $hash_comp{"mature_stop"} - $offset + 1;
			$key_stop = $hash_comp{"mature_stop"};
		}elsif ($hash_comp{'subject_strand'} eq '-'){
			$key_start = $hash_comp{"mature_start"};
            $key_stop = $hash_comp{"mature_start"} + $offset - 1;
		}
	}

    #read into hash
    $hash_comp{"pre_struct"}=$pre_struct;
    $hash_comp{"pre_seq"}=$pre_seq;
	$hash_comp{"key_start"}=$key_start;
	$hash_comp{"key_stop"}=$key_stop;
	
	return 1;
	
	}else{
		return 0;
	}

}

sub do_test_assemble{

#    not currently used, tests if the 'pri_struct' as assembled from the parts (Dicer products, lower flanks)
#    is identical to 'pri_struct' before disassembly into parts

    my $assemble_struct;
	my $key_start;
	my $key_stop;
	my $key_lng = length($hash_comp{"mature_seq"})+length($hash_comp{"loop_seq"})+length($hash_comp{"star_seq"});
    if(defined $hash_comp{"flank_first_struct"} and $hash_comp{"mature_struct"} and $hash_comp{"loop_struct"} and $hash_comp{"star_struct"} and defined $hash_comp{"flank_second_struct"}){
	if($hash_comp{"mature_arm"} eq "first"){
	    $assemble_struct.=$hash_comp{"flank_first_struct"}.$hash_comp{"mature_struct"}.$hash_comp{"loop_struct"}.$hash_comp{"star_struct"}.$hash_comp{"flank_second_struct"};
		if($hash_comp{'subject_strand'} eq '+'){
			$key_start = $hash_comp{"mature_beg"};
			$key_stop = $hash_comp{"mature_start"} + $key_lng - 1;
		}elsif($hash_comp{'subject_strand'} eq '-'){
			$key_start = $hash_comp{"mature_start"} - $key_lng + 1;
			$key_stop = $hash_comp{"mature_start"};
		}
	}else{
	    $assemble_struct.=$hash_comp{"flank_first_struct"}.$hash_comp{"star_struct"}.$hash_comp{"loop_struct"}.$hash_comp{"mature_struct"}.$hash_comp{"flank_second_struct"};
		if ($hash_comp{'subject_strand'} eq '+'){
			$key_start = $hash_comp{"mature_stop"} - $key_lng + 1;
			$key_stop = $hash_comp{"mature_stop"};
		}elsif ($hash_comp{'subject_strand'} eq '-'){
			$key_start = $hash_comp{"mature_stop"};
			$key_stop = $hash_comp{"mature_stop"} + $key_lng - 1;
		}	
	}
	unless($assemble_struct eq $hash_comp{"pri_struct"}){
	    $hash_comp{"test_assemble"}=$assemble_struct;
	    print_hash_comp();
	}
	$hash_comp{"key_start"}=$key_start;
	$hash_comp{"key_stop"}=$key_stop;
    }
    return;
 }

sub fill_structure{

    #reads the dot bracket structure into the 'bp' hash where each key and value are basepaired

    my $struct=$hash_struct{$db_old};
    my $lng=length $struct;

    #local stack for keeping track of basepairings
    my @bps;

    for(my $pos=1;$pos<=$lng;$pos++){
	my $struct_pos=excise_struct($struct,$pos,$pos,"+");

	if($struct_pos eq "("){
	    push(@bps,$pos);
	}

	if($struct_pos eq ")"){
	    my $pos_prev=pop(@bps);
	    $hash_bp{$pos_prev}=$pos;
	    $hash_bp{$pos}=$pos_prev;
	}
    }
    return;
}

sub fill_star{

    #fills specifics on the expected star strand into 'comp' hash ('component' hash)

    #if the mature sequence is not plausible, don't look for the star arm
    my $mature_arm=$hash_comp{"mature_arm"};
    unless($mature_arm){$hash_comp{"star_arm"}=0; return;}

    #if the star sequence is not plausible, don't fill into the hash
    my($star_beg_exp,$star_end_exp)=find_star();
    my $star_arm=arm_star($star_beg_exp,$star_end_exp);
    unless($star_arm){return;}

    #excise expected star sequence and structure
    my $star_seq=excise_seq($hash_comp{"pri_seq"},$star_beg_exp,$star_end_exp,"+");
    my $star_struct=excise_seq($hash_comp{"pri_struct"},$star_beg_exp,$star_end_exp,"+");

    #fill into hash
    $hash_comp{"star_beg_exp"}=$star_beg_exp;
    $hash_comp{"star_end_exp"}=$star_end_exp;
    $hash_comp{"star_seq"}=$star_seq;
    $hash_comp{"star_struct"}=$star_struct;
    $hash_comp{"star_arm"}=$star_arm;

    return;
}

sub find_star{

    #uses the 'bp' hash to find the expected star begin and end positions from the mature positions

    #the -2 is for the overhang
    my $mature_beg=$hash_comp{"mature_beg"};
    my $mature_end=$hash_comp{"mature_end"}-2;
    my $mature_lng=$mature_end-$mature_beg+1;

    #in some cases, the last nucleotide of the mature sequence does not form a base pair,
    #and therefore does not basepair with the first nucleotide of the star sequence.
    #In this case, the algorithm searches for the last nucleotide of the mature sequence
    #to form a base pair. The offset is the number of nucleotides searched through.
    my $offset_star_beg=0;
    my $offset_beg=0;

    #the offset should not be longer than the length of the mature sequence, then it
    #means that the mature sequence does not form any base pairs
    while(not $offset_star_beg and $offset_beg<$mature_lng){
	if($hash_bp{$mature_end-$offset_beg}){
	    $offset_star_beg=$hash_bp{$mature_end-$offset_beg};
	}else{
	    $offset_beg++;
	}
    }
    #when defining the beginning of the star sequence, compensate for the offset
    my $star_beg_exp=$offset_star_beg-$offset_beg;

    #same as above
    my $offset_star_end=0;
    my $offset_end=0;
    while(not $offset_star_end and $offset_end<$mature_lng){
	if($hash_bp{$mature_beg+$offset_end}){
	    $offset_star_end=$hash_bp{$mature_beg+$offset_end};
	}else{
	    $offset_end++;
	}
    }
    #the +2 is for the overhang
    my $star_end_exp=$offset_star_end+$offset_end+2;

    return($star_beg_exp,$star_end_exp);
}

sub fill_pri{

    #fills basic specifics on the precursor into the 'comp' hash

    my $seq=$hash_seq{$db_old};
    my $struct=$hash_struct{$db_old};
    my $mfe=$hash_mfe{$db_old};
    my $length=length $seq;

    $hash_comp{"pri_id"}=$db_old;
    $hash_comp{"pri_seq"}=$seq;
    $hash_comp{"pri_struct"}=$struct;
    $hash_comp{"pri_mfe"}=$mfe;
    $hash_comp{"pri_beg"}=1;
    $hash_comp{"pri_end"}=$length;

    return;
}

sub fill_mature{

    #fills specifics on the mature sequence into the 'comp' hash

    my $mature_query=find_mature_query();
    my($mature_beg,$mature_end,$mature_start,$mature_stop)=find_positions_query($mature_query);
    my $mature_strand=find_strand_query($mature_query);
    my $mature_seq=excise_seq($hash_comp{"pri_seq"},$mature_beg,$mature_end,$mature_strand);
    my $mature_struct=excise_struct($hash_comp{"pri_struct"},$mature_beg,$mature_end,$mature_strand);
    my $mature_arm=arm_mature($mature_beg,$mature_end,$mature_strand);

    $hash_comp{"mature_query"}=$mature_query;
    $hash_comp{"mature_beg"}=$mature_beg;
    $hash_comp{"mature_end"}=$mature_end;
	$hash_comp{"mature_start"}=$mature_start;
	$hash_comp{"mature_stop"}=$mature_stop;
    $hash_comp{"mature_strand"}=$mature_strand;
    $hash_comp{"mature_struct"}=$mature_struct;
    $hash_comp{"mature_seq"}=$mature_seq;
    $hash_comp{"mature_arm"}=$mature_arm;

    return;
}

sub fill_loop{

    #fills specifics on the loop sequence into the 'comp' hash

    #unless both mature and star sequences are plausible, do not look for the loop
    unless($hash_comp{"mature_arm"} and $hash_comp{"star_arm"}){return;}

    my $loop_beg;
    my $loop_end;

    #defining the begin and end positions of the loop from the mature and star positions
    #excision depends on whether the mature or star sequence is 5' of the loop ('first')
    if($hash_comp{"mature_arm"} eq "first"){
	$loop_beg=$hash_comp{"mature_end"}+1;
    }else{
	$loop_end=$hash_comp{"mature_beg"}-1;
    }

    if($hash_comp{"star_arm"} eq "first"){
	$loop_beg=$hash_comp{"star_end_exp"}+1;
    }else{
	$loop_end=$hash_comp{"star_beg_exp"}-1;
    }

    #unless the positions are plausible, do not fill into hash
    unless(test_loop($loop_beg,$loop_end)){return;}

    my $loop_seq=excise_seq($hash_comp{"pri_seq"},$loop_beg,$loop_end,"+");
    my $loop_struct=excise_struct($hash_comp{"pri_struct"},$loop_beg,$loop_end,"+");

    $hash_comp{"loop_beg"}=$loop_beg;
    $hash_comp{"loop_end"}=$loop_end;
    $hash_comp{"loop_seq"}=$loop_seq;
    $hash_comp{"loop_struct"}=$loop_struct;

    return;
}

sub fill_lower_flanks{

    #fills specifics on the lower flanks and unpaired strands into the 'comp' hash

    #unless both mature and star sequences are plausible, do not look for the flanks
    unless($hash_comp{"mature_arm"} and $hash_comp{"star_arm"}){return;}

    my $flank_first_end;
    my $flank_second_beg;

    #defining the begin and end positions of the flanks from the mature and star positions
    #excision depends on whether the mature or star sequence is 5' in the potenitial precursor ('first')
    if($hash_comp{"mature_arm"} eq "first"){
	$flank_first_end=$hash_comp{"mature_beg"}-1;
    }else{
	$flank_second_beg=$hash_comp{"mature_end"}+1;
    }

    if($hash_comp{"star_arm"} eq "first"){
	$flank_first_end=$hash_comp{"star_beg_exp"}-1;
    }else{
	$flank_second_beg=$hash_comp{"star_end_exp"}+1;
    }

    #unless the positions are plausible, do not fill into hash
    unless(test_flanks($flank_first_end,$flank_second_beg)){return;}

    $hash_comp{"flank_first_end"}=$flank_first_end;
    $hash_comp{"flank_second_beg"}=$flank_second_beg;
    unless ($hash_comp{"pri_beg"}==$hash_comp{"flank_first_end"}+1){
        $hash_comp{"flank_first_seq"}=excise_seq($hash_comp{"pri_seq"},$hash_comp{"pri_beg"},$hash_comp{"flank_first_end"},"+");
        $hash_comp{"flank_first_struct"}=excise_struct($hash_comp{"pri_struct"},$hash_comp{"pri_beg"},$hash_comp{"flank_first_end"},"+");
    }else{
        $hash_comp{"flank_first_seq"}="";
        $hash_comp{"flank_first_struct"}="";
    }
    unless ($hash_comp{"pri_end"}==$hash_comp{"flank_second_beg"}-1){
        $hash_comp{"flank_second_seq"}=excise_seq($hash_comp{"pri_seq"},$hash_comp{"flank_second_beg"},$hash_comp{"pri_end"},"+");
        $hash_comp{"flank_second_struct"}=excise_struct($hash_comp{"pri_struct"},$hash_comp{"flank_second_beg"},$hash_comp{"pri_end"},"+");
    }else{
        $hash_comp{"flank_second_seq"}="";
        $hash_comp{"flank_second_struct"}="";
    }

    return;
}

sub arm_mature{

    #tests whether the mature sequence is in the 5' ('first') or 3' ('second') arm of the potential precursor

    my ($beg,$end,$strand)=@_;

    #mature and star sequences should alway be on plus strand
    if($strand eq "-"){return 0;}

    #there should be no bifurcations and minimum one base pairing
    my $struct=excise_seq($hash_comp{"pri_struct"},$beg,$end,$strand);
    if(defined($struct) and $struct=~/^(\(|\.)+$/ and $struct=~/\(/){
	return "first";
    }elsif(defined($struct) and $struct=~/^(\)|\.)+$/ and $struct=~/\)/){
	return "second";
    }
    return 0;
}

sub arm_star{

    #tests whether the star sequence is in the 5' ('first') or 3' ('second') arm of the potential precursor

    my ($beg,$end)=@_;

    #unless the begin and end positions are plausible, test negative
    unless($beg>0 and $beg<=$hash_comp{"pri_end"} and $end>0 and $end<=$hash_comp{"pri_end"} and $beg<=$end){return 0;}

    #no overlap between the mature and the star sequence
    if($hash_comp{"mature_arm"} eq "first"){
	($hash_comp{"mature_end"}<$beg) or return 0;
    }elsif($hash_comp{"mature_arm"} eq "second"){
	($end<$hash_comp{"mature_beg"}) or return 0;
    }

    #there should be no bifurcations and minimum one base pairing
    my $struct=excise_seq($hash_comp{"pri_struct"},$beg,$end,"+");
    if($struct=~/^(\(|\.)+$/ and $struct=~/\(/){
	return "first";
    }elsif($struct=~/^(\)|\.)+$/ and $struct=~/\)/){
	return "second";
    }
    return 0;
}

sub test_loop{

    #tests the loop positions

    my ($beg,$end)=@_;

    #unless the begin and end positions are plausible, test negative
    unless($beg>0 and $beg<=$hash_comp{"pri_end"} and $end>0 and $end<=$hash_comp{"pri_end"} and $beg<=$end){return 0;}

    return 1;
}

sub test_flanks{

    #tests the positions of the lower flanks

    my ($beg,$end)=@_;

    #unless the begin and end positions are plausible, test negative
    unless($beg>0 and $end>0 and $beg<=$end){ return 0; }
    unless($beg<=$hash_comp{"pri_end"} and $end<=$hash_comp{"pri_end"}+1){return 0;}

    return 1;
}

sub comp{

    #subroutine to retrive from the 'comp' hash

    my $type=shift;
    my $component=$hash_comp{$type};
    return $component;
}

sub find_strand_query{

    #subroutine to find the strand for a given query

    my $query=shift;
    my $strand=$hash_query{$query}{"strand"};
    return $strand;
}

sub find_positions_query{

    #subroutine to find the begin and end positions for a given query

    my $query=shift;
    my $beg=$hash_query{$query}{"db_beg"};
    my $end=$hash_query{$query}{"db_end"};
	my $start=$hash_query{$query}{"db_start"};
	my $stop=$hash_query{$query}{"db_stop"};
    return ($beg,$end,$start,$stop);
}

sub find_mature_query{

    #finds the query with the highest frequency of reads and returns it
    #is used to determine the positions of the potential mature sequence

    my @queries=sort {$hash_query{$b}{"freq"} <=> $hash_query{$a}{"freq"}} keys %hash_query;
    my $mature_query=$queries[0];
    return $mature_query;
}

sub reset_variables{

    #resets the hashes for the next potential precursor

    %hash_query=();
    %hash_comp=();
    %hash_bp=();
    %hash_stars=();
	%hash_seq=();
	%hash_struct=();
	%hash_mfe=();
    $lines=();
    $flag=0;

    return;
}

sub excise_seq{

    #excise sub sequence from the potential precursor

    my($seq,$beg,$end,$strand)=@_;

    #begin can be equal to end if only one nucleotide is excised
    unless($beg<=$end){print STDERR "begin can not be smaller than end for $db_old\n";exit;}

    #rarely, permuted combinations of signature and structure cause out of bound excision errors.
    #this happens once appr. every two thousand combinations
    unless($beg<=length($seq)){$out_of_bound++;return 0;}

    #the blast parsed format is 1-indexed, substr is 0-indexed
    my $sub_seq=substr($seq,$beg-1,$end-$beg+1);

    #if on the minus strand, the reverse complement should be returned
    if($strand eq "-"){
	$sub_seq=revcom($sub_seq);
    }

    return $sub_seq;

}

sub excise_struct{

    #excise sub structure

    my($struct,$beg,$end,$strand)=@_;
    my $lng=length $struct;

    #begin can be equal to end if only one nucleotide is excised
    unless($beg<=$end){print STDERR "begin can not be smaller than end for $db_old\n";exit;}

    #rarely, permuted combinations of signature and structure cause out of bound excision errors.
    #this happens once appr. every two thousand combinations
    unless($beg<=length($struct)){return 0;}

    #the blast parsed format is 1-indexed, substr is 0-indexed
    my $sub_struct=substr($struct,$beg-1,$end-$beg+1);

    return $sub_struct;
}

sub create_hash_seeds{

    #parses a fasta file with sequences of known miRNAs considered for conservation purposes
    #reads the seeds into a hash

    my ($file) = @_;
    my ($id, $desc, $sequence, $seed) = ();

    open (FASTA, "<$file") or die "can not open $file\n";
    while (<FASTA>)
    {
        chomp;
        if (/^>(\S+)(.*)/)
	{
	    $id       = $1;
	    $desc     = $2;
	    $sequence = "";
	    $seed  = "";
	    while (<FASTA>){
                chomp;
                if (/^>(\S+)(.*)/){
		    $seed                = substr($sequence,1,$seed_lng);
		    $seed                =~ tr/[T]/[U]/;
		    $hash_mirs{$seed}   .="$id\t";
		    $hash_seeds{$seed} += 1;

		    $id               = $1;
		    $desc             = $2;
		    $sequence         = "";
		    $seed          = "";
		    next;
                }
		$sequence .= $_;
            }
        }
    }
    $seed                = substr($sequence,1,$seed_lng);
    $seed                =~ tr/[T]/[U]/;
    $hash_mirs{$seed}   .="$id\t";
    $hash_seeds{$seed} += 1;
    close FASTA;
}

sub parse_file_struct{

    #parses the output from RNAfold and reads it into hashes

    my($file) = @_;
    my($id,$desc,$seq,$struct,$mfe) = ();

    open (FILE_STRUCT, "<$file") or die "can not open $file\n";
    while (<FILE_STRUCT>)
    {
        chomp;
        if (/^>(\S+)\s*(.*)/)
	{
	    $id          = $1;
	    $desc        = $2;
	    $seq         = "";
	    $struct      = "";
	    $mfe         = "";
	    while (<FILE_STRUCT>){
                chomp;
                if (/^>(\S+)\s*(.*)/){
		    $hash_desc{$id}   = $desc;
		    $hash_seq{$id}    = $seq;
		    $hash_struct{$id} = $struct;
		    $hash_mfe{$id}    = $mfe;

		    $id          = $1;
		    $desc        = $2;
		    $seq         = "";
		    $struct      = "";
		    $mfe         = "";

		    next;
                }
		if(/^\w/){
		    tr/uU/tT/;
		    $seq .= $_;
		}if(/((\.|\(|\))+)/){
		    $struct .=$1;
		}
		if(/\((\s*-\d+\.\d+)\)/){
		    $mfe = $1;
		}

	    }
        }
    }

    $hash_desc{$id}        = $desc;
    $hash_seq{$id}         = $seq;
    $hash_struct{$id}      = $struct;
    $hash_mfe{$id}         = $mfe;

    close FILE_STRUCT;
    return;
}

sub precursor_decision{
	my ($x,$y) = @_;
#if ($x == $y){ $y = $y + 1; }
#	my $m = ($x+$y)/2;
#	my $var = ($x-$m)**2 + ($y-$m)**2;
#	my $value = $var / ( ($x-$y)**2 * $m );
#	return $value;
    my $value = $x + $y;
    return $value;
}

sub find_freq{

    #finds the frequency of a given read query from its id.

    my($query)=@_;

    if($query=~/_x(\d+)/){
	my $freq=$1;
	return $freq;
    }else{
	return 0;
    }
}

sub store_hash_comp{
	unless (defined $hash_value{$hash_comp{'primary_id'}}{$hash_comp{'subject_strand'}}{$hash_comp{'key_start'}}{$hash_comp{'key_stop'}}){
		$hash_value{$hash_comp{'primary_id'}}{$hash_comp{'subject_strand'}}{$hash_comp{'key_start'}}{$hash_comp{'key_stop'}} = 100;
	}
	my $flank1 = length($hash_comp{'flank_first_seq'});
	my $flank2 = length($hash_comp{'flank_second_seq'});
	my $value = precursor_decision($flank1,$flank2);

	my $freq_mature = $hash_comp{'freq_mature'};
    my $freq_star = $hash_comp{'freq_star'};
    my $mature_seq = $hash_comp{'mature_seq'};
    my $star_seq = $hash_comp{'star_seq'};
    my $mature_beg = $hash_comp{'mature_beg'};
    my $mature_end = $hash_comp{'mature_end'};
    my $loop_beg = $hash_comp{'loop_beg'};
    my $loop_end = $hash_comp{'loop_end'};
    my $star_beg = $hash_comp{'star_beg_exp'};
    my $star_end = $hash_comp{'star_end_exp'};
    my $mature_id = $hash_comp{'mature_query'};
    my $star_id = $hash_comp{'mature_query'}.'*';
	
    my $mature_site;
    my $star_site;
	my $key_start_freq;
	my $key_stop_freq;
    if ($mature_beg < $star_beg){
        $mature_site = '5p';
        $star_site = '3p';
		$key_start_freq = $freq_mature;
		$key_stop_freq = $freq_star;
    }else{
        $mature_site = '3p';
        $star_site = '5p';
		$key_start_freq = $freq_star;
		$key_stop_freq = $freq_mature;
    }
	unless (defined $hash_site{$hash_comp{'primary_id'}}{$hash_comp{'subject_strand'}}{$hash_comp{'key_start'}}{$hash_comp{'key_start'}}){
		$hash_site{$hash_comp{'primary_id'}}{$hash_comp{'subject_strand'}}{$hash_comp{'key_start'}}{$hash_comp{'key_start'}}=$key_start_freq;
		$hash_site{$hash_comp{'primary_id'}}{$hash_comp{'subject_strand'}}{$hash_comp{'key_start'}}{$hash_comp{'key_stop'}}=$key_stop_freq;
	}else{
		$hash_site{$hash_comp{'primary_id'}}{$hash_comp{'subject_strand'}}{$hash_comp{'key_start'}}{$hash_comp{'key_stop'}}=$key_stop_freq;
	}
	unless (defined $hash_site{$hash_comp{'primary_id'}}{$hash_comp{'subject_strand'}}{$hash_comp{'key_stop'}}{$hash_comp{'key_stop'}}){
		$hash_site{$hash_comp{'primary_id'}}{$hash_comp{'subject_strand'}}{$hash_comp{'key_stop'}}{$hash_comp{'key_stop'}}=$key_stop_freq;
		$hash_site{$hash_comp{'primary_id'}}{$hash_comp{'subject_strand'}}{$hash_comp{'key_stop'}}{$hash_comp{'key_start'}}=$key_start_freq;
	}else{
		$hash_site{$hash_comp{'primary_id'}}{$hash_comp{'subject_strand'}}{$hash_comp{'key_stop'}}{$hash_comp{'key_start'}}=$key_start_freq;
	}

	my $site = $hash_comp{'primary_id'}." ".$hash_comp{'subject_strand'}." ".$hash_comp{'subject_beg'}." ".$hash_comp{'subject_end'};
    my $key_site = $hash_comp{'primary_id'}." ".$hash_comp{'subject_strand'}." ".$hash_comp{'key_start'}." ".$hash_comp{'key_stop'};
	my $outtable = "$mature_id\t$star_id\t$mature_site\t$star_site\t$freq_mature\t$freq_star\t$mature_seq\t$star_seq\t$site\t$hash_comp{'pre_seq'}";
	my $outresult;
	my @reads;
    my $lr;
    my $bef;

    my $after;
    my $rseq;

    my @sread;
    $hash_comp{"pri_seq"} =~ tr/Tt/Uu/;
    my @pseq = split(//,lc $hash_comp{"pri_seq"});

    my $spacer = " " x ($col1_width - length("pri_seq"));
    my $shift;
    my %reads_hash;

	#print ">$hash_comp{\"pri_id\"}\n";
	$outresult .= ">$hash_comp{\"pri_id\"}\n";
    my $spacer2 = 14;
    my $spacer2s=0;

    if(not $hash_comp{"problem_structure"} and not $hash_comp{"problem_signature"}){
        $outresult .= "microRNA name       $key_site\n";
        $outresult .= "precursor coord     $site\n";
        $outresult .= "mature sequence     $hash_comp{\"mature_seq\"}\n";
        $outresult .= "star sequence       $hash_comp{\"star_seq\"}\n";
        $spacer2s = " " x ($spacer2 - length($hash_comp{"score"}));
        #print "score total\t\t$spacer2s$hash_comp{\"score\"}\n";
		$outresult .= "score total\t\t$spacer2s$hash_comp{\"score\"}\n";
        $spacer2s = " " x ($spacer2 - length($hash_comp{"score_star"}));
        #print "score for star read(s)\t$spacer2s$hash_comp{\"score_star\"}\n";
		$outresult .= "score for star read(s)\t$spacer2s$hash_comp{\"score_star\"}\n";
        $spacer2s = " " x ($spacer2 - length($hash_comp{"score_freq"}));
        #print "score for read counts\t$spacer2s$hash_comp{\"score_freq\"}\n";
		$outresult .= "score for read counts\t$spacer2s$hash_comp{\"score_freq\"}\n";
        $spacer2s = " " x ($spacer2 - length($hash_comp{"score_mfe"}));
        #print "score for mfe\t\t$spacer2s$hash_comp{\"score_mfe\"}\n";
		$outresult .= "score for mfe\t\t$spacer2s$hash_comp{\"score_mfe\"}\n";
	if($file_mature){
        $spacer2s = " " x ($spacer2 - length($hash_comp{"score_seed"}));
        #print "score for cons. seed\t$spacer2s$hash_comp{\"score_seed\"}\n";
		$outresult .= "score for cons. seed\t$spacer2s$hash_comp{\"score_seed\"}\n";
    }
	if($file_mature and defined $hash_comp{"seed_family"}){
        $spacer2s = " " x ($spacer2 - length($hash_comp{"seed_family"}));
        #print "miRNA with same seed\t$spacer2s$hash_comp{\"seed_family\"}\n";
		$outresult .= "miRNA with same seed\t$spacer2s$hash_comp{\"seed_family\"}\n";

    }
        $spacer2s = " " x ($spacer2 - length($hash_comp{"freq_total"}));
        #print "total read count\t$spacer2s$hash_comp{\"freq_total\"}\n";
		$outresult .= "total read count\t$spacer2s$hash_comp{\"freq_total\"}\n";
        $spacer2s = " " x ($spacer2 - length($hash_comp{"freq_mature"}));
        #print "mature read count\t$spacer2s$hash_comp{\"freq_mature\"}\n";
        $outresult .= "mature read count\t$spacer2s$hash_comp{\"freq_mature\"}\n";
		$spacer2s = " " x ($spacer2 - length($hash_comp{"freq_loop"}));
        #print "loop read count\t\t$spacer2s$hash_comp{\"freq_loop\"}\n";
        $outresult .= "loop read count\t\t$spacer2s$hash_comp{\"freq_loop\"}\n";
		$spacer2s = " " x ($spacer2 - length($hash_comp{"freq_star"}));
        #print "star read count\t\t$spacer2s$hash_comp{\"freq_star\"}\n";
		$outresult .= "star read count\t\t$spacer2s$hash_comp{\"freq_star\"}\n";


    }else{
	#structure problems
	if(defined $hash_comp{"problem_structure"}){ print $hash_comp{"problem_structure"};}
	if(defined $hash_comp{"problem_mature_lng"}){ print $hash_comp{"problem_mature_lng"};}
	if(defined $hash_comp{"problem_duplex_bp"}){ print $hash_comp{"problem_duplex_bp"};}
	if(defined $hash_comp{"problem_duplex_lng"}){ print $hash_comp{"problem_duplex_lng"};}
	if(defined $hash_comp{"problem_struct_mature"}){ print $hash_comp{"problem_struct_mature"};}
	if(defined $hash_comp{"problem_struct_star"}){ print $hash_comp{"problem_struct_star"};}
	if(defined $hash_comp{"problem_struct_loop"}){ print $hash_comp{"problem_struct_loop"};}
	if(defined $hash_comp{"problem_struct_us_flanks"}){ print $hash_comp{"problem_struct_us_flanks"};}
	if(defined $hash_comp{"problem_struct_ds_flanks"}){ print $hash_comp{"problem_struct_ds_flanks"};}
	if(defined $hash_comp{"problem_struct_bifurcation"}){ print $hash_comp{"problem_struct_bifurcation"};}

	#signature problems
	if(defined $hash_comp{"problem_signature"}){ print $hash_comp{"problem_signature"};}
	if(defined $hash_comp{"problem_read_cnt"}){ print $hash_comp{"problem_read_cnt"};}
	if(defined $hash_comp{"problem_Dicer_signature"}){ print $hash_comp{"problem_Dicer_signature"};}
    }


    #### printing alignment;
    #print "exp";
	$outresult .= "exp";

    #print " " x ($col1_width-3);
	$outresult .= " " x ($col1_width-3);
    if($hash_comp{"problem_structure"}){
		#print "n" x length($hash_comp{"pri_seq"});
		$outresult .= "n" x length($hash_comp{"pri_seq"});
    }else{
		#print "f" x $hash_comp{"flank_first_end"};
		$outresult .= "f" x $hash_comp{"flank_first_end"};
		if($hash_comp{"mature_beg"}  < $hash_comp{"loop_beg"}){
			$bef = "M" x ($hash_comp{"mature_end"}-$hash_comp{"mature_beg"}+1);
			#print $bef;
			$outresult .= $bef;
			$bef = "l" x ($hash_comp{"loop_end"}-$hash_comp{"loop_beg"}+1);
			#print $bef;
			$outresult .= $bef;
			$bef = "S" x ($hash_comp{"star_end_exp"}-$hash_comp{"star_beg_exp"}+1);
			#print $bef;
			$outresult .= $bef;
		}else{
			$bef = "S" x ($hash_comp{"star_end_exp"}-$hash_comp{"star_beg_exp"}+1);
			#print $bef;
			$outresult .= $bef;
			$bef = "l" x ($hash_comp{"loop_end"}-$hash_comp{"loop_beg"}+1);
			#print $bef;
			$outresult .= $bef;
			$bef = "M" x ($hash_comp{"mature_end"}-$hash_comp{"mature_beg"}+1);
			#print $bef;
			$outresult .= $bef;
		}
		$bef = "f" x ($hash_comp{"pri_end"}-$hash_comp{"flank_second_beg"}+1);
		#print $bef;
		$outresult .= $bef;
    }

    if(defined $hash_comp{"star_beg_obs"}){

	#print "\nobs";
	$outresult .= "\nobs";
	#print " " x ($col1_width-3);
	$outresult .= " " x ($col1_width-3);
	if($hash_comp{"problem_structure"}){
	    #print "n" x length($hash_comp{"pri_seq"});
		$outresult .= "n" x length($hash_comp{"pri_seq"});
	}else{
	    if($hash_comp{"mature_beg"}  < $hash_comp{"loop_beg"}){

		$bef="f" x ($hash_comp{"mature_beg"}-1);
		#print $bef;
		$outresult .= $bef;
		$bef = "M" x ($hash_comp{"mature_end"}-$hash_comp{"mature_beg"}+1);
		#print $bef;
		$outresult .= $bef;
		$bef = "l" x ($hash_comp{"star_beg_obs"}-$hash_comp{"mature_end"}-1);
		#print $bef;
		$outresult .= $bef;
		$bef = "S" x ($hash_comp{"star_end_obs"}-$hash_comp{"star_beg_obs"}+1);
		#print $bef;
		$outresult .= $bef;
		$bef="f" x ($hash_comp{"pri_end"}-$hash_comp{"star_end_obs"});
		#print $bef;
		$outresult .= $bef;

	    }else{

		$bef="f" x ($hash_comp{"star_beg_obs"}-1);
		#print $bef;
		$outresult .= $bef;
		$bef = "S" x ($hash_comp{"star_end_obs"}-$hash_comp{"star_beg_obs"}+1);
		#print $bef;
		$outresult .= $bef;
		$bef = "l" x ($hash_comp{"mature_beg"}-$hash_comp{"star_end_obs"}-1);
		#print $bef;
		$outresult .= $bef;
		$bef = "M" x ($hash_comp{"mature_end"}-$hash_comp{"mature_beg"}+1);
		#print $bef;
		$outresult .= $bef;
		$bef="f" x ($hash_comp{"pri_end"}-$hash_comp{"mature_end"});
		#print $bef;
		$outresult .= $bef;
	    }
	}
    }

    #print "\npri_seq$spacer",lc $hash_comp{"pri_seq"},"\n";
	my $lc_pre_seq = lc $hash_comp{"pri_seq"};
	$outresult .= "\npri_seq$spacer$lc_pre_seq\n";
    $spacer = " " x ($col1_width - length("pri_struct"));
    #print "pri_struct$spacer$hash_comp{\"pri_struct\"}\t#MM\n";
	$outresult .= "pri_struct$spacer$hash_comp{\"pri_struct\"}\t#MM\n";
	my @lines = split /\n/, $lines;
	@reads = map $_->[0], sort {$a->[1] <=> $b->[1] or $a->[2] <=> $b->[2]} map [(split/\s+/)], @lines;
    foreach (@reads){
        my $query_id = $_;
        my $query_seq = $hash_query{$query_id}{"seq"};
        $query_seq =~ tr/t/u/;
        my $subject_beg = $hash_query{$query_id}{"db_beg"};
        my $subject_end = $hash_query{$query_id}{"db_end"};
        my $strand = $hash_query{$query_id}{"strand"};
        my $subject_lng = length($hash_comp{'pre_seq'});
        if ( ($strand eq '+') && ($subject_beg >=1) && ($subject_end <= $subject_lng) ){
            $spacer = ' ' x ($col1_width - length($query_id));
            $bef = '.' x ($subject_beg - 1);
            $after = '.' x ($subject_lng - $subject_end);
            $outresult .= "$query_id$spacer$bef$query_seq$after\t0\n";
        }
    }
	if ($value < $hash_value{$hash_comp{'primary_id'}}{$hash_comp{'subject_strand'}}{$hash_comp{'key_start'}}{$hash_comp{'key_stop'}}){
		$hash_value{$hash_comp{'primary_id'}}{$hash_comp{'subject_strand'}}{$hash_comp{'key_start'}}{$hash_comp{'key_stop'}} = $value;
		#print STDERR "flank1=$flank1; flank2=$flank2; value=$value\n";
		$hash_stored{$hash_comp{'primary_id'}}{$hash_comp{'subject_strand'}}{$hash_comp{'key_start'}}{$hash_comp{'key_stop'}} = $outresult;
		$hash_list{$hash_comp{'primary_id'}}{$hash_comp{'subject_strand'}}{$hash_comp{'key_start'}}{$hash_comp{'key_stop'}} = $outtable;
	    my $mature_region = "$hash_comp{'mature_beg'}-$hash_comp{'mature_end'}";
        my $star_region = "$hash_comp{'star_beg_exp'}-$hash_comp{'star_end_exp'}";
        my $precursor_seq = $lc_pre_seq;
        my $precursor_struct = $hash_comp{'pri_struct'};
        my $struct_id = $hash_comp{'primary_id'};
        $hash_figure{$hash_comp{'primary_id'}}{$hash_comp{'subject_strand'}}{$hash_comp{'key_start'}}{$hash_comp{'key_stop'}} = "$precursor_seq\t$precursor_struct\t$mature_region\t$star_region";
    }
	return;
}

sub error_checking{
	my @subjects = sort keys %hash_site;
	for my $subject (@subjects){
		my @strands = sort keys %{$hash_site{$subject}};
		for my $strand (@strands){
			my @roots = sort {$a<=>$b} keys %{$hash_site{$subject}{$strand}};
			for my $root (@roots){
				my @stems = sort {$a<=>$b} keys %{$hash_site{$subject}{$strand}{$root}};
				next if @stems < 3; #only concerns @stems greater than 3
				my $root_freq = $hash_site{$subject}{$strand}{$root}{$root};
				my (@tmp_stems,$exp_stem,$stem,$i);
				for $stem (@stems){
					next if $stem == $root;
					@tmp_stems=sort keys %{$hash_site{$subject}{$strand}{$stem}};
					unless (@tmp_stems < 3){
						for $i (0..$#stems){
							if ($stems[$i] == $root){
								if ($i==0){
									$exp_stem = $stems[$i+1];
								}elsif ($i==$#stems){
									$exp_stem = $stems[$i-1];
								}else{
									if (($root - $stems[$i-1]) < ($stems[$i+1] - $root)){
										$exp_stem = $stems[$i-1];
									}else{
										$exp_stem = $stems[$i+1];
									}
								}
							}
							last if defined $exp_stem;
						}
					}
					last if defined $exp_stem;
				}
				if (defined $exp_stem){
					for $stem (@stems){
						unless ($stem == $exp_stem){
							$hash_dead{$subject}{$strand}{$root}=$stem;
							$hash_dead{$subject}{$strand}{$stem}=$root;
						}
					}
				}else{
					$exp_stem=(sort {$hash_site{$subject}{$strand}{$root}{$b}<=>$hash_site{$subject}{$strand}{$root}{$a}} keys %{$hash_site{$subject}{$strand}{$root}})[0];
					for $stem (@stems){
						unless ($stem == $exp_stem){
							$hash_dead{$subject}{$strand}{$root}=$stem;
							$hash_dead{$subject}{$strand}{$stem}=$root;
						}
					}
				}
			}
		}
	}
}

sub print_results{
	open RESULT,">$result_out" or die "Can not create file $result_out: $!\n";
	open LIST,">$table_out" or die "Can not create file $table_out: $!\n";
    print LIST "Mature_id\tStar_id\tMature_arm\tStar_arm\tMature_freq\tStar_freq\tMature_seq\tStar_seq\tPrecursor_site\tPrecursor_seq\tKey_coord\n";	
	my @subjects = sort keys %hash_stored;
	for my $subject (@subjects){
		my @strands = sort keys %{$hash_stored{$subject}};
		for my $strand (@strands){
			my @begs = sort {$a<=>$b} keys %{$hash_stored{$subject}{$strand}};
			for my $beg (@begs){
				my @ends = sort {$a<=>$b} keys %{$hash_stored{$subject}{$strand}{$beg}};
				for my $end (@ends){
					next if (defined $hash_dead{$subject}{$strand}{$beg} && $hash_dead{$subject}{$strand}{$beg} == $end);
					print RESULT $hash_stored{$subject}{$strand}{$beg}{$end};
					
                    my $info = $hash_figure{$subject}{$strand}{$beg}{$end};

                    my ($precursor_seq,$precursor_struct,$mature_region,$star_region)=split /\s+/,$info;
                    
                    my $figure = structure_prepare($precursor_seq,$precursor_struct,$mature_region,$star_region);
                    print RESULT "\n", $figure, "\n\n\n\n";

					print LIST $hash_list{$subject}{$strand}{$beg}{$end},"\t","$subject $strand $beg $end","\n";
					
				}
			}

		}
	
	}
	close RESULT;
	close LIST;
	return;
}

sub print_hash_bp{

    #prints the 'bp' hash

    my @keys=sort {$a<=>$b} keys %hash_bp;
    foreach my $key(@keys){
	my $value=$hash_bp{$key};
	print "$key\t$value\n";
    }
    print "\n";
}

sub contained{

    #Is the stretch defined by the first positions contained in the stretch defined by the second?

    my($beg1,$end1,$beg2,$end2)=@_;

    testbeginend($beg1,$end1,$beg2,$end2);

    if($beg2<=$beg1 and $end1<=$end2){
	return 1;
    }else{
	return 0;
    }
}

sub testbeginend{

    #Are the beginposition numerically smaller than the endposition for each pair?

    my($begin1,$end1,$begin2,$end2)=@_;

    unless($begin1<=$end1 and $begin2<=$end2){
	print STDERR "beg can not be larger than end for $db_old\n";
	exit;
    }
}

sub rev_pos{

#   This subroutine reverses the begin and end positions

    my($beg,$end,$lng)=@_;

    my $new_end=$lng-$beg+1;
    my $new_beg=$lng-$end+1;

    return($new_beg,$new_end);
}

sub round {

    #rounds to nearest integer

    my($number) = shift;
    return int($number + .5);

}

sub rev{

    #reverses the order of nucleotides in a sequence

    my($sequence)=@_;

    my $rev=reverse $sequence;

    return $rev;
}

sub com{

    #the complementary of a sequence

    my($sequence)=@_;

    $sequence=~tr/acgtuACGTU/TGCAATGCAA/;

    return $sequence;
}

sub revcom{

    #reverse complement

    my($sequence)=@_;

    my $revcom=rev(com($sequence));

    return $revcom;
}

sub max2 {

    #max of two numbers

    my($a, $b) = @_;
    return ($a>$b ? $a : $b);
}

sub min2  {

    #min of two numbers

    my($a, $b) = @_;
    return ($a<$b ? $a : $b);
}

sub score_freq{

#   scores the count of reads that map to the potential precursor
#   Assumes geometric distribution as described in methods section of manuscript

    my $freq=shift;

    #parameters of known precursors and background hairpins
    my $parameter_test=0.999;
    my $parameter_control=0.6;

    #log_odds calculated directly to avoid underflow
    my $intercept=log((1-$parameter_test)/(1-$parameter_control));
    my $slope=log($parameter_test/$parameter_control);
    my $log_odds=$slope*$freq+$intercept;

    #if no strong evidence for 3' overhangs, limit the score contribution to 0
    unless($hash_comp{"star_read"}){$log_odds=min2($log_odds,0);}

    return $log_odds;
}

sub score_mfe{

#   scores the minimum free energy in kCal/mol of the potential precursor
#   Assumes Gumbel distribution as described in methods section of manuscript

    my ($mfe,$pri_lng)=@_;

    #normalize the MFE by length
    my $mfe_adj=$mfe/$pri_lng;

    #parameters of known precursors and background hairpins, scale and location
    #my $prob_test=prob_gumbel_discretized($mfe_adj,5.5,32);
    #my $prob_background=prob_gumbel_discretized($mfe_adj,4.8,23);
    #instead of finding individual functions for real and bgr, the one for log-odds was directly obtained.
    #its a sigmoid func. with +ve x.f(x)=a/(b+exp(x*c))
    my $log_odds=0;
    $log_odds=$param_a/($param_b+exp($mfe_adj*$param_c));

    return $log_odds;
}

sub prob_gumbel_discretized{

#   discretized Gumbel distribution, probabilities within windows of 1 kCal/mol
#   uses the subroutine that calculates the cdf to find the probabilities

    my ($var,$scale,$location)=@_;

    my $bound_lower=$var-0.5;
    my $bound_upper=$var+0.5;

    my $cdf_lower=cdf_gumbel($bound_lower,$scale,$location);
    my $cdf_upper=cdf_gumbel($bound_upper,$scale,$location);

    my $prob=$cdf_upper-$cdf_lower;

    return $prob;
}

sub cdf_gumbel{

#   calculates the cumulative distribution function of the Gumbel distribution

    my ($var,$scale,$location)=@_;

    my $cdf=$e**(-($e**(-($var-$location)/$scale)));

    return $cdf;
}

# draw hairpin-like figure with dot-bracket notation
# multi-loops are deleted

sub structure_prepare{
	my ($seq,$struct,$mature_site,$star_site) = @_;
	my ($mature_beg,$mature_end) = split /-/,$mature_site;
	my ($star_beg,$star_end) = split /-/,$star_site;
    $mature_beg -= 1;
    $mature_end -= 1;
    $star_beg -= 1;
    $star_end -= 1;
	die "Mature coordinate error!\n" unless ($mature_beg >= 0 && $mature_end >= $mature_beg);
	die "Star coordinate error!\n" unless ($star_beg >= 0 && $star_end >= $star_beg);
	die "Mature overlaps in Star!\n" 
		if ($mature_beg <= $star_beg && $star_beg <= $mature_end) or
			($star_beg <= $mature_beg && $mature_beg <= $star_end);
	my $mature_lng = $mature_end - $mature_beg + 1;
	my $star_lng = $star_end - $star_beg + 1;
	$seq =~ s/^(\w{$mature_beg})(\w{$mature_lng})/$1\U$2/;
	$seq =~ s/^(\w{$star_beg})(\w{$star_lng})/$1\U$2/;
    my $formated_struct = format_multi_loops($struct);
	my ($figure) = draw_hairpin_figure($seq, $formated_struct);
	return $figure;
}

sub format_multi_loops{
	my ($struct) = @_;
	my $lng = length $struct;
	my @bps;
	my %hash_bp;
	my $sum=0;
	my $max=0;
	my $max_site;
	my $formated_struct;
	my @structs = split //, $struct;
	for(my $pos=0;$pos<=$lng-1;$pos++){
		my $struct_pos = $structs[$pos];
		if ($struct_pos eq '('){
			push(@bps,$pos);
			$sum ++;
		}
		if ($struct_pos eq ')'){
			my $pos_prev=pop(@bps);
			$sum --;
			$hash_bp{$pos}=$pos_prev;
		}
		if ($sum > $max){
			$max = $sum;
			$max_site = $pos;
		}
	}
	for my $left (sort {$a<=>$b} keys %hash_bp){
		my $right = $hash_bp{$left};
		if ( ($left < $max_site && $right < $max_site) || ($left > $max_site && $right > $max_site) ){
			$structs[$left] = '.';
			$structs[$right] = '.';
		}
	}
	for my $base (@structs){
		$formated_struct .= $base;
	}
	return $formated_struct;
}


sub draw_hairpin_figure{
    my ($seq, $struct) = @_;

    return 'Missing sequence or structure'
        unless length $seq > 0 and length $struct > 0;
    return 'Missmatch length of sequence and structure'
        unless length $seq == length $struct;
    return 'Illegal character in dot-bracket notation'
        if $struct =~ /[^\.\(\)]/;
    
    my (@l1, @l2, @l3, @l4, @l5);
    
    my $len   = length $struct; 
    my $table = make_pair_table($struct);
    my @left  = sort{ $a <=> $b } keys %$table;
    my @right = map { $$table{$_} } @left;

    #ssRNA
    if ($len - $right[0] >= $left[0] - 1){
        for (1..($len - $right[0] - ($left[0] - 1))){
            push @l1, '-';
            push @l2, ' ';
            push @l3, ' ';
            push @l4, ' ';
            push @l5, substr($seq, $len - $_, 1);
        }
        for (1..($left[0] - 1)){
            push @l1, substr($seq, $_ - 1, 1);
            push @l2, ' ';
            push @l3, ' ';
            push @l4, ' ';
            push @l5, substr($seq, 
                $len - ($len - $right[0] - ($left[0] - 1)) - $_, 1);
        }
    }
    else{
        for (1..($left[0] - 1 - ($len - $right[0]))){
            push @l1, substr($seq, $_ - 1, 1);
            push @l2, ' ';
            push @l3, ' ';
            push @l4, ' ';
            push @l5, '-';
        }
        for (1..($len - $right[0])){
            push @l1, substr($seq, $_, 1);
            push @l2, ' ';
            push @l3, ' ';
            push @l4, ' ';
            push @l5, substr($seq, $len - $_, 1);
        }
    }
    
    # stem region
    my $next5 = $left[0];
    my $next3 = $right[0];
    my ($n5, $n3, $asy);
    while ($next5 <= $left[-1]) {
        # stem
        if ($next5 ~~ @left and $next3 ~~ @right) {
            while ($next5 ~~ @left and $next3 ~~ @right){
                # print " $next5 - $next3\n";
                push @l1, ' ';
                push @l2, substr($seq, $next5 - 1, 1);
                push @l3, '|';
                push @l4, substr($seq, $$table{$next5} - 1, 1);
                push @l5, ' ';
                $next5 ++;
                $next3 --;
            }
        }
        # 5' gap
        elsif ($next5 !~ @left and $next3 ~~ @right) {
            # print "[5' gap],$next5,$next3\n";
            $n5 = 0;
            $n5 ++ until ($next5 + $n5) ~~ @left;
            for(1..$n5){
                push @l1, substr($seq, $next5 + $_ - 2, 1);
                push @l2, ' ';
                push @l3, ' ';
                push @l4, ' ';
                push @l5, '-';
            }
            $next5 += $n5;
        }
        # 3' gap
        elsif ($next5 ~~ @left and $next3 !~ @right) {
            # print "[3' gap], $next5,$next3\n";
            $n3 = 0;
            $n3 ++ until ($next3 - $n3) ~~ @right;
            for(1..$n3){
                push @l1, '-';
                push @l2, ' ';
                push @l3, ' ';
                push @l4, ' ';
                push @l5, substr($seq, $next3 - $_, 1);
            }
            $next3  -= $n3;
        }
        # bulge
        # elsif ($next5 !~ @left and $next3 !~ @right){
        else{
            $n5 = 0;
            $n5 ++ until ($next5 + $n5) ~~ @left;
            $n3 = 0;
            $n3 ++ until ($next3 - $n3) ~~ @right;
            # print "[bulge] $next5,$next3\n";
            if ($n5 > $n3) {
                for(1..$n3){
                    push @l1, substr($seq, $next5 + $n5 - $n3 + $_ - 2, 1);
                    push @l2, ' ';
                    push @l3, ' ';
                    push @l4, ' ';
                    push @l5, substr($seq, $next3 - $_, 1);
                }
                for (1..($n5-$n3)) {
                    push @l1, substr($seq, $next5 + $_ - 2, 1);
                    push @l2, ' ';
                    push @l3, ' ';
                    push @l4, ' ';
                    push @l5, '-';
                }
            }
            elsif ($n5 < $n3) {
                for(1..$n5){
                    push @l1, substr($seq, $next5 + $_ - 2, 1);
                    push @l2, ' ';
                    push @l3, ' ';
                    push @l4, ' ';
                    push @l5, substr($seq, $next3 - ($n3 -$n5) - $_, 1);
                }
                for (1..($n3-$n5)) {
                    push @l1, '-';
                    push @l2, ' ';
                    push @l3, ' ';
                    push @l4, ' ';
                    push @l5, substr($seq, $next3 - $_, 1);
                }
            }
            else {
                for (1..$n5) {
                    push @l1, substr($seq, $next5 + $_ - 2, 1);
                    push @l2, ' ';
                    push @l3, ' ';
                    push @l4, ' ';
                    push @l5, substr($seq, $next3 - $_, 1);
                }
            }
            $next5 += $n5;
            $next3 -= $n3;
        }
    }
    
    # terminal loop
    my $loop = $right[-1] - $left[-1] - 1;
    my $n = int (($loop - 2) / 2);
    # print "loop: $loop,$n\n";
    if ($n > 0) {
        for (1..$n){
            push @l1, substr($seq, $next5 + $_ - 2, 1);
            push @l2, ' ';
            push @l3, ' ';
            push @l4, ' ';
            push @l5, substr($seq, $next3 - $_, 1);
        }
        $next5 += $n;
        $next3 -= $n;
        
        push @l1, ' ';
        push @l2, substr($seq, $next5 - 1, 1);
        push @l3, $loop - 2 * ($n + 1) > 0 
                    ? substr($seq, $next5, 1)
                    : ' ';
        push @l4, substr($seq, $next3 + 1, 1);
        push @l5, ' ';
    }
    elsif ($loop == 3 or $loop == 2) {
        push @l1, ' ';
        push @l2, substr($seq, $next5 - 1, 1);
        push @l3, ' ';
        push @l4, substr($seq, $next3 + 1, 1);
        push @l5, ' ';
        
        if ($loop == 3) {
            push @l1, ' ';
            push @l2, ' ';
            push @l3, substr($seq, $next5, 1);
            push @l4, ' ';
            push @l5, ' ';
        }
    }
    
    
    # out put
    my $s1 = join '', @l1;
    my $s2 = join '', @l2;
    my $s3 = join '', @l3;
    my $s4 = join '', @l4;
    my $s5 = join '', @l5;
    
    my $figure = '';
    $figure .= $s1."\n" . $s2. "\n" . $s3. "\n" . $s4. "\n" . $s5;
    
    # for other using
    my @struct_data  = ($s1, $s2, $s3, $s4, $s5);
    
    # reversed
    my @struct_data_reversed = map { my $tmp = reverse $_ } @struct_data;
    my $graph = join "\n", @struct_data_reversed;
    
    return ($figure, $graph, \@struct_data, \@struct_data_reversed);
}

# returns hash representation of dot-bracket notation.
# table{i} is j if (i.j) pair.
sub make_pair_table{
    my ( $struct ) = @_;
    my ( $i, $j, $length, $table, $stack );
    $length = length $struct;
    my @struct_data = split "", $struct;
    for ( $i = 1; $i <= $length; $i ++ ) {
        if ( $struct_data[$i - 1] eq '(' ){
            unshift @$stack, $i;
        }elsif ($struct_data[$i - 1] eq ')' ){
            if ( @$stack == 0 ) {
                die "unbalanced brackets $struct\n";
                return undef;
            }
            $j = shift @$stack;
            $$table{$j} = $i;
        }
    }
    if ( @$stack != 0 ) {
        die "unbalanced brackets $struct\n";
        return undef;
    }
    undef @$stack;
    return $table;
}

sub var_check {
    if ($options{'v'}){
        die "\nmiR_island_core.pl version = $version\n\n";
        exit 0;
    }
    if ($options{'h'}){ &var_error(); }
    if ($options{'b'}){
        $sorted_bam_file = $options{'b'};
    }else{
        &var_error();
    }
    if ($options{'c'}){
        $cotyledon = $options{'c'};
    }else{
        &var_error();
    }
    if ($options{'g'}){
        $genome = $options{'g'};
    }else{
        &var_error();
    }
    if ($options{'l'}){
        $col1_width = $options{'l'};
    }else{
        $col1_width = 45;
    }
    if ($options{'m'}){
        $file_mature = $options{'m'};
    }
    if ($options{'r'}){
        $result_out = $options{'r'};
    }else{
        &var_error();
    }
    if ($options{'s'}){
        $file_struct = $options{'s'};
    }else{
        &var_error();
    }
    if ($options{'t'}){
        $table_out = $options{'t'};
    }else{
        &var_error();
    }
}

sub var_error {
    print "\n";
    print " miR_island_core.pl is the miR-island package's core algorithm for plant microRNA annotation\n\n";
    print " Version = $version\n\n";
    print " WARNING: You did not provide enough information!\n\n" unless exists $options{'h'};
    print " Usage: miR_island_core.pl -b <reads_vs_genome_sorted.bam> -c <monocot|dicot> -g <genome.fa> \\\n";
    print "           -s <structure_file> -t <table_out> -r <result_out> [-m <miRNA.fa>] [-l <col_width>]\n";
    print "\n";
    print " REQUIRED:\n";
    print " -b <sorted.bam>        profile of reads mapped to the related genome in a formatted, sorted\n";
    print "                        and indexed BAM file\n";
    print " -c <monocot|dicot>     cotyledon type: monocot | dicot\n";
	print " -g <genome.fa>         reference genome file in multi-fasta format\n";
    print " -s <structure_file>    secondary structure file produced by RNAfold\n";
    print " -r <result_out>        an intact output file that annotates miRNA and MIRNA genes\n";
    print " -t <table_out>         a brief tab-delimited text file that annotates miRNA and MIRNA genes\n";
    print "\n";
    print " OPTIONAL:\n";
    print " -l <column_width>      coloum 1 width, if identifier is long, add this value (def: 40)\n";
    print " -m <miRNA.fa>          known miRNAs of the related species in multi-fasta format\n";
    print " -v                     show the current version number\n";
	print " -h                     print intact help message\n";
    print "\n";

    exit 1;
}
