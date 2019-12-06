#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;

use constant DEBUG => 0;
my $version = "1.0";

my ($cotyledon,$mfe,$min_freq,$genome,$flag,$index,$genome_index,$known_miRNA,$hairpin_lng,$map_time,$nb_process,$reads,$time,$coord_file);

my $hints = "Please check the file for the following issues:\n
i).   Identifiers should be in the format of '>identifier';\n
ii).  Identifiers are not allowed to have white spaces;\n
iii). Identifiers should be unique in the whole file;\n
iv).  Input sRNA reads must end with _xNumber, while Number stands for read count;\n
v).   Sequences are not allowed to comprise characters other than [YATCGNatcgn].\n
Note: If the error happens in the input sRNA reads, please run \*format_input_reads.pl\*\n
      to format the file.\n";
my $log = join(" ", $0, @ARGV);

my %opts=();
getopts("c:e:f:g:hi:k:l:m:p:r:x:v",\%opts);
&var_check();

my $current_time = time;
my $folder = "MiR_Island_run";
my $dir = "$folder/run_${current_time}";
my $dir_tmp = "$dir/tmp";
my $dir_struct = "$dir/struct_tmp";
my $log_file="$dir/Run_MiRNA_Island.log";
my $reads_collapsed = "$dir/reads_collapsed.fa";
my $bam = "$dir/reads_mapped_genome.bam";
my $sorted_bam = "$dir/reads_mapped_genome_sorted.bam";
my $structure = "$dir/precursors.struct";
my $output = "$dir/result.out";
my $table = "$dir/result.table";
my $table_flagged = "$dir/result.table.flagged";
my $denovo_diff = "$dir/denovo_miRNA_expression.list";
my $known_diff = "$dir/known_miRNA_expression.list";

$time = localtime;
$log .= "\n\n$time\nmkdir $dir\n\n";

make_dir_tmp() unless DEBUG;

open LOG, ">", $log_file or die "Cannot create file $log_file: $!\n";
record_log($log);

$time = localtime;
$log = "$time\n";
record_log($log);

$log = "Checking installed binaries ... ";
record_log($log);
check_installed_binaries() unless DEBUG;
$log = "done\n\n";
record_log($log);

$time = localtime;
$log = "$time\n";
record_log($log);

$log = "Checking CPU cores ... ";
record_log($log);
check_cpu_cores() unless DEBUG;
$log = "done\n\n";
record_log($log);

$time = localtime;
$log = "$time\n";
record_log($log);

$log = "Checking input files\n";
record_log($log);
check_input_files() unless DEBUG;

$log = "\n";
record_log($log);
$time = localtime;
$log = "$time\n";
record_log($log);

$log = "Analyzing known miRNA differential gene expression ... ";
record_log($log);
diff_mature_analysis() unless DEBUG;

$time = localtime;
$log = "$time\n";
record_log($log);

if (defined $flag){
    $genome = "$dir/pseudo_genome.fa";
    $index = "$dir/pseudo.index";
    $genome_index = "$dir/pseudo_genome.fa";
    $output = "$dir_tmp/result.out";
    $table = "$dir_tmp/result.table";
}

$log = "Checking genome bowtie index file ... ";
record_log($log);
check_index() unless DEBUG;
$log = "\n";
record_log($log);

$time = localtime;
$log = "$time\n";
record_log($log);

$log = "Mapping reads to reference genome ... ";
record_log($log);
reads_map_genome() unless DEBUG;
$log = "done\n\n";
record_log($log);

$time = localtime;
$log = "$time\n";
record_log($log);

$log = "Excising potential precursors ... ";
record_log($log);
excise_potential_precursors() unless DEBUG;
$log = "done\n\n";
record_log($log);

$time = localtime;
$log = "$time\n";
record_log($log);

$log = "Preparing second structure, this step may take several minites ... ";
record_log($log);
prepare_second_structure() unless DEBUG;
$log = "done\n\n";
record_log($log);

$time = localtime;
$log = "$time\n";
record_log($log);

$log = "Final miRNA prediction ... ";
record_log($log);
mirna_prediction() unless DEBUG;
$log = "done\n\n";
record_log($log);

$time = localtime;
$log = "$time\n";
record_log($log);

$log = "Fianl check ... ";
record_log($log);
if (defined $flag){
    final_check() unless DEBUG;
}
$log = "done\n\n";
record_log($log);

$time = localtime;
$log = "Flag genes ... ";
record_log($log);
my $in_table = "$dir/result.table";
my $out_table = "$dir/result.table.flagged";
flag_gene_out($in_table,$out_table) unless DEBUG;
my $in_result = "$dir/result.out";
my $out_result = "$dir/result.out.flagged";
flag_gene_out($in_result,$out_result) unless DEBUG;
if (-r $out_table){
    system "mv $out_table $in_table";
}
if (-r $out_result){
    system "mv $out_result $in_result";
}
$log = "done\n\n";

$time = localtime;
$log = "$time\n";
record_log($log);

$log = "Denovo differential miRNA expression analysis ... ";
record_log($log);
denovo_diff_miRNA_analysis() unless DEBUG;
$log = "done\n\n";
record_log($log);

#remove temporary files and folders
system "rm -rf $dir_tmp";
system "rm -rf $dir_struct";
system "rm -f $bam";

$log = "Mission complete.\n\n";
record_log($log);

$time = localtime;
$log = "$time\n\n";
record_log($log);
close LOG;
exit 0;

########################################## SUBROUTINE ################################################

sub check_installed_binaries{
	my $value;
	
    $value = check_binary("perl -e \'print \"pass\" if $]>=5.010;\'","pass");
    if ($value){
        die "Error: \tThe current perl version $] is too low\nPlease update perl version greater than 5.010\n";
    }
	$value = check_binary("perl -V 2>&1","useithreads=define");
	if ($value){
		die "Error: \tThe current perl does not support multiple threads\nPlease restore perl with multiple threads support.\n";
	}
	
	$value = check_binary("perl -e \'use Bio::DB::Fasta; print \"installed\";\'", "installed");
	if ($value){
		die "Error: \tBioPerl not found\nPlease check if bioperl was installed and its path was set correctly.\n ";
	}
	
	$value = check_binary("perl -e \'use Bio::DB::Sam; print \"installed\";\'", "installed");
	if ($value){
		die "Error: \tBio::DB::Sam module not found\nPlease check if Bio::DB::Sam was installed and its path was set correctly.\n";
	}
	
	$value = check_binary("bowtie --version 2>&1","version");
	if ($value){
		die "Error: \tbowtie not found\nPlease check if bowtie was installed and its path was set correctly.\n";
	}
	
	$value = check_binary("RNAfold --help 2>&1","gamma");
	if ($value){
		die "Error: \tRNAfold not found\nPlease check if RNAfold was installed and its path was set correctly.\n";
	}

	$value = check_binary("samtools 2>&1","faidx");
	if ($value){
		die "Error: \tsamtools not found\nPlease check if samtools was installed and its path was set correctly.\n";
	}
	
	return 0;

}

sub check_binary{
	my ($command,$feature) = @_;
	open (RET, "$command |") || return 0;
	my $value = 1;
	while (<RET>){
		if (/$feature/){
			$value = 0;
		}
	}
	close RET;
	return $value;
}

sub check_input_files{
	
    my @samples = split /,/, $reads;
	my @inputs;

    for my $i (0 .. $#samples){
		$log = "Formatting file $samples[$i] ... ";
		record_log($log);
        reads_sanity_check($samples[$i]);
        my $format = reads_format_check($samples[$i]);
        if ($format eq 'raw'){
            format_raw_reads($samples[$i],"$dir/input_${i}.fa");
        }elsif ($format eq 'std'){
            format_clean_reads($samples[$i],"$dir/input_${i}.fa");
        }else{
            die "Error occured when format file $samples[$i]\n";
        }
        push @inputs, "$dir/input_${i}.fa";
		$log = "done\n";
		record_log($log);
        $i ++;
	}
    $reads = join(",", @inputs);
	$log = "Combine all clean reads ... ";
	record_log($log);
    combine_clean_reads(\@inputs,$reads_collapsed);
	$log = "done\n";
	record_log($log);
	
	if ($known_miRNA){
		$log = "Checking mature miRNAs ... ";
		record_log($log);
		mature_check($known_miRNA);
		$log = "done\n";
		record_log($log);
	}
	
	$log = "Checking genome file ... ";
	record_log($log);
	genome_check($genome);
    my $chr_cnt = `grep -c "^>" $genome`;
    unless ($chr_cnt <= 5000){
        $log = "done\n";
		record_log($log);
        $log = "Transform input genome ...";
		record_log($log);
        my $ret_scaffold = system("transform_genome.pl -d $genome -i $dir/pseudo.index -t T -o $dir/pseudo_genome.fa 2>/dev/null");
        unless ($ret_scaffold == 0){
            die "Error occured when transforming genome file\nPlease check the input files.\n";
        }
        my $ret_faidx_pseudo_genome = system("samtools faidx $genome 2>/dev/null");
        unless ($ret_faidx_pseudo_genome == 0){
            die "Error occured when running samtools faidx on pseudo_genome.fa\nThis may be a bug, please let me [admin\@ncrna.net] know, thank you.\n";
        }
        $flag = 1;
    }else{
        my $ret_faidx_genome = system("samtools faidx $genome 2>/dev/null");
        unless ($ret_faidx_genome == 0){
            die "Error occured when running samtools faidx\nThis may be caused by sequence with different length.\nPlease transform the genome with the script fasta_record_with_same_length.pl.\n";
        }
    }
	$log = "done\n";
	record_log($log);
    return;
}

sub reads_sanity_check{
	my ($input) = @_;
	open FH, "<", $input or die "Cannot open file $input: $!\n";

	while (<FH>){
		chomp;
		if ($.==1){
			my $first = $_;
			if ($first !~ />\S+/){
				die "Error: the first line in file $input does not start with '>identifier'\n$hints";
			}
		}
		if (/^\>.+$/){
		}elsif (!/^[ATCGNatcgn]+$/ or /^\s$/){
				die "Error in line $.:\tthe sequence $_ contains characters other than [ATCGNatcgn]\n$hints";
		}
	}
	close FH;
    return;
}

sub reads_format_check{
    my ($file) = @_;
    my $flag;
    open my $fh, "<", $file or die "Cannot open file: $!\n";
    while (<$fh>){
        chomp;
        if (/^>\S+_x(\d+)/){
            $flag = 'std';
            return $flag;
        }else{
            $flag = 'raw';
            return $flag;
        }
    }
    close $fh;
    
    return;
}


sub format_raw_reads {
        my ($file,$out) = @_;
        my %hash;
        my ($id,$cnt,$desc,$seq,$lng) = ();
        open (FASTA, "<", $file) or die $!;
        while (<FASTA>) {
                chomp;
                if (/^>(\S+).*/){
                        $id     = $1;
                        $cnt    = 1;
                        $seq    = "";
                        while (<FASTA>) {
                                chomp;
                                if (/^>(\S+).*/) {
                                    $lng = length $seq;
                                    if ($lng > 17 && $lng < 31){
                                        if (defined $hash{$seq}){
                                            $hash{$seq} += $cnt;
                                        }else{
                                            $hash{$seq} = $cnt;
                                        }
                                    }
                                    $id     = $1;
                                    $cnt    = 1;
                                    $seq    = "";
                                    next;
                                }
                                $seq .= $_;
                        }
                }
        }
        $lng = length $seq;
        if ($lng > 17 && $lng < 31){
            if (defined $hash{$seq}){
                $hash{$seq} += $cnt;
            }else{
                $hash{$seq} = $cnt;
            }
        }
        close FASTA;

        open my $fh, ">", $out or die "Cannot create file $out: $!\n";
        my $i = 1;
        for my $seq (keys %hash){
            my $cnt = $hash{$seq};
            printf $fh ">t%08s", $i;
            print $fh "_x", $cnt, "\n";
            print $fh "$seq\n";
            $i ++;
        }
        close $fh;

        return;
}

sub format_clean_reads {
        my ($file,$out) = @_;
        my %hash;
        my ($id,$cnt,$desc,$seq,$lng) = ();
        open (FASTA, "<", $file) or die $!;
        while (<FASTA>) {
                chomp;
                if (/^>(\S+)_x(\d+)(.*)/){
                        $id     = $1;
                        $cnt    = $2;
                        $desc   = $3;
                        $seq    = "";
                        while (<FASTA>) {
                                chomp;
                                if (/^>(\S+)_x(\d+)(.*)/) {
                                    $lng = length $seq;
                                    if ($lng > 17 && $lng < 31){
                                        if (defined $hash{$seq}){
                                            $hash{$seq} += $cnt;
                                        }else{
                                            $hash{$seq} = $cnt;
                                        }
                                    }
                                    $id     = $1;
                                    $cnt    = $2;
                                    $desc   = $3;
                                    $seq    = "";
                                    next;
                                }
                                $seq .= $_;
                        }
                }
        }
        $lng = length $seq;
        if ($lng > 17 && $lng < 31){
            if (defined $hash{$seq}){
                $hash{$seq} += $cnt;
            }else{
                $hash{$seq} = $cnt;
            }
        }
        close FASTA;
        
        open my $fh, ">", $out or die "Cannot create file $out: $!\n";
        my $i = 1;
        for my $seq (keys %hash){
            my $cnt = $hash{$seq};
            printf $fh ">t%08s", $i;
            print $fh "_x", $cnt, "\n";
            print $fh "$seq\n";
            $i ++;
        }
        close $fh;
        return;
}

sub combine_clean_reads {
        my ($input,$out) = @_;
        my %hash;
        for my $file (@$input){ 
            my ($id,$cnt,$desc,$seq) = ();
            open my $fh, "<", $file or die $!;
            while (<$fh>) {
                chomp;
                if (/^>(\S+)_x(\d+)(.*)/){
                        $id     = $1;
                        $cnt    = $2;
                        $desc   = $3;
                        $seq    = "";
                        while (<$fh>) {
                                chomp;
                                if (/^>(\S+)_x(\d+)(.*)/) {
                                    if ($seq eq ""){
                                    }elsif (defined $hash{$seq}){
                                        $hash{$seq} += $cnt;
                                    }else{
                                        $hash{$seq} = $cnt;
                                    }
                                    $id     = $1;
                                    $cnt    = $2;
                                    $desc   = $3;
                                    $seq    = "";
                                    next;
                                }
                                $seq .= $_;
                        }
                }
        }
        if ($seq eq ""){
        }elsif (defined $hash{$seq}){
            $hash{$seq} += $cnt;
        }else{
            $hash{$seq} = $cnt;
        }
        close $fh;
        
    }
    open my $fh, ">", $out or die "Cannot create file $out: $!\n";
    my $i = 1;
    for my $seq (keys %hash){
        my $cnt = $hash{$seq};
        printf $fh ">t%08s", $i;
        print $fh "_x", $cnt, "\n";
        print $fh "$seq\n";
        $i ++;
    }
    close $fh;
    return;
}

sub genome_check{
	my ($input) = @_;
	my %hash_id;
	my $id;
	open FH, "<", $input or die "Cannot open file $input\n";

	while (<FH>){
		chomp;
		if ($.==1){
			my $first = $_;
			if ($first !~ />\S+/){
				die "Error: the first line in file $input does not start with '>identifier'\n$hints";
			}
		}
		if (/^>(\S+).*$/){
			$id = $1;
			if (defined $hash_id{$id}){
				die "Error in line $.: the identifier $id is not unique\n$hints";
			}else{
				$hash_id{$id}=1;
			}
		}elsif (!/^[A-Za-z]+$/ or /^\s$/){
				die "Error in line $.:\tthe sequence $_ contains unrecognised character or white space\n$hints";
		}
	}
	close FH;

    return;
}

sub mature_check{
	my ($input) = @_;
	my %hash_id;
	my $id;
	open FH, "<", $input or die "Cannot open file $input\n";

	while (<FH>){
		chomp;
		if ($.==1){
			my $first = $_;
			if ($first !~ />\S+/){
				die "Error: the first line in file $input does not start with '>identifier'\n$hints";
			}
		}
		if (/^>(\S+).*$/){
			$id = $1;
			if (defined $hash_id{$id}){
				die "Error in line $.: the identifier $id is not unique\n$hints";
			}else{
				$hash_id{$id}=1;
			}
		}elsif (!/^[ATCGNatcgn]+$/ or /^\s$/){
				die "Error in line $.:\tthe sequence $_ contains characters other than [ATCGNatcgn]\n$hints";
		}
	}
	close FH;
    return;
}

sub check_index{
    
    my @ebwt_suffix = qw(.rev.1.ebwt .rev.2.ebwt .1.ebwt .2.ebwt .3.ebwt .4.ebwt);
    for my $suffix (@ebwt_suffix){
        my $ebwt_file = $genome_index . $suffix;
        unless (-r $ebwt_file){
            $log = "done\n";
			record_log($log);
            $log = "$ebwt_file is not available ... done\n";
			record_log($log);
            $log = "Rebuilding the bowtie indices, this may take a few minutes ...";
			record_log($log);
            my $ret_build = system("bowtie-build $genome $genome_index 1>/dev/null 2>/dev/null");
            unless ($ret_build == 0){
                die "Error occured when indexing the genome with bowtie-build\nPlease check the genome file.\n";
            }
            $log = "done\n";
			record_log($log);
            return;
        }
    }
    $log = "done\n";
	record_log($log);
    $log = "Bowtie indices file is available.\n";
	record_log($log);
    return;
}

sub reads_map_genome{
    
	my $ret_mapping_bam = system("bowtie -f -a -v 0 -S -m $map_time $genome_index $reads_collapsed 2>/dev/null | samtools view -bhSF 4 - -o $dir/reads_mapped_genome.bam 2>/dev/null");

    unless ($ret_mapping_bam == 0){
        die "Error happened when mapping reads to the reference genome\nPlease check the input files.\n";
    }

    my $ret_sorted_bam = system("samtools sort $dir/reads_mapped_genome.bam $dir/reads_mapped_genome_sorted 2>/dev/null");

    unless ($ret_sorted_bam == 0){
        die "Error happened when sorting bam file with samtools\nPlease check the memory usage.\n";
    }

    my $ret_bam_indexed = system("samtools index $dir/reads_mapped_genome_sorted.bam 2>/dev/null");

    unless ($ret_bam_indexed == 0){
        die "Error happened when indexing the sorted bam file using samtools\nPleans check the sorted bam file.\n";
    }

    return;

}

sub excise_potential_precursors{

    # excise potential precursors
    my $ret_precursors;
    $ret_precursors = system("excise_potential_precursors.pl -b $sorted_bam -m $min_freq -l $hairpin_lng -r $genome -o $dir/precursors.fa 2>/dev/null");
    unless ($ret_precursors == 0){
        die "Error occured when excising potential precursors\nPlease check if the BioPerl and Bio::DB::Sam  modules are installed correctly.\n";
    }
    return;
}

sub prepare_second_structure{

    my $ret_structure = system("RNAfold_with_multi-threads.pl -i $dir/precursors.fa -p $nb_process -d $dir_struct -o $dir/precursors.struct 2>/dev/null");
    unless ($ret_structure == 0){
        die "Error occured when RNAfolding precursors\nPlease check RNAfold -noPS/--noPS.\n";
    }
    return;
}

sub mirna_prediction{

    my $ret_prediction;
	if ($known_miRNA){
		$ret_prediction = system("miR_island_core.pl -b $sorted_bam -c $cotyledon -s $structure -g $genome -r $output -t $table -m $known_miRNA 2>/dev/null");
	}else{
		$ret_prediction = system("miR_island_core.pl -b $sorted_bam -c $cotyledon -s $structure -g $genome -r $output -t $table 2>/dev/null");
	}

    unless ($ret_prediction == 0){
        die "Error occured in the key step when predicting microRNAs\nPlease check the input files.\n";
    }
	
    return;
}

sub final_check{
    
    my $ret_result = system("transform_genome.pl -d $dir_tmp/result.out -i $index -t F -o $dir/result.out 2>/dev/null");
    unless ($ret_result == 0){
        die "Error occured when transforming result file.\nPlease check if the command is correct.\n";
    }

    my $ret_table = system("transform_genome.pl -d $dir_tmp/result.table -i $index -t F -o $dir/result.table 2>/dev/null");
    unless ($ret_result == 0){
        die "Error occured when transforming table file\nPlease check if the command is correct.\n";
    }
    
    

    return;

}

sub diff_mature_analysis{

    unless ($known_miRNA){
        $log = "done\nSorry, mature miRNA is not available.\nKnown miRNA expression analysis is ignored.\n";
		record_log($log);
    }else{
        my $ret_mature_dge = system("diff_miRNA_expression_analysis.pl -m $known_miRNA -c $reads_collapsed -r $reads -o $known_diff 2>/dev/null");
        unless($ret_mature_dge == 0){
            die "Error occured when Calculating known miRNA differential gene expression\nPlease check the privilege and the input files.\n";
        }
        $log = "done\n\n";
        record_log($log);
    }
    return;
}

sub denovo_diff_miRNA_analysis{
    my $ret_denovo;
    if ($known_miRNA){
        $ret_denovo = system("diff_miRNA_expression_analysis.pl -m $known_miRNA -l $dir/result.table -r $reads -o $denovo_diff 2>/dev/null");
    }else{
        $ret_denovo = system("diff_miRNA_expression_analysis.pl -l $dir/result.table -r $reads -o $denovo_diff 2>/dev/null");
    }
    unless ($ret_denovo == 0){
        die "Error occured when doing differential denovo miRNA analysis\nPlease check if $dir/result.table is available.\nand each $reads is in the allowed format.\n";
    }

    return;
}

sub flag_gene_out{
    my ($in_file,$out_file)=@_;
    return unless defined $coord_file;
    my %hash_coord;
    
    open my $fh, "<$coord_file" or die "Can not open file $coord_file: $!\n";
	my ($subject,$strand,$beg,$end,$name);
	while (<$fh>){
		chomp;
		if (/^(\S+)\s+(\+|\-)\s+(\d+)\s+(\d+)\s+(\S+)/){
			$subject = $1;
			$strand = $2;
			$beg = $3;
			$end = $4;
			$name = $5;
			$hash_coord{$subject}{$strand}{$beg}{$end} = $name;
		}
	}
	close $fh;

	open my $fh_in, "<$in_file" or die "Can not open file $in_file: $!\n";
	open my $fh_out, ">$out_file" or die "Can not create file $out_file: $!\n";
	while (<$fh_in>){
		chomp;
		my ($primary,$subject,$strand,$beg,$end,$flag);
		if (/(.*\t)(\S+)\s+(\+|\-)\s+(\d+)\s+(\d+)$/ || /(microRNA name\s+)(\S+)\s+(\+|\-)\s+(\d+)\s+(\d+)$/){
			$primary = $1;
			$subject = $2;
			$strand = $3;
			$beg = $4;
			$end = $5;
			if (exists $hash_coord{$subject} && exists $hash_coord{$subject}{$strand}){
				for my $db_beg (sort {$a<=>$b} keys %{$hash_coord{$subject}{$strand}}){
					for my $db_end (sort {$a<=>$b} keys %{$hash_coord{$subject}{$strand}{$db_beg}}){
						if ($db_beg<=$beg && $end<=$db_end){
							$flag = $hash_coord{$subject}{$strand}{$db_beg}{$db_end};
						}
						last if defined $flag;
					}
					last if defined $flag;
				}			
			}
			if (defined $flag){
				print $fh_out "$primary$flag\n";
			}else{
				print $fh_out "$primary-\n";
			}
		}else{
			print $fh_out "$_\n";
		}
	}
	close $fh_in;
	close $fh_out;
    return;
}

sub make_dir_tmp{

    #make temporary directory
    if (not -d "$folder"){
        my $ret_dir = system("mkdir '$folder' 2>/dev/null");
        unless ($ret_dir==0){
            die "Error occured when creating directory $folder: $!\n";
        }
    }
    
    my $ret_mkdir = system("mkdir $dir 2>/dev/null");
    unless ($ret_mkdir==0){ die "Error occured when making directory $dir: $!\n"; }
    
    my $ret_mkdir_tmp = system("mkdir $dir_tmp 2>/dev/null");
    unless ($ret_mkdir_tmp==0){ die "Error occured when making directory $dir_tmp: $!\n";}
    
    my $ret_dir_tmp = system("mkdir $dir_struct 2>/dev/null");
    unless ($ret_dir_tmp==0){ die "Error occured when making directory $dir_struct: $!\n";}

    return;

}

sub check_cpu_cores{
    my $cpu_cores;
    my $ret_cpu_check = system("cat /proc/cpuinfo 1>/dev/null 2>/dev/null");
    unless ($ret_cpu_check==0){
        $log = "Cannot find CPU info with command 'cat /proc/cpuinfo'\nPlease check cpu cores with command 'top'\n";
		record_log($log);
        return;
    }
    open (CPU_INFO, "cat /proc/cpuinfo |") or die "Cannot find CPU cores with 'cat /proc/cpuinfo'\n";
    while (<CPU_INFO>){
        chomp;
        if (/cpu cores\s+:\s+(\d+)/){
            $cpu_cores = $1;
            last;
        }
    }
    unless (defined $nb_process){
        unless (defined $cpu_cores){
            $nb_process = 1;
        }else{
            $nb_process = $cpu_cores;
        }
    }else{
        if (defined $cpu_cores){
            if ($nb_process > $cpu_cores){
                $nb_process = $cpu_cores;
                $log = "WARNING: Your machine has $cpu_cores cpu cores, the parametre may be incorrect\n";
				record_log($log);
            }
        }
    }
    return;
}

sub record_log{
	my ($sentence) = @_;
	print STDERR $sentence;
	print LOG $sentence;
	return;
}

sub var_check{
    if (exists $opts{'v'}) {
        die "\nmiR-island version = $version\n\n";
        exit 0;
    }
    if (exists $opts{'h'}){
        &var_error();
    }
    if ($opts{'c'}){
        $cotyledon=$opts{'c'};
    }else{
        &var_error();
    }
    if ($opts{'c'}){
        $genome=$opts{'g'};
    }else{
        &var_error();
    }
    if ($opts{'r'}){
        $reads=$opts{'r'};
    }else{
        &var_error();
    }
    unless ($opts{'e'}){
	    $mfe=-17.5;
    }else{
	    $mfe=$opts{'e'};
    }
    unless ($opts{'f'}){
	    $min_freq=15;
    }else{
	    $min_freq=$opts{'f'};
    }
    unless ($opts{'i'}){
    	$genome_index=$genome;
    }else{
    	$genome_index=$opts{'i'};
    }
    if ($opts{'k'}){
        $known_miRNA = $opts{'k'};
    }
    unless ($opts{'m'}){
    	$map_time=10;
    }else{
    	$map_time=$opts{'m'};
    }
    if ($opts{'l'}){
        $hairpin_lng = $opts{'l'};
    }else{
        $hairpin_lng = 277;
    }
    if ($opts{'p'}){
    	$nb_process=$opts{'p'};
    }
    if ($opts{'x'}){
        $coord_file=$opts{'x'};
    }
    return;
}

sub var_error {

    print "\n";
    print " MiRNA_Island: an ultrafast package for annotation and quantification of miRNA and MIRNA\n";
    print " genes with High throughput sequencing\n\n";
    print " Version = $version\n\n";
    print " WARNING: You did not provide enough information!\n\n" unless exists $opts{'h'};
    print " Usage: miRNA_island.pl -c <monocot|dicot> -g <genome.fa> -r <reads1.fa[,...,readsN.fa]> \\\n";
    print "           [-i <bowtie_ebwt>] [-e <mfe>] [-f <min_freq>] [-m <max_hits>] [-k <miRNA.fa>] \\\n";
    print "           [-l <haipin_length>] [-p <num_threads>] [-x <pre-miRNA_coord>]\n";
    print "\n";
    print " REQUIRED:\n";
	print " -c <monocot|dicot>             cotyledon type: monocot | dicot\n";
	print " -g <genome.fa>                 reference genome file in multi-fasta format\n";
	print " -r <reads1.fa[,...,readsN.fa]> a comma-separated list of files with small RNA reads in\n";
    print "                                multi-fasta format";
    print "\n";
    print " OPTIONAL:\n";
    print " -e <min_mfe>     maximum threshold MFE to be annotated as a potential precursor (def: -17.5)\n";
	print " -f <int>         minimum frequency of reads to trigger a precursor excising (def: 15)\n";
	print " -i <bowtie_ebwt> prefix of the reference genome's bowtie indexes file (def: same as option \"g\")\n";
	print " -k <miRNA.fa>    known miRNAs of the related species in multi-fasta format\n";
    print " -l <int>         maximum length of the potential precursor to be excised (def: 277)\n";
    print " -m <int>         maximum number of sites that a read could map to related genome (def: 10)\n";
	print " -p <int>         number of threads to use (def: all threads available)\n";
    print " -x <coord>       known pre-miRNA coordinates in tab-delimited file (only seperate\n";
    print "                  the location and the pre-miRNA name with a TAB (such as):\n";
    print "                  Chr2 - 10676451 10676573 ath-MIR156a\n";
    print "                  Chr4 + 9888982 9889070   ath-MIR160b\n";
	print " -h               print intact help message\n";
    print "\n";
    exit 1;
}
########################################################################################################
