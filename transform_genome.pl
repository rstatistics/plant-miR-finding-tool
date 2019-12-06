#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;

my %opts=();

getopt("d:o:i:t:hv",\%opts);
my $version = "1.0";
my ($data, $out, $index, $trans);
&var_check();

if ($trans eq 'T'){
	transform();
}elsif ($trans eq 'F'){
	restore();
}elsif ($trans eq 'O'){
    fasta_tidy();
}else{
	die "Option 't' $trans is not recognised, only T/F/O is allowed.\n";
}
exit 0;

######################################### SUBROUTINE #############################################

sub fasta_tidy{
    my $len = 60;
    my ($id, $seq, @processed);
    open FA, "<", $data or die "Cannot open file $data: $!\n";
    open FO, ">", $out or die "Cannot create file $out: $!\n";
    while(<FA>){
        chomp;
        if (/^>/){
            $id = $_;
            $seq = "";
        }
        while(<FA>){
            chomp;
            if (/^>/){
                @processed = $seq =~ /\w{$len}/g;
                print FO $id, "\n";
                print FO join ("\n",@processed), "\n";
                print FO $', "\n";
                $id = $_;
                $seq = "";
                next;
            }
            $seq .= $_;
        }
    }
    @processed = $seq =~ /\w{$len}/g;
    print FO $id, "\n";
    print FO join ("\n",@processed), "\n";
    print FO $', "\n";
    return;
}

sub transform{
	
	my $len = 60;
	my $num = 50;
	my $i = 0;
	my $scaffold_cnt = 0;
	open INDEX, ">", $index or die "Cannot create file $index: $!\n";
	open OUT, ">", $out or die "Cannot create file $out: $!\n";
	open DATA, "<", $data or die "Cannot open file $data: $!\n";
	my $scaffold_total = `grep -c "^>" $data`;
	my $scaffold_num = int($scaffold_total/($num-1));
	my ($id, $seq, $scaffold_beg, $scaffold_end, $seq_lng, $odd, @processed);
	while (<DATA>){
		chomp;
		if (/^>/){
			$id = $_;
			$id =~ s/^>(\S+).*/$1/;
			$seq = "";
		}
		while (<DATA>){
			chomp;
			if (/^>/){
				$scaffold_cnt ++;
				if ($i==0){
					$i ++;
					print OUT ">Chr${i}\n";
					$scaffold_beg = 1;
					$scaffold_end = 1;
					$odd = "";
					$seq_lng = 0;
				}
				$seq .= 'N' x 20;
				@processed = $seq =~ /\w{$len}/g;
				if (@processed != 0){
					print OUT join ("\n",@processed),"\n";
				}    
				if ($seq_lng==0){
					$seq_lng += length($seq);    
				}else{
					$seq_lng += length($seq) - length($odd);
				}
				$scaffold_end = $seq_lng;
				print INDEX "$id\tChr${i}\t$scaffold_beg\t$scaffold_end\n";
				$scaffold_beg = $scaffold_end + 1;
				$odd = $';
				if ($scaffold_cnt % $scaffold_num == 0){
					$i ++;
					print OUT $odd, "\n" unless $odd eq "";
					print OUT ">Chr${i}\n";
					$scaffold_beg = 1;
					$scaffold_end = 1;
					$odd = "";
					$seq_lng = 0;
				}			
				$id = $_;
				$id =~ s/^>(\S+).*/$1/;
				$seq = $odd;
				next;
			}
			$seq .= $_;
		}
	}
	$seq .= 'N' x 20;
	@processed = $seq =~ /\w{$len}/g;
	if (@processed != 0){
		print OUT join ("\n",@processed), "\n";
	}
	print OUT $';
	if ($seq_lng==0){
		$seq_lng += length($seq);
	}else{
		$seq_lng += length($seq) - length($odd);
	}
	$scaffold_end = $seq_lng;
	print INDEX "$id\tChr${i}\t$scaffold_beg\t$scaffold_end\n";

	close DATA;
	close INDEX;
	close OUT;
	return;
}

sub restore{
	my %hash_coord;
	open INDEX, "<", $index or die "Cannot open file $index: $!\n";
	while (my $coord=<INDEX>){
		if ($coord =~ /^(\S+)\t(\S+)\t(\d+)\t(\d+)/){
			my $scf_name = $1;
			my $subject = $2;
			my $scf_beg = $3;
			my $scf_end = $4;
			$hash_coord{$subject}{$scf_beg}{$scf_end} = $scf_name;
		}
	}
	close INDEX;

	open OUT, ">", $out or die "Cannot create file $out: $!\n";
	open DATA, "<", $data or die "Cannot open file $data: $!\n";
	my ($subject,$strand,$subject_beg,$subject_end);
	while (my $line=<DATA>){
		chomp $line;
		if ($line =~ /(Chr\d+)\s+(\S)\s+(\d+)\s+(\d+)/){
			$subject = $1;
			$strand = $2;
			$subject_beg = $3;
			$subject_end = $4;
			for my $beg (sort {$a<=>$b} keys %{$hash_coord{$subject}}){
				my $end = (keys %{$hash_coord{$subject}{$beg}})[0];
				if ($subject_beg >= $beg && $subject_end <= $end){
					my $scf_id = $hash_coord{$subject}{$beg}{$end};
					my $scf_beg = $subject_beg - $beg + 1;
					my $scf_end = $subject_end - $beg + 1;
					$line =~ s/Chr\d+\s+\S\s+\d+\s+\d+/$scf_id $strand $scf_beg $scf_end/;
                    if ($line =~ /(Chr\d+)\s+(\S)\s+(\d+)\s+(\d+)/){
            			$subject = $1;
            			$strand = $2;
            			$subject_beg = $3;
            			$subject_end = $4;
            			for my $beg (sort {$a<=>$b} keys %{$hash_coord{$subject}}){
            				my $end = (keys %{$hash_coord{$subject}{$beg}})[0];
            				if ($subject_beg >= $beg && $subject_end <= $end){
            					my $scf_id = $hash_coord{$subject}{$beg}{$end};
            					my $scf_beg = $subject_beg - $beg + 1;
            					my $scf_end = $subject_end - $beg + 1;
            					$line =~ s/Chr\d+\s+\S\s+\d+\s+\d+/$scf_id $strand $scf_beg $scf_end/;
                                last;
                            }
                        }
                    }
					last;
				}
			}
		    print OUT $line, "\n";
		}else{
			print OUT $line, "\n";
		}
	}
	close DATA;
	close OUT;
	return;
}

sub var_check{

    if (exists $opts{v}){
        die "\ntransform_genome.pl version = $version\n\n";
    }
    if (exists $opts{h}){
        &var_error();
    }
    if ($opts{d}){
        $data = $opts{d};
    }else{
        &var_error();
    }
    if ($opts{o}){
        $out = $opts{o};
    }else{
        &var_error();
    }
    if ($opts{t}){
        $trans = $opts{t};
    }else{
        &var_error();
    }
    if ($opts{i}){
        $index = $opts{i};
    }else{
        if ($opts{t} ne "O"){
            &var_error();
        }
    }
    return;
}

sub var_error{
    print "\n";
    print " transform_genome.pl transfoms scaffold-level genome into chromosome-level and vice versa;\n";
    print " it can also be used to transform the genome to multi-fasta with sequences in same length\n\n";
    print " WARNING: You did not provide enough information!\n\n" unless exists $opts{'h'};
    print " Usage: transform_genome.pl -d <genome.fa> -i <index_file> -t <T> -o <pseudo_genome.fa>\n";
    print " or:    transform_genome.pl -d <data> -i <index_file> -t <F> -o <data_processed>\n";
    print " or:    transform_genome.pl -d <genome.fa> -t <O> -o <genome_with_same_length.fa>\n";
    print "\n";
    print " OPTIONS:\n";
	print "\n";
	print " -d      data to transform or restroe\n";
	print " -i      index file that record scaffold info\n";
    print " -o      specify the output file\n";
    print " -t      function type to use: T/F/O\n";
	print "         T to transform scaffold-level genome to chromosome-level\n";
	print "         F to restore chromosome-level info to scaffold-level\n";
    print "         O to transform genome with same sequence length\n";
	print " -h      print intact help information\n";
    print "\n";
    exit 1;
}
