#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use threads;
use Thread::Queue;

my $version = "1.0";
my ($precursor, $nb_process, $tmp_struct_dir, $structure);
my %opts=();
getopts("hi:d:o:p:v",\%opts);
if (exists $opts{v}){ die "\nRNAfold_with_multi-threads.pl version = $version\n\n"; exit 0; }
var_error() if $opts{h};
var_error() unless ($opts{i} && $opts{d} && $opts{o});
&var_check();
########################################  MAIN  #############################################

my @data = split_fa($precursor,$tmp_struct_dir,$nb_process);

my $stream = Thread::Queue->new(@data,undef);
my $nb_mission = scalar @data;

my @running = ();
my @Threads;
while (scalar @Threads < $nb_mission) {
    @running = threads->list(threads::running);

    if (scalar @running < $nb_process) {
	    my $file = $stream->dequeue();
    	my $out = ${file}."_struct";
        my $thread = threads->new(\&secondary_structure_prediction,$file,$out);
        push (@Threads, $thread);
        my $tid = $thread->tid;
    }
    @running = threads->list(threads::running);
    foreach my $thr (@Threads) {
        if ($thr->is_running()) {
            my $tid = $thr->tid;
        }
        elsif ($thr->is_joinable()) {
            my $tid = $thr->tid;
            $thr->join;
        }
    }
 
    @running = threads->list(threads::running);
}

while (scalar @running != 0) {
    foreach my $thr (@Threads) {
        $thr->join if ($thr->is_joinable());
    }
    @running = threads->list(threads::running);
    sleep(3);
}

my $ret_cat = system("cat $tmp_struct_dir/*_struct > $structure 2>/dev/null");

unless ($ret_cat==0){
	die "Error occured when gathering the odd precursor_tmp.struct\nPlease check the previlege.\n";
}

exit 0;

############################################# SUBROUTINE ###############################################

sub split_fa{
	my ($filename,$tmpdir,$number) = @_;

	#check temp folder
	mkdir $tmpdir, 0755 or die "Cannot make temp directory: $!"
		unless -e $tmpdir;

	#open filehandles for output splitted files
	my @out_fhs;
	my @out_files;
	for (0..$number-1) {
		my $split_name = "$tmpdir/precursor.$_.tmp";
		open my $fh, ">", $split_name  or die "Cannot create file $split_name: $!";
		push @out_fhs,$fh;
		push @out_files, $split_name;
	}

	#read and split file
	open my $in_fh, "<", $filename or die "Cannot open file $filename: $!";
	my $seq_num;
	while (<$in_fh>) {
		$seq_num++ if /^>/;
		my $index = $seq_num % $number;
		my $fh = $out_fhs[$index];
		print $fh $_;
	}
	#colse filehandles
	for (0..$number-1) {
		close $out_fhs[$_];
	}
	#change to main directory in case of ...
	return @out_files;
}

sub secondary_structure_prediction{
	my ($file, $out) = @_;
	system("cat $file | RNAfold --noPS > $out");
}

sub var_check{
    if (exists $opts{v}){
        die "\nRNAfold_with_multi-threads.pl version = $version\n\n";
    }
    if (exists $opts{h}){
        &var_error();
    }
    if ($opts{i}){
        $precursor = $opts{i};
    }else{
        &var_error();
    }
    if ($opts{p}){
        $nb_process = $opts{p};
    }else{
        $nb_process = 1;
    }
    if ($opts{d}){
        $tmp_struct_dir = $opts{d};
    }else{
        &var_error();
    }
    if ($opts{o}){
        $structure = $opts{o};
    }else{
        &var_error();
    }
    return;
}

sub var_error {
    
    print "\n";
    print " RNAfold_with_multi-threads.pl computes potential precursors using RNAfold with\n";
    print " multiple threads\n\n";
    print " Version = $version\n\n";
    print " WARNING: You did not provide enough information!\n\n" unless exists $opts{h};
    print " Usage: RNAfold_with_multi-threads.pl -i <precursors.fa> -d <tmp_dir> \\\n";
    print "                                      -o <precursor.struct> [-p <num_threads>]\n";
    print " REQUIRED:\n";
    print " -i <precursors.fa>        a multi-fasta file with each sequence in a line\n";
    print " -d <tmp_directory>        a directory to store the temp files when running\n";
    print " -o <precursors.struct>    an RNAfold output file with sequences and structures\n";
    print "\n";
    print " OPTIONAL:\n";
    print " -p <num_threads>          number of threads to use\n";
    print "\n";
    exit 1;
}

############################################  END  ################################################
