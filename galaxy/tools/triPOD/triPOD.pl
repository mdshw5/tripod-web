#!/usr/bin/perl
use strict;
use warnings;
use Algorithm::Cluster 'kcluster';
use Cwd 'getcwd';
use File::Path qw(make_path remove_tree);
use File::Temp qw(tempfile tempdir);
use Getopt::Long;
use POSIX qw(ceil floor);
use Set::IntervalTree;
use Tree::Interval;
use threads;
use Time::HiRes 'usleep';

my $CMD_LINE = join(" ", $0, @ARGV);
# Detect the total number of cores
my $MAX_CORES = `cat /proc/cpuinfo | grep processor | wc -l`;
my $THE_TIME = cpu_time();

# Get process ID, current working directory, and current time
my $CWD = getcwd();
my $PID_VALUE = $$;

# Default Parameters
my $ALPHA       = 0.05;
my $BATCH;
my $BUILD       = "hg18";
my $CORES       = $MAX_CORES - 1;
my $GENDER      = 0;
my $GRAPHICS    = "T";
my $HD          = "Y";
my $HELP;
my $NC_WARNING  = 0.02;
# Number of stdevs from the heterozygous BAF mean - used as boundaries 
# for detection of informative SNPs. (The intent is to ignore at least 50%
# of the members of the "normal" heterozygous SNP distribution.)
my $HET_SD      = 1.414213562373095; # sqrt(2)
# Number of stdevs from the homozygous BAF means - used as boundaries 
# for mBAF conversions. (The intent is to ignore nearly all 
# of the members of the "normal" homozygous SNP distribution.)
my $HOM_SD      = 4;
my $MI1         = "Y";
my $OUTPUT_DIR;
my $POD         = "Y";
my $PODcr       = "Y";
my $SIZE_REGION = 100;
my $STATS;
my $verbose     = "T";

# Ensure that the number of cores for analysis does not exceed 
# number of available cores
$CORES = $MAX_CORES unless ($MAX_CORES >= $CORES);
# Get optional input parameters (e.g. --cores=2 --bin=100)
GetOptions( 'alpha=f'  => \$ALPHA,		 'hetSD=f'  => \$HET_SD,
			'batch=s'  => \$BATCH,       'homSD=f'  => \$HOM_SD,
            'bin=i'    => \$SIZE_REGION, 'mi1!'     => \$MI1,
            'build=s'  => \$BUILD,       'out=s'    => \$OUTPUT_DIR,
            'cores=i'  => \$CORES,       'pod!'     => \$POD,
			'gender=s' => \$GENDER,		 'podcr!'   => \$PODcr,
            'graph!'   => \$GRAPHICS,    'stats!'   => \$STATS,
            'hd!'      => \$HD,   		 'verbose!' => \$verbose,
			'help!'    => \$HELP       
		  );

help() if $HELP;
my $CORES_PER_SAMPLE = $CORES;
my $AMP = 0.1;
my ($HD_LOW, $HD_MID, $HD_HI) = (-3, -2, -1.5);
my ($DEL_UPPER, $DEL_LOWER)   = (-0.1, -2);
# Array indices
use constant {
	CHR         => 1,	 POS     => 2,
	START       => 1,	 STOP    => 2,
	# Data type
	LRR         => 1,	 BAF     => 2,
	# Gender
	MALE        => 1,	 FEMALE  => 2,
	# Detection method
	POD         => 1,	 HD      => 2,
	MI1         => 3,	 PODcr   => 4,
	# ID and contribution
	FATHER      => 1,	 MOTHER  => 2,
	CHILD       => 3,	 NONE    => 4,
	BOTH    	=> 5,	 UNKNOWN => 6,
	UNKNOWN_MI1 => 7,    OUTLIER => 8,

	CHR_X       => 23,
	NA          => -1000,
	PI          => 4 * atan2(1, 1)
};

# Determine if gender was assigned
if ($GENDER =~ /^m/i) {$GENDER = MALE}
elsif ($GENDER =~ /^f/i) {$GENDER = FEMALE}
elsif ($GENDER) {
    $GENDER = 0; 
    print "The gender specification cannot be recognized. Please use M or F.",
        "\nProceeding without analysis of Chromosome X.\n\n";
}

# Check if batch mode is currectly submitted by user
my $check_batch;
if ($BATCH) {$check_batch = $BATCH}
else {$check_batch = $ARGV[-1]}

if(!open(INPUT_FILE, "<", $check_batch)) { 
	print STDERR "Could not open input file: $check_batch!\n";
	exit;
}
my $line = <INPUT_FILE>;
chomp($line);
my @line =  split("\t", $line);
close INPUT_FILE;
if (scalar(@line) == 1 && -e $line[0] && !$BATCH) {
	$BATCH = $check_batch; 
	print STDERR "The input appears to be a list of files. ",
		"Switching to BATCH mode.\n";
}
elsif (scalar(@line) == 12 && $BATCH) {
	push(@ARGV, $check_batch);
	$BATCH = 0;
	print STDERR "The input does not appear to be a list of files. ",
		"Switching to single sample mode.\n";
}

# Open batch file if supplied
my (@files, @batch, @batches);
if ($BATCH) {
	$verbose = "";
    open(BATCH_FILE, "<", $BATCH) || die("Could not open batch file!\n");
    while (<BATCH_FILE>) {
        my $file = $_;
        chomp($file);
        push(@files, $file);
    }
	close BATCH_FILE;
	if (@files >= $CORES) {
		$CORES_PER_SAMPLE = 1;
		# samples per core - 
		# Scenario - for 21 samples, 20 cores the program will calculate
		# that the max numer of samples to be run on any core is 2, and 
		# instead of distributing 2 samples to 1 core and singles to 19
		# cores, for efficiency of system resources and without losing
		# runtime, it will distribute 2 samples to 10 cores and the 
		# remaining samples to the 11th core.
		my $max_spc = ceil(@files / $CORES);
		for (my $i = 0; $i < @files; $i++) {
			if (!$i) {push(@batch, $files[$i])}
			elsif ($i % $max_spc) {push(@batch, $files[$i])}
			else {
				push(@batches, [@batch]);
				@batch = ();
				push(@batch, $files[$i]);
			}
		}
		push(@batches, [@batch]) if @batch;	
	}
	else {push(@batches, \@files)}
}
else {
    if ($ARGV[-1]) {
		$batches[0] = [$ARGV[-1]];
		push(@files, $ARGV[-1]);
	}
    else {die("Please include an input file.\n")}
}

# Open genome build file
my @CENTRO_REFS;
#my $BUILD_FILE = join("", "genome_build/", $BUILD, "_centromeres.txt");
my $BUILD_FILE = $BUILD;
open(BUILD_FILE, "<", $BUILD_FILE)
    || die("\nCould not open genome build file!\n",
    "Please ensure that the genome_build folder is located in the same ",
    "directory as triPOD.pl.\n");
    
<BUILD_FILE> for 1..2;
while (<BUILD_FILE>) {
    my $line = $_;
    chomp($line);
    my @line_array = split(/\t/, $line);
    $line_array[1] = (split(/chr/, $line_array[1]))[1];
    
    # Check Chromosome column for anything other than acceptable
    # chromosomes
    if ($line_array[1] !~ /^[1-9]$|^[1][0-9]|^[2][0-2]$|x|y|xy|m|mt/i) {
        print "\nUnrecognized character in Chromosome column of genome ",
        "build file!\n";
        print_proper_format();
        exit();
    }
    # Check Position columns for non-numeric characters and negative
    # numbers  
    elsif ($line_array[2] !~ /^-?\d/ || $line_array[2] < 0) {
        print "\nUnrecognized character in Position Start column of genome ",
        "build file!\n";
        print_proper_format();
        exit();
    }
    elsif ($line_array[3] !~ /^-?\d/ || $line_array[3] < 0) {
        print "\nUnrecognized character in Position End column of genome ",
        "build file!\n";
        print_proper_format();
        exit();
    }
    
    if ($line_array[CHR] !~ /^\d/) {
        if ($line_array[CHR] eq "X" && $GENDER) {$line_array[CHR] = CHR_X}
        else {next}
    }
    $CENTRO_REFS[$line_array[1]] = [@line_array[1..3]];
}

my $input_name;
if ($BATCH) {
	$input_name = (split(/[\/]/,(split(/.txt/, $BATCH))[0]))[-1];
}
else {$input_name = (split(/[\/]/,(split(/.txt/, $ARGV[-1]))[0]))[-1]}

if (!$OUTPUT_DIR) {$OUTPUT_DIR = "$CWD/triPOD_Results/$input_name"}
# If user input contains a trailing backslash, remove it
else {$OUTPUT_DIR =~ s|/\z||}
#eval {make_path("$OUTPUT_DIR/log_files")} 
make_path("$OUTPUT_DIR/log_files");
-e "$OUTPUT_DIR/log_files" || die "Could not create output directory: $OUTPUT_DIR\n";
print "Output directory = $OUTPUT_DIR/\n";
print "Number of Samples = ", scalar(@files), "\n";

my $OUTPUT_FILE;
if ($BATCH) {$OUTPUT_FILE = "$OUTPUT_DIR/$input_name\_triPOD_Batch_Results.txt"}
else {$OUTPUT_FILE = "$OUTPUT_DIR/triPOD_Results.txt"}
open(OUTPUT_FILE, '>', $OUTPUT_FILE) || die 
	"Could not create output file - $OUTPUT_FILE!\n"; 
print OUTPUT_FILE join("\n", $THE_TIME, $CMD_LINE), "\n\n";

# Redirect STDERR to a log file
my $PERL_LOG = "$OUTPUT_DIR/log_files/$input_name\_perl_log.txt";
open(STDERR, '>', $PERL_LOG) || die 
	"Could not create perl log file - $PERL_LOG!\n"; 

# Create stats file if requested
my $STATS_FILE;
if ($STATS) {
	$STATS_FILE = "$OUTPUT_DIR/$input_name\_stats.txt";
	open(STATS_FILE, '>', $STATS_FILE) || die 
		"Could not create stats file - $STATS_FILE!\n"; 
	print STATS_FILE "File\tChild\tFather\tMother\tChild_NC_rate\t",
	"Father_NC_rate\tMother_NC_rate\tChild_HD_rate\tChild_Min_HD\t",
	"Father_HD_rate\tFather_Min_HD\tMother_HD_rate\tMother_Min_HD\t",
	"MI1_rate\tMin_MI1\tMI1_Upper_Thresh\tMI1_Lower_Thresh\tMin_POD\t",
	"Prob\t<1per?Genomes\tError_rate\tAcceptable_Errors\tAA_Bound\t",
	"BB_BOUND\tAB_UPPER_BOUND\tAB_LOWER_BOUND\tNum_Detected_Regions\t",
	"Num_POD_regions\tNum_PODcr_regions\tNum_MI1_regions\tNum_HD_regions\n";
}

# Declare global variables
my ($AA_BOUND, $AB_LOWER_BOUND, $AB_UPPER_BOUND, $ACCEPTABLE_ERRORS, 
    $ADJUSTED_MI1, $BB_BOUND, $BOUNDARY_EXTENSION, $CH_LRR_MED, 
	$CH_mBAF_MED, $CH_NAME, $ERROR_RATE, $FILENAME, $HET_MI1_RATE, $INF_SNPS, $INPUT_FILE, 
	$MIN_BAF_OUT, $MIN_POD, $LOCAL_CH_BAF_MEAN, $LOCAL_CH_mBAF_MED, 
	$LOCAL_CH_LRR_MED, $LOCAL_P1_mBAF_MED, $LOCAL_P1_LRR_MED, $LOCAL_P2_mBAF_MED,
	$LOCAL_P2_LRR_MED, $MI1_LOWER_THRESH, $MI1_UPPER_THRESH, $MIN_MI1, $MIN_P1_HD,
	$MIN_P2_HD, $MIN_CH_HD, $P1_AA_BOUND, $P1_BB_BOUND, $P2_AA_BOUND, $P2_BB_BOUND,
	$P1_LRR_MED, $P1_mBAF_MED, $P1_NAME, $P2_LRR_MED, $P2_mBAF_MED, $P2_NAME,
	$PERL_TO_R_FILE, $PERL_TO_R_FILENAME, $POD_ALPHA, $PODCR_LOWER_THRESH, 
	$PODCR_UPPER_THRESH, $R_ATTEMPTS, $R_LOG_FILENAME, $R_PID, $R_PID_FILE, 
	$R_PID_FILENAME, $R_TO_PERL_FILE, $R_TO_PERL_FILENAME, $REFINING);
my (@BAF_OUTLIERS, @CURR_CHR_REFS, @INIT_CALC_REFS, @INIT_CH_STATS, @INIT_P1_STATS, 
	@INIT_P2_STATS, @LOCAL_STATS, @LRR_REFS, @NONOVERLAP_RESULTS,  
	@P1_INF_SNPS, @P2_INF_SNPS, @STREAK_ARRAY);
my (%BAF_BY_POS, %BAF_OUTLIERS_BY_POS, %CH_HD_BY_POS, %LRR_BY_POS, %MI1_BY_POS,  
	%P1_HD_BY_POS, %P1_INF_SNP_BY_POS, %P2_HD_BY_POS, %P2_INF_SNP_BY_POS, %SNP_BY_NUM, 
	%SNP_BY_POS);
my ($BATCH_REF, @batch_threads, $main_results);
# Count of completed samples - [0] = with abn, [1] = no abn
my @completed;
if ($#batches) {
	foreach my $batch_ref (@batches[0..$#batches-1]) {
		$BATCH_REF = $batch_ref;
		my $batch_thread = threads->new(\&main_program);
		push(@batch_threads, $batch_thread);
	}
}
$BATCH_REF = $batches[-1];
my (@output_array, @stats_array, @trio_def);
$main_results = main_program();
push(@output_array, @{$main_results->[0]}) if $main_results->[0];
push(@stats_array, @{$main_results->[1]}) if $STATS;
push(@trio_def, @{$main_results->[2]});
$completed[0] += $main_results->[3]->[0];
$completed[1] += $main_results->[3]->[1];

if (@batches > 1) {
	# Check if threads are finished
	while (threads->list(threads::running)) {sleep(2)}
	for (0..$#batch_threads) { 
		if ($batch_threads[$_]->is_joinable()) {
			$main_results = $batch_threads[$_]->join();
			push(@output_array, @{$main_results->[0]}) if $main_results->[0];
			push(@stats_array, @{$main_results->[1]}) if $STATS;
			push(@trio_def, @{$main_results->[2]});
			$completed[0] += $main_results->[3]->[0];
			$completed[1] += $main_results->[3]->[1];
		}
	}
}

my @sorted_output_array = sort {
	$$a[0] cmp $$b[0] || $$a[1] <=> $$b[1] || $$a[2] <=> $$b[2]
	} @output_array;
	
my @sorted_stats_array = sort {$$a[0] cmp $$b[0]} @stats_array;

print OUTPUT_FILE join("\n", @trio_def), "\n\n" unless $BATCH;
print OUTPUT_FILE "Detected Chromosomal Abnormalities", "\n",
	join("\t", qw(Sample Chr Start Stop Type Parent_of_Origin 
		Inheritance Size(SNPs) Informative_SNPs Size(bp) 
		Detection Median_mBAF Median_LRR Father-Median_mBAF 
		Father-Median_LRR Mother-Median_mBAF Mother-Median_LRR)
	),"\n";
map {print OUTPUT_FILE join("\t", @{$_}), "\n"} @sorted_output_array;
map {print STATS_FILE join("\t", @{$_}), "\n"} @sorted_stats_array;

$GENDER ||= "NA";			
my @tot_param;
push(@tot_param, "alpha=$ALPHA");
push(@tot_param, "batch=$BATCH") if $BATCH;
push(@tot_param, "bin=$SIZE_REGION", "build=$BUILD", "cores=$CORES",
	"gender=$GENDER");
if ($GRAPHICS) {push(@tot_param, "graph")}
else {push(@tot_param, "nograph")}
if ($HD) {push(@tot_param, "hd")}
else {push(@tot_param, "nohd")}
push(@tot_param, "hetSD=$HET_SD", "homSD=$HOM_SD");
if ($MI1) {push(@tot_param, "mi1")}
else {push(@tot_param, "nomi1")}
push(@tot_param, "out=$OUTPUT_DIR/"); 
if ($POD) {push(@tot_param, "pod")}
else {push(@tot_param, "nopod")}
print OUTPUT_FILE "\n\nPARAMETERS:\n--", join(" --", @tot_param);
close OUTPUT_FILE;     
close STATS_FILE if $STATS; 

my $temp = join("_", "$OUTPUT_DIR\/log_files\/temp.txt");
my $cleanup_script = "grep -v 'Wanted\\|Parameter' $PERL_LOG > $temp\n
	mv $temp $PERL_LOG";
system($cleanup_script);
my $successful = $completed[0] + $completed[1];
print "Samples successfully analyzed = $successful\n";
print "Number of failed samples = ", scalar(@files) - $successful, "\n"
	if (scalar(@files) - $successful);
print "Samples with no detectable abnormalities = $completed[1]\n"
	if ($completed[1]);
print STDERR "Samples successfully analyzed = $successful\n";
print STDERR "Number of failed samples = ", 
	scalar(@files) - $successful, "\n";
print STDERR "Samples with no detectable abnormalities = $completed[1]\n";
close STDERR; 

#THE END



sub main_program {
	############ Iterate analysis for number of files supplied #################
	my (@output_array, @stats_array, @trio_def);
	my @completed = (0) x 2;
	FILE:
	foreach my $file (@$BATCH_REF) {  
		print "starting $file\n";
		$INF_SNPS = 0;
#		$sample_count++;
#		print "\nSAMPLE $sample_count: $file\n";
		
		($AA_BOUND, $AB_LOWER_BOUND, $AB_UPPER_BOUND, $ACCEPTABLE_ERRORS, 
		 $ADJUSTED_MI1, $BB_BOUND, $BOUNDARY_EXTENSION, $CH_LRR_MED,
		 $CH_mBAF_MED, $CH_NAME, $ERROR_RATE, $FILENAME, $HET_MI1_RATE, $INF_SNPS, $INPUT_FILE, 
		 $MIN_BAF_OUT, $MIN_POD, $LOCAL_CH_BAF_MEAN, $LOCAL_CH_mBAF_MED, 
		 $LOCAL_CH_LRR_MED, $LOCAL_P1_mBAF_MED, $LOCAL_P1_LRR_MED, 
		 $LOCAL_P2_mBAF_MED, $LOCAL_P2_LRR_MED, $MI1_LOWER_THRESH, 
		 $MI1_UPPER_THRESH, $MIN_MI1, $MIN_P1_HD, $MIN_P2_HD, $MIN_CH_HD, 
		 $P1_AA_BOUND, $P1_BB_BOUND, $P2_AA_BOUND, $P2_BB_BOUND, $P1_LRR_MED,
		 $P1_mBAF_MED, $P1_NAME, $P2_LRR_MED, $P2_mBAF_MED, $P2_NAME, 
		 $PERL_TO_R_FILE, $PERL_TO_R_FILENAME, $POD_ALPHA, $PODCR_LOWER_THRESH, 
		 $PODCR_UPPER_THRESH,$R_ATTEMPTS, $R_LOG_FILENAME,
		 $R_PID, $R_PID_FILE, $R_PID_FILENAME, $R_TO_PERL_FILE, 
		 $R_TO_PERL_FILENAME, $REFINING) = (0) x 53;
		(@BAF_OUTLIERS, @CURR_CHR_REFS, @INIT_CALC_REFS, @INIT_CH_STATS, @INIT_P1_STATS, 
		 @INIT_P2_STATS, @LOCAL_STATS, @LRR_REFS, @NONOVERLAP_RESULTS, 
		 @P1_INF_SNPS, @P2_INF_SNPS, @STREAK_ARRAY) = (()) x 12;
		(%BAF_BY_POS, %BAF_OUTLIERS_BY_POS, %CH_HD_BY_POS, %LRR_BY_POS, %MI1_BY_POS, 
		 %P1_HD_BY_POS, %P1_INF_SNP_BY_POS, %P2_HD_BY_POS, %P2_INF_SNP_BY_POS, 
		 %SNP_BY_NUM, %SNP_BY_POS) = (()) x 11;
		
		
		# Input file, tab delim, sorted by chromosome and position, 
		# format = SNP Name, Chromosome, Position, Father, Father BAF, Father LRR, 
		# Mother, Mother BAF, Mother LRR, Child, Child BAF, Child LRR
		$INPUT_FILE = $file;
		if(!open(INPUT_FILE, "<", $INPUT_FILE)) { 
			print "Could not open input file: $INPUT_FILE!\n";
			print STDERR "Could not open input file: $INPUT_FILE!\n";			
			next FILE;
		}
		
		# Get filename
		$FILENAME = (split(/[\/]/,(split(/.txt/, $file))[0]))[-1];

		# Get sample names
		my $HEADERS = <INPUT_FILE>;
		$P1_NAME = (split(/[.]/, (split(/\t/, $HEADERS))[3]))[0];
		$P2_NAME = (split(/[.]/, (split(/\t/, $HEADERS))[6]))[0];
		$CH_NAME = (split(/[.]/, (split(/\t/, $HEADERS))[9]))[0];

		# Create temp files for communication with R and call Rscript    
		if ($GRAPHICS) {
			$R_LOG_FILENAME = join("", $OUTPUT_DIR, "/log_files/R_log.txt");
			($R_PID_FILE, $R_PID_FILENAME)         = tempfile(UNLINK => 1);
			($R_TO_PERL_FILE, $R_TO_PERL_FILENAME) = tempfile(UNLINK => 1);
			($PERL_TO_R_FILE, $PERL_TO_R_FILENAME) = tempfile(UNLINK => 1);
			print $PERL_TO_R_FILENAME, "\n";
			start_R();
			$R_ATTEMPTS = 0;
		}
		############################## Initial Calculations #########################    
		print "Performing Initial Calculations...\n" if $verbose;
		
		my $check = 0;
		my ($curr_chr, $thread_list, $thread_count, $prev_position) = (0) x 4;
		my ($centro_start, $centro_stop, $gap);
		my (@init_sums, @threads, @chromosomes);
				
		while (<INPUT_FILE>) {
			my $line = $_;
			chomp($line);
			my @line_array = split(/\t/, $line);
			next if (!$line_array[0]);
			
			#Setup
			if ($. > 1) {
				# Convert chromosome identifier to numeric 
				if ($line_array[CHR] !~ /^\d/) {
					if ($line_array[CHR] eq "X" && $GENDER) {
						$line_array[CHR] = CHR_X;
					}
					else {next}
				}
				$curr_chr = $line_array[CHR];
				$chromosomes[$curr_chr] = $curr_chr;
				# Store centromere info
				($centro_start, $centro_stop) = @{$CENTRO_REFS[$curr_chr]}[1,2];
				$gap = 1 if ($line_array[POS] <  $centro_stop);
				push(@INIT_CALC_REFS, \@line_array); 
				last;
			}
		}        

		while (<INPUT_FILE>) {			
			my $line = $_;
			chomp($line);
			my @line_array = split(/\t/, $line);
			
			next if (!$line_array[0]);
			# Convert chromosome identifier to numeric 
			if ($line_array[CHR] !~ /^\d/) {
				if ($line_array[CHR] eq "X" && $GENDER) {
					$line_array[CHR] = CHR_X
				}
				else {next}
			}
			if ($line_array[POS] > $centro_start 
				&& $line_array[POS] < $centro_stop) {
					# Check for correct genome build
					print "\n\nThe genome build does not match the input data!\n",
						"Please indicate the correct UCSC genome build ",
						"(e.g. --build=hg19).\n\n";
					next FILE;
			}
			###################### Input formatting checks #######################
			if ($. <= 1000) {
				# Check Chromosome column for anything other than acceptable
				# chromosomes
				if ($line_array[CHR] 
					!~ /^[1-9]$|^[1][0-9]|^[2][0-2]$|x|y|xy|m|mt/i) {
					print "\nSample $FILENAME : FAILED - Unrecognized ",
						"character ($line_array[CHR]) in Chromosome column!\n";
					print STDERR "Sample $FILENAME : FAILED - Unrecognized ",
						"character ($line_array[CHR]) in Chromosome column!\n";
					print_proper_format();
					next FILE;
				}
				# Check Position column for non-numeric characters and negative
				# numbers  
				elsif ($line_array[POS] !~ /^-?\d/ || $line_array[POS] < 0) {
					print "\nSample $FILENAME : FAILED - Unrecognized number ",
						"($line_array[POS]) in Position column!\n";
					print STDERR "Sample $FILENAME : FAILED - Unrecognized ",
						"number ($line_array[POS]) in Position column!\n";					
					print_proper_format();
					next FILE;
				}
				# Check Genotype columns for proper genotype format
				elsif ($line_array[3] !~ /AA|AB|BB|^N/ || 
					   $line_array[6] !~ /AA|AB|BB|^N/ || 
					   $line_array[9] !~ /AA|AB|BB|^N/) {
					print "\nSample $FILENAME : FAILED - Unrecognized ",
						"character (@line_array[3,6,9]) in a Genotype column!\n";
					print STDERR "Sample $FILENAME : FAILED - Unrecognized ",
						"character (@line_array[3,6,9]) in a Genotype ",
						"column!\n";					
					print_proper_format();
					next FILE;
				}
				# Check BAF columns for numbers < 0 or > 1
				elsif ((($line_array[4]  < 0 || ($line_array[4] > 1
					   && $line_array[4] != 2) || ($line_array[4] == 2
					   && substr($line_array[0], 0, 2) !~ /cn/i)) 
					   && $line_array[4] ne "NaN") ||
					   (($line_array[7]  < 0 || ($line_array[7] > 1
					   && $line_array[7] != 2) || ($line_array[7] == 2
					   && substr($line_array[0], 0, 2) !~ /cn/i)) 
					   && $line_array[7] ne "NaN") ||
					   (($line_array[10]  < 0 || ($line_array[10] > 1
					   && $line_array[10] != 2) || ($line_array[10] == 2
					   && substr($line_array[0], 0, 2) !~ /cn/i)) 
					   && $line_array[10] ne "NaN")) {
					print "\nSample $FILENAME : FAILED - B allele frequency ",
						"is out of range (@line_array[4,7,10]) !\n";
					print STDERR "Sample $FILENAME : FAILED - B allele ",
						"frequency is out of range (@line_array[4,7,10]) !\n";					
					print_proper_format();
					next FILE;
				}
				# Check for input in last expected column
				elsif (!defined($line_array[11])) {
					print "\nSample $FILENAME : FAILED - Too few columns.\n";
					print STDERR "Sample $FILENAME : FAILED - Too few ",
						"columns.\n";					
					print_proper_format();
					next FILE;
				}
				# Check if Position column is properly sorted
				elsif ($prev_position > $line_array[POS] 
						&& $curr_chr == $line_array[CHR]) {
					print "\nSample $FILENAME : FAILED - The input file ",
						"is not properly sorted!\n";
					print STDERR "Sample $FILENAME : FAILED - The input ",
						"file is not properly sorted!\n";					
					print_proper_format();
					next FILE;
				}
				$prev_position = $line_array[POS];
			} 
			#######################################################################    

			# Start new thread when centromere or new chromosome is encountered                
			if ($curr_chr != $line_array[CHR]
				|| ($gap && $line_array[POS] > $centro_start)) { 
				$gap = 0;
				if (@INIT_CALC_REFS >= $SIZE_REGION) {
					if ($CORES_PER_SAMPLE == 1) {
						my $results = initial_calculations();
						if ($results) {
							push(@INIT_CH_STATS, $$results[0]);
							push(@INIT_P1_STATS, $$results[1]);
							push(@INIT_P2_STATS, $$results[2]);
						}
					}
					else {    
						while (threads->list(threads::running)
							>= $CORES_PER_SAMPLE) {usleep(100000)}
						join_init_threads(\@threads, $thread_count);     
						$thread_list = threads->new(\&initial_calculations);
						push(@threads, $thread_list);
						$thread_count++;
					} 
				}
				if ($curr_chr != $line_array[CHR]) {
					$curr_chr = $line_array[CHR];			  
					# Store centromere info
					($centro_start, $centro_stop) =
						@{$CENTRO_REFS[$curr_chr]}[1,2];
					$gap = 1 if ($line_array[POS] <  $centro_stop);
					# Check for proper chromosome sorting
					if ($chromosomes[$curr_chr]) {
						print "\nThe input file is not properly sorted!\n";
						print_proper_format();
						next FILE;
					}
					else {$chromosomes[$curr_chr] = $curr_chr}
				}
				@INIT_CALC_REFS = ();
			}
			push(@INIT_CALC_REFS, \@line_array); 
		}
		close(INPUT_FILE);
		if (@INIT_CALC_REFS >= $SIZE_REGION) {
			if ($CORES_PER_SAMPLE == 1) {
				my $results = initial_calculations();
				if ($results) {
					push(@INIT_CH_STATS, $$results[0]);
					push(@INIT_P1_STATS, $$results[1]);
					push(@INIT_P2_STATS, $$results[2]);
				}
			}
			else {            
				while (threads->list(threads::running) >= $CORES_PER_SAMPLE) {
					usleep(100000);
				}
				join_init_threads(\@threads, $thread_count);   
				$thread_list = threads->new(\&initial_calculations);
				push(@threads, $thread_list);
				$thread_count++;
			} 
		}
		if ($CORES_PER_SAMPLE > 1) {
			# Check if threads are finished
			while (threads->list(threads::running)) {usleep(100000)}
			join_init_threads(\@threads, $thread_count);     
		}
		@INIT_CALC_REFS = ();  

		my @sorted_init_p1 = sort {
		$$a[0] <=> $$b[0] || $$a[1] <=> $$b[1]} @INIT_P1_STATS;
		my @sorted_init_p2 = sort {
		$$a[0] <=> $$b[0] || $$a[1] <=> $$b[1]} @INIT_P2_STATS;
		my @sorted_init_ch = sort {
		$$a[0] <=> $$b[0] || $$a[1] <=> $$b[1]} @INIT_CH_STATS;

		(my $p1_normals, $P1_mBAF_MED, $P1_LRR_MED, my $p1_nc_rate, 
		 $P1_AA_BOUND, $P1_BB_BOUND) = calc_init_stats(\@sorted_init_p1);

		(my $p2_normals, $P2_mBAF_MED, $P2_LRR_MED, my $p2_nc_rate, 
		 $P2_AA_BOUND, $P2_BB_BOUND) = calc_init_stats(\@sorted_init_p2);

		(my $ch_normals, $CH_mBAF_MED, $CH_LRR_MED, my $ch_nc_rate, 
		 $AA_BOUND, $BB_BOUND, $AB_UPPER_BOUND, $AB_LOWER_BOUND)
		 = calc_init_stats(\@sorted_init_ch);
		
		my $rounded_AA_bound     = sprintf("%.3f",$AA_BOUND);
		my $rounded_BB_bound     = sprintf("%.3f",$BB_BOUND);
		my $rounded_AB_up_bound  = sprintf("%.3f",$AB_UPPER_BOUND);
		my $rounded_AB_low_bound = sprintf("%.3f",$AB_LOWER_BOUND);
		
		# Calculate MI1 thresholds
		my ($ab_sums, $ab_sum_sq, $ab_ct) = (0) x 3;
		foreach (@$ch_normals) {
			$ab_sums   += $$_[2]->[3];
			$ab_sum_sq += $$_[2]->[4];
			$ab_ct     += $$_[2]->[5];
		}	
		my @AB_stats = st_dev($ab_sums, $ab_sum_sq, $ab_ct);
		$MI1_UPPER_THRESH = $AB_stats[0] + ($AB_stats[1] * 5);
		$MI1_UPPER_THRESH = 0.95 if $MI1_UPPER_THRESH > 0.95;
		$MI1_LOWER_THRESH = $AB_stats[0] - ($AB_stats[1] * 5);
		$MI1_LOWER_THRESH = 0.05 if $MI1_LOWER_THRESH < 0.05;
		$PODCR_UPPER_THRESH = $AB_stats[0] + ($AB_stats[1] * 2);
		$PODCR_LOWER_THRESH = $AB_stats[0] - ($AB_stats[1] * 2);

		@threads = ();
		($thread_count, $thread_list, $curr_chr) = (0) x 3;
		
		###### Calculate minimum size of region for std POD algorithm ########
		
		if(!open(INPUT_FILE, "<", $INPUT_FILE)) { 
			print "Could not open input file: $INPUT_FILE!\n";
			print STDERR "Could not open input file: $INPUT_FILE!\n";			
			next FILE;
		}

		while (<INPUT_FILE>) {
			my $line = $_;
			chomp($line);
			my @line_array = split(/\t/, $line);
			next if (!$line_array[0]);
			#Setup
			if ($. > 1) {
				# Convert chromosome identifier to numeric 
				if ($line_array[CHR] !~ /^\d/) {
					if ($line_array[CHR] eq "X" && $GENDER) {$line_array[CHR] = CHR_X}
					else {next}
				}
				$curr_chr = $line_array[CHR];
				
				# Store centromere info
				($centro_start, $centro_stop) = @{$CENTRO_REFS[$curr_chr]}[1,2];
				$gap = 1 if ($line_array[POS] <  $centro_stop);
				push(@CURR_CHR_REFS, \@line_array); 
				last;
			}
		}        

		while (<INPUT_FILE>) {
			my $line = $_;
			chomp($line);
			my @line_array = split(/\t/, $line);
			next if (!$line_array[0]);
			
			# Convert chromosome identifier to numeric 
			if ($line_array[CHR] !~ /^\d/) {
				if ($line_array[CHR] eq "X" && $GENDER) {$line_array[CHR] = CHR_X}
				else {next}
			}

			if ($line_array[POS] > $centro_start 
				&& $line_array[POS] < $centro_stop) {
					# Check for correct genome build
					print "\n\nThe genome build does not match the input data!\n",
						"Please indicate the correct UCSC genome build ",
						"(e.g. --build=hg19).\n\n";
					next FILE;
			}

			# Start new thread when centromere or new chromosome is encountered                
			if ($curr_chr != $line_array[CHR]
				|| ($gap && $line_array[POS] > $centro_start)) { 
				$gap = 0;
				
				if (@CURR_CHR_REFS >= $SIZE_REGION) {
					# Create threads only if >1 processor is available                
					if ($CORES_PER_SAMPLE == 1) {
						my $results = count_nonoverlap_regions();
						push(@NONOVERLAP_RESULTS, $results) if $results;
					}
					else {            
						while (threads->list(threads::running)
							>= $CORES_PER_SAMPLE) {usleep(100000)}
						join_stdprob_threads(\@threads, $thread_count);
						$thread_list = threads->new(\&count_nonoverlap_regions);
						push(@threads, $thread_list);
						$thread_count++;
					} 
				}
				if ($curr_chr != $line_array[CHR]) {
					$curr_chr = $line_array[CHR];			  
					# Store centromere info
					($centro_start, $centro_stop) =
						@{$CENTRO_REFS[$curr_chr]}[1,2];
					$gap = 1 if ($line_array[POS] <  $centro_stop);					
				}
				(@CURR_CHR_REFS, @LRR_REFS, @P1_INF_SNPS, @P2_INF_SNPS
				) = (()) x 4;
				(%CH_HD_BY_POS, %P1_HD_BY_POS, %P1_INF_SNP_BY_POS, 
				 %P2_INF_SNP_BY_POS, %P2_HD_BY_POS, %BAF_BY_POS, %LRR_BY_POS
				) = (()) x 7;
	 
			}
			push(@CURR_CHR_REFS, \@line_array); 
		}
		close(INPUT_FILE);
		if (@CURR_CHR_REFS >= $SIZE_REGION) {
			if ($CORES_PER_SAMPLE == 1) {
				my $results = count_nonoverlap_regions();
				push(@NONOVERLAP_RESULTS, $results) if $results;
			}
			else {            
				while (threads->list(threads::running) 
					>= $CORES_PER_SAMPLE) {usleep(100000)}
				join_stdprob_threads(\@threads, $thread_count);
				$thread_list = threads->new(\&count_nonoverlap_regions);
				push(@threads, $thread_list);
				$thread_count++;
			}
		}
		if ($CORES_PER_SAMPLE > 1) {
			# Check if threads are finished
			while (threads->list(threads::running)) {usleep(100000)}
			join_stdprob_threads(\@threads, $thread_count);
		}
		@CURR_CHR_REFS = ();   

		my @sorted_nonoverlap = sort {
			$$a[0] <=> $$b[0] || $$a[1] <=> $$b[1]} @NONOVERLAP_RESULTS;
			
		for (0..$#sorted_nonoverlap) {
			push(@{$sorted_nonoverlap[$_]}, @{$INIT_CH_STATS[$_]}[6,7]);
		}
		
		my @max_sizes = my @normals = ();		
		# Make an array of refs for kmeans
				# map {print OUTPUT_FILE "$$_[5]\t"}@sorted_nonoverlap;
				# print OUTPUT_FILE "\n";
		# map {print OUTPUT_FILE "$$_[6]\t"}@sorted_nonoverlap;
				# print OUTPUT_FILE "\n";

			# map {print OUTPUT_FILE "$$_[10]\t"}@sorted_nonoverlap;
	
		# print OUTPUT_FILE "\n";
					# map {print OUTPUT_FILE $$_[6] / $$_[10], "\t"}@sorted_nonoverlap;
	# print OUTPUT_FILE "\n";

			# map {print OUTPUT_FILE "$$_[3]\t"}@sorted_nonoverlap;
	# print OUTPUT_FILE "\n";
			# map {print OUTPUT_FILE $$_[3] / $$_[6], "\t"}@sorted_nonoverlap;
	# print OUTPUT_FILE "\n";

		map {push(@max_sizes, [$$_[2]])} @sorted_nonoverlap;
			# map {print OUTPUT_FILE "@{$_}\t"}@max_sizes;
	# print OUTPUT_FILE "\n";

		my ($num_clusters, $normal_cluster, $cluster_ref) 
			= cluster(\@max_sizes, 0.5);
		
		if ($num_clusters == 1) {@normals = @sorted_nonoverlap}
		else {
			my @temp = @sorted_nonoverlap;
			until ($num_clusters == 1) {
				@normals = @max_sizes = ();
				
				for (my $i = 0; $i < @temp; $i++) {
					if ($$cluster_ref[$i] == $normal_cluster) {
						push(@normals, $temp[$i]);				
					}
				}
				map {push(@max_sizes, [$$_[2]])} @normals;
				@temp = @normals;
				($num_clusters, $normal_cluster, $cluster_ref) 
					= cluster(\@max_sizes, 0.45);
			}
		}
		# map{print OUTPUT_FILE "@{$$_[4]}\n"} @sorted_nonoverlap;

		my ($baf_snp_ct, $norm_baf_ct) = (0) x 2;
		$baf_snp_ct  += $_ for map {$$_[-2]} @sorted_nonoverlap;
		$norm_baf_ct += $_ for map {$$_[-2]} @normals;
		my @nonoverlap_totals;
		for (0..$SIZE_REGION) {$nonoverlap_totals[$_] = 0}
		for (my $i = 3; $i <= $SIZE_REGION; $i++) {
			my $adj_snp_ct = 0;
			my ($stdev, $mean) = st_dev([map{$$_[4]->[$i]} @normals]);
			for (0..$#normals) {
				my $curr_value = $normals[$_]->[4]->[$i];
				unless ($curr_value > $mean + $stdev * 4) {
					$nonoverlap_totals[$i] += $curr_value;
					$adj_snp_ct += $normals[$_]->[-2];
				}
			}
			my $adjustment_factor = $baf_snp_ct / $adj_snp_ct;
			$nonoverlap_totals[$i] *= $adjustment_factor;
		}
		my $adjusted_mi1_ct = 0;
		$adjusted_mi1_ct += $_ for map {$$_[5]} @normals;
		my $mi1_rate = $adjusted_mi1_ct / $norm_baf_ct;
		#$mi1_rate  += int(estimate_additional_instances($norm_baf_ct, 
		#	$mi1_rate, 1 / $norm_baf_ct));
		$mi1_rate = $adjusted_mi1_ct / $norm_baf_ct;
		$MIN_MI1  = calc_min_adj_snps($norm_baf_ct, $mi1_rate, $ALPHA);
		
		my $total_mi1_ct;
		map {$total_mi1_ct += "$$_[5]\t"}@sorted_nonoverlap;	
		my $total_mi1_rate = $total_mi1_ct / $baf_snp_ct;

		my $het_mi1_ct;
		map {$het_mi1_ct += "$$_[10]\t"} @normals;	
		$HET_MI1_RATE = $het_mi1_ct / $norm_baf_ct;
		#$het_mi1_ct += int(estimate_additional_instances($norm_baf_ct,
		#	$HET_MI1_RATE, 1 / $norm_baf_ct));
		#$HET_MI1_RATE = $het_mi1_ct / $norm_baf_ct;
		#print "het mi1 rate $HET_MI1_RATE $het_mi1_ct\n";
		
		my $baf_outliers;
		$ab_ct = 0;
		map {$baf_outliers += "$$_[11]\t"} @normals;	
		map {$ab_ct += "$$_[12]\t"} @normals;	
		my $baf_outlier_rate = $baf_outliers / $ab_ct;
		$MIN_BAF_OUT = calc_min_adj_snps($ab_ct, $baf_outlier_rate, $ALPHA);
		#print "$baf_outliers $baf_outlier_rate $ab_ct $MIN_BAF_OUT\n";
		# my $min_snps;				
		# my ($n, $k, $p) = (100, 2, $ab_outlier_rate);
		# until ($min_snps) {
			# my $bin_coeff = calculate_bin_coeff($n, $k);
			# if ($bin_coeff * ($p**$k) * ((1 - $p)**($n - $k)) <= $ALPHA) {
				# $min_snps = $k;
			# }
			# $k++;
		# }
		# $MIN_BAF_OUT_RATE = $min_snps / 100;
		
		
		#print "baf outlier rate $MIN_BAF_OUT $ab_outlier_rate $ab_outliers $ab_ct\n";
		
		# Calculate hd rates and minimum size of hd regions
		my ($p1_hd_ct, $p2_hd_ct, $ch_hd_ct, $norm_lrr_snp_ct) = (0) x 4;				
		$p1_hd_ct        += $_ for map {$$_[7]}  @normals;
		$p2_hd_ct        += $_ for map {$$_[8]}  @normals;
		$ch_hd_ct        += $_ for map {$$_[9]}  @normals;		
		$norm_lrr_snp_ct += $_ for map {$$_[-1]} @normals;
		my $p1_hd_rate  = $p1_hd_ct / $norm_lrr_snp_ct;
		my $p2_hd_rate  = $p2_hd_ct / $norm_lrr_snp_ct;
		my $ch_hd_rate  = $ch_hd_ct / $norm_lrr_snp_ct;
		#$p1_hd_ct      += int(estimate_additional_instances($norm_lrr_snp_ct,
		#					$p1_hd_rate, 1 / $norm_lrr_snp_ct));
		$p1_hd_rate     = $p1_hd_ct / $norm_lrr_snp_ct;
		#$p2_hd_ct      += int(estimate_additional_instances($norm_lrr_snp_ct,
		#					$p2_hd_rate, 1 / $norm_lrr_snp_ct));
		$p2_hd_rate     = $p2_hd_ct / $norm_lrr_snp_ct;
		#$ch_hd_ct      += int(estimate_additional_instances($norm_lrr_snp_ct,
		#					$ch_hd_rate, 1 / $norm_lrr_snp_ct));
		$ch_hd_rate     = $ch_hd_ct / $norm_lrr_snp_ct;
		$MIN_P1_HD      = calc_min_adj_snps($norm_lrr_snp_ct, $p1_hd_rate, $ALPHA);
		$MIN_P2_HD      = calc_min_adj_snps($norm_lrr_snp_ct, $p2_hd_rate, $ALPHA);
		$MIN_CH_HD      = calc_min_adj_snps($norm_lrr_snp_ct, $ch_hd_rate, $ALPHA);
		
		# Default assumption of error rate is 3.75 times the rate of single
		# Mendelian errors as determined by the ratio of all single 
		# errors derived from mutating all genotype combinations 
		# of normal inheritance to all detectable errors (24 detectable 
		# changes out of 90 possible changes). The ratio (90/24) is then used
		# to estimate the overall "error" rate (single Mendelian errors and 
		# single SNP biological abnormalities). In a previous step, single 
		# Mendelian errors were counted while discarding any which were 
		# adjacent to reduce influence of large biological abnormalities in
		# error rate estimations.
		$ERROR_RATE = $mi1_rate * 3.75;
		my $rounded_error_rate = sprintf("%.5f", $ERROR_RATE); 
		
		# Calculate minimum SNP size for region detected by standard POD analysis.
		# When the number of nonoverlapping windows >= the incrementing snp size 
		# divided into the 1 / probability of the snp size is >= the 
		# default number of genomes, the minimum SNP size is stored.
		# For a window containing 9 inf snps of which are from the same parent,
		# the probability is calculated as the binomial coefficent of 9 out of 9,
		# plus the bin coeff of 8/9 * the probability that one is an error to the 
		# power of 1 * the number of possible ways to achieve 8/9 and one error 
		# by chance (9), + bin coeff of 7/9 * error rate ^2 * 9 * 8, etc.
		my $ct = 3;
		my $totals;
		until ($MIN_POD || $ct == $SIZE_REGION) {
			my $n = my $k = $ct;
			my $p = 0.5;
			my $bin_coeff = calculate_bin_coeff($n, $k);
			my $pvalue = $bin_coeff * ($p**$k) * ((1 - $p)**($n - $k));

			#my $pvalue = 0; 
			#for (0..$k) {
				# my $bin_coeff = calculate_bin_coeff($n, $k);
				# $pvalue += $bin_coeff * ($p**$k) * ((1 - $p)**($n - $k));
				# print OUTPUT_FILE "$n $_ $pvalue\n";
			#}
			#print "$ct $nonoverlap_totals[$ct] $prob\n";
			if ($nonoverlap_totals[$ct]) {$totals = $nonoverlap_totals[$ct]}
			else {$totals = 1}
			my $sidak_alpha = 1 - (1 - $ALPHA)**(1 / $totals);
			print "$n pvalue $pvalue sidak $sidak_alpha\n";

			if ($pvalue <= $sidak_alpha) {
				$MIN_POD = $ct;
				my ($n, $k, $p) = ($MIN_POD, 0, 0.5);
				my $bin_coeff = calculate_bin_coeff($n, $k);
				my $prob = 2 * ($bin_coeff * ($p**$k) * ((1 - $p)**($n - $k)));				
				$POD_ALPHA = $prob;
			}
			$ct++;
		}
			#print $MIN_POD, "\n";
		my $next = check_quality($total_mi1_rate, $p1_nc_rate, $p2_nc_rate, 
			$ch_nc_rate); 
		next FILE if $next;

		# Calculate Error Threshold as percent of a window, given the significance 
		# threshold and estimated error rate.
		# Binomial Probability formula
		# P(k) = (n choose k) * p^k * (1-p)^(n-k), where
		# n = number of SNPs per window
		# k = number of errors per window
		# n - k = number of correct calls per window
		# p = probability of an error in one trial
		
		# Iterate until the maximum of "errors" is reached before P(k)
		# < the significance threshold. This is the number of "acceptable
		# errors" (informative SNPs contributed by the opposite parent in
		# a true abnormal region) per bin that may happen by chance, without 
		# indicating that a region has ended.
		
		my ($stop, $errors) = (0) x 2;
		until ($stop) {
			my ($n, $k, $p) = ($SIZE_REGION, $errors, $ERROR_RATE);
			my $prob = 0;
			for ($k..$n) {
				my $bin_coeff = calculate_bin_coeff($n, $_);
				$prob += $bin_coeff * ($p**$_) * ((1 - $p)**($n - $_));
			}
			if ($prob <= 0.001) {
				$stop = 1;
				$ACCEPTABLE_ERRORS = $errors - 1;
				$ACCEPTABLE_ERRORS ||= 0;
			}
			$errors++;
		}
		
		###############################Print Parameters##############################
			
		$GENDER ||= "NA";

		if ($verbose) {
			#Print Parameters to screen
			print "PARAMETERS:\n\n",
				"Window size (SNPs)               = $SIZE_REGION\n",
				"Significance threshold           = $ALPHA\n";
			if ($POD) {print 
				"Min POD region (SNPs)            = $MIN_POD\n"}
			if ($MI1) {print 
				"Min MI1 region (SNPs)            = $MIN_MI1\n"}
			if ($HD) {print 
				"Father Min HD region (SNPs)      = $MIN_P1_HD\n",
				"Mother Min HD region (SNPs)      = $MIN_P2_HD\n",
				"Child Min HD region  (SNPs)      = $MIN_CH_HD\n"}
			print 
				"AA boundary                      = $rounded_AA_bound\n",
				"BB boundary                      = $rounded_BB_bound\n",
				"AB upper boundary                = $rounded_AB_up_bound\n",
				"AB lower boundary                = $rounded_AB_low_bound\n",
				"Estimated Error Rate             = $rounded_error_rate\n",
				"Gender                           = $GENDER\n\n";
		}
		
		$GENDER = 0 if ($GENDER eq "NA");

		if ($GRAPHICS) { 
			# Open file from R containing the current R PID
			get_R_PID();
		
			# Check if the R script is still running
			$R_ATTEMPTS = determine_R_status($R_ATTEMPTS);
		}
		
		################################## Analysis #################################
		($thread_count, $thread_list, $curr_chr) = (0) x 3;
		@threads = ();
		($centro_start, $centro_stop, $gap) = (0) x 3;
		my (@detected_regions, @results);
		
		if(!open(INPUT_FILE, "<", $INPUT_FILE)) { 
			print "Could not open input file: $INPUT_FILE!\n";
			print STDERR "Could not open input file: $INPUT_FILE!\n";		
			next FILE;
		}
		
		while (<INPUT_FILE>) {
			my $line = $_;
			chomp($line);
			my @line_array = split(/\t/, $line);
			next if (!$line_array[0]);
			#Setup
			if ($. > 1) {
				# Convert chromosome identifier to numeric 
				if ($line_array[CHR] !~ /^\d/) {
					if ($line_array[CHR] eq "X" 
						&& $GENDER) {$line_array[CHR] = CHR_X}
					else {next}
				}
				$curr_chr = $line_array[CHR];
				print "Analyzing Chromosome ", $curr_chr, "\n" if $verbose;
				
				# Store centromere info
				($centro_start, $centro_stop) = @{$CENTRO_REFS[$curr_chr]}[1,2];
				$gap = 1 if ($line_array[POS] <  $centro_stop);
				push(@CURR_CHR_REFS, \@line_array); 
				last;
			}
		}        

		INPUT: while (<INPUT_FILE>) {			
			my $line = $_;
			chomp($line);
			my @line_array = split(/\t/, $line);
			next INPUT if (!$line_array[0]);
			
			# Convert chromosome identifier to numeric 
			if ($line_array[CHR] !~ /^\d/) {
				if ($line_array[CHR] eq "X" 
					&& $GENDER) {$line_array[CHR] = CHR_X}
				else {next INPUT}
			}

			# Start new thread when centromere or new chromosome is encountered                
			if ($curr_chr != $line_array[CHR]
				|| ($gap && $line_array[POS] > $centro_start)) { 
				$gap = 0;
				
				if (@CURR_CHR_REFS >= $SIZE_REGION) {
					# Create threads only if >1 processor is available                
					if ($CORES_PER_SAMPLE == 1) {
						my $results = manage_chromosome_arm();
						if ($results) {
							push(@detected_regions, @{$$results[0]});
							$INF_SNPS       += $$results[1];
							$ADJUSTED_MI1   += $$results[2];
							if ($$results[3]) {
								if ($LOCAL_STATS[$$results[3]]) {
									$LOCAL_STATS[$$results[3] + 30] = $$results[4];
								}                
								else {$LOCAL_STATS[$$results[3]] = $$results[4]}
							}
						}
					}
					else {            
						while (threads->list(threads::running)
							>= $CORES_PER_SAMPLE) {usleep(100000)}
						join_analysis_threads(\@results,
							\@threads, $thread_count);     
						$thread_list = threads->new(\&manage_chromosome_arm);
						push(@threads, $thread_list);
						$thread_count++;
					} 
				}
				if ($curr_chr != $line_array[CHR]) {
					$curr_chr = $line_array[CHR];
					print "Analyzing Chromosome ", $curr_chr, "\n" if $verbose;
			  
					# Store centromere info
					($centro_start, $centro_stop) 
						= @{$CENTRO_REFS[$curr_chr]}[1,2];
					$gap = 1 if ($line_array[POS] <  $centro_stop);
				}
				(@CURR_CHR_REFS, @LRR_REFS, @P1_INF_SNPS, @P2_INF_SNPS
				) = (()) x 4;
				(%CH_HD_BY_POS, %P1_HD_BY_POS, %P1_INF_SNP_BY_POS, 
				 %P2_INF_SNP_BY_POS, %P2_HD_BY_POS, %BAF_BY_POS, %LRR_BY_POS
				) = (()) x 7;
	 
			}
			push(@CURR_CHR_REFS, \@line_array); 
		}
		close(INPUT_FILE);
		if (@CURR_CHR_REFS >= $SIZE_REGION) {
			if ($CORES_PER_SAMPLE == 1) {
				my $results = manage_chromosome_arm();
				if ($results) {
					push(@detected_regions, @{$$results[0]});
					$INF_SNPS       += $$results[1];
					$ADJUSTED_MI1   += $$results[2];
					if ($$results[3]) {
						if ($LOCAL_STATS[$$results[3]]) {
							$LOCAL_STATS[$$results[3] + 30] = $$results[4];
						}                
						else {$LOCAL_STATS[$$results[3]] = $$results[4]}
					}
				}
			}
			else {            
				while (threads->list(threads::running)
					>= $CORES_PER_SAMPLE) {usleep(100000)}
				$thread_list = threads->new(\&manage_chromosome_arm);
				push(@threads, $thread_list);
				$thread_count++;
				join_analysis_threads(\@results, \@threads, $thread_count);
			}
		}
		if ($CORES_PER_SAMPLE > 1) {
			# Check if threads are finished
			while (threads->list(threads::running)) {usleep(100000)}
			join_analysis_threads(\@results, \@threads, $thread_count);
			map {push(@detected_regions, @{$_}) if $_} @results;
		}
		@CURR_CHR_REFS = ();   
		
		# Refine autosomal medians to reflect only normal regions in child
		
		my @temp;
		map {push(@temp, $LOCAL_STATS[$_]) if $LOCAL_STATS[$_]} (0..52);	
		$P1_mBAF_MED = median([map {$$_[0]} @temp]);
		$P2_mBAF_MED = median([map {$$_[1]} @temp]);
		$CH_mBAF_MED = median([map {$$_[2]} @temp]);
		$P1_LRR_MED  = median([map {$$_[3]} @temp]);
		$P2_LRR_MED  = median([map {$$_[4]} @temp]);
		$CH_LRR_MED  = median([map {$$_[5]} @temp]);
		@temp = ();
		#print "refined medians $P1_mBAF_MED $P2_mBAF_MED $CH_mBAF_MED $P1_LRR_MED $P2_LRR_MED $CH_LRR_MED\n";
		####################### Descriptive Calculations #######################            
		#my $informative_snp_rate = $INF_SNPS / $INIT_SUMS[17];
		# Check if regions were detected
		if ($#detected_regions < 0) {
			print OUTPUT_FILE "\n\nNo regions of abnormal parental contribution ",
				"were detected.\n\n" unless $BATCH;
			print STDERR "Sample $FILENAME : No regions of abnormal parental contribution ",
				"were detected.\n" if $BATCH;	
			print "\n\nNo regions of abnormal parental contribution were detected.",
				"\n\n" if $verbose;
			$completed[1]++;
			next FILE;
		}
	   
		# Sort output array by chromosome and position     
		my @sorted_detected_regions = sort {
			$$a[0] <=> $$b[0] || $$a[1] <=> $$b[1]} @detected_regions;
					
		my ($num_pod, $num_mi1, $num_hd, $num_podcr) = (0) x 4;

		######################## Prepare results for output file ####################
				
		my $count = 0;	
		foreach my $ref (@sorted_detected_regions) {
			
			if ($$ref[17] == POD)      {$num_pod++}
			elsif ($$ref[17] == MI1)   {$num_mi1++}
			elsif ($$ref[17] == HD)    {$num_hd++}
			elsif ($$ref[17] == PODcr) {$num_podcr++}			
			
			# If a region comprised most of a chromosome arm and underwent
			# a single round of analyses, adjust the regions medians by
			# the medians of the normal SNPs from the other arm, or if 
			# most of the chromosome is abnormal, adjust using the autosomal
			# medians.
			my ($p1_mBAF, $p2_mBAF, $ch_mBAF, $p1_LRR, $p2_LRR, $ch_LRR);
			if ($$ref[19] == 1) {
				if ($LOCAL_STATS[$$ref[0]]) {
					# local medians 
					($p1_mBAF, $p2_mBAF, $ch_mBAF, $p1_LRR, $p2_LRR, $ch_LRR) = 
						@{$LOCAL_STATS[$$ref[0]]};
				}
				elsif ($LOCAL_STATS[$$ref[0] + 30]) {
					# local medians
					($p1_mBAF, $p2_mBAF, $ch_mBAF, $p1_LRR, $p2_LRR, $ch_LRR) = 
						@{$LOCAL_STATS[$$ref[0] + 30]};
				}
				else {
					# refined autosomal medians
					($p1_mBAF, $p2_mBAF, $ch_mBAF, $p1_LRR, $p2_LRR, $ch_LRR) = 
						($P1_mBAF_MED, $P2_mBAF_MED, $CH_mBAF_MED,
						$P1_LRR_MED, $P2_LRR_MED, $CH_LRR_MED);
				}			
				$$ref[10] = sprintf("%.4g", $$ref[10] - ($ch_mBAF - 0.5))
					unless $$ref[10] eq "NA";
				$$ref[11] = sprintf("%.4g", $$ref[11] - $ch_LRR)
					unless $$ref[11] eq "NA";
				$$ref[12] = sprintf("%.4g", $$ref[12] - ($p1_mBAF - 0.5))
					unless $$ref[12] eq "NA";
				$$ref[13] = sprintf("%.4g", $$ref[13] - $p1_LRR)
					unless $$ref[13] eq "NA";
				$$ref[14] = sprintf("%.4g", $$ref[14] - ($p2_mBAF - 0.5))
					unless $$ref[14] eq "NA";
				$$ref[15] = sprintf("%.4g", $$ref[15] - $p2_LRR)
					unless $$ref[15] eq "NA";
				if ($$ref[11] > $AMP + $ch_LRR) {$$ref[3] = "AMP"}
				elsif ($$ref[11] < $DEL_UPPER + $ch_LRR && 
					$$ref[11] > $HD_HI + $ch_LRR) {$$ref[3] = "DEL"}
				elsif ($$ref[11] <= $HD_HI + $ch_LRR) {$$ref[3] = "HD"}
				else {$$ref[3] = " "}
				$$ref[20] = call_inheritance($ref);
			}
			if    ($$ref[4] == FATHER)  {$$ref[4] = "Father"}
			elsif ($$ref[4] == MOTHER)  {$$ref[4] = "Mother"}
			elsif ($$ref[4] == NONE)    {$$ref[4] = "None"}
			elsif ($$ref[4] == BOTH)    {$$ref[4] = "Both"}        
			elsif ($$ref[4] == UNKNOWN) {$$ref[4] = "Unknown"}
			if (!$$ref[3] && ($$ref[4] eq "Father" || $$ref[4] eq "Mother")) {
				$$ref[4] = $$ref[4] . " (C)";
			}
			
			my $detection;         
			if ($$ref[17] == HD)       {$detection = "HD"}
			elsif ($$ref[17] == MI1)   {$detection = "MI1"}
			elsif ($$ref[17] == POD)   {$detection = "POD"}
			elsif ($$ref[17] == PODcr) {$detection = "PODcr"}			
			
			push(@output_array, [$CH_NAME, @$ref[0..4,20,5..7], $detection, 
				@$ref[10..15]]);					
		}
		push(@trio_def, "Trio Definition Autosomal NoCall Rate",
			"Father =   $P1_NAME\t$p1_nc_rate",
			"Mother =   $P2_NAME\t$p2_nc_rate",
			"Child  =   $CH_NAME\t$ch_nc_rate") unless ($BATCH);		
			
		############## Create file to pass results to R script ##################
		if ($GRAPHICS) { 
			print $PERL_TO_R_FILE "Chr\tStart\tStop\tType\tParent\tSNPs",
				"\tInformative SNPs\tSize of Region (bp)\tRadius of region for R",
				"\tbp of midpoint for R\tMedian BAF\tMedian LRR\tFather-Median BAF",
				"\tFather-Median LRR\tMother-Median BAF\tMother-Median LRR\n";

			#for (@sorted_detected_regions) {
			#	print $PERL_TO_R_FILE join("\t", @$_
			map {print $PERL_TO_R_FILE join("\t", @{$_}[0..15]), "\n"} @sorted_detected_regions;
#print get_R_PID(), "\n";

			print "Creating Graphics\n" if $verbose;
			
			#$R_ATTEMPTS = determine_R_status($R_ATTEMPTS);
			#next FILE if $R_ATTEMPTS > 1;
			
			# Wait for R script to finish and return results
			my $waiting_for_R = 0;
			until ($waiting_for_R) {
				if (-s $R_TO_PERL_FILENAME) {$waiting_for_R = 1}
				else {
					$R_ATTEMPTS = determine_R_status($R_ATTEMPTS);
					usleep(500000);
				}
			}
		}                        
		my $TIME = time() - $^T;
		print "finished in $TIME seconds\n" if $verbose;
		if ($STATS) {
			push(@stats_array, [$FILENAME, $CH_NAME, $P1_NAME, 
			$P2_NAME, $ch_nc_rate, $p1_nc_rate, $p2_nc_rate, $ch_hd_rate, 
			$MIN_CH_HD, $p1_hd_rate, $MIN_P1_HD, $p2_hd_rate, $MIN_P2_HD,
			$mi1_rate, $MIN_MI1, $MI1_UPPER_THRESH, $MI1_LOWER_THRESH, 
			$MIN_POD, $POD_ALPHA, $ALPHA, $ERROR_RATE, 
			$ACCEPTABLE_ERRORS, $AA_BOUND, $BB_BOUND, $AB_UPPER_BOUND, 
			$AB_LOWER_BOUND, scalar(@sorted_detected_regions), $num_pod, 
			$num_podcr, $num_mi1, $num_hd]);
		}
		$completed[0]++;	
			print "completed $file\n";

	}
	return([\@output_array, \@stats_array, \@trio_def, \@completed]);
}












################################# Sub Functions ##############################
sub cluster {
	# Returns optimal number of clusters (1 or 2) and array of cluster results
	my ($data_ref, $y) = @_;	
	my %param = (
		nclusters => 2,
		data => $data_ref,
		mask => '',
		weight => '',
		transpose => 0,
		npass => 10,
		method => 'a',
		dist => 'e',
		initialid => [],
	);
	my @clusters = Algorithm::Cluster::kcluster(%param);
	my (@cluster1, @cluster2);

	# The jump method is used to determine the optimal number of clusters
	# between 1 and 2.
	for (my $i = 0; $i < @{$clusters[0]}; $i++) {
		if ($clusters[0]->[$i] == 0) {push(@cluster1, $$data_ref[$i]->[0])} 
		if ($clusters[0]->[$i] == 1) {push(@cluster2, $$data_ref[$i]->[0])} 
	}
	
	#print OUTPUT_FILE "cluster1 @cluster1\n";
	#print OUTPUT_FILE "cluster2 @cluster2\n";

	my @clustered = ([@cluster1, @cluster2], \@cluster1, \@cluster2);
	my @sum_dists;
	
	foreach my $cluster (@clustered) {
		my ($sum, $sum_dist) = (0) x 2;
		$sum += $_ for @$cluster;
		my $mean = $sum / @$cluster;
		$sum_dist += ($_ - $mean)**2 for @$cluster;
		push(@sum_dists, $sum_dist);
	}
	
	my @distortions;
	push(@distortions, $sum_dists[0] / @$data_ref); 
	push(@distortions, ($sum_dists[1] + $sum_dists[2]) / @$data_ref);
	
	# Compute transformed distortions
	my @trans_dist;
	for (@distortions) {push(@trans_dist, $_**-$y)} 
	my $optimal_num_clusters = 0;
	# Compute max jump
	if ($trans_dist[0] >= $trans_dist[1] - $trans_dist[0] ) {
		$optimal_num_clusters = 1;
	}
	else {$optimal_num_clusters = 2}
	#print OUTPUT_FILE "$trans_dist[0] ", $trans_dist[1] - $trans_dist[0], "\n"; 
	my $normal_cluster = 0;
	if ($optimal_num_clusters == 2) {
		my ($a_sum, $a_mean, $b_sum, $b_mean) = (0) x 4;
		map {$a_sum += $_} @cluster1;
		$a_mean = $a_sum / @cluster1;
		map {$b_sum += $_} @cluster2;
		$b_mean = $b_sum / @cluster2;
		
		if ($a_mean > $b_mean) {$normal_cluster = 1}
	}
	#print OUTPUT_FILE "optimal num clusters $optimal_num_clusters\n";
	return($optimal_num_clusters, $normal_cluster, $clusters[0]);
}
		
sub calculate_bin_coeff {
    # Estimates and returns the binomial coefficient for n and k.  
    my ($n, $k) = @_;
    return 1 if $k == 0;
    my $r = $n - $k + 1;
    for (2..$k) {$r *= (($n - $k)/$_ + 1)}
    return ($r);
}

sub log10 {
	my $x = $_[0];
	return if ($x == 0);
	return(log($x)/log(10));
}

sub approx_power_pdf {
	my ($n, $k, $p) = @_;
	return 1 if $k == 0;
	my $numerator = ($k * log10($n * $p)) + log10(exp(1)**(-$n * $p));
	my $denominator = (($k + 0.5) * log($k) - $k + log(2 * PI)/2) 
		* log10(exp(1));
	my $log_approx = $numerator - $denominator;
	#print OUTPUT_FILE "in approx $numerator $denominator $log_approx\n";
	return($log_approx);
}

# sub pdf {
	# my ($n, $k, $p) = @_;
	# return 1 if $k == 0;
	# my $k_fact = 1;
	# $k_fact *= $_ for 1..$k;
	# my $numerator = 2.7182818284590452353602874713527**(-$n*$p) * $n*$p**$k;
	# return(0) if ($numerator == 0);
	# my $prob = (2.7182818284590452353602874713527**(-$n*$p) * ($n*$p)**$k) / $k_fact;
	# return($prob);
# }
	
sub calc_init_stats {
	my $init_stats = $_[0];
	# Identify and ignore abnormal chromosome arms
	my @stdevs = my @normals = ();
	map {push(@stdevs, [$$_[9]])} @$init_stats;
	#map {print OUTPUT_FILE "@{$_}\t"}@stdevs;
	#print OUTPUT_FILE "\n";
	my ($num_clusters, $normal_cluster, $cluster_ref) 
		= cluster(\@stdevs, 0.5);
		
	if ($num_clusters == 1) {@normals = @$init_stats}
	else {
		my @temp = @$init_stats;
		until ($num_clusters == 1) {
			@normals = @stdevs = ();
			for (my $i = 0; $i < @temp; $i++) {
				if ($$cluster_ref[$i] == $normal_cluster) {
					push(@normals, $temp[$i]);				
				}
			}
			map {push(@stdevs, [$$_[9]])} @normals;
			@temp = @normals;
			($num_clusters, $normal_cluster, $cluster_ref) 
				= cluster(\@stdevs, 0.45);
		}
	}
	# Calculate median of mBAF medians of "normal" chromosome arms
	my $mBAF_med    = median([map{$$_[3]} @normals]);
	my $lrr_med     = median([map{$$_[4]} @normals]);
	# Calculate NC rate of "normal" chromosome arms				
	my $nc_ct    = my $nc_snp_ct = 0;
	$nc_ct      += $_ for map {$$_[5]} @normals;
	$nc_snp_ct  += $_ for map {$$_[8]} @normals;					
	my $nc_rate  = sprintf("%.4f", $nc_ct / $nc_snp_ct); 
	
	# Calculate BAF thresholds
	my ($aa_sums, $aa_sum_sq, $aa_ct) = (0) x 3;
	foreach (@normals) {
		$aa_sums   += $$_[2]->[0];
		$aa_sum_sq += $$_[2]->[1];
		$aa_ct     += $$_[2]->[2];
	}	
	my @AA_stats = st_dev($aa_sums, $aa_sum_sq, $aa_ct);
	my ($ab_sums, $ab_sum_sq, $ab_ct) = (0) x 3;
	foreach (@normals) {
		$ab_sums   += $$_[2]->[3];
		$ab_sum_sq += $$_[2]->[4];
		$ab_ct     += $$_[2]->[5];
	}	
	my @AB_stats = st_dev($ab_sums, $ab_sum_sq, $ab_ct);
	my ($bb_sums, $bb_sum_sq, $bb_ct) = (0) x 3;
	foreach (@normals) {
		$bb_sums   += $$_[2]->[6];
		$bb_sum_sq += $$_[2]->[7];
		$bb_ct     += $$_[2]->[8];
	}	
	my @BB_stats = st_dev($bb_sums, $bb_sum_sq, $bb_ct);
	my $aa_bound = $AA_stats[0] + ($AA_stats[1] * $HOM_SD);
	my $bb_bound = $BB_stats[0] - ($BB_stats[1] * $HOM_SD);
	my $ab_upper = $AB_stats[0] + ($AB_stats[1] * $HET_SD);
	my $ab_lower = $AB_stats[0] - ($AB_stats[1] * $HET_SD);
	
	return(\@normals, $mBAF_med, $lrr_med, $nc_rate,
			$aa_bound, $bb_bound, $ab_upper, $ab_lower);
}

sub calc_min_adj_snps {
    # Calculate minimum number of SNPs, given a significance 
    # threshold and estimated rate.
    # Iterate until the number SNPs 
    # < the significance threshold. 
    my ($n, $p, $thresh) = @_;
    my $min_SNPs;
    my $k = 2;
    until ($min_SNPs) {
		my $prob = streak_prob($n, $k, $p);
		@STREAK_ARRAY = ();
        if ($prob <= $thresh) {$min_SNPs = $k}
        $k++;
    }
    return($min_SNPs);
}
  
sub overlapping_windows {
	# Counts and stores parental contribution for each window
	my ($array_ref, $window_size, $small) = @_;
	my ($p1_ct, $p2_ct, $index_count) = (0) x 3;
	my $chr = $$array_ref[0]->[CHR];
	#print OUTPUT_FILE "$chr $$array_ref[0]->[POS] $$array_ref[-1]->[POS]\n";
	my $array_size = scalar(@$array_ref);
	my @return_refs;
    for (my $i = 0; $i < $window_size; $i++) {   
        my $curr_ref = $$array_ref[$i];
        if ($$curr_ref[-1] == FATHER) {$p1_ct++} 
        elsif ($$curr_ref[-1] == MOTHER) {$p2_ct++}
    }
	
    while ($index_count + $window_size < $array_size) {	 
        my $start_ref = $$array_ref[$index_count];
        my $end_ref = $$array_ref[$index_count + $window_size];
		
        if ($$end_ref[-1] == FATHER) {$p1_ct++} 
        elsif ($$end_ref[-1] == MOTHER) {$p2_ct++} 
        
		if ($small) {
			#unless ($p1_ct < $ACCEPTABLE_ERRORS && $p2_ct < $ACCEPTABLE_ERRORS
			#	&& $p1_ct + $p2_ct < $MIN_POD) { 
			#print OUTPUT_FILE join("\t", $chr, $$start_ref[POS], $$end_ref[POS], 
			#		$p1_ct, $p2_ct, $outlier_ct), "\n";
				push(@return_refs, [$chr, $$start_ref[POS], $$end_ref[POS], 
					$p1_ct, $p2_ct]);
			#}
		}
		else {
			# my ($n, $k, $p) = ($p1_ct + $p2_ct, $p1_ct, 0.5);
            # my $bin_coeff = calculate_bin_coeff($n, $k);
            # my $prob = ($bin_coeff * $p**$n);

			my $max = $p1_ct >= $p2_ct ? $p1_ct : $p2_ct;
			my $min = $p1_ct <= $p2_ct ? $p1_ct : $p2_ct;
			#print OUTPUT_FILE join("\t", $chr, $$start_ref[POS], $$end_ref[POS], 
			#		$p1_ct, $p2_ct, $outlier_ct), "\n";
			if ($max >= 15 && $min <= 1) { 
			#if ($prob <= $POD_ALPHA) {
				push(@return_refs, [$chr, $$start_ref[POS], $$end_ref[POS], 
					$p1_ct, $p2_ct]);
			}
		}
         
        # Adjust counts     
        if ($$start_ref[-1] == FATHER) {$p1_ct--} 
        elsif ($$start_ref[-1] == MOTHER) {$p2_ct--} 		
        $index_count++;  
    }
	return(\@return_refs)
}

  
sub calculate_region_median {
    my ($region_ref, $input, $index) = @_;    
    my ($array_ref, $hash_ref, @median_array);
    my $count = 0;
    my ($median, $lowest_index, $highest_index, $lowest_value, $highest_value);
	my ($aa_bound, $bb_bound) = ($AA_BOUND, $BB_BOUND);
    if ($input == LRR) {
        $array_ref = \@LRR_REFS;
        $hash_ref  = \%LRR_BY_POS;
    }
    elsif ($input == BAF) {
        $array_ref = \@CURR_CHR_REFS;
        $hash_ref  = \%BAF_BY_POS;
		if ($index == 4) {
			($aa_bound, $bb_bound) = ($P1_AA_BOUND, $P1_BB_BOUND);
		}
		elsif ($index == 7) {
			($aa_bound, $bb_bound) = ($P2_AA_BOUND, $P2_BB_BOUND);
		}
    }
	
	#print STDERR "calcregmed to find\n";
		#				print STDERR "@$region_ref[START,STOP]\n";

    my ($start_index, $stop_index, $present) 
        = find_closest_indices(@$region_ref[START,STOP], $hash_ref);
    if ($present) {
        $lowest_value = $highest_value = 0;
        #$lowest_index = $highest_index = $start_index + 1;
        for ($start_index..$stop_index) {
            if ($input == LRR || $$array_ref[$_]->[$index] >= $aa_bound &&
                $$array_ref[$_]->[$index] <= $bb_bound) {
                if ($$array_ref[$_]->[$index] < $lowest_value 
					&& $_ != $start_index && $_ != $stop_index) {
                    $lowest_value = $$array_ref[$_]->[$index];
                    $lowest_index = $_;
                }
                if ($$array_ref[$_]->[$index] > $highest_value
					&& $_ != $start_index && $_ != $stop_index) {
                    $highest_value = $$array_ref[$_]->[$index];                
                    $highest_index = $_;
                }
                if ($input == BAF && $$array_ref[$_]->[$index] < 0.5) {
                    push(@median_array, (1 - $$array_ref[$_]->[$index]));
                }
                else {push(@median_array, $$array_ref[$_]->[$index])}
            }
        }
        if (@median_array) {$median = median(\@median_array)}
        elsif ($input == BAF) {$median = 1}        
    }
    else {$median = NA}
    return([$median, $start_index, $stop_index, $lowest_index, $highest_index,
    $lowest_value, $highest_value]);
}
   
sub calculate_region_stats {
    # Chr, Start, Stop, Type, Parent, #SNPs, #Informative SNPs, Size(Mb), Radius, Midpoint 
    my $region_ref = $_[0];
    my @stats = @$region_ref;
    # Calculate number of SNPs in region    
    $stats[5] = count_SNPs(@stats[START,STOP], \%SNP_BY_POS);
    # Calculate size(bp) of region
    $stats[7] = $stats[STOP] - $stats[START];
    # Calculate radius for graphics
    $stats[8] = $stats[7] / 2;
    #Calculate midpoint for graphics
    $stats[9] = $stats[STOP] - $stats[8];
             
    # Calculate parental and progeny BAF and LRR medians for a 
    # detected region.            
    my $p1_lrr_median  = calculate_region_median($region_ref, LRR, 3)->[0];
    my $p2_lrr_median  = calculate_region_median($region_ref, LRR, 4)->[0];
    my $ch_lrr_median  = calculate_region_median($region_ref, LRR, 5)->[0];
    my $p1_mbaf_median = calculate_region_median($region_ref, BAF, 4)->[0];
    my $p2_mbaf_median = calculate_region_median($region_ref, BAF, 7)->[0];
    my $ch_mbaf_median = calculate_region_median($region_ref, BAF, 10)->[0];
    
    if($REFINING) {
        if ($ch_mbaf_median == 1) {$stats[10] = sprintf("%.4g", $ch_mbaf_median)}
        else {$stats[10] = sprintf("%.4g", $ch_mbaf_median - ($LOCAL_CH_mBAF_MED - 0.5))}
        $stats[11] = sprintf("%.4g", $ch_lrr_median - $LOCAL_CH_LRR_MED);
        if ($p1_mbaf_median == 1) {$stats[12] = sprintf("%.4g", $p1_mbaf_median)}
        else {$stats[12] = sprintf("%.4g", $p1_mbaf_median - ($LOCAL_P1_mBAF_MED - 0.5))}
        $stats[13] = sprintf("%.4g", $p1_lrr_median - $LOCAL_P1_LRR_MED);
        if ($p2_mbaf_median == 1) {$stats[14] = sprintf("%.4g", $p2_mbaf_median)}
        else {$stats[14] = sprintf("%.4g", $p2_mbaf_median - ($LOCAL_P2_mBAF_MED - 0.5))}
        $stats[15] = sprintf("%.4g", $p2_lrr_median - $LOCAL_P2_LRR_MED);
    }
    else {
        $stats[10] = sprintf("%.4g", $ch_mbaf_median);
        $stats[11] = sprintf("%.4g", $ch_lrr_median);
        $stats[12] = sprintf("%.4g", $p1_mbaf_median);
        $stats[13] = sprintf("%.4g", $p1_lrr_median);
        $stats[14] = sprintf("%.4g", $p2_mbaf_median);
        $stats[15] = sprintf("%.4g", $p2_lrr_median);
    }
    
    if ($ch_lrr_median < $HD_MID || $ch_mbaf_median == NA) {$stats[10] = "NA"}
    if ($p1_lrr_median < $HD_MID || $p1_mbaf_median == NA) {$stats[12] = "NA"}
    if ($p2_lrr_median < $HD_MID || $p2_mbaf_median == NA) {$stats[14] = "NA"}
    if ($ch_lrr_median == NA) {$stats[11] = "NA"}
    if ($p1_lrr_median == NA) {$stats[13] = "NA"}
    if ($p2_lrr_median == NA) {$stats[15] = "NA"}
	
	if ($stats[11] > $AMP) {$stats[3] = "AMP"}
	elsif ($stats[11] < $DEL_UPPER && $stats[11] > $HD_HI) {$stats[3] = "DEL"}
	elsif ($stats[11] <= $HD_HI) {$stats[3] = "HD"}
	else {$stats[3] = " "}
	
    @$region_ref = @stats;
}

sub call_inheritance {
    # Determines whether a likely inheritance pattern can be called.
    # Calls are made for two categories : inherited and 
    # inherited with novel copy number state. The absence of a called
    # inheritance pattern does not imply a de novo abnormality, but 
    # includes de novo, UPD, low-level mosaicism, and regions which
    # the artificial thresholds have failed to determine the true pattern.    
    my $region = $_[0];
	#print OUTPUT_FILE "$CH_NAME inh @$region\n";
    #my $ch_AMP = $$region[11] > $AMP;
	#my $ch_HD  = $$region[11] <= $HD_HI && $$region[11] != NA;
    #my $ch_DEL = $$region[11] < $DEL_UPPER && $$region[11] > $HD_HI;
	#if    ($ch_AMP) {$$region[3] = "AMP"}
	#elsif ($ch_DEL) {$$region[3] = "DEL"}
	#elsif ($ch_HD)  {@$region[3,4,17] = ("HD", NONE, HD)}
	#else {$$region[3] = ""}
	my $ch_AMP = $$region[3] eq "AMP";
	my $ch_DEL = $$region[3] eq "DEL";
	my $ch_HD  = $$region[3] eq "HD";
	
	if ($$region[3] eq "DEL" && !$$region[19]) {
	#print "@$region\n";
        if ($$region[4] == FATHER) {$$region[4] = MOTHER}
        elsif ($$region[4] == MOTHER) {$$region[4] = FATHER}
    }
    my $p1_AMP = $$region[13] > $AMP;
    my $p2_AMP = $$region[15] > $AMP;
    my $p1_HD  = $$region[13] <= $HD_HI && $$region[13] != NA;
    my $p2_HD  = $$region[15] <= $HD_HI && $$region[15] != NA;
    my $p1_DEL = $$region[13] < $DEL_UPPER && $$region[13] > $HD_HI;
    my $p2_DEL = $$region[15] < $DEL_UPPER && $$region[15] > $HD_HI;
	if ($$region[4] == UNKNOWN) {
		if ($ch_AMP) {
			if ($p1_AMP && !$p2_AMP)    {$$region[4] = FATHER}
			elsif ($p2_AMP && !$p1_AMP) {$$region[4] = MOTHER}
		}
		if ($ch_DEL) {
			if ($p1_DEL && $p2_DEL)    {$$region[4] = FATHER}
			elsif ($p2_DEL && $p1_DEL) {$$region[4] = MOTHER}
		}
	}
    my $P1     = $$region[4] == FATHER;
    my $P2     = $$region[4] == MOTHER;
    my $NONE   = $$region[4] == NONE;
	
    # Likely inherited 
    my $inh    = $ch_AMP && ($p1_AMP && $P1 || $p2_AMP && $P2)
              || $ch_DEL && (($p1_DEL || $p1_HD) && $P1
                  || ($p2_DEL || $p2_HD) && $P2)    
              || $ch_HD && (($p1_HD && ($p2_HD || $p2_DEL))  
			      || ($p2_HD && ($p1_HD || $p1_DEL)));
    
    # Likely inherited with unique copy number state
    my $inh_cn = $ch_DEL && (!$p1_DEL && !$p1_HD && $p2_HD && $P2
                          || !$p2_DEL && !$p2_HD && $p1_HD && $P1)
              || $ch_HD && $p1_DEL && $p2_DEL;
              
    if ($ch_HD && (($p1_DEL || $p1_HD)
        && ($p2_DEL || $p2_HD))) {$$region[4] = BOTH}
    if ($$region[4] == NONE) {$$region[4] = UNKNOWN}
	
    my $state = " ";
    if ($inh) {$state = "INH"}
    if ($inh_cn) {$state = "INH-CN"}
    return($state);                        
}              

sub check_quality {
    # Ends the current analysis if quality checks are not met.
    my ($mi1_rate, $p1_nc_rate, $p2_nc_rate, $ch_nc_rate) = @_;    
    my ($next, $alerts) = (0) x 2;
    # Next file if estimated non-adjacent single Mendelian error rate is above 3%
    if ($mi1_rate >= 0.02) {
        $next = 1;
        print "FAILED - The Mendelian error rate (", 
			sprintf("%.3f",$mi1_rate), " %) is above acceptable levels,\n",
			"indicating a relationship annotation or quality problem.",
			"\nThe analysis cannot be completed!\n" if $verbose;
        print STDERR "Sample $FILENAME : FAILED - The Mendelian error rate ",
			"($mi1_rate %)is above acceptable levels, indicating a ",
			"relationship annotation or quality problem.\n";
        clean_up() unless $BATCH;
    }
 
    # Next file if NC rate is above 10%       
    if ($p1_nc_rate >= 0.1) {
        $next = 1;
        print "FAILED - The Father's NC rate is greater than 10%, indicating ",
			"a quality problem!\nThe analysis cannot be completed!\n\n" 
			if $verbose;
        print STDERR "Sample $FILENAME : FAILED - The Father's NC rate is ",
			"greater than 10%, indicating a quality problem!\nThe analysis ",
			"cannot be completed!\n";
        clean_up() unless $BATCH;            
    }   
    if ($p2_nc_rate >= 0.1) {
        $next = 1;
        print "FAILED - The Mother's NC rate is greater than 10%, indicating ",
			"a quality problem!\nThe analysis cannot be completed!\n\n" 
			if $verbose;
        print STDERR "Sample $FILENAME : FAILED - The Mother's NC rate is ",
			"greater than 10%, indicating a quality problem!\nThe analysis ",
			"cannot be completed!\n";
        clean_up() unless $BATCH;
    }   
    if ($ch_nc_rate >= 0.1) {
        $next = 1;
        print "FAILED - The Child's NC rate is greater than 10%, indicating ",
			"a quality problem!\nThe analysis cannot be completed!\n\n" 
			if $verbose;
        print STDERR "Sample $FILENAME : FAILED - The Child's NC rate is ", 
			"greater than 10%, indicating a quality problem!\nThe analysis ",
			"cannot be completed!\n";
        clean_up() unless $BATCH;
    }   

    if (!$next) {
        # Alert user if NC rate is above 2%.
        my $alerts = 0;  
        if ($p1_nc_rate > $NC_WARNING) {
            print "\nWarning - The Father's NC rate is greater than 2%!\n" 
				if $verbose;
            print OUTPUT_FILE "\nWarning - The Father's NC rate is greater ",
			"than 2%!\n" unless $BATCH;
            print STDERR "Sample $FILENAME : Warning - ",
				"The Father's NC rate is greater than 2%!\n" if $BATCH;				
            $alerts++;
        } 
        if ($p2_nc_rate > $NC_WARNING) {
            print "\nWarning - The Mother's NC rate is greater than 2%!\n" 
				if $verbose;
            print OUTPUT_FILE "\nWarning - The Mother's NC rate is greater ",
			"than 2%!\n" unless $BATCH;
            print STDERR "Sample $FILENAME : Warning - ",
				"The Mother's NC rate is greater than 2%!\n" if $BATCH;				
            $alerts++;
        }   
        if ($ch_nc_rate > $NC_WARNING) {
            print "\nWarning - The Child's NC rate is greater than 2%!\n" 
				if $verbose;
            print OUTPUT_FILE "\nWarning - The Child's NC rate is greater ",
			"than 2%!\n" unless $BATCH;
            print STDERR "Sample $FILENAME : Warning - ",
				"The Child's NC rate is greater than 2%!\n" if $BATCH;				
            $alerts++;
        } 
    }
    print "\n" if $verbose;
    print OUTPUT_FILE "\n" x (4 - $alerts) unless $BATCH;
    return($next); 
}

sub clean_up {
    # Removes files if analysis was not completed due to quality checks
    my $cmd = "rm $OUTPUT_FILE $STATS_FILE";
    system($cmd);
}

sub collapse_region_list {
    my $regions_ref  = $_[0];
	#print  "in collapse\n";
	#print STDERR "in collapse $regions_ref\n";
    #map{print "@{$_}\n"} @$regions_ref;
	
    my $done = 0;
	#if(!$regions_ref) {print "collapse empty\n"}
    return unless ($regions_ref && @$regions_ref && scalar(@$regions_ref) > 1);
	#print STDERR "in collapse2\n";
    my @revised;
    until ($done) {
	#print "done\n";
	#print "@$regions_ref\n";
	    $done = 1 if !@$regions_ref;
        OUTER: for (my $i = 0; $i <= $#$regions_ref; $i++) {
		#print "outer\n";
            my $reg1 = $$regions_ref[$i];
			#print  "reg1 @$reg1\n";
            for (my $j = 0; $j <= $#$regions_ref; $j++) {
                if ($i == $j) {next}
                my $reg2 = $$regions_ref[$j];
						#	print  "reg2 @$reg2\n";

				if ($$reg1[17] == HD) {next unless ($$reg2[17] == HD)}
                if (overlap(@$reg1[1,2], @$reg2[1,2])->[0]) {
				#print "in if\n";
                    for (my $k = 0; $k < @$regions_ref; $k++) {
                        push(@revised, $$regions_ref[$k]) if ($k != $i && $k != $j);
                    }
                    @$regions_ref = @revised; 
                    @revised = ();
									#print "to splice\n";

                    my @overlap = @{splice_overlap_regions($reg1, $reg2)};
									#print "back from splice @overlap\n";

                    map {$$_[5] = count_SNPs(@$_[START,STOP], \%SNP_BY_POS)} @overlap;
                    push(@$regions_ref, @overlap);
														#print "here @$regions_ref\n";

                    last OUTER;
                }
            }
			#print "$i scalar(@$regions_ref)\n";
            $done = 1 if (!@$regions_ref || $i == $#$regions_ref);
        }
		#$done = 1 if (@$regions_ref < 2);
        map {for my $x (4..15,18) {$$_[$x] ||= 0}} @$regions_ref if @$regions_ref;
    }
	#map{print  "leaving collapse @{$_}\n"} @$regions_ref if @$regions_ref;
}

sub count_lines {
    # Retrieves and returns the total number of lines in the input file
    my $filename = $_[0];
       my $num_lines = `wc -l $filename | awk \'{print \$1}\'`;
    return($num_lines);
}

sub count_std_mi1_overlap {
	my $regions_ref  = $_[0];
	map{$$_[18] = 0} @$regions_ref;
	for (my $i = 0; $i < @$regions_ref; $i++) {
		my $reg1 = $$regions_ref[$i];
		if ($$reg1[17] == POD) {
			for (my $j = 0; $j < @$regions_ref; $j++) {
				next if ($i == $j);
				my $reg2 = $$regions_ref[$j];
				if ($$reg2[17] == MI1) {
					if (overlap($$reg1[1], $$reg1[2], $$reg2[1], 
						$$reg2[2])->[0]) {$$reg1[18]++}
				}
			}
		}
	}
}

sub count_nonoverlap_regions {
	# Counts and stores number of non-overlapping windows with inf snps >= 3.
	# Also returns refined mi1 and informative SNP counts.
	my ($inf_SNP_ct, $mi1_ct, $snp_ct, $curr_inf_ct, $prev_p1_hd, 
		$prev_p2_hd, $prev_ch_hd, $prev_MI1, $prev_outlier, $last_inf_snp
		) = (0) x 10;
	my $size = scalar(@CURR_CHR_REFS);
	
	# The following are *_stats array indices
	my ($CHR, $START, $MAX_SIZE, $TOTAL, $NON_OVER_REF, $MI1_CT, 
		$INF_SNP_CT, $P1_HD_CT, $P2_HD_CT, $CH_HD_CT, $HET_MI1_CT,
		$BAF_OUTLIERS, $AB_CT) = (0..12);
	my (@rolling_array, @prev_ends, @nonoverlap, @stats);
	@stats[$CHR,$START] = @{$CURR_CHR_REFS[0]}[1,2];
	$stats[$_] = 0 for (2..12); 
	for (0..$SIZE_REGION + 1) {$prev_ends[$_] = $nonoverlap[$_] = 0}
	my $outlier_ct;
	my @overlapping;
	for (0..$SIZE_REGION + 1) {$overlapping[$_] = 0}
	for (my $i = 0; $i < @CURR_CHR_REFS; $i++) {
		next if (!$CURR_CHR_REFS[$i]->[3] || !$CURR_CHR_REFS[$i]->[6] 
			|| !$CURR_CHR_REFS[$i]->[9]);
	    next if ((substr($CURR_CHR_REFS[$i]->[0], 0, 2) =~ /cn/i)
			  || substr($CURR_CHR_REFS[$i]->[4],  0, 1) eq "N" 
			  || substr($CURR_CHR_REFS[$i]->[7],  0, 1) eq "N" 
			  || substr($CURR_CHR_REFS[$i]->[10], 0, 1) eq "N");
		$snp_ct++;
		# if ($CURR_CHR_REFS[$i]->[10] > $AB_UPPER_BOUND 
			# && $CURR_CHR_REFS[$i]->[10] < $BB_BOUND
			# || $CURR_CHR_REFS[$i]->[10] > $AB_LOWER_BOUND 
			# && $CURR_CHR_REFS[$i]->[10] < $AA_BOUND) {
			# $outlier_ct++;
			# }
        my $results = find_informative_snps($CURR_CHR_REFS[$i]);
        $stats[$INF_SNP_CT] += $$results[2];
		# Count nonconsecutive SNPs with MI1 (or imputed)
		if ($$results[3]) {
			$stats[$HET_MI1_CT]++ if ($$results[4] == NONE);
			$stats[$MI1_CT]++ unless $prev_MI1;
			$stats[$MI1_CT]-- if ($prev_MI1 == 1);
			$prev_MI1++;
		}
		else {$prev_MI1 = 0}
		if ($$results[5]) {
			$stats[$BAF_OUTLIERS]++ unless $prev_outlier;
			$stats[$BAF_OUTLIERS]-- if ($prev_outlier == 1);
			$prev_outlier++;
		}
		# else {$prev_outlier = 0 if ($CURR_CHR_REFS[$i]->[9] ne "AA" 
			# && $CURR_CHR_REFS[$i]->[9] ne "BB" 
			# && $CURR_CHR_REFS[$i]->[10] < $BB_BOUND 
			# && $CURR_CHR_REFS[$i]->[10] > $AA_BOUND);
		# }
		else {$prev_outlier = 0 if ($CURR_CHR_REFS[$i]->[9] eq "AB")}


		#if ($CURR_CHR_REFS[$i]->[-1] == OUTLIER) {$stats[$BAF_OUTLIERS]++}
		# $stats[$AB_CT]++ if ($CURR_CHR_REFS[$i]->[9] ne "AA" 
			# && $CURR_CHR_REFS[$i]->[9] ne "BB" 
			# && $CURR_CHR_REFS[$i]->[10] < $BB_BOUND 
			# && $CURR_CHR_REFS[$i]->[10] > $AA_BOUND);
		$stats[$AB_CT]++ if ($CURR_CHR_REFS[$i]->[9] eq "AB");

		# Count nonconsecutive SNPs with LRR < HD threshold		
		if ($CURR_CHR_REFS[$i]->[5] < $HD_HI) {
			$stats[$P1_HD_CT]++ unless ($prev_p1_hd);
			$stats[$P1_HD_CT]-- if ($prev_p1_hd == 1);
			$prev_p1_hd++;
        }
        else {$prev_p1_hd = 0} 
		if ($CURR_CHR_REFS[$i]->[8] < $HD_HI) {
			$stats[$P2_HD_CT]++ unless ($prev_p2_hd);
			$stats[$P2_HD_CT]-- if ($prev_p2_hd == 1);
			$prev_p2_hd++;
        }
        else {$prev_p2_hd = 0} 
		if ($CURR_CHR_REFS[$i]->[11] < $HD_HI) {
			$stats[$CH_HD_CT]++ unless ($prev_ch_hd);
			$stats[$CH_HD_CT]-- if ($prev_ch_hd == 1);
			$prev_ch_hd++;
        }
        else {$prev_ch_hd = 0} 

		# Track nonoverlapping windows with >= 3 informative SNPs 
		if ($snp_ct < $SIZE_REGION) {
			push(@rolling_array, $$results[2]);
			$curr_inf_ct = $stats[$INF_SNP_CT];
			if ($$results[2]) {
				$last_inf_snp = $CURR_CHR_REFS[$i]->[POS];
			}				
		}
		else {
			$curr_inf_ct-- if $rolling_array[0];
			if ($$results[2]) {
				$curr_inf_ct++;
				$last_inf_snp = $CURR_CHR_REFS[$i]->[POS];
			}
			if ($curr_inf_ct >= 3) {
				$overlapping[$curr_inf_ct]++;
				for (3..$curr_inf_ct) {
					if($CURR_CHR_REFS[$i - $SIZE_REGION]->[POS]
						> $prev_ends[$_]) {$nonoverlap[$_]++}
				}
				for (3..$curr_inf_ct) {$prev_ends[$_] = $last_inf_snp}
			}
			shift(@rolling_array);
			push(@rolling_array, $$results[2]);
		}
	}
	#print "het_mi1s $stats[$HET_MI1_CT]\n";
	return() unless ($snp_ct >= $SIZE_REGION); 
	$stats[$NON_OVER_REF] = \@nonoverlap;
	for (0..$#nonoverlap) {$stats[$MAX_SIZE] = $_ if $nonoverlap[$_]};
	for (0..$#nonoverlap) {$stats[$TOTAL] += $nonoverlap[$_] if $_ >= 6};
	#for (0..$#nonoverlap) {$stats[$TOTAL] += $nonoverlap[$_]};
	
	#print "$CURR_CHR_REFS[0]->[CHR] $CURR_CHR_REFS[0]->[POS] $outlier_ct $snp_ct ", $outlier_ct / $snp_ct, "\n";
	print OUTPUT_FILE "@overlapping\n";
	return(\@stats);
}
		
sub count_SNPs {
    my ($position1, $position2, $hash_ref) = @_;
    my $num_snps  = 0;
    if (!$position1 || !$position2) {return($num_snps)}
    my $value1 = $$hash_ref{$position1};
    my $value2 = $$hash_ref{$position2};   
    if (!$value1 || !$value2) {
        foreach my $key (keys %$hash_ref) {
            $num_snps++ if ($key >= $position1 && $key <= $position2);
        }
    }
    else {$num_snps = $value2 - $value1 + 1}
    return($num_snps);
}    

sub cpu_time {
    # Retrieves and returns the time and date as set on the current machine
    my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my @week_days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    my ($second, $minute, $hour, $day_of_month, $month, $year_offset, 
		$day_of_week, $day_of_year, $daylight_savings) = localtime();
    my $year = 1900 + $year_offset;
    $day_of_year = $daylight_savings = 0; # to make -w happy
    my $THE_TIME = join("", $hour, ":", $minute,":", $second, " ",
        $week_days[$day_of_week], " ", $months[$month], $day_of_month,
		" ", $year);
    return($THE_TIME);
}

sub cusum {
# Use cumulative sums of LRR difference from a threshold to detect 
# critical points indicating start and stop of cnv.
# For efficiency, the LRR medians and thresholds for deletions are
# inverted to allow for detection of maximum critical points.
#print "in cusum\n";
    my ($region_ref, $array_ref, $reg_index, $median, $start_index, $stop_index,
        $critical_index) = @_;
    my @data_refs = @$array_ref; 
	my $data;
    my ($prior_median, $post_median, $begin_cusum, $end_cusum, $min,
        $max_index, $min_index, $found_start_index, $found_end, $cusum, $hd, 
        $del, $amp, $other, $thresh, $local_thresh, $start_crit_index, 
        $end_crit_index, $adj_abn, $father, $mother) = (0) x 21;
    my $max = undef;  
	
	# For region of 2 snps begin cusum from both ends
    if ($$region_ref[5] == 2) {
        ($start_crit_index, $end_crit_index)  = ($stop_index, $start_index);
    }
	elsif ($$region_ref[5] == 3) {
		if ($critical_index) {$start_crit_index = $end_crit_index = $critical_index}
		else {($start_crit_index, $end_crit_index)  = ($stop_index, $start_index)}
	}
    else {$start_crit_index = $end_crit_index = $critical_index}
	
	if (!$start_crit_index) {print STDERR "nostartcrit @$region_ref\n"}
    my (@prior, @post);
    my $count = 1;
    my $extend = $$region_ref[5] * 5;
    $extend = 25 if $extend > 25;
	#print "begin cusum $reg_index $start_crit_index $end_crit_index @$region_ref\n";
	if    ($$region_ref[4] == FATHER) {$father = 1}
	elsif ($$region_ref[4] == MOTHER) {$mother = 1}
		
    if ($REFINING) {
		if ($reg_index == 10)   {$local_thresh = $LOCAL_CH_mBAF_MED}
		#elsif ($father) {$local_thresh = -$LOCAL_P1_LRR_MED}        
        #elsif ($mother) {$local_thresh = -$LOCAL_P2_LRR_MED}
        elsif ($reg_index == 3) {$local_thresh = -$LOCAL_P1_LRR_MED}        
        elsif ($reg_index == 4) {$local_thresh = -$LOCAL_P2_LRR_MED}
        else                    {$local_thresh = -$LOCAL_CH_LRR_MED}
    }
    else {
		if ($reg_index == 10)   {$local_thresh = $CH_mBAF_MED}
        #elsif ($father) {$local_thresh = -$P1_LRR_MED}        
        #elsif ($mother) {$local_thresh = -$P2_LRR_MED}
		elsif ($reg_index == 3) {$local_thresh = -$P1_LRR_MED}        
        elsif ($reg_index == 4) {$local_thresh = -$P2_LRR_MED}
        else                    {$local_thresh = -$CH_LRR_MED}
    }
    
    if ($$region_ref[3] eq "HD") {$hd  = 1}
    elsif ($$region_ref[3] eq "DEL")  {$del = 1}
    elsif ($$region_ref[3] eq "AMP")  {$amp = 1}
    else {$other = 1}
    $median = -$median if ($hd || $del);
	
    if ($median > -$HD_HI * 2) {$thresh = -$HD_HI}
    else {$thresh = $median - (($median - $local_thresh) / 2)}
	
    # If a non-hd region is large, start cusum $MIN_POD number of informative SNPs
    # from each end.
    if ($$region_ref[6] > $MIN_POD * 3) {
        if ($$region_ref[17] == POD) {
			my ($hash_ref, $inf_array_ref);
				#print STDERR "cusum to find\n";
				#		print STDERR "@$region_ref[START,STOP]\n";

            my ($start, $stop) 
                = find_closest_indices(@$region_ref[START,STOP,16]);
			if    ($father) {$inf_array_ref = \@P1_INF_SNPS}
			elsif ($mother) {$inf_array_ref = \@P2_INF_SNPS}
			if ($other)  {$hash_ref = \%BAF_BY_POS}
			else {$hash_ref = \%LRR_BY_POS}
            my $start_pos = $$inf_array_ref[$start + $MIN_POD - 1]->[POS];
            my $end_pos   = $$inf_array_ref[$stop  - $MIN_POD]->[POS];	
			#print STDERR "cusum2 to find\n";
			#						print STDERR "$start_pos, $end_pos\n";

			($start_crit_index, $end_crit_index) 
                = find_closest_indices($start_pos, $end_pos, $hash_ref);
        }
    }

	$count = 0;
    if ($start_index - $extend >= 0) {
        for (my $x = 1; @prior < $extend; $x++) {
			$count++;
			my $data = $data_refs[$start_index - $x]->[$reg_index];
			if ($other) {
				next if ($data < $AA_BOUND || $data > $BB_BOUND);
				$data = 1 - $data if ($data < $LOCAL_CH_BAF_MEAN);
			}
            push(@prior, $data);
        }
        $prior_median = median(\@prior);
        $prior_median = -$prior_median if ($hd || $del);
        # Determine if adjacent region is a unique abnormal region
        # with a greater median
        $adj_abn = 1 if ($prior_median >= $median * 2);
        $begin_cusum = $start_index - $count;
        $count  = 0;
		FIND_PRIOR_NORM:
        until ($prior_median < $thresh) {
            for (my $x = $count; $x <= $extend + $count; $x++) {
                if ($start_index - ($extend + $x) < 0) {
                   $$region_ref[START] = $data_refs[0]->[POS];
                   $found_start_index = 1;
                   last FIND_PRIOR_NORM;
                }
                
				shift(@prior);
				$data = $data_refs[$start_index - ($extend + $x)]->[$reg_index];
				if ($other) {
					next if ($data < $AA_BOUND || $data > $BB_BOUND);
					$data = 1 - $data if ($data < $LOCAL_CH_BAF_MEAN);
				}
				push(@prior, $data);
			}
            $prior_median = median(\@prior);
			$prior_median = -$prior_median if ($hd || $del);
			$begin_cusum = $start_index - ($extend + $count);
			$count += $extend;
        }    
    }
    else {
        for (my $x = 1; $start_index - $x >= 0; $x++) {
			$data = $data_refs[$start_index - $x]->[$reg_index];
			if ($other) {
				next if ($data < $AA_BOUND || $data > $BB_BOUND);
				$data = 1 - $data if ($data < $LOCAL_CH_BAF_MEAN);
			}
            push(@prior, $data);
        }
        $prior_median = median(\@prior);
        $prior_median = -$prior_median if ($prior_median != NA && ($hd || $del));       
        if ($prior_median > $thresh) {
            $$region_ref[START] = $data_refs[0]->[POS];
            $found_start_index = 1;
        }
        else {$begin_cusum = 0}
    }
	
    unless ($found_start_index) {
		# If HD region has complete separation from surrounding normal region
		# adjust thresh to ensure accurate change point detection.
		# if ($median > -$HD_HI * 2) {
			# my @norm;
			# for (my $i = $begin_cusum; $i < $start_index; $i++) {
				# push(@norm, $data_refs[$i]->[$reg_index]);
			# }
			# my @sort_norm = sort { $a <=> $b } @norm;
			# my $min_norm = $sort_norm[0];
			# my @abn;
			# for (my $i = $start_index; $i <= $stop_index; $i++) {
				# push(@abn, $data_refs[$i]->[$reg_index]);
			# }
			# my @sort_abn = sort { $a <=> $b } @abn;
			# my $max_abn = $sort_abn[-1];
			# if ($min_norm > $max_abn) {
				# $thresh = -($min_norm + (($max_abn - $min_norm) / 2));
			# }
		# }
        for (my $x = $start_crit_index; $x >= $begin_cusum; $x--) {
            my $value = $data_refs[$x]->[$reg_index];
            # Minimize affect of large outliers in del region
            #if ($del && $median < $HD_HI) {$value = 0 if $value > $HD_HI}
			
			next if $del && $median < $HD_HI && $value > $HD_HI; 
			if ($other) {
				next if ($value < $AA_BOUND || $value > $BB_BOUND);
				$value = 1 - $value if ($value < $LOCAL_CH_BAF_MEAN);
			}
			# Minimize affect of large outliers in upd region
			next if ($other && $value > $median + 0.1); 
			
            $value = -$value if ($hd || $del);
            $cusum += $value - $thresh;
            if (!defined($max)) {
                $max = $min = $cusum;
                $max_index = $min_index = $x;
            }
            if ($cusum > $max) {
                $max = $cusum;
                $max_index = $x;
            }
            elsif ($cusum < $min) {
                $min = $cusum;
                $min_index = $x;
            }

        }
        # If adjacent SNPs are likely to be in another abnormal region keep  
        # initial boundary unless there is a detectable cusum minimum
        if ($adj_abn) {
            if ($min_index != $start_crit_index 
                && $min_index < $begin_cusum + ($extend / 5)) { 
                $$region_ref[START] = $data_refs[$min_index]->[POS];
            }
        }
        elsif ($max_index != $start_crit_index) { 
            $$region_ref[START] = $data_refs[$max_index]->[POS];
        }
        if ($$region_ref[17] == HD && $reg_index != 3 && $reg_index != 4) {
            if ($$region_ref[START] < $data_refs[$start_index]->[POS]) {
                $$region_ref[START] = $data_refs[$start_index]->[POS];
            }
        }
    }  
    ($max_index, $min, $min_index, $cusum, $adj_abn) = (0) x 5;
	$max = undef;
    $count = 0;

    if ($stop_index + $extend <= $#data_refs) { 
        for (my $x = 1; @post < $extend; $x++) {
			$count++;	
			#print $data_refs[$stop_index + $x], "\n";
			$data = $data_refs[$stop_index + $x]->[$reg_index];
			if ($other) {
				last if ($stop_index + $x >= $#data_refs);
				next if ($data < $AA_BOUND || $data > $BB_BOUND);
				$data = 1 - $data if ($data < $LOCAL_CH_BAF_MEAN);
			}
            push(@post, $data); 
			last if ($stop_index + $x >= $#data_refs);
        }
		
        $post_median = median(\@post);
        $post_median = -$post_median if ($hd || $del);
        # Determine if adjacent region is a unique abnormal region
        # with a greater median
        if ($post_median >= $median * 2) {
		$adj_abn = 1;
		}        
        $end_cusum = $stop_index + $count;
		FIND_POST_NORM:
        until ($post_median < $thresh) {
            for (my $x = $count; $x <= $extend + $count; $x++) {
                if ($stop_index + $extend + $x > $#data_refs) {
                   $$region_ref[STOP] = $data_refs[-1]->[POS];
                   $found_end = 1;
                   last FIND_POST_NORM;
                }
                shift(@post);
				$data = $data_refs[$stop_index + ($extend + $x)]->[$reg_index];
				if ($other) {
					next if ($data < $AA_BOUND || $data > $BB_BOUND);
					$data = 1 - $data if ($data < $LOCAL_CH_BAF_MEAN);
				}
				push(@post, $data);
			}
			$post_median = median(\@post);
            $post_median = -$post_median if ($hd || $del);               
            $end_cusum = $stop_index + ($extend + $count);                            
            $count += $extend;
        }
    }
    else {
        for (my $x = 1; $x <= $#data_refs - $stop_index; $x++) {
			$data = $data_refs[$stop_index + $x]->[$reg_index];
			if ($other) {
				next if ($data < $AA_BOUND || $data > $BB_BOUND);
				$data = 1 - $data if ($data < $LOCAL_CH_BAF_MEAN);
			}
            push(@post, $data);                    
        }
        $post_median = median(\@post);
        $post_median = -$post_median if ($prior_median != NA && ($hd || $del));               
        if ($post_median > $thresh) {
            $$region_ref[STOP] = $data_refs[-1]->[POS];
            $found_end = 1;
        }
        else {$end_cusum = $#data_refs}
    }

    unless ($found_end) {
		# If HD region has complete separation from surrounding normal region
		# adjust thresh to ensure accurate change point detection.
		# if ($median > -$HD_HI * 2) {
			# my @norm;
			# for (my $i = $end_cusum; $i > $stop_index; $i--) {
				# push(@norm, $data_refs[$i]->[$reg_index]);
			# }
			# my @sort_norm = sort { $a <=> $b } @norm;
			# my $min_norm = $sort_norm[0];
			# my @abn;
			# for (my $i = $start_index; $i <= $stop_index; $i++) {
				# push(@abn, $data_refs[$i]->[$reg_index]);
			# }
			# my @sort_abn = sort { $a <=> $b } @abn;
			# my $max_abn = $sort_abn[-1];
			# if ($min_norm > $max_abn) {
				# $thresh = -($min_norm + (($max_abn - $min_norm) / 2));
			# }
		# }

        for (my $x = $end_crit_index; $x <= $end_cusum; $x++) {
            my $value = $data_refs[$x]->[$reg_index];
            # Minimize affect of large outliers in del region
            #if ($del && $median < $HD_HI) {$value = 0 if $value > $HD_HI}
			next if $del && $median < $HD_HI && $value > $HD_HI; 
			if ($other) {
				next if ($value < $AA_BOUND || $value > $BB_BOUND);
				$value = 1 - $value if ($value < $LOCAL_CH_BAF_MEAN);
			}

			# Minimize affect of large outliers in upd region
			next if ($other && $value > $median + 0.1);  
			
            $value = -$value if ($hd || $del);
            $cusum += $value - $thresh;
            if (!defined($max)) {
                $max = $min = $cusum;
                $max_index = $min_index = $x;
            }
            if ($cusum > $max) {
                $max = $cusum;
                $max_index = $x;
            }
            elsif ($cusum < $min) {
                $min = $cusum;
                $min_index = $x;
            }
        }
        # If adjacent SNPs are likely to be in another abnormal region keep  
        # initial boundary unless there is a detectable cusum minimum
        if ($adj_abn) {
            if ($min_index != $end_crit_index 
                && $min_index < $end_cusum - ($extend / 5)) { 
                $$region_ref[STOP] = $data_refs[$min_index]->[POS];
            }
			elsif ($max_index != $end_crit_index) {
				$$region_ref[STOP] = $data_refs[$max_index]->[POS];
			}
        }
        elsif ($max_index != $end_crit_index) { 
            $$region_ref[STOP] = $data_refs[$max_index]->[POS];
        }
        if ($$region_ref[17] == HD && $reg_index != 3 && $reg_index != 4) {
            if ($$region_ref[STOP] > $data_refs[$stop_index]->[POS]) {
                $$region_ref[STOP] = $data_refs[$stop_index]->[POS];
            }
        }
    }
	#print "leaving cusum\n";
	#print "end cusum @$region_ref\n";

    return($region_ref);
}

sub define_POD_region {
    # Combines and calculates descriptive info for a detected POD region.
    my ($ref, $region1, $last_sig_SNP) = @_;
    my (@abn_region, $array_ref);
    for (my $i = 0; $i < @$ref; $i++) {
        if ($i == 0) {@abn_region[0,START,STOP] = @{$$ref[$i]}[0,START,STOP]}
        elsif ($i == $last_sig_SNP) {$abn_region[STOP] = $$ref[$i]->[STOP]}
        if ($i == @$ref - 1) {
            if ($region1) {
                @abn_region[3,4,16,17] = (" ",FATHER,\%P1_INF_SNP_BY_POS, POD);
                $array_ref = \@P1_INF_SNPS;
            } 
            else {
                @abn_region[3,4,16,17] = (" ",MOTHER,\%P2_INF_SNP_BY_POS, POD);
                $array_ref = \@P2_INF_SNPS;
            }
        }
    }   
    for (4..17) {$abn_region[$_] ||= 0}
				#print OUTPUT_FILE "definePOD to find\n";
#print OUTPUT_FILE "@abn_region[START,STOP]\n";
    my ($closest_start_ind, $closest_stop_ind) 
		= find_closest_indices(@abn_region[START,STOP,16]);
    $abn_region[START] = $$array_ref[$closest_start_ind]->[POS];
    $abn_region[STOP]  = $$array_ref[$closest_stop_ind]->[POS];
    $abn_region[6] = count_SNPs(@abn_region[START,STOP,16]);
		#print "in define @abn_region\n";

    @abn_region = @{evaluate_region(\@abn_region)};
    $abn_region[6] = count_SNPs(@abn_region[START,STOP,16]);
	#print "leaving define @abn_region\n";
    return(\@abn_region);   
}                
    
sub define_mi1_regions {
    my $mi1_ref = $_[0];
    map {@{$_}[16,17] = (\%MI1_BY_POS, MI1)} @$mi1_ref;
    map {for my $x (7..15,18) {$$_[$x] ||= 0}} @$mi1_ref;
	#print  "indefinemi1\n";
	#map{print  "@{$_}\n"} @$mi1_ref if @$mi1_ref;
    map {$_ = evaluate_region($_)} @$mi1_ref;
		#print  "after eval\n";
	#map{print  "@{$_}\n"} @$mi1_ref if @$mi1_ref;
	#print "mi1tocollapse\n";
    collapse_region_list(\@$mi1_ref);
		#print  "after collapse\n";
	#map{print  "@{$_}\n"} @$mi1_ref if @$mi1_ref;
	
    map {$_ = evaluate_region($_) unless $$_[6]} @$mi1_ref;
	my @mi1_regions;
	map {push(@mi1_regions, $_) if ($$_[6] >= $MIN_MI1)} @$mi1_ref;
    map {calculate_region_stats($_)} @mi1_regions;
    return (\@mi1_regions);
}

sub define_PODcr_regions {
	my $outlier_regs = $_[0];
	my $array_ref;
	map {@{$_}[16,17] = (\%BAF_OUTLIERS_BY_POS, PODcr)}@$outlier_regs;
	
    map {for my $x (7..15,18) {$$_[$x] ||= 0}} @$outlier_regs;
	#print STDERR "indefinePODcr\n";
	map{print STDERR "@{$_}\n"} @$outlier_regs if @$outlier_regs;
	#print STDERR "podcr to collapse\n";
	collapse_region_list(\@$outlier_regs);
		#print STDERR "back from collapse\n";

	for my $ref (@$outlier_regs) {
					#print STDERR "definePODcr to find\n";
					#print STDERR "@$ref[START,STOP]\n";

		my ($closest_start_ind, $closest_stop_ind) 
			= find_closest_indices(@$ref[START,STOP,16]);
			#print "@$ref[START,STOP]\n";
			#print "closest $closest_start_ind $closest_stop_ind\n";
			#print "$$array_ref[$closest_start_ind]->[POS] $$array_ref[$closest_stop_ind]->[POS]\n";
		$$ref[START] = $BAF_OUTLIERS[$closest_start_ind]->[POS];
		$$ref[STOP]  = $BAF_OUTLIERS[$closest_stop_ind]->[POS];
		$$ref[6] = count_SNPs(@$ref[START,STOP,16]);
	}
	
	#print STDERR "back from ref\n";
	#map{print STDERR "@{$_}\n"} @$outlier_regs if @$outlier_regs;
	#map{print "@{$_}\n"} @$outlier_regs if @$outlier_regs;
	
    #print "sending to eval\n";
    map {$_ = evaluate_region($_)} @$outlier_regs;
		#print STDERR "cryptic after eval\n";
	#map{print STDERR "@{$_}\n"} @$outlier_regs if @$outlier_regs;
		#print STDERR "podcr to collapse2\n";

    collapse_region_list(\@$outlier_regs);
		#print STDERR "cryptic after collapse\n";
	#map{print STDERR "@{$_}\n"} @$outlier_regs if @$outlier_regs;
	
    map {$_ = evaluate_region($_) unless $$_[6]} @$outlier_regs;
	#my @return;
	#map {push(@return, $_) if ($$_[6] >= $MIN_POD)} @$outlier_regs;
    #map {calculate_region_stats($_)} @return;
    #return (\@return);
	map {calculate_region_stats($_)} @$outlier_regs;
	
	# To avoid regions due to BAF genomic waves, determine if the distribution
	# of AB calls on is as expected above the upper thresh and below the lower
	my @results;
	foreach my $ref (@$outlier_regs) {
	#print STDERR "PODcr @$ref\n";
		my ($upper, $lower) = (0) x 2;
		my ($reg_start, $reg_stop) 
			= find_closest_indices(@$ref[START,STOP], \%BAF_OUTLIERS_BY_POS);
		#print STDERR "$CURR_CHR_REFS[$reg_start]->[POS] $CURR_CHR_REFS[$reg_stop]->[POS]\n";
		if (defined($reg_start) && defined($reg_stop)) {
			for (my $j = $reg_start; $j <= $reg_stop; $j++) {
				if ($BAF_OUTLIERS[$j]->[10] > $AB_UPPER_BOUND) {$upper++}
				elsif ($BAF_OUTLIERS[$j]->[10] < $AB_LOWER_BOUND) {$lower++}
			}
			my $prob = 0;	
#print STDERR "lower upper $lower $upper\n";	
			if ($lower == $upper) {$prob = 1}
			else {
				my ($n, $k, $p) = ($lower + $upper, $lower, 0.5);
				my ($x, $y);
				if ($k <= $n / 2) {($x, $y) = (0, $k)}
				else {($x, $y) = ($k, $n)}
				for ($x..$y) {
					my $bin_coeff = calculate_bin_coeff($n, $_);
					$prob += $bin_coeff * ($p**$_) * ((1 - $p)**($n - $_));
				}
				$prob = 2 * $prob unless $lower == $upper;
			}
			my $sidak_alpha = 1 - (1 - $ALPHA)**(1 / scalar(@$outlier_regs));
						#print STDERR "$prob $sidak_alpha\n";

			if ($prob > $sidak_alpha) {push(@results, $ref)}
		}
	}
    return (\@results);
}

sub detect_hom_del { 
    my (@p1_hd_indices, @p2_hd_indices, @child_hd_indices); 
    my ($prev_snp, $next_snp, $size);
    # Locate SNPs indicating a homozygous deletion
    # Will be labeled as parent lacking information content and 
    # switched to parent contributing before printing
    for (my $i = 0; $i < @LRR_REFS; $i++) {
        my $line_ref  = $LRR_REFS[$i];
        my ($p1_LRR, $p2_LRR, $child_LRR) = @$line_ref[3..5];
        if ($child_LRR <= $HD_HI) {
            $CH_HD_BY_POS{$$line_ref[POS]} = keys %CH_HD_BY_POS;
            push(@child_hd_indices, [$i, NONE]);
        }
        if ($p1_LRR <= $HD_HI) {
            $P2_HD_BY_POS{$$line_ref[POS]} = keys %P2_HD_BY_POS;
            push(@p2_hd_indices, [$i, MOTHER]);
        }
        if ($p2_LRR <= $HD_HI) {
            $P1_HD_BY_POS{$$line_ref[POS]} = keys %P1_HD_BY_POS;
            push(@p1_hd_indices, [$i, FATHER]);
        }
    }
     
    # Extend each SNP into regions
    my @p1_hd = @{extend_hd(\@p1_hd_indices)};
    my @p2_hd = @{extend_hd(\@p2_hd_indices)};
    my @ch_hd = @{extend_hd(\@child_hd_indices)};
    
    #my @temp;
    if (@p1_hd) {
	    map {for my $x (4..18) {$$_[$x] ||= 0}} @p1_hd;
        map {$_ = evaluate_region($_, FATHER)} @p1_hd;
        #map {push(@temp, evaluate_region($_))} @p1_hd;
        #push(@p1_hd, @temp);
		map {$_ = evaluate_region($_)} @p1_hd;
        #@temp = ();
		#print STDERR "p1hd to collapse\n";
        collapse_region_list(\@p1_hd);
        map {@{$_}[16,17] = (\%P1_HD_BY_POS, HD)} @p1_hd;
        map {$_ = evaluate_region($_) unless $$_[6]} @p1_hd;
    }

    if (@p2_hd) {
	    map {for my $x (4..18) {$$_[$x] ||= 0}} @p2_hd;
        map {$_ = evaluate_region($_, MOTHER)} @p2_hd;
        #map {push(@temp, evaluate_region($_))} @p2_hd;
        #push(@p2_hd, @temp);
		map {$_ = evaluate_region($_)} @p2_hd;	
		#print STDERR "p2hd to collapse\n";
		
        collapse_region_list(\@p2_hd);
        map {@{$_}[16,17] = (\%P2_HD_BY_POS, HD)} @p2_hd;
        map {$_ = evaluate_region($_) unless $$_[6]} @p2_hd;
    }
        
    if (@ch_hd) {
	    map {for my $x (4..18) {$$_[$x] ||= 0}} @ch_hd;
        map {$_ = evaluate_region($_)} @ch_hd;
				#print STDERR "chhd to collapse\n";

        collapse_region_list(\@ch_hd);
        map {@{$_}[16,17] = (\%CH_HD_BY_POS, HD)} @ch_hd;
        map {$_ = evaluate_region($_) unless $$_[6]} @ch_hd;
    }
    # Correct any overlap between trio members
    my (@combine_regions, @return_regions);
    push(@combine_regions, @p1_hd, @p2_hd, @ch_hd);
    if (@combine_regions) { 
			#print STDERR "combinehd to collapse\n";

	    collapse_region_list(\@combine_regions);
        map {$_ = evaluate_region($_) unless $$_[6]} @combine_regions;
		
        # Count SNPs and informative SNPS
        my $hd_thresh;
        foreach my $ref (@combine_regions) {
            $$ref[5] = count_SNPs(@$ref[1,2], \%SNP_BY_POS);
            if ($$ref[4] == FATHER) {
                $$ref[16] = \%P1_HD_BY_POS;
                $$ref[6] = count_SNPs(@$ref[1,2,16]);
                $hd_thresh = $MIN_P2_HD;
            }
            elsif ($$ref[4] == MOTHER) {
                $$ref[16] = \%P2_HD_BY_POS;
                $$ref[6] = count_SNPs(@$ref[1,2,16]);
                $hd_thresh = $MIN_P1_HD;
            }
            elsif ($$ref[4] == NONE) {
                $$ref[16] = \%CH_HD_BY_POS;
                $$ref[6] = count_SNPs(@$ref[1,2,16]);
                $hd_thresh = $MIN_CH_HD;
            }
			elsif ($$ref[4] == BOTH) {
				if ($MIN_P1_HD <= $MIN_P2_HD) {
					$$ref[16] = \%P2_HD_BY_POS;
					$hd_thresh = $MIN_P1_HD;
				}
				else {
					$$ref[16] = \%P1_HD_BY_POS;
					$hd_thresh = $MIN_P2_HD;
				}
            }
            else {$$ref[6] = NA}
            if ($$ref[6] >= $hd_thresh) {push(@return_regions, $ref)}
        }        
    }
    map {@{$_}[17] = HD} @return_regions if @return_regions;

    # Calculate region stats
    map {calculate_region_stats($_)} @return_regions if @return_regions;
    return(\@return_regions);
}            
  
sub determine_R_status {
    # Determines if the R script is still running and makes one attempt to 
    # restart the Rscript.
    my $R_ATTEMPTS = $_[0];
    my $R_status = `ps -a | grep $R_PID`;
    if (!$R_status) { 
        if (!$R_ATTEMPTS) {
            $R_ATTEMPTS++;
            start_R();
            sleep(1);
            ($R_PID_FILE, $R_PID_FILENAME) = tempfile(UNLINK => 1);
            if(!get_R_PID()) {$R_ATTEMPTS++}
        }            
        if ($R_ATTEMPTS > 1) {
            print "A fatal error has occurred in communication with the R",
                "software.\n" if $verbose;
            print STDERR "A fatal error has occurred in communication with ",
				"the R software.\n";				
        }        
    }
    return($R_ATTEMPTS);
} 

# sub estimate_additional_instances {
    # # Estimate number of adjacent SNPs per genome due to chance
    # # as determined by the rate of non-adjacent SNPs.
    # # Iterate until the number SNPs 
    # # < the significance threshold. 
    # my ($n, $p, $thresh) = @_;
	# my $add_SNPs = my $stop = 0;
	# my $k = 2;
    # until ($stop) {
	# ???
		# if ($prob > $thresh) {$add_SNPs += ($prob * (1 / $thresh) * $k)}
		# else {$stop = 1}
        # $k++;
    # }
	# return($add_SNPs);
# }

sub evaluate_region {
    my ($region_ref, $contrib) = @_;
    #print "in eval @$region_ref\n";
    my $count = 1;
	my $index;
    $contrib ||= CHILD;
    my @region_array = @$region_ref;
    my @array = ();
    if (!@$region_ref) {return()}
    my $local_lrr;
    if ($contrib == FATHER) {
        if ($REFINING) {$local_lrr = $LOCAL_P2_LRR_MED}
        else {$local_lrr = $P2_LRR_MED}
        $index = 4;
    }
    elsif ($contrib == MOTHER) {
        if ($REFINING) {$local_lrr = $LOCAL_P1_LRR_MED}
        else {$local_lrr = $P1_LRR_MED}
        $index = 3;
    }
    else {
        if ($REFINING) {$local_lrr = $LOCAL_CH_LRR_MED}
        else {$local_lrr = $CH_LRR_MED}
        $index = 5;
    }
	#print "lrrs meds $P1_LRR_MED $LOCAL_P1_LRR_MED $P2_LRR_MED $LOCAL_P2_LRR_MED $CH_LRR_MED $LOCAL_CH_LRR_MED\n";
    $region_array[5] = count_SNPs(@region_array[START,STOP], \%SNP_BY_POS);
    my ($median, $start_index, $stop_index, $lowest_index, 
        $highest_index) = @{calculate_region_median($region_ref, LRR, $index)};
#		print  "@region_array[0..2] median= $median\n";
#				print OUTPUT_FILE "@region_array[0..2] median= $median\n";
	#print "in eval $median ", $DEL_UPPER + $local_lrr,  " @region_array\n";
	#print "local lrr $local_lrr\n";
    #if ($median < ($DEL_UPPER / 2) + $local_lrr) {
	if ($median < $DEL_UPPER + $local_lrr) {
#print "in del\n";
		# If median is greater than $DEL_UPPER threshold, evaluate at threshold
		#$median = $DEL_UPPER if ($median > $DEL_UPPER);
        $region_array[3] = "DEL" if ($region_array[3] ne "HD");
        @region_array = @{cusum(\@region_array, \@LRR_REFS, $index, $median,
            $start_index, $stop_index, $lowest_index)};
    }
	#elsif ($median >= ($AMP / 2) + $local_lrr) {
    elsif ($median >= $AMP + $local_lrr) {
	#print "in amp\n";
		#$median = $AMP if ($median < $AMP);
        $region_array[3] = "AMP";
        @region_array = @{cusum(\@region_array, \@LRR_REFS, $index, $median,
            $start_index, $stop_index, $highest_index)};
    }
    else {
		#print "in eval else\n";

        $index = 10;
        $region_array[3] = " ";
		my $median_index;
        ($median, $start_index, $stop_index, $lowest_index, 
            $highest_index) = @{calculate_region_median($region_ref, BAF, $index)};
          #      print "@region_array[0..2] median= $median\n";

		#print OUTPUT_FILE "@region_array[0..2] median= $median\n";
		if ($median == NA) {return($region_ref)}
		elsif ($median > 0.1 + $LOCAL_CH_mBAF_MED && $median < $BB_BOUND) {
			#$median = 0.6 if ($median < 0.6);
			my (@inf, $median_index);
			#if ($region_array[4] == FATHER) {
			#	if ($region_array[17] ==  $hash_ref = \%P1_INF_SNP_BY_POS}
			#elsif ($region_array[4] == MOTHER) {$hash_ref = \%P2_INF_SNP_BY_POS}
			my $hash_ref = $region_array[16];
			#else {return(\@region_array)}
			foreach my $key (keys %$hash_ref) {
				push(@inf, $key) if ($key >= $region_array[1] &&
						$key <= $region_array[2]);
			}
			my @sorted_inf = sort { $a <=> $b } @inf;
			my $median_inf = $sorted_inf[int(@sorted_inf / 2)];
			#print "median inf $median_inf\n";
			if ($median_inf) {
				for (my $i = $start_index; $i <= $stop_index; $i++) {
					$median_index = $i if ($CURR_CHR_REFS[$i]->[POS] == $median_inf)
				}

				@region_array = @{cusum(\@region_array, \@CURR_CHR_REFS, $index, $median,
					$start_index, $stop_index, $median_index)};
			}
		}
	    else {$region_array[3] = " "}	
    }
    if ($contrib != CHILD) {$region_array[3] = " "}
    $region_array[5] = count_SNPs(@region_array[START,STOP], \%SNP_BY_POS);
    $region_array[6] = count_SNPs(@region_array[START,STOP,16]);
		#print "leaving eval $median @region_array\n";

    return(\@region_array);  
}   
  
sub extend_hd {
    my $indices_refs = $_[0];
    my $tree = Tree::Interval->new();
    # Extension array - first index = 0 or 1 has been extended, 
    # second index = fallback position, third index = count
    my @regions;
    foreach my $index_ref (@$indices_refs) {
        my ($upper_end, $lower_end) = (0) x 2;
        my ($up_index, $low_index) =  (1) x 2;
        my @region_array = ();

        my $line_ref = $LRR_REFS[$$index_ref[0]];
        @region_array[0..3] = (@$line_ref[CHR,POS], $$line_ref[POS], 
            $$index_ref[1]);
        if (@regions) {next if ($tree->find($$line_ref[POS]))}

        # Extend parental hom del regions
        my $id; 
        if    ($$index_ref[1] == FATHER) {$id = FATHER}
        elsif ($$index_ref[1] == MOTHER) {$id = MOTHER}
        elsif ($$index_ref[1] == NONE)   {$id = NONE}
               
        until ($upper_end) {
            if ($LRR_REFS[$$index_ref[0] + $up_index]) {
                my $new_line_ref = $LRR_REFS[$$index_ref[0] + $up_index];

                if (!$tree->find($$new_line_ref[POS])) {
                    my $LRR;
                    if    ($id == FATHER) {$LRR = $$new_line_ref[4]}
                    elsif ($id == MOTHER) {$LRR = $$new_line_ref[3]}
                    elsif ($id == NONE)   {$LRR = $$new_line_ref[5]}
                    if ($LRR < $HD_HI) {
                        $region_array[STOP] = $$new_line_ref[POS];
                    }
                    else {$upper_end = 1}
                }
                else {$upper_end = 1}
            }
            else {$upper_end = 1}
            $up_index++;
        }
   
        until ($lower_end) {
            if ($LRR_REFS[$$index_ref[0] - $low_index]) {
                my $new_line_ref = $LRR_REFS[$$index_ref[0] - $low_index];

                if (!$tree->find($$new_line_ref[POS])) {                    
                    my $LRR;
                    if    ($id == FATHER) {$LRR = $$new_line_ref[4]}
                    elsif ($id == MOTHER) {$LRR = $$new_line_ref[3]}
                    elsif ($id == NONE)   {$LRR = $$new_line_ref[5]}
    
                    if ($LRR < $HD_HI) {
                        $region_array[START] = $$new_line_ref[POS];
                    }
                    else {$lower_end = 1}
                }
                else {$lower_end = 1}
            }
            else {$lower_end = 1}
            $low_index++;
        }
        if ($$index_ref[1] == NONE) {$region_array[3] = "HD"}
        else {$region_array[3] = " "}
        my $hd_thresh;
        if ($id == FATHER)    {
            $region_array[4] = FATHER;
            $region_array[16] = \%P1_HD_BY_POS;
            $hd_thresh = $MIN_P2_HD;
        }
        elsif ($id == MOTHER) {            
            $region_array[4] = MOTHER;
            $region_array[16] = \%P2_HD_BY_POS;
            $hd_thresh = $MIN_P1_HD;
        }
        elsif ($id == NONE)   {            
            $region_array[4] = NONE;
            $region_array[16] = \%CH_HD_BY_POS;
            $hd_thresh = $MIN_CH_HD;
        }
        $region_array[6] = count_SNPs(@region_array[START,STOP,16]);
        $region_array[17] = HD;

        if ($region_array[6] >= $hd_thresh) {
            for (4..15,18) {$region_array[$_] ||= 0}
            push(@regions, \@region_array);

            $tree->insert(@region_array[START,STOP,START]);
        }
    }
    return(\@regions);
}                            

sub find_large_mosaic_POD_regions {
#print "find_large_mosaic_POD_regions\n";

	my $array_ref = $_[0];
	my @large_refs;
	#map{print "@{$_}\n"} @large_refs if @large_refs;
	my ($chr, $start, $stop, $p1, $p2) = (0) x 5;
	foreach my $ref (@$array_ref) {
		if (!$chr) {($chr, $start, $stop, $p1, $p2) = @$ref[0..4]}
		elsif ($chr == $$ref[0] && $stop >= $$ref[1]) {
			$stop = $$ref[2];
			$p1  += $$ref[3];
			$p2  += $$ref[4];
		}
		else {
			push(@large_refs, [$chr, $start, $stop, $p1, $p2]);
			($chr, $start, $stop, $p1, $p2) = @$ref[0..4];
		}
	}
	#map{print "@{$_}\n"} @large_refs if @large_refs;

	foreach my $ref (@large_refs) {
		$$ref[$_] = 0 for (5..20);
		if ($$ref[3] > $$ref[4]) {
			@$ref[3..6] = (" ", FATHER, 500, $$ref[3]);
			@$ref[16,17] = (\%P1_INF_SNP_BY_POS, POD);
		}
		else {
			@$ref[3..6] = (" ", MOTHER, 500, $$ref[4]);
			@$ref[16,17] = (\%P2_INF_SNP_BY_POS, POD);
		}
	}
	#print STDERR "indefinePODcr\n";
	#map{print STDERR "@{$_}\n"} @large_refs if @large_refs;
	#print STDERR "podcr to collapse\n";
	collapse_region_list(\@large_refs);
		#print STDERR "back from collapse\n";

	for my $ref (@large_refs) {
					#print STDERR "definePODcr to find\n";
					#print STDERR "@$ref[START,STOP]\n";
		my $array_ref;
		if ($$ref[4] == FATHER) {$array_ref = \@P1_INF_SNPS}
		else {$array_ref = \@P2_INF_SNPS}
		
		my ($closest_start_ind, $closest_stop_ind) 
			= find_closest_indices(@$ref[START,STOP,16]);
			#print "@$ref[START,STOP]\n";
			#print "closest $closest_start_ind $closest_stop_ind\n";
			#print "$$array_ref[$closest_start_ind]->[POS] $$array_ref[$closest_stop_ind]->[POS]\n";
		$$ref[START] = $$array_ref[$closest_start_ind]->[POS];
		$$ref[STOP]  = $$array_ref[$closest_stop_ind]->[POS];
		$$ref[6] = count_SNPs(@$ref[START,STOP,16]);
	}

	#print STDERR "back from ref\n";
	#map{print STDERR "@{$_}\n"} @large_refs if @large_refs;
	#map{print "@{$_}\n"} @large_refs if @large_refs;
	
    #print "sending to eval\n";
    map {$_ = evaluate_region($_)} @large_refs;
		#print STDERR "cryptic after eval\n";
	#map{print STDERR "@{$_}\n"} @large_refs if @large_refs;
		#print STDERR "podcr to collapse2\n";

    collapse_region_list(\@large_refs);
		#print STDERR "cryptic after collapse\n";
	#map{print STDERR "@{$_}\n"} @large_refs if @large_refs;

    map {$_ = evaluate_region($_) unless $$_[6]} @large_refs;
	#my @return;
	#map {push(@return, $_) if ($$_[6] >= $MIN_POD)} @large_refs;
    #map {calculate_region_stats($_)} @return;
    #return (\@return);

	map {calculate_region_stats($_)} @large_refs;
	#map{print "@{$_}\n"} @large_refs if @large_refs;
	#print "back from large POD\n";
    return (\@large_refs);
}

sub find_POD_regions {
    # Calculates statistical significance for windows with adequate 
    # informative SNPs and combines overlapping abnormal windows.
    # Calls define_POD_region and returns a reference to the refined regions.
	my (@abnormal_regions, @abn_region, @temp_region_refs, @return_regions, 
		@cryptic_regions, @large_windows, @temp_large_win);
    my ($flag1, $flag2, $in_region1, $in_region2, $last_sig_SNP, 
        $p1_sig, $p2_sig, $chr, $end_of_region, $abn_region, $region_size,
		$last_ext) = (0) x 12;

    my $size = scalar(@CURR_CHR_REFS);
    if($size < $SIZE_REGION) {$region_size = $size}
    else {$region_size = $SIZE_REGION}
	my $small_window = 1;
	my $array_ref = overlapping_windows(\@CURR_CHR_REFS, $region_size, 
		$small_window);
		
    my $large_window_size = $SIZE_REGION * 5;
	$small_window = 0;

    for (my $i = 0; $i < @$array_ref; $i++) {
        my $item_ref = $$array_ref[$i];
        my ($p1_ct, $p2_ct) = @$item_ref[3,4];
        my $not_in_region = (!$in_region1 && !$in_region2);
        $chr = $$item_ref[0];  
        # If the last window overlaps a previous abnormal region, it may 
        # contain another small nonoverlapping region. In order to detect
        # the additional region, the window boundaries are adjusted to be
        # nonoverlapping and the informative SNPs are recalculated.
        if ($$item_ref[START] <= $end_of_region && $not_in_region) { 
            if ($i == @$array_ref - 1) { 
                $$item_ref[START] 
                    = $SNP_BY_NUM{$SNP_BY_POS{$end_of_region} + 1};
                $p1_ct = count_SNPs(@$item_ref[START,STOP],
					\%P1_INF_SNP_BY_POS);
                $p2_ct = count_SNPs(@$item_ref[START,STOP], 
					\%P2_INF_SNP_BY_POS);
            }
            else {
                $p1_sig = $p2_sig = 0;
                next; 
            }
        }
        # Ignore windows with number of informative SNPs less than
        # the previously calculated minimum number of SNPs for 
        # which the probability can possibly be below the significance
        # threshold. 
		my $abnormal = my $next = 0;
        if ($p1_ct + $p2_ct < $MIN_POD) { 
            if ($not_in_region) { 
                $p1_sig = $p2_sig = 0;
				push(@temp_large_win, $item_ref);
                next; 
            }
        }
        else {
            # Calculate probability of parental contributions in a 
            # window (normal or abnormal).
            # P(k) = (n choose k) * p^n, where
            # n = number of informative SNPs per window
            # k = number of paternal contributions per window
            # p = probability of a paternal contribution for an informative SNP
			my $prob = 0;
			if ($p1_ct == $p2_ct) {$prob = 1}
			else {
				my ($n, $k, $p) = ($p1_ct + $p2_ct, $p1_ct, 0.5);
				my ($x, $y);
				if ($k <= $n / 2) {($x, $y) = (0, $k)}
				else {($x, $y) = ($k, $n)}
				for ($x..$y) {
					my $bin_coeff = calculate_bin_coeff($n, $_);
					$prob += $bin_coeff * ($p**$_) * ((1 - $p)**($n - $_));
					#print $bin_coeff * ($p**$_) * ((1 - $p)**($n - $_)), "\n";
					#print "$prob $n $_\n";
				}
				$prob = 2 * $prob;
			}
			#print OUTPUT_FILE "$prob $POD_ALPHA $p1_ct $p2_ct\n";
            # If the probability is significantly abnormal, determine
            # which parent was the major contributor 
            if ($prob <= $POD_ALPHA) {
				$abnormal = 1;
				#print OUTPUT_FILE "abn @$item_ref\n";
                if ($p1_ct > $p2_ct) {($p1_sig, $p2_sig) = (1,0)} 
                else {($p1_sig, $p2_sig) = (0,1)}
				
				if (scalar(@temp_large_win) >= $large_window_size) {
					my ($start, $stop) = find_closest_indices($temp_large_win[0]->[1], 
					$temp_large_win[-1]->[2], \%BAF_BY_POS);
					#print [@CURR_CHR_REFS[$start..$stop]], "\n";
					#print OUTPUT_FILE "$temp_large_win[0]->[0] $temp_large_win[0]->[1] $temp_large_win[-1]->[2]\n";
					push(@large_windows, @{overlapping_windows([@CURR_CHR_REFS[$start..$stop]], 
						$large_window_size, $small_window)});
				}
				@temp_large_win = ();
            }
            elsif ($not_in_region) { 
                $p1_sig = $p2_sig = 0;
				push(@temp_large_win, $item_ref);				
                next; 
            }
            else {$p1_sig = $p2_sig = 0}
        }
 		
		if ($in_region1) {
			#print OUTPUT_FILE "$$item_ref[START] $last_sig_SNP region1 $in_region1 $p1_ct $p2_ct $p1_sig $p2_sig $POD_ALPHA\n";

			if (!$p1_sig && $p2_ct > $ACCEPTABLE_ERRORS) {$flag2 = 1}
			elsif ($$item_ref[START] > $last_ext) {
				$BOUNDARY_EXTENSION++;
				$last_ext = $$item_ref[STOP];
			}
		}
        elsif ($in_region2) {
			if (!$p2_sig && $p1_ct > $ACCEPTABLE_ERRORS) {$flag1 = 1}
			elsif ($$item_ref[START] > $last_ext) {
				$BOUNDARY_EXTENSION++;
				$last_ext = $$item_ref[STOP];
			}
		}
			
        # Define abnormal regions
        if ($not_in_region) {
            if ($p1_sig) {
                push(@temp_region_refs, $item_ref);
                ($last_sig_SNP, $in_region1) = (0,1);
                next;
            }
            elsif ($p2_sig) {
                push(@temp_region_refs, $item_ref);
                ($last_sig_SNP, $in_region2) = (0,1);
                next;
            }
        }
        elsif ($in_region1) {
            #Store current line or finalize calculations for region
            if (!$flag2) {
                push(@temp_region_refs, $item_ref);
                $last_sig_SNP = $#temp_region_refs if ($p1_sig);           
            }
            else {
                $abn_region = define_POD_region(\@temp_region_refs, $in_region1,
                    $last_sig_SNP); 
                push(@abnormal_regions, $abn_region);
                $end_of_region = $$abn_region[2];
                @temp_region_refs = ();
                $in_region1 = $flag2 = $last_sig_SNP = 0;            
                if ($p2_sig) { 
                    push(@temp_region_refs, $item_ref);
                ($last_sig_SNP, $in_region2) = (0,1);
                }
            }      
        }
        elsif ($in_region2) {
            #Store current line or finalize calculations for region
            if (!$flag1) {
                push(@temp_region_refs, $item_ref);            
                $last_sig_SNP = $#temp_region_refs if ($p2_sig);
            }
            else { 
                unless ($chr == CHR_X && $GENDER == MALE) {
                    $abn_region = define_POD_region(\@temp_region_refs, 
                        $in_region1, $last_sig_SNP); 
                    push(@abnormal_regions, $abn_region);
                }
                $end_of_region = $$abn_region[2];                
                @temp_region_refs = ();
                $in_region2 = $flag1 = $last_sig_SNP = 0;
                
                if ($p1_sig) { 
                    push(@temp_region_refs, $item_ref);
                ($last_sig_SNP, $in_region1) = (0,1);
                }
            }
        }   
        $p1_sig = $p2_sig = 0;    
    }

    if ($temp_region_refs[-1]) { 
        unless ($in_region2 && $chr == CHR_X && $GENDER == MALE) {
            $abn_region = define_POD_region(\@temp_region_refs, $in_region1,
                $last_sig_SNP); 
            push(@abnormal_regions, $abn_region);
        }
    }
	
	if (scalar(@temp_large_win) >= $large_window_size) {
		my ($start, $stop) = find_closest_indices($temp_large_win[0]->[1], 
		$temp_large_win[-1]->[2], \%BAF_BY_POS);
		#print [@CURR_CHR_REFS[$start..$stop]], "\n";
		#print OUTPUT_FILE "$temp_large_win[0]->[0] $temp_large_win[0]->[1] $temp_large_win[-1]->[2]\n";
		push(@large_windows, @{overlapping_windows([@CURR_CHR_REFS[$start..$stop]], 
			$large_window_size, $small_window)});
	}
	@temp_large_win = ();

    if (@abnormal_regions) {
		#map{print "in abn @{$_}\n"} @abnormal_regions;
			#print STDERR "find_abn tocollapse\n";
		
        collapse_region_list(\@abnormal_regions);
		#map{print "in abn after collapse @{$_}\n"} @abnormal_regions;

        map {$_ =  evaluate_region($_)} @abnormal_regions;
		
        map {calculate_region_stats($_)} @abnormal_regions;
		#map{print "in abn after eval @{$_}\n"} @abnormal_regions;

    }
    #map{print OUTPUT_FILE "cryptic @{$_}\n"} @cryptic_regions if @cryptic_regions;
	
	return(\@abnormal_regions, \@large_windows);    
}

sub find_closest_indices {
    my ($start, $stop, $hash_ref) = @_;
	return unless ($start && $stop && $hash_ref);
    my ($closest_lower, $closest_higher, $value1, $value2);
	my $present = 0;
    my (@indices, @sorted_indices);
    $value1 = $$hash_ref{$start};
    $value2 = $$hash_ref{$stop};   
    if (!$value1 || !$value2) {
	#print "in find closest neither $start $stop\n";
        foreach my $key (keys %$hash_ref) {
            push(@indices, $$hash_ref{$key}) if ($key >= $start && $key <= $stop);
        }
		#print "indices @indices\n";
        if (@indices) {
            @sorted_indices = sort { $a <=> $b } @indices;
			#print "sortedindices @indices\n";
		
            $closest_lower  = $sorted_indices[0];
            $closest_higher = $sorted_indices[-1];
            $present = 1;
        }
    }
    else {
        $closest_lower  = $value1;
        $closest_higher = $value2;
        $present = 1;
    }
    return($closest_lower, $closest_higher, $present);
}

sub find_informative_snps {
    # Analyzes each SNP for parental contribution. 
    # Results are added to the input array in binary form.
    # Returns a number of variables for genomewide calculations.

    my $array_ref   = $_[0];
    my ($p1, $p2, $p1_mi1, $p2_mi1, $un_mi1, $outlier, $ch_upper_MI1, 
		$ch_lower_MI1, $outlier_ct);
    if ($$array_ref[10] > $MI1_UPPER_THRESH) {$ch_upper_MI1 = 1} 
    elsif ($$array_ref[10] < $MI1_LOWER_THRESH) {$ch_lower_MI1 = 1} 
    
	if ($$array_ref[3] eq "BB") { 
	    if ($$array_ref[6] eq "BB" && $$array_ref[9] eq "AB") {$un_mi1 = 1}
        elsif ($$array_ref[6] eq "AA") {
            if ($$array_ref[10] > $AB_UPPER_BOUND) {
                $p1 = 1;
                $p1_mi1 = 1 if ($ch_upper_MI1);
            }
            elsif ($$array_ref[10] < $AB_LOWER_BOUND) {
                $p2 = 1;
                $p2_mi1 = 1 if ($ch_lower_MI1);
            }
        }
        elsif ($$array_ref[6] eq "AB") {
            if ($$array_ref[10] < $AB_LOWER_BOUND) {
                $p2 = 1;
                $p2_mi1 = 1 if ($ch_lower_MI1);
            }
        }
        elsif (substr($$array_ref[6], 0, 1) eq "N") {
            if (($$array_ref[10] < $AB_LOWER_BOUND || $$array_ref[9] eq "AA") 
				&& $$array_ref[8] > $HD_MID) {$p2 = 1}
        }
    }
    elsif ($$array_ref[3] eq "AA") {    
        if ($$array_ref[6] eq "BB") {
		    if ($$array_ref[10] > $AB_UPPER_BOUND) {
                $p2 = 1;
                $p2_mi1 = 1 if ($ch_upper_MI1);
            }
			elsif ($$array_ref[10] < $AB_LOWER_BOUND) {
                $p1 = 1;
                $p1_mi1 = 1 if ($ch_lower_MI1);
            }
        }
        elsif ($$array_ref[6] eq "AB") {
            if ($$array_ref[10] > $AB_UPPER_BOUND) {
                $p2 = 1;
                $p2_mi1 = 1 if ($ch_upper_MI1);
            }
			elsif ($$array_ref[10] < $AB_LOWER_BOUND 
				&& $$array_ref[9] eq "AB") {$outlier = 1}
        }
        elsif (substr($$array_ref[6], 0, 1) eq "N") {
            if (($$array_ref[10] > $AB_UPPER_BOUND || $$array_ref[9] eq "BB") 
				&& $$array_ref[8] > $HD_MID) {$p2 = 1}
        }
        elsif ($$array_ref[6] eq "AA" && $$array_ref[9] eq "AB") {$un_mi1 = 1}
    }
    elsif ($$array_ref[3] eq "AB") {
        if ($$array_ref[6] eq "BB") {
            if ($$array_ref[10] < $AB_LOWER_BOUND) {
                $p1 = 1;
                $p1_mi1 = 1 if ($ch_lower_MI1);
            }
        }
        elsif ($$array_ref[6] eq "AA") {
            if ($$array_ref[10] > $AB_UPPER_BOUND) {
                $p1 = 1;
                $p1_mi1 = 1 if ($ch_upper_MI1);
            }
        }
    }
    elsif (substr($$array_ref[3], 0, 1) eq "N") {
        if ($$array_ref[6] eq "AA") {
            if (($$array_ref[10] > $AB_UPPER_BOUND || $$array_ref[9] eq "BB") 
				&& $$array_ref[5] > $HD_MID) {$p1 = 1}
        }                            
        elsif ($$array_ref[6] eq "BB") {
            if (($$array_ref[10] < $AB_LOWER_BOUND || $$array_ref[9] eq "AA") 
				&& $$array_ref[5] > $HD_MID) {$p1 = 1}
        }                    
    }
        	
	# if ($$array_ref[9] ne "AA" && $$array_ref[9] ne "BB" 
		# && (($$array_ref[10] > $AB_UPPER_BOUND && $$array_ref[10] < $BB_BOUND) 
		# || $$array_ref[10] < $AB_LOWER_BOUND && $$array_ref[10] > $AA_BOUND)) {
			# $outlier_ct = 1;
	# }
	if ($$array_ref[9] eq "AB" && ($$array_ref[10] > $PODCR_UPPER_THRESH
		|| $$array_ref[10] < $PODCR_LOWER_THRESH)) {$outlier_ct = 1}
			
	my ($inf_snp, $mi1) = (0) x 2;
    if ($$array_ref[5] < $HD_MID || $$array_ref[8] < $HD_MID 
		|| $$array_ref[11] < $HD_MID) {push(@$array_ref, UNKNOWN)}   
    elsif ($p1) {
        push(@$array_ref, FATHER);
        $inf_snp = 1 unless ($$array_ref[CHR] == CHR_X);
    }
    elsif ($p2) {
        push(@$array_ref, MOTHER);
        $inf_snp = 1 unless ($$array_ref[CHR] == CHR_X);
    }
	elsif ($un_mi1) {push(@$array_ref, UNKNOWN_MI1)}
	else {push(@$array_ref, UNKNOWN)}
    my $type = 0;
    if ($p1_mi1) {
        $mi1 = $p1_mi1;
        $type = FATHER;
    }
    elsif ($p2_mi1) {
        $mi1 = $p2_mi1;
        $type = MOTHER;
    }
    if ($un_mi1) {
        $mi1 = $un_mi1;
        $type = NONE;
    }
    return([$p1, $p2, $inf_snp, $mi1, $type, $outlier_ct]);
}

sub get_R_PID {
    $R_PID = <$R_PID_FILE>;
	print $R_PID, "\n";
    if (!$R_PID) {
        print "An error has occurred in communication with the R ", 
        "software.\n"
    }
    return($R_PID);
}

sub help {
    print <<END;

Usage: perl $0 [--options] [INPUT_FILE]
  --alpha   A threshold for which a p-value <= alpha is considered significant.
            Default = --alpha=0.05
  --batch   Submit a file containing a list of file names to be run sequentially
  --bin     Number of SNPs per bin for sliding window analysis
            Default = --bin=100
  --build   UCSC Genome assembly for centro locations (e.g. --build=hg19)
            Default = --build=hg18
  --cores   Number of CPU cores to employ
            Default = maximum cores - 1 (e.g. --cores=8)
  --gender  Gender designation for sample (M or F)
            Default = NA
  --graph   Creates graphic output of results (--graph or --nograph)
            Default = --graph 
  --hd		Abnormality detection by homozygous deletion analysis (PODhd)
			(--hd or --nohd) Default = --hd 
  --help    Prints a help message and exits
  --hetSD   Heterozygous SNP threshold. StDevs from mean BAF
            Default = --hetSD=1.414213562373095 (sqrt of 2)
  --homSD   Homozygous SNP threshold. StDevs from mean BAF
            Default = --homSD=4
  --mi1		Abnormality detection by Mendelian "error" analysis (PODmi1)
			(--mi1 or --nomi1) Default = --mi1 
  --out     Specify an output directory (e.g. --out=results)
            Default = --out=triPOD_Results  
  --pod		Abnormality detection by standard Parent-of-Origin-based Detection
			alogrithm 
			(--pod or --nopod) Default = --pod	
  --podcr	Abnormality detection by cryptic Parent-of-Origin-based Detection
			algorithm (PODcr) 
			(--podcr or --nopodcr) Default = --podcr				
  --verbose Prints progress info to screen. Negated using --noverbose.
			--batch mode will run in --noverbose mode.
			Default = --verbose
  
END
    print_proper_format();
    exit;
}

sub initial_calculations {
	my (@ch_stats, @p1_stats, @p2_stats);
	@p1_stats = @p2_stats = map {$ch_stats[$_] = 0} (0..10);
	# The following are *_stats array indices
	my ($CHR, $START, $BAF_REF, $MBAF_MED, $LRR_MED, $NC_CT,
		$BAF_SNP_CT, $LRR_SNP_CT, $NC_SNP_CT, $STDEV) = (0..9);
    my ($ch_baf_clust_sum, $ch_baf_clust_sum_sq, $ch_baf_clust_ct, 
		$p1_baf_clust_sum, $p1_baf_clust_sum_sq, $p1_baf_clust_ct,
		$p2_baf_clust_sum, $p2_baf_clust_sum_sq, $p2_baf_clust_ct);		
	my (@ch_baf, @ch_mbaf, @ch_lrr, @mi1_count, @p1_baf, @p1_mbaf, @p1_lrr,
		@p2_baf, @p2_mbaf, @p2_lrr);
	for (0..10) {$ch_baf[$_] = 0}
	@p1_baf = @p2_baf = @ch_baf;

    foreach my $line_ref (@INIT_CALC_REFS) {
		my @line_array = @$line_ref;
		next if (!$line_array[3] || !$line_array[6] || !$line_array[9]);
        $p1_stats[$NC_CT]++ if substr($line_array[3], 0, 1) eq "N";
        $p2_stats[$NC_CT]++ if substr($line_array[6], 0, 1) eq "N";
        $ch_stats[$NC_CT]++ if substr($line_array[9], 0, 1) eq "N";
        $p1_stats[$NC_SNP_CT]++;
        next if (substr($line_array[4], 0, 1) eq "N" 
			||  substr($line_array[7],  0, 1) eq "N" 
			||  substr($line_array[10], 0, 1) eq "N");
        push(@p1_lrr, $line_array[5]);
	    push(@p2_lrr, $line_array[8]);
        push(@ch_lrr, $line_array[11]);
		$p1_stats[$LRR_SNP_CT]++;
		
        if (substr($line_array[0], 0, 2) =~ /cn/i) {
            $p1_stats[$NC_CT]-- if substr($line_array[3], 0, 1) eq "N";
            $p2_stats[$NC_CT]-- if substr($line_array[6], 0, 1) eq "N";
            $ch_stats[$NC_CT]-- if substr($line_array[9], 0, 1) eq "N";
            next;
        }
        $p1_stats[$BAF_SNP_CT]++;

        # Store values for BAF StDev calculations
        baf_stats_for_stdev($line_array[4],  0.25, 0.75, \@p1_baf, \@p1_mbaf); 
		baf_stats_for_stdev($line_array[7],  0.25, 0.75, \@p2_baf, \@p2_mbaf); 
		baf_stats_for_stdev($line_array[10], 0.25, 0.75, \@ch_baf, \@ch_mbaf); 
		 		
		# Store data for stdev clustering
		if ($line_array[10] > 0.02 && $line_array[10] < 0.98) {
			$ch_baf_clust_sum += $line_array[10];
			$ch_baf_clust_sum_sq += $line_array[10]**2;
			$ch_baf_clust_ct++;
		}
		if ($line_array[4] > 0.02 && $line_array[4] < 0.98) {
			$p1_baf_clust_sum += $line_array[4];
			$p1_baf_clust_sum_sq += $line_array[4]**2;
			$p1_baf_clust_ct++;
		}
		if ($line_array[7] > 0.02 && $line_array[7] < 0.98) {
			$p2_baf_clust_sum += $line_array[7];
			$p2_baf_clust_sum_sq += $line_array[7]**2;
			$p2_baf_clust_ct++;
		}		
    }
	return() unless ($p1_stats[$BAF_SNP_CT] >= $SIZE_REGION);
	
	(@ch_stats[$CHR,$START], @p2_stats[$CHR,$START], @p1_stats[$CHR,$START]) 
		= (@{$INIT_CALC_REFS[0]}[1,2]) x 3;
	$ch_stats[$BAF_SNP_CT] = $p2_stats[$BAF_SNP_CT] = $p1_stats[$BAF_SNP_CT];
	$ch_stats[$LRR_SNP_CT] = $p2_stats[$LRR_SNP_CT] = $p1_stats[$LRR_SNP_CT];
	$ch_stats[$NC_SNP_CT]  = $p2_stats[$NC_SNP_CT]  = $p1_stats[$NC_SNP_CT];
	$ch_stats[$BAF_REF] = \@ch_baf;
	$p1_stats[$BAF_REF] = \@p1_baf;
	$p2_stats[$BAF_REF] = \@p2_baf;
	
	my $p1_baf_mean = $p1_baf[9] / $p1_baf[10];
	my $p2_baf_mean = $p2_baf[9] / $p2_baf[10];
	my $ch_baf_mean = $ch_baf[9] / $ch_baf[10];
    map {$_ = 1 - $_ if ($_ < $p1_baf_mean)} @p1_mbaf;    
    map {$_ = 1 - $_ if ($_ < $p2_baf_mean)} @p2_mbaf; 
	map {$_ = 1 - $_ if ($_ < $ch_baf_mean)} @ch_mbaf;
	$p1_stats[$MBAF_MED] = median(\@p1_mbaf);
	$p2_stats[$MBAF_MED] = median(\@p2_mbaf);
	$ch_stats[$MBAF_MED] = median(\@ch_mbaf);
	$p1_stats[$LRR_MED]  = median(\@p1_lrr);
	$p2_stats[$LRR_MED]  = median(\@p2_lrr);
	$ch_stats[$LRR_MED]  = median(\@ch_lrr);	
	$ch_stats[$STDEV] = st_dev($ch_baf_clust_sum, $ch_baf_clust_sum_sq,
								$ch_baf_clust_ct);
	$p1_stats[$STDEV] = st_dev($p1_baf_clust_sum, $p1_baf_clust_sum_sq,
								$p1_baf_clust_ct);
	$p2_stats[$STDEV] = st_dev($p2_baf_clust_sum, $p2_baf_clust_sum_sq,
								$p2_baf_clust_ct);
    return([\@ch_stats, \@p1_stats, \@p2_stats]);
}

sub baf_stats_for_stdev {
	my ($baf, $lower, $upper, $baf_ref, $mbaf_ref) = @_;
	my ($aa_sum, $aa_sum_squares, $aa_count, $ab_sum, $ab_sum_squares,
		$ab_count, $bb_sum, $bb_sum_squares, $bb_count, $baf_sum, 
		$baf_ct);
		
	# $baf_ref array indices aa_sum, aa_sum_squares, aa_count, ab_sum, ab_sum_squares,
	#	 ab_count, bb_sum, bb_sum_squares, bb_count, baf sum for mbaf, baf ct for mbaf
	if ($baf < $upper) {
		if ($baf <= $lower) {
			# AA
			$$baf_ref[0] += $baf;
			$$baf_ref[1] += $baf**2;
			$$baf_ref[2]++;
		}
		else {
			# AB
			$$baf_ref[3] += $baf;
			$$baf_ref[4] += $baf**2;
			$$baf_ref[5]++;
			push(@$mbaf_ref, $baf);
			$$baf_ref[9] += $baf;
			$$baf_ref[10]++;
		}
	}
	else {
		# BB
		$$baf_ref[6] += $baf;
		$$baf_ref[7] += $baf**2;
		$$baf_ref[8]++;
	}
}


sub join_analysis_threads {
    my ($results_ref, $threads_ref, $thread_count) = @_;
    for (0..$thread_count - 1) { 
        if ($$threads_ref[$_]->is_joinable()) {
            my $results        = $$threads_ref[$_]->join();
			if($results) {
				$$results_ref[$_]  = $$results[0];
				$INF_SNPS         += $$results[1];
				$ADJUSTED_MI1     += $$results[2];
				if ($$results[3]) {
					if ($LOCAL_STATS[$$results[3]]) {
						$LOCAL_STATS[$$results[3] + 30] = $$results[4];
					}                
					else {$LOCAL_STATS[$$results[3]] = $$results[4]}
				}
			}
			#else {print "no results\n"}
        }               
    }
}

sub join_init_threads {
    my ($threads_ref, $thread_count) = @_;
    for (0..$thread_count - 1) {
        if ($$threads_ref[$_]->is_joinable()) {
			my $results = $$threads_ref[$_]->join();
			if ($results) {
				push(@INIT_CH_STATS, $$results[0]);
				push(@INIT_P1_STATS, $$results[1]);
				push(@INIT_P2_STATS, $$results[2]);
			}
        }               
    }
}

sub join_stdprob_threads {
    my ($threads_ref, $thread_count) = @_;
    for (0..$thread_count - 1) { 
        if ($$threads_ref[$_]->is_joinable()) {
            my $results = $$threads_ref[$_]->join();
			if ($results) {
				push(@NONOVERLAP_RESULTS, $results);
			}
        }               
    }
}

sub manage_chromosome_arm {
    # This is the hub for the analysis. Each thread uses this function to 
    # start a cascade of function calls and to return the detected regions to 
    # the main program. 
    # Descriptive info is collected for later genomewide calculations.
    # Returns a reference to the final detected regions and a number of 
    # descriptive variables.
	($LOCAL_CH_BAF_MEAN, $LOCAL_CH_mBAF_MED, $LOCAL_P1_mBAF_MED, 
	 $LOCAL_P2_mBAF_MED, $LOCAL_CH_LRR_MED, $LOCAL_P1_LRR_MED, 
	 $LOCAL_P2_LRR_MED, $BOUNDARY_EXTENSION) = (0) x 8;
    (@LRR_REFS, @BAF_OUTLIERS, @P1_INF_SNPS, @P2_INF_SNPS) = (()) x 4;
    (%BAF_BY_POS, %BAF_OUTLIERS_BY_POS, %CH_HD_BY_POS, %LRR_BY_POS, %MI1_BY_POS, 
	 %P1_HD_BY_POS, %P1_INF_SNP_BY_POS, %P2_HD_BY_POS, %P2_INF_SNP_BY_POS, 
	 %SNP_BY_NUM, %SNP_BY_POS) = (()) x 11;
    my (@ch_mbaf_array, @ch_lrr_array, @p1_mbaf_array, @p1_lrr_array, 
		@p2_mbaf_array, @p2_lrr_array, @remove);
	my ($ch_baf_sum, $ch_baf_ct, $p1_baf_sum, $p1_baf_ct, $p2_baf_sum, 
	    $p2_baf_ct) = (0) x 7;

    $REFINING = 0;
    for (my $i = 0; $i < @CURR_CHR_REFS; $i++) {
		if (!$CURR_CHR_REFS[$i]->[3] || !$CURR_CHR_REFS[$i]->[6] 
			|| !$CURR_CHR_REFS[$i]->[9]) {
			push(@remove, $i);
			next;
		}
		# Create hashes for SNP number calculations  
        $SNP_BY_NUM{keys %SNP_BY_POS} = $CURR_CHR_REFS[$i]->[POS];        
        $SNP_BY_POS{$CURR_CHR_REFS[$i]->[POS]} = keys %SNP_BY_POS; 

        if (substr($CURR_CHR_REFS[$i]->[4], 0, 1) eq "N" 
            || substr($CURR_CHR_REFS[$i]->[7], 0, 1) eq "N"
            || substr($CURR_CHR_REFS[$i]->[10], 0, 1) eq "N") {
            push(@remove, $i);
            next;
        }
        # Create chromosome array including cnv markers for LRR calculations
        my $lrr_ref = [@{$CURR_CHR_REFS[$i]}[0..2,5,8,11]];
        push(@LRR_REFS, $lrr_ref);
        $LRR_BY_POS{$$lrr_ref[POS]} = $#LRR_REFS;    
                
        # Ignore intensity only markers for baf analyses
        if (substr($CURR_CHR_REFS[$i]->[0], 0, 2) =~ /cn/i) {
			push(@remove, $i);
		}
        else {
			my $ch_baf = $CURR_CHR_REFS[$i]->[10];
			if ($ch_baf > $AA_BOUND && $ch_baf < $BB_BOUND) {
				$ch_baf_sum += $ch_baf;
				$ch_baf_ct++;
			}
        }
    }
   
    delete @CURR_CHR_REFS[@remove] if @remove;
    my @temp;
    map {push(@temp, $_) if $_} @CURR_CHR_REFS;
    @CURR_CHR_REFS = @temp;
	return() unless ($#CURR_CHR_REFS >= $SIZE_REGION);
	$LOCAL_CH_BAF_MEAN = $ch_baf_sum / $ch_baf_ct;

    (@ch_mbaf_array, @ch_lrr_array, @p1_mbaf_array, @p1_lrr_array, 
	 @p1_mbaf_array, @p1_lrr_array, @remove) = (()) x 7;
	$ch_baf_sum = $ch_baf_ct = 0;

    ############################## Process chromosome #############################
    my @result = @{process_chromosome_arm()};
    my (@baf_normal, @lrr_normal) = (()) x 2;

    my @final_regions = @{$result[0]};
    my $total_snps = @CURR_CHR_REFS;
    my $abn_snps = 0;
    my (@baf_abn_indices, @lrr_abn_indices);

    if (@final_regions) {
        foreach my $ref (@final_regions) {$abn_snps += $$ref[5]}
        my (@baf_abn_indices, @lrr_abn_indices);
		#print "$abn_snps $total_snps\n";
        if ($abn_snps < $total_snps * 0.75) {
            foreach my $ref (@final_regions) {
							#print STDERR "manage to find\n";
#print STDERR "@$ref[START,STOP]\n";
                my ($baf_start_ind, $baf_stop_ind) 
                    = find_closest_indices(@$ref[START,STOP], \%BAF_BY_POS);
                if ($baf_start_ind && $baf_stop_ind) {
					push(@baf_abn_indices, $baf_start_ind..$baf_stop_ind);
				}
				my ($lrr_start_ind, $lrr_stop_ind) 
                    = find_closest_indices(@$ref[START,STOP], \%BAF_BY_POS);
				if ($lrr_start_ind && $lrr_stop_ind) {
					push(@lrr_abn_indices, $lrr_start_ind..$lrr_stop_ind);		
				}
            }
        }
        else {
            # If most of the chromosome arm is abnormal analyses are not 
            # repeated using local medians.
            # Indicate single pass through analyses and return results.
			#print "returning\n";
            map {$$_[19] = 1} @final_regions;
            return(\@result);
		}

        if (@baf_abn_indices) {
            my @baf_temp = @CURR_CHR_REFS;
            delete @baf_temp[@baf_abn_indices];
            map {push(@baf_normal, $_) if $_} @baf_temp;

            @baf_temp = (); 
        }
        else {@baf_normal = @CURR_CHR_REFS}
                    
        my @lrr_temp = @LRR_REFS;
        delete @lrr_temp[@lrr_abn_indices];
        map {push(@lrr_normal, $_) if $_} @lrr_temp;
        @lrr_temp = (); 
    }
    else {
        @baf_normal = @CURR_CHR_REFS;
        @lrr_normal = @LRR_REFS;
    }       
    my ($aa_sum, $aa_sum_squares, $aa_count, $ab_sum, $ab_sum_squares,
        $ab_count, $bb_sum, $bb_sum_squares, $bb_count);
    foreach my $ref (@baf_normal) {
        # Store values for BAF StDev calculations
        my $ch_baf = $$ref[10];   
        if ($ch_baf < 0.75) {
            if ($ch_baf <= 0.25) {
                $aa_sum += $ch_baf;
                $aa_sum_squares += $ch_baf**2;
                $aa_count++;
            }
            else {
                $ab_sum += $ch_baf;
                $ab_sum_squares += $ch_baf**2;
                $ab_count++;
            }
			if ($ch_baf > $AA_BOUND) {
				push(@ch_mbaf_array, $ch_baf);
				$ch_baf_sum += $ch_baf;
				$ch_baf_ct++;
			}
        }
        else {
            $bb_sum += $ch_baf;
            $bb_sum_squares += $ch_baf**2;
            $bb_count++;
			if ($ch_baf < $BB_BOUND) {
				push(@ch_mbaf_array, $ch_baf) ;
				$ch_baf_sum += $ch_baf;
				$ch_baf_ct++;
			}				
        }
        my $baf = $$ref[4];
        if ($baf > $P1_AA_BOUND && $baf < $P1_BB_BOUND) {
            push(@p1_mbaf_array, $baf); 
			$p1_baf_sum += $baf;
			$p1_baf_ct++;			
        }
        $baf = $$ref[7];
        if ($baf > $P2_AA_BOUND && $baf < $P2_BB_BOUND) {
            push(@p2_mbaf_array, $baf);                
			$p2_baf_sum += $baf;
			$p2_baf_ct++;			
        }
    }
    
    #Store initial thresholds
    my @temp_bounds = ($AA_BOUND, $BB_BOUND, $AB_UPPER_BOUND, $AB_LOWER_BOUND,
		$MI1_UPPER_THRESH, $MI1_LOWER_THRESH, $ACCEPTABLE_ERRORS);
    # Calculate standard deviation for BAFs and BAF thresholds for analysis
    my @AA_stats = st_dev($aa_sum, $aa_sum_squares, $aa_count);
    my @BB_stats = st_dev($bb_sum, $bb_sum_squares, $bb_count);
    my @AB_stats = st_dev($ab_sum, $ab_sum_squares, $ab_count);
    $AA_BOUND         = $AA_stats[0] + ($AA_stats[1] * $HOM_SD);
    $BB_BOUND         = $BB_stats[0] - ($BB_stats[1] * $HOM_SD);
    $AB_UPPER_BOUND   = $AB_stats[0] + ($AB_stats[1] * $HET_SD);
    $AB_LOWER_BOUND   = $AB_stats[0] - ($AB_stats[1] * $HET_SD);
    $MI1_UPPER_THRESH = $AB_stats[0] + ($AB_stats[1] * 5);
    $MI1_LOWER_THRESH = $AB_stats[0] - ($AB_stats[1] * 5);
    @AA_stats = @AB_stats = @BB_stats = ();
 #print "@temp_bounds $AA_BOUND $BB_BOUND $AB_UPPER_BOUND $AB_LOWER_BOUND $MI1_UPPER_THRESH $MI1_LOWER_THRESH\n";
    foreach my $lrr_ref (@lrr_normal) {
        push(@p1_lrr_array, $$lrr_ref[3]);
        push(@p2_lrr_array, $$lrr_ref[4]);
        push(@ch_lrr_array, $$lrr_ref[5]);
    }

    $LOCAL_P1_LRR_MED     = median(\@p1_lrr_array);
    $LOCAL_P2_LRR_MED     = median(\@p2_lrr_array);
    $LOCAL_CH_LRR_MED     = median(\@ch_lrr_array);    
    $LOCAL_CH_BAF_MEAN    = $ch_baf_sum / $ch_baf_ct;
    my $local_p1_baf_mean = $p1_baf_sum / $p1_baf_ct;
    my $local_p2_baf_mean = $p2_baf_sum / $p2_baf_ct;
    map {$_ = 1 - $_ if ($_ < $LOCAL_CH_BAF_MEAN)} @ch_mbaf_array;
    map {$_ = 1 - $_ if ($_ < $local_p1_baf_mean)} @p1_mbaf_array;    
    map {$_ = 1 - $_ if ($_ < $local_p2_baf_mean)} @p2_mbaf_array; 
    $LOCAL_CH_mBAF_MED = median(\@ch_mbaf_array);
    $LOCAL_P1_mBAF_MED = median(\@p1_mbaf_array);
    $LOCAL_P2_mBAF_MED = median(\@p2_mbaf_array);
    
    my $chr = $baf_normal[0]->[CHR];
    
	if ($BOUNDARY_EXTENSION) {
		my $stop   = 0;
		my $errors = 1;
		my $sidak_alpha = 1 - (1 - $ALPHA)**(1 / $BOUNDARY_EXTENSION);
		#print OUTPUT_FILE "alpha $sidak_alpha $ERROR_RATE\n";
		until ($stop) {
			my ($n, $k, $p) = ($SIZE_REGION, $errors, $ERROR_RATE);
			my $prob = 0;
			for ($k..$n) {
				my $bin_coeff = calculate_bin_coeff($n, $_);
				#print OUTPUT_FILE $bin_coeff, "\n";
				#print OUTPUT_FILE "dollar1 $_\n";

				$prob += $bin_coeff * ($p**$_) * ((1 - $p)**($n - $_));
				#print OUTPUT_FILE $bin_coeff * ($p**$_) * ((1 - $p)**($n - $_)), "\n";
				#print OUTPUT_FILE "prob $prob $n $_\n";
				#print OUTPUT_FILE "dollar1 $_\n";

			}
			if ($prob <= $sidak_alpha) {
				$stop = 1;
				$ACCEPTABLE_ERRORS = $errors - 1;
				$ACCEPTABLE_ERRORS ||= 0;
			}
			$errors++;
		}
	}
	
    (%BAF_OUTLIERS_BY_POS, %CH_HD_BY_POS, %MI1_BY_POS, %P1_HD_BY_POS, 
	 %P1_INF_SNP_BY_POS, %P2_INF_SNP_BY_POS, %P2_HD_BY_POS) = (()) x 7; 
    (@baf_normal, @lrr_normal, @p1_lrr_array, @p2_lrr_array, @ch_lrr_array,
     @ch_mbaf_array, @BAF_OUTLIERS, @P1_INF_SNPS, @P2_INF_SNPS) = (()) x 9; 

    $REFINING = 1;
    # Reprocess chromosome
    my $return = process_chromosome_arm();
    map {$$_[19] = 2} @{$$return[0]};
    
    # Store local stats for each chromosome arm for which they are calculated
    push (@$return, $chr, [($LOCAL_P1_mBAF_MED, $LOCAL_P2_mBAF_MED,
        $LOCAL_CH_mBAF_MED, $LOCAL_P1_LRR_MED, $LOCAL_P2_LRR_MED, 
        $LOCAL_CH_LRR_MED)]);

    # Restore initial thresholds     
    ($AA_BOUND, $BB_BOUND, $AB_UPPER_BOUND, $AB_LOWER_BOUND, 
		$MI1_UPPER_THRESH, $MI1_LOWER_THRESH, $ACCEPTABLE_ERRORS) = @temp_bounds;
    return($return);
}

sub median {
    my $array_ref = $_[0];
    my $median = 0;
    if (!$$array_ref[0] || @$array_ref == 1) {$median = NA}
    elsif (@$array_ref == 2) {$median = ($$array_ref[0] + $$array_ref[1]) / 2}
    else {
        my @sorted_array = sort { $a <=> $b } @$array_ref;
        my $mid_index = @sorted_array / 2;
        if (@sorted_array % 2) {
            $median = $sorted_array[int(($mid_index) - 0.5)];
        }
        else {
            $median = ($sorted_array[$mid_index] 
                + $sorted_array[$mid_index - 1]) / 2;
        }
    }
    return($median);
}

sub overlap {
    my ($start1, $stop1, $start2, $stop2) = @_;
    my ($overlap_start, $overlap_stop, $end1_start, $end1_stop, 
        $end2_start, $end2_stop, $prev_snp, $next_snp, $type) = (0) x 9;
    # -----
    # -----       
    if ($start1 == $start2 && $stop1 == $stop2) {
        ($overlap_start, $overlap_stop)  = ($start1, $stop1);
        $type = 1;
    }
    #   -----
    # -----
    elsif ($start1 > $start2 && $stop2 >= $start1 && $stop2 < $stop1) {
        $prev_snp = $SNP_BY_NUM{$SNP_BY_POS{$start1} - 1};
        $next_snp = $SNP_BY_NUM{$SNP_BY_POS{$stop2} + 1};
        ($end1_start, $end1_stop)       = ($start2, $prev_snp);
        ($overlap_start, $overlap_stop) = ($start1, $stop2);
        ($end2_start, $end2_stop)       = ($next_snp, $stop1);
        $type = 2;
    }
    # -----
    # ---
    elsif ($start1 == $start2 && $stop2 > $start1 && $stop2 < $stop1) {
        $next_snp = $SNP_BY_NUM{$SNP_BY_POS{$stop2} + 1};
        ($overlap_start, $overlap_stop) = ($start2, $stop2);
        ($end2_start, $end2_stop)       = ($next_snp, $stop1);
        $type = 3;
    }
    #   ---
    # -----
    elsif ($start1 > $start2 && $stop2 > $start1 && $stop2 == $stop1) {
        $prev_snp = $SNP_BY_NUM{$SNP_BY_POS{$start1} - 1};
        ($end1_start, $end1_stop)       = ($start2, $prev_snp);
        ($overlap_start, $overlap_stop) = ($start1, $stop1);
        $type = 4;
    }
    #  ---
    # -----
    elsif ($start1 > $start2 && $stop1 < $stop2) { 
        $prev_snp = $SNP_BY_NUM{$SNP_BY_POS{$start1} - 1};
        $next_snp = $SNP_BY_NUM{$SNP_BY_POS{$stop1}  + 1};        
        ($end1_start, $end1_stop)       = ($start2, $prev_snp);
        ($overlap_start, $overlap_stop) = ($start1, $stop1);
        ($end2_start, $end2_stop)       = ($next_snp, $stop2);
        $type = 5;        
    }
    # -----
    #   -----
    elsif ($start2 > $start1 && $stop1 >= $start2 && $stop1 < $stop2) {
        $prev_snp = $SNP_BY_NUM{$SNP_BY_POS{$start2} - 1};
        $next_snp = $SNP_BY_NUM{$SNP_BY_POS{$stop1}  + 1};        
        ($end1_start, $end1_stop)       = ($start1, $prev_snp);
        ($overlap_start, $overlap_stop) = ($start2, $stop1);
        ($end2_start, $end2_stop)       = ($next_snp, $stop2);
        $type = 6;        
    }
    # ---
    # -----
    elsif ($start2 == $start1 && $stop1 > $start2 && $stop1 < $stop2) {
        $next_snp = $SNP_BY_NUM{$SNP_BY_POS{$stop1}  + 1};        
        ($overlap_start, $overlap_stop) = ($start1, $stop1);
        ($end2_start, $end2_stop)       = ($next_snp, $stop2);
        $type = 7;        
    }
    # -----
    #   ---
    elsif ($start2 > $start1 && $stop1 > $start2 && $stop1 == $stop2) {
        $prev_snp = $SNP_BY_NUM{$SNP_BY_POS{$start2} - 1};
        ($end1_start, $end1_stop)       = ($start1, $prev_snp);
        ($overlap_start, $overlap_stop) = ($start2, $stop2);
        $type = 8;        
    }
    # -----
    #  ---
    elsif ($start2 > $start1 && $stop2 < $stop1) {
        $prev_snp = $SNP_BY_NUM{$SNP_BY_POS{$start2} - 1};
        $next_snp = $SNP_BY_NUM{$SNP_BY_POS{$stop2}  + 1};        
        ($end1_start, $end1_stop)       = ($start1, $prev_snp);
        ($overlap_start, $overlap_stop) = ($start2, $stop2);
        ($end2_start, $end2_stop)       = ($next_snp, $stop1);
        $type = 9;
    }
    return([$type, $end1_start, $end1_stop, $overlap_start, $overlap_stop,
          $end2_start, $end2_stop]);
}


sub print_proper_format {
    # Prints a message describing the proper format of the input file.
    print <<END;
    
The input file must be tab delimited, sorted by chromosome and position, 
and in the following order: SNP Name, Chromosome, Position, 
Father GType, Father BAF, Father LRR, Mother GType, Mother BAF, Mother LRR, 
Child GType, Child BAF, Child LRR.
The genotypes must be AA, AB, BB, NC or NoCall.
B allele frequencies must be >= 0 and <= 1 for polymorphic markers.
A header line is expected but is not used to determine column identity.
    
END
}    

sub process_chromosome_arm {
    my ($adjusted_mi1_ct, $chr, $end, $hom_reg_ref, $inf_SNP_count, $mi1, 
        $mi1_contrib, $mi1_ct, $mi1_reg_ref, $POD_ref, $PODcr_ref, $podcr_regs,
		$start, $het_mi1_ct, $out_start, $out_end, $outlier_ct) = (0) x 17;
    my (@abnormal_regions, @adj_mi1_regions, @return, @large_regions, @outlier_regions);
    $chr = $CURR_CHR_REFS[0]->[CHR];
    
    for (my $i = 0; $i < @CURR_CHR_REFS; $i++) {
        $BAF_BY_POS{$CURR_CHR_REFS[$i]->[POS]} = $i;
		
        my $results = find_informative_snps($CURR_CHR_REFS[$i]);
		if ($$results[0]) {
		    push(@P1_INF_SNPS, $CURR_CHR_REFS[$i]);
			$P1_INF_SNP_BY_POS{$CURR_CHR_REFS[$i]->[POS]} = keys %P1_INF_SNP_BY_POS;
		}
		elsif ($$results[1]) {
		    push(@P2_INF_SNPS, $CURR_CHR_REFS[$i]);
			$P2_INF_SNP_BY_POS{$CURR_CHR_REFS[$i]->[POS]} = keys %P2_INF_SNP_BY_POS;
		}
        $inf_SNP_count += $$results[2];
        if ($$results[3] && $CURR_CHR_REFS[$i]->[-1] != UNKNOWN_MI1) {
            if ($mi1_ct && $mi1_contrib != $$results[4]) {
                if ($mi1_ct >= $MIN_MI1) {
                    push(@adj_mi1_regions, [$chr, $start, $CURR_CHR_REFS[$i]->[POS], 
						" ", $mi1_contrib, $mi1_ct, $mi1_ct]);
                }
                $mi1_ct = $start = 0;
            }
            $mi1_ct++;
            $mi1_contrib = $$results[2];
            $MI1_BY_POS{$CURR_CHR_REFS[$i]->[POS]} = keys %MI1_BY_POS;
            if ($mi1_ct == 1) {$start = $end = $CURR_CHR_REFS[$i]->[POS]}
            else {$end = $CURR_CHR_REFS[$i]->[POS]}
            $adjusted_mi1_ct++ if !$mi1_ct;
        }
        else {
            if ($mi1_ct >= $MIN_MI1) {
                push(@adj_mi1_regions, [$chr, $start, $end, " ", $mi1_contrib, 
					$mi1_ct, $mi1_ct]);
            }
            $mi1_ct = $start = $end = $mi1_contrib = 0;
        }
		if ($$results[5]) {
            $outlier_ct++;
            $BAF_OUTLIERS_BY_POS{$CURR_CHR_REFS[$i]->[POS]} 
				= keys %BAF_OUTLIERS_BY_POS;
			push(@BAF_OUTLIERS, $CURR_CHR_REFS[$i]);	
            if ($outlier_ct == 1) {$out_start = $out_end = $CURR_CHR_REFS[$i]->[POS]}
			else {$out_end = $CURR_CHR_REFS[$i]->[POS]}
        }
        elsif ($CURR_CHR_REFS[$i]->[9] eq "AB") {
            if ($outlier_ct >= $MIN_BAF_OUT) {
                push(@outlier_regions, [$chr, $out_start, $out_end, 
					" ", UNKNOWN, $outlier_ct, $outlier_ct]);
            }
            $outlier_ct = $out_start = $out_end = 0;
        }
    }
	
	if ($mi1_ct >= $MIN_MI1) {
        push(@adj_mi1_regions, [$chr, $start, $end, " ", $mi1_contrib, $mi1_ct, 
			$mi1_ct]);
    }
	if ($outlier_ct >= $MIN_BAF_OUT) {
        push(@outlier_regions, [$chr, $out_start, $out_end, 
			" ", UNKNOWN, $outlier_ct, $outlier_ct]);
    }


    #map{print OUTPUT_FILE "@{$_}\n"} @adj_mi1_regions if @adj_mi1_regions;
    if ($POD) {
	#print "starting POD\n";
        ($POD_ref, my $large_windows) = find_POD_regions();
        push(@abnormal_regions, @$POD_ref) if @$POD_ref;
				#print "ending POD\n";
		#map {print OUTPUT_FILE "POD @{$_}\n"} @abnormal_regions;

		my $large_regions = find_large_mosaic_POD_regions($large_windows);
		#map {print OUTPUT_FILE "large @{$_}\n"} @$large_regions;
        push(@abnormal_regions, @$large_regions) if @$large_regions;
		#collapse_region_list(@abnormal_regions);
		#map {print OUTPUT_FILE  "large and POD @{$_}\n"} @abnormal_regions;
		#print STDERR "abn regs1 @abnormal_regions\n";
					#print "ending POD\n";

    }
	if ($PODcr && @outlier_regions) {
		#print "starting PODcr\n";
		$podcr_regs = define_PODcr_regions(\@outlier_regions);
        push(@abnormal_regions, @$podcr_regs) if @$podcr_regs;
		   # map{print STDERR "@{$_}\n"} @abnormal_regions if (@abnormal_regions);

		#print "ending PODcr\n";
	}		
    if ($MI1 && @adj_mi1_regions) {
			#print "starting mi1n";

        $mi1_reg_ref = define_mi1_regions(\@adj_mi1_regions);
        push(@abnormal_regions, @$mi1_reg_ref) if @$mi1_reg_ref;
		#print STDERR "abn regs2 @abnormal_regions\n";
				#print "ending mi1n";

    }  
	#map{print "mi1 @{$_}\n"} @$mi1_reg_ref if (@$mi1_reg_ref);
    if ($HD) {
		#print "starting HD\n";

        $hom_reg_ref = detect_hom_del();
        push(@abnormal_regions, @$hom_reg_ref) if @$hom_reg_ref;
		#print STDERR "abn regs3 @abnormal_regions\n";
			#	print "ending HD\n";

    }
    if (@abnormal_regions) {		
		count_std_mi1_overlap(\@abnormal_regions);
		   # map{print OUTPUT_FILE "@{$_}\n"} @abnormal_regions if (@abnormal_regions);
    # map{print STDERR "@{$_}\n"} @abnormal_regions if (@abnormal_regions);
					#print STDERR "process1 tocollapse\n";

        collapse_region_list(\@abnormal_regions);
		   # map{print OUTPUT_FILE "1st @{$_}\n"} @abnormal_regions if (@abnormal_regions);

        map {$_ = evaluate_region($_) unless $$_[6]} @abnormal_regions;
        map {calculate_region_stats($_) unless $$_[7]} @abnormal_regions;
				  # map{print OUTPUT_FILE "2nd @{$_}\n"} @abnormal_regions if (@abnormal_regions);

		count_std_mi1_overlap(\@abnormal_regions);	
					#print STDERR "process2 tocollapse\n";

		collapse_region_list(\@abnormal_regions);
				  # map{print OUTPUT_FILE "3rd @{$_}\n"} @abnormal_regions if (@abnormal_regions);

		map{$$_[6] = count_SNPs(@$_[START,STOP,16]) unless $$_[6]} @abnormal_regions;
        map {calculate_region_stats($_) unless $$_[7]} @abnormal_regions;
		
		for (my $i = 0; $i <= $#abnormal_regions; $i++) {
			my $ref = $abnormal_regions[$i];
			#print OUTPUT_FILE "checking for het mi1 @$ref\n";
			my $mi1_het_ct = 0;
			#Check for MI1s indicative of abn region in parent instead of child
						#print STDERR "process to find\n";
#print STDERR "@$ref[START,STOP]\n";
			my ($reg_start, $reg_stop) 
				= find_closest_indices(@$ref[START,STOP], \%BAF_BY_POS);
				#print OUTPUT_FILE "$SNP_BY_POS{$$ref[START]} $CURR_CHR_REFS[$SNP_BY_POS{$$ref[START]}]->[POS] $SNP_BY_POS{$$ref[STOP]} $CURR_CHR_REFS[$SNP_BY_POS{$$ref[STOP]}]->[POS]\n"; 
 				#print OUTPUT_FILE "reg start stop @$ref[START,STOP] $CURR_CHR_REFS[$reg_start]->[POS] $CURR_CHR_REFS[$reg_stop]->[POS]\n";
				
			if ($reg_start && $reg_stop && $reg_stop - $reg_start > 0) {
				for (my $j = $reg_start; $j <= $reg_stop; $j++) {
					$mi1_het_ct++ if ($CURR_CHR_REFS[$j]->[-1] == UNKNOWN_MI1);
				}
				#print OUTPUT_FILE "het ct $mi1_het_ct @$ref\n";
				if ($mi1_het_ct) {
					my $thresh = 0.001;
					my ($n, $k, $p) = ($reg_stop - $reg_start, $mi1_het_ct, $HET_MI1_RATE);
					#print OUTPUT_FILE "here $mi1_het_ct $n $k $p @$ref\n";		

					my $approx = approx_power_pdf($n, $k, $p);
					#print OUTPUT_FILE "approx $approx\n";
					if (log10($thresh) - $approx + 1 < 15) {
						my $bin_coeff = calculate_bin_coeff($n, $k);
						unless ($bin_coeff eq "nan") {
							my $prob = $bin_coeff * ($p**$k) * ((1 - $p)**($n - $k));
							#my $prob = pdf($n, $k, $p);	
							next if ($prob < $thresh);
							#if ($prob < $thresh) {print OUTPUT_FILE "if next\n"; next;}
							
							#print OUTPUT_FILE "prob $mi1_het_ct $n $k $p $bin_coeff $prob @$ref\n";		
							#print OUTPUT_FILE "prob $mi1_het_ct $n $k $p $prob @$ref\n";
						}
					}
					else {next}
				}
			}

			if ($$ref[17] == HD) {
				my $hd_thresh;
				if    ($$ref[4] == FATHER) {$hd_thresh = $MIN_P2_HD}
				elsif ($$ref[4] == MOTHER) {$hd_thresh = $MIN_P1_HD}
				elsif ($$ref[4] == NONE)   {$hd_thresh = $MIN_CH_HD}
				elsif ($$ref[4] == BOTH)   {
					if ($MIN_P1_HD <= $MIN_P2_HD) {$hd_thresh = $MIN_P1_HD}
					else {$hd_thresh = $MIN_P2_HD}
				}
				push(@return, $ref) if ($$ref[6] >= $hd_thresh);
			}
			elsif ($$ref[17] == MI1) {
				push(@return, $ref) if ($$ref[6] >= $MIN_MI1);
			}
			elsif ($$ref[17] == POD || $$ref[17] == PODcr) {
			#	push(@return, $ref) if ($$ref[6] >= $MIN_POD);
				push(@return, $ref);
			}
		}
		map {$$_[20] = call_inheritance($_)} @return;
    }
    my $results = [(\@return, $inf_SNP_count, $adjusted_mi1_ct)];
    return($results);
}

sub splice_overlap_regions {
	# Splices two overlapping regions depending on type of overlap 
    my ($reg1, $reg2) = @_;
    my (@results, @overlap);
    
    if (!@$reg1 && !@$reg2) {return(\@results)}
    elsif (!@$reg1) {return($reg2)}
    elsif (!@$reg2) {return($reg1)}

	my ($chr, $start1, $stop1, $type1, $contrib1, $size1, $inf1, $detect1
		) = @$reg1[0..6,17];
#print STDERR "reg1 @$reg1\n";
#print STDERR "reg2 @$reg2\n";
	my ($start2, $stop2, $type2, $contrib2, $size2, $inf2, $detect2
		) = @$reg2[1..6,17]; 
	my ($orient, $end1_start, $end1_stop, $overlap_start, $overlap_stop,
		$end2_start, $end2_stop, $no_split) 
		= @{overlap($start1, $stop1, $start2, $stop2)};
		#print "$chr $start1, $stop1 $start2 $stop2 $$reg2[18] $$reg2[18]\n";
	# Combinations of contribution and type 
	my $eq = ($contrib1 == $contrib2 && $size1 == $size2 
		&& $inf1 == $inf2 && $detect1 == $detect2);
	my $mi1_1_trumps = $detect2 == POD && $detect1 == MI1 
		&& $$reg2[18] <= 2 && $inf1 >= 5 && (($inf2 - $inf1) < $MIN_POD);
	my $mi1_2_trumps = $detect1 == POD && $detect2 == MI1 
		&& $$reg1[18] <= 2 && $inf2 >= 5 && (($inf1 - $inf2) < $MIN_POD);
	my $combine1 = !$mi1_2_trumps && $contrib1 == $contrib2;
	my $combine2 = !$mi1_1_trumps && $contrib1 == $contrib2;			
	my $reg1_trumps = ($detect1 == POD && !$mi1_2_trumps) 
		|| ($contrib1 == $contrib2 && $detect1 == $detect2 
			&& ((!$inf2 || $size1 > $size2) 
			|| ($size1 == $size2 && $inf1 > $inf2)))
		|| $contrib1 == NONE && $contrib2 != NONE 
		|| $contrib2 == UNKNOWN
		|| $mi1_1_trumps;
	my $reg2_trumps = ($detect2 == POD && !$mi1_1_trumps) 
		|| ($contrib1 == $contrib2 && $detect1 == $detect2
			&& ((!$inf1 || $size2 > $size1) 
			|| ($size1 == $size2 && $inf2 > $inf1)))
		|| $contrib2 == NONE && $contrib1 != NONE 
		|| $contrib1 == UNKNOWN
		|| $mi1_2_trumps;    
	
	my %empty;
	# Output combinations
	my $combo111 = [$chr, $end1_start, $end1_stop, $type1, $contrib1,
					(0) x 11, $$reg1[16], $detect1];
	my $combo112 = [$chr, $end1_start, $end1_stop, $type2, $contrib2,
					(0) x 11, $$reg2[16], $detect2];
	my $combo121 = [$chr, $end1_start, $end2_stop, $type1, $contrib1,
					(0) x 11, $$reg1[16], $detect1];
	my $combo122 = [$chr, $end1_start, $end2_stop, $type2, $contrib2,
					(0) x 11, $$reg2[16], $detect2];
	my $combo221 = [$chr, $end2_start, $end2_stop, $type1, $contrib1,
					(0) x 11, $$reg1[16], $detect1];
	my $combo222 = [$chr, $end2_start, $end2_stop, $type2, $contrib2,
					(0) x 11, $$reg2[16], $detect2];
	my @both_contrib;
	if ($detect2 == POD || $inf2 > $inf1) {@both_contrib = ($$reg2[16], $detect2)}
	else {@both_contrib = ($$reg1[16], $detect1)}
	my $comboOOB = [$chr, $overlap_start, $overlap_stop, " ", BOTH,
					(0) x 11, @both_contrib[0,1]];

	# -----
	# -----
	if ($orient == 1) {
		$no_split = 1;
		if ($eq || $reg1_trumps) {push(@results, $reg1)}
		elsif ($reg2_trumps)     {push(@results, $reg2)}                
		else {push(@results, $comboOOB)}
	}
	#   -----
	# -----
	elsif ($orient == 2) {
		if ($eq || $combine1) {push(@results, $combo121); $no_split = 1}
		elsif ($combine2)     {push(@results, $combo122); $no_split = 1}
		elsif ($reg1_trumps)  {push(@results, $reg1, $combo112)}
		elsif ($reg2_trumps)  {push(@results, $reg2, $combo221)}
		else {push(@results, $combo112, $comboOOB, $combo221)} 
	}
	# -----
	# ---
	elsif ($orient == 3) {
		if ($eq || $reg1_trumps) {push(@results, $reg1); $no_split = 1}
		elsif ($reg2_trumps) {push(@results, $reg2, $combo221)}
		else {push(@results, $comboOOB, $combo221)}
	}
	#   ---
	# -----
	elsif ($orient == 4) {
		if ($eq || $reg2_trumps) {push(@results, $reg2); $no_split = 1}
		elsif ($reg1_trumps) {push(@results, $reg1, $combo112)}
		else {push(@results, $combo112, $comboOOB)}
	}
	#  ---
	# -----
	elsif ($orient == 5) { 
		if ($eq || $reg2_trumps) {push(@results, $reg2); $no_split = 1}
		elsif ($reg1_trumps) {push(@results, $reg1, $combo112, $combo222)}
		else {push(@results, $combo112, $comboOOB, $combo222)}
	}
	# -----
	#   -----
	elsif ($orient == 6) {
		if ($eq || $combine1) {push(@results, $combo121); $no_split = 1}
		elsif ($combine2)     {push(@results, $combo122); $no_split = 1}
		elsif ($reg1_trumps)  {push(@results, $reg1, $combo222)}
		elsif ($reg2_trumps)  {push(@results, $reg2, $combo111)}
		else {push(@results, $combo111, $comboOOB, $combo222)}
	}
	# ---
	# -----
	elsif ($orient == 7) {
		if ($eq || $reg2_trumps) {push(@results, $reg2); $no_split = 1}
		elsif ($reg1_trumps) {push(@results, $reg1, $combo222)}
		else {push(@results, $comboOOB, $combo222)}
	}
	# -----
	#   ---
	elsif ($orient == 8) {
		if ($eq || $reg1_trumps) {push(@results, $reg1); $no_split = 1}
		elsif ($reg2_trumps) {push(@results, $reg2, $combo111)}
		else {push(@results, $combo111, $comboOOB)}
	}
	# -----
	#  ---
	elsif ($orient == 9) {
		if ($eq || $reg1_trumps) {push(@results, $reg1); $no_split = 1}
		elsif ($reg2_trumps) {push(@results, $reg2, $combo111, $combo221)}
		else {push(@results, $combo111, $comboOOB, $combo221)}
	}
	else {push(@results, $reg1, $reg2); $no_split = 1}
	
	map {$$_[6] = count_SNPs(@$_[START,STOP,16])} @results;
	#map {print "results @{$_}\n"} @results;

	my @return;
	
	if ($no_split) {@return = @results}
	else {
		for my $ref (@results) {
			if ($$ref[17] == HD) {
				my $hd_thresh;
				if    ($$ref[4] == FATHER) {$hd_thresh = $MIN_P2_HD}
				elsif ($$ref[4] == MOTHER) {$hd_thresh = $MIN_P1_HD}
				elsif ($$ref[4] == NONE)   {$hd_thresh = $MIN_CH_HD}
				elsif ($$ref[4] == BOTH)   {
					if ($MIN_P1_HD <= $MIN_P2_HD) {$hd_thresh = $MIN_P1_HD}
					else {$hd_thresh = $MIN_P2_HD}
				}
				push(@return, $ref) if ($$ref[6] >= $hd_thresh);
			}
			elsif ($$ref[17] == MI1) {
				push(@return, $ref) if ($$ref[6] >= $MIN_MI1);
			}
			elsif ($$ref[17] == POD || $$ref[17] == PODcr) {
				#print "$MIN_POD $$ref[6]\n";
				push(@return, $ref) if ($$ref[6] >= $MIN_POD);
			}
			else {push(@return, $ref)}
		}
	}
	
	#map {print "results @{$_}\n"} @return;
    return(\@return);
}

sub st_dev {
    # Calculates standard deviation.
	# Takes an array ref or sum, sum_squares, count
    # Returns the mean, standard deviation, variance, and 
    # coefficient of variance
	my ($in1, $in2, $in3) = @_;
    my ($sum, $sum_squares, $ct, $mean, $variance, $coef_of_var, 
		$st_dev) = (0) x 7;    
	
	if ($in2) {($sum, $sum_squares, $ct) = @_}
	else {
		$sum += $_ for @$in1;
		$sum_squares += $_**2 for @$in1;
		$ct = @$in1;
	}
    if ($ct > 1) { 
        $mean = $sum / $ct;
        if ($mean) {
            $variance = ($sum_squares - ($sum**2 / $ct)) / ($ct - 1);      
            if ($variance < 0 && $variance > -0.000001) {
                $variance = $st_dev = $coef_of_var = 0; 
            }
            elsif ($variance < -0.000001) {
                print STDERR "SAMPLE $FILENAME : A negative variance value ",
					"has been calculated for global BAF!\n"; 
            }
            else {
                $st_dev = sqrt($variance); 
                $coef_of_var = $st_dev / $mean;
            }
        }
    } 
    return($mean, $st_dev, $variance, $coef_of_var);
}

sub start_R {
    # Makes a system call to kick off the Rscript.
    my $r_script = join(" ", "Rscript --slave --no-save --no-restore",
        "--no-environ --silent",
        "/data/Joe/Mosaicism/SNP-POA/SNP-POA_Rscript_cleanv5.R",
        $OUTPUT_DIR, $INPUT_FILE, $FILENAME, $P1_NAME, 
        $P2_NAME, $CH_NAME, $PID_VALUE, $R_PID_FILENAME,
        $PERL_TO_R_FILENAME, $R_TO_PERL_FILENAME, "> $R_LOG_FILENAME 2>&1 &");
    system($r_script);
}

sub streak_prob {
	# Calculate the probability of a streak size >=k in n trials
	# P[n,k] = p^k + sum(from j=1,k) of p^(j-1) (1-p) P[n-j,k]
	# P[n,k] = 0 if n=0 or k>n
	# P[n,k] = p^k for n=k
	my($n, $minK, $p) = @_;
	my (@array, $prevS);
	if (!$n || $n < $minK) {return(0)}
	elsif ($n == $minK) {return($p**$minK)}
	my $result = $p**$minK;

	for (my $i = 0; $i <= $n - 1; $i++) {
		if ($i < $minK) {$STREAK_ARRAY[$i] = 0}
		elsif ($i == $minK) {$STREAK_ARRAY[$i] = $p**$minK}
		else {
			$STREAK_ARRAY[$i] = $p**$minK;
			$STREAK_ARRAY[$i] += $p**($_ - 1) * (1 - $p) * $STREAK_ARRAY[$i - $_] for(1..$minK);
		}
	}
	$result += $p**($_-1) * (1 - $p) * $STREAK_ARRAY[$n - $_] for(1..$minK);		
	return($result);
}
