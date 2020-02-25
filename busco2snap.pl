#!/usr/bin/env perl

use strict;
use warnings;
use Cwd 'abs_path';
use IPC::Cmd qw[can_run run];
use Parallel::Loops;

my $version = "0.1";

my $busco_dir = "";
my $full_table = "";
my $retraining_parameters = "";
my $species = "";
my $duplicated = 0;
my $fragmented = 0;
my $num_s = 0;
my $num_d = 0;
my $num_f = 0;
my $threads = 1;
my $out_dir = abs_path("./");
my $prefix = "";
my $keep_tmp = 0;
my $verbose = 0;
my $dry = 0;
my $cmd;
my @arg_cp = @ARGV;
my @error_args = ();
my $length = 0;
my %genome_ann;

my $input_error = 0;

sub print_help{
	print STDOUT "\n";
	print STDOUT "busco2snap.pl v$version\n";
	print STDOUT "\n";
	print STDOUT "Description:\n";
	print STDOUT "\tCreation of a SNAP model from BUSCOs by predicting the gene structure with the retrained\n\tAugustus model, converting to zff format. The tools augustus and blastp need to be in\n\tyour \$PATH. If fathom, forge and hmm-assembler.pl are in your \$PATH the SNAP model will be\n\tcreated automatically.\n";
	print STDOUT "\tIf Augustus predicts more than one gene at a BUSCO locus the predicted proteins are blasted\n\tagainst the ancestral variants of the corresponding BUSCO lineage. The gene(s) where the best\n\thit (defined by highest bit score) is on the searched BUSCO are used only.\n";
	print STDOUT "\n";
	print STDOUT "Usage:\n";
	print STDOUT "\tbusco2snap.pl [-b <busco_dir> | -ft <full_table> -rp <retraining_parameters>]\n";
	print STDOUT "\n";
	print STDOUT "Mandatory:\n";
	print STDOUT "\t-b STR\t\tBUSCO output directory (run_*)\n";
	print STDOUT "\tOR\n";
	print STDOUT "\t-ft STR\t\tBUSCOs full table file (full_table_*)\n";
	print STDOUT "\t-rp STR\t\tDirectory containing retrining parameters from Augustus\n\t\t\t(augustus_output/retraining_parameters/)\n";
	print STDOUT "\n";
	print STDOUT "Options: [default]\n";
	print STDOUT "\t-dup\t\tInclude duplicated BUSCOs [off]\n";
	print STDOUT "\t-frag\t\tInclude fragmented BUSCOs [off]\n";
	print STDOUT "\t-o STR\t\tOutput directory [.]\n";
	print STDOUT "\t-pre STR\tPrefix for file name of the SNAP model file [species model name]\n";
	print STDOUT "\t-t INT\t\tNumber of parallel processed BUSCOs [1]\n";
	print STDOUT "\t-kt\t\tKeep temporary files [off]\n";
	print STDOUT "\t-v\t\tPrint executed commands to STDERR [off]\n";
	print STDOUT "\t-dry-run\tOnly print commands to STDERR instead of executing [off]\n";
	print STDOUT "\n";
	print STDOUT "\t-h or -help\tPrint this help and exit\n";
	print STDOUT "\t-version\tPrint version number and exit\n";
	exit;
}

sub exe_cmd{
	my ($cmd,$verbose,$dry) = @_;
	if($verbose == 1){
		print STDERR "CMD\t$cmd\n";
	}
	if($dry == 0){
		system("$cmd") == 0 or die "ERROR\tsystem $cmd failed: $?";
	}
}

sub empty_element{
	my ($i,$l) = @_;
	for (my $x = 0; $x <= $l-1; $x++){
		$arg_cp[$i+$x] = "";
	}
}

sub check_rp{
	my ($rp_dir) = @_;
	opendir (DIR, $rp_dir) or die "ERROR\tcould not open directory $rp_dir\n";
	
	my $base = "";
	while (my $file = readdir(DIR)){
		if($file =~ m/^\./){
			next;
		}
		my @model_part = split(/_/,$file);
		if($model_part[-1] eq "probs.pbl"){
			if($model_part[-2] eq "intron" or
			   $model_part[-2] eq "exon" or
			   $model_part[-2] eq "igenic"){
				my @b = @model_part;
				splice(@b,-2);
				if($base eq ""){
					$base = join("_",@b);
				}
				else{
					if($base ne join("_",@b)){
						print STDERR "ERROR\tAugustus species name differs from $base for file $rp_dir/$file\n";
						$input_error = 1;
					}
				}
				next;
			}
			else{
				print STDERR "ERROR\tUnrecognized Augustus model part $file in $rp_dir\n";
				$input_error = 1;
			}
		}
		if($model_part[-1] eq "metapars.cfg" or
		   $model_part[-1] eq "metapars.utr.cfg" or
		   $model_part[-1] eq "parameters.cfg" or
		   $model_part[-1] eq "weightmatrix.txt" or
		   $model_part[-1] eq "metapars.cgp.cfg" or
		   $model_part[-1] eq "parameters.cfg.orig1"){
			my @b = @model_part;
			splice(@b,-1);
			if($base eq ""){
				$base = join("_",@b);
			}
			else{
				if($base ne join("_",@b)){
					print STDERR "ERROR\tAugustus species name differs from $base for file $rp_dir/$file\n";
					$input_error = 1;
				}
			}
			$base = join("_",@b);
			next;
		}
		else{
			print STDERR "ERROR\tUnrecognized Augustus model part $file in $rp_dir\n";
			$input_error = 1;
		}
	}
	return $base;
}

if(scalar(@ARGV==0)){
	print_help;
}

for (my $i = 0; $i < scalar(@ARGV);$i++){
	if ($ARGV[$i] eq "-b"){
		$busco_dir = abs_path($ARGV[$i+1]);
		empty_element($i,2);
	}
	if ($ARGV[$i] eq "-ft"){
		$full_table = abs_path($ARGV[$i+1]);
		empty_element($i,2);
	}
	if ($ARGV[$i] eq "-rp"){
		$retraining_parameters = abs_path($ARGV[$i+1]);
		empty_element($i,2);
	}
	if ($ARGV[$i] eq "-dup"){
		$duplicated = 1;
		empty_element($i,1);
	}
	if ($ARGV[$i] eq "-frag"){
		$fragmented = 1;
		empty_element($i,1);
	}
	if ($ARGV[$i] eq "-t"){
		$threads = $ARGV[$i+1];
		empty_element($i,2);
	}
	if ($ARGV[$i] eq "-o"){
		$out_dir = abs_path($ARGV[$i+1]);
		empty_element($i,2);
	}
	if ($ARGV[$i] eq "-pre"){
		$prefix = $ARGV[$i+1];
		empty_element($i,2);
	}
	if ($ARGV[$i] eq "-kt"){
		$keep_tmp = 1;
		empty_element($i,1);
	}
	if ($ARGV[$i] eq "-v"){
		$verbose = 1;
		empty_element($i,1);
	}
	if ($ARGV[$i] eq "-dry-run"){
		$verbose = 1;
		$dry = 1;
		empty_element($i,1);
	}
	if ($ARGV[$i] eq "-h" or $ARGV[$i] eq "-help"){
		print_help;
	}
	if ($ARGV[$i] eq "-version"){
		print STDERR $version . "\n";
		exit;
	}
}

for(my $i = 0; $i < scalar(@arg_cp); $i++){
	if($arg_cp[$i] ne ""){
		push(@error_args,$arg_cp[$i]);
	}
}

if(scalar(@error_args) > 0){
	print STDERR "ERROR\tunrecognized prarameters: " . join(" ",@error_args) . "\n";
	$input_error = 1;
}

if($full_table eq "" and $retraining_parameters eq "" and $busco_dir ne ""){

	$full_table = `find $busco_dir -mindepth 1 -maxdepth 1 -type f -name "full_table_*"`;
	chomp $full_table;
	if($full_table eq ""){
		print STDERR "ERROR\tfull_table not found\n";
		$input_error = 1;
	}
	
	$retraining_parameters = `find $busco_dir -mindepth 2 -maxdepth 2 -type d -name "retraining_parameters"`;
	chomp $retraining_parameters;
	if($retraining_parameters eq ""){
		print STDERR "ERROR\tretraining_parameters not found\n";
		$input_error = 1;
	}
	else{
		$species = check_rp($retraining_parameters);
	}
}
else{
	if($full_table eq ""){
		print STDERR "ERROR\t-ft not specified\n";
		$input_error = 1;
	}
	else{
		if(not -f $full_table){
			print STDERR "ERROR\t$full_table is not a file\n";
			$input_error = 1;
		}
	}
	if($retraining_parameters eq ""){
		print STDERR "ERROR\t-rp not specified\n";
		$input_error = 1;
	}
	else{
		if(not -d $retraining_parameters){
			print STDERR "ERROR\t$retraining_parameters is not a directory\n";
			$input_error = 1;
		}
		else{
			$species = check_rp($retraining_parameters);
		}
	}
}

if($threads !~ m/^\d+$/ or $threads <= 0){
	print STDERR "ERROR\t-t is not an integer >=0\n";
	$input_error = 1;
}

if(-f $out_dir){
	print STDERR "ERROR\t$out_dir is already a file\n";
	$input_error = 1;
}

if(not defined(can_run("augustus"))){
	print STDERR "ERROR\taugustus is not in your \$PATH\n";
	$input_error = 1;
}

if(not defined(can_run("blastp"))){
	print STDERR "ERROR\tblastp is not in your \$PATH\n";
	$input_error = 1;
}

if(not defined($ENV{'AUGUSTUS_CONFIG_PATH'})){
	print STDERR "ERROR\t\$AUGUSTUS_CONFIG_PATH needs to be set\n";
	$input_error = 1;
}

if($input_error == 1){
	print STDERR "ERROR\tinput error detected\n";
	exit 1;
}


if($prefix eq ""){
	$prefix = $species;
	print STDERR "INFO\tChanging prefix to $prefix\n";
}

if($prefix =~ m/[\/><|:&\\,;?*]/){
	$prefix =~ tr/[\/><|:&\\,;?*]/_/;
	print STDERR "INFO\tChanging prefix to $prefix\n";
}

if(not -e $ENV{'AUGUSTUS_CONFIG_PATH'} . "/species/" . $species){
	$cmd = "ln -s $retraining_parameters $ENV{'AUGUSTUS_CONFIG_PATH'}/species/$species";
	exe_cmd($cmd,$verbose,$dry);
}
else{
	if(-l $ENV{'AUGUSTUS_CONFIG_PATH'} . "/species/" . $species){
		my $target = readlink($ENV{'AUGUSTUS_CONFIG_PATH'} . "/species/" . $species);
		$target =~ s/\/$//;
		if($target ne $retraining_parameters){
			print STDERR "ERROR\texisting symlink $ENV{'AUGUSTUS_CONFIG_PATH'}/species/$species with target different from $retraining_parameters\n";
			exit 1;
		}
	}
	else{
		if(-d $ENV{'AUGUSTUS_CONFIG_PATH'} . "/species/" . $species){
			print STDERR "INFO\t$ENV{'AUGUSTUS_CONFIG_PATH'}/species/$species directory exists. Checking if files are the same...\n";
			my $md5sum = `md5sum $ENV{'AUGUSTUS_CONFIG_PATH'}/species/$species/* | awk '{print \$1}' | md5sum | awk '{print \$1}'`;
			chomp $md5sum;
			my $rp_md5sum = `md5sum $retraining_parameters/* | awk '{print \$1}' | md5sum | awk '{print \$1}'`;
			chomp $rp_md5sum;
			if($md5sum ne $rp_md5sum){
				print STDERR "ERROR\tAugustus model in $ENV{'AUGUSTUS_CONFIG_PATH'}/species/$species is different from $retraining_parameters\n";
				exit 1;
			}
			else{
				print STDERR "INFO\t...Passed\n";
			}
		}
		else{
			print STDERR "ERROR\t$ENV{'AUGUSTUS_CONFIG_PATH'}/species/$species exists and is no symlink nor a directory.\n";
			exit 1;
		}
	}
}

if(not -d $out_dir){
	$cmd = "mkdir -p $out_dir";
	exe_cmd($cmd,$verbose,$dry);
}

my %buscos;
my %target_ids;
my $target_fasta = "";
my $lineage_path = "";

print STDERR "INFO\tExtracting BUSCOs and their locations from $full_table\n";

open (FT, '<', $full_table) or die "ERROR\tcould not open $full_table\n";

while (my $line = <FT>){
	if($line =~ m/^#/){
		if($line =~ m/^# To reproduce this run: /){
			my @run = split(/ /,$line);
			for(my $i = 0; $i < scalar(@run); $i++){
				if($run[$i] eq "-i"){
					$target_fasta = $run[$i+1];
				}
				if($run[$i] eq "-l"){
					$lineage_path = $run[$i+1];
					if(substr($lineage_path,-1) eq "/"){
						$lineage_path = substr($lineage_path,0,-1);
					}
				}
			}
			next;
		}
		else{
			next;
		}
	}
	chomp $line;
	my @full_table_line = split(/\t/,$line);
	if($full_table_line[1] eq "Complete"){
		$num_s++;
		$buscos{$full_table_line[0]} = join("\t",@full_table_line[1..4]);
		$target_ids{$full_table_line[2]} = 1;
	}
	if($duplicated == 1){
		if($full_table_line[1] eq "Duplicated"){
			$num_d++;
			if(exists($buscos{$full_table_line[0]})){
				$buscos{$full_table_line[0]} = $buscos{$full_table_line[0]} . "\n" . join("\t",@full_table_line[1..4]);
			}
			else{
				$buscos{$full_table_line[0]} = join("\t",@full_table_line[1..4]);
			}
			$target_ids{$full_table_line[2]} = 1;
		}
	}
	if($fragmented == 1){
		if($full_table_line[1] eq "Fragmented"){
			$num_f++;
			if(exists($buscos{$full_table_line[0]})){
				$buscos{$full_table_line[0]} = $buscos{$full_table_line[0]} . "\n" .  join("\t",@full_table_line[1..4]);
			}
			else{
				$buscos{$full_table_line[0]} = join("\t",@full_table_line[1..4]);
			}
			$target_ids{$full_table_line[2]} = 1;
		}
	}
}

close FT;

my $num_target_ids = scalar(keys(%target_ids));
my $num_buscos = scalar(keys(%buscos));
my $num_c = $num_s + $num_d;

print STDERR "INFO\tExtracting $num_target_ids sequences with C:$num_c\[S:$num_s,D:$num_d],F:$num_f,n:$num_buscos BUSCOs from $target_fasta\n";

my @out_fa = ();
my $fa_progress = 0;
my $fa_frac = 0;

if($dry == 0){
	open (FA, '<', $target_fasta) or die "ERROR\tcould not open $target_fasta\n";

	my $head = "";
	my $out_file = "";

	while (my $line = <FA>){
		chomp $line;
		if($line =~ m/^>/){
			if($out_file ne ""){
				close OUT;
			}
			$head = $line;
			$head =~ s/^>//;
			$head =~ s/ .*//;
			if(exists($target_ids{$head})){
				if($verbose == 0){
					$fa_progress++;
					$fa_frac = sprintf "%.1f", ($fa_progress / $num_target_ids) * 100;
					print STDERR $fa_progress . " / " . $num_target_ids . " [" . $fa_frac . "%]\r";
				}
				
				$out_file = $head . ".fa";
				$out_file =~ tr/[\/><|:&\\,;?*]/_/;
				open (OUT, '>', "$out_dir/$out_file") or die "ERROR\tcould not open $out_dir/$out_file\n";
				push(@out_fa,"$out_dir/$out_file");
			}
		}
		if(exists($target_ids{$head})){
			print OUT $line . "\n";
		}
	}
	
	close OUT;
	close FA;
}

$length = length($fa_progress . " / " . $num_target_ids . " [" . $fa_frac . "%]");
print STDERR " " x $length . "\r";
print STDERR "INFO\tPredicting genestructure for $num_buscos BUSCOs and converting to zff format\n";

if(not -l "$out_dir/ancestral_variants"){
	$cmd = "ln -s $lineage_path/ancestral_variants $out_dir/ancestral_variants";
	exe_cmd($cmd,$verbose,$dry);
}
if(not -f "$out_dir/ancestral_variants.pog" or not -f "$out_dir/ancestral_variants.psd" or not -f "$out_dir/ancestral_variants.psi" or not -f "$out_dir/ancestral_variants.phr"     or not -f "$out_dir/ancestral_variants.psq" or not -f "$out_dir/ancestral_variants.pin"){
	$cmd = "makeblastdb -in $out_dir/ancestral_variants -dbtype prot -parse_seqids -out $out_dir/ancestral_variants > $out_dir/makeblastdb.log 2> $out_dir/makeblastdb.err";
	exe_cmd($cmd,$verbose,$dry);
}

my @cegma_gffs = ();
my @buscos_keys = keys(%buscos);

my $b_progress = 0;
my $b_frac = 0;
my $parent_pid = $$;
if(-f "$out_dir/.b2cp$parent_pid"){
	system("rm $out_dir/.b2cp$parent_pid");
}

my $multiple_buscos = Parallel::Loops->new($threads);
$multiple_buscos->share(\@cegma_gffs);
$multiple_buscos->share(\%genome_ann);

$multiple_buscos->foreach( \@buscos_keys, sub {
	if($verbose == 0){
		system("echo 1 >> $out_dir/.b2cp$parent_pid");
		$b_progress = `grep -c "1" $out_dir/.b2cp$parent_pid`;
		chomp $b_progress;
		$b_frac = sprintf "%.1f", ($b_progress / $num_buscos) * 100;
		print STDERR $b_progress . " / " . $num_buscos . " [" . $b_frac . "%]\r";
	}
	
	my $busco_id = $_;
	my @sites = split(/\n/,$buscos{$_});
	for(my $i = 0; $i < scalar(@sites); $i++){
		my ($status,$fa_id,$start,$end) = split(/\t/,$sites[$i]);
		my $fasta_in = $fa_id . ".fa";
		$fasta_in =~ tr/[\/><|:&\\,;?*]/_/;
		my $cmd = "augustus --noInFrameStop=true --species=$species --gff3=on --predictionStart=$start --predictionEnd=$end $out_dir/$fasta_in";
		if($status eq "Complete" or $status eq "Duplicated"){
			$cmd = $cmd . " --genemodel=complete";
			if($status eq "Duplicated"){
				$busco_id = (split(/_/,$busco_id))[0];
				$busco_id = $busco_id . "_" . $i;
			}
		}
		$cmd = $cmd . " > $out_dir/$busco_id.gff 2> $out_dir/$busco_id.err";
		exe_cmd($cmd,$verbose,$dry);
		
		if($dry == 0){
			open (GFF, '<', "$out_dir/$busco_id.gff") or die "ERROR\tCould not open $out_dir/$busco_id.gff\n";
			
			my %cds;
			my %strand;
			my %faa;
			my $faa_swithch = 0;
			my $gene_id = "";
			
			while (my $line = <GFF>){
				chomp $line;
				if($line =~ m/^# protein sequence = \[/){
					$faa_swithch = 1;
					my $aa = $line;
					$aa =~ s/^# protein sequence = \[//;
					$faa{$busco_id . "_" . $gene_id} = $aa;
				}
				else{
					if($faa_swithch == 1){
						my $aa = $line;
						$aa =~ s/^# //;
						if($aa =~ m/\]$/){
							$faa_swithch = 0;
							$aa =~ s/\]$//;
						}
						$faa{$busco_id . "_" . $gene_id} = $faa{$busco_id . "_" . $gene_id} . $aa;
					}
					else{
						if($line =~ m/^#/){
							next;
						}
						else{
							my @augustus = split(/\t/,$line);
							if($augustus[2] eq "CDS"){
								$gene_id = (split(/;/,$augustus[8]))[0];
								$gene_id =~ s/^ID=//;
								if(exists($cds{$gene_id})){
									$cds{$gene_id} = $cds{$gene_id} . "\n" . $line;
								}
								else{
									$cds{$gene_id} = $line;
								}
								if(not exists($strand{$gene_id})){
									$strand{$gene_id} = $augustus[6];
								}
							}
						}
					}
				}
			}
			
			close GFF;
			
			my @cds_keys = keys(%cds);
			
			if(scalar(keys(%faa)) > 1){
				@cds_keys = ();
				open (FAA, '>', "$out_dir/$busco_id.faa") or die "ERROR\tCould not open $out_dir/$busco_id.faa\n";
				foreach(keys(%faa)){
					print FAA ">" . $_ . "\n";
					print FAA $faa{$_} . "\n";
				}
				close FAA;
				$cmd = "blastp -db $out_dir/ancestral_variants -query $out_dir/$busco_id.faa -outfmt 6 > $out_dir/$busco_id.blastp 2> $out_dir/$busco_id.blastp.err";
				exe_cmd($cmd,$verbose,$dry);
				
				my %bit;
				my %target;
				open (BLAST, '<', "$out_dir/$busco_id.blastp") or die "ERROR\tCould not open $out_dir/$busco_id.blastp\n";
				while (my $line = <BLAST>){
					chomp $line;
					my @blast = split(/\t/,$line);
					if(exists($bit{$blast[0]})){
						if($blast[11] > $bit{$blast[0]}){
							$bit{$blast[0]} = $blast[11];
							$target{$blast[0]} = $blast[1];
						}
					}
					else{
						$bit{$blast[0]} = $blast[11];
						$target{$blast[0]} = $blast[1];
					}
				}
				close BLAST;
				
				foreach(keys(%target)){
					my $qid = (split(/_/,$_))[0];
					my $gid = (split(/_/,$_))[-1];
					my $tid = (split(/_/,$target{$_}))[0];
					if($qid eq $tid){
						push(@cds_keys,$gid);
					}
				}
				
				if($keep_tmp == 0){
					$cmd = "rm $out_dir/$busco_id.faa $out_dir/$busco_id.blastp $out_dir/$busco_id.blastp.err";
					exe_cmd($cmd,$verbose,$dry);
				}
			}
			
			for(my $g = 0; $g < scalar(@cds_keys); $g++){
				my @gene = split(/\n/,$cds{$cds_keys[$g]});
				
				if($strand{$cds_keys[$g]} eq "-"){
					my @reorder = ();
					for (my $i = scalar(@gene)-1; $i > -1; $i--){
						push(@reorder,$gene[$i]);
					}
					@gene = @reorder;
				}
				
				for(my $i = 0; $i < scalar(@gene); $i++){
					my @arr = split(/\t/,$gene[$i]);
					my @zff = ();
					$zff[3] = $busco_id;
					if(scalar(@cds_keys) > 1){
						$zff[3] = $zff[3] . "_" . $cds_keys[$g];
					}
					if($strand{$cds_keys[$g]} eq "-"){
						$zff[1] = $arr[4];
						$zff[2] = $arr[3];
					}
					else{
						$zff[1] = $arr[3];
						$zff[2] = $arr[4];
					}
					
					if(scalar(@gene) == 1){
						$zff[0] = "Esngl";
					}
					else{
						if($i == 0){
							$zff[0] = "Einit";
						}
						if($i > 0 and $i < scalar(@gene)-1){
							$zff[0] = "Exon";
						}
						if($i == scalar(@gene)-1){
							$zff[0] = "Eterm";
						}
					}
					
					if(exists($genome_ann{$fa_id})){
						$genome_ann{$fa_id} = $genome_ann{$fa_id} . "\n" . join("\t",@zff);
					}
					else{
						$genome_ann{$fa_id} = join("\t",@zff);
					}
				}
			}
				
			if($keep_tmp == 0){
				$cmd = "rm $out_dir/$busco_id.gff $out_dir/$busco_id.err";
				exe_cmd($cmd,$verbose,$dry);
			}
		}
		else{
			print STDERR "INFO\tI would open $out_dir/$busco_id.gff and convert to zff format\n";
		}
	}
});

$length = length($b_progress . " / " . $num_buscos . " [" . $b_frac . "%]");
print STDERR " " x $length . "\r";

if($keep_tmp == 0){
	my @db = ("ancestral_variants","ancestral_variants.pog","ancestral_variants.psd","ancestral_variants.psi","ancestral_variants.phr","ancestral_variants.psq","ancestral_variants.pin");
	foreach(@db){
		if(-f "$out_dir/$_"){
			$cmd = "rm $out_dir/$_";
			exe_cmd($cmd,$verbose,$dry);
		}
	}
}

if($dry == 0){
	print STDERR "INFO\tWriting zff files to $out_dir/genome.ann and $out_dir/genome.dna\n";
	
	my @genome_ann_keys = sort(keys(%genome_ann));
	if(-f "$out_dir/genome.dna"){
		$cmd = "rm $out_dir/genome.dna";
		exe_cmd($cmd,$verbose,$dry);
	}
	if(-f "$out_dir/genome.ann"){
		$cmd = "rm $out_dir/genome.ann";
		exe_cmd($cmd,$verbose,$dry);
	}
	open (GA, '>', "$out_dir/genome.ann") or die "ERROR\tCould not open $out_dir/genome.ann\n";
	for(my $i = 0; $i < scalar(@genome_ann_keys); $i++){
		print GA ">" . $genome_ann_keys[$i] . "\n";
		print GA $genome_ann{$genome_ann_keys[$i]} . "\n";
		my $fa = $genome_ann_keys[$i] . ".fa";
		$fa =~ tr/[\/><|:&\\,;?*]/_/;
		system("cat $out_dir/$fa >> $out_dir/genome.dna");
	}
	close GA;
}
else{
	print STDERR "INFO\tI would write zff files to $out_dir/genome.ann and $out_dir/genome.dna\n";
}

if($keep_tmp == 0){
	foreach(@out_fa){
		$cmd = "rm $_";
		exe_cmd($cmd,$verbose,$dry);
	}
}

if(defined(can_run("fathom")) and defined(can_run("forge")) and defined(can_run("hmm-assembler.pl"))){
	print STDERR "INFO\tCreating snap model $out_dir/$prefix.snap.hmm\n";
	if($verbose == 1){
		print STDERR "CMD\tcd $out_dir\n";
	}
	chdir "$out_dir";
	$cmd = "fathom genome.ann genome.dna -categorize 1000 > fathom.log 2> fathom.err";
	exe_cmd($cmd,$verbose,$dry);
	$cmd = "fathom -export 1000 -plus uni.ann uni.dna >> fathom.log 2>> fathom.err";
	exe_cmd($cmd,$verbose,$dry);
	$cmd = "forge export.ann export.dna > forge.log 2> forge.err";
	exe_cmd($cmd,$verbose,$dry);
	$cmd = "hmm-assembler.pl $prefix . > $prefix.snap.hmm 2> hmm-assembler.err";
	exe_cmd($cmd,$verbose,$dry);
	if($keep_tmp == 0){
		my @patterns = ("*count","*model","*duration","transitions","phaseprefs","alt.ann","alt.dna","err.ann","err.dna","export.aa","export.ann","export.dna","export.tx","olp.ann","olp.dna","uni.ann","uni.dna","wrn.ann","wrn.dna");
		$cmd = "rm \$(find -mindepth 1 -maxdepth 1 -type f -name \"" . join("\" -or -name \"",@patterns) . "\")";
		exe_cmd($cmd,$verbose,$dry);
	}
}

if(-f "$out_dir/.b2cp$parent_pid"){
	system("rm $out_dir/.b2cp$parent_pid");
}

exit;
