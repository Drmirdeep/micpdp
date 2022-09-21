#!/usr/bin/perl -w

use strict;
use Getopt::Std;
use File::Basename;
use Storable;

use Cwd 'abs_path';

my $script_path;

BEGIN {
	$script_path = dirname(abs_path($0));
}

use lib $script_path;
use My_micPDP qw(&get_codon_pos_index dnds);

my $version='0.99.11';

my $usage= qq{

micPDP v$version 

\nUsage: $0 -m multiple_alignment_file -s species_list -l min_orf_len -c chromosome -r annotation_files

[options]:
    m     multiple species alignment file in fasta format
    r     annotation files in ucsc table format
    n     non-coding gene annotations
    c     chromosome that is currently processed
    l     mininum orf length; default 27 nt = 9 AA
    s     list of species in maf according to phylogenetic tree
    t     output dir in subfolder projects
    S     ignore strand when checking for overlap with annotation files
    b     bed file with 12 coloumns so a maf file is not needed anymore, option -B and -b are mutually exclusive
 

[additional options]:
    F     mount point to folder 'ucsc' at which all maf files and indices are located
    E     exon mode - supply this option if you want only look for conservation in your exons bed file supplied by option -b
    M     maximum allowed length of an orf; default 300 nt = 100 AA
    W     set window size which is used when calculating entropies in neighboring regions; default ORF length
    A     output all orfs found and abort
    R     skip reapeat masked regions, only works if lower case letters in MAF file are repeats
    N     File with snps if these files are present
    P     path to reference sequence edits so we can filter out species that have frameshift mutations
    G     if given then gaps in the beginning and end of an alignment are corrected, default off since buggy still
    T     sequence is allowed to have this fraction of gaps, default is 0
    B     bed file with transcripts informing us about codon start positions --- must be used with option Y, otherwise there will be no effect !
    K     circular mode on genome, needs a 12 column bed file where CDS (c6 and c7) correspond to the circle in question 
    k     print full exon alignment , only useful for circmode !
    Y     for Ka/Ks calc we ignore codons with gaps, more fair I think - On by default
    y     we return 1 for each codon substitution and not the hamming distance
    X     score also with phyloCSF
    x     mincodon number to have for scoring with phylocsf
    v     verbose mode 
    F     disable flank scoring
    q     look also for ORFs on antisense strand

	to get all orfs use -A -S

};

die $usage if(not $ARGV[0]);

my %options=();
getopts("F:c:m:i:j:a:s:z:l:r:EM:W:n:o:I:t:SRO:AP:GT:b:Kk:YyB:X:x:vFDC:Q:Jq",\%options);

checkopts(%options);

## change this value to whereever BIO2 is mounted
## can be done on command line by specifying option 'F'
my $ucsc_mount_point='/scratch/local/ucsc/';
if($options{'F'}){
	$ucsc_mount_point=$options{'F'};
}


my %ii;
my @l;
my @iib;

if($options{'P'}){
open IN,"$options{'P'}/$options{'c'}.maf.index.sorted.i2" or die "Index file $options{'P'}/$options{'c'}.maf.index.sorted.i2 for index not given\n";

while(<IN>){
	chomp;
	if(/\#/){
		$ii{0}=1;
		push(@iib,0);
		next;
	}
	@l=split();
	$ii{$l[0]}=$l[1];
	push(@iib,$l[0]);
}
close IN;
}
# open my $iif,"$options{'P'}/$options{'c'}.maf.index.sorted" or die "index file could not be opened\n";

# my $res=binarysearch2(\@iib,$pos_to_go);
# my $hi=$iib[$res];
# my $pos=$ii{$hi};
# seek $iif,$pos,0;


print STDERR "
#########################################################################
# micPep - microPeptide detection tool 
# Version: $version
# Date: 16-12-2019
# Author: Sebastian Mackowiak 
# Email: sd.mackowiak\@gmail.com                            
#                                                                       
#########################################################################
";


## phylocsf parameters
my $currentgene;
my $phylo="$ucsc_mount_point/PhyloCSF";
my $param_phylo;
if($options{'s'} =~ /mm10/){
	$param_phylo="mouse60";
}elsif($options{'s'} =~ /mm9/){
	$param_phylo="mouse29";
}elsif($options{'s'} =~ /ce6/){
	$param_phylo="elegans";
}elsif($options{'s'} =~ /dr7/){
	$param_phylo="zebrafish";
}elsif($options{'s'} =~ /dm3/){
	$param_phylo="15flies";
}elsif($options{'s'} =~ /mm10/){
	$param_phylo="mouse60";
}else{
	$param_phylo="46vertebrates";
}

#$options{'Y'} =1;
#$options{'y'} =1;
#die $usage if(not $options{'m'} and not $options{'b'});


use POSIX qw/strftime/;
print STDERR "#\n#start\t", strftime("%d-%m-%Y\t%T",localtime),"\n#\n###########################################\n\n";


## fraction of gaps allowed in sequence
my $gapthreshold=0;
$gapthreshold=$options{'T'} if($options{'T'});
my $mincodons=8;
$mincodons=$options{'x'} if($options{'x'});

## get information about CDS files if options{'C'} is given
my %cdsinfo=();
my $id;
if($options{'Q'}){ ## change this back to C when bugfixed
	if(-f $options{'C'}){
	open IN,$options{'C'} or die "File $options{'C'} could not be opened\n";
	while(<IN>){
		chomp;
		@l=split();
		next if($l[0] ne $options{'c'});
		## this pattern match is for ensemble genes only, adapt when using refseq or sth else
		$id=$l[1];
		$cdsinfo{$id}{'start'}=$1 if($l[2] =~ /(\d)/);
		$cdsinfo{$id}{'end'}=$1 if($l[3] =~ /(\d)/);
		$cdsinfo{$id}{'len'}=$1 if($l[4] =~ /(\d)/);
		$cdsinfo{$id}{'offset'}=$1 if($l[5] =~ /(\d)/);
		$cdsinfo{$id}{'startX'}=$1 if($l[6] =~ /(\d)/);
		$cdsinfo{$id}{'startshift'}=$1 if($l[7] =~ /(\d)/);
		$cdsinfo{$id}{'endshift'}=$1 if($l[8] =~ /(\d)/);
		$cdsinfo{$id}{'no-cdsf'}=$1 if($l[9] =~ /(\d)/);
		$cdsinfo{$id}{'bonafide'}=$1 if($l[10] =~ /(\d)/);
	}
	close IN;
}
}

my $kaks_table=retrieve("$script_path/kaks_synmut.perlhash_storable");
my $kaks=$kaks_table;

## read in order of species in MAF
my @species;
lines_to_array($options{'s'},\@species);


my ($chr,$start,$len,$strand,$ms); 


$ms=$species[0]; ## this is our reference species


$chr=$options{'c'};


my @selected;

my $hg19="$ucsc_mount_point/hg19/maf/";
my $hg19e=".maf.stitched.cmpl.repeats_lc";

my $ce6="$ucsc_mount_point/ucsc/celegans/stitched_mafs/";
my $ce6e=".stitched.cmpl.repeats_lc";

my $ce10="$ucsc_mount_point/ucsc/celegans/ws220_stitched_mafs/";
my $ce10e=".maf.stitched.cmpl.r";

my $dm3="$ucsc_mount_point/ucsc/dm3/stitched_mafs/";
my $dm3e=".stiched.cmpl.repeats_lc";

my $dr7="$ucsc_mount_point/ucsc/zebrafish/maf/splitup_blocks/";
my $dr7e=".maf.stitched.cmpl.repeats_lc";

my $mm9="$ucsc_mount_point";
my $mm9e=".maf.stitched.cmpl.repeats_lc.r";

my $mm10="$ucsc_mount_point/ucsc/mm10/mafs_ucsc/";
my $mm10e=".maf.stitched.cmpl.repeats_lc.r";


my ($msb,$mse,$sf,$sf2);

#if($ARGV[5]){$msb="$ARGV[5]/"}else{
if($ms eq 'hg19' or $ms eq 'hsa'){
  $msb=$hg19;
}elsif($ms eq 'ce6' or $ms eq 'cel'){
  $msb=$ce6;
}elsif($ms eq 'ce10'){
  $msb=$ce10;
}elsif($ms eq 'dm3' or $ms eq 'dme'){
  $msb=$dm3;
}elsif($ms eq 'dr7' or $ms eq 'danRer7'){
  $msb=$dr7;
}elsif($ms eq 'mm9'){
  $msb=$mm9;
}elsif($ms eq 'mm10'){
  $msb=$mm10;
}else{
  die "species not found\n";
}



if ($ms eq 'hg19') {
	$mse=$hg19e;
	$sf="$msb/hg19_species";
} elsif ($ms eq 'ce6') {
	$mse=$ce6e;
	$sf="$msb/ce6_species";
} elsif ($ms eq 'ce10') {
	$mse=$ce10e;
	$sf="$msb/ce10_species_ucsc_order";
} elsif ($ms eq 'dm3') {
	$mse=$dm3e;
	$sf="$msb/insect_species";
} elsif ($ms eq 'dr7' or $ms eq 'danRer7') {
	$mse=$dr7e;
	$sf="$msb/dr7_species";
} elsif ($ms eq 'mm9') {
	$mse=$mm9e;
	$sf="$msb/mm9_spec.reordered";
} elsif ($ms eq 'mm10') {
	$mse=$mm10e;
	$sf="$msb/mm10_species.reordered";
} else {
	die "Species $ms is not valid\n";
}

if(-f $options{'m'}){
	open MAF,'<',$options{'m'};
}else{
	if(not -f "$msb$chr$mse"){ print STDERR "file $msb$chr$mse maf file not found\nskipping\n";exit;}
	open MAF,"$msb$chr$mse" or die "maf file $msb$chr$mse not found\n";
}

my $fh=*MAF;

my ($refl,$end,$offset,$lenseq);

my $dmp;

while (<MAF>) {
	if (/>[a-zA-Z0-9]+.$chr+\(\S\):(\d+)-(\d+)/) {
		$offset=$1;
		$end=$2;
		$lenseq=$end-$offset;
		$refl=length($_);
		last;
	} else {
		die "pattern not matched $_\n";
	}
}

my $length_chr=$end;

my %spec=();

open IN,$sf or die "species file $sf not found\n";

my @dmp;
my $mult=-1;
my $specl;
my $seql=0;
my $first=1;
my $s='';

my @species2;
my $origs;
my $real;

while (<IN>) {
	if (/(\S+)/) {
		if($ms =~ /mm\d/){
			my $d=$_;
			chomp $d;
			if($d =~ /^(\S+)\s*(\S*)/){
				$origs=$1;
				$real=$2;
			}
			if ($ms ne 'mm10') { ## so it is mm9
				$d=~ s/\s+/:/g;
			} else {
				if($real){
					$d.=":";
				}
			}
			$mult++;
			$s="$origs:$real";
			push(@species2,$origs);
			if ($first) {
				$first=0;
				$specl+=$refl;
			} else {
				$specl+=length($s)+2;
				$seql+=$lenseq+1;
			}
			$spec{$origs}=$specl+$seql;
		}else{
			@dmp=split(",",$1);
			foreach my $d(@dmp){
				$mult++;
				if(/^>*(\S+)/){
					$s=$1;
					push(@species2,$s);
					if($first){
						$first=0;
						$specl+=$refl;
					}else{
						$specl+=length($s)+2;
						$seql+=$lenseq+1;
					}
					
					$spec{$s}=$specl+$seql;
				}
			}
		}
	}
}
close IN;

my %gapindex;
my %gapindexr; ## this has position and then as second key the species to delete
my $href=\%gapindex;
## if we want to retrieve hash with gap information and so on

my $line;

my $itergap;
if ($options{'P'}) {
	#$href = retrieve("$options{'P'}/$options{'c'}.maf.perlhash_storable");
	open GI,"$options{'P'}/$options{'c'}.maf.index.sorted" or die "Index.sorted file not found\n";
	$itergap=filehandle_iterator(*GI);
	$line=$itergap->(); 
	$line=$itergap->() if($line =~ /\#/); ## only if we have a header line
	@l=split(/\s+/,$line);
}


if(not -d "projects"){
   mkdir "projects";
}


my $tag="default";
$tag=$options{'t'} if($options{'t'});
if (not -d "projects/$tag/") {
	mkdir "projects/$tag/";
}

my $window=30;
$window=$options{'W'} if($options{'W'});

my %tss; ## holds all tss
my %tse; ## holds all tse
my @tssap;
my @tesap;
my @tssam;
my @tesam;
my ($dtss,$dtes);

my %anno=();
if ($options{'a'}) {
	readin_annotation_files($options{'a'},\%anno);
}

my $min_orf_len=27;
$min_orf_len=$options{'l'} if($options{'l'});

my $max_size=300; ## max orf len
$max_size=$options{'M'} if $options{'M'};
$max_size++;

my $igs=0;
$igs=1 if($options{'S'});


my %ignores;
$ignores{'+'}='-';
$ignores{'-'}='+';

## read protein alphabet to hashes codons and rcodons
my %codons;
my %rcodons;
initProt(\%codons,\%rcodons);


my %refseq;
my %nc_refseq=();
my %oa_files=();
my ($pi,$mi);

my %fhh=();

my @files;



## read in file with transcript info
# my ($fhb,%codonpos);
my %cp;
if($options{'B'}){
 	open my $optionsB ,$options{'B'} or die "File $options{'B'} not found\n";
	get_codon_pos_index($options{'B'},\%cp,$options{'c'});
	
	print STDERR "codon pos read in for all transcripts on $options{'c'}\n";
	my @j=sort {$b <=> $a} keys %cp;
#	die join(",",@j),"\n";

	## this can be optimized by using iterators on the file

# 	$fhb=filehandle_iterator(*$optionsB);
# 	$fhl=$fhb->(); #iteration1
# 	die $fhl;
}

if ($options{'r'}) {
	@files=split(",",$options{'r'}) or die "$options{'r'} not given\n";
	foreach my $f (@files) {
#		if ($options{'Z'}) {
			$f.=".sorted";
#		}

		open my $fh ,$f or die "File $f not found\n";

		$fhh{$f}{'type'} = 'c';
		$fhh{$f}{'fh'} = *$fh;
		$fhh{$f}{'it'} =  filehandle_iterator(*$fh);
		$fhh{$f}{'l'} = $fhh{$f}{'it'} ->(); ## iteration 1
		if ($fhh{$f}{'l'} =~ /\#/) {
			readin_refseq(\%refseq,$fhh{$f}{'l'},$f);
			$fhh{$f}{'l'} = $fhh{$f}{'it'} ->(); ## iteration2
		}
		get_tss_tse(\%fhh,$f);
	}
}

if ($options{'n'}) {
  @files=split(",",$options{'n'}) or die "$options{'n'} not given\n";
  foreach my $f (@files) {
	  $f.=".sorted";

    open my $fh ,$f or die "File $f not found\n";
    $fhh{$f}{'type'} = 'n';
    $fhh{$f}{'fh'} = *$fh;
    $fhh{$f}{'it'} =  filehandle_iterator(*$fh);
    $fhh{$f}{'l'} = $fhh{$f}{'it'} ->(); ## iteration 1
    if ($fhh{$f}{'l'} =~ /\#/) {
	  readin_refseq(\%nc_refseq,$fhh{$f}{'l'},$f);
	  $fhh{$f}{'l'} = $fhh{$f}{'it'} ->();
    }
    get_tss_tse(\%fhh,$f);
  }
}


@tssap=sort {$a <=> $b} keys %{$tss{$chr}{'+'}};
@tssam=sort {$a <=> $b} keys %{$tss{$chr}{'-'}};
@tesap=sort {$a <=> $b} keys %{$tse{$chr}{'+'}};
@tesam=sort {$a <=> $b} keys %{$tse{$chr}{'-'}};


## read in maf with species alignments
my %maf=();
print STDERR "reading maf file now\n";
my $finalpos='';

my %positions=();
## if we use a bed file we ignore the given global maf file cause we extract everything directly from genome

my %ph=();

my @bl;
my @bs;

 ## open output files

my $what="genomic-strand\tmp-orientation";
my $whatc='#chr';
$what="tx-orientation_genomic\tmp-orientation" if($options{'b'});
$whatc='#transcript' if($options{'b'});
my $add='';

if($options{'K'}){
	$add.="\tcirc_len\tstart.frag\tend.frag\t#junctionc";
}else{
	$add="\ttype\tscoretxregion";
}


my %seen=();
if($options{'X'}){
	mkdir $options{'X'} if(not -d $options{'X'});
	$add.="\tomega\tomega_flankl\tomega_flankr\tomega_cds\tcodons-scored";
}

my $header="$whatc\t$what\tstart\tend\tlen\tframe\t#species\tsyn_mut\tAA\tDNA\tdnaH(orf)\taaH(orf)\taaH(pre)\taaH(suf)\tspecies\tkozak\tshine-d\tdistNTSS\tdistNTSE\tdistNstart\tdistNend\tother_annotations\twsum,codMut,specMut\tsyn_mut,nsyn_mut,uniq_syn_mut,uniq_nsyn_mut,ka,ks\tomega\t#startC:val\t#stopC:val\tucsc_coords\tmedian(dndsORF)\tmedian(dndsORF_iss)\tmedian(dndsleftF)\tmedian(dndsrightF)\tmean(dndsORF)\tmean(dndsleftF)\tmean(dndsrightF)\tsd(dndsORF)\tsd(dndsleftF)\tsd(dndsrightF)\t#frameshiftSpecies\tframeshiftSpecies$add";

$header.="\t#nonsense";
$header.="\n";

if(not $options{'b'}){
  open COD,">projects/$tag/${chr}_coding_overlap" or die "File ${chr}_coding_overlap could not be created\n";
  print COD $header;
  
  open NC,">projects/$tag/${chr}_noncoding_overlap" or die "File ${chr}_noncoding_overlap could not be created\n";
  print NC $header;
  
  open IG,">projects/$tag/${chr}_intergenic" or die "File ${chr}_intergenic could not be created\n";
  print IG $header;
  
  open ASN,">projects/$tag/${chr}_asn" or die "File ${chr}_asn could not be created\n";
  print ASN $header;
  
  open ASC,">projects/$tag/${chr}_asc" or die "File ${chr}_asc could not be created\n";
  print ASC $header;
  open ERR ,">projects/$tag/${chr}.err" or die "Could not open file projects/$tag/${chr}.err\n";

}else{

  ## we need to put the filename in here
	my $file;
	if($options{'K'}){
		if(not -d "$chr"){
			mkdir "projects/$tag/${chr}";
			my @sp=split("/",$options{'b'});
			
			$file=$sp[$#sp];
		}
		#open IG,">projects/$tag/${chr}/${file}_transcripts" or die "Could not open file projects/$tag/${chr}_transcripts\n";

	}else{
		open IG,">projects/$tag/${chr}_transcripts" or die "Could not open file projects/$tag/${chr}_transcripts\n";
	}
  print IG $header;
  open ERR ,">projects/$tag/${chr}_transcripts.err" or die "Could not open file projects/$tag/${chr}_transcripts.err\n";
}

if($options{'b'}){
	if(not $options{'K'}){
		open BEDI,">projects/$tag/${chr}_transcripts.bed" or die "Could not open Bed file for writing\n";
		#open BEDIE,">projects/$tag/${chr}_transcripts_exons.bed" or die "Could not open Bed file for writing\n";
	}
	my $bs="$options{'b'}.sorted";
	if(not -f $bs){
		print STDERR "File $bs does not exist, checking local dir\n";
		my ( $name2, $path2, $extension2 ) = fileparse ( $options{'b'}, '\..*' );
		$bs="$name2$extension2.sorted";
		if(not -f $bs){
			print STDERR "File $bs does not exist, trying to generate it in local dir\n";
			system("sort -nk2 $options{'b'} > $bs");
		}else{
			print STDERR "File $bs exists, using this one\n";
		}
	}

	open IN,"$bs" or die "Bed file $bs not found for processing\n";
	
	%seen=();
	my @line;
	my @iline;
	my $pold;
	my $tc=0;
	my $rel;
	my $curbp;
	my $type;
	my $i;
	my $s;
	my $il;
	my %tmaf;
	while(<IN>){
		next if(/^\s*$/);
		next if(/\#/);
		@iline=split();
		next if($iline[0] ne $chr);
		$tc++;
#		print STDERR "$tc\r";
		%maf=();
		%ph=();
		chomp;
		$il=join("\t",@iline);
		if(not defined $iline[6]){			
			$il.="\t$iline[2]\t$iline[2]\t0\t1\t".($iline[2]-$iline[1])."\t0";
		}
		@line=split(/\s+/,$il);

		print STDERR "$line[3]\n" if($options{'v'});

		$currentgene=$line[3];

		getseq_bed($fh,$offset,$chr,$il,\%maf,\@species);
		$finalpos=length($maf{$species[0]});
		
		$rel=-1; ## relative 0 based postion
		$curbp=0;
		$type='';
	
		$ph{'strand'} = $line[5];
		$ph{'name'} = $line[3];
        if($line[5] ne '+' and $line[5] ne '-'){
            $line[5] = '+';
            $ph{'name'}.='nostrand';
        }
		$ph{'txs'} = $line[1];
		$ph{'txe'} = $line[2];
		if(defined $line[6]){
			$ph{'cdss'} = $line[6];
			$ph{'cdse'} = $line[7];	
			$ph{'bs'} = $line[11];
			$ph{'bl'} = $line[10];
			$ph{'bc'} = $line[9];
			@bl=split(",",$line[10]);
			@bs=split(",",$line[11]);
		}else{
			$ph{'cdss'} = $line[1];
			$ph{'cdse'} = $line[1];	
			$ph{'bs'} = $line[1];
			$ph{'bl'} = $line[2]-$line[1];
			$ph{'bc'} = $line[1];
			@bl=split(",",$ph{'bl'});
			@bs=split(",",0);
		}

   		my %newindex=();
		## delete from gap index here  until we reach next tx start so we always follow this guys here, maybe we only check postions where we have indeed gaps ??
		## this may create huge computing overhead
		if($pold and $options{'P'}){
			if($pold != $ph{'txs'}){
			for ($i = $pold;$i < $ph{'txs'}; $i++){
				if($gapindexr{$i}){
					foreach $s(keys %{$gapindexr{$i}}){
						delete $gapindex{$s}{$i} if(exists $gapindex{$s}{$i});
					}
					delete $gapindexr{$i};
				}
			}
		}
			$pold=$ph{'txs'}; ## set new pold
		}


		## set the initial $pold variable
		if(not $pold){
			$pold = $ph{'txs'};
		}

		my $llen=0;
		for($i= 0;$i<scalar @bs;$i++){
			$llen+=$bl[$i];
		}
		my $pa=-1;
		$llen--;
		
		
		for($i= 0;$i<scalar @bs;$i++){
			$curbp=-1;
			$ph{'len'}+=$bl[$i];
			for($s=$bs[$i]; $s < $bs[$i]+$bl[$i];$s++){
				$curbp++;
				$pa++;
				if($ph{'strand'} eq '+'){
					$rel=$pa;
				}else{
					$rel=$llen-$pa;
				}
				
				$ph{$rel}{'abs'}=$line[1]+$s;  ## absolute genome position
				

				if($ph{'strand'} eq '-'){
					$ph{'cdsse'} = $rel if($ph{$rel}{'abs'} == $line[6]);
					$ph{'cdssr'} = $rel if($ph{$rel}{'abs'} == $line[7]-1);

				}else{
					$ph{'cdssr'} = $rel if($ph{$rel}{'abs'} == $line[6]);
					$ph{'cdsse'} = $rel if($ph{$rel}{'abs'} == $line[7]-1);
				}

				## where are we, utr or cds
				if($line[6] != $line[7]){
					if($ph{$rel}{'abs'} >= $line[6] and $ph{$rel}{'abs'} < $line[7]){$type='cds';} ## we also need a counter for the frames
					if($ph{$rel}{'abs'} < $line[6]){$type='utr5';$type='utr3' if($line[5] eq '-')}
					if($ph{$rel}{'abs'} >= $line[7]){$type='utr3';$type='utr5' if($line[5] eq '-')}
				}else{
					$type='nc';
				}
					
				$ph{$rel}{'t'}=$type;
				$ph{$rel}{'pos'}=$bs[$i];
				$ph{$rel}{'cpos'}=$curbp;
			}
		}

		if($options{'K'}){ ## so we use circular mode , then extract first the circ from the full bed and then concat it 4 times
			## if cds could not be identified cause start or end are not in extracted sequence then skip. This happens if one of them is not part of the exon structure given
			
			if(not defined $ph{'cdsse'} or not defined $ph{'cdssr'}){
				print STDERR "$ph{'name'} has inconsistent cds start or end and will be skipped\n";
				next;
			}

			$ph{'circl'} = $ph{'cdsse'}-$ph{'cdssr'}+1;
			if(not $ph{'circl'}){
				print STDERR "here circl $ph{'name'}\n";
				next;
			}

			foreach my $s(keys %maf){
				$tmaf{$s}=substr($maf{$s},$ph{'cdssr'},$ph{'circl'});
#				if($ph{'circl'} % 3 != 0){
					$tmaf{$s}=$tmaf{$s}.$tmaf{$s}.$tmaf{$s}.$tmaf{$s};   ## we put it 4 times behind each other since we need to have a circle 
#				}
			}

			%maf=%tmaf;
			$finalpos=length($maf{$species[0]});
			if($options{'k'}){
				if(not -d $options{'k'}){
					mkdir $options{'k'};
				}

				open OUT,">$options{'k'}/${line[3]}";
				foreach my $s(@species){
					print OUT ">$s";
					if($s eq $species[0]){
						if($ph{'strand'} eq '+'){
							print OUT "|$line[0]:$line[6]-$line[7]:$line[5]:0-$ph{'circl'}:0-",4*$ph{'circl'},":$ph{'cdssr'}-",$ph{'cdsse'}+1;
						}else{
							print OUT "|$line[0]:$line[6]-$line[7]:$line[5]:0-$ph{'circl'}:0-",4*$ph{'circl'},":",$ph{'cdssr'}+1,"-$ph{'cdsse'}";
						}
					}
						print OUT "\n$maf{$s}\n";
				}
				close OUT;
			}
		}


		## in case of a noncoding RNA we set this to 0 end length
		$ph{'cdssr'} = 0 if(not defined $ph{'cdssr'});   ### in case we have a noncoding rna in here
		$ph{'cdsse'} = $ph{'len'}-1 if(not defined $ph{'cdsse'});   ### in case we have a noncoding rna in here


		## if we check for coords now we first have the offset of $ph{'cdssr'} and then we need to keep in mind, that we can be in circle round 1 2 3 or 4 when checking for gaps and so on
		if($options{'K'}){
			$ph{'txs'} = $line[6];
			$ph{'txe'} = $line[7];
		}			

#		die $cp{$ph{1543}{'abs'}};

		## this seems to work, I need to track the coords for the circles if there is something in and we need to print the alignemnts in the end
		my $skipnt=0;
		if($options{'E'}){ ## we score only exons for conservation
			detect_orfs_exons(\%maf,\%positions,\@species,\%refseq,\%nc_refseq,\%oa_files,$igs,$window,$kaks,$finalpos,\%ph,$itergap,\%cp,*GI,$skipnt);
		}else{
			detect_orfs(\%maf,\%positions,\@species,\%refseq,\%nc_refseq,\%oa_files,$igs,$window,$kaks,$finalpos,\%ph,$itergap,\%cp,*GI,$skipnt);
		}
	}
	
	close IN;
	close BEDI;
	#close BEDIE;
	close IG;
  }else{
	my $skipnt;
	read_fasta($options{'m'},\%maf,\$finalpos,$species[0],\$skipnt); 
	print STDERR "Fasta file $options{'m'} read\n";
	## now do the real job and detect micropeptidegs
	detect_orfs(\%maf,\%positions,\@species,\%refseq,\%nc_refseq,\%oa_files,$igs,$window,$kaks,$finalpos,\%ph,$itergap,\%cp,*GI,$skipnt);
  }
close MAF;


print STDERR "#\n#end\t", strftime("%d-%m-%Y\t%T",localtime),"\n#\n###########################################\n\n";
exit;


## cat files together and make bed files and html files
my $oo;
if (not $options{'O'}){$oo = "-o ";}else{ $oo = "-o $options{'O'} ";}
print STDERR "now do
cd projects/$tag
cat chr*intergenic > all_orfs_intergenic
perl ../../score.pl -i all_orfs_intergenic -o $oo > all_orfs_intergenic.tmp
head -n1  all_orfs_intergenic.tmp >  all_orfs_intergenic.mssm
grep -v \\# all_orfs_intergenic.tmp |sort -rnk 2  >> all_orfs_intergenic.mssm
perl ../../make_bed_file.pl all_orfs_intergenic.mssm
";

print STDERR "Now you should upload your bedfiles to UCSC and then 
create an html file with the correct hgsid
#mirpep2table.pl s_6_m_6.mssm.cand dme 1 hgsid  30
../../mp_make_html_red.pl -t all_orfs_intergenic.mssm -s 3 -e 4 -c 1 -i 0 -S cel -w 30 -I hgsid`
";


#`../../mp_make_html_red.pl -t all_orfs_intergenic.mssm -s 3 -e 4 -c 1 -i 0 -S cel -w 30`

exit;

my @F;
for(my $i=1;$i<=12;$i++){
  open IIN,"all_orfs_intergenic_.mssm" or die "File all_orfs_intergenic.mssm not found\n";
  open OUTI,">s_${i}_m_${i}.mssm.cand" or die "File s_${i}_m_${i}.mssm.cand could not be created\n";
  while(<IIN>){
    @F=split();
    if(/\#/){print OUTI;}elsif($F[7]>$i && $F[8] > $i){print OUTI;}else{}
  }
  close OUTI;
  close IIN;
  `perl ../../mirpep_make_bed_file.pl s_${i}_m_${i}.mssm.cand $tag${i}${i}`;
}

print STDERR "Now you should upload your bedfiles to UCSC and then 
create an html file with the correct hgsid
../../mp_make_html_red.pl -t all_orfs_intergenic.mssm -s 3 -e 4 -c 1 -i 0 -S cel -w 30 -a 6 -b 6`
";

#######################################################
#######################################################
#######################################################
sub getseq_bed{
    ## fuction needs to get a filehandle
    ## must be called like check(*IN,1234)
    ## the star is really important
    ## this is random access to index file, saves a lot of time

    my ($fh,$offset,$chrr,$line,$maf,$selected)=@_;
	
	my ($tseq,$seq,$pos,$first);
    my ($chr,$start,$end,$name,$score,$strand,$dmp1,$dmp2,$color,$blockcount,$blocksizes,$blockstarts)=split(/\s+/,$line);
    $blocksizes=$1 if($blocksizes =~ /(\S+),$/);
    $blockstarts=$1 if($blockstarts =~ /(\S+),$/);

    my @bs=split(",",$blocksizes);
    my @bb=split(",",$blockstarts);


    ## get coords
    my @cc;
    my $coords='';
    my ($b,$e);
    for(my $i=0;$i<scalar @bs;$i++){
        $b=$start+$bb[$i];
        $e=$start+$bb[$i]+$bs[$i];
        $coords.="$b,$e,";
        push(@cc,$b);
        push(@cc,$e);
    }


    # my @cc=split(",",$coords);
    if($cc[0] == $cc[1]){
        print STDERR "$name\tlen0\n";
        next;
    }
    $first=0;
    foreach my $s(@$selected){
        $seq='';
        if(not $first){
            $first = 1;
            #print ">$s.$chr($strand):$coords:$name\n";
        }else{
            #print ">$s\n";
        }

        for (my $i=0;$i<scalar @cc;$i+=2){
            $len=$cc[$i+1]-$cc[$i]; ## get the length correctely
            $start=$cc[$i];
            $pos=$spec{$s}+$start-$offset;
            seek $fh,$pos,0;
            $tseq='';
            read($fh,$tseq,$len);
            $seq.=$tseq;
        }
        if($strand eq '-'){
            $seq=rc($seq);
		}
        #print $seq,"\n";
		$$maf{$s} = $seq;
    }
}


sub getseq{
    ## fuction needs to get a filehandle
    ## must be called like check(*IN,1234)
    ## the star is really important
    ## this is random access to index file, saves a lot of time

    my ($fh,$start,$len,$offset,$strand,$tmphash,$species)=@_;

    my ($seq,$pos,$first);
    foreach my $s (@$species) {

        $pos=$spec{$s}+$start-$offset;

        seek $fh,$pos,0;
        read($fh,$seq,$len);
        if ($strand eq '-') {
            $seq=rc($seq);
        }
		$$tmphash{$s} = $seq;
    }
}


## some iterator subs 
sub Iterator (&) { return $_[0];}

sub filehandle_iterator {
	my $fh = shift;
	return Iterator { <$fh> };
}

sub get_tss_tse{
    my ($fhh,$f)= @_;
	my ($cname,$cchr,$cstrand,$ctxs,$ctxe,$ccds,$ccde,$ces,$cee)=(-1,-1,-1,-1,-1,-1,-1,-1,-1);
    my @line;
    ($cchr,$cstrand,$ctxs,$ctxe,$ccds,$ccde,$cname)= ($fhh{$f}{'token'}{'chr'},$fhh{$f}{'token'}{'strand'}, $fhh{$f}{'token'}{'txstart'},$fhh{$f}{'token'}{'txend'},$fhh{$f}{'token'}{'cdsstart'},$fhh{$f}{'token'}{'cdsen\
d'},  $fhh{$f}{'token'}{'name'});
    open IN,$f or die "File $f could not be opened\n";
    while(<IN>){
        next if(/\#/);
        @line=split(/\s+/);
        $tss{$line[$cchr]}{$line[$cstrand]}{$line[$ctxs]}=1;
        $tse{$line[$cchr]}{$line[$cstrand]}{$line[$ctxe]}=1;
    }
    close IN;
}



sub check{
    my ($i,$h,$line,$f,$maxlen,$fhh,$g2,$chr,$startpos) = @_;
    return if(not defined $$fhh{$f}{'l'});
	return if($startpos > $i+$maxlen+1);

	#die $i if($f eq 'ce6_wormbase_RNAs_12_02_13' and $i == 895);


    my ($utr,$cds,%ee,%es,$ie,$id,$pseudo,$extra)=();
    ## adapt this here
    delete $$h{$chr}{'+'}{($i-3*$maxlen)} if (exists $$h{$chr}{'+'}{($i-3*$maxlen)}); ## delete previous positions if existed in hash
    delete $$h{$chr}{'-'}{($i-3*$maxlen)} if (exists $$h{$chr}{'-'}{($i-3*$maxlen)}); ## delete previous positions if existed in hash

    ## get variables here for ctxstart,ctxend,cname and chromosome
	my ($cname,$cchr,$cstrand,$ctxs,$ctxe,$ccds,$ccde,$ces,$cee)=(-1,-1,-1,-1,-1,-1,-1,-1,-1);
    ($cchr,$cstrand,$ctxs,$ctxe,$ccds,$ccde,$cname,$pseudo)= ($$fhh{$f}{'token'}{'chr'},$$fhh{$f}{'token'}{'strand'}, $$fhh{$f}{'token'}{'txstart'},$$fhh{$f}{'token'}{'txend'},$$fhh{$f}{'token'}{'cdsstart'},$$fhh{$f}{'token'}{'cdsend'},  $$fhh{$f}{'token'}{'name'},$$fhh{$f}{'token'}{'pseudo'});

	if (not $ccde) {
		die $f;
	}
	#	die "here" if($f eq 'ce6_wormbase_RNAs_12_02_13' and $i == 895);

    ## goto next line if current does not match chromosome
    while ($$line[$cchr] ne $chr) {
        $$fhh{$f}{'l'} = $fhh{$f}{'it'} ->();
        last if(not defined $$fhh{$f}{'l'});
		@$line=split(/\s+/,$fhh{$f}{'l'}); ## write to same reference
    }
    return if(not defined $$fhh{$f}{'l'}); ## exit from here

	## read in until we reach this limit, eh?
	while ($$line[$ctxs] < $i+$max_size+1) {
		
		#	print "$i $$line[$ctxs] $$line[$cname] $f \n";

		#		die "$$line[$ctxs]  $f here \n" if($i == 11315 and $$line[$ctxs] == 11315);

		
		#		if($$line[$ctxs] > $i){
		#			return;  ## in this case we wait until i == $l[4]
		#		}elsif($$line[$ctxs] == $i){
		#		print "inserting to hash $i $$line[$cname] \t  $f \n";# if($i == 11315);
		#now we read in the line and ask for next line until l[4] > $i;
		## read in until we have all coords in;
		
		($utr,$cds,%ee,%es,$ie,$id)=();
		$id=$$line[$cname];
		## skip if not correct chromosome and get next line
			
		## newly added 30.09.13
		if (exists $anno{$id}) {
			$extra="$anno{$id}{'type'}_";
			#print $extra,"\n";
		} else {
			$extra='';
		}
		
		## for each position between txStart and txEnd
		if ($ces != -1) {
			@es{split(",",$$line[$ces])}=undef; ## put array as keys in hash
			@ee{split(",",$$line[$cee])}=undef;
		}
		
		my ($utr,$cds,$exon)=("","","");
		for (my $i=$$line[$ctxs];$i<=$$line[$ctxe];$i++) {
			if (exists $es{$i}) {
				$exon=",exon;";
			} elsif (exists $ee{$i+1}) { ## bed files have end coord not included
				$exon=",intron;";
			} else {
				$exon='';
			}
				
			
			## add utr info if necessary  ## would also need to check if in utr exon or not ... -> better check spliced 5utrs alone -> did this already ...
			if ($ccds != -1) {
				if ($$line[$cstrand] eq '+') {
					if ($i < $$line[$ccds]) {
						$utr=',utr5';
					} elsif ($i > $$line[$ccde]) {
						$utr=',utr3';
					} else {
						$utr='';
					}
				} else {
					if ($i < $$line[$ccds]) {
						$utr=',utr3';
					} elsif ($i > $$line[$ccde]) {
						$utr=',utr5';
					} else {
						$utr='';
					}
				}
				
					
				if ($i >= $$line[$ccds] and $i < $$line[$ccde]) {
					$cds=',cds';
				} else {
					$cds = '';
				}
				## if noncoding then no utr present and not cds present but exons may be present
				if ($$line[$ccds] == $$line[$ccde]) {
					$cds='';
					$utr='';
				}
			}
			if ($$line[$cname] =~ /NR_/ and $g2 != 0) {
				$$g2{$$line[$cchr]}{$$line[$cstrand]}{$i} .= "$extra$pseudo$$line[$cname]$utr$cds$exon";
			} else {			## this is for protein coding regions
				$$h{$$line[$cchr]}{$$line[$cstrand]}{$i} .= "$extra$pseudo$$line[$cname]$utr$cds$exon";
			}
		}
		  
		$$fhh{$f}{'l'} = $fhh{$f}{'it'} ->();
		last if(not defined $$fhh{$f}{'l'});
		@$line=split(/\s+/,$fhh{$f}{'l'}); ## write to same reference
	}
}
			

sub check_gaps{
	my ($href,$start,$end,$hs,$species,$strand)=@_;
	my ($keep,$s);
	my $i=0;
	#print "$start\t@$species\n";
	my $beg=$start;
	$beg++; ## since ucsc index starts at one and we at 0 correct
	#$beg++ if($strand eq '-');
	
	for(my $pos=$beg;$pos<$end+1;$pos++){
		#print "$species[0]\t$pos\t$$href{$species[0]}{$pos}{'gaps'}\n";
		if($$href{$$species[0]}{$pos}{'gaps'}){ ## we found gaps in alignmet, (insertions in other species), now check which and if we need to kick them out ...
			foreach $s(@$species){
				next if($s eq $$species[0]); ## next if refspecies
				$keep=1;
				if($$href{$s}{$pos}{'ins'}){ ## ok we have insertions in the species
					#if($$href{$s}{$pos}{'ins'} %3 ==0 and $i%3==0){ ## so if multiple of 3 for insertions and we are in frame 0 for ref then insertions only add amino acids and can be kept
				#	print "$s\t$pos\t$$href{$s}{$pos}{'ins'}\n";
					if(($$href{$s}{$pos}{'ins'} %3 ==0) and ($$href{$s}{$pos}{'ins'} < ($end-$start))){ ## must simply be a multiple of 3, so no frameshift mutation and length of insert is less than orflength, should maybe just half of orflength
						$keep = 1;
			#			print "$pos $s keep $$href{$s}{$pos}{'ins'}\n" if($start == 4220);
					}else{
						$keep = 0;
             #           print "$pos $s kickout $$href{$s}{$pos}{'ins'}\n" if($start == 4220);
						## print STDERR "key $s must be deleted due to $$href{$s}{$pos}{'ins'} frameshift insertions at pos $pos\n";
						}
				}
				if($keep == 0){
					if(exists $$hs{$s}){
						delete $$hs{$s};
							## print STDERR "key $s delete from species list due to frameshift insertions\n";
					}
				}
			}
		}
		$i++;
	}
	#print "$end\t@$species\n";
}

sub check_gaps_tx{
	my ($href,$start,$end,$hs,$species,$strand,$ph)=@_;
	my ($keep,$s);
	my $i=0;
	my $frame=0;
	my $pos; ## we changed pos to rpos and get for pos the genomic pos from utr hash, thats it finish
## sth wrong here

	my $rels;
	my $rele;
	my $ostr='+';
	my $as='';
	$pos=$start;
	
	#my $id2="\t$options{'c'}:".($$ph{$rels}{'abs'})."-".($$ph{$rele}{'abs'}+1).":$ostr";
	#$dndsoutstr.=$id2;
	if($$ph{'strand'} eq '+'){
		$rels=$start;
		$rele=$end;
	}else{
		$rels=$$ph{'len'}-$end-1; ## is correct
		$rele=$$ph{'len'}-$start-1; ## is correct
	}
	
	for(my $rpos=$rels;$rpos<$rele;$rpos++){
		$i++; ## tells us the relative pos in ORF, important to get correct frame later
		## now we recalc the relative coords to genomic coords
		$pos=$$ph{$rpos}{'abs'};
		

		$pos++; # our gap index is one based but we are 0 based in utrs, fault of the u if (necessary?)
		#print "$strand\t$i\t$pos\t$rpos\t";
		$pos-- if($strand eq '-' and $$ph{'strand'} eq '+'); ## if we are on minus strand then gap start is on pos downstream
		$pos-- if($strand eq '+' and $$ph{'strand'} eq '-'); ## if we are on minus strand then gap start is on pos downstream
		
		if($$href{$$species[0]}{$pos}{'gaps'}){ ## we found gaps in alignmet, (insertions in other species), now check which and if we need to kick them out ...
			foreach $s(@$species){
				next if($s eq $$species[0]); ## next if refspecies
				$keep=1;
				if($$href{$s}{$pos}{'ins'}){ ## ok we have insertions in the species
					#print "$s\t$pos\t$$href{$s}{$pos}{'ins'}\n";
					if(($$href{$s}{$pos}{'ins'} %3 ==0) and ($$href{$s}{$pos}{'ins'} < $end-$start)){ ## length of insert < orflength and multiple of 3
						$keep =1;
					}else{
				#		print STDERR "$start  $s  $pos\n"; ## delme
						$keep=0;
						}
				}
				if($keep ==0){
					
					if(exists $$hs{$s}){
						delete $$hs{$s};
					}
				}
			}
		}
	}
}



sub check_gaps_tx2{
	my ($href,$start,$end,$hs,$species,$strand,$ph)=@_;
	my ($keep,$s);
	my $i=0;
	my $frame=0;
	my $pos;
# 	print "====== $start\t$end\t$$ph{$start}{'abs'} $$ph{$end}{'abs'}\n";
	my $pcheck;
	for(my $rpos=$start;$rpos<$end;$rpos++){
		$pcheck=$rpos;
		$i++; ## tells us the relative pos in ORF, important to get correct frame later
		
		## now we recalc the relative coords to genomic coords
		if($options{'K'}){ ## in circ mode
			while($pcheck -$ph{'circl'} >= 0){
				$pcheck-=$ph{'circl'};
			}
			$pcheck+=$ph{'cdssr'}; ## add the offset, now we are at correct genomic position again
		}
		
		

		$pos=$$ph{$pcheck}{'abs'};
		$pos++; # our gap index is one based but we are 0 based in utrs, fault of the u if (necessary?)
		$pos-- if($$ph{'strand'} eq '-'); ## if we are on minus strand then gap start is on pos downstream
		
		if($$href{$$species[0]}{$pos}{'gaps'}){ ## we found gaps in alignmet, (insertions in other species), now check which and if we need to kick them out ...
			foreach $s(@$species){
				next if($s eq $$species[0]); ## next if refspecies
				$keep=1;
				if($$href{$s}{$pos}{'ins'}){ ## ok we have insertions in the species

					if(($$href{$s}{$pos}{'ins'} %3 ==0) and ($$href{$s}{$pos}{'ins'} < $end-$start)){ ## length of insert < orflength and multiple of 3
						$keep =1;
					}else{
						$keep=0;
					}
				}
				if($keep ==0){
					if(exists $$hs{$s}){
						delete $$hs{$s};
					}
				}
			}
		}
	}
}


sub check_codon_pos{
	my ($cp,$start,$end,$strandtx,$strand,$kick,$ph)=@_;
	my ($pos,$pos2);
	my $pcheck;

	my $gstrand='+';
	if($strand eq '+'){
		$gstrand=$strandtx;
	}else{
		if($strandtx eq '+'){
			$gstrand='-';
		}else{
			$gstrand='+';
		}
	}
#	print "$$ph{400}{'abs'}   $$cp{142264327}\n" if($start == 400);

	for(my $rpos=$start;$rpos<$end;$rpos+=3){
		$pcheck=$rpos;
		## now we recalc the relative coords to genomic coords
		if($options{'K'}){ ## in circ mode
			while($pcheck -$ph{'circl'} >= 0){
				$pcheck-=$ph{'circl'};
			}
			$pcheck+=$ph{'cdssr'}; ## add the offset, now we are at correct genomic position again
		}
		$pos=$$ph{$pcheck}{'abs'};
#		print "$rpos\t$pos\n";
	   


		if($strand eq '+'){ ## so 
			## this is for same direction only on same strand 
			if(exists $$cp{$pos}){
				if($$cp{$pos} =~ /[$strandtx]/){
					$$kick{$pos}=1;
					
				}
			}			
		}else{
			if($strandtx eq '+'){ ### and strand eq '-'
				## this is for antisense to plus transcript
                $pos2=$pos;
				$pos2++;       ## this is correct since we go from left to right with increasing genomic positions
#				$pos2--;
				if(exists $$cp{$pos2}){
					if($$cp{$pos2} =~ /[$strandtx]/){
						$$kick{$pos}=1;
					}
				}
			}else{#strandtx eq '-' and strand eq '+' this is for antisense to minus trancscript
				$pos2=$pos;
				$pos2--;     ## this is correct since we go from right to left with decreasing genomic positions
				if(exists $$cp{$pos2}){
					if($$cp{$pos2} =~ /[$strandtx]/){
						$$kick{$pos}=1;
					}
				}
			}
		}
	}	
}



sub readin_annotation_files{
	my ($a,$h)= @_;
	my @l;

	my @files=split(",",$a);
	foreach my $f(@files){
		open IN,$f or die "Could not open annotation files\n";
		while(<IN>){
			@l=split();
			$$h{$l[0]}{'type'} = $l[1];
			$$h{$l[0]}{'transcript_id'} = $l[4];
		}
		close IN
	}
}

## sub-routines
sub readin_refseq{
	my ($fhh,$ltp,$f)=@_;
	my $pseudo='';
	if($f =~ /pseudo/i){
		$pseudo = $pseudo;
	}

	my $extra='';
	my ($utr,$cds,%ee,%es,$ie);
	my @line;
	my $id;

	my ($exon,$pos)=();
	my ($cname,$cchr,$cstrand,$ctxs,$ctxe,$ccds,$ccde,$ces,$cee)=(-1,-1,-1,-1,-1,-1,-1,-1,-1);
	if($ltp =~ /^#/){
		$ltp =~ s/chromStart/txStart/;
		$ltp =~ s/chromEnd/txEnd/;
		$ltp =~ s/chromStarts/blockStarts/;
		
		chomp $ltp;
		@line = split(/\s+/,$ltp);

		for(my $c=0;$c < scalar @line;$c++){
			if($line[$c] =~ /name$/i){
				if($cname == -1){
					$cname=$c;
				}else{
					die "error in cname assignment from f $f, check header line please\n";
				}
			}elsif($line[$c] =~ /chr/i){
				if($cchr== -1){
					$cchr=$c;
				}else{
					die "error in chr assignment from file $f, check header line please\n";
				}
			}elsif($line[$c] =~ /strand/i){
				if($cstrand == -1){
					$cstrand=$c;
				}else{
					die "error in cstrand assignment from file $f, check header line please\n";
				}
			}elsif($line[$c] =~ /txStart$/i){
				if($ctxs == -1){
					$ctxs=$c;
				}else{
					die "error in ctxs assignment from file $f, check header line please\n";
				}

			}elsif($line[$c] =~ /txEnd$/i){
				if($ctxe == -1){
					$ctxe=$c;
				}else{
					die "error in ctxe assignment from file $f, check header line please\n";
				}
			}elsif($line[$c] =~ /cdsStart$/i){
				if($ccds == -1){
					$ccds=$c;
				}else{
					die "error in ctxe assignment from file $f, check header line please\n";
				}
			}elsif($line[$c] =~ /cdsEnd$/i){
				if($ccde == -1){
					$ccde=$c;
				}else{
					die "error in ctxe assignment from file $f, check header line please\n";
				}
			}elsif($line[$c] =~ /exonStarts/){
				$ces = $c;
			}elsif($line[$c] =~ /exonEnds/){
				$cee = $c;
			}
		}
	}

	$$fhh{$f}{'token'}{'txstart'} = $ctxs;
	$$fhh{$f}{'token'}{'txend'} = $ctxe;
	$$fhh{$f}{'token'}{'cdsstart'} = $ccds;
	$$fhh{$f}{'token'}{'cdsend'} = $ccde;
	$$fhh{$f}{'token'}{'exonstarts'} = $ces;
	$$fhh{$f}{'token'}{'exonends'} = $cee;
	$$fhh{$f}{'token'}{'strand'} = $cstrand;
	$$fhh{$f}{'token'}{'chr'} = $cchr;
	$$fhh{$f}{'token'}{'name'} = $cname;
	$$fhh{$f}{'token'}{'pseudo'} = $pseudo;
	return;
}


## sub-routines
sub readin_refseq_old{
  my ($genome,$file,$g2)=@_;
  my $pseudo='';
  my $extra='';
  if($file =~ /Pseudo/i){
	  $pseudo='pseudo-';
  }

  my ($utr,$cds,%ee,%es,$ie);
  my @line;
  my $id;
  print STDERR "reading file $file\n";
  open IN,"<$file" or die "No file with refseq gene Annotations given\n";
  my ($exon,$pos)=();
  my ($cname,$cchr,$cstrand,$ctxs,$ctxe,$ccds,$ccde,$ces,$cee)=(-1,-1,-1,-1,-1,-1,-1,-1,-1);
  while(<IN>){
	  if(/^#/){
		  $_ =~ s/chromStart/txStart/;
		  $_ =~ s/chromEnd/txEnd/;
		  chomp;
		  @line = split(/\s+/);
			for(my $c=0;$c < scalar @line;$c++){
				if($line[$c] =~ /name$/i){
					if($cname == -1){
						$cname=$c;
					}else{
						die "error in cname assignment from file $file, check header line please\n";
					}
				}elsif($line[$c] =~ /chr/i){
					if($cchr== -1){
						$cchr=$c;
					}else{
						die "error in cname assignment from file $file, check header line please\n";
					}
				}elsif($line[$c] =~ /strand/i){
					if($cstrand == -1){
						$cstrand=$c;
					}else{
						die "error in cstrand assignment from file $file, check header line please\n";
					}
				}elsif($line[$c] =~ /txStart$/i){
					if($ctxs == -1){
						$ctxs=$c;
					}else{
						die "error in ctxs assignment from file $file, check header line please\n";
					}
					
				}elsif($line[$c] =~ /txEnd$/i){
					if($ctxe == -1){
						$ctxe=$c;
					}else{
						die "error in ctxe assignment from file $file, check header line please\n";
					}
				}elsif($line[$c] =~ /cdsStart$/i){
					if($ccds == -1){
						$ccds=$c;
					}else{
						die "error in ctxe assignment from file $file, check header line please\n";
					}
				}elsif($line[$c] =~ /cdsEnd$/i){
					if($ccde == -1){
						$ccde=$c;
					}else{
						die "error in ctxe assignment from file $file, check header line please\n";
					}
				}elsif($line[$c] =~ /exonStarts/){
					$ces = $c;
				}elsif($line[$c] =~ /exonEnds/){
					$cee = $c;
				}

			}
			next;
	  }
	  
	  chomp;
	  ($utr,$cds,%ee,%es,$ie)=();
	  @line = split(/\s+/);
	  ## skip if not correct chromosome
	  $id=$line[$cname];
	  next if($line[$cchr] ne $chr);
	  
	  ## newly added 30.09.13
	  if(exists $anno{$id}){
		  $extra="$anno{$id}{'type'}_";
		  print $extra,"\n";
	  }else{
		  $extra='';
	  }

	  $tss{$line[$cchr]}{$line[$cstrand]}{$line[$ctxs]}.="$extra$pseudo$id,";
	  $tse{$line[$cchr]}{$line[$cstrand]}{$line[$ctxe]}.="$extra$pseudo$id,";
	  ## for each position between txStart and txEnd
	  if($ces != -1){
	     @es{split(",",$line[$ces])}=undef; ## put array as keys in hash
	     @ee{split(",",$line[$cee])}=undef;
	  }
	  
	  ($utr,$cds,$exon)=("","","");
	  for(my $i=$line[$ctxs];$i<=$line[$ctxe];$i++){
		  if(exists $es{$i}){
			  $exon=",exon;";
		  }elsif(exists $ee{$i+1}){ ## bed files have end coord not included
			  $exon=",intron;";
		  }else{
			  $exon='';
		  } 

		  ## add utr info if necessary     ## would also need to check if in utr exon or not ... -> better check spliced 5utrs alone -> did this already ...
		  if($ccds != -1){
			  if($line[$cstrand] eq '+'){
				  if($i < $line[$ccds]){
					  $utr=',utr5';
				  }elsif($i > $line[$ccde]){
					  $utr=',utr3';
				  }else{
				  $utr='';
				  }
			  }else{
				  if($i < $line[$ccds]){
					  $utr=',utr3';
				  }elsif($i > $line[$ccde]){
					  $utr=',utr5';
				  }else{
					  $utr='';
				  }
			  }
			  
			  if($i >= $line[$ccds] and $i < $line[$ccde]){
				  $cds=',cds';
			  }else{
				  $cds = '';
			  }
		  
		  

			  ## if noncoding then no utr present and not cds present but exons may be present
			  if($line[$ccds] == $line[$ccde]){
				  $cds='';
				  $utr='';
			  }
		  }

		  if($line[12]){
			  if($line[$cname] =~ /NR_/ and $g2 != 0){
				  $$g2{$line[$cchr]}{$line[$cstrand]}{$i} .= $line[12].",$pseudo".$line[$cname]."$utr$cds$exon";
			  }else{ ## this is for protein coding regions
				  $$genome{$line[$cchr]}{$line[$cstrand]}{$i} .= $line[12].",$pseudo".$line[$cname]."$utr$cds$exon";
			  }
		  }else{
			  if($line[$cname] =~ /NR_/ and $g2 != 0){
				  $$g2{$line[$cchr]}{$line[$cstrand]}{$i} .= "$extra$pseudo$line[$cname]$utr$cds$exon";
			  }else{	## this is for protein coding regions			
					  $$genome{$line[$cchr]}{$line[$cstrand]}{$i} .= "$extra$pseudo$line[$cname]$utr$cds$exon";				
			  }
		  }
	  }
  }
  close IN;
  print STDERR "reading file $file finished\n";
}



## check for overlap in genes
				
sub overlap_wrapper{
	my ($chr,$strand,$pos,$e,$hashc,$hashnc,$hasho,$igs,$cons,$orf_len,$syn_mut,$aa,$dna,$dnaH,$oe,$pe,$se,$synmutstat,$stringATG,$dndsout,$fsios,$tmpseq) = @_;

	my @l;
	my $coding = check_overlap($chr,$strand,$pos,$e,$hashc,$igs);
	my $ncoding= check_overlap($chr,$strand,$pos,$e,$hashnc,$igs);

	my ($overlap_with,$res) =('','');

	if($options{'o'}){
	 	foreach my $entry(keys %$hasho){
	 		$res=check_overlap($chr,$strand,$pos,$e,$$hasho{$entry},$igs);
	 		$overlap_with.="$entry," if($res);
	 	}
	}

	## this could be refined so that we now if something is antisense to coding or a noncoding RNA
	my $outstr="$chr\t$strand\t$pos\t$e\t$orf_len\t0\t".(scalar @$cons)."\t$syn_mut\t$aa\t$dna\t$dnaH\t$oe\t$pe\t$se\t".(join(",",@$cons))."\t-\t-\t-\t-\t-\t-\t$ncoding,$coding,$overlap_with".(getCodonMut($tmpseq,$cons,0,$orf_len,$strand)).$synmutstat.$stringATG.$dndsout.$fsios."\n";

	if($coding){
#		print STDERR "$dna,$pos,$e\n" if($pos ==5080121);
		## get new AS category
		if($coding =~ /__AS__/){
			print ASC $outstr;
			return 1;
		}
			

		if(substr($coding,0,1) eq 1){
			print COD $outstr;
			return 1;
		}
	}
	
	if($ncoding){
		if($ncoding =~ /__AS__/){
			print ASN $outstr;
			return 1;
		}
#		print STDERR "noncoding $dna,$pos,$e",substr($ncoding,0,1),"\n" if($pos ==5080121);
		if(substr($ncoding,0,1) eq 1){
			
			print NC $outstr;
			return 1;
		}
	}
	$overlap_with.=";$coding;$ncoding;";
	return $overlap_with;
}



sub check_overlap{
	
	my ($chr,$strand,$start,$end,$hash,$ignore)=@_;
	my $ret='';

	if($$hash{$chr}{$strand}{$start} ){
		$ret.= "1:b:$$hash{$chr}{$strand}{$start};";
	}
	if($$hash{$chr}{$strand}{$end} ){
		$ret.= "1:e:$$hash{$chr}{$strand}{$end};";
	}
	
#	if($ignore){
	$strand=$ignores{$strand};
		
	if($$hash{$chr}{$strand}{$start} ){
		$ret.= "1:b:__AS__$$hash{$chr}{$strand}{$start};";
	}
	if($$hash{$chr}{$strand}{$end} ){
		$ret.= "1:e:__AS__$$hash{$chr}{$strand}{$end};";
	}
	return $ret;
}




sub check_overlap_ra{
	## fuction needs to get a filehandle
	## must be called like check(*IN,1234)
	## the star is really important
	## this is random access to index file, saves a lot of time

	my ($fh,$position,$end)=@_;
	my $l;
#	open IN,$file or die "No file given\n";
	#my $pos=3*$position;
	my $pos;
	for(my $i=$position;$i<$end;$i++){
		$pos=3*$i;
		seek $fh,$pos,0;
		read($fh,$l,2);
		chomp $l;
		if($l =~ /^_+$/){
			
		}else{
			return($l);
		}
	}
	return 0;
}


sub detect_orfs{
    my ($maf,$positions,$species,$refseq,$nc_refseq,$oa_files,$igs,$window,$kaks,$finalpos,$ph,$itergap,$cp,$iif,$skipnt) = @_; 
    my $start=0;
    my ($message,$cchrom,$pos);
	my $overlap_with;
    my @cons;

    my $tri="";
    my ($b_f,$b_r,$e,$e_r,$frame,$frame_r,$s1,$s1_r,$s2,$s2_r,$s3,$s3_r) = 0;
    my ($f1,$f2,$f3,$f1_r,$f2_r,$f3_r) = (-1,-1,-1,-1,-1,-1);
    
    my $orf_len;
	my $shined;
	my $dnaH;

    ## make reference to ref_sequence
	
    my $ref_seq=\$$maf{$$species[0]};
    my $lref = length($$ref_seq);

    my ($pe,$oe,$se);
   
    my $syn_mut;
    my $res;
	my $kozak;
	my $tssa;
	my $tesa;
	my ($coding,$nccoding,$strand);

	my $line;
	my @annol;

	my $totalC=-1;
	my $pcheck;
	my %gotcha;

	

	if($options{'P'}){
		## here we set the current file pos and get the iter!!!
		#		  open my $iif,"$options{'P'}/$options{'c'}.maf.index.sorted" or die "index file could not be opened\n";
		my $pos_to_go= $$ph{0}{'abs'} < $$ph{$lref-1}{'abs'} ? $$ph{0}{'abs'} : $$ph{$lref-1}{'abs'}; 
		
		my $res=binarysearch2(\@iib,$pos_to_go);
		
		my $hi=$iib[$res];
		my $pos=$ii{$hi};
		while($hi > $pos_to_go){
			$res--;
			$hi=$iib[$res];
			$pos=$ii{$hi};
		}
		seek $iif,$pos,0;
	}
	
	


    for(my $i=$skipnt;$i<$lref;$i++){
		$pcheck=$i;
		## now we recalc the relative coords to genomic coords
		if($options{'K'}){ ## in circ mode
			$pcheck=$i;
			while($pcheck -$ph{'circl'} >= 0){
				$pcheck-=$ph{'circl'};
			}
			$pcheck+=$ph{'cdssr'}; ## add the offset, now we are at correct genomic position again
			$gotcha{$pcheck}=1;
		}
		
		$pos=$$ph{$pcheck}{'abs'};

	  if($i%1000000 == 0 and $i > 0){
	    print STDERR Nicenumber($i)," bases from ",Nicenumber($lref)," processed\n";  
	  }

	  
	  ## here we do the actual annotation file readin with an iterator
	  foreach my $f(keys %fhh){
        next if(not defined $fhh{$f}{'l'});
        @annol=split(/\s+/,$fhh{$f}{'l'});
		#	print "$f  $fhh{$f}{'type'} \n";
		
		if($fhh{$f}{'type'} eq 'c'){
		  #			print "$f\t$i\n" if($i == 11336);
		  check($i,$refseq,\@annol,$f,$max_size,\%fhh,$nc_refseq,$chr,$i);
		}else{
		  check($i,$nc_refseq,\@annol,$f,$max_size,\%fhh,0,$chr,$i);
		}
    }

	  $tri = uc(substr($$ref_seq,$i,3)); 
	  
	  next if($tri =~ /-/); ## skip if gap is encountered in ref_sequence ## this is necessary,added 25.01.2013
	   

	  if($options{'P'}){
## here we set the current file pos and get the iter!!!
#		  open my $iif,"$options{'P'}/$options{'c'}.maf.index.sorted" or die "index file could not be opened\n";
#		  my $res=binarysearch2(\@iib,$pos_to_go);
#		  my $hi=$iib[$res];
#		  my $pos=$ii{$hi};
#		  seek $iif,$pos,0;
		if(not $options{'b'}){
		  if(defined $l[1] and $l[1] =~ /\d/){
			
			#			if(not $options{'b'}){
			if($l[1] > $i+$max_size+2){ ## if we are above the max size delete this
			}elsif($l[1] < $i -$max_size+2 and $l[1] > $i - (3*$max_size)){ ## as long as we are smaller + window than current pos take next line
			  while($l[1] < $i - $max_size+2){
				$line=$itergap->();
				if(defined $line){
				  @l=split(/\s+/,$line);
				}else{
				  @l=();
				}
				last if(not $l[1]);
			  }
			}else{
			  while($l[1] > $i -($max_size+2) and $l[1] < ($i+$max_size+2)){
				$gapindex{$l[0]}{$l[1]}{'gaps'}=$l[2];
				$gapindex{$l[0]}{$l[1]}{'ins'}=$l[3];
				$line=$itergap->();
				if(defined $line){
				  @l=split(/\s+/,$line);
				}else{
				  @l=();
				}
				last if(not $l[1]);
				#				print STDERR "insert $l[0] $l[1]\n" if($i == 842679);
				#				print STDERR "insert $line at $i\n" if($l[1] == 842679);
			  }
			}
		  }

		  if($i >= 2*$max_size){
			foreach my $s(@species){
			  delete $gapindex{$s}{$i-$max_size-2} if(exists $gapindex{$s}{$i-$max_size-2});
			}
		  }
		}else{
		  ## we have options b we need to check if the abs position is gapped or not, so how do we do this
		  ## what to read in from gap index now ???
		  ## tx start to tx end and for next transcript delete everything until next txstart so we should be fine, no ?
			
		  if(defined $l[1] and $l[1]=~ /\d/){
			## we crossed tx end already 
			  if($l[1] >= $$ph{'txe'}){
			  }elsif($l[1] < $$ph{'txs'}){

				  ## as long as gap pos is smaller then our tx start get next lines ....
				  while($l[1] < $$ph{'txs'}){
					  #print "$l[1] $$ph{'txs'} skipping\n";
					  $line=$itergap->();
					  if(defined $line){
						  @l=split(/\s+/,$line);
					  }else{
						  @l=();
					  }
					  last if(not $l[1]);
#					  print $l[1],"\n";
				  }
			  }
			  
			  
			  
			  

			  ## now we are at txs or above
			  if($l[1] >= $$ph{'txs'} and $l[1] < $$ph{'txe'}){ ## we are in the transcript right now so we read until we reach what we need to get to txe
				  while($l[1] >= $$ph{'txs'} and $l[1] < $$ph{'txe'}){ ## insert into gap hash everything in transcript
					  $gapindexr{$l[1]}{$l[0]}++;
					  $gapindex{$l[0]}{$l[1]}{'gaps'}=$l[2];
					  $gapindex{$l[0]}{$l[1]}{'ins'}=$l[3];
					  $line=$itergap->();
					  if($line =~ /spec/){$line=$itergap->();}
					  if(defined $line){
						  @l=split(/\s+/,$line);
					  }else{
						  @l=();
					  }
					  last if(not $l[1]);
				  }
			  }
		  }else{
			  print STDERR "some error in checking gaps here ...\n";
		  }

		  ## when and how to delete stuff
	  }
	}
		if( $tri eq "ATG"){
			$b_f = $i+$start;
			$frame = $i%3;
			if($frame == 0){ $f1 = $b_f; $s1=$i; } 
			if($frame == 1){ $f2 = $b_f; $s2=$i; }
			if($frame == 2){ $f3 = $b_f; $s3=$i; }
		}

	  
	  
	  if($tri eq "TAA" or $tri eq "TGA" or $tri eq "TAG"){
		$strand = '+';
		$b_f = $i+$start;
		$e = $b_f+2;      ## define correct end, coord included
		$frame = $i%3;
#		print "$tri\t$i\t $f2 $e $frame $f1 $f2 $f3\n"; ## orf goes from 601-739(741)
#		die;



		my $tr2;
		$tssa=\@tssap;
		$tesa=\@tesap;
		
		if($frame == 0 and $f1 != -1){
		  check_orf($strand,$frame,\$f1,$e,$start,$maf,$species,$refseq,$nc_refseq,$oa_files,$igs,$window,$ref_seq,$tssa,$tesa,$kaks,$finalpos,$ph,$totalC,$cp);
	
        }
		if($frame == 1 and $f2 != -1){
		  check_orf($strand,$frame,\$f2,$e,$start,$maf,$species,$refseq,$nc_refseq,$oa_files,$igs,$window,$ref_seq,$tssa,$tesa,$kaks,$finalpos,$ph,$totalC,$cp);
	
        }
		if($frame == 2 and $f3 != -1){
		  check_orf($strand,$frame,\$f3,$e,$start,$maf,$species,$refseq,$nc_refseq,$oa_files,$igs,$window,$ref_seq,$tssa,$tesa,$kaks,$finalpos,$ph,$totalC,$cp);
	
    	}
		
	  }
	  
	 
	## now for the reverse strand
		if($options{'q'}){
		if($tri eq "TTA" or $tri eq "CTA" or $tri eq "TCA"){
		$b_r = $i+$start;
		#print "start codon found in frame $frame at pos $i >> $b\n";
		$frame_r = $i%3;
		if($frame_r == 0){ $f1_r = $b_r;$s1_r=$i;}
		if($frame_r == 1){ $f2_r = $b_r;$s2_r=$i;}
		if($frame_r == 2){ $f3_r = $b_r;$s3_r=$i;}
	  }
	 
	  
	  if( $tri eq "CAT"){                 
        $strand = '-';
		$b_r = $i+$start;
		$e_r = $b_r+2;
		$frame_r = $i%3;
		
		
		$tssa=\@tssam;
		$tesa=\@tesam;
    	my $tr2;
		
	    ## now we run it for the minus strand
    	if($frame_r == 0 and $f1_r != -1){
		  check_orf($strand,$frame_r,\$f1_r,$e_r,$start,$maf,$species,$refseq,$nc_refseq,$oa_files,$igs,$window,$ref_seq,$tssa,$tesa,$kaks,$finalpos,$ph,$totalC,$cp);
        }
		if($frame_r == 1 and $f2_r != -1){
		  check_orf($strand,$frame_r,\$f2_r,$e_r,$start,$maf,$species,$refseq,$nc_refseq,$oa_files,$igs,$window,$ref_seq,$tssa,$tesa,$kaks,$finalpos,$ph,$totalC,$cp);
        }
		if($frame_r == 2 and $f3_r != -1){
		  check_orf($strand,$frame_r,\$f3_r,$e_r,$start,$maf,$species,$refseq,$nc_refseq,$oa_files,$igs,$window,$ref_seq,$tssa,$tesa,$kaks,$finalpos,$ph,$totalC,$cp);
    	}
      } 
    }
	}
}

sub detect_orfs_exons{
    my ($maf,$positions,$species,$refseq,$nc_refseq,$oa_files,$igs,$window,$kaks,$finalpos,$ph,$itergap,$cp,$iif) = @_; 
    my $start=0;
    my ($message,$cchrom,$pos);
	my $overlap_with;
    my @cons;

    my $tri="";
    my ($b_f,$b_r,$e,$e_r,$frame,$frame_r,$s1,$s1_r,$s2,$s2_r,$s3,$s3_r) = 0;
    my ($f1,$f2,$f3,$f1_r,$f2_r,$f3_r) = (-1,-1,-1,-1,-1,-1);
    
    my $orf_len;
	my $shined;
	my $dnaH;

    ## make reference to ref_sequence
    my $ref_seq=\$$maf{$$species[0]};
    my $lref = length($$ref_seq);

    my ($pe,$oe,$se);
   
    my $syn_mut;
    my $res;
	my $kozak;
	my $tssa;
	my $tesa;
	my ($coding,$nccoding,$strand);

	my $line;
	my @annol;

	my $totalC=-1;
	my $pcheck;
	my %gotcha;
	my $i=0; ## we take whole exon



	$pcheck=$i;
	## now we recalc the relative coords to genomic coords
		
	$pos=$$ph{$pcheck}{'abs'};


	  if($i%1000000 == 0 and $i > 0){
	    print STDERR Nicenumber($i)," bases from ",Nicenumber($lref)," processed\n";  
	  }

	  
	  ## here we do the actual annotation file readin with an iterator
	  foreach my $f(keys %fhh){
        next if(not defined $fhh{$f}{'l'});
        @annol=split(/\s+/,$fhh{$f}{'l'});
		#	print "$f  $fhh{$f}{'type'} \n";
		
		if($fhh{$f}{'type'} eq 'c'){
		  #			print "$f\t$i\n" if($i == 11336);
		  check($i,$refseq,\@annol,$f,$max_size,\%fhh,$nc_refseq,$chr,$i);
		}else{
		  check($i,$nc_refseq,\@annol,$f,$max_size,\%fhh,0,$chr,$i);
		}
    }

	  $tri = uc(substr($$ref_seq,$i,3)); 
	  
	  next if($tri =~ /-/); ## skip if gap is encountered in ref_sequence ## this is necessary,added 25.01.2013
	  

	  if($options{'P'}){
		if(not $options{'b'}){
		  if(defined $l[1] and $l[1] =~ /\d/){
			

			if($l[1] > $i+$max_size+2){ ## if we are above the max size delete this
			}elsif($l[1] < $i -$max_size+2 and $l[1] > $i - (3*$max_size)){ ## as long as we are smaller + window than current pos take next line
			  while($l[1] < $i - $max_size+2){
				$line=$itergap->();
				if(defined $line){
				  @l=split(/\s+/,$line);
				}else{
				  @l=();
				}
				last if(not $l[1]);
			  }
			}else{
			  while($l[1] > $i -($max_size+2) and $l[1] < ($i+$max_size+2)){
				$gapindex{$l[0]}{$l[1]}{'gaps'}=$l[2];
				$gapindex{$l[0]}{$l[1]}{'ins'}=$l[3];
				$line=$itergap->();
				if(defined $line){
				  @l=split(/\s+/,$line);
				}else{
				  @l=();
				}
				last if(not $l[1]);
			}
		  }
		}

		  if($i >= 2*$max_size){
			  foreach my $s(@species){
				  delete $gapindex{$s}{$i-$max_size-2} if(exists $gapindex{$s}{$i-$max_size-2});
			  }
		  }
		}else{
		  ## we have options b we need to check if the abs position is gapped or not, so how do we do this
		  ## what to read in from gap index now ???
		  ## tx start to tx end and for next transcript delete everything until next txstart so we should be fine, no ?
			
		  if(defined $l[1] and $l[1]=~ /\d/){
			## we crossed tx end already 
			  if($l[1] >= $$ph{'txe'}){
			  }elsif($l[1] < $$ph{'txs'}){

				  ## as long as gap pos is smaller then our tx start get next lines ....
				  while($l[1] < $$ph{'txs'}){
					  #print "$l[1] $$ph{'txs'} skipping\n";
					  $line=$itergap->();
					  if(defined $line){
						  @l=split(/\s+/,$line);
					  }else{
						  @l=();
					  }
					  last if($l[1] !~ /\^\d+$/);
				  }
			  }
			  
			  
			  

			  ## now we are at txs or above
			  if($l[1] =~ /\^\d+$/){
			  if($l[1] >= $$ph{'txs'} and $l[1] < $$ph{'txe'}){ ## we are in the transcript right now so we read until we reach what we need to get to txe
				  while($l[1] >= $$ph{'txs'} and $l[1] < $$ph{'txe'}){ ## insert into gap hash everything in transcript
					  $gapindexr{$l[1]}{$l[0]}++;
					  $gapindex{$l[0]}{$l[1]}{'gaps'}=$l[2];
					  $gapindex{$l[0]}{$l[1]}{'ins'}=$l[3];
					  $line=$itergap->();
					  if(defined $line){
						  @l=split(/\s+/,$line);
					  }else{
						  @l=();
					  }
					  last if(not $l[1]);
				  }
			  }
			}
		  }
		  ## when and how to delete stuff
	  }
	}

		$b_f = $i+$start;
		$frame = $i%3;
		if($frame == 0){ $f1 = $b_f; $s1=$i; } 
		if($frame == 1){ $f2 = $b_f; $s2=$i; }
		if($frame == 2){ $f3 = $b_f; $s3=$i; }


	  
		$strand = '+';
		$b_f = $i+$start;
		$e = $b_f+2;      ## define correct end, coord included
	    $e = $lref - 1 ;
		$frame = $i%3;
	



		my $tr2;
		$tssa=\@tssap;
		$tesa=\@tesap;
		
		if($frame == 0 and $f1 != -1){
		  check_orf($strand,$frame,\$f1,$e,$start,$maf,$species,$refseq,$nc_refseq,$oa_files,$igs,$window,$ref_seq,$tssa,$tesa,$kaks,$finalpos,$ph,$totalC,$cp);
        }
	  return;
}




sub check_orf{
	my ($strand,$frame,$fs,$e,$start,$maf,$species,$refseq,$nc_refseq,$oa_files,$igs,$window,$ref_seq,$tssa,$tesa,$kaks,$finalpos,$ph,$totalC,$cp) = @_;
	my ($oe,$pe,$se,$tr2,$dnaH,$overlap_with,$kozak,$shined)=("-","-","-","-","-",'','','','','');
	my $codons_scored=0;
	my @cons;
	my $orf_len = ($e-$$fs+1);
	my $pos = $$fs-$start;
 
	my %hs;
	my @nspec=();
	my ($syn_mut,$nsyn_mut,$uniq_syn_mut,$uniq_nsyn_mut,$ka,$ks,$omega);
	my %kaks_pairwise;

	## R is set if we want to skip repeat regions completely so we skip orf if in repeat region
	if($options{'R'}){
		if(substr($$ref_seq,$$fs,$orf_len) =~ /[acgtn]+/){
			print ERR substr($$ref_seq,$$fs,$orf_len),"\tcontains repeat elements\n";
			return;
		}
	}

	my $fst=$pos;



	$$fs = -1;
	## if this is the case try to find a met before without passing the maxsize limitation
	my $new=$pos;return if($orf_len >= $max_size);
	if($strand eq '+'){
		if($orf_len < $max_size){
			for(my $i=$pos;$i > $e-$max_size;$i-=3){
				last if($i > $finalpos); ## if i bigger than that stop
				$tr2=substr($$ref_seq,$i,3);
				## R is set if we want to skip repeat regions completely so we dont extend orfs into repeat regions 
				if($options{'R'}){
					last if($tr2 =~ /[actgn]+/);
				}
				$tr2=uc($tr2);
				
				
				last if($i < 0);
				last if($tr2 eq "TAA" or $tr2 eq "TGA" or $tr2 eq "TAG");
				$new=$i if($tr2 eq 'ATG');
			}
		}
		$pos=$new;
	}else{
		
		$new = $e;
		if($orf_len < $max_size){
			
			for(my $i=$e+1;$i < $pos+$max_size-2;$i+=3){
				last if($i > $finalpos); ## if i bigger than that stop
				$tr2=substr($$ref_seq,$i,3);
				if($options{'R'}){
					last if($tr2 =~ /[actgn]+/);
				}
				$tr2=uc($tr2);
				
				
				last if($i < 0);
				last if($tr2 eq "TTA" or $tr2 eq "CTA" or $tr2 eq "TCA");
				if( $tr2 eq 'CAT'){
					$new=$i;
				}
			}
			    }
		$new+=2 if($new != $e);
		if($e != $new){
			$e=$new;
		}
	}
	

	$orf_len = ($e-$pos+1);


	
	return if($e > $finalpos); ## if we left the chromosome already return
	return if($orf_len < $min_orf_len);



	if($options{'K'}){
		return if($pos > 2*$ph{'circl'} and $totalC == -1); ## we abort searching here if not circle was found in the first to round of searching
	}
	#if($fst==601 ){die "$pos  $e $strand  $finalpos here";}		



	my $dna=uc(substr($$ref_seq,$pos,$orf_len));
	
	return if($dna =~ /N/); ## why is he not returning here ???	
	if($options{'K'}){
		return if(exists $seen{$dna});
		$seen{$dna}=();
	}
	
	if($strand eq '-'){ $dna = rc($dna);}
	
	
	my $aa=peptide($dna);

	if($options{'A'}){
		if($options{'b'}){
			print IG "$$ph{'name'}\t$$ph{'strand'}";
		}else{
			print IG "$chr\t$strand";
		}
		print IG "\t+\t$pos\t$e\t$orf_len\t$frame\t-\t-\t$aa\t$dna\n";
		return;
	}
	$$ph{'strand'} = 'na' if(not defined $$ph{'strand'});
	$ph{$pos}{'abs'} = 'na' if(not defined $ph{$pos}{'abs'});
	$ph{$e}{'abs'} = 'na' if(not defined $ph{$e}{'abs'});
	$$ph{'name'} = 'na' if(not defined $$ph{'name'});
	
	my $erro="$$ph{'name'}\t$$ph{'strand'}\t$strand\t$pos\t$e\t$orf_len\t$aa\t$dna\t$options{'c'}\t$ph{$pos}{'abs'}\t$ph{$e}{'abs'}";
	#print "here $$ph{'name'}\t$$ph{'strand'}\t$strand\t$pos\t$e\t$orf_len\t$aa\t$dna\t$options{'c'}\t$ph{$pos}{'abs'}\t$ph{$e}{'abs'}\n";
	
	my %tmpseqm;
	my $tmpseq=\%tmpseqm;
	if($options{'b'}){
		foreach my $s(@$species){
			$$tmpseq{$s}=substr($$maf{$s},$pos,$orf_len);
		}
	}else{
		getseq($fh,$pos,$orf_len,$offset,$strand,$tmpseq,$species);
	}

	## the tmpseq hash holds only the orf sequence of the maf from the full transcript maf
	
	# print "$pos,$orf_len,$strand\n";
	# foreach my $s(@$species){print ">$s\n$tmpseq{$s}\n";}
	# exit;
									 
	
	@cons=();		    
	my @frameshiftI=();
	get_cons_species_relaxed(\@cons,$tmpseq,$species,0,$orf_len,$strand,\@frameshiftI);

	# foreach my $s(@cons){
	# 	print ">$s\n$$tmpseq{$s}\n";
	# }
	# exit;

	if(scalar @cons ==1){
		print ERR "$erro\tget_cons_species_relaxed\tno_cons_species_left\n";
		return;
	}


	my %atgstopcons=();
	#			print @cons,"\n";
	check_stop_and_start_codon_conservation($tmpseq,\@cons,\%atgstopcons,0,$orf_len,$strand);

	## here could further filter but for now we just report the numbers
	#			print "$transcript\tspecies\t",scalar @cons,"\nATG = $atgstopcons{$$species[0]}{'start'}\nSTOP cons = $atgstopcons{$$species[0]}{'stop'}{'cons'}\nSTOP synmut = $atgstopcons{$$species[0]}{'stop'}{'syn'}\n";			
		
	
	## here we filter out species with frameshift insertions 
	if ($options{'P'} and not $options{'b'}) {
		@nspec=();
		@hs{@cons} = undef;
		check_gaps($href,$pos,$pos+$orf_len-1,\%hs,\@cons,$strand);
		if(scalar keys %hs < 2){ ## if only the refspecies is left
			print ERR "$erro\tcheck_gaps\tonly ref_species left\n";
			return;
		}
		## now we fill up the nspec array
		foreach my $s (@cons) {
			if (exists $hs{$s}) {
				push(@nspec,$s); 
			} else {
				push(@frameshiftI,$s); ## also get here ones that have frameshift indels or huge insertions
			}
		}
		@cons=@nspec;			## now we set the new cons array
	}
	if(scalar @cons == 1){
		print ERR "$erro\tcheck_gaps;end of block\tonly ref_species left\n";
		return;
	}

	if ($options{'P'} and $options{'b'}) {
		@nspec=();
		@hs{@cons} = undef;
	  
		#		my          ($href,$start,$end,$hs,$species,$id,$strand,$ph)=@_;
#		check_gaps_tx($href,$pos,$pos+$orf_len-1,\%hs,\@cons,$strand,$ph);

		check_gaps_tx2($href,$pos,$pos+$orf_len-1,\%hs,\@cons,$strand,$ph);

		#		die $pos if(scalar keys %hs < 2); ## delme
        #   	die "$pos $e",keys %hs,"\n" if($e == 288);
		if(scalar keys %hs < 2){ ## if only the refspecies is left
			print ERR "$erro\tcheck_gaps_tx2\tno_cons_species_left\n";
			return;
		}


		## now we fill up the nspec array
		foreach my $s (@cons) {
			if (exists $hs{$s}) {
				push(@nspec,$s); 
			} else {
				push(@frameshiftI,$s); ## also get here ones that have frameshift indels or huge insertions
			}
		}
		@cons=@nspec;			## now we set the new cons array
	}
	
#	return if(scalar @cons == 1);
	if(scalar @cons ==1){
		print ERR "$erro\tcheck_gaps_tx2_end_routine\tno_cons_species_left\n";
		return;
	}
#	print IG "$erro\n";

#	return; ##leave routine


	## evtl. say that all species must have stop codon conserved
	my ($startC,$stopC,$stopS)=($atgstopcons{$$species[0]}{'start'}{'cons'},$atgstopcons{$$species[0]}{'stop'}{'cons'},$atgstopcons{$$species[0]}{'stop'}{'syn'});
	my $stringATG="\t$startC:";

	#,$stopC,$stopS";
	foreach my $s(@cons){
	  if($s eq $cons[0]){
		$stringATG.="$atgstopcons{$$species[0]}{'start'}{$s}";
	  }else{
		$stringATG.=",$atgstopcons{$$species[0]}{'start'}{$s}";
	  }
	}
	$stringATG.="\t$stopC:";
	foreach my $s(@cons){
	  if($s eq $cons[0]){
		$stringATG.="$atgstopcons{$$species[0]}{'stop'}{$s}";
	  }else{
		$stringATG.=",$atgstopcons{$$species[0]}{'stop'}{$s}";
	  }
	}
	
	if(not $options{'A'}){
	  return if(scalar @cons == 1);
	}
	
	
	## get number of synonymous mutations in orfs and ka ks
	%kaks_pairwise=();
	my $nonsense;
	($syn_mut,$nsyn_mut,$uniq_syn_mut,$uniq_nsyn_mut,$ka,$ks,$nonsense)=split(",",syn_mut($tmpseq,\@cons,0,$orf_len,$strand,$kaks,\%kaks_pairwise));
	
	if($ks == 0 or $ka ==0){
	  $omega="ka=$ka/ks=$ks";
	}else{
	  $omega=$ka/$ks;
	}
	my $synmutstat="\t$syn_mut,$nsyn_mut,$uniq_syn_mut,$uniq_nsyn_mut,$ka,$ks\t$omega";
	
	if($strand eq '+'){	
	  $kozak=uc(substr($$ref_seq,$pos-6,6));
	  $shined=0;
	  if(uc(substr($$ref_seq,$pos-20,20)) =~ /AGGAGG/){
		$shined=1;
	  }
	}else{
	  $kozak=rc(uc(substr($$ref_seq,$e+1,6)));
	  $shined=0;
	  if(rc(uc(substr($$ref_seq,$e+1,20))) =~ /AGGAGG/){
		$shined=1;
	  }   
	}
	
	## get dist do next TSS and next TSE;
	($dtss,$dtes)=('na','na'); 
	($dtss,$dtes)=split(",",get_dtss_dtes($tssa,$e,$tesa,$pos)) if(scalar @$tssa and scalar @$tesa);
	
	
	### here comes the whole dnds story
	my %mafp=();
	my %rflank=();
	my %lflank=();

	my %dnds=();
	my %kaks=();
	my %lambda=();
	
	my $flankl=$orf_len;
	
	my $id="${chr}_${pos}_${orf_len}_$strand";
	my $name=$id;
	
	#die $pos,$flankl,substr($$maf{'ce6'},$pos,$flankl),"\n";
	if($pos-$flankl >=0){ ## this is for genome mode only true
	}else{
	  $flankl=$pos;
#	  $flankl-=($flankl%3);## make it a multiple of 3 if necessary				
	  $flankl-=($flankl%3) if(not $options{'b'}); ## only in genome mode important and if we go really at chromosome borders, otherwise we are fine here
	}

	my $ts;
	my $diff = $orf_len-$flankl;
	
	## this also depends on plus or minus strand in general
	my $es='+';
	
	## the issue with the flank is that if a 5 or 3 utr is missing we dont get enoungh statistics
	if(not $options{'b'}){
	  getseq($fh,$pos-$flankl,$flankl,$offset,$strand,\%lflank,\@cons);
	  ## on chrM in flies we start at genomic position 1 so no left flank available
	  

	}else{
		if ($$ph{'txs'} < $orf_len or $$ph{'txe'}+$orf_len > $length_chr) {
		} else {
			## if diff is > 0 then the transcript is too short on the flank so diff is the number of nt we need to get from intergenic region
			if ($diff > 0) {
				if ($strand eq '+' and $$ph{'strand'} eq '-') {
					$es='-';
				} elsif ($strand eq '-' and $$ph{'strand'} eq '+') {
					$es='-';
				}
		
				## we get what is genomically upstream of the transcript
				if ($$ph{'strand'} eq '+') {
					getseq($fh,$$ph{'txs'}-$diff,$diff,$offset,$es,\%lflank,\@cons);
		  
					foreach my $s (@cons) {
						$ts=substr($$maf{$s},$pos-$flankl,$flankl); ## we start from 0 basically 
						if ($strand eq '-') { ## we are in antisense orientation from given maf extraction
							$ts=rc($ts);
							$lflank{$s} = $ts.$lflank{$s}; ## this is not correct yet since we are on the other side actually so it should be right flank than
							#die $lflank{$s}; ## works this is then the right flank actually, but we fix this in the output so no problems, so next codon right to it must be a stop
			  
						} else {
							$lflank{$s} = $lflank{$s}.$ts;
							#die $lflank{$s}; # works, so next codon right to it must be ATG
							#die "$pos  $orf_len  $flankl $diff  here $lflank{$s}\n";
						}
					}
				} else {		## transcript is on genomic minus
					## we get what is genomically downstream so we really get the left flank for things in sense orientation
					## we need to check here that we dont cross the end of the genome ###
		  
					if ($$ph{'txe'} + $diff > $length_chr) {
						$diff = $length_chr-$$ph{'txe'};
						$diff-=($diff %3); ## make seq a multiple of 3 , really important otherwise it will crash
					}
					#						die "here $$ph{'txe'},$diff, ";
					getseq($fh,$$ph{'txe'},$diff,$offset,$es,\%lflank,\@cons);

					foreach my $s (@cons) {
						$ts=substr($$maf{$s},$pos-$flankl,$flankl);
						if ($strand eq '-') { ## we are in antisense orientation from given maf extraction
							$ts=rc($ts);
							$lflank{$s} = $ts.$lflank{$s};
							#die $lflank{$s}; # works this is actually the right flank of the as transcript on minus strand, so next codon left to it must be stop
						} else {
							$lflank{$s} = $lflank{$s}.$ts;
						}
					}
				}
				## I think it is correct now, lets go to the right flank sequence now
			} else {
				$flankl-=($flankl%3);
				foreach my $s (@cons) {	## also check here for the strands	
					$lflank{$s} = substr($$maf{$s},$pos-$flankl,$flankl);
				}
			}
		}
	}
	
	

	#print "lflank $$ph{'name'},",%lflank,"\n";
	# print "$id\n"; ## deletehere
	#get_rate(\%lflank,\%lambda,\@cons,'l',$id) if(%lflank);
	#print "$lflank{$cons[0]}\n";  ## deletehere
	#			print "\n\nnext\n\n";
######################### right flank ## this is only in case we are in genome mode or genomic plus, otherwise we need to check if we get smaller 0, when transcript is on minus
	my $flankr=$orf_len;
	if($pos+$orf_len+$flankr > $finalpos){
	  $flankr = $finalpos-$pos-$orf_len;
	  $flankr-=($flankr %3); ## make seq a multiple of 3 , really important otherwise it will crash
	}
	
	$diff = $orf_len-$flankr;
	
	if(not $options{'b'}){
	  getseq($fh,$pos+$orf_len,$flankr,$offset,$strand,\%rflank,\@cons);
	}else{
		if ($$ph{'txs'} < $orf_len or $$ph{'txe'}+$orf_len > $length_chr) {
		} else {
			if($diff > 0){
				#die "$strand  $$ph{'strand'}\n";
				if($strand eq '+' and $$ph{'strand'} eq '-'){
					$es='-';
				}elsif($strand eq '-' and $$ph{'strand'} eq '+'){
					$es='-';
				}
		
		## now we get what is genomically downstream of the stop codon
		if($$ph{'strand'} eq '+'){
		  getseq($fh,$$ph{'txe'},$diff,$offset,$es,\%rflank,\@cons);
		  foreach my $s(@cons){
			$ts=substr($$maf{$s},$pos+$orf_len,$flankr);
			if($strand eq '-'){ ## we are in antisense orientation from given maf extraction
			  $ts=rc($ts);
			  $rflank{$s} = $rflank{$s}.$ts; ## this is not correct yet since we are on the other side actually so it should be right flank than
			  #die $rflank{$s}; #workds ## we are genomic minus, so next codon left is ATG
			}else{
			  $rflank{$s} = $ts.$rflank{$s};
			  #die $rflank{$s}; #works we are genomic plus so next codon left is a stop
			  #die "$pos  $orf_len  $flankl $diff  here $lflank{$s}\n";
			}
		  }
		  
		}else{ ## transcript is on genomic minus
		  # if($$ph{'txe'} + $diff > $length_chr){
		  # 	$diff = $length_chr-$$ph{'txe'};
		  # 	$diff-=($diff %3); ## make seq a multiple of 3 , really important otherwise it will crash
		  # }
		  
		  ## we can get the problem that the right flank goes below 0 so we need to check this first
		  if($$ph{'txs'} - $diff < 0){
			$diff=$$ph{'txs'};
			## now we need to make sure that the sequence is still a multiple of 3 
			my $mod=(($flankr+$diff) %3);
			$diff-=$mod;
		  }
		  
		  ## now we get what is genomically upstream of the stop codon since we are on minus
		  getseq($fh,$$ph{'txs'}-$diff,$diff,$offset,$es,\%rflank,\@cons);
		  foreach my $s(@cons){
			$ts=substr($$maf{$s},$pos+$orf_len,$flankr);
			if($strand eq '-'){ ## we are in antisense orientation from given maf extraction
			  $ts=rc($ts);
			  $rflank{$s} = $rflank{$s}.$ts;
			  #die $rflank{$s}; #works  ## here we are on genomic plus
			}else{
			  $rflank{$s} = $ts.$rflank{$s};
			  #die $rflank{$s}; #works we are on genomic minus
			  #								die "$pos  $orf_len  $flankl $diff  here $lflank{$s}\n";
			}
		  }
		}
	  }else{	
		$flankr-=($flankr %3); ## make it a multiple of 3
		foreach my $s(@cons){
		  $rflank{$s} = substr($$maf{$s},$pos+$orf_len,$flankr);
		  #$rflank{$s} -= ($rflank{$s} %3); ## make it a multiple of 3
		}
	  }
		}
	}
	
	
	#print "$diff rflank $$ph{'name'},",%rflank,"\n";
	#get_rate(\%rflank,\%lambda,\@cons,'r',$id) if(%rflank);
	#			print "\n\nnext\n\n";			
	
	
	

	#%mafp; ## we need to fill this with the sequences we have to have
	foreach my $s(@cons){
	  #$mafp{$s} = substr($$maf{$s},$pos,$orf_len);
	  $mafp{$s}=$$tmpseq{$s};
	  ## on genome run we have the correct orientation already
	  ## on bed files we are always on sense by default, so if sth is antisense we need to reverse it here for the maf, flanks are automatically correct! 
	  $mafp{$s} = rc($mafp{$s}) if($strand eq '-' and $options{'b'}); ## rc if strand eq minus so we are anti sense, 
	}
			# foreach my $k(keys %mafp){
			# 	print "$k\n$mafp{$k}\n";
			# }

			## dnds for orf

#	print "$pos\t$e\tchr8:$$ph{$pos}{'abs'}-$$ph{$e}{'abs'}\t$strand\n";
	my %kick=();
	if($options{'K'}){
#		die "$cp{$ph{1543}{'abs'}} $pos $ph{'cdssr'}\n" if($pos == 40);
		check_codon_pos($cp,$pos,$e,$$ph{'strand'},$strand,\%kick,$ph);
		# if($pos == 400){
		# 	for my $k(sort {$b <=> $a} keys %kick){
		# 		print "$k,";
		# 	}
		# 	print "\n";
		# }
	}
	
	my %phylocsf_score=();
	$phylocsf_score{'orf'}="na";
	$phylocsf_score{'flankl'}="na";
	$phylocsf_score{'flankr'}="na";
	$phylocsf_score{'cds'}="na";
	
	$phylocsf_score{'orf'}=kaks(\%mafp,\@cons,$orf_len,$kaks_table,\%kaks,$id,'orf',$ph,\%kick,$pos,1,\$codons_scored);
#	print $codons_scored,"\n";
	dnds(\%lambda,\%kaks,\%dnds,\@cons,$id,'orf',$orf_len);   # getseq($fh,$start,$len,$offset,$strand);

	## hack - not very elegant but works
	my %kaks_iss;
	my %dnds_iss;
	my %mafp_iss;
	## remove start and end from score
	foreach my $s(keys %mafp){
		$mafp_iss{$s} = substr($mafp{$s},3,$orf_len-6);
	}

	my $jk1;
	if($orf_len > 8){
		$jk1=kaks(\%mafp_iss,\@cons,$orf_len-6,$kaks_table,\%kaks_iss,$id,'orf',$ph,\%kick,$pos,1,\$codons_scored);
    #	print $codons_scored,"\n";
		dnds(\%lambda,\%kaks_iss,\%dnds_iss,\@cons,$id,'orf',$orf_len-6);   # getseq($fh,$start,$len,$offset,$strand);
	}else{
		$dnds_iss{$name}{'orf'}{'median'} = -1000;
	}


	my $dmp;
	if(%lflank){
		$phylocsf_score{'flankl'}=kaks(\%lflank,\@cons,$flankl,$kaks_table,\%kaks,$id,'flankl',$ph,\%kick,$pos,	$phylocsf_score{'orf'},\$dmp);
		dnds(\%lambda,\%kaks,\%dnds,\@cons,$id,'flankl',$flankl);   # getseq($fh,$start,$len,$offset,$strand);
	}else{
		$dnds{$name}{'flankl'}{'median'}="na";
		$dnds{$name}{'flankl'}{'mean'}="na";
		$dnds{$name}{'flankl'}{'sd'}="na";
	}
	
	## this is the right flank
	if(%rflank){
		$phylocsf_score{'flankr'}=kaks(\%rflank,\@cons,$flankr,$kaks_table,\%kaks,$id,'flankr',$ph,\%kick,$pos,	$phylocsf_score{'orf'},\$dmp);
		dnds(\%lambda,\%kaks,\%dnds,\@cons,$id,'flankr',$flankr);   # getseq($fh,$start,$len,$offset,$strand);
	}else{
		$dnds{$name}{'flankr'}{'median'}="na";
		$dnds{$name}{'flankr'}{'mean'}="na";
		$dnds{$name}{'flankr'}{'sd'}="na";
	}


	my $dndsoutstr;
	my $rels;
	my $rele;
	my $ostr='+';
	my $as='';
	
	my $fragS=1;
	my $fragE=1;
	my $junctionC=0;
	if($options{'b'}){
        ## this we need for the bed file downstream but we need to get the correct coords for the id2		
		$rels=$pos;
		$rele=$e;
				

		
		if($options{'K'}){ ## in circ mode
			while($rels -$ph{'circl'} > 0){
				$fragS++;
				$rels-=$ph{'circl'};
				$rele-=$ph{'circl'};
				$fragE++;
			}
			while($rele -$ph{'circl'} > 0){
				$rele-=$ph{'circl'};
				$fragE++;
			}

			$rels+=$ph{'cdssr'}; ## add the offset, now we are at correct genomic position again
			$rele+=$ph{'cdssr'};
			$junctionC=$fragE-$fragS;
			$totalC++ if($junctionC);
		}
		


		if($$ph{'strand'} eq '+'){
			if($strand eq '-'){$ostr='-'; $as='.as'}
			
		}else{
#			$rels=$$ph{'len'}-$e-1; ## is correct
#			$rele=$$ph{'len'}-$pos-1; ## is correct
			if($strand eq '+'){$ostr='-';}else{$as='.as';}
			#die "$e $pos $rels $rele";
		}
#		print "$pos $e\n";
		my $id2;


		if($options{'K'}){
			my ($jk1,$jk2);
			$jk1 = $$ph{$rels}{'abs'};
			$jk2 = $$ph{$rele}{'abs'};
			if($jk1 < $jk2){
				$id2="\t$options{'c'}:".($jk1+1)."-".($jk2+1).":$ostr";
			}else{
				$id2="\t$options{'c'}:".($jk2+1)."-".($jk1+1).":$ostr";
			}
		}else{
			my ($jk1,$jk2);
			$jk1 = $$ph{$pos}{'abs'};
			$jk2 = $$ph{$e}{'abs'};
			if($jk1 < $jk2){
				$id2="\t$options{'c'}:".($jk1+1)."-".($jk2+1).":$ostr";
			}else{
				$id2="\t$options{'c'}:".($jk2+1)."-".($jk1+1).":$ostr";
			}
		}

		$dndsoutstr.=$id2;

	}else{
		my $id2="${chr}:".($pos+1)."-".($pos+$orf_len).":$strand";
		$dndsoutstr.="\t$id2";
	}

	$dndsoutstr.="\t$dnds{$name}{'orf'}{'median'}\t$dnds_iss{$name}{'orf'}{'median'}";
	my $dnm=$dndsoutstr."\t$dnds{$name}{'flankr'}{'median'}\t$dnds{$name}{'flankl'}{'median'}\t$dnds{$name}{'orf'}{'mean'}\t$dnds{$name}{'flankr'}{'mean'}\t$dnds{$name}{'flankl'}{'mean'}\t$dnds{$name}{'orf'}{'sd'}\t$dnds{$name}{'flankr'}{'sd'}\t$dnds{$name}{'flankl'}{'sd'}";
	my $dnp=$dndsoutstr."\t$dnds{$name}{'flankl'}{'median'}\t$dnds{$name}{'flankr'}{'median'}\t$dnds{$name}{'orf'}{'mean'}\t$dnds{$name}{'flankl'}{'mean'}\t$dnds{$name}{'flankr'}{'mean'}\t$dnds{$name}{'orf'}{'sd'}\t$dnds{$name}{'flankl'}{'sd'}\t$dnds{$name}{'flankr'}{'sd'}";
	
	
	## get correct outstring for dnds
	if(not $options{'b'}){
		if($strand eq '-'){
#			$dndsoutstr=$dnm;
			$dndsoutstr=$dnm;
		}else{
			$dndsoutstr=$dnp;
		}
#            print "\t$name\t$dnds{$name}{'orf'}{'median'}\t$dnds{$name}{'flankl'}{'median'}\t$dnds{$name}{'flankr'}{'median'}\n";
	}else{
		$dndsoutstr=$dnp;  ## corrected above already when getting flanks
		# if($$ph{'strand'} eq '+' and $strand eq '-'){
		# 	$dndsoutstr=$dnm;
			
		# }elsif($$ph{'strand'} eq '-' and $strand eq '+'){
		# 	$dndsoutstr=$dnm;
		# }else{
		# 	$dndsoutstr=$dnp;
		# }
	}
	

	
	
	
	
	
	
	my $frameshiftIoutstr="\t".scalar @frameshiftI;
	if(scalar @frameshiftI == 0){
		$frameshiftIoutstr.="\t-";
	}else{
		$frameshiftIoutstr.="\t".join(",",@frameshiftI);
	}

	$overlap_with=overlap_wrapper($chr,$strand,$pos,$e,$refseq,$nc_refseq,$oa_files,$igs,\@cons,$orf_len,$syn_mut,$aa,$dna,$dnaH,$oe,$pe,$se,$synmutstat,$stringATG,$dndsoutstr,$frameshiftIoutstr,$tmpseq);
	return if($overlap_with eq 1); ## that means we printed already something
	
	if($options{'b'}){
		print IG "$$ph{'name'}\t$$ph{'strand'}\t$strand";
	}else{
		print IG "$chr\t$strand\t+";
	}
	
	print IG "\t$pos\t$e\t$orf_len\t$frame\t",scalar @cons,"\t$syn_mut\t$aa\t$dna\t$dnaH\t$oe\t$pe\t$se\t",join(",",@cons),"\t$kozak\t$shined\t$dtss\t$dtes\tdsc\tdec\t$overlap_with",getCodonMut($tmpseq,\@cons,0,$orf_len,$strand),$synmutstat,$stringATG,$dndsoutstr,$frameshiftIoutstr;
	
	## determine correct relative position in orf
	if($options{'b'} and not $options{'K'}){ ## only if we are not in circle mode, otherwise we get stupid results for cds
		## get exons now in which our micpep sits .... 
		my %eh=();
		for(my $i = $rels;$i<=$rele;$i++){
			$eh{$$ph{$i}{'pos'}}++;
		}
	   
			# my @bsa;
			# my $bs;
		# my $bl;
		# my $bc;
		# for my $k(sort {$a<=>$b} keys %eh){
		# 	push(@bsa,$k);
		# 	$bl.="$eh{$k},";
		# 	$bc++;
			# }
		# $bsa[0]+=$eh{$$ph{$rels}{'pos'}};
		# 			 #			chop $bs;
		# 			 #			chop $bl;
		# 			 $bs=join(",",@bsa).",";
		
			## add block information if necessary
#			print  "$options{'c'}\t",$$ph{$rels}{'abs'},"\t",$$ph{$rele}{'abs'}+1,"\t$$ph{'name'}$as\t$dnds{$name}{'orf'}{'median'}\t$ostr\t$$ph{$rels}{'abs'}\t",$$ph{$rele}{'abs'}+1,"\t$bc\t$bc\t$bl\t$bs\n";
#			print  "$options{'c'}\t",$$ph{'txs'},"\t",$$ph{'txe'},"\t$$ph{'name'}$as\t$dnds{$name}{'orf'}{'median'}\t$ostr\t$$ph{$rels}{'abs'}\t",$$ph{$rele}{'abs'}+1,"\t$bc\t$bc\t$bl\t$bs\n";
#			print  BEDI "$options{'c'}\t",$$ph{'txs'},"\t",$$ph{'txe'},"\t$$ph{'name'}$as\t$dnds{$name}{'orf'}{'median'}\t$ostr\t$$ph{$rels}{'abs'}\t",$$ph{$rele}{'abs'}+1,"\t$bc\t$bc\t$bl\t$bs\n";
#			print BEDI "$options{'c'}\t",$$ph{$rels}{'abs'},"\t",$$ph{$rele}{'abs'}+1,"\t$$ph{'name'}$as\t$dnds{$name}{'orf'}{'median'}\t$ostr\t$$ph{$rels}{'abs'}\t",$$ph{$rele}{'abs'}+1,"\n";
			## plus one since our end coord is included but in bed it must be excluded, we are 0 based already
		
		my $type;
		#print "$rels:$rele:$$ph{$rels}{'t'}:$$ph{$rele}{'t'}\n";
		$type="$$ph{$rels}{'t'}:$$ph{$rele}{'t'}";
		

		if($strand eq '-'){
			$type="$type:AS";
		}else{
			$type="$type:S";
		}

## old type def 
		# $type='utr5' if($$ph{$rels}{'abs'} < $$ph{'cdss'} and $strand eq '+' and $$ph{'strand'} eq '+');
		# $type='utr5as' if($$ph{$rels}{'abs'} < $$ph{'cdss'} and $strand eq '-' and $$ph{'strand'} eq '+');
		# $type='utr3' if($$ph{$rels}{'abs'} >= $$ph{'cdse'} and $strand eq '+' and $$ph{'strand'} eq '+');
		# $type='utr3as' if($$ph{$rels}{'abs'} >= $$ph{'cdse'} and $strand eq '-' and $$ph{'strand'} eq '+');
		
		# $type='utr5' if($$ph{$rels}{'abs'} >= $$ph{'cdse'} and $strand eq '+' and $$ph{'strand'} eq '-');
		# $type='utr5as' if($$ph{$rels}{'abs'} >= $$ph{'cdse'} and $strand eq '-' and $$ph{'strand'} eq '-');   
		# $type='utr3' if($$ph{$rels}{'abs'} < $$ph{'cdss'} and $strand eq '+' and $$ph{'strand'} eq '-');
		# $type='utr3as' if($$ph{$rels}{'abs'} < $$ph{'cdss'} and $strand eq '-' and $$ph{'strand'} eq '-'); 
		
		
			
		# ## this is a true cds containing thing
		# if($$ph{$rels}{'abs'} >= $$ph{'cdss'} and $$ph{$rels}{'abs'} < $$ph{'cdse'}){
		# 	if($strand eq '+'){
		# 		$type='cds';
		# 	}else{
		# 		$type='cdsas';
		# 	}
		# }

		# if($$ph{'cdss'} == $$ph{'cdse'}){
		# 	$type='nc' if($strand eq '+'); 
		# 	$type='ncas' if($strand eq '-');
		# }
####################		
		## now we get the frame information here
		my ($frame,$framex);
		
		$frame = abs($rels-$$ph{'cdssr'}) %3; ### since transcript cssr is always relative in correct orientation the frames are the same for plus and minus as transcripts 
#		print "$rels $rels-$$ph{'cdssr'}\n";
#			die "$$ph{'name'}\t$pos\t$e\t$$ph{'len'}\t$strand\t$$ph{'strand'}\n" if(not defined $$ph{'cdssr'});
#
		if($type =~ /cds:S/){
			if($frame == 0){
				$framex=0;
			}else{
				$framex=$frame;
			}

		}elsif($type =~ /cds:AS/){
			#if ($$ph{'strand'} eq '+') {
			if ($frame == 2) { ## if orig is on plus and we are as then modulo needs to be two 
				$framex=2;
			}else{
				$framex=$frame;
			}
			
		}else{
			$framex ='';
		}
#			die "$$ph{'name'}\t$frame\t($rels-$$ph{'cdssr'})\n" if(not defined $framex);
		
		if($framex ne ''){
			$type.=":$framex";
		}
			 
		#print "$frame $pos $mafp{$cons[0]}\n";
		## we want to check if the dnds we have is better than the one in the coding sequence
		%kaks=();
		my %dndsc=();
		%lambda=();
		%mafp=();
		
		### old we dont need it like this
		## we get a problem here if we get out of range, so we reduce by three which should not make a big difference in dnds calc!!! but we need to check this carefully
		my $st=0;
		
		## we add the frame instead of subtracting it if we are above the finalpos;
		my $posr=$pos;
		if($$ph{'strand'} eq '-' and ($pos-$frame+$orf_len) > $finalpos){
			$frame=abs($frame-3);
			$st=3;
			
		}
		
		
		## if we are below 0 lets add 3 to pos so we miss one codon in beginning but get one from the end
		if($pos-$frame < 0){
			$posr+=3;
			$st=3;
		}


		## we always get sequence from sense of transcript !

		foreach my $s(@cons){ ## orig maf			
			if($strand eq '+' and $$ph{'strand'} eq '+'){                   ## sense
				$mafp{$s} = substr($$maf{$s},$posr-$frame,$orf_len);
			}elsif($strand eq '+' and $$ph{'strand'} eq '-'){               ## sense
				#die "$pos $frame $$maf{$s}";
				$mafp{$s} = substr($$maf{$s},$posr-$frame,$orf_len);


			}elsif($strand eq '-' and $$ph{'strand'} eq '+'){      ## antisense
				$mafp{$s} = substr($$maf{$s},$posr-$frame,$orf_len-$st);
			}elsif($strand eq '-' and $$ph{'strand'} eq '-'){      ## antisense
				$mafp{$s} = substr($$maf{$s},$posr-$frame,$orf_len-$st);
			}
		}
		
#		print $mafp{$species[0]},"\n";
		$phylocsf_score{'cds'}=kaks(\%mafp,\@cons,$orf_len-$st,$kaks_table,\%kaks,$id,'cds',$ph,\%kick,$pos,"na-cod",\$dmp);
		dnds(\%lambda,\%kaks,\%dndsc,\@cons,$id,'cds',$orf_len-$st);   # getseq($fh,$start,$len,$offset,$strand);
		
		
		#		print  "$options{'c'}\t",$$ph{'txs'},"\t",$$ph{'txe'},"\t$$ph{'name'}.$type\t$dnds{$name}{'cds'}{'median'}\t$ostr\t$$ph{$rels}{'abs'}\t",$$ph{$rele}{'abs'}+1,"\t$dndsc{$name}{'cds'}{'median'}\t$$ph{'bc'}\t$$ph{'bl'}\t$$ph{'bs'}\n";
		#print  BEDI "$options{'c'}\t",$$ph{'txs'},"\t",$$ph{'txe'},"\t$$ph{'name'}.$type\t$dnds{$name}{'cds'}{'median'}\t$ostr\t$$ph{$rels}{'abs'}\t",$$ph{$rele}{'abs'}+1,"\t$dndsc{$name}{'cds'}{'median'}\t$$ph{'bc'}\t$$ph{'bl'}\t$$ph{'bs'}\n";
		print  BEDI "$options{'c'}\t",$$ph{'txs'},"\t",$$ph{'txe'},"\t$$ph{'name'}.$type\t$dnds{$name}{'orf'}{'median'}\t$ostr\t$$ph{$rels}{'abs'}\t",$$ph{$rele}{'abs'}+1,"\t$dndsc{$name}{'cds'}{'median'}\t$$ph{'bc'}\t$$ph{'bl'}\t$$ph{'bs'}\n";
		
		## print to exons of peptide only, necessary for pdf and pvalue calculation
#		print BEDIE "$options{'c'}\t",$$ph{'txs'},"\t",$$ph{'txe'},"\t$$ph{'name'}\t$dnds{$name}{'cds'}{'median'}\t$ostr\t$$ph{$rels}{'abs'}\t",$$ph{$rele}{'abs'}+1,"\t$dndsc{$name}{'cds'}{'median'}\t$$ph{'bc'}\t$$ph{'bl'}\t$$ph{'bs'}\n";
				
		print IG "\t$type\t$dndsc{$name}{'cds'}{'median'}";
	}

	if($options{'K'}){
		print IG "\t$$ph{'circl'}\t$fragS\t$fragE\t$junctionC";
	}

	if($options{'X'}){
		print IG "\t$phylocsf_score{'orf'}\t$phylocsf_score{'flankl'}\t$phylocsf_score{'flankr'}\t$phylocsf_score{'cds'}\t$codons_scored";
	}

	
	print IG "\t$nonsense";
	print IG "\n";
	return;		    
} ## here it should end

sub dnds_old{
    my ($lambda,$kaks,$dnds,$species,$name,$type,$len)= @_;

    my $sum=0;
    my $L;

    my $sum1=0;
    my $sum2=0;
    my $sumt=0;
    my $sumt2=0;

    my @sorted;
    my $mid;

    $L=$len;
    my %e=();

    foreach my $s(@$species){
        next if($s eq $$species[0]); ## without refspec
#           $L=$ids{$name}{'l'};
         ## size of flanks used to calc lamdas

		$sum1= $$kaks{$name}{$s}{$type}{'ka'};
		$sum2= $$kaks{$name}{$s}{$type}{'ks'};

        $e{$s} = $sum1/$sum2;
        $sumt+= $e{$s};
    }

    @sorted=sort { $a <=> $b }values %e;

    $mid = $#sorted/2;  ## this is the number of elements -1 so if its mod 2 ==0 then we have an odd number of elements otherwise not
    if($mid % 2 == 0){
        $$dnds{$name}{$type}{'median'}=$sorted[ $mid ];
    }else{ ## tied ranks
        $mid=int($mid);
        $$dnds{$name}{$type}{'median'}=(($sorted[$mid] + $sorted[$mid+1])/2);
    }


    ## without refspec
	if(scalar @$species ==1){ ## why does this happen ????
		$$dnds{$name}{$type}{'median'} =-1;
		$$dnds{$name}{$type}{'mean'} =-1;
		$$dnds{$name}{$type}{'var'} = -1;
        $$dnds{$name}{$type}{'sd'} =-1;
		print STDERR  @$species,"\t$name\n";
		return;
	}
    $$dnds{$name}{$type}{'mean'} = $sumt/(scalar @$species -1);
    my $var2=0;
    foreach my $s(@$species){
		next if($s eq $$species[0]); ## without refspec
		$var2+= ($e{$s}-$$dnds{$name}{$type}{'mean'})**2
    }

    if(scalar @$species ==2){
        $$dnds{$name}{$type}{'var'} = 0;
        $$dnds{$name}{$type}{'sd'} =0;
        return;
    }
    ## we need at least 3 species otherwise var cannot be calculated, one spec is the reference ........

    $var2/=(scalar @$species -2); ## corresponds to n-1 since we used already the mean and thus loose one degree of freedom

    #$dnds{$name}{$type}{'var'} = ($sumt2/(scalar @$species)) - ($sumt/(scalar @$species))**2 ;
    $$dnds{$name}{$type}{'var'} = $var2;

    ## somehow we got ausloeschung since numbers are small for machine formula, other
    #print "$dnds{$name}{$type}{'var'}\t$dnds{$name}{$type}{'var2'}\n";
    #if($dnds{$name}{$type}{'var'} <0){
    #   die "$name\t$dnds{$name}{$type}{'var'}\t$var2\n";
    #}
    $$dnds{$name}{$type}{'sd'} = sqrt( $var2 );# if($dnds{$name}{$type}{'var'} > 0);
}


sub getrate{
    my ($maf,$lambda,$species,$side,$name)=@_;
    my $mut=0;

    my $refseq;
    my $seq;
    my $codref;
    my $codspec;

	
	my $len=1;
	

    foreach my $s(@$species){
        $mut=0;
        if($s eq $$species[0]){
            $refseq=uc($$maf{$s});
			$len=length($refseq);
			$len=1 if(not $len);
		}else{
            $seq=uc($$maf{$s});
            for(my $i = 0;$i < length($refseq);$i++){
                if(substr($refseq,$i,1) ne substr($seq,$i,1)){
                    $mut++;
                }
            }
        }
        ## now we get the lamdas per side
        if(exists  $$lambda{$name}{$s}){
            $$lambda{$name}{$s}{'m'}= ($$lambda{$name}{$s}{'l'}+($mut/$len) )/2;
        }
        $$lambda{$name}{$s}{$side}+=$mut/$len;
    }
}


sub hamming {
    ## strings are compared in binary format with xor, then we simply count the remaining zeros
	return 1 if($options{'y'});
    return ($_[0] ^ $_[1]) =~ tr/\0//c;
}



sub kaks{
    my ($maf,$species,$len,$kaks_table,$pairwise_kaks_result,$name,$type,$ph,$kick,$rels,$orfscore,$codons_scored) = @_;
#   my %delete;

    my ($syn,$nsyn,$usyn,$unsyn,$ka,$ks)= (0,0,0,0,0,0);
    my $refs=uc($$maf{$$species[0]});

	## sanity check for the length
	$len=length($refs); ## set the length here again since we have sometimes bad behavior

	## this can only happen for boundary regions
	my $diff=length($refs) %3;

	if($diff != 0){
		$len-=$diff;

	}
	die "length < 0 for $refs   $name \n" if($len <0);

    my %uniq=();
    my %pairs=();
    my ($cod1,$cod2);
    my $hd=0;

	my $subtractlen=0;
	my %scodignore=();

	my $kicks;
	my $codk=0;
	my %ignore;

    foreach my $s (@$species){
        $pairs{$s}{$type}{'s'}{'h'}=0;
        $pairs{$s}{$type}{'n'}{'h'}=0;
		#print "$refs,$s,@$species\n";
		$subtractlen=0; ## how many codons to discard ?
		$scodignore{$s}{'l'}=0;

		next if($s eq $species[0]);

        for (my $i=0+$len%3; $i<$len;$i+=3){
            $cod1=uc(substr($refs,$i,3));
            $cod2=uc(substr($$maf{$s},$i,3));

			## we ignore codons with gaps
			$codk=0;
			if($type eq 'orf' and $options{'B'}){
				$kicks=$rels+$i;
				while($kicks >= $$ph{'circl'}){
					$kicks -= $$ph{'circl'};
				}	

				$kicks+=$ph{'cdssr'};

#				die $$ph{$kicks+3}{'abs'};

				if(exists $$ph{$kicks}){
					if(exists $$kick{$$ph{$kicks}{'abs'}}){
						$ignore{$i}=1;
						$scodignore{$s}{'ignorepos'}{$i}=1; ## put codon position to ignore in hash
						$scodignore{$s}{'l'}+=3;            ## remove 3 from length in the end
						$codk=1;
					}
				}
				
				if($codk == 0){
					my $k3=$kicks+3;
					if(exists $$ph{$k3}){
						if(exists $$kick{$$ph{$k3}{'abs'}}){
							$ignore{$i}=1;
							$scodignore{$s}{'ignorepos'}{$i}=1; ## put codon position to ignore in hash
							$scodignore{$s}{'l'}+=3;            ## remove 3 from length in the end
							$codk=1;
						}
					}
				}
			}

			
			## we ignore codons with gaps			
			if($options{'Y'}){
				if($cod2 =~ /\-/){
					if(not $codk){
						$scodignore{$s}{'ignorepos'}{$i}=1; ## put codon position to ignore in hash
						$scodignore{$s}{'l'}+=3;            ## remove 3 from length in the end
					}
					$codk=1;
				}
			}
			next if($codk ==1 ); 

            ## if codons are different
            if($cod2 ne $cod1){
                ## if AA are the same then we have a synonymous mutation
                $hd=hamming($cod1,$cod2);

                if(peptide($cod1) eq peptide($cod2)){
                    $syn++;
                    $uniq{$i}{'s'}{$cod2}++;
                    $pairs{$s}{$type}{'s'}{$i}=$hd;
                    $pairs{$s}{$type}{'s'}{'h'}+=$hd;  ## number of changed nt in codon for substitution
                    $pairs{$s}{$type}{'s'}{'n'}++;     # number of substitutions
                }else{
                    $nsyn++;
                    $uniq{$i}{'n'}{$cod2}++;
                    $pairs{$s}{$type}{'n'}{$i}=$hd;
                    $pairs{$s}{$type}{'n'}{'h'}+=$hd; ## number of changed nt in codon for substitution
                    $pairs{$s}{$type}{'n'}{'n'}++;    ## number of substitutions
                }
            }
        }

		if($options{'Y'}){
			## init sites
			$scodignore{$s}{'synsites'}=0;
			$scodignore{$s}{'nonsynsites'}=0;
			$scodignore{$s}{'len'}=$len-$scodignore{$s}{'l'};    ## remove 3 from length in the end
		#	print "$s $scodignore{$s}{'len'}   $refs\n";
		}
		
    }

    ## here we do ka and ks calculation
    my $synsites=0;
    my $nonsynsites=0;

    my $cod;
    ## first we get the number of synonymous and nonsynonymous sites
    for(my $i=0+$len%3;$i<$len;$i+=3){
        $cod=substr($refs,$i,3);

#		next if($cod !~ /[ACGT][ACGT][ACGT]/);
		next if($cod =~ /N/);
		if(not defined $$kaks_table{$cod}{1} or not defined $$kaks_table{$cod}{2} or not defined $$kaks_table{$cod}{3}){
			print STDERR  "error missing sequence $refs,$i,$cod\t$name\t$$ph{'name'}\n";
			return 1000;
		}
        $synsites+=$$kaks_table{$cod}{1}+$$kaks_table{$cod}{2}+$$kaks_table{$cod}{3};

		## foreach species we count the synsites separately	
		if($options{'Y'}){
			foreach my $s(@$species){
				next if($s eq $species[0]);
				if(not $scodignore{$s}{'ignorepos'}{$i}){
					$scodignore{$s}{'synsites'}+=$$kaks_table{$cod}{1}+$$kaks_table{$cod}{2}+$$kaks_table{$cod}{3};
				}
			}
		}
    }
	
    $nonsynsites=$len-$synsites;
	if($options{'Y'}){
		## now species specific nonsynsites
		foreach my $s(@$species){
			next if($s eq $species[0]);
			if($scodignore{$s}{'len'} >0){
				$scodignore{$s}{'nonsynsites'}=$scodignore{$s}{'len'}-$scodignore{$s}{'synsites'};						
				$scodignore{$s}{'addsyn'}=$scodignore{$s}{'synsites'}/$scodignore{$s}{'len'};
				$scodignore{$s}{'addnsyn'}=$scodignore{$s}{'nonsynsites'}/$scodignore{$s}{'len'};

				$scodignore{$s}{'synsites'}+=$scodignore{$s}{'addsyn'};
				$scodignore{$s}{'nonsynsites'}+=$scodignore{$s}{'addnsyn'};
			}else{ ## so we have 0 sequence to score then we set it to 1 
				$scodignore{$s}{'addnsyn'}=1; ## else we set it to -1 for undefined
				$scodignore{$s}{'addsyn'}=1;  ## why does it become negative with my new scoring
			}
		}
	}


	my $addsyn=$synsites/$len;
	my $addnsyn=$nonsynsites/$len;
	
	$synsites+=$addsyn;
	$nonsynsites+=$addnsyn;


    ## here comes now the tricky part
    ## this is Ka Ks with hamming distance
    foreach my $s (@$species){
        next if($s eq $$species[0]);
        $$pairwise_kaks_result{$name}{$s}{$type}{'ka'} = $pairs{$s}{$type}{'n'}{'h'};

		if($options{'Y'}){
			$$pairwise_kaks_result{$name}{$s}{$type}{'ka'}+=$scodignore{$s}{'addnsyn'};
			if($scodignore{$s}{'nonsynsites'}){
				$$pairwise_kaks_result{$name}{$s}{$type}{'ka'} /= $scodignore{$s}{'nonsynsites'};
			}else{
				$$pairwise_kaks_result{$name}{$s}{$type}{'ka'} = 1; ## we have 0 mutations
			}
		}else{
			$$pairwise_kaks_result{$name}{$s}{$type}{'ka'}+=$addnsyn;
			if($nonsynsites){
				$$pairwise_kaks_result{$name}{$s}{$type}{'ka'} /= $nonsynsites;
			}else{  ## this should never happen since we added a pseudocount already above
				$$pairwise_kaks_result{$name}{$s}{$type}{'ka'} =1;  ## we have 0 mutations
			}
		}
		
        $$pairwise_kaks_result{$name}{$s}{$type}{'ks'} = $pairs{$s}{$type}{'s'}{'h'};
		

		## $$pairwise_kaks_result{$name}{$s}{$type}{'ks'}++ if($options{'C'});
		## we just add one to the counts
		#if($options{'C'}){
		if($options{'Y'}){
			$$pairwise_kaks_result{$name}{$s}{$type}{'ks'}+=$scodignore{$s}{'addsyn'};
			if($scodignore{$s}{'synsites'}){ ## can be 0 if length is 0
				$$pairwise_kaks_result{$name}{$s}{$type}{'ks'} /= $scodignore{$s}{'synsites'};
			}else{
				$$pairwise_kaks_result{$name}{$s}{$type}{'ks'} =1;
			}
		}else{
			$$pairwise_kaks_result{$name}{$s}{$type}{'ks'}+=$addsyn;
			if($synsites){
				$$pairwise_kaks_result{$name}{$s}{$type}{'ks'} /= $synsites;
			}else{ ## this should never happen since we added a pseudocount already above
				$$pairwise_kaks_result{$name}{$s}{$type}{'ks'} =1;
			}
		}	
	}

	## if we have options X we also score with PhyloCSF omega mode
	if($options{'X'}){
		my $ocod=0;
		print STDERR "scoring phylocsf $name $type \n" if($options{'v'});
		open PHY,">$options{'X'}/$currentgene${name}_$type.maf" or die "File $options{'X'}/$currentgene${name}_$type.maf could not be created\n";
		foreach my $s (@$species){
			print PHY ">$s";
			if($s eq $$species[0]){print PHY "|$currentgene${name}";}
			print PHY "\n";
			for (my $i=0+$len%3; $i<$len;$i+=3){
#				print "$i\t",uc(substr($$maf{$s},$i,3)),"\n";
				next if($ignore{$i});	
				$ocod++ if($s eq $$species[0]);
				print PHY uc(substr($$maf{$s},$i,3));
			}
			$$codons_scored=$ocod if($s eq $$species[0]);
			print PHY "\n";
		}
		close PHY;
		my $sss=join(",",@$species);
		if($ocod == 0){
			print STDERR "$currentgene${name}_$type.maf aborted. NO codons output for alignment since all positions were skipped\n";
			return "na-cod";
		}
		
		if($ocod < $mincodons){
			my $skipped=scalar keys %ignore;
			print STDERR "$currentgene${name}_$type.maf aborted. Not enough codons output for alignment to score. $mincodons needed, $ocod available, $skipped skipped to coding gene ambiguity\n";
			return "na-min";
		}

		## we skip the flank scoring with phyloCSF when option F is given
		if($options{'F'} and $type =~ /flank/){
			return "na";
		}

		if($options{'D'}){
			return "na";
		}else{

			if($orfscore !~ /cod/){
				my $res=`$phylo $param_phylo $options{'X'}/$currentgene${name}_$type.maf --strategy=omega --species=$sss`;
				print STDERR " done \n" if($options{'v'});
				if($res =~ /(\S+)$/){
					return $1;
				}else{
					return "na";
				}
			}else{
				print STDERR "aborted cause ORF had no codons output for alignment since all positions were skipped\n";
				return "na";
			}
		}
	}
	return "na";
}


sub check_stop_and_start_codon_conservation{
	my ($maf,$species,$atgstopcons,$pos,$len,$strand) = @_;
	my ($start,$stop,$refstop);

	foreach my $s(@$species){
		#if($strand eq '+'){
			$start=substr($$maf{$s},$pos,3);
			$stop=substr($$maf{$s},$pos+$len-3,3);
		#}else{
		#	$stop=rc(substr($$maf{$s},$pos,3));
	#		$start=rc(substr($$maf{$s},$pos+$len-3,3));
	#	}
		$stop=uc($stop);
		$start=uc($start);

		if($s eq $$species[0]){
			$refstop=$stop;
			$$atgstopcons{$$species[0]}{'stop'}{'syn'}=0;
			$$atgstopcons{$$species[0]}{'stop'}{'cons'}=0;
			$$atgstopcons{$$species[0]}{'start'}{'cons'}=0;
		}


		$$atgstopcons{$$species[0]}{'start'}{$s}=-1;
		$$atgstopcons{$$species[0]}{'stop'}{$s}=-1;
		if($start eq 'ATG'){
			$$atgstopcons{$$species[0]}{'start'}{'cons'}++; ## count how often we have ATG in start
			$$atgstopcons{$$species[0]}{'start'}{$s}=1;
		}elsif(substr($start,0,1) eq '-'){ ## if first is a gap
			if(substr($start,1,1) eq '-'){ ## if second is a gap
				if(substr($start,2,1) =~ /[G-]/){ ## if third is a gap or G
					$$atgstopcons{$$species[0]}{'start'}{$s}=2;      # set if codon is --- or --G
				}

			}elsif(substr($start,1,1) eq 'T'){
				if(substr($start,2,1) eq 'G'){
					$$atgstopcons{$$species[0]}{'start'}{$s}=2;     # set if codon is -TG
				}
			}
		}
			
		
		
		if($stop =~ /TGA|TAG|TAA/){
			$$atgstopcons{$$species[0]}{'stop'}{'cons'}++; ## count for species 0 how often it is conserved
			$$atgstopcons{$$species[0]}{'stop'}{$s}=0;
			if($stop ne $refstop){
				$$atgstopcons{$$species[0]}{'stop'}{'syn'}++; ## how often stop has synmut
				$$atgstopcons{$$species[0]}{'stop'}{$s}=1; ## if species has mutation in stop then we set it to 
			}
		}elsif(substr($stop,2,1) eq '-'){ ## if first is a gap
			if(substr($stop,1,1) eq '-'){ ## if second is a gap
				if(substr($stop,0,1) =~ /[T-]/){ ## if third is a gap or T
					$$atgstopcons{$$species[0]}{'stop'}{$s}=2;      # set if codon is --- or T--
				}

			}elsif(substr($stop,1,1) =~ /[GA]/){ #TA- or TG-
				if(substr($stop,2,1) eq '-'){
					$$atgstopcons{$$species[0]}{'stop'}{$s}=2;     # set if codon is TG- or TA-      ##this is covering basically what we want
				}
			}
		}
	}
}


sub find_prev_met{
	my ($res,$orf_len,$max_size,$pos,$e,$tr,$ref_seq)=@_;
	my $tr2;
	if($orf_len < $max_size){
		for(my $i=$pos;$i > $e-$max_size;$i-=3){
			$tr2=uc(substr($$ref_seq,$i,3));
			last if($i < 0);
			last if($tr2 eq "TAA" or $tr2 eq "TGA" or $tr2 eq "TAG");
			if(uc(substr($$ref_seq,$i,3)) eq 'ATG'){
				$res=$i;
			}
		}
	}
	return $res;
}


sub syn_mut{
    my ($maf,$species,$pos,$len,$strand,$kaks_table,$pairwise_kaks_result) = @_;

#	my %delete;

	my ($syn,$nsyn,$usyn,$unsyn,$ka,$ks,$nonsense)= (0,0,0,0,0,0,0);
    my $refs=substr($$maf{$$species[0]},$pos,$len);
 
	#$refs=rc($refs) if($strand eq '-');

	my %uniq=();
	my %pairs=();
	my ($cod1,$cod2);
	my %nsc;
	$nsc{'TGA'}=1;
	$nsc{'TAA'}=1;
	$nsc{'TAG'}=1;	

	for (my $i=0; $i<$len;$i+=3){
		foreach my $s (@$species){
			$cod1=uc(substr($refs,$i,3));
			next if($s eq $$species[0]);
			
			$cod2=uc(substr($$maf{$s},$pos+$i,3));
		#	$cod1=rc($cod1) if($strand eq '-');
		#	$cod2=rc($cod2) if($strand eq '-');
#			print "\t$cod2\t" if($pos == 10408);
			
			## if codons are different
			if($cod2 ne $cod1){
				## if AA are the same then we have a synonymous mutation
				if(peptide($cod1) eq peptide($cod2)){
					$syn++;
					$uniq{$i}{'s'}{$cod2}++;
					$pairs{$s}{'s'}{$i}=undef; 
				}else{
					$nsyn++;
					$uniq{$i}{'n'}{$cod2}++;
					$pairs{$s}{'n'}{$i}=undef;
				}
			}
			$nonsense++ if($nsc{$cod2} and $i < $len-3);
		}
    }
	

	## here we do ka and ks calculation
	my $synsites=0;
	my $nonsynsites=0;
	my $cod;
	## first we get the number of synonymous and nonsynonymous sites
	# for(my $i=0;$i<$len-3;$i+=3){
	# 	$cod=substr($refs,$i,3);
	# 	#print "$cod  $$kaks{$cod}{1}+$$kaks{$cod}{2}+$$kaks{$cod}{3}\n";
	# 	$synsites+=$$kaks_table{$cod}{1}+$$kaks_table{$cod}{2}+$$kaks_table{$cod}{3};
	# }
	# $nonsynsites=$len-$synsites;

	# ## now we do the pairwise ka and ks calculation and sum up later
	# foreach my $s (@$species){
	# 	next if($s eq $$species[0]);
	# 	$$pairwise_kaks_result{$s}{'ka'} = scalar keys %{$pairs{$s}{'n'}};
	# 	$$pairwise_kaks_result{$s}{'ka'} /= $nonsynsites;
		
	# 	$$pairwise_kaks_result{$s}{'ks'} = scalar keys %{$pairs{$s}{'s'}};
	# 	$$pairwise_kaks_result{$s}{'ks'} /= $synsites;
		
	# 	if(!$$pairwise_kaks_result{$s}{'ka'} or !$$pairwise_kaks_result{$s}{'ks'}){
	# 		$$pairwise_kaks_result{$s}{'omega'} = -1;
	# 	}else{
	# 		$$pairwise_kaks_result{$s}{'omega'} = ($$pairwise_kaks_result{$s}{'ka'} / $$pairwise_kaks_result{$s}{'ks'});
	# 	}
	# }

	# ## we can now weigh it with the phylotree pairwise distances and then take the mean values, nice eh ?

    # ## load pairwise distances and weigh is

	# ## ok, this is somehow experimental but could maybe also work, needs some benchmarking
	# ## now sum up uniq codon mutation for syn mut and non syn mut
	 for my $k(keys %uniq){
	 	$usyn+= scalar keys %{$uniq{$k}{'s'}};
	 	$unsyn+= scalar keys %{$uniq{$k}{'n'}};
	 }
	# $ka= $unsyn/$nonsynsites;
	# $ks= $usyn/$synsites;

    return "$syn,$nsyn,$usyn,$unsyn,$ka,$ks,$nonsense";
}


sub getCodonMut{
	my ($maf,$species,$pos,$len,$strand) = @_;
	my (%pos,%mut,%spec,$seq2,$cod);
	my $seq=substr($$maf{$$species[0]},$pos,$len);
	if($strand eq '-'){
		#$seq=rc($seq);
	}

	for(my $i=0;$i<length($seq);$i+=3){
		for my $s(@$species){
			$seq2=substr($$maf{$s},$pos,$len);
			
			if($strand eq '-'){
		#		$seq2=rc($seq2);
			}
			
			$cod=lc(substr($seq2,$i,3));
			next if(lc(substr($seq,$i,3)) eq $cod);
			$mut{$i}{$cod}++;
			$spec{$s}++;
		}
	}

    ## get number of postitions for different species
	my $sum=0;
	for(my $i=0;$i<$len;$i+=3){
		if($mut{$i}){
			$sum+=scalar keys %{$mut{$i}}; ## here we score number of different syn mutations per site
		} 
	}

	my $spos=scalar keys %mut; ## how many pos are mutated
	my $sspec=scalar keys %spec;
	return "\t$sum,$spos,$sspec";
}

sub entropyAA{
    my ($maf,$species,$pos,$len,$strand) = @_;
    my $H=0;
    my %eh=();
    my $sn = scalar @$species;
    my $prob;
    my $h;
	my $cod;

   #

    for (my $i=$pos;$i< $pos+$len;$i+=3){
		foreach my $s (@$species){
			$cod='';
			if($strand eq '-'){
				$cod=$rcodons{uc(substr($$maf{$s},$i,3))};
			}else{
				$cod=$codons{uc(substr($$maf{$s},$i,3))};
			}
			$cod = '-' if(not $cod);
			$eh{$i}{$cod}++;
			
		}

		for my $nt(keys %{$eh{$i}}){

			$prob = $eh{$i}{$nt}/$sn;
			$h = ($prob*log2($prob));
			$H -= $h;
		}
    }

	$len/=3;

    return $H/$len;
}

sub entropy{
    my ($maf,$species,$pos,$len,$strand) = @_;
    my $H=0;
    my %eh=();
    my $sn = scalar @$species;
    my $prob;
    my $h;
	my $cod;

   #

    for (my $i=$pos;$i< $pos+$len;$i++){
		foreach my $s (@$species){
			$cod=uc(substr($$maf{$s},$i,1));
			$eh{$i}{$cod}++;			
		}

		for my $nt(keys %{$eh{$i}}){

			$prob = $eh{$i}{$nt}/$sn;
			$h = ($prob*log2($prob));
			$H -= $h;
		}
    }

    return $H/$len;
}


sub get_cons_species{
    my ($cons,$maf,$species,$aa,$pos,$orf_len,$strand) = @_;
    my $oseq='';
    
    
    ## check all species
    ##put refspec as first species
    push(@$cons,$$species[0]);
    for(my $i=1;$i < scalar @$species;$i++){
	## check if peptide is the same and put species to cons hash if inside

	## correct for wrong gaps in the end and beginning ...
	## not done yet
		
		next if(not $$maf{$$species[$i]});

		if($strand eq '-'){
			if(peptide(rc(substr($$maf{$$species[$i]},$pos,$orf_len))) eq $$aa){
				push(@$cons,$$species[$i]);
			}

		}else{
#			print "$$species[$i]\t",peptide(substr($$maf{$$species[$i]},$pos,$orf_len)),"\n";
			if(peptide(substr($$maf{$$species[$i]},$pos,$orf_len)) eq $$aa){
				push(@$cons,$$species[$i]);
			}
		}
    }
}

## here we check for dna species availability
sub get_cons_species_relaxed{
    my ($cons,$maf,$species,$pos,$orf_len,$strand,$frameshiftI) = @_;
    my $oseq='';
    my ($count,$checkseq,$kickout,$lp,$lg); 
    
    ## check all species
    ##put refspec as first species
    push(@$cons,$$species[0]);




	my $seq;
	my ($exchange,@border1,@border2,$len1,$len2,%fill,$need);
    my $maxg=4;
    for(my $i=1;$i < scalar @$species;$i++){
		next if(not $$maf{$$species[$i]});
		
		if($options{'G'}){
		## correct for wrong gaps in the end and beginning ...
		$seq = substr($$maf{$$species[$i]},$pos,$orf_len);
        if($seq =~ /^(-*)[acgtuACGTUnN]+(-*)$/){
			$len1 = length($1);
			$len2 = length($2);
 			
		#	print  " $$species[$i] $seq  $len1,$len2 \n" if($pos == 2192);

			if($len1 > 0 and $len1 < $orf_len/$maxg){
				## get index numbers of letters
				%fill=();
				# how many letters we need
				$need=$len1;
				
				## get border regions
				@border1 = split("",substr($$maf{$$species[$i]},$pos-$orf_len,$orf_len));
				for(my $i=$orf_len-1; $i >=0; $i--){
					if($border1[$i] ne '-'){
						$fill{$i}=undef;
						$need--;
						last if($need == 0);
					}
				}
				$exchange='';
				for my $k(sort {$a <=> $b} keys %fill){
					$exchange.=$border1[$k];   ## we add the key from the border1 letters
                    $border1[$k]='-'; ## no  we set the right border key
				}
                substr($$maf{$$species[$i]},$pos,$len1,$exchange); ## set new maf letters
                substr($$maf{$$species[$i]},$pos-$orf_len,$orf_len,join("",@border1)); ## set new maf letters in border region
			}
			if($len2 > 0 and $len2 < $orf_len /$maxg){
				## get index numbers of letters
				%fill=();
				# how many letters we need
				$need=$len2;
				## get border regions
				@border2 = split("",substr($$maf{$$species[$i]},$pos+$orf_len,$orf_len));
			
				for(my $i=0; $i < $orf_len; $i++){
					if($border2[$i] ne '-'){
						$fill{$i}=undef;
						$need--;
						last if($need == 0);
					}
				}
				$exchange='';
				for my $k(sort {$a <=> $b} keys %fill){
					$exchange.=$border2[$k];   ## we add the key from the border1 letters
                    $border2[$k]='-'; ## no  we set the right border key
				}
				

                substr($$maf{$$species[$i]},$pos+$orf_len-$len2,$len2,$exchange); ## set new maf letters

				substr($$maf{$$species[$i]},$pos+$orf_len,$orf_len,join("",@border2)); ## set new maf letters in border region
			}
			
		}
	    }
		$checkseq=substr($$maf{$$species[$i]},$pos,$orf_len);
#		print "$$species[$i] $checkseq\n" if($orf_len == 342);


		$count = ($checkseq =~ s/-/-/g); ## count gaps in alignments
		if(($count / $orf_len) <= $gapthreshold){
			## no frameshift gaps please
			
			#if($count % 3 == 0){
			## now we also check the gaps separately
			
				if(!$count){
					push(@$cons,$$species[$i]);
					
				}else{
#					for(my $i=0;$i<$orf_len;$i+=3){
					$kickout=0;
					$lp=();
					$lg=(); 
					while($checkseq =~ /(\-+)/g){ ## 
						$lp=$-[0];
						$lg=length($1);
						
						next if($lp == 0);              #only gap within sequence, ignore from beginning 
						next if($orf_len-$lg == $lp);   #only gap within sequence, ignore from end
						
					   

						if($lg % 3 == 0){
						}else{
							$kickout=1;
#							print "$lp $lg $count $pos $orf_len $$species[$i]\n" if($orf_len > 3600);
						}
					}
					if(not $kickout){
						push(@$cons,$$species[$i]);
					}else{
						push(@$frameshiftI,$$species[$i]); ## get species that have frameshift deletions
					}
				}
		}
    }
}




sub lines_to_array{
    my ($file,$arr) = @_;
    my $i=0;
    open IN,"<$file" or die "file $file could not be openend in lines_to_array\n";
    while(<IN>){
	last if(/^\s*$/);
	if(/^>*(\S+)/){
	    push(@$arr,$1);
	}
    }
}

sub initProt {
    my ($codons,$rcodons) = @_;
    my ($cod,$aa,$rcod);
    while(<DATA>){
        next if(/^\s*$/);
	chomp;
	($cod,$aa) = split();
	$$codons{uc($cod)} = $aa;
	$rcod = reverse($cod);
	$rcod =~ tr/acgtACGT/TGCATGCA/;
	$$rcodons{uc($rcod)} = $aa;
    }
}

sub peptide{
    my ($seq) = @_;
    return '-' if($seq =~ /-/);
    my $p='';
    my $cc;
    for(my $i=0;$i<length($seq);$i+=3){
        $cc =$codons{uc(substr($seq,$i,3))};
    $p.=$cc if($cc);
    }
    return $p;
}

sub read_fasta{
    my ($file,$hashabc,$finalpos,$speciesref,$skipnt)=@_;
    my $id;
	my $offset=0;
	my $add;
    open IN,"$file" or die "file not found\n";
#    open IN,"<$options{'m'}" or die "no maf file given with options m\n";
    while(<IN>){
	if(/>([a-zA-Z0-9]*)/){ ## read only id, no dots or other things allowed
	    $id = $1;
		return if($id ne $speciesref);
		## add offset characters
		if(/(\d+)-(\d+)/){
			$offset = $1;
			$$finalpos=$2;
		}
		$$skipnt = $offset;
		$$hashabc{$id}='-'x$offset;
	}elsif(/(\S+)/){
	    $$hashabc{$id}.=$1;
	}else{
	}
    }
    close IN;
}


sub log2{
    return log($_[0])/log(2);
}

## subroutine to insert each 3 digits a dot 
sub Nicenumber{
    my @numarr=split(/[.,]/,shift);    
    my $number = $numarr[0];

    my @n = split(//,reverse($number));
    my $res="";
    for(my $j=0; $j < length($number); $j++){
        if($j%3 eq 0 and $j ne 0){
            $res.=".$n[$j]";
        }else{
            $res.="$n[$j]";
        }
    }
    $res=reverse($res);
    $res.=",$numarr[1]" if($numarr[1]);
    return($res);
}

sub rc{
    my ($seq) = @_;
    $seq = reverse($seq);
    $seq =~ tr/acgtACGTuU/TGCATGCAaA/;
    return $seq;
}

## does a binary search for the closest element in an array given a query and an array reference, returns the index of the array with the closest value
## this returns already the distance not the index
## does a binary search for the closest elemnt in an array given a query and an array reference, returns the index of the array with the closest value
sub binarysearch2{
	my ($array,$query) = @_;
	my $i1=0;
	my $i2=scalar @$array;
	$i2--;
	my $mid='';


	while($i1<$i2){
		$mid=int(($i1+$i2)/2);
#		print @$array,"\n";
		if($$array[$mid] < $query){
			$i1=$mid+1;
		}else{
			$i2=$mid;
		}
	}
	$a=abs($$array[$i1]-$query);
	$b=abs($$array[$i1-1]-$query);

#	die $i1,"\n";

######## return index lowest value in array
	return $i1-1;# if($a < $b);
#	return ($i1-1) if($b <= $a);
#

}


sub binarysearch{
	my ($array,$query) = @_;
	my $i1=0;
	my $i2=scalar @$array;
	$i2--;
	my $mid='';

	while($i1<$i2){
		$mid=int(($i1+$i2)/2);
		if($$array[$mid] < $query){
			$i1=$mid+1;
		}else{
			$i2=$mid;
		}
	}
	$a=abs($$array[$i1]-$query);
	$b=abs($$array[$i1-1]-$query);


######## return index
	return $i1 if($a <= $b);
	return ($i1-1) if($b < $a);
######## return distance not index	
	return $a if($a < $b);
	return $b if($b < $a);

}

sub get_dtss_dtes{
	my ($tssa,$e,$tesa,$pos) = @_;
	my ($dtss,$dtes);

	$dtss=binarysearch($tssa,$e);
	
	while($$tssa[$dtss] < $e){
		$dtss++;
		
		if($dtss >= scalar @$tssa){
			$dtss="na";
			last;
		}
	}
	$dtss=abs($$tssa[$dtss]-$e) if($dtss ne 'na');
	
	$dtes=binarysearch($tesa,$pos);
	while($$tesa[$dtes] > $pos){
		$dtes--;
		$dtes='na' if($dtes < 0);
		last if($dtes eq 'na');
	}
	$dtes=abs($$tesa[$dtes]-$pos) if($dtes ne 'na');
	
	return "$dtss,$dtes";
}

sub checkopts{
	my (%options) = @_;
	if(defined $options{'m'} and not -f $options{'m'}){ die "File $options{'m'} supplied by option -m does not exist\n";}
	if(defined $options{'s'} and not -f $options{'s'}){ die "File $options{'s'} supplied by option -s does not exist\n";}
	if(defined $options{'b'} and not -f $options{'b'}){ die "File $options{'b'} supplied by option -b does not exist\n";}
	if(not defined $options{'c'} or $options{'c'} !~ /^\w+/){ die "No chromosome supplied by $options{'c'} or invalid chromosome name\n";}
	if(not defined $options{'t'} or $options{'t'} !~ /^\w+/){ die "No project name supplied by $options{'t'} or invalid project name\n";}
}






__DATA__
aaa                 K
aac                 N
aag                 K
aat                 N
aca                 T
acc                 T
acg                 T
act                 T
aga                 R
agc                 S
agg                 R
agt                 S
ata                 I
atc                 I
atg                 M
att                 I
caa                 Q
cac                 H
cag                 Q
cat                 H
cca                 P
ccc                 P
ccg                 P
cct                 P
cga                 R
cgc                 R
cgg                 R
cgt                 R
cta                 L
ctc                 L
ctg                 L
ctt                 L
gaa                 E
gac                 D
gag                 E
gat                 D
gca                 A
gcc                 A
gcg                 A
gct                 A
gga                 G
ggc                 G
ggg                 G
ggt                 G
gta                 V
gtc                 V
gtg                 V
gtt                 V
taa                 *
tac                 Y
tag                 *
tat                 Y
tca                 S
tcc                 S
tcg                 S
tct                 S
tga                 *
tgc                 C
tgg                 W
tgt                 C
tta                 L
ttc                 F
ttg                 L
ttt                 F

