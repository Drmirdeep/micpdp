package My_micPDP; 

use strict; 
use Exporter; 
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK); 
$VERSION = 1.00; 
@ISA = qw(Exporter); 
@EXPORT = (); 
@EXPORT_OK = qw(get_codon_pos_index dnds rc cmpl round log2 kaks peptide); 


sub rc{
    my ($rc) = @_;
	$rc=reverse($rc);
    cmpl(\$rc);
    return $rc;
}

sub cmpl{
    my ($seq)=(@_);
    $$seq =~ tr/acgtuACGTU/tgcaaTGCAA/;
}



sub log2{
    return log($_[0])/log(2);
}

sub round{
	return int(1000*$_[0])/1000;
}

sub get_codon_pos_index{
	## reads in bed file and notes first codon position of cds for each codon in hash cp!
	my ($f,$cp,$chr)= @_;
	open IN,$f or die "File with exons not given\n";
	
	my (@line,$rel,$curbp,%ph,@bl,@bs);

	while(<IN>){
		@line=split();
		next if($line[0] ne $chr);
		$ph{'strand'} = $line[5];
		$ph{'name'} = $line[3];
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
			die "no 12 coloumn bed format given\n";
		}
		
		my $llen=0;
		for(my $i= 0;$i<scalar @bs;$i++){
			$llen+=$bl[$i];
		}
		my $pa=-1;
		$llen--;
	
		my $codoncount=-1;

		for(my $i= 0;$i<scalar @bs;$i++){
			$curbp=-1;
			$ph{'len'}+=$bl[$i];
			for(my $s=$bs[$i]; $s < $bs[$i]+$bl[$i];$s++){
				$curbp++;
				$pa++;
				if($ph{'strand'} eq '+'){
					$rel=$pa;
				}else{
					$rel=$llen-$pa;
				}
				
				$ph{$rel}{'abs'}=$line[1]+$s;  ## absolute genome position
								
			

				if($ph{'strand'} eq '-'){
					if($ph{$rel}{'abs'} == $line[6]){
						$codoncount=1;
#						print "$rel $ph{$rel}{'abs'} $line[6]\n";
#						die;
					}
					
					if($ph{$rel}{'abs'} == $line[7]){ 
						$codoncount=-1;
#						print "$rel $ph{$rel}{'abs'} $line[7]\n";
#						die;
					}
					
					if($codoncount!=-1){
						$codoncount++ if( $ph{$rel}{'abs'} != $line[6]);
						#print "$rel\t$ph{$rel}{'abs'}\n";
						if($codoncount %3 == 0){
							$$cp{$ph{$rel}{'abs'}}.=$line[5];
#							die "$rel\t$s\t$ph{$rel}{'abs'} $codoncount\n";
						}
					}
				}else{
					if($ph{$rel}{'abs'} == $line[6]){ 
						$codoncount++;
					}
					if($ph{$rel}{'abs'} == $line[7]-1){
						$codoncount=-1;
					}
					if($codoncount>-1){
						$codoncount++ if( $ph{$rel}{'abs'} != $line[6]);
						$$cp{$ph{$rel}{'abs'}}.=$line[5] if($codoncount %3 == 0);
#						print $ph{$rel}{'abs'},"\n" if($codoncount %3 == 0);
#						die if( $ph{$rel}{'abs'} >  17051730);
					}	
				}
		#		print "$rel\t$s\t$ph{$rel}{'abs'} $codoncount\n";
		#		die if($ph{$rel}{'abs'} =~ /0645$/); 

			}
		}
#		die $$cp{17051719};
#		exit;
	}
return 1;	
}
# for my $k(sort {$a <=> $b} keys %cp){
# 	print "$k\n";
# }


sub dnds{
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

	my $ms=0;
	my %is=();

    foreach my $s(@$species){
        next if($s eq $$species[0]); ## without refspec
		$sum1= $$kaks{$name}{$s}{$type}{'ka'};
		$sum2= $$kaks{$name}{$s}{$type}{'ks'};

		if($sum1 ==1 and $sum2 == 1){            #if(ka and ks ==1 then we only have gaps in sequence and take it out from calculus
			$ms++;
			$is{$s}=1;
			next;
		}

        $e{$s} = $sum1/$sum2;
		#print "$s\t$sum1\t$sum2\t$e{$s}\n";
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
	if((scalar @$species) -$ms ==1){ ## why does this happen ????
		$$dnds{$name}{$type}{'median'} =-1;
		$$dnds{$name}{$type}{'mean'} =-1;
		$$dnds{$name}{$type}{'var'} = -1;
        $$dnds{$name}{$type}{'sd'} =-1;
		#print STDERR  @$species,"\t$name\n";
		return;
	}

	

    $$dnds{$name}{$type}{'mean'} = $sumt/(scalar @$species -1 -$ms);
    my $var2=0;
    foreach my $s(@$species){
		next if($s eq $$species[0]); ## without refspec
		next if($is{$s});
		$var2+= ($e{$s}-$$dnds{$name}{$type}{'mean'})**2
    }

    if((scalar @$species) -$ms ==2){
        $$dnds{$name}{$type}{'var'} = 0;
        $$dnds{$name}{$type}{'sd'} =0;
        return;
    }
    ## we need at least 3 species otherwise var cannot be calculated, one spec is the reference ........

    $var2/=(scalar @$species -2 -$ms); ## corresponds to n-1 since we used already the mean and thus loose one degree of freedom

    #$dnds{$name}{$type}{'var'} = ($sumt2/(scalar @$species)) - ($sumt/(scalar @$species))**2 ;
    $$dnds{$name}{$type}{'var'} = $var2;

    ## somehow we got ausloeschung since numbers are small for machine formula, other
    #print "$dnds{$name}{$type}{'var'}\t$dnds{$name}{$type}{'var2'}\n";
    #if($dnds{$name}{$type}{'var'} <0){
    #   die "$name\t$dnds{$name}{$type}{'var'}\t$var2\n";
    #}
    $$dnds{$name}{$type}{'sd'} = sqrt( $var2 );# if($dnds{$name}{$type}{'var'} > 0);
}


sub hamming {
    ## strings are compared in binary format with xor, then we simply count the remaining zeros
	my ($cod1,$cod2,$optionsY)= @_;
	return 1 if($optionsY);
    return ($_[0] ^ $_[1]) =~ tr/\0//c;
}


sub kaks{
    my ($maf,$species,$len,$kaks_table,$pairwise_kaks_result,$name,$type,$ph,$kick,$rels,$codons,$optionsB,$optionsY,$optionsy) = @_;

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
	my $codk;

    foreach my $s (@$species){
        $pairs{$s}{$type}{'s'}{'h'}=0;
        $pairs{$s}{$type}{'n'}{'h'}=0;
		#print "$refs,$s,@$species\n";
		$subtractlen=0; ## how many codons to discard ?
		$scodignore{$s}{'l'}=0;
#		print $s;
		next if($s eq $$species[0]);
        for (my $i=0+$len%3; $i<$len;$i+=3){
            $cod1=uc(substr($refs,$i,3));

#            next if($s eq $$species[0]);   #### skip refspec
            $cod2=uc(substr($$maf{$s},$i,3));

		
			## here we ignore things which overlap with codons from coding genes
			$codk=0;
			if($type eq 'orf' and $optionsB){
				$kicks=$rels+$i;
				while($kicks > $$ph{'circl'}){
					$kicks -= $$ph{'circl'};
				}
				
				
				if(exists $$kick{$$ph{$kicks}{'abs'}}){
					$scodignore{$s}{'ignorepos'}{$i}=1; ## put codon position to ignore in hash
					$scodignore{$s}{'l'}+=3;    ## remove 3 from length in the end
					# if($rels == 400 and $s eq 'mm9'){
					# 	print "$i ";
					# }
					$codk=1;
				}
			}

			## we ignore codons with gaps
			if($optionsY){
				if($cod2 =~ /\-/ and $codk == 0){
					$scodignore{$s}{'ignorepos'}{$i}=1; ## put codon position to ignore in hash
					$scodignore{$s}{'l'}+=3;    ## remove 3 from length in the end
				}elsif($cod2 =~ /-/){
					next;
				}
			}




            ## if codons are different
            if($cod2 ne $cod1){
                ## if AA are the same then we have a synonymous mutation
                $hd=hamming($cod1,$cod2,$optionsy);

                if(peptide($cod1,$codons) eq peptide($cod2,$codons)){
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
#		if($optionsC){
			## init sites
			$scodignore{$s}{'synsites'}=0;
			$scodignore{$s}{'nonsynsites'}=0;

			$scodignore{$s}{'len'}=$len-$scodignore{$s}{'l'};    ## remove 3 from length in the end
		#	print "$s $scodignore{$s}{'len'}   $refs\n";
#		}
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
		#if($optionsC){
			foreach my $s(@$species){
				next if($s eq $$species[0]);
				if(not $scodignore{$s}{'ignorepos'}{$i}){
					$scodignore{$s}{'synsites'}+=$$kaks_table{$cod}{1}+$$kaks_table{$cod}{2}+$$kaks_table{$cod}{3};
				}
			}
		#}


    }
	
    $nonsynsites=$len-$synsites;
	#if($optionsC){
		## now species specific nonsynsites
		foreach my $s(@$species){
			next if($s eq $$species[0]);
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
	#}

#	print "\n$refs\n";
#	print "$nonsynsites=$len-$synsites\n";
#	foreach my $s(@$species){
#		next if($s eq $species[0]);
#		print "$s $scodignore{$s}{'len'} \n";
#		print "$s ",keys %{$scodignore{$s}{'ignorepos'}},"\n";
#		print "$s $scodignore{$s}{'synsites'}  $scodignore{$s}{'nonsynsites'} \n";
#		print "$s $scodignore{$s}{'nonsynsites'}=$scodignore{$s}{'len'}-$scodignore{$s}{'synsites'}\n";	
#	}

	my $addsyn=$synsites/$len;
	my $addnsyn=$nonsynsites/$len;
	
	$synsites+=$addsyn;
	$nonsynsites+=$addnsyn;


    ## here comes now the tricky part
    ## this is Ka Ks with hamming distance
    foreach my $s (@$species){
        next if($s eq $$species[0]);
        $$pairwise_kaks_result{$name}{$s}{$type}{'ka'} = $pairs{$s}{$type}{'n'}{'h'};

		#if($optionsC){
			$$pairwise_kaks_result{$name}{$s}{$type}{'ka'}+=$scodignore{$s}{'addnsyn'};
			if($scodignore{$s}{'nonsynsites'}){
				$$pairwise_kaks_result{$name}{$s}{$type}{'ka'} /= $scodignore{$s}{'nonsynsites'};
			}else{
				$$pairwise_kaks_result{$name}{$s}{$type}{'ka'} = 1; ## we have 0 mutations
			}
		# }else{
		# 	$$pairwise_kaks_result{$name}{$s}{$type}{'ka'}+=$addnsyn;
		# 	if($nonsynsites){
		# 		$$pairwise_kaks_result{$name}{$s}{$type}{'ka'} /= $nonsynsites;
		# 	}else{  ## this should never happen since we added a pseudocount already above
		# 		$$pairwise_kaks_result{$name}{$s}{$type}{'ka'} =1;  ## we have 0 mutations
		# 	}
		# }
		
        $$pairwise_kaks_result{$name}{$s}{$type}{'ks'} = $pairs{$s}{$type}{'s'}{'h'};
		

		## $$pairwise_kaks_result{$name}{$s}{$type}{'ks'}++ if($options{'C'});
		## we just add one to the counts
		#if($options{'C'}){
		#if($optionsC){
			$$pairwise_kaks_result{$name}{$s}{$type}{'ks'}+=$scodignore{$s}{'addsyn'};
			if($scodignore{$s}{'synsites'}){ ## can be 0 if length is 0
				$$pairwise_kaks_result{$name}{$s}{$type}{'ks'} /= $scodignore{$s}{'synsites'};
			}else{
				$$pairwise_kaks_result{$name}{$s}{$type}{'ks'} =1;
			}
		# }else{
		# 	$$pairwise_kaks_result{$name}{$s}{$type}{'ks'}+=$addsyn;
		# 	if($synsites){
		# 		$$pairwise_kaks_result{$name}{$s}{$type}{'ks'} /= $synsites;
		# 	}else{ ## this should never happen since we added a pseudocount already above
		# 		$$pairwise_kaks_result{$name}{$s}{$type}{'ks'} =1;
		# 	}
		# }
	
#		print "$s\t$$pairwise_kaks_result{$name}{$s}{$type}{'ka'}\t$$pairwise_kaks_result{$name}{$s}{$type}{'ks'}\t",$$pairwise_kaks_result{$name}{$s}{$type}{'ka'}/$$pairwise_kaks_result{$name}{$s}{$type}{'ks'},"\n";
	}
}



sub peptide{
    my ($seq,$codons) = @_;
    my $p='';
    for(my $i=0;$i<length($seq);$i+=3){
	$p.=$$codons{uc(substr($seq,$i,3))};
    }
    return $p;
}
