#!/usr/bin/perl

open result, ">Target_combineXJtumblerParlorRoller60_fele20otherTumbler6.lsbl.W10s5";
            
open vcf, "<Target_combine60_fele20otherTumbler6_SNP.weir.replace0.LSBL";
while (<vcf>) {
        chomp $_;
        @eachline=split (/\t/, $_);
        push @{$chr_position{$eachline[0]}}, $eachline[1];
		$fst{$eachline[0]}{$eachline[1]} = $eachline[5];
}
close vcf;


open fai,"<fai.list";



#for ($ii=1;$ii<=38;$ii++) {
while(<fai>){
    chomp;
    $ii=$_; 
       @sort_position= @{$chr_position{$ii}};

        @sort_position= sort {$a <=> $b} @sort_position;
                print "$sort_position[0]\n";
                print "$sort_position[-1]\n";

        for ($j=$sort_position[0];$j<=$sort_position[-1];$j=$j+5000) {
                $start=$j;
                $end= $j+10000;
				
				$fst_window=0;
				$snp_number= 0;

                for ($x=$start;$x<= $end;$x++) {
                        $t= exists $fst{$ii}{$x};
                        if ($t==1) {
							$fst_window= $fst_window + $fst{$ii}{$x};
							$snp_number++;
						}
				}
				if ($snp_number>3) {
				$mean_fst_window = $fst_window/$snp_number;

				print result "$ii\t$start\t$end\t$fst_window\t$snp_number\t$mean_fst_window\n";
				}
		}
}



