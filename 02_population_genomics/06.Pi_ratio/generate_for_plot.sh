join -j 2 combineXJtumblerParlorRoller60_pi_w20k_s5k.windowed.scaffold60.pi others55_pi_w20k_s5k.windowed.scaffold60.pi | awk '{print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $8 "\t" $9}'> combineXJtumblerParlorRoller60_vs_others55_w20k_s5k.windowed.pi
sed -i '1d' combineXJtumblerParlorRoller60_vs_others55_w20k_s5k.windowed.pi
awk '{print $1 ":" $2 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $7/$5}' combineXJtumblerParlorRoller60_vs_others55_w20k_s5k.windowed.pi > combineXJtumblerParlorRoller60_vs_others55_w20k_s5k.windowed.piRatio.txt
sed -i "1iSNP\tCHR\tBP\tBP_END\tvariants_combineXJtumblerParlorRoller\tpi_combineXJtumblerParlorRoller\tvariants_other\tpi_other\tP" combineXJtumblerParlorRoller60_vs_others55_w20k_s5k.windowed.piRatio.txt
awk '$5>=50 && $7>=50 {print $0}' combineXJtumblerParlorRoller60_vs_others55_w20k_s5k.windowed.piRatio.txt > combineXJtumblerParlorRoller60_vs_others55_w20k_s5k.windowed.piRatio.filterN10.txt


join -j 2 XJtumbler38_pi_w20k_s5k.windowed.scaffold60.pi others55_pi_w20k_s5k.windowed.scaffold60.pi | awk '{print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $8 "\t" $9}'> XJtumbler38_vs_others55_w20k_s5k.windowed.pi
sed -i '1d' XJtumbler38_vs_others55_w20k_s5k.windowed.pi
awk '{print $1 ":" $2 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $7/$5}' XJtumbler38_vs_others55_w20k_s5k.windowed.pi > XJtumbler38_vs_others55_w20k_s5k.windowed.piRatio.txt
sed -i "1iSNP\tCHR\tBP\tBP_END\tvariants_XJtumbler\tpi_XJtumbler\tvariants_other\tpi_other\tP" XJtumbler38_vs_others55_w20k_s5k.windowed.piRatio.txt
awk '$5>=50 && $7>=50 {print $0}' XJtumbler38_vs_others55_w20k_s5k.windowed.piRatio.txt > XJtumbler38_vs_others55_w20k_s5k.windowed.piRatio.filterN10.txt


join -j 2 parlorRoller22_pi_w20k_s5k.windowed.scaffold60.pi others55_pi_w20k_s5k.windowed.scaffold60.pi | awk '{print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $8 "\t" $9}'> parlorRoller22_vs_others55_w20k_s5k.windowed.pi
sed -i '1d' parlorRoller22_vs_others55_w20k_s5k.windowed.pi
awk '{print $1 ":" $2 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $7/$5}' parlorRoller22_vs_others55_w20k_s5k.windowed.pi > parlorRoller22_vs_others55_w20k_s5k.windowed.piRatio.txt
sed -i "1iSNP\tCHR\tBP\tBP_END\tvariants_parlorRoller\tpi_parlorRoller\tvariants_other\tpi_other\tP" parlorRoller22_vs_others55_w20k_s5k.windowed.piRatio.txt
awk '$5>=50 && $7>=50 {print $0}' parlorRoller22_vs_others55_w20k_s5k.windowed.piRatio.txt > parlorRoller22_vs_others55_w20k_s5k.windowed.piRatio.filterN10.txt
