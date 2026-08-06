open f1,"<fst_XJtumbler38_vs_fele20otherTumbler6_SNP.weir.fst.edit2";
open f2,"<fst_XJtumbler38_vs_other147_SNP.weir.fst.edit2";
open f3,"<fst_other147_vs_fele20otherTumbler6_SNP.weir.fst.edit2";
open f4,">Target_xinjiang38_fele20otherTumbler6_SNP.weir.replace0.LSBL";

while (<f1>){
chomp;
@aa=split(/\s+/,$_);
$hash1{$aa[0]}{$aa[1]}=$aa[2];
}

while (<f2>){
chomp;
@ab=split(/\s+/,$_);
$hash2{$ab[0]}{$ab[1]}=$ab[2];
}

print f4 "#chr\tpos\txinjiang38_vs_fele20otherTumbler6\txinjiang38_vs_other147\tothers147_vs_fele20otherTumbler6\tTarget_xinjiang38_LSBL\n";

while (<f3>){
chomp;
@ac=split(/\s+/,$_);
if ((exists $hash1{$ac[0]}{$ac[1]})and (exists $hash2{$ac[0]}{$ac[1]})){
   $lsbl=($hash1{$ac[0]}{$ac[1]}+ $hash2{$ac[0]}{$ac[1]}-$ac[2])/2;
#   if ($lsbl<0){$lsbl=0};
   print f4"$ac[0]\t$ac[1]\t$hash1{$ac[0]}{$ac[1]}\t$hash2{$ac[0]}{$ac[1]}\t$ac[2]\t$lsbl\n";
  }
#$hash{$aa[0]}{$aa[1]}=$aa[2];
}
