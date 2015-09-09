#!/usr/bin/perl


# open bed file: normalized_filtered_hits.bed
while(<STDIN>)
{
  @array = split;
  $start = $array[1];
  $end = $array[2];
  $length = $end - $start +1;

## start to shuffle
  $current = $length;
   #print "in loop";
   $randChrom = int(rand(6));

   #print "random" .  $randChrom;


   if($randChrom eq 0)
   {
     $chrom = "I";
     #print "CHRPM" . $chrom . "\t";
     $randStart = int(rand(15323679 - $current));
    # print "START" . $randStart . "\t";
   }
   if($randChrom eq 1)
   {
     $chrom = "II";
     $randStart = int(rand(15534032 - $current));
   }
   if($randChrom eq 2)
   {
     $chrom = "III";
     $randStart = int(rand(14013465 - $current));
   }
   if($randChrom eq 3)
   {
     $chrom = "IV";
     $randStart = int(rand(17785400- $current));
   }
   if($randChrom eq 4)
   {
     $chrom = "V";
     $randStart = int(rand(21272929 - $current));
   }
   if($randChrom eq 5) 
   {
     $chrom = "X";
     $randStart = int(rand(18014219 - $current));
   }


$end = $randStart + $current;

print $chrom . "\t" . $randStart . "\t" . $end;
for ($i=3;$i<8;$i++) {	print "\t".$array[$i];		}
print "\n";
}



