############ Perl Script to generate rnk files for GSEPreranked

use POSIX;
my $markerfile = $ARGV[0]; chomp $markerfile; ###### Seurat Marker File
my $prefix=$ARGV[1]; chomp $prefix; ######### For example: Cluster 
my $outfolder=$ARGV[2]; chomp $outfolder; ###### Output directory
my %hash;
my $counter=0;
open(F,"$markerfile");
while(my $data = <F>)
{
 $data =~ s/^\s+|\s+$//g;
 chomp $data;
 $counter++;
 if($counter == 1)
 {}
 else
 {
  my @arr=split(/\t/,$data);
  my $clusterId=$prefix.$arr[6];
  my $logFC=$arr[2];
  my $gene=$arr[7];
  my $value=$gene."\t".$logFC;
  $hash{$clusterId}{$value}=0;
 }
}
close F;
if(-d 'Cluster_Rnk_Files')
{}
else
{
 mkdir($outfolder);
}
foreach my $m(sort keys %hash)
{
 my @logfc; my @genes;
 open(OUT,">$outfolder/$m.rnk");
 print OUT "\#\#\n";
 foreach my $n(keys %{$hash{$m}})
 {
  my @arr=split(/\t/,$n);
  push(@genes,$arr[0]);
  push(@logfc,$arr[1]);
 }
 my @logfcsort=sort { $b <=> $a } @logfc;
 my @indexes;
 for(my $i=0;$i<@logfcsort;$i++)
 {
  for(my $j=0;$j<@logfc;$j++)
  {
   if($logfcsort[$i] == $logfc[$j])
   {
    push(@indexes,$j);
   }
  }
 }
 my @genes_use=@genes[@indexes];
 for(my $i=0;$i<@genes_use;$i++)
 {
  print OUT "$genes_use[$i]\t$logfcsort[$i]\n";
 }
 close OUT;
}