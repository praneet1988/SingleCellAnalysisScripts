########### Perl Script to automate GSEAPreranked  analysis

my $dir=$ARGV[0]; chomp $dir; ########Directory with Ranked genes to test for enrichment
my $gmt=$ARGV[1]; chomp $gmt; ####GMT Gene Set file
my $outdir=$ARGV[2]; chomp $outdir;  ##### Output Drectory
opendir(F,"$dir");
foreach my $f(readdir(F))
{
 if($f =~ /.rnk/)
 {
  my @arrf=split(/\.rnk/,$f);
  my $file=$dir."\/".$f;
  system('/Applications/GSEA_4.0.3/./gsea-cli.sh', 'GSEAPreranked', '-gmx', $gmt, '-rnk', $file, '-set_min', '2', '-set_max', '1000', '-out', $outdir, '-rpt_label', $arrf[0]);
 }
}
closedir F;