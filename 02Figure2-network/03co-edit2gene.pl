use 5.24.1;
my %hash;
open(IN0,"D:/HFU/RNA_editing_combine_regulation/data/RADAR_edit_type.txt")||die;
while(<IN0>){
chomp;
my @a = split/\t/,$_;
$hash{$a[0]} = $a[5];
 }

open(IN,"D:/HFU/RNA_editing_combine_regulation/res/RADAR_CO_network.txt")||die;
open(OUT,">D:/HFU/RNA_editing_combine_regulation/res/RADAR_CO_network_gene.txt")||die;
while(<IN>){
	chomp;
	my @b = split/\t/,$_;
my $tmp =$b[0]."\t".$b[1]."\t". $hash{$b[0]}."\t".$hash{$b[1]};
print OUT "$tmp\n";
	}
