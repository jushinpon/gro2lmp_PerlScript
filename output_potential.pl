use warnings;
use strict;

my %pairsty;
   $pairsty{1} = "lj/gromacs/coul/gromacs 9.0 11.0  0.000001 11.0"; 
my %bondsty;
   $bondsty{1} = "harmonic";                      
my %angsty;
   $angsty{1} = "harmonic";                      
   $angsty{2} = "cosine/squared";                      
my %dihsty;
   $dihsty{1} = "fourier";                      
   $dihsty{9} = "fourier";                      
   $dihsty{2} = "quadratic";                      
my $dielectric = 15.0;

sub output_potential{

my ($atomobj_a,$atom_lookup_h,$atomtype_a,$atomtypeID_h,$coorobj_a,$grobox_a,
$bondobj_a,$bond_lookup_h,$bondtype2name_a,$bondname2type_h,$bondID2name_h,
$cbond_a,$ebond_a,
$angleobj_a,$angle_lookup_h,$angletype2name_a,$anglename2type_h,$angleID2name_h,
$Bdihobj_a,$Bdih_lookup_h,$Bdihtype2name_a,$Bdihname2type_h,$BdihID2name_h,
$Sdihobj_a,$Sdih_lookup_h,$Sdihtype2name_a,$Sdihname2type_h,$SdihID2name_h,
$nonbond_lookup_h) = @_;
`rm -rf ./output_potential/`;
`mkdir output_potential`;
############make potential file
open my $df, "> ./output_potential/potential.in"; 
#pair_style
print $df "pair_style $pairsty{1}\n";
for (sort keys %{$nonbond_lookup_h}){	
	print $df "pair_coeff $_ @{$nonbond_lookup_h->{$_}}[1..3]\n";
}

print $df "\n";
print $df "bond_style $bondsty{1}\n";
for (1..@{$bondtype2name_a}){	
	my $name = $bondtype2name_a->[$_ - 1];
	my $ID1 = 
	        # @{$bond_lookup_h{$name}}[1,2];
	print $df "bond_coeff $_ @{$bond_lookup_h->{$name}}[1,2] #$name\n";
}
#fake bond
my $type = scalar @{$bondtype2name_a} + 1;	
#for (@{$ebond_a}){	
	print $df "bond_coeff $type 0.0 2.5 #for fake bond only\n";
#}

print $df "\n";
print $df "angle_style hybrid $angsty{1} $angsty{2}\n";
for (1..@{$angletype2name_a}){	
	my $name = $angletype2name_a->[$_ - 1];
	my $stlyeID = ${$angle_lookup_h->{$name}}[0];	
	print $df "angle_coeff $_ $angsty{$stlyeID} @{$angle_lookup_h->{$name}}[1,2] #$name\n";
}

print $df "\n";
print $df "dihedral_style hybrid $dihsty{1} $dihsty{2}\n";
for (1..@{$Bdihtype2name_a}){	
	my $name = $Bdihtype2name_a->[$_ - 1];
	my $stlyeID = ${$Bdih_lookup_h->{$name}}[0][0];# hash-> two D array
	if($dihsty{$stlyeID} eq "fourier"){	
		my $mult = scalar @{$Bdih_lookup_h->{"$name"}};
		print $df "dihedral_coeff $_ $dihsty{$stlyeID} $mult ";
		for my $ml (0..$mult-1){			
			print $df " @{$Bdih_lookup_h->{$name}->[$ml]}[1..3]";
		}
		print $df " #$name\n";
	}
}
#side chain dihedral
for (@{$Bdihtype2name_a}+1..@{$Bdihtype2name_a}+@{$Sdihtype2name_a}){	
	my $name = $Sdihtype2name_a->[$_ - @{$Bdihtype2name_a} - 1];
	my $stlyeID = ${$Sdih_lookup_h->{$name}}[0][0];# hash-> two D array
	if($dihsty{$stlyeID} eq "quadratic"){	
		print $df "dihedral_coeff $_ $dihsty{$stlyeID} @{$Sdih_lookup_h->{$name}->[0]}[1,2] #$name\n";
		
	}
}
print $df "\n";
print $df "dielectric $dielectric\n";
close($df);

##make restrain files
open my $rs1, "> ./output_potential/restrain1.in";
open my $rs2, "> ./output_potential/restrain2.in";
open my $rl, "> ./output_potential/restrain_leng.in";

print $rs1 "variable firstID index "; 
print $rs2 "variable secondID index "; 
print $rl "variable bondleng index "; 
for (@{$cbond_a}){
	print $rs1 "&\n\"$_->[0]\" "; 
	print $rs2 "&\n\"$_->[1]\" "; 
	print $rl "&\n\"$_->[2]\" ";
}
close($rs1);
close($rs2);
close($rl);

}# end sub
1;
