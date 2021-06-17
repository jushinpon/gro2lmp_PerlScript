use warnings;
use strict;

sub output_datafile{

my ($atomobj_a,$atom_lookup_h,$atomtype_a,$atomtypeID_h,$coorobj_a,$grobox_a,
$bondobj_a,$bond_lookup_h,$bondtype2name_a,$bondname2type_h,$bondID2name_h,
$cbond_a,$ebond_a,
$angleobj_a,$angle_lookup_h,$angletype2name_a,$anglename2type_h,$angleID2name_h,
$Bdihobj_a,$Bdih_lookup_h,$Bdihtype2name_a,$Bdihname2type_h,$BdihID2name_h,
$Sdihobj_a,$Sdih_lookup_h,$Sdihtype2name_a,$Sdihname2type_h,$SdihID2name_h,$link_arraya,
$nonbond_lookup_h) = @_;
#`rm -rf ./output_datafile/`;
#`mkdir output_datafile`;
############make data file
open my $df, "> output.data"; 
print $df "#lammps data file from Perl script\n";
print $df "\n";
print $df scalar(@{$atomobj_a}) . " atoms\n";
print $df scalar(@{$bondobj_a}) + scalar(@{$ebond_a}) + scalar(@{$cbond_a}) . " bonds\n";# real + fake bonds 
print $df scalar(@{$angleobj_a}) . " angles\n";
print $df scalar(@{$Bdihobj_a}) . " dihedrals\n";
print $df scalar(@{$Sdihobj_a}) . " impropers\n";


print $df "\n";
print $df scalar(@{$atomtype_a}) . " atom types\n";
print $df scalar(@{$bondtype2name_a}) + scalar(@{$cbond_a}) + 1 . " bond types\n";# add one for fake bond
print $df scalar(@{$angletype2name_a}) . " angle types\n";
print $df scalar(@{$Bdihtype2name_a}) . " dihedral types\n";
print $df scalar(@{$Sdihtype2name_a}) . " improper types\n";
print $df "\n";
print $df "0.0 ". ${$grobox_a}[0][0]." xlo xhi\n";
print $df "0.0 ". ${$grobox_a}[0][1]." ylo yhi\n";
print $df "0.0 ". ${$grobox_a}[0][2]." zlo zhi\n";
print $df "\n";
print $df "Masses\n";
print $df "\n";
for (1..@{$atomtype_a}){
	my $atomtype = ${$atomtype_a}[$_ - 1];
	print $df "$_ " . ${$atom_lookup_h}{$atomtype} . " #$atomtype\n";
}
print $df "\n";
print $df "Atoms\n";
print $df "\n";
for (1..@{$atomobj_a}){
    my $atomname = ${$atomobj_a}[$_ - 1][1];
    my $atomtype = ${$atomtypeID_h}{$atomname};# name to type ID
    my $mID = ${$atomobj_a}[$_ - 1][2];
    my $charge = ${$atomobj_a}[$_ - 1][3];
    my @coorxyz = @{$coorobj_a->[$_ - 1]}[0..2];
	print $df "$_ $mID $atomtype $charge @coorxyz" . " #$atomname\n";
}
print $df "\n";
print $df "Bonds\n";
print $df "\n";

for (1..@{$bondobj_a}){

    my $name = ${$bondID2name_h}{$_};
    my $type = ${$bondname2type_h}{$name};# name to type ID
    my @bondID = @{$bondobj_a->[$_ - 1]}[0..1];
	print $df "$_ $type @bondID" . " #$name\n";
}
#constrained bond
my $extra_bondtype = scalar(@{$bondtype2name_a}); # this number so far from the read bond above
for (0..$#{$cbond_a}){
    $extra_bondtype++;
    my $bondID = scalar(@{$bondobj_a}) + $_ + 1;
    my @pairID = @{$cbond_a->[$_]}[0..1];
	print $df "$bondID $extra_bondtype @pairID" . " #constrained bond\n";
}

## fake bond (the last bond type)
$extra_bondtype++;

for (0..$#{$ebond_a}){
    my $bondID = scalar(@{$bondobj_a}) + scalar(@{$cbond_a})+ $_ + 1;
    my @pairID = @{$ebond_a->[$_]}[0..1];
	print $df "$bondID $extra_bondtype @pairID" . " #fake bond\n";
}
print $df "\n";
print $df "Angles\n";
print $df "\n";
for (1..@{$angleobj_a}){

    my $name = ${$angleID2name_h}{$_};
    my $type = ${$anglename2type_h}{$name};# name to type ID
    my @angleID = @{$angleobj_a->[$_ - 1]}[0..2];
	print $df "$_ $type @angleID" . " #$name\n";
}
print $df "\n";
print $df "Dihedrals\n";
print $df "\n";
for (1..@{$Bdihobj_a}){
    my $name = ${$BdihID2name_h}{$_};
    my $type = ${$Bdihname2type_h}{$name};# name to type ID
    my @BdihID = @{$Bdihobj_a->[$_ - 1]}[0..3];
	print $df "$_ $type @BdihID" . " #$name\n";
}
## side chain dihedral
print $df "\n";
print $df "Impropers\n";
print $df "\n";
for (1..@{$Sdihobj_a}){#@{$Sdihobj_a}
	my $Sdihcounter = $_ - 1;
	print "${$link_arraya}[$Sdihcounter]\n";
    my $type = scalar ${$Sdihname2type_h}{${$link_arraya}[$Sdihcounter]};# name to type ID
    my @SdihID =  @{$Sdihobj_a->[$Sdihcounter]}[0..3];#from 0
	print $df "$_ $type @SdihID" . " #${$link_arraya}[$Sdihcounter] for sidechain\n";
	#}
}

close($df);

}# end sub
1;
