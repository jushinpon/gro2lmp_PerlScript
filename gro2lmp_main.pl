use warnings;
use strict;
use Data::Dumper;
require './output_datafile.pl';
require './output_potential.pl';
my $molecule_itp = "Nucleic_A.itp";
my $forcefield_itp = "martini_v2.1-dna.itp";# for non-bond parameters
my $input = "box.gro"; # pdb or gro

my $filenamein = "ssDNA";
my $filenamedata = "ssDNA";
my $forceconstant = "200";
my $kj2kcal = "0.2390057361";
my $pairstyle = "lj/gromacs/coul/gromacs 9.0 11.0  0.000001 11.0";
my $water = "no";
############################################################################
open my $itp1,"< $molecule_itp" or die("Can't open $molecule_itp");
my @mol_itp = <$itp1>;
close($itp1);
chomp @mol_itp;
########################################################################## ATOMS #############################################
#  (1)   (SN0)     (1)    DC   BB2     2  (0.0000) (45.0000) ; C
my @atomobj = grep {if(m/^\s*(\d+)\s+(\w+)\s+(\d+)\s+\w+\s+\w+\s+\d+\s+(-?\d*\.?\d*)\s+(-?\d*\.*\d*)\s*;?.+/)
	{$_ = [$1,$2,$3,$4,$5];}} @mol_itp;
if(! @atomobj) {die "No atom information\n";}
#print Dumper(\@atomobj);
my %atom_lookup;
for (@atomobj){
	$atom_lookup{$_->[1]} = $_->[4];# filter out the duplicate one and keep mass
}
#print Dumper(\%atom_lookup);
my @atomtype = sort keys %atom_lookup;# define the type order for lammps
my %atomtypeID = map {$atomtype[$_ - 1] => $_ } 1..@atomtype;
#print Dumper(\%atomtypeID);

open my $groinput,"< $input" or die("Can't open $input");
my @gro = <$groinput>;
close($groinput);
chomp @gro;
#1DC     BB2    1  (34.660)  (35.277)  (35.017)
my @coorobj = grep {if(m/^\s*\w+\s+\w+\s+\w+\s+(-?\d*\.?\d*)\s+(-?\d*\.?\d*)\s+(-?\d*\.?\d*)\s*$/){$_ = [$1*10.,$2*10.,$3*10.];}} @gro;
if(! @coorobj) {die "No coorobj information\n";}
#print Dumper(\@coorobj);
# gro box information
my @grobox = grep {if(m/^\s*(\d*\.\d*)\s+(\d*\.\d*)\s+(\d*\.\d*)\s*$/){$_ = [$1*10.,$2*10.,$3*10.];}} @gro;
if(!@grobox) {die "No grobox information\n";}
#print Dumper(\@grobox);
######################################################## BONDS ##########################################
#    (1)     (2)      (1)   (0.19800) (80000) ; DC(C)-DC(C)
my @bondobj = grep {if(m/^\s*(\d+)\s+(\d+)\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.?\d+)\s+;.+/ )
	{$_ = [$1,$2,$3,$4*10.,$5*$kj2kcal*0.5*0.01];}} @mol_itp;
if(! @bondobj) {die "No bondobj information\n";}
#print Dumper(\@bondobj);
my %bond_lookup; #look up table for bonding parameters
my %bondIDpair; # filter the duplicate pairs for exclusion pairs
my $bc = 1; #bond counter
my %bondID2name;
for (@bondobj){
	my $a = $_->[0] - 1;# the last one after push, 
	my $b = $_->[1] - 1;# array ID of the second bead
	$bondIDpair{"$_->[0]-$_->[1]"} = 1;
	$bondIDpair{"$_->[1]-$_->[0]"} = 1;
	my @temp1 = sort ($atomobj[$a][1],$atomobj[$b][1]);
	my @temp = (@temp1,$_->[4],$_->[3]);
	my $link = join ("-",@temp);
	$bond_lookup{"$link"} = [$_->[2],$_->[4],$_->[3]]; #hash -> array, for keeping function type and parameters
	$bondID2name{"$bc"} = "$link";
	$bc++;
}
#print Dumper(\%bond_lookup);
#print Dumper(\%bondID2name);
my @bondtype2name = sort keys %bond_lookup;
my %bondname2type = map {$bondtype2name[$_ - 1] => $_ } 1..@bondtype2name;

#constrained bond (need to use restrain in lammps, 
#outout variable xxx index xxx)
#3     4      1   0.22000 ; DC
my @cbond = grep {if(m/^\s*(\d+)\s+(\d+)\s+\d+\s+(-?\d*\.?\d*)\s*;.+/){$_ = [$1,$2,$3*10.0];}} @mol_itp;
if(! @cbond) {die "No constrained bond information\n";}

#exclusion pair (make fake bonds in lammps, but need to remove the duplicate ones from bond types!)
#    1     6   ; DC
my @temp = grep {if(m/^\s+(\d+)\s+(\d+)\s+;?\s*\w{0,2}$/){$_ = [$1,$2];}} @mol_itp;
if(! @temp) {die "No exclusion bond information\n";}
my @ebond;
for (@temp){
	my $a = $_->[0];
	my $b = $_->[1];
	chomp ($a,$b);
	my $link = join ("-",sort ($atomobj[$a-1][1],$atomobj[$b-1][1]));
	if(!$bondIDpair{"$a-$b"}){push @ebond,[$a,$b,"0.00 2.5","#$link for fake bond"];}
	else {print "$a-$b pair exists in both bond and ebond => skipped in ebond!\n"}
}
#print Dumper(\@ebond);

########################################################### ANGLES ##############################################
##    (1)     (2)     (3)      (2)  (95.00000)   (210) ; DC                         Sidechain angles
my @angleobj = grep {if(m/^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(-?\d*\.?\d*)\s+(-?\d*\.?\d*)\s*;.*$/)
	{$_ = [$1,$2,$3,$4,$5,$6*$kj2kcal*0.5];}} @mol_itp;
if(! @angleobj) {die "No angleobj information\n";}
#print Dumper(\@angleobj);
my $ac = 1; #angle counter
my %angleID2name; 
my %angle_lookup; #look up table for bending parameters
for (@angleobj){
	my $a = $_->[0] - 1;# array ID of the first bead 
	my $b = $_->[1] - 1;# array ID of the second bead
	my $c = $_->[2] - 1;# array ID of the third bead
	my @temp = ($atomobj[$a][1],$atomobj[$b][1],$atomobj[$c][1]);
	my $link;
	if($atomobj[$a][1] le $atomobj[$c][1]){
		
		my @tempa = (@temp,$_->[3],$_->[5],$_->[4]);
		$link = join ("-",@tempa);		
		}
		
	else {
		my @tempa = ((reverse @temp),$_->[3],$_->[5],$_->[4]);
		$link = join ("-",@tempa);
	}
	$angle_lookup{"$link"} = [$_->[3],$_->[5],$_->[4]]; #hash -> array, for keeping parameters
	$angleID2name{"$ac"} = "$link";
	$ac++;
}
#print Dumper(\%angle_lookup);
my @angletype2name = sort keys %angle_lookup;
my %anglename2type = map {$angletype2name[$_ - 1] => $_ } 1..@angletype2name;

############################################################Backbone DIHEDRALS ##########################################################
#   (1)     (2)     (6)     (7)      (1)  (180.00000)   (2)  (3) ; DC(C)-DC(C)-DG(C)-DG(C)       Backbone dihedrals
# need to consider Fourier type
my @Bdihobj4all = grep {if(m/^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(-?\d*\.\d*)\s+(-?\d+.?\d*)\s+(-?\d+.?\d*)\s*;.+/)
	{$_ = [$1,$2,$3,$4,$5,$6,$7*$kj2kcal,$8];}} @mol_itp;
if(! @Bdihobj4all) {die "No Backbone dihedral information\n";}
my @Bdihobj;# only keep one dihedral IDs for fourier type with multiplicity higher than 1
my %Bdih_lookup; #look up table for dihedral parameters
my $Bdc = 1; #backbone dih counter
my %BdihID2name;

for (@Bdihobj4all){
	my $a = $_->[0] - 1;# array ID of the first bead 
	my $b = $_->[1] - 1;# array ID of the second bead
	my $c = $_->[2] - 1;# array ID of the third bead
	my $d = $_->[3] - 1;# array ID of the fourth bead
	my @set_temp = ($a + 1,$b + 1 ,$c + 1,$d + 1);#dihedral set
	my @temp = ($atomobj[$a][1],$atomobj[$b][1],$atomobj[$c][1],$atomobj[$d][1]);
	my $link;
	my @temp_para = ($_->[4],$_->[6],$_->[7],$_->[5]);#current parameters from Bdihobj

	if($atomobj[$a][1] le $atomobj[$d][1]){
		if($_->[4] eq 9){
			$link = join ("-",@temp,@set_temp);
		}else {
			$link = join ("-",@temp,@temp_para);
		}
	}
	else {
		if($_->[4] eq 9){
			$link = join ("-",(reverse @temp),@set_temp);
		}else {
			$link = join ("-",(reverse @temp,@temp_para));
		}
	}
# if fourier type exists. lmp format m (by array No.), k, n,degress
	if($Bdih_lookup{"$link"}){# already exists and real Fourier
		    my $check = 0;
		for my $temp (@{$Bdih_lookup{"$link"}}){	
			if(@temp_para ~~ @$temp){$check = 1;}# if two arrays are equal
		}
		if($check == 0){push @{$Bdih_lookup{"$link"}}, [@temp_para];}# no duplicate parameters for the same key
	}
	else{
		$Bdih_lookup{"$link"} = [[@temp_para]]; #hash -> 2D array, for keeping parameters possible with Fourier
		push @Bdihobj,[@{$_}[0..3]];
		$BdihID2name{"$Bdc"} = "$link";
	    $Bdc++;
	}
}
#print Dumper(\@Bdihobj4all);
#print Dumper(\@Bdihobj);
#print Dumper(\%Bdih_lookup);
my @Bdihtype2name = sort keys %Bdih_lookup;
my %Bdihname2type = map {$Bdihtype2name[$_ - 1] => $_ } 1..@Bdihtype2name;

#for (@Bdihtype){my $multi = @{$Bdih_lookup{$_}}; print "Bdih Multiplicity for $_: $multi\n"}

############################################################ Sidechain DIHEDRALS ##########################################################
##   (1)     (2)     (3)     (4)      (2)  (-90.00000)    (20) ; DC                                  Sidechain dihedrals
# need to consider Fourier type
my @Sdihobj4all = grep {if(m/^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(-?\d*\.\d*)\s+(-?\d+.?\d*)\s+\s*;.+/)
	{$_ = [$1,$2,$3,$4,$5,abs($6),$7*$kj2kcal*0.5];}} @mol_itp;
if(! @Sdihobj4all) {die "No Sidechain dihedral information\n";}
my @Sdihobj;
my %Sdih_lookup; #look up table for dihedral parameters
my $Sdc = 1; #side chain dih counter
my %SdihID2name;
for (@Sdihobj4all){
	my $a = $_->[0] - 1;# array ID of the first bead 
	my $b = $_->[1] - 1;# array ID of the second bead
	my $c = $_->[2] - 1;# array ID of the third bead
	my $d = $_->[3] - 1;# array ID of the fourth bead
	my @set_temp = ($a + 1,$b + 1 ,$c + 1,$d + 1);#dihedral set

	my @temp = ($atomobj[$a][1],$atomobj[$b][1],$atomobj[$c][1],$atomobj[$d][1]);
	my $link;
	my @temp_para = ($_->[4],$_->[6],$_->[5]);#current parameters from Bdihobj

	if($atomobj[$a][1] le $atomobj[$d][1]){
		if($_->[4] eq 9){
			$link = join ("-",@temp,@set_temp);
		}else {
			$link = join ("-",@temp,@temp_para);
		}
	}
	else {
		if($_->[4] eq 9){
			$link = join ("-",(reverse @temp),@set_temp);
		}else {
			$link = join ("-",(reverse @temp,@temp_para));
		}
	}
	# if fourier type exists.
	if($Sdih_lookup{"$link"}){# already exists
		    my $check = 0;
		for my $temp (@{$Sdih_lookup{"$link"}}){	
			if(@temp_para ~~ @$temp){$check = 1;}# if two arrays are equal
		}
		if($check == 0){push @{$Sdih_lookup{"$link"}}, [@temp_para];}# no duplicate parameters for the same key
	}
	else{
		$Sdih_lookup{"$link"} = [[@temp_para]]; #hash -> 2D array, for keeping parameters possible with Fourier
	    push @Sdihobj,[@{$_}[0..3]];
	    $SdihID2name{"$Sdc"} = "$link";
	    $Sdc++;
	}
}
#print Dumper(\@Sdihobj);
#print Dumper(\%SdihID2name);
#print Dumper(\%Sdih_lookup);
my @Sdihtype2name = sort keys %Sdih_lookup;
my %Sdihname2type = map {$Sdihtype2name[$_ - 1] => $_ } 1..@Sdihtype2name;
#for (@Sdihtype){my $multi = @{$Sdih_lookup{$_}}; print "Sdih Multiplicity for $_: $multi\n"}

# non-bonded
open my $itp2,"< $forcefield_itp" or die("Can't open $forcefield_itp");
my @forcefield_itp = <$itp2>;
close($itp2);
#(vSNa)     (vSQa)  (1)  (4.300000e-01) (3.000000e+00)
#[-+]?([0-9]*[.])?[0-9]+([eE][-+]?\d+)? 
my @nonbondobj4all = grep {if(m/^\s*(\w+)\s+(\w+)\s+(\d+)\s+([-+]?\d*\.?\d*[eE][-+]?\d+)\s+([-+]?\d*\.?\d*[eE][-+]?\d+)/)
	{$_ = [$1,$2,$3,$4*10.,$5*$kj2kcal];}} @forcefield_itp;
if(! @nonbondobj4all) {die "No nonbondobj4all information\n";}
#print Dumper(\@nonbondobj4all);
my %nonbond4all_lookup; #look up table for non-bond parameters
for (@nonbondobj4all){
	my $a = $_->[0];# the name of the first atom 
	my $b = $_->[1];# the name of the second atom
	my $link = join ("-",sort ($a,$b));
	$nonbond4all_lookup{"$link"} = [$_->[2],$_->[4],$_->[3],"\#$link"]; #hash -> array, for keeping function type and parameters
}
#print Dumper(\%bond_lookup);
my @nonbond4alltype = sort keys %nonbond4all_lookup;

## the following is for the current system:
my %nonbond_lookup;

for my $t1 (0..$#atomtype)
{
	for my $t2 ($t1..$#atomtype)
	{
			my $a = $atomtype[$t1];# the name of the first atom 
			my $b = $atomtype[$t2];# the name of the second atom
			my $type1 = $t1 + 1;
			my $type2 = $t2 + 1;			
			my $link = join ("-",sort ($a,$b));
			if($nonbond4all_lookup{"$link"}){
				$nonbond_lookup{"$type1 $type2"} = $nonbond4all_lookup{"$link"};
			}
			else{ die "No $link pair coeff in the database.\n";}  
	}
}

#print Dumper(\%nonbond_lookup);
# make data file
&output_datafile(\@atomobj,\%atom_lookup,\@atomtype,\%atomtypeID,\@coorobj,\@grobox,
\@bondobj,\%bond_lookup,\@bondtype2name,\%bondname2type,\%bondID2name,
\@cbond,\@ebond,# other bond treatments
\@angleobj,\%angle_lookup,\@angletype2name,\%anglename2type,\%angleID2name,# information about bending
\@Bdihobj,\%Bdih_lookup,\@Bdihtype2name,\%Bdihname2type,\%BdihID2name,# information about backbone dihedral
\@Sdihobj,\%Sdih_lookup,\@Sdihtype2name,\%Sdihname2type,\%SdihID2name,#information about side chain dihedral
\%nonbond_lookup # non-bond information
);
#make potential file
&output_potential(\@atomobj,\%atom_lookup,\@atomtype,\%atomtypeID,\@coorobj,\@grobox,
\@bondobj,\%bond_lookup,\@bondtype2name,\%bondname2type,\%bondID2name,
\@cbond,\@ebond,# other bond treatments
\@angleobj,\%angle_lookup,\@angletype2name,\%anglename2type,\%angleID2name,# information about bending
\@Bdihobj,\%Bdih_lookup,\@Bdihtype2name,\%Bdihname2type,\%BdihID2name,# information about backbone dihedral
\@Sdihobj,\%Sdih_lookup,\@Sdihtype2name,\%Sdihname2type,\%SdihID2name,#information about side chain dihedral
\%nonbond_lookup # non-bond information
);
