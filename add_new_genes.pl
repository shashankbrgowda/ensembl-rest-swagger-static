#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;
use Bio::EnsEMBL::Registry;

my $registry = "Bio::EnsEMBL::Registry";
$registry->load_all();

my $slice_adaptor = $registry->get_adaptor('Human', 'Core', 'Slice');
my $attr_adaptor = $registry->get_adaptor('Human', 'Core', 'Attribute');
my $gene_adaptor = $registry->get_adaptor('Human', 'Core', 'Gene');

# ex: Chromosome
print "Enter coordinate system name:\n";  
my $csName = <STDIN>;
chomp $csName;

# ex: 21
print "Enter sequence region name:\n";  
my $seqRegionName = <STDIN>;
chomp $seqRegionName;

# ex: 36948925
print "Enter sequence region start:\n";  
my $seqRegionStart = <STDIN>;
chomp $seqRegionStart;

# ex: 37058704
print "Enter sequence region end:\n";  
my $seqRegionEnd = <STDIN>;
chomp $seqRegionEnd;

# ex: 1
print "Enter coordinate system version:\n";
my $csVersion = <STDIN>;
chomp $csVersion;

my $slice = $slice_adaptor->fetch_by_region($csName, $seqRegionName, $seqRegionStart, $seqRegionEnd, $csVersion);
my @genes = @{$slice->get_all_Genes};

my %attrib_gene;

# ex:
# {
# 	HCLS => [Gene1(554079), Gene2(572275)],
# 	HCLS1 => [Gene3]
# }
foreach my $gene (@genes) {
	my ($attrib) = @{ $gene->get_all_Attributes('name') };

	if(exists $attrib_gene{$attrib->value()}) {
		push(@{$attrib_gene{$attrib->value()}}, $gene);
	} else {
		$attrib_gene{$attrib->value()} = [$gene];
	}
}

foreach my $key ( keys %attrib_gene ) { 
   my @dup_genes = @{$attrib_gene{$key}};
   my $size = scalar(@dup_genes);
   my $c = 1;

   if($size > 1) {
   		foreach my $dup_gene (@dup_genes) {
   			$dup_gene->is_current(0);
   			$dup_gene->modified_date(time());
   			$gene_adaptor->update($dup_gene);

   			my $gene_id = $dup_gene->dbID();
   			$dup_gene->is_current(1);
   			my ($attrib) = @{ $dup_gene->get_all_Attributes('name') };
   			$attrib->value($attrib->value() . "." . $c);

				foreach my $t (@{$dup_gene->get_all_Transcripts}) {
				  my $tt = $t->translation();
				  if ($tt) {
				    foreach my $pf (@{$tt->get_all_ProteinFeatures}){
				      $pf->dbID(undef);
				      $pf->adaptor(undef);
				    }
				    $tt->dbID(undef);
				    $tt->adaptor(undef);
				  }
				  foreach my $e (@{$t->get_all_Exons}) {
				    $e->dbID(undef);
				    $e->adaptor(undef);
				  }
				  $t->dbID(undef);
				  $t->adaptor(undef);
				}

   			$dup_gene->dbID(undef);
   			$dup_gene->adaptor(undef);
   			$dup_gene->created_date(time());
   			$dup_gene->modified_date(time());

   			my $dbId = $gene_adaptor->store($dup_gene);
   			print "Gene ", $gene_id, " has duplicate name and new gene ", $dbId, " is created.\n";

   			$c = $c + 1;
   		}
   }
}