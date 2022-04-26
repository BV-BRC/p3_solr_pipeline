#!/usr/bin/env perl

###########################################################
#
# rast2solr.pl: 
#
# Script to parse RASTtk genome object and input genbank file
# and prepare separate JSON files for each of the Solr cores.  
#
# Input: 
#	- JSON file containing RASTtk genome object
#	- Optional Original GenBank file used as input to RASTtk job 
# 
# Output: Eight JSON files each correspodning to a Solr core
#	- genome.json
#	- genome_sequence.json
#	- genome_feature.json
#	- feature_sequence.json
#	- pathway.json
#	- subsystem.json
#	- sp_gene.json
# - genome_amr.json  
#
###########################################################

use strict;
use FindBin qw($Bin);
use Getopt::Long::Descriptive;
use POSIX;
use JSON;
use Data::Dumper;
use Digest::MD5 qw(md5 md5_hex md5_base64);
use Bio::SeqIO;
use Bio::PrimarySeq;
use Bio::SeqFeature::Generic;
use Date::Parse;
use XML::Simple;
use File::Slurp;

use lib "$Bin/../lib";
use SolrAPI;
use P3RateLimitedUserAgent;
use Bio::P3::NCBILookup::NCBILookupClient;

our $have_config;
eval
{
    require Bio::KBase::AppService::AppConfig;
    $have_config = 1;
};

our $user_agent = P3RateLimitedUserAgent->new(3);

my ($data_api_url, $reference_data_dir);

my $data_api_url = $ENV{PATRIC_DATA_API};
my $reference_data_dir = $ENV{PATRIC_REFERENCE_DATA};

if (!defined($data_api_url) && $have_config)
{
    $data_api_url = Bio::KBase::AppService::AppConfig->data_api_url;
}

if (!defined($reference_data_dir) && $have_config)
{
    $reference_data_dir = Bio::KBase::AppService::AppConfig->reference_data_dir;
}

my $json = JSON->new->allow_nonref;

my ($opt, $usage) = describe_options("%c %o",
				     [],
				     ["write-reference-data", "Query reference data from Solr and write reference files"],
				     ["genomeobj-file=s", "RASTtk annotations as GenomeObj.json file"],
				     ["genbank-file=s", "Original GenBank file that was used as input to RASTtk"],
				     ["public", "public, default is private"],
				     ["data-api-url=s", "Data API URL", { default => $data_api_url }],
				     ["reference-data-dir=s", "Data API URL", { default => $reference_data_dir }],
				     ["lookup-client-url=s", "NCBI lookup client url"],
				     [],
				     ["help|h", "Print usage message and exit"] );

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless $opt->write_reference_data || $opt->genomeobj_file;
my $lookup_client = Bio::P3::NCBILookup::NCBILookupClient->new($opt->lookup_client_url) if $opt->lookup_client_url;
my $solrh = SolrAPI->new($opt->data_api_url, $opt->reference_data_dir);

if ($opt->write_reference_data)
{
    $solrh->getECRef();
    $solrh->getPathwayRef();
    $solrh->getSpGeneRef();
    exit 0;
}

my $genomeobj_file = $opt->genomeobj_file;
my $genbank_file = $opt->genbank_file;
my $outfile = $genomeobj_file;
$outfile=~s/(.gb|.gbf|.json)$//;

print "Processing $genomeobj_file\n";

# Read GenomeObj
open GOF, $genomeobj_file or die "Can't open input RASTtk GenomeObj JSON file: $genomeobj_file\n";
my @json_array = <GOF>;
my $json_str = join "", @json_array;
close GOF;

my $genomeObj = $json->decode($json_str);

# Set global parameters
my $public = $opt->public? 1 : 0;
my $annotation = "PATRIC"; 


# Get EC, pathway, and spgene reference dictionaries
my $ecRef = $solrh->getECRef();
my $pathwayRef = $solrh->getPathwayRef();
my $spgeneRef = $solrh->getSpGeneRef();


# Initialize global arrays to hold genome data
my %seq=();
my $genome;
my @sequences = ();
my @features = ();
my @feature_sequences = ();
my @pathwaymap = ();
my @subsystemmap = ();
my @spgenemap = (); 
my @taxonomy = ();
my @genome_amr = ();
my $featureIndex;
my %subsystem_assignments=();
my %md5=();


# Process GenomeObj
getGenomeInfo();
getGenomeQuality();

# Get additional genome metadata 
getMetadataFromGenBankFile($genbank_file) if -f $genbank_file;

getMetadataFromBioSample($genome->{biosample_accession}) if $genome->{biosample_accession} && $genome->{superkingdom} ne "Viruses";

# Auto curate genome metadata
curateMetadata();

# Get predicted AMR phenotypes 
getAMRPhenotypes();

# Get genome sequences
getGenomeSequences();

# Get subsystems
getSubsystems();

# Get Genome features from genome obj
getGenomeFeatures();

# Get additional features from the GenBank file
getGenomeFeaturesFromGenBankFile($genbank_file) if -f $genbank_file;

# delete $genome->{superkingdom};

# write to json files
writeJson();


sub writeJson {

	print "Preparing JSON files ...\n";

	my $genome_json = $json->pretty->encode($genome);
	my $sequence_json = $json->pretty->encode(\@sequences);
	my $feature_json = $json->pretty->encode(\@features);
	my $feature_sequence_json = $json->pretty->encode(\@feature_sequences);
	my $pathwaymap_json = $json->pretty->encode(\@pathwaymap);
	my $subsystemmap_json = $json->pretty->encode(\@subsystemmap);
	my $spgenemap_json = $json->pretty->encode(\@spgenemap);
	my $taxonomy_json = $json->pretty->encode(\@taxonomy);
	my $genome_amr_json = $json->pretty->encode(\@genome_amr);
	
	open FH, ">genome.json" or die "Cannot write genome.json: $!"; 
	print FH "[".$genome_json."]";
	close FH;

	open FH, ">genome_sequence.json" or die "Cannot write genome_sequence.json: $!"; 
	print FH $sequence_json;
	close FH;

	open FH, ">genome_feature.json" or die "Cannot write genome_feature.json: $!"; 
	print FH $feature_json;
	close FH;
	
	open FH, ">feature_sequence.json" or die "Cannot write feature_sequence.json: $!"; 
	print FH $feature_sequence_json;
	close FH;

	open FH, ">pathway.json" or die "Cannot write pathway.json: $!"; 
	print FH $pathwaymap_json;
	close FH;
	
	open FH, ">subsystem.json" or die "Cannot write subsystem.json: $!"; 
	print FH $subsystemmap_json;
	close FH;

	open FH, ">sp_gene.json" or die "Cannot write sp_gene.json: $!"; 
	print FH $spgenemap_json;
	close FH;

	open FH, ">taxonomy.json" or die "Cannot write taxonomy.json: $!"; 
	print FH $taxonomy_json;
	close FH;

	open FH, ">genome_amr.json" or die "Cannot write genome_amr.json: $!"; 
	print FH $genome_amr_json;
	close FH;

}

sub getGenomeInfo {

	print "Getting genome metadata ...\n";

	my ($chromosomes, $plasmids, $contigs, $sequences, $cds, $genome_length, $gc_count, $taxon_lineage_ranks);

	$genome->{owner} = $genomeObj->{owner}? $genomeObj->{owner} : "PATRIC";
	$genome->{public} = $public;

	$genome->{genome_id} = $genomeObj->{id};
	$genome->{genome_name} = $genomeObj->{scientific_name};
	$genome->{common_name} = $genomeObj->{scientific_name};
	$genome->{common_name}=~s/\W+/_/g;
	$genome->{common_name}=~s/_*$//g;

	$genome->{taxon_id}    =  $genomeObj->{ncbi_taxonomy_id};
	($genome->{taxon_lineage_ids}, $genome->{taxon_lineage_names}, $taxon_lineage_ranks)  = $solrh->getTaxonLineage($genome->{taxon_id});

	prepareTaxonomy($genome->{taxon_lineage_ids}) if $public;

	my %desired_ranks = map { $_ => 1 } qw(superkingdom kingdom phylum class order family genus species);
	if (ref $taxon_lineage_ranks eq "ARRAY"){
		for (my $i = 0; $i < @$taxon_lineage_ranks; $i++){
			my $rank = $taxon_lineage_ranks->[$i];
			my $name = $genome->{taxon_lineage_names}->[$i];
			if ($desired_ranks{lc($rank)}){
				$genome->{$rank} = $name;
			}
		}
	}

	foreach my $type (@{$genomeObj->{typing}}){
		$genome->{mlst} .= "," if $genome->{mlst};
		$genome->{mlst} .= $type->{typing_method}.".".$type->{database}.".".$type->{tag};
	}

	foreach my $seqObj (@{$genomeObj->{contigs}}) {

		$sequences++;

		if ($sequences == 1){
			foreach my $dblink (@{$seqObj->{genbank_locus}->{dblink}}){
				$genome->{bioproject_accession} = $1 if $dblink=~/BioProject:\s*(.*)/;
				$genome->{biosample_accession} = $1 if $dblink=~/BioSample:\s*(.*)/;
			}

			foreach my $reference (@{$seqObj->{genbank_locus}->{references}}){
				$genome->{publication} .= $reference->{PUBMED}."," unless $genome->{publication}=~/$reference->{PUBMED}/;
			}
			$genome->{publication}=~s/,*$//g;
			$genome->{completion_date} = strftime "%Y-%m-%dT%H:%M:%SZ", localtime str2time($seqObj->{genbank_locus}->{date}) if $seqObj->{genbank_locus}->{date};
		}

		if ($seqObj->{genbank_locus}->{accession}[1]=~/^([A-Z]{4})\d+00000$/){ # wgs, capture only master accession
			$genome->{genbank_accessions} = $seqObj->{genbank_locus}->{accession}[1];
		}else{
			$genome->{genbank_accessions} .= $seqObj->{genbank_locus}->{accession}[0]."," unless length($genome->{genbank_accessions}) > 100 ;
		}

	}
	$genome->{genbank_accessions}=~s/,*$//g;
	
	$genome->{assembly_accession}=$1 if $genbank_file=~/(GCA_\d+|GCF_\d+)/;

	
}


sub getGenomeQuality {

	print "Getting genome quality ...\n";

	my $qc = $genomeObj->{quality}; 

	$genome->{chromosomes} = $qc->{chromosomes};
	$genome->{plasmids} = $qc->{plasmids};
	$genome->{contigs} = $qc->{contigs};
	$genome->{genome_length} = $qc->{genome_length}; # Missing for viruses
	$genome->{gc_content} = $qc->{gc_content}; # Missing for viruses
	$genome->{contig_l50} = $qc->{genome_metrics}->{L50};
	$genome->{contig_n50} = $qc->{genome_metrics}->{N50};

	$genome->{genome_status} = $qc->{genome_status};
	$genome->{genome_status} = "Partial" if $genome->{superkingdom} eq "Viruses"; 

	$genome->{trna} = $qc->{feature_summary}->{tRNA};
	$genome->{rrna} = $qc->{feature_summary}->{rRNA};
	$genome->{cds} = $qc->{feature_summary}->{cds};

	$genome->{cds_ratio} = $qc->{cds_ratio};
	$genome->{hypothetical_cds} = $qc->{protein_summary}->{hypothetical};
	$genome->{hypothetical_cds_ratio} = $qc->{hypothetical_cds_ratio};
	$genome->{partial_cds} = $qc->{feature_summary}->{partial_cds};
	$genome->{partial_cds_ratio} = $qc->{partial_cds_ratio};
	$genome->{plfam_cds} = $qc->{protein_summary}->{plfam_assignment};
	$genome->{plfam_cds_ratio} = $qc->{plfam_cds_ratio};

	$genome->{coarse_consistency} = $qc->{coarse_consistency};
	$genome->{fine_consistency} = $qc->{fine_consistency};
	$genome->{checkm_completeness} = $qc->{completeness};
	$genome->{checkm_contamination} = $qc->{contamination};
	#$genome->{checkm_completeness} = $qc->{checkm_data}->{Completeness};
	#$genome->{checkm_contamination} = $qc->{checkm_data}->{Contamination};

	$genome->{genome_quality_flags} = $qc->{genome_quality_flags};
	$genome->{genome_quality} = $qc->{genome_quality};

}


sub getAMRPhenotypes {

	foreach my $amr1 (@{$genomeObj->{classifications}}){
	
		my $amr;
		next if $amr1->{name}=~/combined/;

		$amr->{owner} = $genome->{owner};
		$amr->{public} = $public;
		$amr->{genome_id} = $genome->{genome_id};
		$amr->{genome_name} = $genome->{genome_name};
		$amr->{taxon_id} = $genome->{taxon_id};

		$amr->{antibiotic} = lc $amr1->{name}; 	
		$amr->{resistant_phenotype} = ucfirst $amr1->{sensitivity};
		$amr->{resistant_phenotype} = "Susceptible" if $amr->{resistant_phenotype}=~/sensitive/i;	
		$amr->{evidence} = "Computational Prediction"; 	
		$amr->{computational_method} = "AdaBoost Classifier"; 	
		$amr->{computational_method_performance} = "Accuracy:$amr1->{accuracy}, F1 score:$amr1->{f1_score}, AUC:$amr1->{area_under_roc_curve}";
		$amr->{vendor} = "PATRIC"; 	

		push @{$genome->{antimicrobial_resistance}}, ucfirst $amr->{resistant_phenotype} unless (grep {$_ eq ucfirst $amr->{resistant_phenotype}} @{$genome->{antimicrobial_resistance}});
		$genome->{antimicrobial_resistance_evidence} = "Computational Prediction";

		push @genome_amr, $amr;

		}

}


sub getGenomeSequences {

	print "Getting genome sequences ...\n";

	my $count=0;

	foreach my $seqObj (@{$genomeObj->{contigs}}){
		
		$count++;
		my $sequence;
		
		$sequence->{owner} = $genome->{owner};
		$sequence->{public} = $public;
	
		$sequence->{genome_id} = $genome->{genome_id};
		$sequence->{genome_name} = $genome->{genome_name};
		$sequence->{taxon_id} = $genome->{taxon_id};
		$sequence->{sequence_id} = $sequence->{genome_id}.".con.".sprintf("%04d", $count);

		my $seq_id = $seqObj->{id};	
		$sequence->{gi} = $seqObj->{genbank_locus}->{gi};

		if($seqObj->{genbank_locus}->{accession}[0]=~/\S+/){
			$sequence->{accession} = $seqObj->{genbank_locus}->{accession}[0];
		}elsif($seq_id=~/gi\|(\d+)\|(ref|gb)\|([\w\.]+)/){
			$sequence->{gi} = $1;
			$sequence->{accession} = $3;
		}elsif($seq_id=~/accn\|([\w\.]+)/){
			$sequence->{accession} = $1;
		}else{
			$sequence->{accession} = $sequence->{sequence_id};
		}

		$sequence->{description} = $seqObj->{genbank_locus}->{definition}? $seqObj->{genbank_locus}->{definition} : $seqObj->{id};
		$sequence->{topology} =	$seqObj->{genbank_locus}->{geometry};
		#$sequence->{mol_type} =	$seqObj->{genbank_locus}->{mol_type};

		$sequence->{sequence_type} = $1 if $sequence->{description}=~/(chromosome|plasmid|segment|contig|scaffold)/i;
		$sequence->{sequence_status} = $1 if $sequence->{description}=~/(complete|partial)/i;

		$genome->{genome_status} = "Complete" 
			if $sequence->{description}=~/complete genome|complete segment|assembly|complete sequence/i && $genome->{superkingdom} eq "Viruses";

		$sequence->{chromosome} = $1 if $sequence->{description}=~/chromosome (\S*)\s*,/i;
		$sequence->{plasmid} = $1 if $sequence->{description}=~/plasmid (\S*)\s*,/i;
		$sequence->{segment} = $1	if $seqObj->{description}=~/segment (\S*)\s*,/i;

		$sequence->{gc_content} = sprintf("%.2f", ($seqObj->{dna}=~tr/GCgc//)*100/length($seqObj->{dna})) if $seqObj->{dna};
		$sequence->{length} = length($seqObj->{dna});
		$sequence->{sequence} = lc($seqObj->{dna});
		$sequence->{sequence_md5} = md5_hex(lc $seqObj->{dna});
	
		$sequence->{version} = $1 if $seqObj->{genbank_locus}->{version}[0]=~/^.*\.(\d+)$/;
		$sequence->{release_date} = strftime "%Y-%m-%dT%H:%M:%SZ", localtime str2time($seqObj->{genbank_locus}->{date});
	
		# look up hash
		$seq{$seq_id}{sequence_id} = $sequence->{sequence_id};
		$seq{$seq_id}{accession} = $sequence->{accession};
		$seq{$seq_id}{sequence} = $sequence->{sequence};

		push @sequences, $sequence;

	}

}


sub getSubsystems{

	print "Getting subsystems ...\n";

	foreach my $subsystem (@{$genomeObj->{subsystems}}){
		my $subsystem_name = $subsystem->{name};
		my $active = $subsystem->{variant_code};
		my ($superclass, $class, $subclass) = @{$subsystem->{classification}};
		foreach my $role (@{$subsystem->{role_bindings}}){
			my $role_id = $role->{role_id};
			foreach my $patric_id (@{$role->{features}}){
				push @{$subsystem_assignments{$patric_id}}, "$patric_id\t$role_id\t$subsystem_name\t$superclass\t$class\t$subclass\t$active\n";
			}
		}
	}

}


sub getGenomeFeatures{
	
	print "Getting genome features ...\n";

	foreach my $featObj (@{$genomeObj->{features}}){
			
		my ($feature, $sequence, $aa_sequence, $pathways, $ecpathways);
		my (@segments, @go, @ec_no, @ec, @pathways, @ecpathways, @spgenes, @uniprotkb_accns, @ids);

		$feature->{owner} = $genome->{owner};
		$feature->{public} = $public;
		$feature->{annotation} = $annotation;
			
		$feature->{genome_id} = $genome->{genome_id};

		$feature->{genome_name} = $genome->{genome_name};
		$feature->{taxon_id} = $genome->{taxon_id};

		$feature->{patric_id} = $featObj->{id};

		$feature->{feature_type} = $featObj->{type};
		$feature->{feature_type} = $1 if ($featObj->{type} eq "rna" && $featObj->{function}=~/(rRNA|tRNA)/);
		$feature->{feature_type} = 'misc_RNA' if ($featObj->{type} eq "rna" && !($featObj->{function}=~/rRNA|tRNA/));
		$feature->{feature_type} = 'repeat_region' if ($featObj->{type} eq "repeat");

		$feature->{product} = $featObj->{function};
		$feature->{product} = "hypothetical protein" if ($feature->{feature_type} eq 'CDS' && !$feature->{product});
		$feature->{product}=~s/\"/''/g;
		$feature->{product}=~s/^'* *| *'*$//g;
		
		foreach my $locObj (@{$featObj->{location}}){
			
			my ($seq_id, $pstart, $start, $end, $strand, $length);

			($seq_id, $pstart, $strand, $length) = @{$locObj};
			
			if ($strand eq "+"){
				$start = $pstart;
				$end = $start+$length-1;
			}else{
				$end = $pstart;
				$start = $pstart-$length+1;
			}
			push @segments, $start."..". $end;

			$feature->{sequence_id} = $seq{$seq_id}{sequence_id};
			$feature->{accession} = $seq{$seq_id}{accession};

			$feature->{start} = $start if ($start < $feature->{start} || !$feature->{start});
			$feature->{end} = $end if ($end > $feature->{end} || !$feature->{end});
			$feature->{strand} = $strand;

			$sequence .= lc substr($seq{$seq_id}{sequence}, $start-1, $length);

		}

		$sequence =~tr/ACGTacgt/TGCAtgca/ if $feature->{strand} eq "-";
		$sequence = reverse($sequence) if $feature->{strand} eq "-";

		$feature->{segments} = \@segments;
		$feature->{location} = $feature->{strand} eq "+"? join(",", @segments): "complement(".join(",", @segments).")"; 

		if ($sequence && $feature->{feature_type} ne "source"){
			$feature->{na_length} = length($sequence);
			$feature->{na_sequence_md5} = md5_hex(lc $sequence);
			if ($md5{$feature->{na_sequence_md5}} !=1){
				$md5{$feature->{na_sequence_md5}} = 1;
				push @feature_sequences, { "md5" => $feature->{na_sequence_md5}, "sequence_type" =>"NA", "sequence" => $sequence };
			}
		}

		if ($featObj->{protein_translation}){
			$aa_sequence = $featObj->{protein_translation};
		}elsif($feature->{feature_type} eq "mat_peptide"){
			my $seqObj = Bio::PrimarySeq->new (-seq => $sequence);
			my $translatedSeqObj = $seqObj->translate();
			$aa_sequence = $translatedSeqObj->seq();	
		}else{
			# No AA sequence
		}

		if ($aa_sequence){
			$feature->{aa_length} 	= length($aa_sequence);
			$feature->{aa_sequence_md5} = md5_hex($aa_sequence);	
			if ($md5{$feature->{aa_sequence_md5}} !=1){
				$md5{$feature->{aa_sequence_md5}} = 1;
				push @feature_sequences, { "md5" => $feature->{aa_sequence_md5}, "sequence_type" =>"AA", "sequence" => $aa_sequence };
			}
		}

		my $strand = ($feature->{strand} eq '+')? 'fwd':'rev';
		$feature->{feature_id}		=	"$annotation.$feature->{genome_id}.$feature->{accession}.".
																	"$feature->{feature_type}.$feature->{start}.$feature->{end}.$strand";

		if ($feature->{feature_type}=~/classifier_predicted_region/){
			my $evidence =  $featObj->{annotations}->[0]->[0];
			($feature->{classifier_score}, $feature->{classifier_round}) = $evidence=~/alpha=(\S*) round=(\S*)/;
		}


		foreach my $alias (@{$featObj->{alias_pairs}}){
			my ($alias_type, $alias_value) = @{$alias};
			$feature->{refseq_locus_tag} = $alias_value if ($alias_type eq "locus_tag");
			$feature->{protein_id} = $alias_value if ($alias_type eq "");
			$feature->{gene_id} = $alias_value if ($alias_type eq "GeneID");
			$feature->{gene} = $alias_value if ($alias_type eq "gene");
		}

		foreach my $family (@{$featObj->{family_assignments}}){
			my ($family_type, $family_id, $family_function) = @{$family};
			$feature->{figfam_id} = $family_id if ($family_id=~/^FIG/);
			$feature->{plfam_id} = $family_id if ($family_id=~/^PLF/);
			$feature->{pgfam_id} = $family_id if ($family_id=~/^PGF/);
		}

		@ec_no = $feature->{product}=~/\( *EC[: ]([\d-\.]+) *\)/g if $feature->{product}=~/EC[: ]/;
		push @{$feature->{property}}, "EC number" if scalar @ec_no > 0;
		
		foreach my $ec_number (@ec_no){

			my $ec_description = $ecRef->{$ec_number}->{ec_description};
			push @ec, $ec_number.'|'.$ec_description unless (grep {$_ eq $ec_number.'|'.$ec_description} @ec);
			
			foreach my $go_term (@{$ecRef->{$ec_number}->{go}}){
				push @go, $go_term unless (grep {$_ eq $go_term} @go);
			}
			
			foreach my $pathway (@{$pathwayRef->{$ec_number}->{pathway}}){
				my ($pathway_id, $pathway_name, $pathway_class) = split(/\t/, $pathway);
				push @pathways, $pathway_id.'|'.$pathway_name unless (grep {$_ eq $pathway_id.'|'.$pathway_name} @pathways);
				my $ecpathway = "$ec_number\t$ec_description\t$pathway_id\t$pathway_name\t$pathway_class";
				push @ecpathways, $ecpathway unless (grep {$_ eq $ecpathway} @ecpathways);
			}

		}

		$feature->{go} = \@go if scalar @go;
		push @pathwaymap, preparePathways($feature, \@ecpathways);
		push @{$feature->{property}}, "Pathway" if scalar @ecpathways > 0;

		@spgenes = @{$featObj->{similarity_associations}} if $featObj->{similarity_associations};			
		push @spgenemap, prepareSpGene($feature, $_) foreach(@spgenes);

		# Prepare PATRIC AMR genes for matching functions 
		push @spgenemap, prepareSpGene($feature, ()) if $spgeneRef->{$feature->{product}};

		push @{$feature->{property}}, "Subsystem" if ref $subsystem_assignments{$feature->{patric_id}} eq "ARRAY";	
		foreach my $subsystem_assignment (@{$subsystem_assignments{$feature->{patric_id}}}){
			my @values = split /\t/, $subsystem_assignment;
			my $subsystem_name = $values[2];
			$subsystem_name =~s/_/ /g;
		
			push @subsystemmap, prepareSubsystem($feature, $subsystem_assignment);	
		}

		push @features, $feature  unless $feature->{feature_type} eq 'gene'; 

		$genome->{lc($annotation).'_cds'}++ if $feature->{feature_type} eq 'CDS';

	}

}


sub prepareSpGene {

		my ($feature, $spgene_match) = @_;

		my $spgene;
		my ($property, $gene_name, $locus_tag, $organism, $function, $classification, $antibiotics_class, $antibiotics, $pmid, $assertion); 
		my ($source, $source_id, $qcov, $scov, $identity, $evalue);

		if($spgene_match){ # All specialty genes from external sources
			($source, $source_id, $qcov, $scov, $identity, $evalue) = @$spgene_match;
			$source_id=~s/^\S*\|//;
			($property, $gene_name, $locus_tag, $organism, $function, $classification, $antibiotics_class, $antibiotics, $pmid, $assertion) 
			= split /\t/, $spgeneRef->{$source.'_'.$source_id} if ($source && $source_id);
		}elsif($spgeneRef->{$feature->{product}}){ # PATRIC AMR genes, match by functional role
			($property, $gene_name, $locus_tag, $organism, $function, $classification, $antibiotics_class, $antibiotics, $pmid, $assertion) 
			= split /\t/, $spgeneRef->{$feature->{product}};
			$source = "PATRIC";	
		}

		return unless $property && $source;

		my ($qgenus) = $feature->{genome_name}=~/^(\S+)/;
		my ($qspecies) = $feature->{genome_name}=~/^(\S+ +\S+)/;
		my ($sgenus) = $organism=~/^(\S+)/;
		my ($sspecies) = $organism=~/^(\S+ +\S+)/;

		my ($same_genus, $same_species, $same_genome, $evidence); 

		$same_genus = 1 if ($qgenus eq $sgenus && $sgenus ne "");
		$same_species = 1 if ($qspecies eq $sspecies && $sspecies ne ""); 
		$same_genome = 1 if ($feature->{genome} eq $organism && $organism ne "") ;

		$evidence = ($source && $source_id)? 'BLAT' : "K-mer Search";

		$spgene->{owner} = $feature->{owner};
		$spgene->{public} = $public;
		
		$spgene->{genome_id} = $feature->{genome_id};	
		$spgene->{genome_name} = $feature->{genome_name};	
		$spgene->{taxon_id} = $feature->{taxon_id};	
		
		$spgene->{feature_id} = $feature->{feature_id};
		$spgene->{patric_id} = $feature->{patric_id};	
		$spgene->{alt_locus_tag} = $feature->{alt_locus_tag};
		$spgene->{refseq_locus_tag} = $feature->{refseq_locus_tag};
		
		$spgene->{gene} = $gene_name? $gene_name : $feature->{gene};
		$spgene->{product} = $feature->{product};

		$spgene->{property} = $property;
		$spgene->{source} = $source;
		$spgene->{property_source} = $property.': '.$source;
		
		$spgene->{source_id} = $source_id;
		$spgene->{organism} = $organism;
		$spgene->{function} = $function;
		$spgene->{classification} = $classification; 
		$spgene->{antibiotics_class} = $antibiotics_class if $antibiotics_class;
		$spgene->{antibiotics} = [split /[,;]/, $antibiotics] if $antibiotics;
		$spgene->{pmid} = [split /[,;]/, $pmid] if $pmid;
		$spgene->{assertion} = $assertion;

		$spgene->{query_coverage} = $qcov; 
		$spgene->{subject_coverage} =  $scov;
		$spgene->{identity} = $identity;
		$spgene->{e_value} = $evalue;

		$spgene->{same_genus} = $same_genus;
		$spgene->{same_species} = $same_species;
		$spgene->{same_genome} = $same_genome;
	  $spgene->{evidence} = $evidence;	

		return $spgene;
}


sub preparePathways {

	my ($feature, $ecpathways) = @_;
	my @pathways = ();

	foreach my $ecpathway (@$ecpathways){
		my $pathway;
		my ($ec_number, $ec_description, $pathway_id, $pathway_name, $pathway_class) = split /\t/, $ecpathway;
		
		$pathway->{owner} = $feature->{owner};
		$pathway->{public} = $public;

		$pathway->{genome_id} = $feature->{genome_id};
		$pathway->{genome_name} = $feature->{genome_name};
		$pathway->{taxon_id} = $feature->{taxon_id};

		$pathway->{sequence_id} = $feature->{sequence_id};
		$pathway->{accession} = $feature->{accession};
		
		$pathway->{annotation} = $feature->{annotation};
		
		$pathway->{feature_id} = $feature->{feature_id};
		$pathway->{patric_id} = $feature->{patric_id};
		$pathway->{alt_locus_tag} = $feature->{alt_locus_tag};
		$pathway->{refseq_locus_tag} = $feature->{refseq_locus_tag};
		
		$pathway->{gene} = $feature->{gene};
		$pathway->{product} = $feature->{product};
		
		$pathway->{ec_number} = $ec_number;
		$pathway->{ec_description} = $ec_description;
		
		$pathway->{pathway_id} = $pathway_id;
		$pathway->{pathway_name} = $pathway_name;
		$pathway->{pathway_class} = $pathway_class;
		
		$pathway->{genome_ec} = $feature->{genome_id}.'_'.$ec_number;
		$pathway->{pathway_ec} = $pathway_id.'_'.$ec_number;

		push @pathways, $pathway;

	}

	return @pathways;

}

sub prepareSubsystem {

	my ($feature, $subsystem_assignment) = @_;

	my $subsystem;

	chomp $subsystem_assignment;
	my ($patric_id, $role_id, $subsystem_id, $superclass, $class, $subclass, $active) = split /\t/, $subsystem_assignment;
	my $subsystem_name = $subsystem_id;
	$subsystem_id =~s/ /_/g;
	$subsystem_name =~s/_/ /g;
	my $role_name = $role_id;
	$role_id =~s/ /_/g;
	$role_name =~s/_/ /g;

	$subsystem->{owner} = $feature->{owner};
	$subsystem->{public} = $public;

	$subsystem->{genome_id} = $feature->{genome_id};
	$subsystem->{genome_name} = $feature->{genome_name};
	$subsystem->{taxon_id} = $feature->{taxon_id};

	$subsystem->{feature_id} = $feature->{feature_id};
	$subsystem->{patric_id} = $feature->{patric_id};
	$subsystem->{refseq_locus_tag} = $feature->{refseq_locus_tag};
		
	$subsystem->{gene} = $feature->{gene};
	$subsystem->{product} = $feature->{product};
		
	$subsystem->{role_id} = $role_id;
	$subsystem->{role_name} = $role_name;
		
	$subsystem->{subsystem_id} = $subsystem_id;
	$subsystem->{subsystem_name} = $subsystem_name;
		
	$subsystem->{superclass} = $superclass;
	$subsystem->{class} = $class;
	$subsystem->{subclass} = $subclass;
	
	$subsystem->{active} = $active;
		
	return $subsystem;

}


sub getMetadataFromGenBankFile {

	my ($genbak_file) = @_;

	print "Getting genome metadata from genbank file: $genbank_file ...\n";

	open GB, "<$genbank_file" || return "Can't open genbank file: $genbank_file\n";
	my $gb;
	my $in_contig = 0;
	while (<GB>){
		if ($in_contig){
			if (!/^\s*\d/){
		    $in_contig = 0;
			}
		}else{
			$gb .= $_;
			$in_contig = /^ORIGIN/;
	  }
	}
	close GB;
	
	$genome->{host_name} = "Patent" if $gb=~/LOCUS.*PAT/;
	if($gb=~/JOURNAL\s+Patent: *(.*)\s(\S+);/){
		my $patent_number = $1;
		my $patent_date = $2;
		$patent_number=~s/ (\d+)$/\/$1/;
	
		push @{$genome->{additional_metadata}}, "Patent Number: $patent_number";
		push @{$genome->{additional_metadata}}, "Patent Date: $patent_date";
	}

	# parse metadata from genabnk comment 
	$genome->{assembly_method} = $1 if $gb=~/Assembly Method\s*:: (.*)/;
	$genome->{sequencing_depth} = $1 if $gb=~/Genome Coverage\s*:: (.*)/;
	$genome->{sequencing_platform} = $1 if $gb=~/Sequencing Technology\s*:: (.*)/;

	# Parse flu metadata from GenBank comment 
	my $flu_data=$1 if $gb=~/##FluData-START##(.+?)##FluData-END##/s;
	foreach my $entry (split /\n/, $flu_data){
		next unless $entry=~/::/;
		my ($attrib,$value) = $entry=~/^\s*(\S+)\s+::\s+(.*)\s*$/;
		$attrib = lc $attrib;
		if ($attrib=~/host_age|host_gender|passage/i){
			$genome->{$attrib} = $value;
		}elsif($attrib=~/^collect_date$/i){
			$genome->{collection_date}=$value;
		}elsif($attrib=~/^type$/i){
			$genome->{serovar}=$value;
		}else{
			push @{$genome->{other_clinical}}, "$attrib:$value"; 
		}
	} 


	# Parse metadata from source feature 
	$gb=~s/\n */ /g;
	
	my $strain = ""; 
	if ($gb=~/\/strain="([^"]*)"/){
		$strain = $1;
	}elsif($gb=~/\/isolate="([^"]*)"/){
		$strain = $1;
	}else{
		# no strain or isolate
	}

	$genome->{strain} = $strain unless $strain=~/^ *(-|missing|na|n\/a|not available|not provided|not determined|nd|unknown) *$/i;
	
	$genome->{genome_name} .= " $genome->{strain}" if ($genome->{strain} && (not $genome->{genome_name}=~/$genome->{strain}/i));
	$genome->{genome_name}=~s/\($genome->{strain}\)/$genome->{strain}/;

	$genome->{segment} = $1 if $gb=~/\/segment="([^"]*)"/ && (grep {$_=~/Bunyavirales|Reoviridae|Orthomyxoviridae/} @{$genome->{taxon_lineage_names}});
	$genome->{serovar} = $1 if $gb=~/\/serotype="([^"]*)"/;
	$genome->{geographic_location} = $1 if $gb=~/\/country="([^"]*)"/;
	$genome->{host_name} = $1 if $gb=~/\/host="([^"]*)"/;
	$genome->{lab_host} = $1 if $gb=~/\/lab_host="([^"]*)"/;
	$genome->{isolation_source} = $1 if $gb=~/\/isolation_source="([^"]*)"/;
	$genome->{collection_date} = $1 if $gb=~/\/collection_date="([^"]*)"/;
	$genome->{culture_collection} = $1 if $gb=~/\/culture_collection="([^"]*)"/;

	if ($gb=~/\/note="(passage.details|passage.history) *: *([^"]*) *"/){
		$genome->{passage} = $1;
	}elsif($genome->{lab_host}=~/passage/i){
		$genome->{passage} = $genome->{lab_host};
	}else{
		# No passage info
	} 

	# For Influenza, get missing metadata from strain name	
	if ($genome->{genome_name}=~/Influenza (A|B|C|D)/){	
		if ($genome->{strain}=~/\s*\(([HN0-9-x]*)\)\s*$/){
			$genome->{serovar}=$1;
			$genome->{strain}=~s/\s*\([HN0-9-x]*\)\s*$//;
		}
		if ($genome->{strain}=~/(A|B|C|D)\/(.*?)\/(.*?)\/(.*?)\/(.*)/){ # type/host/location/identifier/year
			$genome->{host_name} = $2 unless $genome->{host_name};
			$genome->{geographic_location} = $3 unless $genome->{geographic_location};
			$genome->{collection_date} = $5 unless $genome->{collection_date};
		}elsif($genome->{strain}=~/(A|B|C|D)\/(.*?)\/(.*?)\/(.*?)/){ # type/location/identifier/year
			$genome->{host_name} = "Homo sapiens" unless $genome->{host_name};
			$genome->{geographic_location} = $2 unless $genome->{geographic_location};
			$genome->{collection_date} = $4 unless $genome->{collection_date};	
		}else{
			# strain is not expected format
		}
	}


}


sub getGenomeFeaturesFromGenBankFile {

	my ($genbank_file) = @_;

	print "Getting genome features from the genbank file: $genbank_file ...\n";
	
	my $annotation = "RefSeq";

	my $genomeObj = Bio::SeqIO->new( -file   => "<$genbank_file", -format => 'GenBank');

	while (my $seqObj = $genomeObj->next_seq){

		my ($accession, $sequence_id);

		$accession = $seqObj->accession_number;

		for (my $i=0; $i < scalar @sequences; $i++){
			next unless $sequences[$i]->{accession} eq $accession;
			$sequence_id = $sequences[$i]->{sequence_id}; 
			last;	
		}
	
		for my $featObj ($seqObj->get_SeqFeatures){

			my ($accn, $feature, $na_sequence, $aa_sequence, $pathways, $ecpathways);
			my (@go, @ec_no, @ec, @pathway, @ecpathways, @spgenes, @uniprotkb_accns, @ids);

			$feature->{owner} = $genome->{owner};
			$feature->{public} = $public;
			$feature->{annotation} = $annotation;
			
			$feature->{genome_id} = $genome->{genome_id};
			$feature->{genome_name} = $genome->{genome_name};
			$feature->{taxon_id} = $genome->{taxon_id};

			$feature->{sequence_id} = $sequence_id;
			$feature->{accession} = $accession;

			$feature->{feature_type} = $featObj->primary_tag;
			next if ($feature->{feature_type} eq 'gene' && (grep {$_=~/Bacteria|Archaea/} @{$genome->{taxon_lineage_names}}));
			$featureIndex->{$feature->{feature_type}}++;

			$feature->{start} = $featObj->start;
			$feature->{end} = $featObj->end;
			$feature->{strand} = ($featObj->strand==1)? '+':'-';
			$feature->{location} = $featObj->location->to_FTstring;

			my @segments;
			if ($featObj->location->isa('Bio::Location::SplitLocationI')){
				for my $location ($featObj->location->sub_Location){
					push @segments, $location->start."..". $location->end;
				}
			}else{
				push @segments, $featObj->start."..". $featObj->end;	
			}
			$feature->{segments} = \@segments;

			$na_sequence = lc $featObj->spliced_seq->seq;

			my $strand = ($feature->{strand} eq '+')? 'fwd':'rev';
			$feature->{feature_id}		=	"$annotation.$feature->{genome_id}.$feature->{accession}.".
																	"$feature->{feature_type}.$feature->{start}.$feature->{end}.$strand";
			
			for my $tag ($featObj->get_all_tags){

				for my $value ($featObj->get_tag_values($tag)){
					
					$feature->{codon_start} = $value if $tag eq 'codon_start';

					$feature->{feature_type} 	= 'pseudogene' if ($tag eq 'pseudo' && $feature->{feature_type} eq 'gene');

					$feature->{patric_id} = $1 if ($tag eq 'db_xref' && $value=~/SEED:(fig.*)/);					
					$feature->{refseq_locus_tag} 	= $value if ($tag eq 'locus_tag' && $annotation eq "RefSeq");
					$feature->{refseq_locus_tag} 	= $1 if ($tag eq 'db_xref' && $value=~/Refseq_locus_tag:(.*)/i);
					$feature->{protein_id} 	= $value if ($tag eq 'protein_id');
					$feature->{protein_id} 	= $1 if ($tag eq 'db_xref' && $value=~/protein_id:(.*)/);
					$feature->{gene_id} = $1 if ($tag eq 'db_xref' && $value=~/^GeneID:(\d+)/);

					$aa_sequence = $value if ($tag eq "translation");

					$feature->{gene} = $value if ($tag eq 'gene');
					$feature->{product} = $value if ($tag eq 'product');

					$feature->{figfam_id} 	= $value if ($tag eq 'FIGfam');
					
					if ($tag eq 'EC_number'){
						my $ec_number = $value;	
						my $ec_description = $ecRef->{$ec_number}->{ec_description};
						push @ec, $ec_number.'|'.$ec_description unless (grep {$_ eq $ec_number.'|'.$ec_description} @ec);
					}

					push @ids, $value if ($tag eq 'db_xref');
				
				}

			}

			next if ($feature->{feature_type} eq 'gene' && (grep {$_=~/Bacteria|Archaea/} @{$genome->{taxon_lineage_names}}));

			if ($na_sequence && $feature->{feature_type} ne "source"){
				$feature->{na_length} = length($na_sequence);
				$feature->{na_sequence_md5} = md5_hex(lc $na_sequence);
				if ($md5{$feature->{na_sequence_md5}} !=1){
					$md5{$feature->{na_sequence_md5}} = 1;
					push @feature_sequences, { "md5" => $feature->{na_sequence_md5}, "sequence_type" =>"NA", "sequence" => $na_sequence };				
				}
			}

			if($feature->{feature_type} eq "mat_peptide" && not $aa_sequence){
				my $seqObj = Bio::PrimarySeq->new (-seq => $na_sequence);
				my $translatedSeqObj = $seqObj->translate();
				$aa_sequence = $translatedSeqObj->seq();
			}

			if ($aa_sequence){
				$feature->{aa_length} 	= length($aa_sequence);
				$feature->{aa_sequence_md5} = md5_hex($aa_sequence);
				if ($md5{$feature->{aa_sequence_md5}} !=1){
					$md5{$feature->{aa_sequence_md5}} = 1;
					push @feature_sequences, { "md5" => $feature->{aa_sequence_md5}, "sequence_type" =>"AA", "sequence" => $aa_sequence };
				}
			}					
			
			push @features, $feature;
	
			$genome->{lc($annotation).'_cds'}++ if $feature->{feature_type} eq 'CDS';	

		}

	}

	$genomeObj->close();

}

sub getMetadataFromBioSample {

	my($biosample_accn) = @_;

	print "Getting genome metadata from BioSample: $biosample_accn ...\n";

	my($biosample_id, $biosample_xml);

	if ($lookup_client)
	{
	    ($biosample_id, $biosample_xml) = $lookup_client->get_metadata_from_biosample($biosample_accn);
	}
	else
	{	

	    my $xml = get_xml("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=biosample&term=$biosample_accn");
	    $xml =~ s/\n//;
	    ($biosample_id) = $xml=~/<Id>(\d+)<\/Id>/;

	    $biosample_xml = get_xml("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&retmode=xml&id=$biosample_id");

	    return unless $biosample_xml;
	}
	
	my $xml = XMLin($biosample_xml, ForceArray => ["Row"]);

	return unless ref $xml->{BioSample}->{Attributes}->{Attribute} eq 'ARRAY';

	foreach my $attrib (@{$xml->{BioSample}->{Attributes}->{Attribute}}){
	
		my $attrib_name = $attrib->{harmonized_name};
		my $attrib_value = $attrib->{content};

		if ($attrib_name=~/lat_lon/i){ 
			my ($latitude, $longitude) = $attrib=~/(\d+\.*\d* [NS])\s+(\d+\.*\d* [EW])/;
			$genome->{latitude} = $latitude;
			$genome->{longitude} = $longitude;	
		}else{
			my $patric_attrib_name = biosample2patricAttrib($attrib_name);
			next unless ($patric_attrib_name && $attrib_value && not $attrib_value=~/^ *(-|missing|na|n\/a|not available|not provided|not determined|nd|unknown) *$/i);
		
			if ($patric_attrib_name=~/other|additional/i){
				my ($attrib1, $attrib2) = $patric_attrib_name=~/^([^:]*):(.*)/;
				push @{$genome->{$attrib1}}, "$attrib2:$attrib_value";
			}elsif ($patric_attrib_name=~/comments/i){
				push @{$genome->{$patric_attrib_name}}, $attrib_value;
			}else{
				$genome->{$patric_attrib_name} = $attrib_value unless $genome->{$patric_attrib_name}; # get metadata from biosample, only if not available from genbank file
			} 
		}

	}

	$genome->{genome_name} .= " strain $genome->{strain}" if ($genome->{strain} && (not $genome->{genome_name}=~/$genome->{strain}/) ) ;
	
	# parse AMR metadata

	return unless $xml->{BioSample}->{Description}->{Comment}->{Table}->{class}=~/Antibiogram/i;

	foreach my $row (@{$xml->{BioSample}->{Description}->{Comment}->{Table}->{Body}->{Row}}){
	
		my @amr1 = @{$row->{Cell}};

		my $amr;

		$amr->{owner} = $genome->{owner};
		$amr->{public} = $public;
		$amr->{genome_id} = $genome->{genome_id};
		$amr->{genome_name} = $genome->{genome_name};
		$amr->{taxon_id} = $genome->{taxon_id};

		$amr->{antibiotic} = lc $amr1[0]; 	
		$amr->{resistant_phenotype} = ucfirst $amr1[1];	
		$amr->{measurement_sign} = $amr1[2] unless ref $amr1[2] eq ref {};	
		$amr->{measurement_value} = $amr1[3] unless ref $amr1[3] eq ref {};
		$amr->{measurement} = $amr->{measurement_sign}.$amr->{measurement_value}; 	
		$amr->{measurement_unit} = $amr1[4] unless ref $amr1[4] eq ref {}; 	
		$amr->{laboratory_typing_method} = $amr1[5] unless ref $amr1[5] eq ref {}; 	
		$amr->{laboratory_typing_platform} = $amr1[6] unless ref $amr1[6] eq ref {}; 	
		$amr->{vendor} = $amr1[7] unless ref $amr1[7] eq ref {}; 	
		$amr->{laboratory_typing_method_version} = $amr1[8] unless ref $amr1[8] eq ref {}; 	
		$amr->{testing_standard} = $amr1[9] unless ref $amr1[9] eq ref {}; 	

		push @{$genome->{antimicrobial_resistance}}, ucfirst $amr->{resistant_phenotype} unless (grep {$_ eq ucfirst $amr->{resistant_phenotype}} @{$genome->{antimicrobial_resistance}}); 
		$genome->{antimicrobial_resistance_evidence} = "AMR Panel";
	
		push @genome_amr, $amr;

	}

}


sub curateMetadata {

	print "Auto curate genome metadata\n";

	my ($href1, $href2, $href3) = readMetadataRefs();

	my %host_map = %{$href1};
	my %country_map = %{$href2};
	my %country_flu = %{$href3};

	# clean host name
	my $host_name_orig = $genome->{host_name}; 
	$genome->{host_name}=~s/^\s*|\s*$|"//g;
	$genome->{host_name}=~s/ *;.*//g;
	$genome->{host_name}=~s/.*, *//;
	$genome->{host_name} = ucfirst lc $genome->{host_name};
	$genome->{host_name} = "Human" 
		if $genome->{host_name}=~/\b(homosap\W+|homo|human|patient|boy|girl|man|woman)\b/ && not $genome->{host_name}=~/non[-\s]*human/;
	
	# host group and common name
	if ($host_map{$genome->{host_name}}){ 
		($genome->{host_common_name}, $genome->{host_group}) = split /\t/, $host_map{$genome->{host_name}};
		$genome->{host_common_name} = "" if $genome->{host_common_name}=~/null/i; 
		$genome->{host_group} = "" if $genome->{host_group}=~/null/i; 
	}else{
		#$host_common_name = "Unknown";
		#$host_group = "Unknown";
	} 

	# host gender
	if ($host_name_orig=~/\b(male|man|boy|M)\b/i){
		$genome->{host_gender} = "male";
	}elsif($host_name_orig=~/\b(female|woman|girl|F)\b/i){
		$genome->{host_gender} = "female";
	}else{
		# no host gender  
	}
	
	# host age 
	if ($host_name_orig=~/([<>\d]+)[\s-]*(year|y)/i){
		$genome->{host_age} = "$1 years";
	}elsif($host_name_orig=~/([<>\d]+)[\s-]*(month|m)/i){
		$genome->{host_age} = "$1 months";
	}elsif($host_name_orig=~/([<>\d]+)[\s-]*(week|w)/i){
		$genome->{host_age} = "$1 weeks";
	}elsif($host_name_orig=~/([<>\d]+)[\s-]*(day|d)/i){
		$genome->{host_age} = "$1 days";
	}elsif($host_name_orig=~/age[\s]*([<>\d]+)$/i){
		$genome->{host_age} = "$1 years";
	}else{
		# no host age 
	}
		
	# lab host
	if ($genome->{lab_host}){
		#$genome->{host_common_name} = "Lab host";
		#$genome->{host_group} = "Lab";
	}elsif($host_name_orig=~/passage|clone|laboratory|culture|select|cell|mdck/i){
		#$genome->{host_common_name} = "Lab host";
		#$genome->{host_group} = "Lab";
		#$genome->{lab_host} = $host_name_orig;
		#$genome->{host_name} = "";
	}else{
		# No lab host
	}

	# lab reassortment
	if ($genome->{genome_name}=~/recomb|reassort/i){
		$genome->{host_common_name} = "Lab reassortment";
		$genome->{host_group} = "Lab";
		$genome->{lab_host} = $host_name_orig unless $genome->{lab_host};
		$genome->{host_name} = "";
	}
	
	# passage 
	if ($genome->{passage}){
		$genome->{passage}=~s/.*(passage.details|passage.history|passage level|passaged in)([ :]*)//ig;
		$genome->{passage}=~s/.*passage: *//i;
		$genome->{passage}=~s/ *passage$| *passage\(s\)$//ig;
		$genome->{passage}=~s/ *\([0-9\-\/]*\)//;
		$genome->{passage}=~s/;([^;]*)(antigenic|lineage|sensitive|resist)(.*)$//;
		$genome->{passage}=~s/mdck\s*(\d+)/MDCK$1/i;
		$genome->{passage}=~s/egg\s*(\d+)/Egg$1/i;
		$genome->{passage}=~s/egg/Egg/i;
		$genome->{passage}=~s/([A-Z]+)([- ]+)(\d+)/$1$3/g;
	}
	
	# isolation country and geographic group
	if ($genome->{geographic_location}){
		my $country = $genome->{geographic_location};
		$country=~s/"//g;
		$country=~s/\s*[:,;].*//;
		if ($country_map{$country}){
			$genome->{isolation_country} = $country;
			$genome->{geographic_group} = $country_map{$country};
		}else{
			# not on the country list
		}
	} 
	
	# collection year
	if ($genome->{collection_date}=~/(\d\d\d\d).*(\d\d\d\d)/){
		$genome->{collection_year} = $1 if $1 eq $2;
	}elsif($genome->{collection_date}=~/(\d\d\d\d)/){
		$genome->{collection_year} = $1;
	}else{
		# no collection year
	}
	
	# flu season
	my ($month, $year, $y1, $y2);
	if ($genome->{genome_name}=~/Influenza/i && $country_flu{$genome->{isolation_country}}){	
		if ($genome->{collection_date}=~/(\d\d)-([A-Za-z]+)-\d\d(\d\d)/){
			$month = $2;
			$year = $3;	
		}elsif($genome->{collection_date}=~/\d\d(\d\d)-(\d\d)-(\d\d)/){
			$month = $2;
			$year = $1;
		}else{
			# unexpected date format
		}

		if ($month=~/oct|nov|dec|10|11|12/i){
			$y1 = $year;
			$y2 = $y1==99? "00": $y1+1;	
		}elsif($month=~/jan|feb|mar|apr|may|01|02|03|04|05/i){
			$y1 = $year=="00"? 99: $year-1;
			$y2 = $y1==99? "00": $y1+1;	
		}else{
		}		

		$genome->{season} = sprintf("%02d", $y1)."-".sprintf("%02d", $y2) if $y1 && $y2;			
	}
	
	# serovar / serotype 
	if ($genome->{genome_name}=~/Influenza/i){
		if ($genome->{serovar}){
			$genome->{subtype} = $genome->{serovar};
		}elsif($genome->{genome_name}=~/(H\d+N\d+)/){
			$genome->{subtype} = $1;
		}else{
		}	
		$genome->{subtype}=~s/\(untyped\)|\?/x/g;
  	$genome->{subtype}=~s/"|\(|\)//g;
  	$genome->{subtype}= "" if $genome->{subtype}=~/^(unknown\|unidentified)$/i;
  	$genome->{subtype} = "Mixed" if $genome->{subtype}==~/mixed/i;
		($genome->{h_type}, $genome->{n_type}) = $genome->{subtype}=~/H0*([\d]+)N0*([\d]+)/;
	}

	# Segments labels for segmented viruses
	$genome->{segment}=~s/^(segment|rna)\s*//i;
	if (grep {$_=~/Influenza/} @{$genome->{taxon_lineage_names}}){		
		#
	}elsif(grep {$_=~/Bunyavirales/} @{$genome->{taxon_lineage_names}}){
		#	
	}elsif(grep {$_=~/Reoviridae/} @{$genome->{taxon_lineage_names}}){
		#
	}else{

	}

}


sub biosample2patricAttrib{

	my ($attribute) = @_;

	my %biosample2patric = (

		"altitude" => "altitude",
		"biomaterial_provider" => "additional_metadata:biomaterial_provider",
		"collected_by" => "additional_metadata:collected_by",
		"collection_date" => "collection_date",			
		"culture_collection" => "culture_collection", 
		"depth" => "depth", 
		"description" => "comments", 
		"env_biome" => "other_environmental:env_biome", 
		"genotype" => "other_typing:genotype", 
		"geo_loc_name" => "geographic_location", 
		"host" => "host_name", 
		"host_age" => "host_age", 
		"host_description" => "other_clinical:host_description", 
		"host_disease" => "host_health",
		"host_disease_outcome" => "other_clinical:host_disease_outcome",
		"host_disease_stage" => "other_clinical:host_disease_stage",
		"host_health_state" => "other_clinical:host_health_state",
		"host_sex" => "host_gender",
		"host_subject_id" => "other_clinical:host_subject_id", 
		"host_tissue_sampled" => "isolation_source",
		"identified_by" => "additional_metadata:identified_by",
		"isolation_source" => "isolation_source", 
		"lab_host" => "lab_host", 
		"lat_lon" => "other_environmental:lat_lon", 
		"mating_type" => "additional_metadata:mating_type", 
		"organism" => "",
		"passage_history" => "passage",
		"passage" => "passage",
		"pathotype" => "pathovar",
		"sample_name" => "",
		"sample_title" => "",
		"sample_type" => "additional_metadata:sample_type",
		"samp_size" => "",
		"serotype" => "serovar",
		"serovar" => "serovar",
		"specimen_voucher" => "additional_metadata:specimen_voucher", 
		"strain" => "strain",
		"isolate" => "strain",
		"subgroup" => "",
		"subtype" => "",
		"temp" => "other_environmental:temperature"
	);

	return $biosample2patric{$attribute} if $attribute;	

}


sub prepareTaxonomy {
	
	"Preparing taxonomy update ...\n";

	my ($taxon_lineage_ids) = @_;

	foreach my $taxon_id (@{$taxon_lineage_ids}){
	
		my $taxon;	
		$taxon->{taxon_id} = $taxon_id;
		$taxon->{genomes}->{inc} = 1;		

		push @taxonomy, $taxon

	}

}

sub get_xml
{
    my($url) = @_;
    my $res = $user_agent->get($url);
    if (!$res->is_success)
    {
	warn "Failure retrieving $url: " . $res->content;
	return;
    }
    my $xml = $res->content;
    return $xml;
}

sub get_xml_to_file
{
    my($file, $url) = @_;
    my $xml = get_xml($url);
    write_file($file, $xml);
}

sub readMetadataRefs {

    my $refs_dir = "$ENV{KB_TOP}/lib/autocuration-metadata";
    $refs_dir = "$Bin/refs" unless -d $refs_dir;

    my(%host_map, %country_map, %country_flu);

    # process host mappings
    if (open FH, "$refs_dir/host_mapping")
    {
	%host_map = ();
	while (my $entry = <FH>){
	    chomp $entry;
	    my ($host_name, $host_common_name, $host_group) = $entry=~/(.*)\t(.*)\t(.*)/;
	    $host_map{$host_name} = "$host_common_name\t$host_group";
			$host_map{lc $host_name} = "$host_common_name\t$host_group";
			$host_map{ucfirst lc $host_name} = "$host_common_name\t$host_group";
	}
	close FH;
    }
    else
    {
	warn "Cannot open host mapping file $refs_dir/host_mapping: $!";
    }

    # process country mapping for flu season
    if (open FH, "$refs_dir/country_mapping")
    {
	%country_map = ();
	%country_flu = ();
	while (my $entry = <FH>) {
	    
	    chomp $entry;
	    my ($country, $geographic_group, $flu_season) = $entry=~/(.*)\t(.*)\t(.*)/;
	    $country_map{$country} = "$geographic_group";
	    $country_flu{$country} = 1 if $flu_season=~/yes/i;
	}
	close FH;
    }
    else
    {
	warn "Cannot open country mapping file $refs_dir/country_mapping: $!";
    }
    return (\%host_map, \%country_map, \%country_flu);
}
