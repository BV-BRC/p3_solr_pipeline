package Bio::P3::NCBILookup::NCBILookupImpl;
use strict;
# Use Semantic Versioning (2.0.0-rc.1)
# http://semver.org 
our $VERSION = "0.1.0";

=head1 NAME

NCBIAnnotation

=head1 DESCRIPTION



=cut

#BEGIN_HEADER

use SolrAPI;
use P3RateLimitedUserAgent;
use Data::Dumper;
use JSON::XS;

our $have_config;

eval {
    require Bio::KBase::AppService::AppConfig;
    $have_config = 1;
};
use XML::Simple;

sub get_xml
{
    my($self, $url) = @_;
    print "retr: $url\n";
    my $res = $self->{user_agent}->get($url);
    if (!$res->is_success)
    {
	die "Failure retrieving $url: " . $res->content;
    }
    my $xml = $res->content;
    return $xml;
}

#END_HEADER

sub new
{
    my($class, @args) = @_;
    my $self = {
    };
    bless $self, $class;
    #BEGIN_CONSTRUCTOR

    my $limit = 3;

    if (my $k = $ENV{NCBI_API_KEY})
    {
	$self->{api_key} = $k;
	$limit = 7;
    }
    $self->{user_agent} = P3RateLimitedUserAgent->new($limit);

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

    $self->{data_api_url} = $data_api_url;
    $self->{reference_data_dir} = $reference_data_dir;
    
    $self->{json} = JSON->new->allow_nonref;


    #END_CONSTRUCTOR

    if ($self->can('_init_instance'))
    {
	$self->_init_instance();
    }
    return $self;
}
=head1 METHODS
=head2 get_metadata_from_bioproject

  $bioproject_id, $bioproject_xml = $obj->get_metadata_from_bioproject($accn)

=over 4


=item Parameter and return types

=begin html

<pre>
$accn is a string
$bioproject_id is a string
$bioproject_xml is a string
</pre>

=end html

=begin text

$accn is a string
$bioproject_id is a string
$bioproject_xml is a string

=end text



=item Description


=back

=cut

sub get_metadata_from_bioproject
{
    my $self = shift;
    my($accn) = @_;

    my @_bad_arguments;
    (!ref($accn)) or push(@_bad_arguments, "Invalid type for argument \"accn\" (value was \"$accn\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to get_metadata_from_bioproject:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::P3::NCBILookup::Service::CallContext;
    my($bioproject_id, $bioproject_xml);
    #BEGIN get_metadata_from_bioproject

    if (exists $self->{bioproject_cache}->{$accn})
    {
	($bioproject_id, $bioproject_xml) = @{ $self->{bioproject_cache}->{$accn}};
	print "Cache hit $accn\n";
    }
    else
    {
	my $url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=bioproject&term=$accn";
	$url .= "&api_key=$self->{api_key}" if $self->{api_key};
	my $xml = $self->get_xml($url);
	$xml =~ s/\n//;
	
	
	($bioproject_id) = $xml=~/<Id>(\d+)<\/Id>/;
	
	$bioproject_xml = $self->get_xml("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=bioproject&retmode=xml&id=$bioproject_id");
	$self->{bioproject_cache}->{$accn} = [$bioproject_id, $bioproject_xml];
    }

    #END get_metadata_from_bioproject
    my @_bad_returns;
    (!ref($bioproject_id)) or push(@_bad_returns, "Invalid type for return variable \"bioproject_id\" (value was \"$bioproject_id\")");
    (!ref($bioproject_xml)) or push(@_bad_returns, "Invalid type for return variable \"bioproject_xml\" (value was \"$bioproject_xml\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to get_metadata_from_bioproject:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($bioproject_id, $bioproject_xml);
}


=head2 get_metadata_from_biosample

  $biosample_id, $biosample_xml = $obj->get_metadata_from_biosample($accn)

=over 4


=item Parameter and return types

=begin html

<pre>
$accn is a string
$biosample_id is a string
$biosample_xml is a string
</pre>

=end html

=begin text

$accn is a string
$biosample_id is a string
$biosample_xml is a string

=end text



=item Description


=back

=cut

sub get_metadata_from_biosample
{
    my $self = shift;
    my($accn) = @_;

    my @_bad_arguments;
    (!ref($accn)) or push(@_bad_arguments, "Invalid type for argument \"accn\" (value was \"$accn\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to get_metadata_from_biosample:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::P3::NCBILookup::Service::CallContext;
    my($biosample_id, $biosample_xml);
    #BEGIN get_metadata_from_biosample

    if (exists $self->{biosample_cache}->{$accn})
    {
	($biosample_id, $biosample_xml) = @{$self->{biosample_cache}->{$accn}};
	print "Got sample cache $accn\n";
    }
    else
    {
	print "Getting genome metadata from BioSample: $accn ...\n";
	my $url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=biosample&term=$accn";  
	$url .= "&api_key=$self->{api_key}" if $self->{api_key};
	my $xml = $self->get_xml($url);
	$xml =~ s/\n//;
	($biosample_id) = $xml=~/<Id>(\d+)<\/Id>/;
	
	$biosample_xml = $self->get_xml("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&retmode=xml&id=$biosample_id");
	$self->{biosample_cache}->{$accn} = [$biosample_id, $biosample_xml];
    }


    #END get_metadata_from_biosample
    my @_bad_returns;
    (!ref($biosample_id)) or push(@_bad_returns, "Invalid type for return variable \"biosample_id\" (value was \"$biosample_id\")");
    (!ref($biosample_xml)) or push(@_bad_returns, "Invalid type for return variable \"biosample_xml\" (value was \"$biosample_xml\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to get_metadata_from_biosample:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($biosample_id, $biosample_xml);
}


=head2 get_assembly_accession

  $assembly_accession = $obj->get_assembly_accession($accn)

=over 4


=item Parameter and return types

=begin html

<pre>
$accn is a string
$assembly_accession is a string
</pre>

=end html

=begin text

$accn is a string
$assembly_accession is a string

=end text



=item Description


=back

=cut

sub get_assembly_accession
{
    my $self = shift;
    my($accn) = @_;

    my @_bad_arguments;
    (!ref($accn)) or push(@_bad_arguments, "Invalid type for argument \"accn\" (value was \"$accn\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to get_assembly_accession:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::P3::NCBILookup::Service::CallContext;
    my($assembly_accession);
    #BEGIN get_assembly_accession

    if (exists $self->{assembly_cache}->{$accn})
    {
	$assembly_accession = $self->{assembly_cache}->{$accn};
    }
    else
    {
	my $url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=$accn";
	
	$url .= "&api_key=$self->{api_key}" if $self->{api_key};
	my $xml = $self->get_xml($url);
	
	$xml =~ s/\n//;
	my ($gi) = $xml=~/<Id>(\d+)<\/Id>/;
	
	my $xml = $self->get_xml("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=assembly&id=$gi");
	$xml =~ s/\n//;
	my($assembly_id) = $xml=~/<Link>\s*<Id>(\d+)<\/Id>/;
	
	my $xml = $self->get_xml("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id=$assembly_id");
	
	($assembly_accession) = $xml=~/<Genbank>(\S*)<\/Genbank>/;
	$self->{assembly_cache}->{$accn} = $assembly_accession;
    }
    #END get_assembly_accession
    my @_bad_returns;
    (!ref($assembly_accession)) or push(@_bad_returns, "Invalid type for return variable \"assembly_accession\" (value was \"$assembly_accession\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to get_assembly_accession:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($assembly_accession);
}





=head2 version 

  $return = $obj->version()

=over 4

=item Parameter and return types

=begin html

<pre>
$return is a string
</pre>

=end html

=begin text

$return is a string

=end text

=item Description

Return the module version. This is a Semantic Versioning number.

=back

=cut

sub version {
    return $VERSION;
}



=head1 TYPES


=cut

1;
