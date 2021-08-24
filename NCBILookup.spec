module NCBIAnnotation {

    funcdef get_metadata_from_bioproject(string accn) returns (string bioproject_id,
							       string bioproject_xml);
    funcdef get_metadata_from_biosample(string accn) returns (string biosample_id,
							       string biosample_xml);

    funcdef get_assembly_accession(string accn) returns (string assembly_accession);

};
