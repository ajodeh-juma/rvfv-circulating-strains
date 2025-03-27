////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

params.summary_params = [:]

// Validate input parameters

// Check input path parameters to see if they exist
checkPathParamList = [
    params.fasta,
    params.metadata,
    params.outliers,
    params.lineages,
    params.recombinants
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


// check given segment
def segmentsList = ['S-NSS', 'S-NP', 'M', 'L']
if (!segmentsList.contains(params.segment)) {
    exit 1, "Invalid segment option: ${params.segment}. Valid options: ${segmentsList.join(', ')}"
}

def referencesList = ["DQ375404", "DQ375417", "DQ375430", "DQ380208", "DQ380213", "DQ380193", "DQ380154", "DQ380182", "DQ380157"]
if (!referencesList.contains(params.vaccine_reference)) {
    exit 1, "Invalid reference option: ${params.vaccine_reference}. Valid options: ${referencesList.join(', ')}"
}

// def filterColumnsList = ['country', 'location', 'date', 'host']
// if (!filterColumnsList.contains(params.filter_column)) {
//     exit 1, "Invalid column option: ${params.filter_column}. Valid options: ${filterColumnsList.join(', ')}"
// }

// def groupingColumnsList = ['lineage', 'host', 'country', 'year', 'strain', 'strain_type']
// if (!groupingColumnsList.contains(params.grouping_column)) {
//     exit 1, "Invalid grouping column option: ${params.grouping_column}. Valid options: ${groupingColumnsList.join(', ')}"
// }

// def seqTypeColumnsList = ['dna', 'protein']
// if (!seqTypeColumnsList.contains(params.seq_type)) {
//     exit 1, "Invalid sequence type option: ${params.seq_type}. Valid options: ${seqTypeColumnsList.join(', ')}"
// }

// def snpTypeColumnsList = ['singleton', 'multiple', 'conserved', 'all']
// if (!snpTypeColumnsList.contains(params.snp_type)) {
//     exit 1, "Invalid snp type option: ${params.snp_type}. Valid options: ${snpTypeColumnsList.join(', ')}"
// }



ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
ch_metadata = Channel.fromPath(params.metadata, checkIfExists: true)
ch_outliers = Channel.fromPath(params.outliers, checkIfExists: true)
ch_lineages = Channel.fromPath(params.lineages, checkIfExists: true)
ch_vaccine_reference = Channel.of(params.vaccine_reference)
// ch_filter_column = Channel.of(params.filter_column)
ch_grouping_column = Channel.of(params.grouping_column)
ch_recombinants = Channel.fromPath(params.recombinants, checkIfExists: true)


ch_prefix = Channel.of(params.prefix)
                    .map { row ->
                        def meta = [:]
                        def id = params.prefix
                        meta.id = "$id"
                        [ meta ,  params.prefix ]
                        }


////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()


////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

include { GET_SOFTWARE_VERSIONS } from '../modules/local/process/get_software_versions' addParams( options: [publish_files : ['csv':'']]                 )
include { SUBSET_BY_COUNTRY }     from '../modules/local/process/subset_by_country' addParams( options: modules['subset_by_country']                            )
include { REFORMAT_HEADERS }     from '../modules/local/process/reformat_headers' addParams( options: modules['reformat_headers']                            )
include { MAFFT_ALIGN }     from '../modules/local/process/mafft_align' addParams( options: modules['alignment']                            )
include { SUBSET_COORDINATES }     from '../modules/local/process/subset_coordinates' addParams( options: modules['subset_coordinates']                            )
include { SEQKIT_RMDUP }     from '../modules/local/process/seqkit_rmdup' addParams( options: modules['deduplicate']                            )
include { GET_DUPLICATED_TAXA }     from '../modules/local/process/get_duplicated_taxa' addParams( options: modules['get_duplicated_taxa']                            )
include { FILTER_TAXA }     from '../modules/local/process/filter_taxa' addParams( options: modules['filter_taxa']                            )
include { MODEL_TEST }     from '../modules/local/process/model_test' addParams( options: modules['model_test']                            )
include { IQTREE_TREE }     from '../modules/local/process/iqtree_tree' addParams( options: modules['tree']                            )
include { PREPARE_GEOLOCATIONS }     from '../modules/local/process/prepare_geolocations' addParams( options: modules['prepare_geolocations']                            )
include { TREETIME_TREE }     from '../modules/local/process/treetime_tree' addParams( options: modules['treetime']                            )
include { GEOCODE_LOCATIONS }     from '../modules/local/process/geocode_locations' addParams( options: modules['prepare_geolocations']                            )
// include { HYPHY_DEDUP }     from '../modules/local/process/hyphy_dedup' addParams( options: modules['hyphy_dedup']                            )
include { HYPHY_FUBAR }     from '../modules/local/process/hyphy_fubar' addParams( options: modules['hyphy']                            )
include { HYPHY_FEL }     from '../modules/local/process/hyphy_fel' addParams( options: modules['hyphy']                            )
include { HYPHY_SLAC }     from '../modules/local/process/hyphy_slac' addParams( options: modules['hyphy']                            )

include { STRAIN_TYPES }     from '../modules/local/process/strain_types' addParams( options: modules['strain_types']                            )
include { EXTRACT_ACCESSIONS }     from '../modules/local/process/extract_accessions' addParams( options: modules['extract_accessions']                            )
include { SNPSITES }     from '../modules/local/process/snpsites' addParams( options: modules['snp_sites']                            )
include { CONVERT_DNA_TO_AMINOACID }     from '../modules/local/process/convert_dna_to_aminoacid' addParams( options: modules['convert_dna_to_aminoacid']                            )
include { GET_SNPS_RECORDS }     from '../modules/local/process/get_snps_records' addParams( options: modules['get_snps_records']                            )
include { EXTRACT_ABUNDANT_SNPS }     from '../modules/local/process/extract_abundant_snps' addParams( options: modules['extract_abundant_snps']                            )
include { SNIPIT }     from '../modules/local/process/snipit' addParams( options: modules['snipit']                            )
include { COMPUTE_TSTV_STATS }     from '../modules/local/process/compute_tstv_stats' addParams( options: modules['compute_tstv_stats']                            )
include { SUMMARIZE_SNPS }     from '../modules/local/process/summarize_snps' addParams( options: modules['summarize_snps']                            )
include { SUBSET_BY_HOST }     from '../modules/local/process/subset_by_host' addParams( options: modules['subset_by_host']                            )
include { REVERSE_COMPLEMENT }     from '../modules/local/process/reverse_complement' addParams( options: modules['reverse_complement']                            )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////


workflow RVFVCIRCULATINGSTRAINS {

    // MODULE: subset by country
    SUBSET_BY_COUNTRY(ch_prefix, ch_fasta, ch_metadata)
    ch_country_fasta = SUBSET_BY_COUNTRY.out.fasta
    ch_country_csv = SUBSET_BY_COUNTRY.out.csv

    // MODULE: reformat sequence headers
    ch_input_reformat_headers = ch_country_fasta.join(ch_country_csv)
    ch_input_reformat_headers
        .map { row -> 
        def meta = [:]
        meta.id = row[0].id
        [ meta , [ file(row[1], checkIfExists: true), file(row[2], checkIfExists: true) ] ]
    }.set { ch_input_reformat_headers }

    REFORMAT_HEADERS(ch_input_reformat_headers)

    // MODULE: multiple sequence alignment
    MAFFT_ALIGN(REFORMAT_HEADERS.out.fasta)

    // MODULE: subset coordinates
    SUBSET_COORDINATES(MAFFT_ALIGN.out.alignment)

    // MODULE: deduplicate sequences by composition
    SEQKIT_RMDUP(SUBSET_COORDINATES.out.fasta)
    GET_DUPLICATED_TAXA(SEQKIT_RMDUP.out.txt, ch_outliers)

    // MODULE: filter out duplicated and outlier taxa
    FILTER_TAXA(SUBSET_COORDINATES.out.fasta,
                REFORMAT_HEADERS.out.txt,
                GET_DUPLICATED_TAXA.out.dup_taxa
                )
    
    if ( params.segment == 'S-NP' ) {
        REVERSE_COMPLEMENT(ch_prefix, FILTER_TAXA.out.fasta)
        ch_filtered_fasta = REVERSE_COMPLEMENT.out.fasta
    } else {
        ch_filtered_fasta = FILTER_TAXA.out.fasta
    }

    // MODULE: test evolutionary models
    MODEL_TEST(ch_filtered_fasta)

    // MODULE: maximum likelihood tree
    IQTREE_TREE(ch_filtered_fasta, MODEL_TEST.out.log)

    // MODULE: prepare metadata for molecular clock analysis 
    PREPARE_GEOLOCATIONS(FILTER_TAXA.out.txt,
                        SUBSET_BY_COUNTRY.out.csv
                        )
    
    // MODULE: treetime 
    TREETIME_TREE(IQTREE_TREE.out.treefile,
                    ch_filtered_fasta,
                    PREPARE_GEOLOCATIONS.out.dates
                    )
    // MODULE: geocode locations
    GEOCODE_LOCATIONS(PREPARE_GEOLOCATIONS.out.geolocations)

    // MODULE: add strain types
    STRAIN_TYPES(ch_filtered_fasta, ch_metadata, ch_lineages)
    EXTRACT_ACCESSIONS(ch_filtered_fasta, ch_vaccine_reference, STRAIN_TYPES.out.txt, ch_grouping_column)
    SNPSITES(EXTRACT_ACCESSIONS.out.fasta, ch_grouping_column, ch_vaccine_reference)
    CONVERT_DNA_TO_AMINOACID(ch_filtered_fasta, ch_grouping_column, ch_vaccine_reference)
    GET_SNPS_RECORDS(CONVERT_DNA_TO_AMINOACID.out.amino_acid, STRAIN_TYPES.out.txt, ch_vaccine_reference, ch_grouping_column, ch_lineages)
    EXTRACT_ABUNDANT_SNPS(SNPSITES.out.vcf, ch_vaccine_reference)
    SNIPIT(EXTRACT_ACCESSIONS.out.fasta, ch_vaccine_reference, EXTRACT_ACCESSIONS.out.csv, EXTRACT_ABUNDANT_SNPS.out.txt)
    COMPUTE_TSTV_STATS(SNPSITES.out.vcf, STRAIN_TYPES.out.txt, ch_vaccine_reference, ch_recombinants, ch_grouping_column)
    SUMMARIZE_SNPS(SNIPIT.out.csv, STRAIN_TYPES.out.txt, ch_recombinants, ch_vaccine_reference)

    SUBSET_BY_HOST(ch_prefix, CONVERT_DNA_TO_AMINOACID.out.amino_acid, FILTER_TAXA.out.txt)


    // // MODULE: selection analysis
    // // HYPHY_DEDUP(FILTER_TAXA.out.fasta, IQTREE_TREE.out.treefile)
    HYPHY_FUBAR(ch_filtered_fasta, IQTREE_TREE.out.treefile)
    HYPHY_FEL(ch_filtered_fasta, IQTREE_TREE.out.treefile)
    HYPHY_SLAC(ch_filtered_fasta, IQTREE_TREE.out.treefile)

    


    // channel to collect software versions
    ch_software_versions = Channel.empty()

    /*
     * MODULE: Pipeline reporting
     */
    
    GET_SOFTWARE_VERSIONS ( 
        ch_software_versions.map { it }.collect()
    )
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    // Completion.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    Completion.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/