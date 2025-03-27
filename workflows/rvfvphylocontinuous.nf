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
    params.lineages
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


// check given segment
def segmentsList = ['S-NSS', 'S-NP', 'M', 'L']
if (!segmentsList.contains(params.segment)) {
    exit 1, "Invalid segment option: ${params.segment}. Valid options: ${segmentsList.join(', ')}"
}

//def filterColumnsList = ['country', 'location', 'date', 'host']
//if (!filterColumnsList.contains(params.filter_columns)) {
//    exit 1, "Invalid column option: ${params.filter_columns}. Valid options: ${filterColumnsList.join(', ')}"
//}

ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
ch_metadata = Channel.fromPath(params.metadata, checkIfExists: true)
ch_outliers = Channel.fromPath(params.outliers, checkIfExists: true)
ch_lineages = Channel.fromPath(params.lineages, checkIfExists: true)
ch_filter_columns = Channel.of(params.filter_columns)



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
include { FILTER_ON_COLUMN } from '../modules/local/process/filter_on_column' addParams(options: modules['filter_on_column'])

// include { SUBSET_BY_COUNTRY }     from '../modules/local/process/subset_by_country' addParams( options: modules['subset_by_country']                            )
include { REFORMAT_HEADERS }     from '../modules/local/process/reformat_headers' addParams( options: modules['reformat_headers']                            )
include { SEQKIT_RMDUP }     from '../modules/local/process/seqkit_rmdup' addParams( options: modules['deduplicate']                            )
include { GET_DUPLICATED_TAXA }     from '../modules/local/process/get_duplicated_taxa' addParams( options: modules['get_duplicated_taxa']                            )
include { FILTER_TAXA }     from '../modules/local/process/filter_taxa' addParams( options: modules['alignment']                            )
include { MAFFT_ALIGN }     from '../modules/local/process/mafft_align' addParams( options: modules['alignment']                            )
include { SUBSET_COORDINATES }     from '../modules/local/process/subset_coordinates' addParams( options: modules['subset_coordinates']                            )
include { MODEL_TEST }     from '../modules/local/process/model_test' addParams( options: modules['model_test']                            )
include { IQTREE_TREE }     from '../modules/local/process/iqtree_tree' addParams( options: modules['tree']                            )
include { IQTREE_TREE_NO_MODEL }     from '../modules/local/process/iqtree_nomodel' addParams( options: modules['tree']                            )
include { PREPARE_GEOLOCATIONS }     from '../modules/local/process/prepare_geolocations' addParams( options: modules['prepare_geolocations']                            )
include { TREETIME_TREE }     from '../modules/local/process/treetime_tree' addParams( options: modules['treetime']                            )
include { GEOCODE_LOCATIONS }     from '../modules/local/process/geocode_locations' addParams( options: modules['prepare_geolocations']                            )
include { REVERSE_COMPLEMENT }     from '../modules/local/process/reverse_complement' addParams( options: modules['reverse_complement']                            )


////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////


workflow RVFVPHYLOCONTINUOUS {

    // MODULE: filter on specified columns
    FILTER_ON_COLUMN(ch_prefix, ch_fasta, ch_metadata, ch_filter_columns)

    

    // // MODULE: subset by country
    // SUBSET_BY_COUNTRY(ch_prefix, ch_fasta, ch_metadata)
    // ch_country_fasta = SUBSET_BY_COUNTRY.out.fasta
    // ch_country_csv = SUBSET_BY_COUNTRY.out.csv

    // MODULE: reformat sequence headers
    ch_input_reformat_headers = FILTER_ON_COLUMN.out.fasta.join(FILTER_ON_COLUMN.out.csv)
    ch_input_reformat_headers
        .map { row -> 
        def meta = [:]
        meta.id = row[0].id
        [ meta , [ file(row[1], checkIfExists: true), file(row[2], checkIfExists: true) ] ]
    }.set { ch_input_reformat_headers }

    REFORMAT_HEADERS(ch_input_reformat_headers)

    // MODULE: deduplicate sequences by composition
    SEQKIT_RMDUP(REFORMAT_HEADERS.out.fasta)
    GET_DUPLICATED_TAXA(SEQKIT_RMDUP.out.txt, ch_outliers)

    // MODULE: filter out duplicated and outlier taxa
    FILTER_TAXA(REFORMAT_HEADERS.out.fasta,
                REFORMAT_HEADERS.out.txt,
                GET_DUPLICATED_TAXA.out.dup_taxa
                )

    // MODULE: multiple sequence alignment
    MAFFT_ALIGN(FILTER_TAXA.out.fasta)

    // MODULE: subset coordinates
    SUBSET_COORDINATES(MAFFT_ALIGN.out.alignment)


    if ( params.segment == 'S-NP' ) {
        REVERSE_COMPLEMENT(ch_prefix, SUBSET_COORDINATES.out.fasta)
        ch_filtered_fasta = REVERSE_COMPLEMENT.out.fasta
    } else {
        ch_filtered_fasta = SUBSET_COORDINATES.out.fasta
    }


    if (params.skip_modeltesting) {
        IQTREE_TREE_NO_MODEL(SUBSET_COORDINATES.out.fasta)
        ch_treefile = IQTREE_TREE_NO_MODEL.out.treefile

    } else {
        // MODULE: test evolutionary models
        // MODEL_TEST(SUBSET_COORDINATES.out.fasta)
        MODEL_TEST(ch_filtered_fasta)
        
        // MODULE: maximum likelihood tree
        IQTREE_TREE(ch_filtered_fasta, MODEL_TEST.out.log)
        // IQTREE_TREE(SUBSET_COORDINATES.out.fasta, MODEL_TEST.out.log)
        ch_treefile = IQTREE_TREE.out.treefile
    }

    // MODULE: prepare metadata for molecular clock analysis 
    PREPARE_GEOLOCATIONS(FILTER_TAXA.out.txt,
                        FILTER_ON_COLUMN.out.csv
                        )
    
    // MODULE: treetime 
    TREETIME_TREE(ch_treefile,
                    ch_filtered_fasta,
                    PREPARE_GEOLOCATIONS.out.dates
                    )
    // SUBSET_COORDINATES.out.fasta

    // MODULE: geocode locations
    GEOCODE_LOCATIONS(PREPARE_GEOLOCATIONS.out.geolocations)


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