/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPROCESS_REPERTOIRE   } from '../modules/local/preprocess_repertoire/main'
include { SIMULATE_CONVERGENCE    } from '../modules/local/simulation/main'
include { SELECT_SIMULATED_SEQUENCES } from '../modules/local/select_simulated_sequences/main'

// nf-core subworkflows
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_stereotyper_pipeline'


// nf-core modules
include { AMULETY_ANTIBERTA2 } from '../modules/nf-core/amulety/antiberta2/main'
include { AMULETY_ANTIBERTY } from '../modules/nf-core/amulety/antiberty/main'
include { AMULETY_BALMPAIRED } from '../modules/nf-core/amulety/balmpaired/main'
include { AMULETY_ESM2 } from '../modules/nf-core/amulety/esm2/main'
include { AMULETY_ANTIBERTA2 as AMULETY_ANTIBERTA2_SIM } from '../modules/nf-core/amulety/antiberta2/main'
include { AMULETY_ANTIBERTY as AMULETY_ANTIBERTY_SIM } from '../modules/nf-core/amulety/antiberty/main'
include { AMULETY_BALMPAIRED as AMULETY_BALMPAIRED_SIM } from '../modules/nf-core/amulety/balmpaired/main'
include { AMULETY_ESM2 as AMULETY_ESM2_SIM } from '../modules/nf-core/amulety/esm2/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow STEREOTYPER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Process simulation parameters
    //

    def fuzziness = params.fuzziness_param.toString().split(',').collect { it.trim().toInteger() }
    ch_fuzziness = Channel.from(fuzziness)
                         .map { it -> it as Float } // Convert to Float for consistency
                         .dump(tag: 'fuzziness') // Debugging
    def abundance = params.clonal_abundance.toString().split(',').collect { it.trim().toFloat() }
    ch_abundance = Channel.from(abundance)
                            .map { it -> it as Float } // Convert to Float for consistency
                            .dump(tag: 'abundance') // Debugging
    def repertoire_sample = params.subsample_size.toString().split(',').collect { it.trim().toInteger() }
    ch_repertoire_sample = Channel.from(repertoire_sample)
                            .map { it -> it as Integer } // Convert to Integer for consistency
                            .dump(tag: 'repertoire_sample') // Debugging
    def witness = params.witness.toString().split(',').collect { it.trim() }
    ch_witness = Channel.from(witness)
                        .dump(tag: 'witness') // Debugging


    //
    // MODULE: Preprocess repertoire
    //
    PREPROCESS_REPERTOIRE (
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(PREPROCESS_REPERTOIRE.out.versions.first())

    ch_repertoire_embeddings = Channel.empty()

    // Get repertoire embeddings
    if (params.embeddings && params.embeddings.split(',').contains('antiberty') ){
        AMULETY_ANTIBERTY(
            PREPROCESS_REPERTOIRE.out.repertoire,
            params.embedding_chain
        )
        ch_versions = ch_versions.mix(AMULETY_ANTIBERTY.out.versions.first())
        ch_repertoire_embeddings = ch_repertoire_embeddings
                        .mix(AMULETY_ANTIBERTY.out.embedding
                        .map{ it -> ["antiberty", it[0].id, it[0], it[1]] })
    }

    if (params.embeddings && params.embeddings.split(',').contains('antiberta2') ){
        AMULETY_ANTIBERTA2(
            PREPROCESS_REPERTOIRE.out.repertoire,
            params.embedding_chain
        )
        ch_versions = ch_versions.mix(AMULETY_ANTIBERTA2.out.versions.first())
        ch_repertoire_embeddings = ch_repertoire_embeddings
                        .mix(AMULETY_ANTIBERTA2.out.embedding
                        .map{ it -> ["antiberta2", it[0].id, it[0], it[1]] })
    }

    if (params.embeddings && params.embeddings.split(',').contains('esm2') ){
        AMULETY_ESM2(
            PREPROCESS_REPERTOIRE.out.repertoire,
            params.embedding_chain
        )
        ch_versions = ch_versions.mix(AMULETY_ESM2.out.versions.first())
        ch_repertoire_embeddings = ch_repertoire_embeddings
                        .mix(AMULETY_ESM2.out.embedding
                        .map{ it -> ["esm2", it[0].id, it[0], it[1]] })
    }

    if (params.embeddings && params.embeddings.split(',').contains('balmpaired') ){
        AMULETY_BALMPAIRED(
            PREPROCESS_REPERTOIRE.out.repertoire,
            params.embedding_chain
        )
        ch_versions = ch_versions.mix(AMULETY_BALMPAIRED.out.versions.first())
        ch_repertoire_embeddings = ch_repertoire_embeddings
                        .mix(AMULETY_BALMPAIRED.out.embedding
                        .map{ it -> ["balmpaired", it[0].id, it[0], it[1]] })
    }

    // Simulate convergence
    SIMULATE_CONVERGENCE (
        PREPROCESS_REPERTOIRE.out.repertoire
    )
    ch_versions = ch_versions.mix(SIMULATE_CONVERGENCE.out.versions.first())

    ch_sim_embeddings = Channel.empty()
    // Get simulated sequences embeddings
    if (params.embeddings && params.embeddings.split(',').contains('antiberty') ){
        AMULETY_ANTIBERTY_SIM(
            SIMULATE_CONVERGENCE.out.sequences,
            params.embedding_chain
        )
        ch_sim_embeddings = ch_sim_embeddings
                        .mix(AMULETY_ANTIBERTY_SIM.out.embedding
                        .map{ it -> ["antiberty", it[0].id, it[0], it[1]] })
    }

    if (params.embeddings && params.embeddings.split(',').contains('antiberta2') ){
        AMULETY_ANTIBERTA2_SIM(
            SIMULATE_CONVERGENCE.out.sequences,
            params.embedding_chain
        )
        ch_sim_embeddings = ch_sim_embeddings
                        .mix(AMULETY_ANTIBERTA2_SIM.out.embedding
                        .map{ it -> ["antiberta2", it[0].id, it[0], it[1]] })
    }

    if (params.embeddings && params.embeddings.split(',').contains('esm2') ){
        AMULETY_ESM2_SIM(
            SIMULATE_CONVERGENCE.out.sequences,
            params.embedding_chain
        )
        ch_sim_embeddings = ch_sim_embeddings
                        .mix(AMULETY_ESM2_SIM.out.embedding
                        .map{ it -> ["esm2", it[0].id, it[0], it[1]] })
    }

    if (params.embeddings && params.embeddings.split(',').contains('balmpaired') ){
        AMULETY_BALMPAIRED_SIM(
            SIMULATE_CONVERGENCE.out.sequences,
            params.embedding_chain
        )
        ch_sim_embeddings = ch_sim_embeddings
                        .mix(AMULETY_BALMPAIRED_SIM.out.embedding
                        .map{ it -> ["balmpaired", it[0].id, it[0], it[1]] })
    }

    ch_repertoire_embeddings.dump(tag: 'repertoire_embeddings') // Debugging
    ch_sim_embeddings.dump(tag: 'simulated_embeddings') // Debugging
    ch_repertoire_meta = PREPROCESS_REPERTOIRE.out.repertoire
                        .map { it -> [it[0].id, it[0], it[1]] } // channel: [ [meta.id, meta, repertoire] ]
                        .dump(tag: 'repertoire') // Debugging
    ch_simulation_meta = SIMULATE_CONVERGENCE.out.sequences
                        .map { it -> [it[0].id, it[0], it[1]] } // channel: [ [meta.id, meta, simulated_seqs] ]
                        .dump(tag: 'simulated_sequences') // Debugging

    ch_rep_sim_meta = ch_repertoire_meta.join(ch_simulation_meta, by: 0)
        .map { it -> [it[1].id, it[1], it[2], it[4]] } // channel: [ [meta.id, meta, repertoire, simulated_seqs] ]
        .dump(tag: 'repertoire_sim_meta') // Debugging

    // Join repertoire embedding and simulated embeddings from the same embedding type and sample id
    ch_embeddings = ch_repertoire_embeddings.join(ch_sim_embeddings, by: [0,1], failOnMismatch: true)
                .dump(tag: 'embeddings joined')
                .map { it -> [it[1], it[0], it[2], it[3], it[5]] } // channel: [ [meta.id, embedding_type, meta, repertoire_embedding, simulation_embedding] ]
                .dump(tag: 'embeddings joined mapped')


    // Combine all pairs of rep_sim_meta and embeddings for the same sample
    ch_rep_sim_meta_embeddings = ch_rep_sim_meta.cross(ch_embeddings)
        .map { it ->
        def fmeta = [:] // create new meta map with id containing sample id and model type
        fmeta.id = it[0][0] + '_' + it[1][1]
        fmeta.model = it[1][1]
        fmeta.sample_id = it[0][0]
        [fmeta, it[0][2], it[0][3], it[1][3], it[1][4]] } // channel: [ [meta, repertoire, simulated_seqs, repertoire_embedding, simulation_embedding] ]
        .dump(tag: 'rep_sim_meta_embeddings') // Debugging

    ch_f_a = ch_fuzziness.combine(ch_abundance)
    ch_f_a_s = ch_f_a.combine(ch_repertoire_sample)
    ch_simulation_params = ch_f_a_s.combine(ch_witness)
                                    .dump(tag: 'simulation_params') // Debugging

    ch_rep_sim_meta_embeddings = ch_rep_sim_meta_embeddings
        .combine(ch_simulation_params)
        .map { it ->
        def fmeta = [:]
        fmeta.id = it[0].id + '_f' + it[5].toString() + '_a' + it[6].toString() + '_s' + it[7].toString() + '_w' + it[8].toString()
        fmeta.model = it[0].model
        fmeta.sample_id = it[0].sample_id
        fmeta.fuzziness = it[5]
        fmeta.abundance = it[6]
        fmeta.repertoire_sample = it[7]
        fmeta.witness = it[8]
        [fmeta, it[1], it[2], it[3], it[4]] } // channel: [ [meta, repertoire, simulated_seqs, repertoire_embedding, simulation_embedding] ]
        .dump(tag: 'rep_sim_meta_embeddings final') // Debugging

    SELECT_SIMULATED_SEQUENCES (
        ch_rep_sim_meta_embeddings
    )
    ch_versions = ch_versions.mix(SELECT_SIMULATED_SEQUENCES.out.versions.first())


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
