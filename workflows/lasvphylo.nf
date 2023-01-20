/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check input path parameters to see if they exist
def checkPathParamList = [params.input_S,params.input_L,params.alignment_S,params.alignment_L,params.tree_L,params.tree_S,params.input_POL,params.input_Z,params.input_NP,params.input_GPC ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: created locally
//

include { SEQKIT_CONCAT as SEQKIT_CONCAT_L } from '../modules/local/seqkit/concat.nf'
include { SEQKIT_CONCAT as SEQKIT_CONCAT_S } from '../modules/local/seqkit/concat.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules (and slightly modified in this case)
//

include { MAFFT as MAFFT_ORIENT_POL } from '../modules/nf-core/mafft/main.nf'
include { MAFFT as MAFFT_ORIENT_Z } from '../modules/nf-core/mafft/main.nf'
include { MAFFT as MAFFT_ORIENT_GPC } from '../modules/nf-core/mafft/main.nf'
include { MAFFT as MAFFT_ORIENT_NP} from '../modules/nf-core/mafft/main.nf'
include { MAFFT as ALIGN } from '../modules/nf-core/mafft/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

// TODO: ALIGN with MAFT -- keep length and try and have only the L and GPC genes. & Create a tree.
workflow LASVPHYLO {

    ch_versions = Channel.empty()
    ch_data_L = Channel.of(
        [
        meta: [ id :params.input_id_L ],
        seq: params.input_L,
        alignment: params.alignment_L,
        tree: params.tree_L
        ])
    ch_data_S = Channel.of(
        [
        meta: [ id :params.input_id_S ],
        seq: params.input_S,
        alignment: params.alignment_S,
        tree: params.tree_S
        ]
        )
    ch_L_ref= Channel.of(params.input_POL, params.input_Z)
    ch_S_ref= Channel.of(params.input_NP, params.input_GPC)

    // orient & isolate genes
    MAFFT_ORIENT_POL(ch_data_L, params.input_POL, "TRUE")
    MAFFT_ORIENT_Z(ch_data_L, params.input_Z, "TRUE")
    MAFFT_ORIENT_NP(ch_data_S, params.input_NP, "TRUE")
    MAFFT_ORIENT_GPC(ch_data_S, params.input_GPC, "TRUE")


    //Concat the genes again
    SEQKIT_CONCAT_L(MAFFT_ORIENT_POL.out.fasta,MAFFT_ORIENT_Z.out.fasta)
    SEQKIT_CONCAT_S(MAFFT_ORIENT_NP.out.fasta,MAFFT_ORIENT_GPC.out.fasta)

    ch_data = Channel.of(
        [meta: [ id :params.input_id_L ], seq: SEQKIT_CONCAT_L.out.fasta, alignment: params.alignment_L, tree: params.tree_L],
        [ meta: [ id :params.input_id_S ], seq: SEQKIT_CONCAT_S.out.fasta, alignment: params.alignment_S, tree: params.tree_S]
        )


    // MAFFT align add the additional sequences to the old alignment
    MAFFT_ALIGN(ch_data,[],"FALSE")

    //IQTREE use the previous alignment and make the tree using the previous tree as a constrain to speed it up
    IQTREE(MAFFT_ALIGN.out.addseq)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
