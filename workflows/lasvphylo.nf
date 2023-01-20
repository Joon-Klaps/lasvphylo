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

include { MAFFT as MAFFT_ORIENT_POL         } from '../modules/nf-core/mafft/main.nf'
include { MAFFT as MAFFT_ORIENT_Z           } from '../modules/nf-core/mafft/main.nf'
include { MAFFT as MAFFT_ORIENT_GPC         } from '../modules/nf-core/mafft/main.nf'
include { MAFFT as MAFFT_ORIENT_NP          } from '../modules/nf-core/mafft/main.nf'
//include { MAFFT as MAFFT_ALIGN } from '../modules/nf-core/mafft/main.nf'
include { MUSCLE                            } from '../modules/nf-core/muscle/main.nf'
include { SEQTK_SUBSEQ as SEQTK_SUBSEQ_POL  } from '../modules/nf-core/seqtk/subseq/main.nf'
include { SEQTK_SUBSEQ as SEQTK_SUBSEQ_Z    } from '../modules/nf-core/seqtk/subseq/main.nf'
include { SEQTK_SUBSEQ as SEQTK_SUBSEQ_GPC  } from '../modules/nf-core/seqtk/subseq/main.nf'
include { SEQTK_SUBSEQ as SEQTK_SUBSEQ_NP   } from '../modules/nf-core/seqtk/subseq/main.nf'

include { IQTREE } from '../modules/nf-core/iqtree/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

// TODO: ALIGN with MAFT -- keep length and try and have only the L and GPC genes. & Create a tree.
workflow LASVPHYLO {

    // orient & isolate genes difficult to make them into a single channel as we need to make a distinction the correct order of genes
    ch_data_pol= Channel.of(
        [
        meta: [ id :params.input_id_L ],
        seq: params.input_L,
        alignment: params.input_POL,
        ]
        )
    ch_data_z= Channel.of(
        [
        meta: [ id :params.input_id_L ],
        seq: params.input_L,
        alignment: params.input_Z,
        ]
        )
    ch_data_np= Channel.of(
        [
        meta: [ id :params.input_id_S ],
        seq: params.input_S,
        alignment: params.input_NP,
        ]
        )
    ch_data_gpc= Channel.of(
        [
        meta: [ id :params.input_id_S ],
        seq: params.input_S,
        alignment: params.input_GPC,
        ]
        )

    MAFFT_ORIENT_POL(ch_data_pol,"pol")
    MAFFT_ORIENT_Z(ch_data_z,"z")
    MAFFT_ORIENT_NP(ch_data_np,"np")
    MAFFT_ORIENT_GPC(ch_data_gpc,"gpc")

    //isolate ony the new sequence
    SEQTK_SUBSEQ_POL(MAFFT_ORIENT_POL.out.fasta, MAFFT_ORIENT_POL.out.pattern, "pol")
    SEQTK_SUBSEQ_Z(MAFFT_ORIENT_Z.out.fasta, MAFFT_ORIENT_Z.out.pattern, "z")
    SEQTK_SUBSEQ_NP(MAFFT_ORIENT_NP.out.fasta, MAFFT_ORIENT_NP.out.pattern, "np")
    SEQTK_SUBSEQ_GPC(MAFFT_ORIENT_GPC.out.fasta, MAFFT_ORIENT_GPC.out.pattern, "gpc")

    //Concat the gene segments again
    // L: RC_Pol then Z
    // S: RC_NP then GPC
    SEQKIT_CONCAT_L(SEQTK_SUBSEQ_POL.out.sequences,SEQTK_SUBSEQ_Z.out.sequences)
    SEQKIT_CONCAT_S(SEQTK_SUBSEQ_NP.out.sequences,SEQTK_SUBSEQ_GPC.out.sequences)

    ch_cutted_genes= SEQKIT_CONCAT_L.out.fasta.mix(SEQKIT_CONCAT_S.out.fasta)


    ch_data = Channel.of(
        [[ id :params.input_id_L ], params.alignment_L],
        [[ id :params.input_id_S ], params.alignment_S]
        )

    ch_modSeqs= ch_cutted_genes.join(ch_data)

    // MAFFT align add the additional sequences to the old alignment
    MUSCLE(ch_modSeqs)

    // merge with previous tree as a guidance
    ch_data_tree = MUSCLE.out.aligned_fasta.join(
        Channel.of(
            [[ id :params.input_id_L ], params.tree_L],
            [[ id :params.input_id_S ], params.tree_S]
            )
    )

    //IQTREE use the previous alignment and make the tree using the previous tree as a constrain to speed it up
    IQTREE(ch_data_tree, [])

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
