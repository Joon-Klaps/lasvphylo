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

//
// MODULES
//

include { MAFFT as MAFFT_ORIENT } from '../modules/nf-core/mafft/main.nf'
include { SUBSEQ                } from '../modules/local/subseq/main.nf'

include { SEQKIT_CONCAT as SEQKIT_CONCAT_L } from '../modules/local/seqkit/concat.nf'
include { SEQKIT_CONCAT as SEQKIT_CONCAT_S } from '../modules/local/seqkit/concat.nf'

include { MUSCLE                } from '../modules/nf-core/muscle/main.nf'
include { MAFFT as MAFFT_ALIGN  } from '../modules/nf-core/mafft/main.nf'

include { IQTREE                } from '../modules/nf-core/iqtree/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//ALIGN with MAFT -- keep length and try and have only the L and GPC genes. & Create a tree.
workflow LASVPHYLO {

    ch_cutted_genes= Channel.empty()

    ch_base_alignment = Channel.of(
        [[ id :params.input_id_L ], params.alignment_L],
        [[ id :params.input_id_S ], params.alignment_S]
    )


    if (!params.skip_indiv_gene_extraction) {
            // orient & isolate genes difficult to make them into a single channel as we need to make a distinction the correct order of genes
    ch_data= Channel.of(
        [
        meta: [ id :params.input_id_L ],
        seq: params.input_L,
        alignment: params.input_POL,
        gene: "pol"
        ],
        [
        meta: [ id :params.input_id_L ],
        seq: params.input_L,
        alignment: params.input_Z,
        gene: "z"
        ],
        [
        meta: [ id :params.input_id_S ],
        seq: params.input_S,
        alignment: params.input_NP,
        gene:"np"
        ],
        [
        meta: [ id :params.input_id_S],
        seq: params.input_S,
        alignment: params.input_GPC,
        gene:"gpc"
        ]
    )
    mafft_orient_in = ch_data.branch{meta, seq, alignment, gene ->
        data: [meta + [gene: gene], seq, alignment]
        gene: gene
    }
    MAFFT_ORIENT(mafft_orient_in.data, mafft_orient_in.gene)


    //isolate ony the new sequence
    yml_file = Channel.of(file(params.yml_file, checkIfExists: true) ?: [])
    SUBSEQ(MAFFT_ORIENT.out.fasta, MAFFT_ORIENT.out.gene, MAFFT_ORIENT.out.pattern, yml_file)

    //Concat the gene segments again
    // L: RC_Pol then Z
    // S: RC_NP then GPC
    genes_isolated = SUBSEQ.out.sequences.branch{ meta, gene, sequences ->
        pol: gene == "pol"
        z: gene == "z"
        np: gene == "np"
        gpc: gene == "gpc"
    }

    SEQKIT_CONCAT_L(genes_isolated.pol,genes_isolated.z)
    SEQKIT_CONCAT_S(genes_isolated.np,genes_isolated.gpc)

    ch_cutted_genes= SEQKIT_CONCAT_L.out.fasta.mix(SEQKIT_CONCAT_S.out.fasta)
    ch_modSeqs= ch_cutted_genes.join(ch_base_alignment)

    ch_added_alignment = MUSCLE(ch_modSeqs).aligned_fasta

    }else {

        // if we skip the individual gene extraction, we just use the input files as they are
        ch_input = Channel.of(
            [[ id :params.input_id_L ], params.input_L],
            [[ id :params.input_id_S ], params.input_S]
        )
        ch_modSeqs = ch_input.join(ch_base_alignment)
        ch_added_alignment = MAFFT_ALIGN(ch_modSeqs, "all").fasta

    }


    // MAFFT align add the additional sequences to the old alignment

    // merge with previous tree as a guidance
    ch_data_tree = ch_added_alignment.join(
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
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
