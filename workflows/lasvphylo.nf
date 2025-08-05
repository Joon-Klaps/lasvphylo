/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

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
        ch_newseq_gene = Channel.of(
            tuple([id: params.input_id_L], file(params.input_L, checkIfExists: true), file(params.input_POL, checkIfExists: true), "pol"),
            tuple([id: params.input_id_L], file(params.input_L, checkIfExists: true), file(params.input_Z, checkIfExists: true), "z"),
            tuple([id: params.input_id_S], file(params.input_S, checkIfExists: true), file(params.input_NP, checkIfExists: true), "np"),
            tuple([id: params.input_id_S], file(params.input_S, checkIfExists: true), file(params.input_GPC, checkIfExists: true), "gpc")
        ).multiMap{meta, seq, alignment, gene ->
            data: [meta + ['gene': gene], seq, alignment]
            gene: gene
        }
        MAFFT_ORIENT(ch_newseq_gene.data, ch_newseq_gene.gene)

        ch_newseq_gene.data.view{ meta,seq, alignment -> "Gene: ${meta.gene} - ID: ${meta.id} - File: ${seq} - Alignment: ${alignment}" }

        //isolate ony the new sequence
        yml_file = Channel.of(file(params.modify_list, checkIfExists: true) ?: [])
        SUBSEQ(MAFFT_ORIENT.out.fasta, MAFFT_ORIENT.out.gene, MAFFT_ORIENT.out.pattern, yml_file)

        //Concat the gene segments again
        // L: RC_Pol then Z
        // S: RC_NP then GPC
        genes_isolated = SUBSEQ.out.fasta.branch{ meta, gene, sequences ->
            pol: gene == "pol"
                return [meta, sequences]
            z: gene == "z"
                return [meta, sequences]
            np: gene == "np"
                return [meta, sequences]
            gpc: gene == "gpc"
                return [meta, sequences]
        }

        SEQKIT_CONCAT_L(genes_isolated.pol,genes_isolated.z)
        SEQKIT_CONCAT_S(genes_isolated.np,genes_isolated.gpc)

        ch_cutted_genes= SEQKIT_CONCAT_L.out.fasta.mix(SEQKIT_CONCAT_S.out.fasta)
        ch_modSeqs= ch_cutted_genes.join(ch_base_alignment)

        ch_added_alignment = MUSCLE(ch_modSeqs).aligned_fasta

    }else {

        // if we skip the individual gene extraction, we just use the input files as they are
        ch_new_seq = Channel.of(
            [[ id :params.input_id_L ], file(params.input_L, checkIfExists: true)],
            [[ id :params.input_id_S ], file(params.input_S, checkIfExists: true)]
        )
        ch_modSeqs = ch_new_seq.join(ch_base_alignment)
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

    ch_data_tree.view{ it -> "ID: ${it[0].id} - Alignment: ${it[1]} - Tree: ${it[2]}" }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
