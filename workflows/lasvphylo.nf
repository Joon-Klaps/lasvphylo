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

    input_id_L = file(params.input_L, checkIfExists: true).basename ?: "L"
    input_id_S = file(params.input_S, checkIfExists: true).basename ?: "S"

    // previous alingmnet channel
    ch_base_alignment = Channel.of(
        [[ id :input_id_L ], file(params.alignment_L) ?: []],
        [[ id :input_id_S ], file(params.alignment_S) ?: []]
    )

    if (!params.skip_indiv_gene_extraction) {
        // orient & isolate genes difficult to make them into a single channel as we need to make a distinction the correct order of genes
        ch_newseq_gene = Channel.of(
            tuple([id: input_id_L], file(params.input_L, checkIfExists: true), file(params.input_POL, checkIfExists: true), "pol"),
            tuple([id: input_id_L], file(params.input_L, checkIfExists: true), file(params.input_Z, checkIfExists: true), "z"),
            tuple([id: input_id_S], file(params.input_S, checkIfExists: true), file(params.input_NP, checkIfExists: true), "np"),
            tuple([id: input_id_S], file(params.input_S, checkIfExists: true), file(params.input_GPC, checkIfExists: true), "gpc")
        ).multiMap{meta, seq, alignment, gene ->
            data: [meta, seq, alignment]
            gene: gene
        }
        MAFFT_ORIENT(ch_newseq_gene.data, ch_newseq_gene.gene)

        // ch_newseq_gene.data.view{ meta,seq, alignment -> "Gene: ${meta.gene} - ID: ${meta.id} - File: ${seq} - Alignment: ${alignment}" }

        // Isolate ony the new sequence (optionally) and remove regions from sequences that have only a single genome & contaminate the data.
        yml_file = Channel.value(file(params.modify_list, checkIfExists: true) ?: [])
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

        // Concatenate genomes
        SEQKIT_CONCAT_L(genes_isolated.pol,genes_isolated.z)
        SEQKIT_CONCAT_S(genes_isolated.np,genes_isolated.gpc)

        ch_cutted_genes= SEQKIT_CONCAT_L.out.fasta.mix(SEQKIT_CONCAT_S.out.fasta)
        ch_modSeqs= ch_cutted_genes.join(ch_base_alignment)

        ch_added_alignment = ch_cutted_genes

    }else {

        // if we skip the individual gene extraction, we just use the input files as they are
        ch_new_seq = Channel.of(
            [[ id :input_id_L ], file(params.input_L, checkIfExists: true)],
            [[ id :input_id_S ], file(params.input_S, checkIfExists: true)]
        )
        ch_modSeqs = ch_new_seq.join(ch_base_alignment)
        ch_added_alignment = MAFFT_ALIGN(ch_modSeqs, "all").fasta
    }

    // make a ML tree
    if (!params.skip_tree) {
        // merge with previous tree as a guidance
        ch_iqtree = ch_added_alignment.join(
            Channel.of(
                [[ id :input_id_L ], file(params.tree_L) ?: []],
                [[ id :input_id_S ], file(params.tree_S) ?: []]
                )
        ).multiMap{ meta, seq, tree ->
            seq: [meta, seq]
            tree: [meta, tree]
        }

        //Make ML tree of new seq (optional) and previous alignment where previous sequences are under a constraint tree
        IQTREE(ch_iqtree.seq, ch_iqtree.tree, [])

    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
