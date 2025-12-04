/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES
//

/*--- ALIGN TOOLS ---*/
include { MAFFT_ALIGN as MAFFT_GENE        } from '../modules/nf-core/mafft/align/main.nf'
include { MAFFT_ALIGN as MAFFT_SEGMENT     } from '../modules/nf-core/mafft/align/main.nf'
include { SUBSEQ                           } from '../modules/local/subseq/main.nf'

/*--- CONCATENATE TOOLS ---*/
include { SEQKIT_CONCAT as SEQKIT_CONCAT_L } from '../modules/local/seqkit/concat.nf'
include { SEQKIT_CONCAT as SEQKIT_CONCAT_S } from '../modules/local/seqkit/concat.nf'

/*--- TREE TOOLS ---*/
include { IQTREE                           } from '../modules/nf-core/iqtree/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow LASVPHYLO {
    take:
    input_L
    input_S
    outdir
    alignment_L
    tree_L
    alignment_S
    tree_S
    input_POL
    input_Z
    input_NP
    input_GPC
    modify_list

    main:

    ch_cutted_genes   = channel.empty()
    ch_base_alignment = channel.empty()
    ch_trees          = channel.empty()
    ch_versions       = channel.empty()

    input_id_L =  "L"
    input_id_S =  "S"

    // previous alignment channel - handle optional parameters
    ch_alignment_inputs = channel.empty()

    if (alignment_L) {
        ch_alignment_inputs = ch_alignment_inputs.mix(
            channel.fromPath(alignment_L).map{it -> [[id: input_id_L], it]}
        )
    }

    if (alignment_S) {
        ch_alignment_inputs = ch_alignment_inputs.mix(
            channel.fromPath(alignment_S).map{it -> [[id: input_id_S], it]}
        )
    }

    ch_base_alignment = ch_alignment_inputs

    if (!params.skip_indiv_gene_extraction) {
        // orient & isolate genes difficult to make them into a single channel as we need to make a distinction the correct order of genes
        ch_newseq = channel.of(
            tuple([id: input_id_L, gene: "pol"], file(input_L, checkIfExists: true), file(input_POL, checkIfExists: true)),
            tuple([id: input_id_L, gene: "z"], file(input_L, checkIfExists: true), file(input_Z, checkIfExists: true)),
            tuple([id: input_id_S, gene: "np"], file(input_S, checkIfExists: true), file(input_NP, checkIfExists: true)),
            tuple([id: input_id_S, gene: "gpc"], file(input_S, checkIfExists: true), file(input_GPC, checkIfExists: true))
        ).multiMap{meta, seq, alignment ->
            fasta: [meta, alignment]
            add: [meta, seq]
        }

        MAFFT_GENE(ch_newseq.fasta, ch_newseq.add,[[:],[]], [[:],[]], [[:],[]], [[:],[]], false)


        // ch_newseq_gene.data.view{ meta,seq, alignment -> "Gene: ${meta.gene} - ID: ${meta.id} - File: ${seq} - Alignment: ${alignment}" }

        // Isolate only the new sequence (optionally) and remove regions from sequences that have only a single genome & contaminate the data.
        yml_file = modify_list ? channel.fromPath(modify_list, checkIfExists: true).collect() : []
        ch_gene = MAFFT_GENE.out.fasta.map {meta, _fasta -> meta.gene }
        SUBSEQ(MAFFT_GENE.out.fasta, ch_gene, MAFFT_GENE.out.pattern, yml_file)

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

        ch_added_alignment = ch_cutted_genes
            .join(ch_base_alignment, remainder: true)
            .map{ meta, seq, alignment ->
                alignment ? [ meta, [seq, alignment ] ] : [ meta, [seq] ]
            }
            .transpose()
            .collectFile(newLine: true, storeDir: "${outdir}/combined_alignment") {
                meta, seq ->
                [ "${meta.id}_combined.fa", seq ]
            }.map{ it ->
                def id = it.baseName.replace("_combined", "")
                [ [id: id], it ]
            }

        ch_versions = ch_versions.mix(MAFFT_GENE.out.versions)
            .mix(SUBSEQ.out.versions)
            .mix(SEQKIT_CONCAT_L.out.versions)
            .mix(SEQKIT_CONCAT_S.out.versions)

    }else {

        // if we skip the individual gene extraction, we just use the input files as they are
        ch_new_seq = channel.of(
            [[ id :input_id_L ], file(input_L, checkIfExists: true)],
            [[ id :input_id_S ], file(input_S, checkIfExists: true)]
        )
        ch_modSeqs = ch_new_seq.join(ch_base_alignment)
        ch_added_alignment = MAFFT_SEGMENT(ch_modSeqs, "all").out.fasta
        ch_versions = ch_versions.mix(MAFFT_SEGMENT.out.versions)
    }

    // make a ML tree
    if (!params.skip_tree) {

        // merge with previous tree as a guidance - handle optional parameters
        ch_tree_inputs = channel.empty()

        if (tree_L) {
            ch_tree_inputs = ch_tree_inputs.mix(
                channel.fromPath(tree_L).map{it -> [[id: input_id_L], it]}
            )
        }

        if (tree_S) {
            ch_tree_inputs = ch_tree_inputs.mix(
                channel.fromPath(tree_S).map{it -> [[id: input_id_S], it]}
            )
        }

        ch_trees = ch_tree_inputs

        ch_iqtree = ch_added_alignment.join(ch_trees, remainder:true)
            .multiMap{ meta, seq, tree ->
                seq: [meta, seq, []]
                tree: tree?: []
            }

        //Make ML tree of new seq (optional) and previous alignment where previous sequences are under a constraint tree
        IQTREE(
            ch_iqtree.seq,
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            ch_iqtree.tree,
            [],
            [],
            []
        )
        ch_versions = ch_versions.mix(IQTREE.out.versions)
    }


    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'lasvphylo_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
