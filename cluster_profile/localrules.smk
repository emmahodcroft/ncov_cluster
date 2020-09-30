localrules: filter, partition_sequences, aggregate_alignments, mask, adjust_metadata_regions, clades, colors, recency, export, incorporate_travel_history,   finalize,  rename_legacy_clades


ruleorder: finalize_swiss > finalize
ruleorder: extract_cluster > subsample

rule add_labels:
    message: "Remove extraneous colorings for main build and move frequencies"
    input:
        auspice_json = rules.incorporate_travel_history.output.auspice_json,
        tree = rules.refine.output.tree,
        clades = rules.clades.output.clade_data,
        mutations = rules.ancestral.output.node_data
    output:
        auspice_json = "results/{build_name}/ncov_with_accessions_and_travel_branches_and_labels.json",
    log:
        "logs/add_labels_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/add_labels.py \
            --input {input.auspice_json} \
            --tree {input.tree} \
            --mutations {input.mutations} \
            --clades {input.clades} \
            --output {output.auspice_json} 2>&1 | tee {log}
        """

rule finalize_swiss:
    message: "Remove extraneous colorings for main build and move frequencies"
    input:
        auspice_json = rules.add_labels.output.auspice_json,
        frequencies = rules.tip_frequencies.output.tip_frequencies_json
    output:
        auspice_json = "auspice/ncov_{build_name}.json",
        tip_frequency_json = "auspice/ncov_{build_name}_tip-frequencies.json"
    log:
        "logs/fix_colorings_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/fix-colorings.py \
            --input {input.auspice_json} \
            --output {output.auspice_json} 2>&1 | tee {log} &&
        cp {input.frequencies} {output.tip_frequency_json}
        """

rule extract_cluster:
    input:
        cluster = "cluster_profile/clusters/cluster_{build_name}.txt",
        alignment = rules.mask.output.alignment
    output:
        cluster_sample = "results/{build_name}/sample-cluster.fasta"
    run:
        from Bio import SeqIO

        with open(input.cluster) as fh:
            cluster = set([x.strip() for x in fh.readlines()])

        seq_out = open(output.cluster_sample, 'w')
        for s in SeqIO.parse(input.alignment, 'fasta'):
            if s.id in cluster:
                SeqIO.write(s, seq_out, 'fasta')

        seq_out.close()

