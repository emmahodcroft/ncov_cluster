localrules: clade_files, colors, download

ruleorder: adjust_metadata_regions_ecdc > adjust_metadata_regions

rule reassign_metadata:
    input:
        metadata = rules.download.output.metadata
    output:
        metadata = "results/metadata_with_reassignments.tsv"
    run:
        import pandas as pd

        d = pd.read_csv(input.metadata, sep='\t')
        ind = d['country'] == 'Israel'
        d.loc[ind,'region'] = 'Europe'
        d.loc[ind,'region_exposure'] = 'Europe'
        d.to_csv(output.metadata, sep='\t')

rule adjust_metadata_regions_ecdc:
    message:
        """
        Adjusting metadata for build '{wildcards.build_name}'
        """
    input:
        metadata = rules.reassign_metadata.output.metadata
    output:
        metadata = "results/{build_name}/metadata_adjusted.tsv"
    params:
        region = lambda wildcards: config["builds"][wildcards.build_name]["region"]
    log:
        "logs/adjust_metadata_regions_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/adjust_regional_meta.py \
            --region {params.region:q} \
            --metadata {input.metadata} \
            --output {output.metadata} 2>&1 | tee {log}
        """

ruleorder: finalize_ecdc > finalize

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

rule finalize_ecdc:
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
