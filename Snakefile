# --- config & env names ---
configfile: "config.yaml"

# Database paths (configurable)
MMSEQS_DB = config.get("mmseqs_db", "alphafold3/database/uniref50DB")
AF3_MODELS = config.get("af3_models", "af3_models")
AF3_DB = config.get("af3_db", "alphafold3/database")

# Clear problematic LD_PRELOAD on HPC
import os
if 'LD_PRELOAD' in os.environ:
    del os.environ['LD_PRELOAD']

MUSITE_ENV   = "musite_env"       # TF1/Keras env used for MusiteDeep
PIPE_ENV    = "new_thesis_env"    # your general pipeline/python env
SCRAPE_ENV  = "scraping_env"      # env for dbPTM scraping (selenium/requests/bs4)

POIs    = config["samples"]
TARGETS = config["targets"]

# Create pairs mapping for cleaner access
POI_TARGET_PAIRS = [(poi, target.replace('_Target', '')) for poi, target in zip(POIs, TARGETS)]

# Expand POI-Target pairs with variants - FIXED to include v8
POI_TARGET_VARIANT_COMBINATIONS = []
for poi, target in POI_TARGET_PAIRS:
    for variant in ['v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7', 'v8']:
        POI_TARGET_VARIANT_COMBINATIONS.append(f"{poi}_{target}_{variant}")

VARIANTS = ['v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7', 'v8']

OUTPUT_DIR = config.get("output_dir", ".")

# Prepend OUTPUT_DIR to all output directory paths
ptm_dir         = f"{OUTPUT_DIR}/{config['output_dirs']['ptm_predictions']}"
rsa_dir         = f"{OUTPUT_DIR}/{config['output_dirs']['rsa']}"
dssp_dir        = f"{OUTPUT_DIR}/{config['output_dirs']['dssp']}"
filtered_dir    = f"{OUTPUT_DIR}/{config['output_dirs']['filtered']}"
json_dir        = f"{OUTPUT_DIR}/{config['output_dirs']['jsons']}"
target_ptm_dir  = f"{OUTPUT_DIR}/{config['output_dirs']['target_ptms']}"

# Define additional output directories (hardcoded in original rules)
conservation_dir         = f"{OUTPUT_DIR}/conservation_outputs"
disorder_dir            = f"{OUTPUT_DIR}/disordered_regions_outputs"
master_tables_dir       = f"{OUTPUT_DIR}/master_tables_poi"
poi_variants_dir        = f"{OUTPUT_DIR}/poi_variants"
alphafold_outputs_dir   = f"{OUTPUT_DIR}/alphafold_outputs"
baseline_control_dir    = f"{OUTPUT_DIR}/baseline_control_case"
positive_control_dir    = f"{OUTPUT_DIR}/positive_control_case"
poi_validated_ptms_dir  = f"{OUTPUT_DIR}/poi_validated_ptms"
results_dir             = OUTPUT_DIR  # Final results go in output root

# PTM model names present in MusiteDeep_web/MusiteDeep/models/
PTM_MODELS = [
    "Phosphoserine_Phosphothreonine",
    "Phosphotyrosine",
    "N6-acetyllysine",
    "Methyllysine",
    "Methylarginine",
    "Ubiquitination",
    "SUMOylation",
    "N-linked_glycosylation",
    "O-linked_glycosylation",
    "Hydroxyproline",
    "Hydroxylysine",
    "S-palmitoyl_cysteine",
    "Pyrrolidone_carboxylic_acid",
]

# Helper function to convert short variant code to full variant name
def get_variant_full_name(variant_code):
    """
    Convert short variant code (v1, v2, ...) to full variant name.

    Args:
        variant_code: Short code like 'v1', 'v2', etc.

    Returns:
        Full variant name like 'v1_buried_conserved'
    """
    variant_map = {
        'v1': 'v1_buried_conserved',
        'v2': 'v2_buried_variable',
        'v3': 'v3_interface_conserved',
        'v4': 'v4_interface_variable',
        'v5': 'v5_exposed_ordered',
        'v6': 'v6_disordered_conserved',
        'v7': 'v7_disordered_variable',
        'v8': 'v8_high_confidence'
    }
    return variant_map.get(variant_code, variant_code)

VARIANTS = ['v1_buried_conserved', 'v2_buried_variable',
            'v3_interface_conserved', 'v4_interface_variable',
            'v5_exposed_ordered', 'v6_disordered_conserved',
            'v7_disordered_variable', 'v8_high_confidence']

# Note: Only v8_high_confidence is kept, v8_all_ptms_filtered is excluded

rule all:
    input:
        # PTM predictions (already uses ptm_dir variable)
        expand(f"{ptm_dir}/{{sample}}.tsv", sample=POIs),
        # DSSP structural predictions (already uses dssp_dir)
        expand(f"{dssp_dir}/raw_dssp_{{sample}}.tsv", sample=POIs),
        # Target PTM data (already uses target_ptm_dir)
        expand(f"{target_ptm_dir}/{{target}}_dbptm.csv", target=TARGETS),
        expand(f"{target_ptm_dir}/{{target}}_processed.csv", target=TARGETS),
        # Conservation and disorder outputs - NOW USES VARIABLES
        expand(f"{conservation_dir}/conservation_scores_{{sample}}.tsv", sample=POIs),
        expand(f"{disorder_dir}/disorder_predictions_{{sample}}.tsv", sample=POIs),
        # Merged PTM+DSSP data (already uses dssp_dir)
        expand(f"{dssp_dir}/merged_ptm_dssp_{{sample}}.tsv", sample=POIs),
        # Master tables - NOW USES VARIABLE
        expand(f"{master_tables_dir}/master_table_{{sample}}.tsv", sample=POIs),
        # Variants - NOW USES VARIABLE
        expand(f"{poi_variants_dir}/{{sample}}_{{variant}}.tsv", sample=POIs, variant=VARIANTS),
        # JSONs (already uses json_dir)
        expand(f"{json_dir}/{{combo}}.json", combo=POI_TARGET_VARIANT_COMBINATIONS),
        # AlphaFold outputs - NOW USES VARIABLE
        expand(f"{alphafold_outputs_dir}/{{combo}}", combo=POI_TARGET_VARIANT_COMBINATIONS),
        # CONTROL CASES - NOW USE VARIABLES
        expand(f"{baseline_control_dir}/{{poi}}_{{target}}_baseline.json",
               zip, poi=POIs, target=[t.replace('_Target', '') for t in TARGETS]),
        expand(f"{poi_validated_ptms_dir}/{{sample}}_dbptm.csv", sample=POIs),
        expand(f"{poi_validated_ptms_dir}/{{sample}}_processed.csv", sample=POIs),
        expand(f"{positive_control_dir}/{{poi}}_{{target}}_positive_control.json",
               zip, poi=POIs, target=[t.replace('_Target', '') for t in TARGETS]),
        expand(f"{alphafold_outputs_dir}/{{poi}}_{{target}}_baseline",
               zip, poi=POIs, target=[t.replace('_Target', '') for t in TARGETS]),
        expand(f"{alphafold_outputs_dir}/{{poi}}_{{target}}_positive_control",
               zip, poi=POIs, target=[t.replace('_Target', '') for t in TARGETS]),
        # Alphafold3 results - NOW USES VARIABLE
        expand(f"{results_dir}/{{poi}}_{{target}}_ranking.tsv",
               zip, poi=POIs, target=[t.replace('_Target', '') for t in TARGETS]),
        expand(f"{results_dir}/{{poi}}_{{target}}_summary.html",
               zip, poi=POIs, target=[t.replace('_Target', '') for t in TARGETS]),
        expand(f"{results_dir}/{{poi}}_{{target}}_recommendations.txt",
               zip, poi=POIs, target=[t.replace('_Target', '') for t in TARGETS])

# --- PTM prediction using MusiteDeep (POIs only) ---
rule predict_ptms:
    input:
        fasta=lambda wc: config["poi_fasta"]  # Always POI for this rule
    output:
        tsv_file=f"{ptm_dir}/{{protein}}.tsv"
    params:
        protein="{protein}"
    shell:
        """
        set -euo pipefail

        # Clear LD_PRELOAD to reduce warnings
        unset LD_PRELOAD 2>/dev/null || true

        # Ensure output directory exists
        mkdir -p {ptm_dir}
        mkdir -p {ptm_dir}/{params.protein}_ptm

        echo "🚀 Starting PTM prediction for {params.protein}"
        echo "   Input: {input.fasta}"
        echo "   Working dir: {ptm_dir}/{params.protein}_ptm"
        echo "   Final output: {output.tsv_file}"

        # Run the updated MusiteDeep script
        conda run -n {MUSITE_ENV} python {config[scripts_dir]}/run_musitedeep_all.py {input.fasta} {ptm_dir}/{params.protein}_ptm

        # Debug: Show what was created
        echo "Contents of output directory:"
        ls -la {ptm_dir}/{params.protein}_ptm/

        # Find the results file (based on your output pattern)
        RESULTS_FILE=""
        if [ -f "{ptm_dir}/{params.protein}_ptm/{params.protein}_all_results.txt" ]; then
            RESULTS_FILE="{ptm_dir}/{params.protein}_ptm/{params.protein}_all_results.txt"
            echo "Found MusiteDeep results file: $RESULTS_FILE"
        elif [ -f "{ptm_dir}/{params.protein}_ptm/{params.protein}_results.txt" ]; then
            RESULTS_FILE="{ptm_dir}/{params.protein}_ptm/{params.protein}_results.txt"
            echo "Found alternative results file: $RESULTS_FILE"
        else
            echo "❌ No results file found! Looking for any result files:"
            find {ptm_dir}/{params.protein}_ptm/ -name "*result*" -type f
            exit 1
        fi

        # Convert the results file to TSV format compatible with downstream scripts
        echo "Converting results to TSV format..."
        conda run -n {PIPE_ENV} python {config[scripts_dir]}/convert_musitedeep_to_tsv.py "$RESULTS_FILE" {output.tsv_file}

        # Verify the output file was created and has content
        if [ -s {output.tsv_file} ]; then
            echo "✅ PTM prediction completed successfully for {params.protein}"
            echo "   Output file: {output.tsv_file}"
            echo "   File size: $(wc -l < {output.tsv_file}) lines"
            echo "   First 5 lines:"
            head -5 {output.tsv_file}
        else
            echo "❌ Output file {output.tsv_file} is empty or was not created!"
            exit 1
        fi
        """

# --- Conservation score prediction using MMseqs2 (POIs only) ---
rule predict_conservation:
    input:
        fasta=lambda wc: config["poi_fasta"]  # Always POI
    output:
        conservation=f"{conservation_dir}/conservation_scores_{{sample}}.tsv"
    params:
        db_path=config.get("mmseqs_db", "alphafold3/database/uniref50DB"),
        threads=4
    shell:
        """
        unset LD_PRELOAD 2>/dev/null || true

        echo "🧬 Predicting conservation scores for {wildcards.sample}"
        echo "   Input: {input.fasta}"
        echo "   Output: {output.conservation}"
        echo "   Database: {params.db_path}"

        # Check if database exists and has index
        if [ ! -f "{params.db_path}" ]; then
            echo "❌ MMseqs2 database not found: {params.db_path}"
            exit 1
        fi

        if [ ! -f "{params.db_path}.index" ]; then
            echo "❌ MMseqs2 database index not found: {params.db_path}.index"
            echo "Please create the index:"
            echo "  cd alphafold3/database/"
            echo "  mmseqs createindex uniref50DB tmp --threads 4"
            exit 1
        fi

        # Run conservation prediction for SINGLE protein
        conda run -n {PIPE_ENV} python {config[scripts_dir]}/predict_conservation.py \
            --fasta {input.fasta} \
            --protein_id {wildcards.sample} \
            --output_dir {conservation_dir} \
            --db_path {params.db_path} \
            --threads {params.threads} \
            --verbose

        # Verify output
        if [ -s {output.conservation} ]; then
            echo "✅ Conservation prediction completed for {wildcards.sample}"
            echo "   Output lines: $(wc -l < {output.conservation})"
            echo "   First few lines:"
            head -5 {output.conservation}
        else
            echo "❌ Conservation prediction failed for {wildcards.sample}"
            exit 1
        fi
        """

# --- Disorder prediction using IUPred2A (POIs only) ---
rule predict_disorder:
    input:
        fasta=lambda wc: config["poi_fasta"]  # Always POI
    output:
        disorder=f"{disorder_dir}/disorder_predictions_{{sample}}.tsv"
    params:
        iupred_path=config.get("iupred2a_dir", "iupred2a") + "/iupred2a.py",
        mode="long",
        threshold=0.5
    shell:
        """
        unset LD_PRELOAD 2>/dev/null || true

        echo "🔮 Predicting disorder for {wildcards.sample}"
        echo "   Input: {input.fasta}"
        echo "   Output: {output.disorder}"
        echo "   Mode: {params.mode}"
        echo "   Threshold: {params.threshold}"

        # Check if IUPred2A exists
        if [ ! -f "{params.iupred_path}" ]; then
            echo "❌ IUPred2A script not found: {params.iupred_path}"
            exit 1
        fi
        # Run disorder prediction for SINGLE protein
        conda run -n {PIPE_ENV} python {config[scripts_dir]}/predict_disorder.py \
        --fasta {input.fasta} \
        --protein_id {wildcards.sample} \
        --output_dir {disorder_dir} \
        --iupred_path {params.iupred_path} \
        --mode {params.mode} \
        --threshold {params.threshold} \
        --verbose

        # Verify output
        if [ -s {output.disorder} ]; then
            echo "✅ Disorder prediction completed for {wildcards.sample}"
            echo "   Output lines: $(wc -l < {output.disorder})"
            echo "   First few lines:"
            head -5 {output.disorder}
        else
            echo "❌ Disorder prediction failed for {wildcards.sample}"
            exit 1
        fi
        """

# --- Predict DSSP structural features (RSA and SS) from PDB files ---
rule predict_dssp:
    input:
        pdb_file=lambda wc: config["pdb_file"]  # From config
    output:
        dssp_file=f"{dssp_dir}/raw_dssp_{{sample}}.tsv"
    shell:
        """
        unset LD_PRELOAD 2>/dev/null || true
        echo "🧬 Predicting DSSP structural features for {wildcards.sample}"
        echo "   Input PDB: {input.pdb_file}"
        echo "   Output DSSP: {output.dssp_file}"
        # Check if PDB file exists
        if [ ! -f "{input.pdb_file}" ]; then
            echo "❌ PDB file not found: {input.pdb_file}"
            echo "Please place PDB files in pdb_files_poi/ directory"
            echo "   Expected filename: {wildcards.sample}.pdb"
            echo ""
            echo "Download from AlphaFold Database:"
            echo "   wget https://alphafold.ebi.ac.uk/files/AF-{wildcards.sample}-F1-model_v4.pdb -O {input.pdb_file}"
            exit 1
        fi
        # Create output directory if it doesn't exist
        mkdir -p {dssp_dir}
        # Run Biotite DSSP prediction script
        conda run -n {PIPE_ENV} python {config[scripts_dir]}/predict_dssp.py \
            --pdb_file {input.pdb_file} \
            --output_file {output.dssp_file} \
            --verbose
        # Verify output
        if [ -s {output.dssp_file} ]; then
            echo "✅ DSSP prediction completed for {wildcards.sample}"
            echo "   Output lines: $(wc -l < {output.dssp_file})"
            echo "   First few lines:"
            head -5 {output.dssp_file}
        else
            echo "❌ DSSP prediction failed for {wildcards.sample}"
            exit 1
        fi
        """

# --- Merge PTM predictions with DSSP structural data (NO FILTERING) ---
rule merge_ptm_dssp:
    input:
        ptm_file=f"{ptm_dir}/{{sample}}.tsv",
        dssp_file=f"{dssp_dir}/raw_dssp_{{sample}}.tsv"
    output:
        merged_file=f"{dssp_dir}/merged_ptm_dssp_{{sample}}.tsv"
    shell:
        """
        unset LD_PRELOAD 2>/dev/null || true
        echo "🔬 Merging PTM predictions with DSSP data for {wildcards.sample}"
        echo "   Input PTM file: {input.ptm_file}"
        echo "   Input DSSP file: {input.dssp_file}"
        echo "   Output merged: {output.merged_file}"
        # Check if input files exist
        if [ ! -f "{input.ptm_file}" ]; then
            echo "❌ PTM file not found: {input.ptm_file}"
            exit 1
        fi
        if [ ! -f "{input.dssp_file}" ]; then
            echo "❌ DSSP file not found: {input.dssp_file}"
            exit 1
        fi
        # Merge PTM with DSSP data (no filtering)
        conda run -n {PIPE_ENV} python {config[scripts_dir]}/merge_dssp_ptm.py \
            {input.ptm_file} \
            {input.dssp_file} \
            {output.merged_file} \
            --verbose
        # Check if merge succeeded
        if [ -s {output.merged_file} ]; then
            echo "✅ PTM+DSSP merge completed for {wildcards.sample}"
            echo "   Merged output lines: $(wc -l < {output.merged_file})"
            echo "   First few lines:"
            head -3 {output.merged_file}
        else
            echo "❌ PTM+DSSP merge failed for {wildcards.sample}"
            exit 1
        fi
        """

# --- Create master table by merging all data sources ---
rule create_master_table:
    input:
        ptm_dssp=f"{dssp_dir}/merged_ptm_dssp_{{sample}}.tsv",
        conservation=f"{conservation_dir}/conservation_scores_{{sample}}.tsv",
        disorder=f"{disorder_dir}/disorder_predictions_{{sample}}.tsv"
    output:
        master_table=f"{master_tables_dir}/master_table_{{sample}}.tsv"
    shell:
        """
        unset LD_PRELOAD 2>/dev/null || true

        echo "🔗 Creating master table for {wildcards.sample}"
        echo "   PTM+DSSP: {input.ptm_dssp}"
        echo "   Conservation: {input.conservation}"
        echo "   Disorder: {input.disorder}"
        echo "   Output: {output.master_table}"

        # Create master table by merging all data
        conda run -n {PIPE_ENV} python {config[scripts_dir]}/create_master_table.py \
            --protein_id {wildcards.sample} \
            --ptm_dssp_file {input.ptm_dssp} \
            --conservation_file {input.conservation} \
            --disorder_file {input.disorder} \
            --output_dir {master_tables_dir} \
            --verbose

        # Verify output
        if [ -s {output.master_table} ]; then
            echo "✅ Master table created for {wildcards.sample}"
            echo "   Output lines: $(wc -l < {output.master_table})"
        else
            echo "❌ Master table creation failed for {wildcards.sample}"
            exit 1
        fi
        """

# --- Generate 8 variants using parallel filtering ---
rule generate_variants:
    input:
        master_table=f"{master_tables_dir}/master_table_{{sample}}.tsv"
    output:
        v1=f"{poi_variants_dir}/{{sample}}_v1_buried_conserved.tsv",
        v2=f"{poi_variants_dir}/{{sample}}_v2_buried_variable.tsv",
        v3=f"{poi_variants_dir}/{{sample}}_v3_interface_conserved.tsv",
        v4=f"{poi_variants_dir}/{{sample}}_v4_interface_variable.tsv",
        v5=f"{poi_variants_dir}/{{sample}}_v5_exposed_ordered.tsv",
        v6=f"{poi_variants_dir}/{{sample}}_v6_disordered_conserved.tsv",
        v7=f"{poi_variants_dir}/{{sample}}_v7_disordered_variable.tsv",
        v8=f"{poi_variants_dir}/{{sample}}_v8_high_confidence.tsv"
    shell:
        """
        unset LD_PRELOAD 2>/dev/null || true

        echo "🔬 Generating variants for {wildcards.sample}"
        echo "   Master table: {input.master_table}"
        echo "   Output directory: {poi_variants_dir}"

        # Generate all 8 variants using parallel filtering
        conda run -n {PIPE_ENV} python {config[scripts_dir]}/generate_variants.py \
            --protein_id {wildcards.sample} \
            --master_table {input.master_table} \
            --output_dir {poi_variants_dir} \
            --verbose

        # Verify all outputs were created
        ALL_CREATED=true
        for variant_file in {output.v1} {output.v2} {output.v3} {output.v4} {output.v5} {output.v6} {output.v7} {output.v8}; do
            if [ ! -f "$variant_file" ]; then
                echo "❌ Variant file not created: $variant_file"
                ALL_CREATED=false
            fi
        done

        if [ "$ALL_CREATED" = true ]; then
            echo "✅ All 8 variants generated successfully for {wildcards.sample}"
            echo "   Files created in poi_variants/"
            ls -lh {poi_variants_dir}/{wildcards.sample}_v*.tsv
        else
            echo "❌ Some variant files were not created"
            exit 1
        fi
        """

# --- Target PTM fetching via dbPTM (no filtering) ---
rule fetch_dbptm_ptms:
    input:
        fasta=lambda wc: config["target_fasta"]  # From config
    output:
        f"{target_ptm_dir}/{{target}}_dbptm.csv"
    retries: 2
    shell:
        """
        set -euo pipefail;
        unset LD_PRELOAD 2>/dev/null || true
    
        echo "DEBUG: Checking for file: {output}"
        echo "DEBUG: Current working directory:"
        pwd
        echo "DEBUG: Listing current directory:"
        ls -lh . | head -10
        ls -lh "{output}" || echo "File not found at: {output}"
        # === CHECK IF OUTPUT ALREADY EXISTS ===
        if [ -f "target_ptm_dir/{wildcards.target}_dbptm.csv" ] && [ -s "target_ptm_dir/{wildcards.target}_dbptm.csv" ]; then
            echo "✅ Target PTM file already exists: {output}"
            line_count=$(wc -l < "target_ptm_dir/{wildcards.target}_dbptm.csv")
            echo "   File has $line_count lines (manually created or from previous run)"
            echo "   Skipping dbPTM fetch"
            # Ensure output exists at path Snakemake expects (we are in output_dir already)
            # So {output} path is relative from HERE, but file is in target_ptm_dir/
            # Just touch {output} to mark it as complete
            touch "{output}"
            exit 0
        fi
        
        echo "🌐 Fetching dbPTM data for {wildcards.target}"
        echo "   Input: {input.fasta}"
        echo "   Output: {output}"

        # Show the FASTA header to verify UniProt ID extraction
        echo "FASTA header:"
        head -1 {input.fasta}

        # Extract the expected UniProt ID for verification
        HEADER=$(head -1 {input.fasta})
        if [[ $HEADER == *"|"*"|"* ]]; then
            UNIPROT_ID=$(echo "$HEADER" | cut -d'|' -f3 | cut -d' ' -f1)
            echo "Expected UniProt ID: $UNIPROT_ID"
            echo "Will try URL: https://biomics.lab.nycu.edu.tw/dbPTM/info.php?id=$UNIPROT_ID"
        fi

        # Use the robust Selenium script with proper Chrome isolation
        conda run -n {PIPE_ENV} python {config[scripts_dir]}/fetch_dbptm_ptms.py {input.fasta} {output}

        # Check if output file exists and has content beyond header
        if [ -s {output} ]; then
            line_count=$(wc -l < {output})
            echo "Output file created with $line_count lines"

            if [ $line_count -gt 1 ]; then
                echo "✅ dbPTM fetching completed for {wildcards.target}"
                echo "   Output lines: $line_count"
                echo "   Column headers:"
                head -1 {output}
                echo "   Sample data:"
                head -3 {output} | tail -2
            else
                echo "⚠️  dbPTM fetching completed but only found header for {wildcards.target}"
                echo "   This may be normal if no PTMs are available in dbPTM"
                echo "   File content:"
                cat {output}
            fi
        else
            echo "❌ dbPTM fetching failed - creating minimal output file"
            echo "UniProt_ID,Position,PTM_Type,Source,Substrate,Location_Raw" > {output}
            echo "   Created empty output file to continue pipeline"
        fi

        echo "Final file status:"
        ls -la {output}
        """

# --- Process dbPTM files by adding AA column and removing unnecessary columns ---
rule process_dbptm_ptms:
    input:
        fasta=lambda wc: config["target_fasta"],
        dbptm_file=f"{target_ptm_dir}/{{target}}_dbptm.csv"
    output:
        processed_file=f"{target_ptm_dir}/{{target}}_processed.csv"
    shell:
        """
        unset LD_PRELOAD 2>/dev/null || true
    
        # === CHECK IF OUTPUT ALREADY EXISTS ===
        if [ -f "target_ptm_dir/{wildcards.target}_processed.csv" ] && [ -s "target_ptm_dir/{wildcards.target}_processed.csv" ]; then
            echo "✅ Processed PTM file already exists: {output.processed_file}"
            line_count=$(wc -l < "target_ptm_dir/{wildcards.target}_processed.csv")
            echo "   File has $line_count lines (manually created or from previous run)"
            echo "   Skipping processing"
            exit 0
        fi

        echo "🔬 Processing dbPTM data for {wildcards.target}"
        echo "   Target FASTA: {input.fasta}"
        echo "   dbPTM file: {input.dbptm_file}"
        echo "   Output: {output.processed_file}"

        # Check if input files exist
        if [ ! -f "{input.fasta}" ]; then
            echo "❌ Target FASTA file not found: {input.fasta}"
            exit 1
        fi

        if [ ! -f "{input.dbptm_file}" ]; then
            echo "❌ dbPTM file not found: {input.dbptm_file}"
            exit 1
        fi

        # Create a copy of the original file to work with
        cp {input.dbptm_file} {output.processed_file}

        # Process the dbPTM file
        conda run -n {PIPE_ENV} python {config[scripts_dir]}/process_dbptm_ptms.py \
            {input.fasta} \
            {output.processed_file} \
            --verbose

        # Check if processing succeeded
        if [ -s {output.processed_file} ]; then
            echo "✅ dbPTM processing completed for {wildcards.target}"
            echo "   Processed file lines: $(wc -l < {output.processed_file})"
            echo "   First few lines:"
            head -3 {output.processed_file}

            # Show column structure
            echo "   Columns in processed file:"
            head -1 {output.processed_file}
        else
            echo "❌ dbPTM processing failed or produced empty output"
            exit 1
        fi
        """

# --- JSON generation for AF3 (one JSON per variant) ---
rule generate_af3_json:
    input:
        poi_fasta = lambda wc: config["poi_fasta"],
        target_fasta = lambda wc: config["target_fasta"],
        poi_ptms     = lambda wc: f"{poi_variants_dir}/{wc.poi}_{get_variant_full_name(wc.variant)}.tsv",
        target_ptms  = lambda wc: f"{target_ptm_dir}/{wc.target}_Target_processed.csv",
        ptm_mapping=lambda wc: config["ptm_mapping"]
    output:
        json = f"{json_dir}/{{poi}}_{{target}}_{{variant}}.json"
    shell:
        """
        unset LD_PRELOAD 2>/dev/null || true

        echo "📄 Generating AF3 JSON for {wildcards.poi} + {wildcards.target} (variant: {wildcards.variant})"
        echo "   POI FASTA: {input.poi_fasta}"
        echo "   Target FASTA: {input.target_fasta}"
        echo "   POI PTMs (variant): {input.poi_ptms}"
        echo "   Target PTMs: {input.target_ptms}"
        echo "   Output JSON: {output.json}"

        # Copy PTM mapping to current directory
        #cp {input.ptm_mapping} PTM_to_CCD_mapping.csv

        # Generate JSON
        conda run -n {PIPE_ENV} python {config[scripts_dir]}/generate_af3_json.py \
            --poi_fasta {input.poi_fasta} \
            --target_fasta {input.target_fasta} \
            --poi_ptms {input.poi_ptms} \
            --target_ptms {input.target_ptms} \
            --output {output.json} \
            --verbose

        if [ -s {output.json} ]; then
            echo "✅ JSON generation completed: {output.json}"
        else
            echo "❌ JSON generation failed"
            exit 1
        fi
        """

# --- CONTROL CASES ---

# --- Fetch validated PTMs for POI proteins from dbPTM ---
rule fetch_poi_validated_ptms:
    input:
        poi_fasta=lambda wc: config["poi_fasta"],
        target_fasta=lambda wc: config["target_fasta"]
    output:
        f"{poi_validated_ptms_dir}/{{sample}}_dbptm.csv"
    retries: 2
    shell:
        """
        unset LD_PRELOAD 2>/dev/null || true

        echo "🌐 Fetching validated PTMs for POI: {wildcards.sample}"
        echo "   Input: {input.poi_fasta}"
        echo "   Output: {output}"

        # Use the existing dbPTM fetching script
        conda run -n {PIPE_ENV} python {config[scripts_dir]}/fetch_dbptm_ptms.py {input.poi_fasta} {output}

        # Check if output file exists
        if [ -s {output} ]; then
            line_count=$(wc -l < {output})
            echo "✅ POI validated PTMs fetched: {wildcards.sample}"
            echo "   Output lines: $line_count"
        else
            echo "⚠️  POI validated PTMs fetching completed with minimal data for {wildcards.sample}"
        fi
        """

# --- Process POI validated PTMs (add AA column) ---
rule process_poi_validated_ptms:
    input:
        fasta=lambda wc: config["poi_fasta"],
        dbptm_file=f"{poi_validated_ptms_dir}/{{sample}}_dbptm.csv"
    output:
        processed_file=f"{poi_validated_ptms_dir}/{{sample}}_processed.csv"
    shell:
        """
        unset LD_PRELOAD 2>/dev/null || true

        echo "🔬 Processing POI validated PTMs for {wildcards.sample}"
        echo "   FASTA: {input.fasta}"
        echo "   dbPTM file: {input.dbptm_file}"
        echo "   Output: {output.processed_file}"

        # Create a copy to work with
        cp {input.dbptm_file} {output.processed_file}

        # Process the dbPTM file
        conda run -n {PIPE_ENV} python {config[scripts_dir]}/process_dbptm_ptms.py \
            {input.fasta} \
            {output.processed_file} \
            --verbose

        if [ -s {output.processed_file} ]; then
            echo "✅ POI validated PTMs processed: {wildcards.sample}"
            echo "   Output lines: $(wc -l < {output.processed_file})"
        else
            echo "❌ POI validated PTM processing failed"
            exit 1
        fi
        """

# --- Generate baseline control JSON (no PTMs) ---
rule generate_baseline_control:
    input:
        poi_fasta=lambda wc: config["poi_fasta"],
        target_fasta=lambda wc: config["target_fasta"]
    output:
        json=f"{baseline_control_dir}/{{poi}}_{{target}}_baseline.json"
    shell:
        """
        unset LD_PRELOAD 2>/dev/null || true

        echo "📄 Generating BASELINE control JSON (no PTMs)"
        echo "   POI: {wildcards.poi}"
        echo "   Target: {wildcards.target}"
        echo "   Output: {output.json}"

        mkdir -p {baseline_control_dir}

        # Generate baseline JSON (unmodified proteins)
        conda run -n {PIPE_ENV} python {config[scripts_dir]}/generate_baseline_control.py \
            --poi_fasta {input.poi_fasta} \
            --target_fasta {input.target_fasta} \
            --output {output.json} \
            --verbose

        if [ -s {output.json} ]; then
            echo "✅ Baseline control JSON generated: {output.json}"
        else
            echo "❌ Baseline control JSON generation failed"
            exit 1
        fi
        """

# --- Generate positive control JSON (all validated PTMs) ---
rule generate_positive_control:
    input:
        poi_fasta=lambda wc: config["poi_fasta"],
        target_fasta=lambda wc: config["target_fasta"],
        poi_ptms=f"{poi_validated_ptms_dir}/{{poi}}_processed.csv",
        target_ptms=lambda wc: f"{target_ptm_dir}/{wc.target}_Target_processed.csv",
        ptm_mapping=lambda wc: config["ptm_mapping"]
    output:
        json=f"{positive_control_dir}/{{poi}}_{{target}}_positive_control.json"
    shell:
        """
        unset LD_PRELOAD 2>/dev/null || true

        echo "📄 Generating POSITIVE control JSON (all validated PTMs)"
        echo "   POI: {wildcards.poi}"
        echo "   Target: {wildcards.target}"
        echo "   POI PTMs: {input.poi_ptms}"
        echo "   Target PTMs: {input.target_ptms}"
        echo "   Output: {output.json}"

        mkdir -p {positive_control_dir}

        #copying PTM mapping file to output_dir
        #cp {input.ptm_mapping} PTM_to_CCD_mapping.csv

        # Generate positive control JSON (all validated PTMs)
        conda run -n {PIPE_ENV} python {config[scripts_dir]}/generate_positive_control.py \
            --poi_fasta {input.poi_fasta} \
            --target_fasta {input.target_fasta} \
            --poi_ptms {input.poi_ptms} \
            --target_ptms {input.target_ptms} \
            --output {output.json} \
            --verbose

        #if [ -s {output.json} ]; then
        #    echo "✅ Positive control JSON generated: {output.json}"
        #else
        #    echo "❌ Positive control JSON generation failed"
        #    exit 1
        #fi
        """

# --- Run AlphaFold3 for baseline control ---
rule run_alphafold3_baseline:
    input:
        json_file=f"{baseline_control_dir}/{{poi}}_{{target}}_baseline.json"
    output:
        output_dir=directory(f"{alphafold_outputs_dir}/{{poi}}_{{target}}_baseline")
    params:
        output_prefix=lambda wc: f"{wc.poi}_{wc.target}_baseline"
    resources:
        gpu=4
    threads: 8
    shell:
        """
        echo "🧬 Running AlphaFold3 BASELINE control for {params.output_prefix}"

        mkdir -p {alphafold_outputs_dir}

        export CUDA_VISIBLE_DEVICES=0
        export XLA_PYTHON_CLIENT_PREALLOCATE=false
        export XLA_PYTHON_CLIENT_MEM_FRACTION=0.50
        export XLA_PYTHON_CLIENT_ALLOCATOR=platform
        export TF_FORCE_GPU_ALLOW_GROWTH=true
        export JAX_PLATFORMS=cuda

        conda run -n af3 python {config[alphafold3_dir]}/run_alphafold.py \
            --json_path={input.json_file} \
            --model_dir={config[af3_models]} \
            --db_dir={config[af3_db]} \
            --output_dir={alphafold_outputs_dir} \
            --num_diffusion_samples=1

        echo "✅ AlphaFold3 baseline execution completed"

        # Find the MOST RECENTLY created directory
        ACTUAL_DIR=$(find {alphafold_outputs_dir} -maxdepth 1 -type d ! -name "alphafold_outputs" -printf '%T@ %p\n' | sort -n | tail -1 | cut -d' ' -f2-)

        if [ -z "$ACTUAL_DIR" ]; then
            echo "❌ No new directory found after AF3 run"
            exit 1
        fi

        echo "Found AF3 output: $ACTUAL_DIR"

        if [ "$ACTUAL_DIR" != "{output}" ]; then
            echo "Renaming $ACTUAL_DIR -> {output}"
            mv "$ACTUAL_DIR" {output}
        fi
        """

# --- Run AlphaFold3 for positive control ---
rule run_alphafold3_positive:
    input:
        json_file=f"{positive_control_dir}/{{poi}}_{{target}}_positive_control.json"
    output:
        output_dir=directory(f"{alphafold_outputs_dir}/{{poi}}_{{target}}_positive_control")
    params:
        output_prefix=lambda wc: f"{wc.poi}_{wc.target}_positive_control"
    resources:
        gpu=4
    threads: 8
    shell:
        """
        set -euo pipefail
        unset LD_PRELOAD 2>/dev/null || true
        echo "🧬 Running AlphaFold3 POSITIVE control for {wildcards.poi}_{wildcards.target}_positive_control"
    
        mkdir -p {alphafold_outputs_dir}
    
        export CUDA_VISIBLE_DEVICES=0
        export XLA_PYTHON_CLIENT_PREALLOCATE=false
        export XLA_PYTHON_CLIENT_MEM_FRACTION=0.50
        export XLA_PYTHON_CLIENT_ALLOCATOR=platform
        export TF_FORCE_GPU_ALLOW_GROWTH=true
        export JAX_PLATFORMS=cuda
    
        conda run -n af3 python {config[alphafold3_dir]}/run_alphafold.py \
        --json_path={input.json_file} \
        --model_dir={config[af3_models]} \
        --db_dir={config[alphafold3_dir]}/database \
        --output_dir={alphafold_outputs_dir} \
        --num_diffusion_samples=1
    
        echo "✅ AlphaFold3 positive control execution completed"
    
        # Find the MOST RECENTLY created directory in alphafold_outputs
        ACTUAL_DIR=$(find {alphafold_outputs_dir} -maxdepth 1 -type d ! -name "alphafold_outputs" -printf '%T@ %p\n' | sort -n | tail -1 | cut -d' ' -f2-)
    
        if [ -z "$ACTUAL_DIR" ]; then
            echo "❌ No new directory found after AF3 run"
            exit 1
        fi
    
        echo "Found AF3 output: $ACTUAL_DIR"
    
        # Rename to expected output name
        if [ "$ACTUAL_DIR" != "{output}" ]; then
            echo "Renaming $ACTUAL_DIR -> {output}"
            mv "$ACTUAL_DIR" {output}
        fi
    
        echo "✅ Positive control output ready: {output}"
        """

rule run_alphafold3:
    input:
        json_file = f"{json_dir}/{{poi}}_{{target}}_{{variant}}.json"
    output:
        output_dir = directory(f"{alphafold_outputs_dir}/{{poi}}_{{target}}_{{variant}}")
    params:
        output_prefix = lambda wc: f"{wc.poi}_{wc.target}_{wc.variant}"
    resources:
        gpu=4
    threads: 8
    shell:
        """
        echo "🧬 Running AlphaFold3 for {params.output_prefix}"

        mkdir -p {alphafold_outputs_dir}

        # Aggressive JAX memory management
        export CUDA_VISIBLE_DEVICES=0
        export XLA_PYTHON_CLIENT_PREALLOCATE=false
        export XLA_PYTHON_CLIENT_MEM_FRACTION=0.50
        export XLA_PYTHON_CLIENT_ALLOCATOR=platform
        export TF_FORCE_GPU_ALLOW_GROWTH=true
        export JAX_PLATFORMS=cuda
        
        conda run -n af3 python {config[alphafold3_dir]}/run_alphafold.py \
            --json_path={input.json_file} \
            --model_dir={config[af3_models]} \
            --db_dir={config[af3_db]} \
            --output_dir={alphafold_outputs_dir} \
            --num_diffusion_samples=1

        echo "✅ AlphaFold3 execution completed for {params.output_prefix}"

        # Find the MOST RECENTLY created directory
        ACTUAL_DIR=$(find {alphafold_outputs_dir} -maxdepth 1 -type d ! -name "alphafold_outputs" -printf '%T@ %p\n' | sort -n | tail -1 | cut -d' ' -f2-)

        if [ -z "$ACTUAL_DIR" ]; then
            echo "❌ No new directory found after AF3 run"
            exit 1
        fi

        echo "Found AF3 output: $ACTUAL_DIR"

        if [ "$ACTUAL_DIR" != "{output}" ]; then
            echo "Renaming $ACTUAL_DIR -> {output}"
            mv "$ACTUAL_DIR" {output}
        fi

        if [ -d "{output.output_dir}" ]; then
            echo "✅ Completed: {output.output_dir}"
        else
            echo "❌ Output directory not found!"
            ls -la alphafold_outputs/
            exit 1
        fi
        """

# --- Analysis and Ranking Rule (Final Step) ---
rule analyze_and_rank_variants:
    input:
        # Wait for all variant AF3 runs to complete for this POI-Target pair
        variant_dirs=lambda wc: expand(
            f"{alphafold_outputs_dir}/{wc.poi}_{wc.target}_{{variant}}",
            variant=["v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8"]
        ),
        # Also wait for control runs (baseline and positive control)
        baseline_dir=f"{alphafold_outputs_dir}/{{poi}}_{{target}}_baseline",
        positive_dir=f"{alphafold_outputs_dir}/{{poi}}_{{target}}_positive_control"
    output:
        ranking=f"{results_dir}/{{poi}}_{{target}}_ranking.tsv",
        summary=f"{results_dir}/{{poi}}_{{target}}_summary.html",
        recommendations=f"{results_dir}/{{poi}}_{{target}}_recommendations.txt"
    params:
        af3_dir=alphafold_outputs_dir
    shell:
        """
        conda run -n {PIPE_ENV} python {config[scripts_dir]}/analyze_af3_results.py \
            --poi {wildcards.poi} \
            --target {wildcards.target} \
            --af3_dir {params.af3_dir} \
            --output {results_dir}/{wildcards.poi}_{wildcards.target}
        """
