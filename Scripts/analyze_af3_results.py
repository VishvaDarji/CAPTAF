#!/usr/bin/env python3
"""
CAPTAF: Analyze AlphaFold3 Results and Rank PTM Variants

Extracts confidence metrics from AlphaFold3 outputs and ranks variants
based on predicted binding strength.

Metrics used:
- ipTM (interface predicted TM-score): Primary ranking metric
- pTM (predicted TM-score): Overall fold quality
- pLDDT (per-atom local distance difference test): Local confidence
- Interface PAE: Position uncertainty at interface
- Interface contacts: Number of residue-residue contacts
"""

import argparse
import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import warnings

warnings.filterwarnings('ignore')

###############################################################################
# AlphaFold3 Output Parsing
###############################################################################

def extract_af3_metrics(af3_output_dir: Path) -> Dict:
    """Extract confidence metrics from AlphaFold3 output directory"""
    
    metrics = {
        'ipTM': np.nan,
        'pTM': np.nan,
        'pLDDT_mean': np.nan,
        'interface_PAE': np.nan,
        'contacts': 0,
        'ranking_score': np.nan
    }

    try:
        # Find summary confidences file
        summary_files = list(af3_output_dir.glob("*_summary_confidences.json"))
        if not summary_files:
            print(f"  WARNING: No summary_confidences.json found in {af3_output_dir}")
            return metrics

        summary_file = summary_files[0]

        # Load summary confidences
        with open(summary_file, 'r') as f:
            summary = json.load(f)

        # Extract main metrics
        metrics['ipTM'] = summary.get('iptm', np.nan)
        metrics['pTM'] = summary.get('ptm', np.nan)
        metrics['ranking_score'] = summary.get('ranking_score', np.nan)

        # Get pairwise ipTM
        chain_pair_iptm = summary.get('chain_pair_iptm', [[]])
        if len(chain_pair_iptm) >= 2 and len(chain_pair_iptm[0]) >= 2:
            metrics['ipTM_pairwise'] = chain_pair_iptm[0][1]

        # Get interface PAE
        chain_pair_pae_min = summary.get('chain_pair_pae_min', [[]])
        if len(chain_pair_pae_min) >= 2 and len(chain_pair_pae_min[0]) >= 2:
            metrics['interface_PAE'] = chain_pair_pae_min[0][1]

        # Find full confidences file for detailed metrics
        conf_files = [f for f in af3_output_dir.glob("*_confidences.json") if "summary" not in f.name]
        
        if conf_files:
            with open(conf_files[0], 'r') as f:
                full_data = json.load(f)

            # Calculate mean pLDDT
            atom_plddts = full_data.get('atom_plddts', [])
            print(f"  DEBUG: Found {len(atom_plddts)} atom_plddts values")
            
            if atom_plddts and len(atom_plddts) > 0:
                # Convert to numpy array explicitly
                plddts_array = np.array(atom_plddts, dtype=float)
                metrics['pLDDT_mean'] = float(np.mean(plddts_array))
                print(f"  DEBUG: Calculated pLDDT mean = {metrics['pLDDT_mean']:.2f}")
            else:
                print(f"  DEBUG: atom_plddts is empty or missing!")

            # Count interface contacts
            contact_probs = full_data.get('contact_probs', [])
            print(f"  DEBUG: Found contact_probs with length {len(contact_probs)}")
            
            if contact_probs and len(contact_probs) > 0:
                contact_probs = np.array(contact_probs)
                metrics['contacts'] = int(np.sum(contact_probs > 0.5))
                print(f"  DEBUG: Calculated {metrics['contacts']} contacts")
            else:
                print(f"  DEBUG: contact_probs is empty or missing!")
        else:
            print(f"  DEBUG: No *_confidences.json files found!")

    except Exception as e:
        print(f"  ERROR parsing {af3_output_dir}: {e}")
        import traceback
        traceback.print_exc()

    return metrics

def calculate_composite_score(metrics: Dict) -> float:
    """
    Calculate composite ranking score from individual metrics
    
    Weights:
    - ipTM: 50% (most important - measures interface quality)
    - Interface PAE: 30% (validates interface positions)
    - Contacts: 20% (confirms binding surface size)
    
    Quality filters:
    - pTM > 0.5 (overall fold must be reasonable)
    - pLDDT > 70 (structure quality must be good)
    """
    
    # Check quality thresholds
    if metrics['pTM'] < 0.5:
        return 0.0  # Overall fold is likely wrong
    
    if metrics['pLDDT_mean'] < 70:
        return 0.0  # Structure quality too low
    
    # Calculate normalized scores
    iptm_score = metrics['ipTM']
    
    # PAE score (normalize: lower PAE = better, target ~8Å)
    pae_score = max(0, 1 - metrics['interface_PAE'] / 15) if not np.isnan(metrics['interface_PAE']) else 0.5
    
    # Contact score (normalize to ~50 contacts as good)
    contact_score = min(1, metrics['contacts'] / 50) if metrics['contacts'] > 0 else 0
    
    # Weighted composite
    composite = (
        0.5 * iptm_score +
        0.3 * pae_score +
        0.2 * contact_score
    )
    
    return composite

###############################################################################
# Variant Analysis and Ranking
###############################################################################

def rank_variants(poi: str, target: str, variant_dirs: List[Path],
                 use_controls: bool = False,
                 baseline_dir: Optional[Path] = None,
                 positive_dir: Optional[Path] = None) -> pd.DataFrame:
    """
    Analyze all variants and create ranking table
    """
    results = []
    
    print(f"\n{'='*60}")
    print(f"Analyzing {poi} + {target}")
    print(f"{'='*60}\n")
    
    # Analyze variants
    for variant_dir in variant_dirs:
        variant_name = variant_dir.name.split('_')[-1]
        
        print(f"Processing {variant_name}...")
        try:
            metrics = extract_af3_metrics(variant_dir)
            metrics['variant'] = variant_name
            metrics['type'] = 'variant'
            metrics['composite_score'] = calculate_composite_score(metrics)
            results.append(metrics)
            
            # Print extracted values
            print(f"  ipTM: {metrics['ipTM']:.3f}, pTM: {metrics['pTM']:.3f}, "
                  f"pLDDT: {metrics['pLDDT_mean']:.1f}, PAE: {metrics['interface_PAE']:.2f}Å")
            
        except Exception as e:
            print(f"  WARNING: Could not analyze {variant_name}: {e}")
            continue
    
    # Analyze controls if requested
    if use_controls:
        print("\nProcessing control structures...")
        
        if baseline_dir and baseline_dir.exists():
            print("Processing baseline...")
            try:
                metrics = extract_af3_metrics(baseline_dir)
                metrics['variant'] = 'baseline'
                metrics['type'] = 'control'
                metrics['composite_score'] = calculate_composite_score(metrics)
                results.append(metrics)
                print(f"  ipTM: {metrics['ipTM']:.3f}, pTM: {metrics['pTM']:.3f}")
            except Exception as e:
                print(f"  WARNING: Could not analyze baseline: {e}")
        
        if positive_dir and positive_dir.exists():
            print("Processing positive control...")
            try:
                metrics = extract_af3_metrics(positive_dir)
                metrics['variant'] = 'positive'
                metrics['type'] = 'control'
                metrics['composite_score'] = calculate_composite_score(metrics)
                results.append(metrics)
                print(f"  ipTM: {metrics['ipTM']:.3f}, pTM: {metrics['pTM']:.3f}")
            except Exception as e:
                print(f"  WARNING: Could not analyze positive: {e}")
    
    # Create DataFrame and sort by ipTM (primary) and composite score (secondary)
    df = pd.DataFrame(results)
    df = df.sort_values(['ipTM', 'composite_score'], ascending=[False, False])
    df = df.reset_index(drop=True)
    
    return df

def generate_recommendations(df: pd.DataFrame, poi: str, target: str,
                            use_controls: bool) -> str:
    """
    Generate human-readable recommendations
    """
    lines = []
    
    lines.append("="*60)
    lines.append(f"CAPTAF Analysis Results: {poi} + {target}")
    lines.append("="*60)
    lines.append("")
    
    # Get top variant
    variants_only = df[df['type'] == 'variant'].copy()
    if len(variants_only) == 0:
        lines.append("ERROR: No valid variant results found")
        return "\n".join(lines)
    
    top_variant = variants_only.iloc[0]
    
    lines.append(f"🏆 BEST VARIANT: {top_variant['variant']}")
    lines.append(f"   ipTM score: {top_variant['ipTM']:.3f}")
    lines.append(f"   pTM score: {top_variant['pTM']:.3f}")
    lines.append(f"   Mean pLDDT: {top_variant['pLDDT_mean']:.1f}")
    lines.append(f"   Interface PAE: {top_variant['interface_PAE']:.2f} Å")
    lines.append(f"   Interface contacts: {top_variant['contacts']:.0f}")
    lines.append("")
    
    # Runner-up
    if len(variants_only) > 1:
        second = variants_only.iloc[1]
        lines.append(f"🥈 Runner-up: {second['variant']}")
        lines.append(f"   ipTM score: {second['ipTM']:.3f}")
        lines.append("")
    
    # Interpretation based on ipTM
    if top_variant['ipTM'] > 0.8:
        lines.append("✓ HIGH CONFIDENCE: Strong binding predicted (ipTM > 0.8)")
        lines.append("  Recommended for experimental validation")
    elif top_variant['ipTM'] > 0.6:
        lines.append("~ MODERATE CONFIDENCE: Likely binding (ipTM 0.6-0.8)")
        lines.append("  Consider validation, but interpret cautiously")
    else:
        lines.append("⚠ LOW CONFIDENCE: Uncertain binding (ipTM < 0.6)")
        lines.append("  May be a failed prediction, validate carefully")
    
    lines.append("")
    
    # Quality checks
    if top_variant['pTM'] < 0.5:
        lines.append("⚠ WARNING: pTM < 0.5 suggests overall fold may be incorrect")
    
    if top_variant['pLDDT_mean'] < 70:
        lines.append("⚠ WARNING: pLDDT < 70 indicates low structure quality")
    
    if top_variant['interface_PAE'] > 8:
        lines.append("⚠ WARNING: Interface PAE > 8Å suggests uncertain interface geometry")
    
    lines.append("")
    
    # Compare to controls if available
    if use_controls:
        baseline = df[df['variant'] == 'baseline']
        positive = df[df['variant'] == 'positive']
        
        lines.append("Comparison to Controls:")
        lines.append("-" * 40)
        
        if len(baseline) > 0:
            baseline_iptm = baseline.iloc[0]['ipTM']
            improvement = ((top_variant['ipTM'] - baseline_iptm) / baseline_iptm) * 100
            lines.append(f"  Baseline (no PTMs): ipTM = {baseline_iptm:.3f}")
            lines.append(f"  Improvement: {improvement:+.1f}%")
        
        if len(positive) > 0:
            positive_iptm = positive.iloc[0]['ipTM']
            comparison = (top_variant['ipTM'] / positive_iptm) * 100
            lines.append(f"  Positive (all PTMs): ipTM = {positive_iptm:.3f}")
            lines.append(f"  Performance: {comparison:.1f}% of positive control")
        
        lines.append("")
    
    # Variant context analysis
    variant_contexts = {
        'v1': 'Buried + Conserved PTMs',
        'v2': 'Buried + Variable PTMs',
        'v3': 'Interface + Conserved PTMs',
        'v4': 'Interface + Variable PTMs',
        'v5': 'Exposed + Ordered PTMs',
        'v6': 'Disordered + Conserved PTMs',
        'v7': 'Disordered + Variable PTMs',
        'v8': 'High Confidence (all PTMs)'
    }
    
    best_context = variant_contexts.get(top_variant['variant'], 'Unknown')
    lines.append(f"Optimal PTM context: {best_context}")
    lines.append("")
    
    # Summary
    lines.append("="*60)
    lines.append("RECOMMENDATION:")
    lines.append("="*60)
    
    if top_variant['ipTM'] > 0.8 and top_variant['pTM'] > 0.5:
        lines.append(f"Prioritize {top_variant['variant']} for experimental validation.")
        lines.append("This PTM pattern is predicted to optimize binding.")
    elif top_variant['ipTM'] > 0.6:
        lines.append(f"Consider testing {top_variant['variant']}, but validate carefully.")
        lines.append("Moderate confidence prediction.")
    else:
        lines.append("Low confidence in all variants.")
        lines.append("Consider:")
        lines.append("  - Alternative protein conformations")
        lines.append("  - Different binding modes")
        lines.append("  - Experimental verification of interaction")
    
    return "\n".join(lines)

def create_html_summary(df: pd.DataFrame, poi: str, target: str,
                       output_file: Path):
    """
    Create HTML visualization of results
    """
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>CAPTAF Results: {poi} + {target}</title>
        <style>
            body {{ 
                font-family: 'Segoe UI', Arial, sans-serif; 
                margin: 40px;
                background-color: #f5f5f5;
            }}
            .container {{
                max-width: 1200px;
                margin: 0 auto;
                background-color: white;
                padding: 30px;
                border-radius: 8px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }}
            h1 {{ 
                color: #2c3e50; 
                border-bottom: 3px solid #3498db;
                padding-bottom: 10px;
            }}
            h2 {{ color: #34495e; }}
            table {{ 
                border-collapse: collapse; 
                width: 100%; 
                margin: 20px 0;
                font-size: 14px;
            }}
            th {{ 
                background-color: #3498db; 
                color: white; 
                padding: 12px;
                text-align: left;
            }}
            td {{ 
                padding: 10px; 
                border-bottom: 1px solid #ddd;
            }}
            .winner {{ 
                background-color: #d4edda; 
                font-weight: bold;
            }}
            .control {{ 
                background-color: #f8f9fa;
                font-style: italic;
            }}
            .high-conf {{ color: #27ae60; font-weight: bold; }}
            .med-conf {{ color: #f39c12; font-weight: bold; }}
            .low-conf {{ color: #e74c3c; font-weight: bold; }}
            .metric-box {{
                display: inline-block;
                padding: 15px;
                margin: 10px;
                border-radius: 5px;
                background-color: #ecf0f1;
            }}
            .metric-value {{
                font-size: 24px;
                font-weight: bold;
                color: #2c3e50;
            }}
            .metric-label {{
                font-size: 12px;
                color: #7f8c8d;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>CAPTAF Analysis Results</h1>
            <h2>{poi} binding to {target}</h2>
    """
    
    # Top variant metrics box
    top_variant = df[df['type'] == 'variant'].iloc[0]
    html += f"""
            <h3>Best Variant: {top_variant['variant']}</h3>
            <div>
                <div class="metric-box">
                    <div class="metric-value">{top_variant['ipTM']:.3f}</div>
                    <div class="metric-label">ipTM Score</div>
                </div>
                <div class="metric-box">
                    <div class="metric-value">{top_variant['pTM']:.3f}</div>
                    <div class="metric-label">pTM Score</div>
                </div>
                <div class="metric-box">
                    <div class="metric-value">{top_variant['pLDDT_mean']:.1f}</div>
                    <div class="metric-label">Mean pLDDT</div>
                </div>
                <div class="metric-box">
                    <div class="metric-value">{top_variant['interface_PAE']:.2f}Å</div>
                    <div class="metric-label">Interface PAE</div>
                </div>
            </div>
            
            <h3>All Results</h3>
            <table>
                <tr>
                    <th>Rank</th><th>Variant</th><th>Type</th>
                    <th>ipTM</th><th>pTM</th><th>pLDDT</th>
                    <th>PAE (Å)</th><th>Contacts</th><th>Score</th>
                </tr>
    """
    
    for idx, row in df.iterrows():
        rank = idx + 1
        row_class = ''
        if rank == 1 and row['type'] == 'variant':
            row_class = 'winner'
        elif row['type'] == 'control':
            row_class = 'control'
        
        # Color code ipTM
        iptm_class = 'high-conf' if row['ipTM'] > 0.8 else ('med-conf' if row['ipTM'] > 0.6 else 'low-conf')
        
        html += f"""
            <tr class="{row_class}">
                <td>{rank}</td>
                <td>{row['variant']}</td>
                <td>{row['type']}</td>
                <td class="{iptm_class}">{row['ipTM']:.3f}</td>
                <td>{row['pTM']:.3f}</td>
                <td>{row['pLDDT_mean']:.1f}</td>
                <td>{row['interface_PAE']:.2f}</td>
                <td>{row['contacts']:.0f}</td>
                <td>{row['composite_score']:.3f}</td>
            </tr>
        """
    
    html += """
            </table>
            
            <h3>Interpretation Guide</h3>
            <ul>
                <li><strong>ipTM > 0.8:</strong> <span class="high-conf">High confidence</span> - Strong binding predicted</li>
                <li><strong>ipTM 0.6-0.8:</strong> <span class="med-conf">Moderate confidence</span> - Likely binding</li>
                <li><strong>ipTM < 0.6:</strong> <span class="low-conf">Low confidence</span> - Uncertain binding</li>
            </ul>
            
            <p><em>Generated by CAPTAF - Context-Aware Snakemake Pipeline for PTM-Driven Binding Optimization with AlphaFold3</em></p>
        </div>
    </body>
    </html>
    """
    
    with open(output_file, 'w') as f:
        f.write(html)

###############################################################################
# Auto-Discovery Functions
###############################################################################

def auto_discover_folders(af3_output_dir: Path, poi: str, target: str) -> Dict[str, List[Path]]:
    """
    Automatically discover variant and control folders
    
    Returns:
        Dict with 'variants', 'baseline', 'positive' keys
    """
    af3_output_dir = Path(af3_output_dir)
    
    discovered = {
        'variants': [],
        'baseline': None,
        'positive': None
    }
    
    # Pattern for this POI-Target pair
    pattern = f"{poi}_{target}_"
    
    print(f"Searching for folders matching: {pattern}*")
    
    for folder in af3_output_dir.iterdir():
        if not folder.is_dir():
            continue
        
        folder_name = folder.name
        
        if not folder_name.startswith(pattern):
            continue
        
        # Determine type
        if 'baseline' in folder_name.lower():
            discovered['baseline'] = folder
            print(f"  Found baseline: {folder_name}")
        elif 'positive' in folder_name.lower():
            discovered['positive'] = folder
            print(f"  Found positive control: {folder_name}")
        elif folder_name.endswith(tuple([f'_v{i}' for i in range(1, 9)])):
            discovered['variants'].append(folder)
            print(f"  Found variant: {folder_name}")
    
    print(f"\nDiscovered: {len(discovered['variants'])} variants")
    if discovered['baseline']:
        print(f"            1 baseline control")
    if discovered['positive']:
        print(f"            1 positive control")
    
    return discovered

###############################################################################
# Main
###############################################################################

def main():
    parser = argparse.ArgumentParser(
        description="CAPTAF: Analyze AlphaFold3 results and rank PTM variants",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example (auto-discover mode):
  python analyze_af3_results.py \\
      --poi P27348 --target Q53ET0 \\
      --af3_dir alphafold_outputs \\
      --output results/P27348_Q53ET0

Example (manual mode):
  python analyze_af3_results.py \\
      --poi P27348 --target Q53ET0 \\
      --variant_dirs alphafold_outputs/P27348_Q53ET0_v*/ \\
      --output results/P27348_Q53ET0
        """
    )
    parser.add_argument("--poi", required=True, help="POI ID")
    parser.add_argument("--target", required=True, help="Target ID")
    parser.add_argument("--af3_dir", type=Path, 
                       help="AlphaFold3 outputs directory (auto-discover mode)")
    parser.add_argument("--variant_dirs", nargs='+',
                       help="Paths to variant directories (manual mode)")
    parser.add_argument("--use_controls", type=lambda x: x.lower() == 'true',
                       default=True, help="Include controls if found (default: true)")
    parser.add_argument("--baseline_dir", type=Path, help="Baseline control directory (manual)")
    parser.add_argument("--positive_dir", type=Path, help="Positive control directory (manual)")
    parser.add_argument("--output", required=True, 
                       help="Output prefix (e.g., results/POI_TARGET)")
    
    args = parser.parse_args()
    
    # Determine mode
    if args.af3_dir:
        # Auto-discover mode
        print("\n" + "="*60)
        print("CAPTAF - Auto-Discovery Mode")
        print("="*60 + "\n")
        
        discovered = auto_discover_folders(args.af3_dir, args.poi, args.target)
        
        if len(discovered['variants']) == 0:
            print(f"\nERROR: No variant folders found for {args.poi}_{args.target}")
            print(f"Searched in: {args.af3_dir}")
            return 1
        
        variant_dirs = discovered['variants']
        baseline_dir = discovered['baseline']
        positive_dir = discovered['positive']
        
        # Auto-enable controls if found
        use_controls = args.use_controls and (baseline_dir or positive_dir)
        
    elif args.variant_dirs:
        # Manual mode
        print("\n" + "="*60)
        print("CAPTAF - Manual Mode")
        print("="*60 + "\n")
        
        variant_dirs = [Path(d) for d in args.variant_dirs]
        baseline_dir = args.baseline_dir
        positive_dir = args.positive_dir
        use_controls = args.use_controls
    else:
        print("ERROR: Must provide either --af3_dir or --variant_dirs")
        parser.print_help()
        return 1
    
    # Analyze and rank
    print("\n" + "="*60)
    print("CAPTAF - AlphaFold3 Results Analysis")
    print("="*60)
    
    df = rank_variants(
        args.poi,
        args.target,
        variant_dirs,
        use_controls,
        baseline_dir,
        positive_dir
    )
    
    # Generate output filenames
    output_ranking = f"{args.output}_ranking.tsv"
    output_summary = f"{args.output}_summary.html"
    output_recommendations = f"{args.output}_recommendations.txt"
    
    # Create output directory if needed
    output_dir = Path(args.output).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save outputs
    print("\nSaving results...")
    df.to_csv(output_ranking, sep='\t', index=False)
    print(f"✓ Ranking table: {output_ranking}")
    
    recommendations = generate_recommendations(df, args.poi, args.target,
                                              use_controls)
    with open(output_recommendations, 'w') as f:
        f.write(recommendations)
    print(f"✓ Recommendations: {output_recommendations}")
    
    create_html_summary(df, args.poi, args.target, Path(output_summary))
    print(f"✓ HTML summary: {output_summary}")
    
    # Print top result
    print("\n" + recommendations)
    print("\n" + "="*60)
    print("Analysis complete!")
    print("="*60 + "\n")
    
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
