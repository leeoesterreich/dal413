#!/usr/bin/env python

from __future__ import print_function
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg')  # Use Agg backend for better compatibility
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import argparse
import os
import sys
import numpy as np
import json
import math
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
import plotly.express as px

class ReactomePathwayAnalyzer:
    def __init__(self, gene_pathway_file="data/NCBI2Reactome_All_Levels.txt", 
                 pathway_hierarchy_file="data/ReactomePathwaysRelation.txt"):
        """Initialize with local Reactome data files"""
        print("Loading Reactome data files...")
        
        # Load only essential columns with optimized dtypes
        dtype_dict = {
            'Gene': str,
            'Pathway_ID': str,
            'Pathway_Name': str,
            'Species': str
        }
        
        # Use chunked reading for better memory management
        chunks = []
        for chunk in pd.read_csv(gene_pathway_file, sep='\t', header=None,
                               names=['Gene', 'Pathway_ID', 'Pathway_URL', 
                                     'Pathway_Name', 'Evidence', 'Species'],
                               usecols=['Gene', 'Pathway_ID', 'Pathway_Name', 'Species'],
                               dtype=dtype_dict, chunksize=50000):
            # Filter human pathways in each chunk
            human_chunk = chunk[chunk['Species'].str.contains('Homo sapiens', na=False)]
            chunks.append(human_chunk)
        
        # Combine filtered chunks
        self.gene_pathway_df = pd.concat(chunks, ignore_index=True)
        
        # Create indices for faster lookups
        self.gene_to_pathways = defaultdict(list)
        self.pathway_to_id = {}
        
        # Build indices in a single pass
        for _, row in self.gene_pathway_df.iterrows():
            self.gene_to_pathways[row['Gene']].append((row['Pathway_Name'], row['Gene']))
            self.pathway_to_id[row['Pathway_Name']] = row['Pathway_ID']
        
        # Load pathway hierarchy with minimal memory usage
        self.pathway_hierarchy = defaultdict(set)
        with open(pathway_hierarchy_file, 'r') as f:
            for line in f:
                parent, child = line.strip().split('\t')
                self.pathway_hierarchy[parent].add(child)
        
        print(f"Loaded {len(self.gene_pathway_df)} gene-pathway associations")
        print(f"Loaded {len(self.pathway_hierarchy)} pathway relationships")
        sys.stdout.flush()
        
    def get_affected_pathways(self, gene_list):
        """Find pathways containing input genes using pre-built indices"""
        affected_pathways = defaultdict(set)
        total_genes = len(gene_list)
        
        print(f"\nAnalyzing {total_genes} genes...")
        sys.stdout.flush()
        
        # Process genes in batches for better progress reporting
        batch_size = max(1, total_genes // 20)
        for i in range(0, total_genes, batch_size):
            batch = gene_list[i:i + batch_size]
            
            # Process each gene in the batch
            for gene in batch:
                gene = str(gene).upper()
                pathways = self.gene_to_pathways.get(gene, [])
                if pathways:
                    for pathway_name, gene_id in pathways:
                        affected_pathways[pathway_name].add(gene_id)
            
            # Update progress
            processed = min(i + batch_size, total_genes)
            sys.stdout.write(f"\rProcessed {processed}/{total_genes} genes")
            sys.stdout.flush()
        
        # Convert sets to lists for compatibility
        return {k: list(v) for k, v in affected_pathways.items()}
    
    def get_pathway_hierarchy(self, pathway_list):
        """Get pathway relationships using pre-built hierarchy index"""
        edges = []
        pathway_ids = {self.pathway_to_id[name]: name for name in pathway_list 
                      if name in self.pathway_to_id}
        
        # Check relationships using the pre-built hierarchy
        for parent_id, children in self.pathway_hierarchy.items():
            if parent_id in pathway_ids:
                parent_name = pathway_ids[parent_id]
                for child_id in children:
                    if child_id in pathway_ids:
                        child_name = pathway_ids[child_id]
                        edges.append((parent_name, child_name))
        
        return edges

    def create_static_pathway_network(self, affected_pathways, output_file="output/pathway_impact_static.html"):
        """Create a static network visualization using Plotly"""
        if not affected_pathways:
            print("No pathways found to visualize.")
            return None
        
        try:
            # Ensure output is HTML
            if not output_file.endswith('.html'):
                output_file = os.path.splitext(output_file)[0] + '.html'
            
            # Filter to keep only the most significant pathways (top 100 by gene count)
            significant_pathways = dict(sorted(affected_pathways.items(), 
                                            key=lambda x: len(x[1]), 
                                            reverse=True)[:100])
            
            # Create a directed graph
            G = nx.DiGraph()
            
            # Add nodes with properties
            for pathway, genes in significant_pathways.items():
                G.add_node(pathway, size=len(genes), genes=genes)
            
            # Add hierarchical relationships
            edges = self.get_pathway_hierarchy(list(significant_pathways.keys()))
            G.add_edges_from(edges)
            
            # Define major pathway categories with detailed subcategories
            categories = {
                'Signal Transduction': {
                    'keywords': ['signal', 'cascade', 'pathway'],
                    'subcategories': {
                        'RTK Signaling': ['EGFR', 'FGFR', 'VEGF', 'PDGF', 'RTK'],
                        'MAPK Cascade': ['MAPK', 'ERK', 'RAF', 'RAS', 'MEK', 'JNK', 'p38'],
                        'WNT Signaling': ['WNT', 'beta-catenin', 'frizzled'],
                        'PI3K-AKT Signaling': ['PI3K', 'AKT', 'mTOR', 'PTEN'],
                        'JAK-STAT Signaling': ['JAK', 'STAT', 'cytokine'],
                        'NF-kB Signaling': ['NF-kB', 'IKK', 'TNF'],
                        'TGF-beta Signaling': ['TGF', 'SMAD', 'BMP'],
                        'Notch Signaling': ['NOTCH', 'delta', 'jagged'],
                        'Hedgehog Signaling': ['hedgehog', 'patched', 'smoothened', 'GLI'],
                        'G-protein Signaling': ['G-protein', 'GPCR', 'adenylate cyclase', 'cAMP']
                    }
                },
                'Cell Communication': {
                    'keywords': ['adhesion', 'junction', 'communication'],
                    'subcategories': {
                        'Cell Adhesion': ['adhesion', 'cadherin', 'integrin'],
                        'Cell Junctions': ['junction', 'gap junction'],
                        'ECM Interaction': ['matrix', 'ECM']
                    }
                },
                'Disease': {
                    'keywords': ['disease', 'cancer', 'disorder'],
                    'subcategories': {
                        'Cancer': ['cancer', 'tumor', 'oncogenic'],
                        'Immune Disorders': ['autoimmune', 'inflammation'],
                        'Metabolic Disease': ['diabetes', 'obesity']
                    }
                },
                'DNA Processes': {
                    'keywords': ['DNA', 'repair', 'replication'],
                    'subcategories': {
                        'DNA Repair': ['repair', 'damage'],
                        'DNA Replication': ['replication', 'synthesis'],
                        'Chromatin': ['chromatin', 'histone']
                    }
                },
                'Cell Cycle': {
                    'keywords': ['cycle', 'division', 'mitotic'],
                    'subcategories': {
                        'Mitosis': ['mitotic', 'spindle'],
                        'Checkpoints': ['checkpoint', 'arrest'],
                        'Cell Division': ['cytokinesis', 'division']
                    }
                },
                'Metabolism': {
                    'keywords': ['metabolism', 'metabolic', 'biosynthesis'],
                    'subcategories': {
                        'Lipid Metabolism': ['lipid', 'fatty acid'],
                        'Protein Metabolism': ['protein', 'proteolysis'],
                        'Energy Metabolism': ['glycolysis', 'TCA']
                    }
                },
                'Immune System': {
                    'keywords': ['immune', 'cytokine', 'interferon', 'interleukin'],
                    'subcategories': {
                        'Innate Immunity': ['innate', 'toll', 'TLR', 'inflammasome'],
                        'Adaptive Immunity': ['adaptive', 'antibody', 'T cell', 'B cell'],
                        'Cytokine Signaling': ['cytokine', 'chemokine']
                    }
                },
                'Developmental Biology': {
                    'keywords': ['development', 'differentiation', 'morphogenesis'],
                    'subcategories': {
                        'Embryonic Development': ['embryo', 'embryonic'],
                        'Tissue Development': ['tissue', 'organogenesis'],
                        'Stem Cell': ['stem cell', 'pluripotent']
                    }
                }
            }
            
            # Create color scheme with distinct colors for better visibility
            # Using a more vibrant color palette
            color_palette = [
                '#FF5733', '#33FF57', '#3357FF', '#FF33A8', '#33A8FF', 
                '#A833FF', '#FFD733', '#33FFD7', '#D733FF', '#FF3333',
                '#33FFA8', '#A8FF33', '#3333FF', '#FF33D7', '#33D7FF'
            ]
            
            color_map = {}
            for i, cat in enumerate(categories.keys()):
                color_map[cat] = color_palette[i % len(color_palette)]
                
                # Create slightly different shades for subcategories
                subcats = len(categories[cat]['subcategories'])
                if subcats > 0:
                    base_color = np.array(matplotlib.colors.to_rgb(color_palette[i % len(color_palette)]))
                    for j, subcat in enumerate(categories[cat]['subcategories'].keys()):
                        # Create a lighter version of the base color
                        factor = 0.7 + (j % 3) * 0.1  # Vary the lightness
                        rgb = base_color * factor
                        rgb = np.clip(rgb, 0, 1)  # Ensure values are in [0,1]
                        hex_color = matplotlib.colors.rgb2hex(rgb)
                        color_map[subcat] = hex_color
            
            # Add 'Other' category
            color_map['Other'] = '#CCCCCC'
            
            # Assign node colors and categories
            node_categories = {}
            for node in G.nodes():
                assigned = False
                # Check main categories first
                for category, props in categories.items():
                    if any(kw.lower() in node.lower() for kw in props['keywords']):
                        node_categories[node] = category
                        assigned = True
                        break
                    # Check subcategories
                    if not assigned:
                        for subcat, keywords in props['subcategories'].items():
                            if any(kw.lower() in node.lower() for kw in keywords):
                                node_categories[node] = f"{category} - {subcat}"
                                assigned = True
                                break
                if not assigned:
                    node_categories[node] = 'Other'
            
            # Smart abbreviation function
            def smart_abbreviate(name, max_words=5):
                words = name.split()
                if len(words) <= max_words:
                    return name
                return ' '.join(words[:max_words]) + ' ...'
            
            # Calculate node sizes
            max_gene_count = max(len(G.nodes[node]['genes']) for node in G.nodes())
            min_gene_count = min(len(G.nodes[node]['genes']) for node in G.nodes())
            
            # Get top 20 nodes by gene count
            top_20_nodes = sorted([(node, len(G.nodes[node]['genes'])) 
                                 for node in G.nodes()],
                                key=lambda x: x[1], reverse=True)[:20]
            top_20_set = set(node for node, _ in top_20_nodes)
            
            # Use a force-directed layout algorithm to position nodes
            pos = nx.spring_layout(G, k=0.3, iterations=50, seed=42)
            
            # Prepare data for Plotly
            edge_x = []
            edge_y = []
            edge_trace = []
            
            # Create curved edges
            for edge in G.edges():
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                
                # Calculate control point for the curve (perpendicular to the line)
                mid_x = (x0 + x1) / 2
                mid_y = (y0 + y1) / 2
                
                # Calculate perpendicular direction
                dx = x1 - x0
                dy = y1 - y0
                length = math.sqrt(dx*dx + dy*dy)
                
                # Normalize and rotate 90 degrees
                if length > 0:
                    dx, dy = -dy/length, dx/length
                else:
                    dx, dy = 0, 0
                
                # Control point offset (adjust for curvature)
                control_x = mid_x + dx * 0.2
                control_y = mid_y + dy * 0.2
                
                # Create a Bezier curve with 20 points
                t = np.linspace(0, 1, 20)
                bezier_x = (1-t)**2 * x0 + 2*(1-t)*t * control_x + t**2 * x1
                bezier_y = (1-t)**2 * y0 + 2*(1-t)*t * control_y + t**2 * y1
                
                # Add the edge trace
                edge_trace.append(
                    go.Scatter(
                        x=bezier_x, y=bezier_y,
                        line=dict(width=1, color='rgba(150,150,150,0.8)'),
                        hoverinfo='none',
                        mode='lines',
                        showlegend=False
                    )
                )
                
                # Add arrow at the end of the edge
                arrow_x = [bezier_x[-2], bezier_x[-1], bezier_x[-2]]
                arrow_y = [bezier_y[-2], bezier_y[-1], bezier_y[-2]]
                
                # Calculate arrow head points
                arrow_length = 0.03
                dx = bezier_x[-1] - bezier_x[-2]
                dy = bezier_y[-1] - bezier_y[-2]
                length = math.sqrt(dx*dx + dy*dy)
                
                if length > 0:
                    dx, dy = dx/length, dy/length
                else:
                    dx, dy = 0, 0
                
                # Rotate 30 degrees clockwise and counter-clockwise
                arrow_x1 = bezier_x[-1] - arrow_length * (dx * math.cos(math.pi/6) - dy * math.sin(math.pi/6))
                arrow_y1 = bezier_y[-1] - arrow_length * (dx * math.sin(math.pi/6) + dy * math.cos(math.pi/6))
                
                arrow_x2 = bezier_x[-1] - arrow_length * (dx * math.cos(math.pi/6) + dy * math.sin(math.pi/6))
                arrow_y2 = bezier_y[-1] - arrow_length * (-dx * math.sin(math.pi/6) + dy * math.cos(math.pi/6))
                
                # Add arrow head
                edge_trace.append(
                    go.Scatter(
                        x=[arrow_x1, bezier_x[-1], arrow_x2],
                        y=[arrow_y1, bezier_y[-1], arrow_y2],
                        line=dict(width=1, color='rgba(150,150,150,0.8)'),
                        hoverinfo='none',
                        mode='lines',
                        showlegend=False
                    )
                )
            
            # Create node trace
            node_trace = go.Scatter(
                x=[pos[node][0] for node in G.nodes()],
                y=[pos[node][1] for node in G.nodes()],
                mode='markers+text',
                hoverinfo='text',
                marker=dict(
                    showscale=False,
                    colorscale='YlGnBu',
                    reversescale=True,
                    color=[],
                    size=[],
                    line=dict(width=2, color='black'),
                    opacity=0.8
                ),
                text=[],
                textposition="top center",
                textfont=dict(
                    family="Arial",
                    size=[],
                    color='black'
                )
            )
            
            # Add node properties
            node_colors = []
            node_sizes = []
            node_texts = []
            node_hovers = []
            font_sizes = []
            
            for node in G.nodes():
                gene_count = len(G.nodes[node]['genes'])
                category = node_categories[node]
                main_category = category.split(' - ')[0] if ' - ' in category else category
                
                # Get color from category
                color = color_map.get(main_category, '#CCCCCC')
                node_colors.append(color)
                
                # Calculate size based on gene count
                size = 15 + (gene_count - min_gene_count) / (max_gene_count - min_gene_count) * 40
                node_sizes.append(size)
                
                # Create node label
                node_texts.append(smart_abbreviate(node))
                
                # Create hover text
                hover_text = f"<b>{node}</b><br>Category: {category}<br>Genes: {gene_count}<br>"
                hover_text += ", ".join(G.nodes[node]['genes'][:10])
                if len(G.nodes[node]['genes']) > 10:
                    hover_text += f"<br>...and {len(G.nodes[node]['genes'])-10} more"
                node_hovers.append(hover_text)
                
                # Set font size based on importance
                font_size = 22 if node in top_20_set else 18
                font_sizes.append(font_size)
            
            node_trace.marker.color = node_colors
            node_trace.marker.size = node_sizes
            node_trace.text = node_texts
            node_trace.hovertext = node_hovers
            node_trace.textfont.size = font_sizes
            
            # Create figure with more space for legends
            fig = go.Figure(
                data=edge_trace + [node_trace],
                layout=go.Layout(
                    title='Pathway Impact Network',
                    titlefont=dict(size=24),
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=20, l=5, r=5, t=40),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    width=1400,  # Increased width to accommodate legends
                    height=900,
                    plot_bgcolor='white'
                )
            )
            
            # Create a separate legend for categories
            # Position the legend on the right side with more space
            legend_x = 1.05  # Move further right
            legend_y = 1.0
            
            # Add title for category legend with larger font
            fig.add_annotation(
                x=legend_x,
                y=legend_y + 0.02,  # Position slightly above
                xref="paper",
                yref="paper",
                text="<b>Pathway Categories</b>",  # Bold text
                showarrow=False,
                font=dict(size=18, color="black"),  # Larger font
                align="left"
            )
            
            # Add category items with more spacing and larger color boxes
            for i, (category, color) in enumerate([(cat, color_map[cat]) for cat in categories.keys()]):
                # Add colored rectangle
                fig.add_shape(
                    type="rect",
                    x0=legend_x,
                    y0=legend_y - 0.05 - i * 0.05,  # More vertical spacing
                    x1=legend_x + 0.03,  # Larger rectangle
                    y1=legend_y - 0.02 - i * 0.05,
                    fillcolor=color,
                    line=dict(color="black", width=1),
                    xref="paper",
                    yref="paper"
                )
                
                # Add category text
                fig.add_annotation(
                    x=legend_x + 0.04,  # Position text further from rectangle
                    y=legend_y - 0.035 - i * 0.05,
                    xref="paper",
                    yref="paper",
                    text=category,
                    showarrow=False,
                    font=dict(size=14, color="black"),  # Larger font
                    align="left"
                )
            
            # Add size legend with more spacing
            # Position below the category legend
            size_legend_y = legend_y - 0.05 - len(categories) * 0.05 - 0.1
            
            # Add title for size legend with larger font
            fig.add_annotation(
                x=legend_x,
                y=size_legend_y + 0.02,  # Position slightly above
                xref="paper",
                yref="paper",
                text="<b>Node Sizes</b>",  # Bold text
                showarrow=False,
                font=dict(size=18, color="black"),  # Larger font
                align="left"
            )
            
            # Add size legend items with more spacing and larger circles
            for i, count in enumerate([5, 10, 15, 20]):
                # Calculate size using the same formula as for nodes
                size = 15 + (count - min_gene_count) / (max_gene_count - min_gene_count) * 40
                scale_factor = 1.5  # Make circles larger in the legend
                
                # Add circle
                fig.add_shape(
                    type="circle",
                    x0=legend_x + 0.015 - (size * scale_factor)/1000,
                    y0=size_legend_y - 0.05 - i * 0.07 - (size * scale_factor)/1000,  # More vertical spacing
                    x1=legend_x + 0.015 + (size * scale_factor)/1000,
                    y1=size_legend_y - 0.05 - i * 0.07 + (size * scale_factor)/1000,
                    fillcolor="#4287f5",  # Use a consistent color for size legend
                    line=dict(color="black", width=1),
                    xref="paper",
                    yref="paper"
                )
                
                # Add size text
                fig.add_annotation(
                    x=legend_x + 0.04,  # Position text further from circle
                    y=size_legend_y - 0.05 - i * 0.07,
                    xref="paper",
                    yref="paper",
                    text=f"{count} genes",
                    showarrow=False,
                    font=dict(size=14, color="black"),  # Larger font
                    align="left"
                )
            
            # Add annotation for top 20 pathways
            fig.add_annotation(
                x=0.5,
                y=0.02,
                xref="paper",
                yref="paper",
                text="Top 20 pathways are highlighted with larger font",
                showarrow=False,
                font=dict(size=14, color="black"),
                align="center"
            )
            
            # Save figure with higher quality
            config = {
                'toImageButtonOptions': {
                    'format': 'png',
                    'filename': 'pathway_network',
                    'height': 900,
                    'width': 1400,
                    'scale': 2  # Higher resolution
                }
            }
            
            pio.write_html(fig, file=output_file, auto_open=False, config=config)
            
            # Save detailed analysis
            txt_file = output_file.rsplit('.', 1)[0] + '_details.txt'
            with open(txt_file, 'w') as f:
                f.write("Pathway Analysis Details\n" + "="*50 + "\n\n")
                
                # Write summary statistics
                f.write(f"Total pathways analyzed: {len(significant_pathways)}\n")
                f.write(f"Total genes affecting pathways: {len(set.union(*[set(genes) for genes in significant_pathways.values()]))}\n\n")
                
                # Group by categories and write detailed information
                for category, props in categories.items():
                    f.write(f"\n{category}\n" + "-"*50 + "\n")
                    category_pathways = []
                    
                    # Find pathways in this category
                    for pathway, genes in significant_pathways.items():
                        if node_categories.get(pathway, '').startswith(category):
                            category_pathways.append((pathway, genes))
                    
                    # Sort by gene count and write details
                    for pathway, genes in sorted(category_pathways, key=lambda x: len(x[1]), reverse=True):
                        f.write(f"  {pathway} ({len(genes)} genes):\n")
                        f.write(f"    Genes: {', '.join(genes)}\n")
                        f.write(f"    Category: {node_categories.get(pathway, 'Other')}\n\n")
            
            print(f"\nStatic visualization saved as: {output_file}")
            print(f"Detailed results saved as: {txt_file}")
            return True
            
        except Exception as e:
            print(f"Error in visualization: {e}")
            import traceback
            traceback.print_exc()
            return None

def main():
    parser = argparse.ArgumentParser(description='Create static pathway impact visualization')
    parser.add_argument('--genes', help='Comma-separated list of genes')
    parser.add_argument('--gene-file', help='File containing list of genes (one per line)')
    parser.add_argument('--output', default="output/pathway_impact_static_v2.html", help='Output file name')
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = ReactomePathwayAnalyzer()
    
    # Process gene list
    genes = []
    if args.genes:
        genes = [g.strip() for g in args.genes.split(',')]
    elif args.gene_file:
        with open(args.gene_file, 'r') as f:
            genes = [line.strip() for line in f if line.strip()]
    else:
        print("Error: Either --genes or --gene-file must be provided")
        return
    
    print(f"\nAnalyzing pathways for {len(genes)} genes...")
    
    # Get affected pathways
    affected_pathways = analyzer.get_affected_pathways(genes)
    
    # Create visualization
    network = analyzer.create_static_pathway_network(affected_pathways, args.output)
    
    # Print summary
    print("\nAnalysis Summary:")
    print("Number of input genes: {}".format(len(genes)))
    print("Number of affected pathways: {}".format(len(affected_pathways)))
    print("\nTop affected pathways (by number of genes):")
    for pathway, affected_genes in sorted(affected_pathways.items(), 
                               key=lambda x: len(x[1]), reverse=True)[:10]:
        print("- {}: {} genes".format(pathway, len(affected_genes)))
        print("  Genes: {}".format(", ".join(affected_genes)))
    
    if network:
        print("\nVisualization saved to: {}".format(args.output))
    else:
        print("\nNo visualization was created due to lack of pathway data.")

if __name__ == "__main__":
    main() 