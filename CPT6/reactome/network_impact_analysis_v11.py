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

class ReactomePathwayAnalyzer:
    def __init__(self, gene_pathway_file="data/NCBI2Reactome_All_Levels.txt", 
                 pathway_hierarchy_file="data/ReactomePathwaysRelation.txt",
                 gene_name_file=None):
        """Initialize with local Reactome data files"""
        print("Loading Reactome data files...")
        
        # Load only essential columns with optimized dtypes
        dtype_dict = {
            'Gene': str,
            'Pathway_ID': str,
            'Pathway_Name': str,
            'Species': str
        }
        
        # Load gene-pathway associations
        try:
            self.gene_pathway_df = pd.read_csv(gene_pathway_file, sep='\t', header=None, 
                                              names=['Gene', 'Pathway_ID', 'URL', 'Pathway_Name', 'Evidence_Code', 'Species'],
                                              usecols=['Gene', 'Pathway_ID', 'Pathway_Name', 'Species'],
                                              dtype=dtype_dict)
            
            # Filter for human pathways
            self.gene_pathway_df = self.gene_pathway_df[self.gene_pathway_df['Species'].str.contains('Homo sapiens')]
            print(f"Loaded {len(self.gene_pathway_df)} gene-pathway associations")
        except Exception as e:
            print(f"Error loading gene-pathway file: {e}")
            self.gene_pathway_df = pd.DataFrame(columns=['Gene', 'Pathway_ID', 'Pathway_Name', 'Species'])
        
        # Load pathway hierarchy
        try:
            hierarchy_df = pd.read_csv(pathway_hierarchy_file, sep='\t', header=None,
                                      names=['Parent_ID', 'Child_ID'])
            
            # Create a dictionary for faster lookup
            self.pathway_hierarchy = defaultdict(list)
            for _, row in hierarchy_df.iterrows():
                self.pathway_hierarchy[row['Parent_ID']].append(row['Child_ID'])
            
            print(f"Loaded {len(hierarchy_df)} pathway relationships")
        except Exception as e:
            print(f"Error loading pathway hierarchy file: {e}")
            self.pathway_hierarchy = defaultdict(list)
        
        # Load gene name to NCBI ID mapping if provided
        self.gene_name_map = {}
        if gene_name_file:
            try:
                with open(gene_name_file, 'r') as f:
                    for line in f:
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            gene_name = parts[0]
                            ncbi_id = parts[1]
                            self.gene_name_map[ncbi_id] = gene_name
                print(f"Loaded {len(self.gene_name_map)} gene name mappings")
            except Exception as e:
                print(f"Error loading gene name file: {e}")
    
    def get_affected_pathways(self, genes):
        """Get pathways affected by the input genes"""
        affected_pathways = defaultdict(list)
        
        # Convert all genes to strings for comparison
        genes = [str(g) for g in genes]
        
        print(f"\nAnalyzing {len(genes)} genes...")
        
        # Process genes in batches for better performance
        batch_size = 50
        for i in range(0, len(genes), batch_size):
            batch = genes[i:i+batch_size]
            # Filter dataframe for current batch of genes
            batch_df = self.gene_pathway_df[self.gene_pathway_df['Gene'].isin(batch)]
            
            # Group by pathway and collect genes
            for _, row in batch_df.iterrows():
                affected_pathways[row['Pathway_Name']].append(row['Gene'])
            
            print(f"Processed {min(i+batch_size, len(genes))}/{len(genes)} genes", end='\r')
        
        return affected_pathways
    
    def get_pathway_hierarchy(self, pathway_names):
        """Get hierarchical relationships between pathways"""
        edges = []
        
        # Create a mapping from pathway name to ID
        pathway_ids = {}
        for _, row in self.gene_pathway_df.iterrows():
            if row['Pathway_Name'] in pathway_names:
                pathway_ids[row['Pathway_ID']] = row['Pathway_Name']
        
        # Check relationships using the pre-built hierarchy
        for parent_id, children in self.pathway_hierarchy.items():
            if parent_id in pathway_ids:
                parent_name = pathway_ids[parent_id]
                for child_id in children:
                    if child_id in pathway_ids:
                        child_name = pathway_ids[child_id]
                        edges.append((parent_name, child_name))
        
        return edges

    def create_pathway_network(self, affected_pathways, output_file="output/pathway_impact.html"):
        """Create an interactive network visualization using vis.js"""
        if not affected_pathways:
            print("No pathways found to visualize.")
            return None
        
        try:
            print("Starting visualization creation...")
            # Ensure output is HTML
            if not output_file.endswith('.html'):
                output_file = os.path.splitext(output_file)[0] + '.html'
            
            # Filter to keep only the most significant pathways (top 100 by gene count)
            significant_pathways = dict(sorted(affected_pathways.items(), 
                                            key=lambda x: len(x[1]), 
                                            reverse=True)[:100])
            print(f"Selected {len(significant_pathways)} significant pathways")
            
            # Create a directed graph
            G = nx.DiGraph()
            
            # Add nodes with properties
            for pathway, genes in significant_pathways.items():
                # Map NCBI IDs to gene names if available
                gene_names = []
                for gene_id in genes:
                    if gene_id in self.gene_name_map:
                        gene_names.append(f"{self.gene_name_map[gene_id]} ({gene_id})")
                    else:
                        gene_names.append(gene_id)
                
                G.add_node(pathway, size=len(genes), genes=gene_names)
            print("Added nodes to graph")
            
            # Add hierarchical relationships
            edges = self.get_pathway_hierarchy(list(significant_pathways.keys()))
            G.add_edges_from(edges)
            print(f"Added {len(edges)} edges to graph")
            
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
                'Gene Expression': {
                    'keywords': ['expression', 'transcription', 'translation'],
                    'subcategories': {
                        'Transcription': ['transcription', 'RNA polymerase'],
                        'Translation': ['translation', 'ribosome'],
                        'RNA Processing': ['splicing', 'mRNA']
                    }
                },
                'Developmental Biology': {
                    'keywords': ['development', 'differentiation', 'morphogenesis'],
                    'subcategories': {
                        'Neuronal Development': ['neuronal', 'axon', 'dendrite'],
                        'Embryonic Development': ['embryonic', 'embryo'],
                        'Tissue Development': ['tissue', 'organogenesis']
                    }
                }
            }
            
            # Create color scheme
            main_colors = plt.cm.Set3(np.linspace(0, 1, len(categories)))
            color_map = {}
            
            for i, (category, _) in enumerate(categories.items()):
                # Convert RGB to hex
                rgb = main_colors[i][:3]  # Extract RGB values
                hex_color = '#{:02x}{:02x}{:02x}'.format(int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
                color_map[category] = hex_color
            
            # Assign categories and colors to nodes
            for node in G.nodes():
                # Default category and color
                G.nodes[node]['category'] = 'Other'
                G.nodes[node]['color'] = '#CCCCCC'
                
                # Check if node matches any category
                for category, info in categories.items():
                    keywords = info['keywords']
                    if any(keyword.lower() in node.lower() for keyword in keywords):
                        G.nodes[node]['category'] = category
                        G.nodes[node]['color'] = color_map[category]
                        break
            print("Assigned categories and colors to nodes")
            
            # Prepare data for vis.js
            nodes_data = []
            min_gene_count = min(G.nodes[node]['size'] for node in G.nodes())
            max_gene_count = max(G.nodes[node]['size'] for node in G.nodes())
            print(f"Node size range: {min_gene_count} to {max_gene_count} genes")
            
            # Scale factor for node sizes - make them more proportional to gene count
            scale_factor = 40  # Increased from the original value
            
            for node in G.nodes():
                size = G.nodes[node]['size']
                # More proportional scaling
                scaled_size = 15 + (size - min_gene_count) / (max_gene_count - min_gene_count) * scale_factor
                
                # Create a clean tooltip without HTML tags
                genes_list = G.nodes[node]['genes']
                tooltip_text = f"{node}\nCategory: {G.nodes[node]['category']}\nGenes: {size}\n" + "\n".join(genes_list)
                
                nodes_data.append({
                    'id': node,
                    'label': f"{node}\n({size} genes)",
                    'title': tooltip_text,
                    'value': size,  # This will be used for sizing
                    'color': G.nodes[node]['color'],
                    'shape': 'dot',
                    'font': {
                        'size': 14,  # Increased from 12
                        'face': 'Arial',
                        'color': 'black'
                    },
                    'size': scaled_size  # Use the scaled size
                })
            print("Prepared node data for vis.js")
            
            edges_data = []
            for source, target in G.edges():
                edges_data.append({
                    'from': source,
                    'to': target,
                    'arrows': 'to',
                    'width': 2.5,  # Increased arrow width for bolder arrows
                    'color': {
                        'color': '#848484',
                        'opacity': 0.8
                    },
                    'smooth': {
                        'type': 'curvedCW',
                        'roundness': 0.2
                    }
                })
            print("Prepared edge data for vis.js")
            
            # Create vis.js network options
            options = json.dumps({
                'nodes': {
                    'scaling': {
                        'min': 15,
                        'max': 55,  # Increased max size
                        'label': {
                            'enabled': True,
                            'min': 14,  # Increased from 12
                            'max': 20   # Increased from 18
                        }
                    },
                    'font': {
                        'size': 14,  # Increased from 12
                        'face': 'Arial'
                    }
                },
                'edges': {
                    'arrows': {
                        'to': {
                            'enabled': True,
                            'scaleFactor': 2.0  # Increased for bolder arrows
                        }
                    },
                    'color': {
                        'color': '#848484',
                        'opacity': 0.8
                    },
                    'width': 2.5,  # Increased width for bolder arrows
                    'smooth': {
                        'type': 'curvedCW',
                        'roundness': 0.2,
                        'forceDirection': 'none'
                    }
                },
                'physics': {
                    'barnesHut': {
                        'gravitationalConstant': -2000,
                        'centralGravity': 0.1,
                        'springLength': 150,
                        'springConstant': 0.05,
                        'damping': 0.09
                    },
                    'maxVelocity': 50,
                    'minVelocity': 0.1,
                    'solver': 'barnesHut',
                    'stabilization': {
                        'enabled': True,
                        'iterations': 1000,
                        'updateInterval': 100
                    },
                    'timestep': 0.5,
                    'adaptiveTimestep': True
                },
                'interaction': {
                    'hover': True,
                    'navigationButtons': True,
                    'keyboard': {
                        'enabled': True
                    }
                }
            })
            print("Created network options")
            
            # Create legend data
            legend_items = []
            for category, color in color_map.items():
                if category in [G.nodes[node]['category'] for node in G.nodes()]:  # Only show categories that are used
                    legend_items.append({
                        'label': category,
                        'color': color
                    })
            
            # Add "Other" category to legend if there are any uncategorized nodes
            if 'Other' in [G.nodes[node]['category'] for node in G.nodes()]:
                legend_items.append({
                    'label': 'Other',
                    'color': '#CCCCCC'  # Default gray color for "Other" category
                })
            print("Created category legend data")
            
            # Create size legend data - more representative sizes
            size_legend_items = []
            size_steps = [min_gene_count]
            if max_gene_count > min_gene_count:
                step = (max_gene_count - min_gene_count) / 3
                for i in range(1, 4):
                    size_steps.append(min_gene_count + i * step)
            
            for count in size_steps:
                scaled_size = 15 + (count - min_gene_count) / (max_gene_count - min_gene_count) * scale_factor
                size_legend_items.append({
                    'label': f'{int(count)} genes',
                    'size': scaled_size
                })
            print("Created size legend data")
            
            # Create HTML template with vis.js
            html_template = """
            <!DOCTYPE html>
            <html>
            <head>
                <title>Reactome representation of CPT6</title>
                <meta charset="utf-8">
                <script type="text/javascript" src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>
                <style type="text/css">
                    body, html {
                        margin: 0;
                        padding: 0;
                        font-family: Arial, sans-serif;
                        height: 100%;
                        overflow: hidden;
                    }
                    #network-container {
                        width: 100%;
                        height: 100%;
                        border: 1px solid lightgray;
                    }
                    #legend {
                        position: absolute;
                        top: 10px;
                        right: 10px;
                        background-color: rgba(255, 255, 255, 0.9);
                        border: 1px solid #ccc;
                        border-radius: 5px;
                        padding: 10px;
                        z-index: 1000;
                        max-width: 250px;
                    }
                    .legend-title {
                        font-weight: bold;
                        margin-bottom: 5px;
                        font-size: 16px;  /* Increased from 14 */
                    }
                    .legend-item {
                        display: flex;
                        align-items: center;
                        margin-bottom: 5px;
                    }
                    .legend-color {
                        width: 15px;
                        height: 15px;
                        margin-right: 5px;
                        border: 1px solid #000;
                    }
                    .legend-size {
                        border-radius: 50%;
                        background-color: #97C2FC;
                        border: 1px solid #2B7CE9;
                        margin-right: 5px;
                    }
                    .legend-label {
                        font-size: 14px;  /* Increased from 12 */
                    }
                    .legend-section {
                        margin-bottom: 15px;
                    }
                    #controls {
                        position: absolute;
                        bottom: 10px;
                        left: 10px;
                        background-color: rgba(255, 255, 255, 0.9);
                        border: 1px solid #ccc;
                        border-radius: 5px;
                        padding: 10px;
                        z-index: 1000;
                    }
                    .control-button {
                        margin: 5px;
                        padding: 5px 10px;
                        background-color: #f8f8f8;
                        border: 1px solid #ddd;
                        border-radius: 3px;
                        cursor: pointer;
                        font-size: 14px;  /* Increased font size */
                    }
                    .control-button:hover {
                        background-color: #e8e8e8;
                    }
                    .title-bar {
                        position: absolute;
                        top: 0;
                        left: 0;
                        right: 0;
                        background-color: rgba(255, 255, 255, 0.9);
                        border-bottom: 1px solid #ccc;
                        padding: 10px;
                        text-align: center;
                        font-size: 20px;  /* Increased from 18 */
                        font-weight: bold;
                        z-index: 1000;
                    }
                    .subtitle {
                        font-size: 14px;  /* Increased from 12 */
                        font-weight: normal;
                        margin-top: 5px;
                    }
                </style>
            </head>
            <body>
                <div class="title-bar">
                    Reactome representation of CPT6
                    <div class="subtitle">
                        Node size represents number of affected genes<br>
                        Colors indicate pathway categories<br>
                        Edges show hierarchical relationships
                    </div>
                </div>
                <div id="network-container"></div>
                <div id="legend">
                    <div class="legend-section">
                        <div class="legend-title">Pathway Categories</div>
                        CATEGORY_LEGEND
                    </div>
                    <div class="legend-section">
                        <div class="legend-title">Node Size (Gene Count)</div>
                        SIZE_LEGEND
                    </div>
                </div>
                <div id="controls">
                    <button class="control-button" onclick="network.fit()">Fit View</button>
                    <button class="control-button" onclick="togglePhysics()">Toggle Physics</button>
                    <button class="control-button" onclick="saveImage()">Save Image</button>
                </div>
                
                <script type="text/javascript">
                    // Create network
                    var container = document.getElementById('network-container');
                    var data = {
                        nodes: new vis.DataSet(NODES_DATA),
                        edges: new vis.DataSet(EDGES_DATA)
                    };
                    var options = OPTIONS;
                    var network = new vis.Network(container, data, options);
                    
                    // Toggle physics
                    var physicsEnabled = true;
                    function togglePhysics() {
                        physicsEnabled = !physicsEnabled;
                        network.setOptions({physics: {enabled: physicsEnabled}});
                    }
                    
                    // Save image
                    function saveImage() {
                        var canvas = network.canvas.frame.canvas;
                        var link = document.createElement('a');
                        link.href = canvas.toDataURL('image/png');
                        link.download = 'pathway_network.png';
                        link.click();
                    }
                    
                    // Stabilize and fit
                    network.once('stabilizationIterationsDone', function() {
                        network.fit();
                    });
                    
                    // Override the default tooltip behavior to show plain text
                    network.on("hoverNode", function(params) {
                        var nodeId = params.node;
                        var node = data.nodes.get(nodeId);
                        var tooltip = document.createElement('div');
                        tooltip.innerHTML = node.title.replace(/\\n/g, '<br>');
                        tooltip.style.position = 'absolute';
                        tooltip.style.backgroundColor = 'white';
                        tooltip.style.padding = '10px';
                        tooltip.style.border = '1px solid #ccc';
                        tooltip.style.borderRadius = '5px';
                        tooltip.style.boxShadow = '0 0 10px rgba(0,0,0,0.2)';
                        tooltip.style.fontSize = '14px';
                        tooltip.style.zIndex = '1001';
                        tooltip.style.maxWidth = '300px';
                        tooltip.style.whiteSpace = 'pre-wrap';
                        document.body.appendChild(tooltip);
                        
                        // Position the tooltip near the mouse
                        tooltip.style.left = (params.event.center.x + 10) + 'px';
                        tooltip.style.top = (params.event.center.y + 10) + 'px';
                        
                        // Remove the tooltip when mouse leaves the node
                        network.once("blurNode", function() {
                            document.body.removeChild(tooltip);
                        });
                    });
                </script>
            </body>
            </html>
            """
            print("Created HTML template")
            
            # Create category legend HTML
            category_legend_html = ""
            for item in legend_items:
                category_legend_html += f"""
                <div class="legend-item">
                    <div class="legend-color" style="background-color: {item['color']};"></div>
                    <div class="legend-label">{item['label']}</div>
                </div>
                """
            
            # Create size legend HTML
            size_legend_html = ""
            for item in size_legend_items:
                size_legend_html += f"""
                <div class="legend-item">
                    <div class="legend-size" style="width: {item['size']}px; height: {item['size']}px;"></div>
                    <div class="legend-label">{item['label']}</div>
                </div>
                """
            print("Created legend HTML")
            
            # Replace placeholders in template
            html_content = html_template.replace('NODES_DATA', json.dumps(nodes_data))
            html_content = html_content.replace('EDGES_DATA', json.dumps(edges_data))
            html_content = html_content.replace('OPTIONS', options)
            html_content = html_content.replace('CATEGORY_LEGEND', category_legend_html)
            html_content = html_content.replace('SIZE_LEGEND', size_legend_html)
            print("Replaced placeholders in template")
            
            # Create output directory if it doesn't exist
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            
            # Write HTML file
            with open(output_file, 'w') as f:
                f.write(html_content)
            print(f"Wrote HTML file to {output_file}")
            
            # Save detailed results to a text file
            details_file = os.path.splitext(output_file)[0] + '_details.txt'
            with open(details_file, 'w') as f:
                f.write("Pathway Impact Analysis Results\n")
                f.write("==============================\n\n")
                f.write(f"Total input genes: {len(set(gene for genes in affected_pathways.values() for gene in genes))}\n")
                f.write(f"Total affected pathways: {len(affected_pathways)}\n\n")
                f.write("Top affected pathways (by number of genes):\n")
                for pathway, genes in sorted(affected_pathways.items(), key=lambda x: len(x[1]), reverse=True):
                    gene_names = []
                    for gene_id in genes:
                        if gene_id in self.gene_name_map:
                            gene_names.append(f"{self.gene_name_map[gene_id]} ({gene_id})")
                        else:
                            gene_names.append(gene_id)
                    f.write(f"- {pathway}: {len(genes)} genes\n")
                    f.write(f"  Genes: {', '.join(gene_names)}\n\n")
            print(f"Wrote details file to {details_file}")
            
            return G
        
        except Exception as e:
            print(f"Error creating visualization: {e}")
            import traceback
            traceback.print_exc()
            return None

def main():
    parser = argparse.ArgumentParser(description='Analyze pathway impact of somatic mutations')
    parser.add_argument('--genes', help='Comma-separated list of genes')
    parser.add_argument('--gene-file', help='File containing list of genes (one per line)')
    parser.add_argument('--gene-name-file', help='File mapping NCBI IDs to gene names')
    parser.add_argument('--output', default="output/pathway_impact.html", help='Output file name')
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = ReactomePathwayAnalyzer(gene_name_file=args.gene_name_file)
    
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
    network = analyzer.create_pathway_network(affected_pathways, args.output)
    
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