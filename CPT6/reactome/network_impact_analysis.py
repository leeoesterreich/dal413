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

    def create_pathway_network(self, affected_pathways, output_file="output/pathway_impact.html"):
        """Create an interactive network visualization using vis.js"""
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
                }
            }
            
            # Create color scheme
            main_colors = plt.cm.Set3(np.linspace(0, 1, len(categories)))
            color_map = {}
            for i, (cat, props) in enumerate(categories.items()):
                # Convert RGB to hex
                rgb = main_colors[i][:3]  # Get RGB values (ignore alpha)
                hex_color = "#{:02x}{:02x}{:02x}".format(
                    int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
                color_map[cat] = hex_color
                
                # Create slightly different shades for subcategories
                subcats = len(props['subcategories'])
                if subcats > 0:
                    sub_colors = plt.cm.Pastel1(np.linspace(0, 1, subcats))
                    for j, subcat in enumerate(props['subcategories'].keys()):
                        rgb = sub_colors[j][:3]
                        hex_color = "#{:02x}{:02x}{:02x}".format(
                            int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
                        color_map[subcat] = hex_color
            
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
            
            # Prepare data for vis.js
            nodes_data = []
            edges_data = []
            
            # Get top 20 nodes by gene count
            top_20_nodes = sorted([(node, len(G.nodes[node]['genes'])) 
                                 for node in G.nodes()],
                                key=lambda x: x[1], reverse=True)[:20]
            top_20_set = set(node for node, _ in top_20_nodes)
            
            # Add nodes
            for node in G.nodes():
                gene_count = len(G.nodes[node]['genes'])
                category = node_categories[node]
                main_category = category.split(' - ')[0] if ' - ' in category else category
                
                # Calculate size based on gene count (scaled for vis.js) - EVEN SMALLER NODES
                size_factor = 15 + (gene_count - min_gene_count) / (max_gene_count - min_gene_count) * 40
                
                # Create node data
                node_data = {
                    'id': node,
                    'label': f"{smart_abbreviate(node)}\n({gene_count} genes)",
                    'title': f"<b>{node}</b><br>Genes: {', '.join(G.nodes[node]['genes'][:10])}" + 
                            (f"<br>...and {len(G.nodes[node]['genes'])-10} more" if len(G.nodes[node]['genes']) > 10 else ""),
                    'value': gene_count,  # Used for size
                    'color': {
                        'background': color_map.get(main_category, color_map['Signal Transduction']),
                        'border': '#000000',
                        'highlight': {
                            'background': color_map.get(main_category, color_map['Signal Transduction']),
                            'border': '#2B7CE9'
                        }
                    },
                    'font': {
                        'size': 22 if node in top_20_set else 18,  # EVEN BIGGER FONTS
                        'face': 'arial',
                        'color': '#000000',
                        'bold': node in top_20_set
                    },
                    'category': category,
                    'level': 1 if node in top_20_set else 2  # Set hierarchy level - top 20 at higher level
                }
                nodes_data.append(node_data)
            
            # Add edges
            for source, target in G.edges():
                edges_data.append({
                    'from': source,
                    'to': target,
                    'arrows': {
                        'to': {
                            'enabled': True,
                            'scaleFactor': 1.5  # LONGER ARROWS
                        }
                    },
                    'color': {
                        'color': '#848484',
                        'opacity': 0.8
                    },
                    'width': 1.5,
                    'smooth': {
                        'type': 'curvedCW',
                        'roundness': 0.2,
                        'forceDirection': 'none'
                    }
                })
            
            # Create legend data
            legend_items = []
            for category, color in color_map.items():
                if category in categories:  # Only show main categories
                    legend_items.append({
                        'label': category,
                        'color': color
                    })
            
            # Create size legend data
            size_legend_items = []
            for count in [5, 10, 15, 20]:
                size_legend_items.append({
                    'label': f'{count} genes',
                    'size': 15 + (count - min_gene_count) / (max_gene_count - min_gene_count) * 40  # MATCH SMALLER NODE SIZES
                })
            
            # Create HTML template with vis.js
            html_template = """
            <!DOCTYPE html>
            <html>
            <head>
                <title>Interactive Pathway Network</title>
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
                        font-size: 14px;
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
                        font-size: 12px;
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
                    button {
                        margin-right: 5px;
                        padding: 5px 10px;
                        cursor: pointer;
                    }
                    #search-container {
                        position: absolute;
                        top: 10px;
                        left: 10px;
                        background-color: rgba(255, 255, 255, 0.9);
                        border: 1px solid #ccc;
                        border-radius: 5px;
                        padding: 10px;
                        z-index: 1000;
                    }
                    #search-input {
                        width: 200px;
                        padding: 5px;
                        margin-right: 5px;
                    }
                    #title {
                        position: absolute;
                        top: 10px;
                        left: 50%;
                        transform: translateX(-50%);
                        background-color: rgba(255, 255, 255, 0.9);
                        border: 1px solid #ccc;
                        border-radius: 5px;
                        padding: 10px;
                        z-index: 1000;
                        text-align: center;
                    }
                    #title h1 {
                        margin: 0;
                        font-size: 18px;
                    }
                    #title p {
                        margin: 5px 0 0 0;
                        font-size: 12px;
                    }
                </style>
            </head>
            <body>
                <div id="network-container"></div>
                
                <div id="title">
                    <h1>Interactive Pathway Impact Network</h1>
                    <p>Drag nodes to reposition. Hover for details. Scroll to zoom.</p>
                </div>
                
                <div id="search-container">
                    <input type="text" id="search-input" placeholder="Search pathways...">
                    <button id="search-button">Search</button>
                </div>
                
                <div id="legend">
                    <div class="legend-section">
                        <div class="legend-title">Categories</div>
                        <div id="category-legend"></div>
                    </div>
                    <div class="legend-section">
                        <div class="legend-title">Node Sizes</div>
                        <div id="size-legend"></div>
                    </div>
                </div>
                
                <div id="controls">
                    <button id="save-positions">Save Layout</button>
                    <button id="load-positions">Load Layout</button>
                    <button id="reset-positions">Reset Layout</button>
                    <button id="export-svg">Export as PNG</button>
                </div>
                
                <script type="text/javascript">
                    // Network data
                    const nodes = new vis.DataSet(NODES_DATA);
                    const edges = new vis.DataSet(EDGES_DATA);
                    
                    // Network options
                    const options = {
                        nodes: {
                            shape: 'circle',
                            scaling: {
                                min: 6,  /* EVEN SMALLER MIN SIZE */
                                max: 55, /* EVEN SMALLER MAX SIZE */
                                label: {
                                    enabled: true,
                                    min: 18,  /* EVEN BIGGER MIN FONT */
                                    max: 32   /* EVEN BIGGER MAX FONT */
                                }
                            },
                            font: {
                                multi: true
                            }
                        },
                        edges: {
                            smooth: {
                                type: 'curvedCW',
                                roundness: 0.2,
                                forceDirection: 'none'
                            },
                            arrows: {
                                to: {
                                    enabled: true,
                                    scaleFactor: 1.5
                                }
                            }
                        },
                        physics: {
                            stabilization: {
                                iterations: 200
                            },
                            barnesHut: {
                                gravitationalConstant: -2000,
                                centralGravity: 0.1,
                                springLength: 150,
                                springConstant: 0.05,
                                damping: 0.09
                            },
                            hierarchicalRepulsion: {
                                nodeDistance: 150
                            }
                        },
                        layout: {
                            improvedLayout: true,
                            hierarchical: {
                                enabled: false,
                                levelSeparation: 150,
                                direction: 'UD',
                                sortMethod: 'directed'
                            }
                        },
                        interaction: {
                            navigationButtons: true,
                            keyboard: true,
                            tooltipDelay: 200,
                            hover: true
                        }
                    };
                    
                    // Create network
                    const container = document.getElementById('network-container');
                    const data = {
                        nodes: nodes,
                        edges: edges
                    };
                    const network = new vis.Network(container, data, options);
                    
                    // Create category legend
                    const categoryLegend = document.getElementById('category-legend');
                    const legendItems = LEGEND_ITEMS;
                    
                    legendItems.forEach(item => {
                        const legendItem = document.createElement('div');
                        legendItem.className = 'legend-item';
                        
                        const colorBox = document.createElement('div');
                        colorBox.className = 'legend-color';
                        colorBox.style.backgroundColor = item.color;
                        
                        const label = document.createElement('div');
                        label.className = 'legend-label';
                        label.textContent = item.label;
                        
                        legendItem.appendChild(colorBox);
                        legendItem.appendChild(label);
                        categoryLegend.appendChild(legendItem);
                    });
                    
                    // Create size legend
                    const sizeLegend = document.getElementById('size-legend');
                    const sizeLegendItems = SIZE_LEGEND_ITEMS;
                    
                    sizeLegendItems.forEach(item => {
                        const legendItem = document.createElement('div');
                        legendItem.className = 'legend-item';
                        
                        const sizeCircle = document.createElement('div');
                        sizeCircle.className = 'legend-size';
                        const size = Math.sqrt(item.size) * 1.5;
                        sizeCircle.style.width = size + 'px';
                        sizeCircle.style.height = size + 'px';
                        
                        const label = document.createElement('div');
                        label.className = 'legend-label';
                        label.textContent = item.label;
                        
                        legendItem.appendChild(sizeCircle);
                        legendItem.appendChild(label);
                        sizeLegend.appendChild(legendItem);
                    });
                    
                    // Search functionality
                    document.getElementById('search-button').addEventListener('click', function() {
                        const searchTerm = document.getElementById('search-input').value.toLowerCase();
                        if (searchTerm) {
                            const nodeIds = nodes.getIds();
                            const matchingNodes = nodeIds.filter(id => 
                                id.toLowerCase().includes(searchTerm)
                            );
                            
                            if (matchingNodes.length > 0) {
                                network.selectNodes(matchingNodes);
                                network.focus(matchingNodes[0], {
                                    scale: 1.2,
                                    animation: true
                                });
                            } else {
                                alert('No matching pathways found');
                            }
                        }
                    });
                    
                    // Enter key for search
                    document.getElementById('search-input').addEventListener('keyup', function(event) {
                        if (event.key === 'Enter') {
                            document.getElementById('search-button').click();
                        }
                    });
                    
                    // Save positions
                    document.getElementById('save-positions').addEventListener('click', function() {
                        const positions = network.getPositions();
                        localStorage.setItem('pathwayNetworkPositions', JSON.stringify(positions));
                        alert('Layout saved');
                    });
                    
                    // Load positions
                    document.getElementById('load-positions').addEventListener('click', function() {
                        const savedPositions = localStorage.getItem('pathwayNetworkPositions');
                        if (savedPositions) {
                            const positions = JSON.parse(savedPositions);
                            network.setOptions({ physics: false });
                            
                            // Update positions for existing nodes
                            Object.keys(positions).forEach(nodeId => {
                                if (nodes.get(nodeId)) {
                                    network.moveNode(nodeId, positions[nodeId].x, positions[nodeId].y);
                                }
                            });
                            
                            alert('Layout loaded');
                        } else {
                            alert('No saved layout found');
                        }
                    });
                    
                    // Reset positions
                    document.getElementById('reset-positions').addEventListener('click', function() {
                        network.setOptions({ physics: { enabled: true } });
                        network.stabilize();
                    });
                    
                    // Export as SVG
                    document.getElementById('export-svg').addEventListener('click', function() {
                        // Get network canvas
                        const canvas = network.canvas.frame.canvas;
                        
                        // Create a temporary canvas to draw the network
                        const tempCanvas = document.createElement('canvas');
                        tempCanvas.width = canvas.width;
                        tempCanvas.height = canvas.height;
                        const ctx = tempCanvas.getContext('2d');
                        
                        // Draw white background
                        ctx.fillStyle = 'white';
                        ctx.fillRect(0, 0, canvas.width, canvas.height);
                        
                        // Draw the network canvas on top
                        ctx.drawImage(canvas, 0, 0);
                        
                        // Convert to image
                        const imageURL = tempCanvas.toDataURL('image/png');
                        
                        // Create download link
                        const downloadLink = document.createElement('a');
                        downloadLink.href = imageURL;
                        downloadLink.download = 'pathway_network.png';
                        document.body.appendChild(downloadLink);
                        downloadLink.click();
                        document.body.removeChild(downloadLink);
                    });
                </script>
            </body>
            </html>
            """
            
            # Replace placeholders with actual data
            html_content = html_template.replace('NODES_DATA', json.dumps(nodes_data))
            html_content = html_content.replace('EDGES_DATA', json.dumps(edges_data))
            html_content = html_content.replace('LEGEND_ITEMS', json.dumps(legend_items))
            html_content = html_content.replace('SIZE_LEGEND_ITEMS', json.dumps(size_legend_items))
            
            # Save HTML file
            with open(output_file, 'w') as f:
                f.write(html_content)
            
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
                    for subcategory, keywords in props['subcategories'].items():
                        f.write(f"\n{subcategory}:\n")
                        for pathway, genes in significant_pathways.items():
                            if any(kw.lower() in pathway.lower() for kw in keywords):
                                f.write(f"  {pathway} ({len(genes)} genes):\n")
                                f.write(f"    Genes: {', '.join(genes)}\n")
                                f.write(f"    Category: {node_categories.get(pathway, 'Other')}\n")
            
            print(f"\nInteractive visualization saved as: {output_file}")
            print(f"Detailed results saved as: {txt_file}")
            return True
            
        except Exception as e:
            print(f"Error in visualization: {e}")
            import traceback
            traceback.print_exc()
            return None

def main():
    parser = argparse.ArgumentParser(description='Analyze pathway impact of somatic mutations')
    parser.add_argument('--genes', help='Comma-separated list of genes')
    parser.add_argument('--gene-file', help='File containing list of genes (one per line)')
    parser.add_argument('--output', default="output/pathway_impact.html", help='Output file name')
    
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
