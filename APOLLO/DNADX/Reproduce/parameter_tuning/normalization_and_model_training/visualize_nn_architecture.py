import graphviz

def generate_xgboost_conceptual_diagram(input_feature_example_count=3, num_trees_to_show=3, filename="xgboost_conceptual_diagram"):
    """
    Generates a conceptual diagram of an XGBoost model structure.

    Args:
        input_feature_example_count (int): Number of example input features to display.
        num_trees_to_show (int): Number of example decision trees to display in the ensemble.
        filename (str): The name of the output file (without extension).
    """
    dot = graphviz.Digraph('XGBoostConcept', comment='Conceptual XGBoost Model for RNA Signature Prediction')
    dot.attr(rankdir='TB', splines='line', nodesep='0.5', ranksep='1.0', labeljust='c', labelloc='t')
    dot.graph_attr['fontsize'] = '18'
    dot.graph_attr['label'] = 'Conceptual XGBoost Model Structure'

    node_fontsize = '10'
    cluster_label_fontsize = '14'

    # Input Layer (CNA Segments)
    with dot.subgraph(name='cluster_input') as c:
        c.attr(label='Input Features\n(CNA Segments)', color='#E6E6FA', style='filled', fontsize=cluster_label_fontsize)
        # c.attr(rank='same') # Keep inputs at the top for TB rankdir
        for i in range(input_feature_example_count):
            c.node(f'feat_{i}', f'Feature {i+1}', shape='ellipse', style='filled', color='lightpink', fontsize=node_fontsize)
        if input_feature_example_count > 2:
             c.node('feat_dots', '...', shape='plaintext', fontsize='20')

    # Ensemble of Decision Trees
    with dot.subgraph(name='cluster_trees') as c:
        c.attr(label=f'Ensemble of Decision Trees\n(e.g., {num_trees_to_show} of N trees shown)', color='#FFFACD', style='filled', fontsize=cluster_label_fontsize)
        # c.attr(rank='same') # Trees are the next stage
        for i in range(num_trees_to_show):
            c.node(f'tree_{i}', f'Tree {i+1}', shape='component', style='filled', color='palegreen', fontsize=node_fontsize, peripheries='2')
        if num_trees_to_show > 1: # Or based on a general N_trees > num_trees_to_show
            c.node('tree_dots', '...', shape='plaintext', fontsize='20')

    # Aggregation Step
    with dot.subgraph(name='cluster_aggregation') as c:
        c.attr(label='Aggregation\n(e.g., Weighted Sum of Tree Predictions)', color='#ADD8E6', style='filled', fontsize=cluster_label_fontsize)
        # c.attr(rank='same')
        c.node('agg', 'Combine Predictions', shape='box', style='filled', color='lightblue', fontsize=node_fontsize)

    # Output Layer (RNA Signature Prediction)
    with dot.subgraph(name='cluster_output') as c:
        c.attr(label='Final Prediction\n(RNA Signature Score)', color='#FFDAB9', style='filled', fontsize=cluster_label_fontsize)
        # c.attr(rank='same')
        c.node('out_0', 'Predicted Score', shape='ellipse', style='filled', color='peachpuff', fontsize=node_fontsize)

    # Edges
    # Connect input features to the concept of the tree ensemble
    input_nodes_for_edges = [f'feat_{i}' for i in range(input_feature_example_count)]
    if input_feature_example_count > 2: input_nodes_for_edges.append('feat_dots')
    
    tree_nodes_for_edges = [f'tree_{i}' for i in range(num_trees_to_show)]
    if num_trees_to_show > 1: tree_nodes_for_edges.append('tree_dots')

    for in_node in input_nodes_for_edges:
        for tree_node_target in tree_nodes_for_edges: # Each feature conceptually goes to each tree processing block
             # To avoid clutter, connect inputs to the cluster or a representative point for the trees
             # For simplicity here, let's draw to each representative tree shown
             if in_node == 'feat_dots' and tree_node_target == 'tree_dots':
                 dot.edge(in_node, tree_node_target, style='dotted', arrowhead='vee')
             elif in_node == 'feat_dots' or tree_node_target == 'tree_dots':
                 dot.edge(in_node, tree_node_target, style='dotted', arrowhead='vee')
             else:
                 dot.edge(in_node, tree_node_target, arrowhead='vee')

    # Connect trees to aggregation
    for tree_node_source in tree_nodes_for_edges:
        if tree_node_source == 'tree_dots':
            dot.edge(tree_node_source, 'agg', style='dotted', arrowhead='vee')
        else:
            dot.edge(tree_node_source, 'agg', arrowhead='vee')
        
    # Connect aggregation to output
    dot.edge('agg', 'out_0', arrowhead='vee')

    try:
        dot.render(filename, view=False, format='png')
        print(f"XGBoost conceptual diagram saved as {filename}.png")
    except graphviz.backend.execute.ExecutableNotFound:
        print("ERROR: Graphviz executable not found. Please ensure Graphviz is installed and in your system's PATH.")
    except Exception as e:
        print(f"An error occurred during rendering: {e}")

if __name__ == '__main__':
    # You can adjust these parameters for the generated image
    generate_xgboost_conceptual_diagram(
        input_feature_example_count=3, 
        num_trees_to_show=3, 
        filename="xgboost_model_conceptual_v1"
    ) 