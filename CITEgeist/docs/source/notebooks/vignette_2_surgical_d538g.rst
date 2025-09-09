Vignette 2: ESR1 D538G Mutations in Surgical Samples
====================================================

This vignette demonstrates the analysis of spatial transcriptomics data from surgical samples containing ESR1 D538G mutations, a common mutation in breast cancer. The analysis combines CITEgeist deconvolution with pathway analysis and mutation signature validation.

Overview
--------

This vignette covers:

* Analysis of surgical specimens with known mutations
* Integration of mutation data with spatial transcriptomics
* Pathway analysis and signature validation
* Comparison of mutation-positive and mutation-negative regions
* Advanced visualization techniques

Key Learning Objectives
-----------------------

After completing this notebook, you will understand:

* How to analyze surgical specimens with CITEgeist
* Methods for integrating mutation data with spatial analysis
* Approaches to pathway analysis in spatial context
* Techniques for validating mutation signatures spatially

Download and Run
----------------

.. raw:: html

    <div class="admonition note">
    <p class="admonition-title">Interactive Notebook</p>
    <p>Download this notebook to run it locally with your own data.</p>
    </div>

**Download the notebook**: `vignette_2_surgical_d538g.ipynb <vignette_2_surgical_d538g.ipynb>`_

**Prerequisites**:
* CITEgeist installed (see :doc:`../installation`)
* Gurobi license configured
* Surgical sample spatial transcriptomics data
* Mutation annotation data
* Required Python packages: scanpy, pandas, matplotlib, seaborn

**Data Requirements**:
* Spatial transcriptomics data from surgical specimens
* Mutation annotation (ESR1 D538G status)
* Gene expression and antibody capture data
* Spatial coordinates and tissue morphology

Code Example
------------

Here's a preview of the key steps from the notebook:

.. code-block:: python

   # Load surgical sample data
   surgical_data_path = os.path.join(DATA_FOLDER, "surgical_sample/outs")
   
   # Initialize model for surgical analysis
   model = CitegeistModel(
       sample_name="surgical_d538g",
       adata=surgical_adata,
       output_folder="./surgical_results"
   )
   
   # Process with lower min_count for surgical samples
   model.filter_gex(min_counts=5)  # Lower threshold for surgical samples
   model.preprocess_gex()
   model.preprocess_antibody()
   
   # Run deconvolution
   props, finetuned_props = model.run_cell_proportion_model(radius=75)

Expected Outputs
----------------

This notebook will generate:

* Cell type proportions in mutation-positive and negative regions
* Pathway enrichment analysis results
* Spatial maps of mutation-associated gene expression
* Validation plots for mutation signatures
* Comparative analysis between regions

Next Steps
----------

After completing this vignette:

1. Apply the analysis to your own surgical samples
2. Explore different mutation types and their spatial patterns
3. Integrate with clinical outcome data
4. Try the macrophage analysis in :doc:`vignette_3_responder_macrophages`

Related Documentation
---------------------

* :doc:`../tutorial` - Comprehensive workflow guide
* :doc:`vignette_1_biopsy_heterogeneity` - Basic analysis workflow
* :doc:`../api` - Complete API reference
