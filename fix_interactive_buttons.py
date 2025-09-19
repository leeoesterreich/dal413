#!/usr/bin/env python3
"""
Script to manually add interactive buttons to HTML files
since Sphinx build environment is broken
"""

import re
import os

def add_interactive_buttons_to_html(html_file_path, notebook_name):
    """Add Binder and Colab buttons to HTML file"""
    
    # Read the HTML file
    with open(html_file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Interactive buttons HTML
    buttons_html = f'''<div style="margin: 20px 0;">
    <a href="https://mybinder.org/v2/gh/leeoesterreich/dal413/HEAD?labpath=CITEgeist%2FCITEgeist%2FJupyter%2F{notebook_name}.ipynb" target="_blank">
        <img src="https://mybinder.org/badge_logo.svg" alt="Launch Binder" style="margin-right: 10px;">
    </a>
    <a href="https://colab.research.google.com/github/leeoesterreich/dal413/blob/main/CITEgeist/CITEgeist/Jupyter/{notebook_name}.ipynb" target="_blank">
        <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab">
    </a>
</div>

<div class="admonition tip">
<p class="admonition-title">🚀 Run Interactively</p>
<p><strong>Binder</strong>: Click above to launch this notebook in a live, executable environment with all dependencies pre-installed. No setup required!</p>
<p><strong>Colab</strong>: Run in Google Colab with free GPU access. You may need to install some packages.</p>
<p><strong>Download</strong>: <a href="{notebook_name}.ipynb">Download the notebook</a> to run locally.</p>
</div>'''
    
    # Find the section after the title and insert buttons
    # Look for the first <p> tag after the main heading
    pattern = r'(<h1>.*?</h1>.*?<p>.*?</p>)'
    
    if re.search(pattern, content, re.DOTALL):
        # Insert buttons after the first paragraph
        content = re.sub(
            r'(<h1>.*?</h1>.*?<p>.*?</p>)',
            r'\1\n' + buttons_html + '\n',
            content,
            count=1,
            flags=re.DOTALL
        )
    else:
        # Fallback: insert after the first <h1> tag
        content = re.sub(
            r'(<h1>.*?</h1>)',
            r'\1\n' + buttons_html + '\n',
            content,
            count=1,
            flags=re.DOTALL
        )
    
    # Write back to file
    with open(html_file_path, 'w', encoding='utf-8') as f:
        f.write(content)
    
    print(f"✅ Added interactive buttons to {html_file_path}")

def main():
    """Add interactive buttons to all vignette HTML files"""
    
    base_dir = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/Github/dal413/docs/notebooks"
    
    vignettes = [
        "vignette_1_biopsy_heterogeneity",
        "vignette_2_surgical_d538g", 
        "vignette_3_responder_macrophages",
        "vignette_4_external_tools"
    ]
    
    for vignette in vignettes:
        html_file = os.path.join(base_dir, f"{vignette}.html")
        if os.path.exists(html_file):
            add_interactive_buttons_to_html(html_file, vignette)
        else:
            print(f"❌ File not found: {html_file}")

if __name__ == "__main__":
    main()

