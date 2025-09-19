# Updated Migration Structure

## New Repository Structure

The migration has been updated to organize everything under the `docs/` folder:

```
leeoesterreich/CITEgeist/
├── docs/
│   ├── core_scripts/           # Renamed CITEgeist code
│   │   ├── model/              # Python modules
│   │   └── Jupyter/            # Jupyter notebooks
│   ├── source/                 # Sphinx documentation source
│   │   ├── notebooks/          # Vignette RST files
│   │   ├── api.rst             # API documentation
│   │   └── ...                 # Other documentation
│   ├── build/                  # Built documentation
│   └── (HTML files)            # Final GitHub Pages files
```

## Updated Files

### Configuration Files
- ✅ `conf.py` - Updated paths: `../` instead of `../../`
- ✅ `requirements.txt` - Binder dependencies
- ✅ `postBuild` - Binder setup script

### Vignette Files (All Updated)
- ✅ `vignette_1_biopsy_heterogeneity.rst`
- ✅ `vignette_2_surgical_d538g.rst` 
- ✅ `vignette_3_responder_macrophages.rst`
- ✅ `vignette_4_external_tools.rst`
- ✅ `notebooks.rst` - Main vignettes page

### Binder Links Updated
All Binder/Colab links now point to:
- Repository: `leeoesterreich/CITEgeist`
- Branch: `dev`
- Path: `docs/core_scripts/Jupyter/`

### Migration Commands
- ✅ Updated `MIGRATION_COMMANDS.md` with new folder structure
- ✅ All copy commands point to `docs/core_scripts/`
- ✅ Verification commands updated

## Key Benefits

1. **Clean Organization**: Everything documentation-related is under `docs/`
2. **Proper Sphinx Paths**: Configuration correctly points to `docs/core_scripts/`
3. **Working Binder Links**: All interactive notebooks accessible via correct paths
4. **GitHub Pages Ready**: Structure supports proper deployment

## Next Steps

1. Follow the updated `MIGRATION_COMMANDS.md`
2. Use the `PULL_REQUEST_GUIDE.md` 
3. Configure GitHub Pages to serve from `dev` branch `/docs` folder
4. New documentation URL will be: `https://leeoesterreich.github.io/CITEgeist/`

All files are ready for migration! 🚀

