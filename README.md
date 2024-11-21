# cmsDS
Interact-tools for Stereo-seq in Jupyter


# Installation

Create the conda environment

```bash
conda env create -f environment.yml
```

# Usage

Here we use SAW default output h5ad.  ([Stereo-seq V1.3 Mouse Tongue DEMO data](https://en.stomics.tech/resources/demo-data/list.html))

Please check that your h5ad anndata has gone through `scanpy.pp.calculate_qc_metrics(adata)`

```python
# Stereo-seq V1.3, Mouse Tongue DEMO data
adata = sc.read_h5ad("B04001D2.cellbin_1.0.h5ad")

# Initialization
interactor = InteractH5ad(adata)
```

## Analyzing

### Data preprocessing

![preprocess](https://github.com/user-attachments/assets/6858e3f1-7516-4dc4-9bc2-c2b4df1b2562)

### Reducing dimensions

![reduceDimension](https://github.com/user-attachments/assets/b33a26a9-31e1-4034-9458-949686b68f4d)

### Clustering

![cluster](https://github.com/user-attachments/assets/9d8f0ec7-d00b-4952-b49c-38799b1e70f1)

### Exporting DEGs

![getDEG](https://github.com/user-attachments/assets/15e8bfce-48a2-4968-be53-e70f7d0a66e1)

### Extracting data

You can extract the cells or genes you wanted, for further analyzing.

Please note that you should check the cell query command yourself.

- Filter by `.obs.columns` or list of `genes`
  
![operateData_1](https://github.com/user-attachments/assets/58a0a713-3a99-49f0-a0a9-45ba0ed77d4b)

- Extract subclusters
  
![operateData_2](https://github.com/user-attachments/assets/87e09945-765c-4c31-9404-317cf16be171)

## Visualization

### Spatial

![plot_spatial](https://github.com/user-attachments/assets/b4244180-60e8-4b65-bb1d-1a96d4b28c69)

### UMAP

![plot_umap](https://github.com/user-attachments/assets/4e1794df-24f0-4d7a-af31-dabf3934a35b)

### Violin

![plot_violin](https://github.com/user-attachments/assets/8f14493b-41ed-4f9e-860b-e1998597434b)

### Dotplot

![plot_dotplot](https://github.com/user-attachments/assets/32b6f006-3f7d-46ea-acc2-90fbbae7cfa3)

# Save your data in .h5ad

```python
final_h5ad = interactor.fun_output() # anndata in scanpy flavor

# Save by:
# final_h5ad("path_to_your_file.h5ad")
```
