
import scanpy as sc

# constants
FILE_DATA = "/home/javaprog/Data/Broad/GeneticsML/Gears202512/Data/norman/perturb_processed.h5ad"

map_hda5 = {
  'norman': '/home/javaprog/Data/Broad/GeneticsML/Gears202512/Data/norman/perturb_processed.h5ad',
  'obsity_orig': '/home/javaprog/Data/Broad/GeneticsML/Gears202512/Data/broad_obesity/obesity_challenge_1.h5ad',
  'obesity_gears': '/home/javaprog/Data/Broad/GeneticsML/Gears202512/Data/broad_obesity/perturb_processed.h5ad'
}

for key, value in map_hda5.items():
  print("reading: {}".format(value))
  adata = sc.read_h5ad(value)

  print("dataset: {}, columns: {}".format(key, adata.obs.columns))
  print(adata.obs.head(10))

  if "condition" in adata.obs.columns:
      print("Found 'condition' column. Values:")
      print(adata.obs["condition"].unique())
  else:
      print("No 'condition' column found.")

  if "cell_type" in adata.obs.columns:
      print("Found 'cell_type' column. Values:")
      print(adata.obs["cell_type"].unique())
  else:
      print("No 'cell_type' column found.")

  if "control" in adata.obs.columns:
      print("Found 'control' column. Values:")
      print(adata.obs["control"].value_counts().head(20))
  else:
      print("No 'control' column found.")
  
  print()


