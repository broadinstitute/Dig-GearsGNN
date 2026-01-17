
import scanpy as sc

# constants
DIR_GEARS = "/home/javaprog/Data/Broad/GeneticsML/Gears202512"
FILE_DATA = "{}/Data/norman/perturb_processed.h5ad".format(DIR_GEARS)

map_hda5 = {
  'norman': '{}/Data/norman/perturb_processed.h5ad'.format(DIR_GEARS),
  'obsity_orig': '{}/Data/broad_obesity/obesity_challenge_1.h5ad'.format(DIR_GEARS),
  'obesity_gears': '{}/Data/broad_obesity/perturb_processed.h5ad'.format(DIR_GEARS)
}


def dict_schema(d, indent=0):
    for k, v in d.items():
        print(" " * indent + f"{k}: {type(v).__name__}")
        if isinstance(v, dict):
            dict_schema(v, indent + 2)

if __name__ == "__main__":

    for key, value in map_hda5.items():
        print("reading: {}".format(value))
        adata = sc.read_h5ad(value)
        print(adata)

        print("dataset: {}, columns: {}".format(key, adata.obs.columns))
        print(adata.obs.head(10))

        if "condition" in adata.obs.columns:
            print("Found 'condition' column. Values:")
            print(adata.obs["condition"].value_counts().head(20))
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
        
        # Top feature_call values
        if "feature_call" in adata.obs.columns:
            print("Found 'feature_call' column. Values:")
            print(adata.obs["feature_call"].value_counts().head(50))
        else:
            print("No 'feature_call' column found.")

        print()


