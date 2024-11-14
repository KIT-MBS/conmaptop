import pandas as pd
import os.path
from BioHelpers_FABER.cmap import Cmap
from itertools import chain


def main() -> None:
    df = pd.read_csv("MolToFam_Validation_Set.csv", comment="#")
    df["Unmodelled_Regions"] = ""
    df["Unmodelled_Regions"] = df["Unmodelled_Regions"].astype("object")
    for i, e in df.iterrows():
        if os.path.isfile(f"PDB_Validation_Set/{e.Filename}.pdb"):
            cm = Cmap()
            cm.load_native_pdb(
                f"PDB_Validation_Set/{e.Filename}.pdb",
                id=e.PDB,
                start_from_one=e.From_One,
                added_tail=e.Add_Tail,
            )
            print(
                f"Loaded {e.Filename}.pdb - unmodelled regions {cm.get_unmodelled_regions()} - L={cm.l}"
            )
            df.at[i, "Unmodelled_Regions"] = cm.get_unmodelled_regions()
        else:
            print(f"Couldn't find {e.Filename}.pdb")
    df.to_csv("MolToFam_Validation_Set_regions.csv", index=False)


if __name__ == "__main__":
    main()
