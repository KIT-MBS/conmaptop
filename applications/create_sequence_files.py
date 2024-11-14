import os
import pandas as pd


def main() -> None:
    df = pd.read_csv("MolToFam_Validation_Set.csv", comment="#")
    for i, row in df.iterrows():
        with open(f"MSA_Validation_Set/{row.Family}.faclean", "r") as msa_file:
            msa_file.readline()
            with open(f"Sequence_Validation_Set/{row.PDB}.fa", "w") as seq_file:
                seq_file.write(msa_file.readline().strip())


if __name__ == "__main__":
    main()
