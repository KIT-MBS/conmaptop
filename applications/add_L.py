import pandas as pd


def main() -> None:
    df = pd.read_csv("MolToFam_Validation_Set.csv", comment="#")
    df["From_One"] = True
    df["From_One"] = df["From_One"].astype("bool")
    df["Add_Tail"] = 0
    df["Add_Tail"] = df["Add_Tail"].astype("int")
    df["L"] = 0
    df["L"].astype("int")
    for e in df.iloc:
        with open(f"Sequence_Validation_python Set/{e.Filename}.fa", "r") as file:
            lines = file.readlines()
            df.loc[df["Filename"] == e.Filename, "L"] = len(lines[1])
    df.to_csv("MolToFam_Validation_Set_withL.csv", index=False)


if __name__ == "__main__":
    main()
