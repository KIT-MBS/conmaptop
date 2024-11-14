import pandas as pd


def main():
    df = pd.read_csv("MolToFam_Validation_Set.csv", comment="#")
    df_m = pd.read_csv("effective_validation_set.csv")

    df = pd.merge(df_m, df, on="Family")
    df.to_csv("MolToFam_Validation_Set.csv", index=False)


if __name__ == "__main__":
    main()
