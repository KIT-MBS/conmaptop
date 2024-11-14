import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import ast

OFFDIAGONAL_THRESHOLD = 4

PATH_SUFFIX = "_Validation_Set"
# PATH_SUFFIX = ""
CSV_FILE = f"MolToFam{PATH_SUFFIX}.csv"


def write_res_file(filename: str, contacts: list):
    with open(filename, "w") as f:
        for con in contacts:
            l = " A/" + str(con[0] + 1) + "/N A/" + str(con[1] + 1) + "/N "
            f.write("WELL" + l + "3.5 9.5 1.0 \n")
            f.write("SLOPE" + l + "3.5 9.5 1.0 \n")
            f.write("SLOPE" + l + "0 25 -1.0 \n")


def is_in_modelled_region(i: int, unmodelled_regions: list[tuple[int, int]]) -> bool:
    for reg in unmodelled_regions:
        if reg[0] <= i + 1 <= reg[1]:
            return False
    return True


# CAVEAT: DCA uses site labelling starting by 0
def write_DCA(df: pd.DataFrame, noc_fact: float = 1.0) -> None:
    for e in df.iloc:
        noc = math.floor(e.L * noc_fact)
        contacts_df = pd.read_csv(f"Restraints{PATH_SUFFIX}/DCA/{e.Family}_raw.txt")
        contacts_list = [
            (int(f["i"]), int(f[" j"]))
            for f in contacts_df.iloc
            if (abs(f["i"] - f[" j"]) > OFFDIAGONAL_THRESHOLD)
            and (
                is_in_modelled_region(
                    int(f["i"]), ast.literal_eval(e.Unmodelled_Regions)
                )
            )
            and (
                is_in_modelled_region(
                    int(f[" j"]), ast.literal_eval(e.Unmodelled_Regions)
                )
            )
        ][:noc]
        write_res_file(
            f"Restraints{PATH_SUFFIX}/DCA/simrna_{e.Filename}_L{noc_fact}.res",
            contacts_list,
        )


# CAVEAT: COCONET uses site labelling starting by 1
def write_COCONET(df: pd.DataFrame, noc_fact: float = 1.0) -> None:
    for e in df.iloc:
        noc = math.floor(e.L * noc_fact)
        contacts_df = pd.read_csv(
            f"Restraints{PATH_SUFFIX}/COCONET/COCONET3x3_{e.Family}.txt",
            names=["i", "j", "val"],
            sep="\t",
            comment="#",
        )
        contacts_list = [
            (int(f["i"]) - 1, int(f["j"]) - 1)
            for f in contacts_df.iloc
            if (abs(f["i"] - f["j"]) > OFFDIAGONAL_THRESHOLD)
            and (
                is_in_modelled_region(
                    int(f["i"]) - 1, ast.literal_eval(e.Unmodelled_Regions)
                )
            )
            and (
                is_in_modelled_region(
                    int(f["j"]) - 1, ast.literal_eval(e.Unmodelled_Regions)
                )
            )
        ][:noc]
        write_res_file(
            f"Restraints{PATH_SUFFIX}/COCONET/simrna_{e.Filename}_L{noc_fact}.res",
            contacts_list,
        )


# CAVEAT: Barnacle uses site labelling starting by 0
def write_Barnacle_vanilla(df: pd.DataFrame, noc_fact: float = 1.0) -> None:
    for e in df.iloc:
        noc = math.floor(e.L * noc_fact)
        contacts_list = read_nparray_to_list(
            f"Restraints{PATH_SUFFIX}/BARNACLE_VANILLA/{e.Family}.npy"
        )
        contacts_list = [
            c
            for c in contacts_list
            if is_in_modelled_region(c[0], ast.literal_eval(e.Unmodelled_Regions))
            and is_in_modelled_region(c[1], ast.literal_eval(e.Unmodelled_Regions))
        ]
        write_res_file(
            f"Restraints{PATH_SUFFIX}/BARNACLE_VANILLA/simrna_{e.Filename}_L{noc_fact}.res",
            contacts_list[:noc],
        )


# CAVEAT: Barnacle uses site labelling starting by 0
def write_Barnacle_gauss(df: pd.DataFrame, noc_fact: float = 1.0) -> None:
    for e in df.iloc:
        noc = math.floor(e.L * noc_fact)
        contacts_list = read_nparray_to_list(
            f"Restraints{PATH_SUFFIX}/BARNACLE_GAUSS/{e.Family}.npy"
        )
        contacts_list = [
            c
            for c in contacts_list
            if is_in_modelled_region(c[0], ast.literal_eval(e.Unmodelled_Regions))
            and is_in_modelled_region(c[1], ast.literal_eval(e.Unmodelled_Regions))
        ]
        write_res_file(
            f"Restraints{PATH_SUFFIX}/BARNACLE_GAUSS/simrna_{e.Filename}_L{noc_fact}.res",
            contacts_list[:noc],
        )


def read_nparray_to_list(filename: str) -> list[list]:
    matrix = np.load(filename)
    for i in range(-OFFDIAGONAL_THRESHOLD, OFFDIAGONAL_THRESHOLD + 1):
        ind = np.where(np.eye(matrix.shape[0], k=i, dtype=bool))
        matrix[ind] = 0
    matrix = np.triu(matrix)
    ind = np.unravel_index(np.argsort(matrix, axis=None), matrix.shape)
    out = [[i, j, matrix[i, j]] for i, j in zip(ind[0], ind[1])]
    return out[::-1]  # Reverse list to have the most important contact on top!


def main():
    df = pd.read_csv(f"MolToFam{PATH_SUFFIX}.csv", comment="#")
    noc_fact = 0.5
    # write_DCA(df, noc_fact)
    # write_COCONET(df, noc_fact)
    write_Barnacle_vanilla(df, noc_fact)
    write_Barnacle_gauss(df, noc_fact)
    print("Done")


if __name__ == "__main__":
    main()
