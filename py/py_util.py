import pandas as pd


def align_on_index(df1: pd.DataFrame, df2: pd.DataFrame):
    # sécurise le type + trim
    df1 = df1.copy(); df2 = df2.copy()
    df1.index = df1.index.astype(str).str.strip()
    df2.index = df2.index.astype(str).str.strip()

    # intersection en conservant l'ORDRE de df1
    common = df1.index[df1.index.isin(df2.index)]
    if len(common) == 0:
        raise ValueError("Aucun ID commun entre df1 et df2.")

    df1_al = df1.loc[common]
    df2_al = df2.loc[common]
    assert (df1_al.index == df2_al.index).all()
    print(f"{len(common)} IDs communs conservés (ordre synchronisé).")
    return df1_al, df2_al

