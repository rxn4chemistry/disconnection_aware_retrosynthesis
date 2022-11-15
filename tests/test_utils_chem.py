from collections import OrderedDict

from dar.chem import (
    AtomEnvironment,
    get_all_atom_indices,
    get_atom_list,
    get_atomic_neighbourhoods,
    remove_mapping,
    remove_rxn_mapping,
    remove_unmapped_components,
    standardise_reaction_component,
)


def test_remove_mapping():
    assert (
        remove_mapping(
            "CN(C)C=O.I[CH:28]([CH3:29])[CH3:30].[CH3:1][O:2][c:3]1[cH:4][c:5]([CH3:6])"
            "[c:7]([CH3:8])[cH:9][c:10]1[NH:11][C:12](=[O:13])[N:14]1[CH2:15][CH2:16]"
            "[N:17]([c:18]2[cH:19][c:20]([F:21])[cH:22][c:23]([F:24])[cH:25]2)[CH2:26]"
            "[CH2:27]1.[H-].[Na+]"
        )
        == "CC(C)I.CN(C)C=O.COc1cc(C)c(C)cc1NC(=O)N1CCN(c2cc(F)cc(F)c2)CC1.[H-].[Na+]"
    )
    assert (
        remove_mapping(
            "CC(C)(C)OC(=O)[N:18]1[CH2:17][CH2:16][c:14]2[c:13]([cH:12][n:11][c:10]"
            "([NH:9][C:7]([c:6]3[cH:5][cH:4][c:3]([CH2:2][CH3:1])[cH:29][cH:28]3)"
            "=[O:8])[n:15]2)[CH2:27]1.CCN(CC)CC.ClCCl.Cl[S:19](=[O:20])(=[O:21])"
            "[c:22]1[cH:23][cH:24][s:25][cH:26]1.O=C(O)C(F)(F)F"
        )
        == "CCN(CC)CC.CCc1ccc(C(=O)Nc2ncc3c(n2)CCN(C(=O)OC(C)(C)C)C3)cc1.ClCCl."
        "O=C(O)C(F)(F)F.O=S(=O)(Cl)c1ccsc1"
    )
    assert (
        remove_mapping(
            "C[OH:10].N#[C:8][c:7]1[cH:6][cH:5][cH:4][c:3]([CH2:2][CH3:1])[cH:11]1."
            "[Na+].[OH-:9]"
        )
        == "CCc1cccc(C#N)c1.CO.[Na+].[OH-]"
    )
    assert (
        remove_mapping(
            "C1CCOC1.CO.COC(C)(C)C.Cl[Ni]Cl.[BH4-].[CH3:1][O:2][C:3](=[O:4])[CH:5]="
            "[CH:6][c:7]1[cH:8][cH:9][c:10]([CH2:11][O:12][c:13]2[cH:14][c:15]"
            "([F:16])[cH:17][cH:18][c:19]2[F:20])[cH:21][c:22]1[OH:23].[Na+]"
        )
        == "C1CCOC1.CO.COC(=O)C=Cc1ccc(COc2cc(F)ccc2F)cc1O.COC(C)(C)C.Cl[Ni]Cl."
        "[BH4-].[Na+]"
    )


def test_remove_rxn_mapping():
    test_reactions = {
        "CCN(CC)CC.CCOCC.Cl[S:3]([CH2:2][CH3:1])(=[O:4])=[O:5].[OH:6][CH2:7][CH2:8]"
        "[Br:9]>>[CH3:1][CH2:2][S:3](=[O:4])(=[O:5])[O:6][CH2:7][CH2:8][Br:9]": "CCN"
        "(CC)CC.CCOCC.[CH3][CH2]S(=O)(=O)Cl.[OH][CH2][CH2]Br>>[CH3][CH2]S(=O)(=O)O"
        "[CH2][CH2]Br",
        "C(N([CH2:2][CH3:3])[CH2:6][CH3:5])[CH3:4].Cl[CH2:8]Cl.O=[C:7](Cl)[O:15]"
        "[c:14]1[cH:9][cH:10][cH:11][cH:12][cH:13]1.[Na+].[OH-].[OH2:1]>>[OH:1]"
        "[c:2]1[cH:3][cH:4][cH:5][cH:6][c:7]1[CH2:8][c:9]1[cH:10][cH:11][cH:12]"
        "[cH:13][c:14]1[OH:15]": "Cl[CH2]Cl.O.O=C(Cl)Oc1[cH][cH][cH][cH][cH]1."
        "[CH3]CN([CH2][CH3])[CH2][CH3].[Na+].[OH-]>>[OH]c1[cH][cH][cH][cH]c1[CH2]"
        "c1[cH][cH][cH][cH]c1[OH]",
        "N[c:2]1[cH:3][cH:4][c:5]2[n:6][c:7]3[n:8]([c:9]2[cH:10]1)[CH2:11][CH2:12]"
        "[CH2:13]3.O.O=N[O-].[ClH:1].[Na+]>>[Cl:1][c:2]1[cH:3][cH:4][c:5]2[n:6]"
        "[c:7]3[n:8]([c:9]2[cH:10]1)[CH2:11][CH2:12][CH2:13]3": "Cl.Nc1[cH][cH]"
        "c2nc3n(c2[cH]1)[CH2][CH2][CH2]3.O.O=N[O-].[Na+]>>Clc1[cH][cH]c2nc3n(c2"
        "[cH]1)[CH2][CH2][CH2]3",
        "C1[CH2:14][CH:15]1[C:7]([NH:6][NH:5][C:2]([CH3:1])([CH3:3])[CH3:4])("
        "[C:8]#[N:9])[CH3:10].CC(=O)C1CC1.O=C1CC[CH2:13][CH2:12][CH2:11]C1>>"
        "[CH3:1][C:2]([CH3:3])([CH3:4])[NH:5][NH:6][C:7]1([C:8]#[N:9])[CH2:10]"
        "[CH2:11][CH2:12][CH2:13][CH2:14][CH2:15]1": "CC(=O)C1CC1.O=C1CC[CH2]"
        "[CH2][CH2]C1.[CH3]C([CH3])([CH3])[NH][NH]C([CH3])(C#N)[CH]1C[CH2]1>>"
        "[CH3]C([CH3])([CH3])[NH][NH]C1(C#N)[CH2][CH2][CH2][CH2][CH2][CH2]1",
        "Cl[C:9](=[O:10])[CH:11]1[CH2:12][CH2:13][CH2:14]1.[O:1]=[C:2]([OH:3])"
        "[c:4]1[cH:5][cH:6][c:7]([NH2:8])[cH:15][c:16]1[N+:17](=[O:18])[O-:19]>>"
        "[O:1]=[C:2]([OH:3])[c:4]1[cH:5][cH:6][c:7]([NH:8][C:9](=[O:10])[CH:11]2"
        "[CH2:12][CH2:13][CH2:14]2)[cH:15][c:16]1[N+:17](=[O:18])[O-:19]": "O=C"
        "(Cl)[CH]1[CH2][CH2][CH2]1.[NH2]c1[cH][cH]c(C(=O)[OH])c([N+](=O)[O-])[cH]"
        "1>>O=C([OH])c1[cH][cH]c([NH]C(=O)[CH]2[CH2][CH2][CH2]2)[cH]c1[N+](=O)[O-]",
    }

    for mapped_rxn, unmapped_rxn in test_reactions.items():
        assert remove_rxn_mapping(mapped_rxn) == unmapped_rxn


def test_standardise_reaction_compnent():
    assert (
        standardise_reaction_component(
            "CO.[Ni].O[N:6]=[C:5]1[CH:4]=[C:3]([CH2:2][CH3:1])[C:9](=[O:10])[CH:8]="
            "[CH:7]1"
        )
        == "CO.O[N:6]=[C:5]1[CH:4]=[C:3]([CH2:2][CH3:1])[C:9](=[O:10])[CH:8]=[CH:7]1."
        "[Ni]"
    )
    assert (
        standardise_reaction_component(
            "O=[C:2]([OH:1])[c:3]1[cH:4][c:5]([Br:6])[cH:7][cH:8][c:9]1[Cl:10].B.Cl."
            "C1CCOC1"
        )
        == "B.C1CCOC1.Cl.O=[C:2]([OH:1])[c:3]1[cH:4][c:5]([Br:6])[cH:7][cH:8][c:9]1"
        "[Cl:10]"
    )
    assert (
        standardise_reaction_component(
            "O.[NH2:15][NH2:16].CCO.CCO[C:13]([c:12]1[cH:11][cH:10][c:9]([O:8][CH2:7]"
            "[c:6]2[cH:5][cH:4][c:3]([O:2][CH3:1])[cH:20][cH:19]2)[cH:18][cH:17]1)="
            "[O:14]"
        )
        == "CCO.CCO[C:13]([c:12]1[cH:11][cH:10][c:9]([O:8][CH2:7][c:6]2[cH:5][cH:4]"
        "[c:3]([O:2][CH3:1])[cH:20][cH:19]2)[cH:18][cH:17]1)=[O:14].O.[NH2:15][NH2:16]"
    )


def test_remove_unmapped_components():
    assert (
        remove_unmapped_components(
            "C1CCC2=NCCCN2CC1.C1CCOC1.CO.O.O[CH2:2][c:3]1[cH:4][c:5]([Cl:6])[cH:7]"
            "[cH:8][c:9]1[CH:10]1[CH2:11][CH2:12][N:13]1[C:14](=[O:15])[O:16][CH2:17]"
            "[C:18]([Cl:19])([Cl:20])[Cl:21].[N-]=[N+]=[N:1]P(=O)(c1ccccc1)c1ccccc1."
            "[Na+].[OH-].c1ccc(P(c2ccccc2)c2ccccc2)cc1"
        )
        == "O[CH2:2][c:3]1[cH:4][c:5]([Cl:6])[cH:7][cH:8][c:9]1[CH:10]1[CH2:11]"
        "[CH2:12][N:13]1[C:14](=[O:15])[O:16][CH2:17][C:18]([Cl:19])([Cl:20])[Cl:21]."
        "[N-]=[N+]=[N:1]P(=O)(c1ccccc1)c1ccccc1"
    )
    assert (
        remove_unmapped_components(
            "CO.O=C(OCc1ccccc1)[NH:9][C@@H:7]([C:5]([NH:4][CH:2]([CH3:1])[CH3:3])="
            "[O:6])[CH3:8].[H][H].[Pd]"
        )
        == "O=C(OCc1ccccc1)[NH:9][C@@H:7]([C:5]([NH:4][CH:2]([CH3:1])[CH3:3])=[O:6])"
        "[CH3:8]"
    )
    assert (
        remove_unmapped_components(
            "O=C(Cl)c1ccccc1F.O=[C:8](Cl)[c:7]1[c:2]([CH3:1])[cH:3][cH:4][cH:5][cH:6]"
            "1.[NH2:9][c:10]1[c:11]([Cl:12])[cH:13][cH:14][cH:15][c:16]1[C:17](="
            "[O:18])[OH:19]"
        )
        == "O=[C:8](Cl)[c:7]1[c:2]([CH3:1])[cH:3][cH:4][cH:5][cH:6]1.[NH2:9][c:10]1"
        "[c:11]([Cl:12])[cH:13][cH:14][cH:15][c:16]1[C:17](=[O:18])[OH:19]"
    )


def test_get_atomic_neighbourhoods():
    example_1 = (
        "C1CCOC1.CO.COC(C)(C)C.Cl[Ni]Cl.[BH4-].[CH3:1][O:2][C:3](=[O:4])[CH:5]"
        "=[CH:6][c:7]1[cH:8][cH:9][c:10]([CH2:11][O:12][c:13]2[cH:14][c:15]([F:16])[cH:17]"
        "[cH:18][c:19]2[F:20])[cH:21][c:22]1[OH:23].[Na+]"
    )
    answer_1 = OrderedDict(
        [
            (1, ["1_2_SINGLE"]),
            (2, ["1_2_SINGLE", "2_3_SINGLE"]),
            (3, ["2_3_SINGLE", "3_4_DOUBLE", "3_5_SINGLE"]),
            (4, ["3_4_DOUBLE"]),
            (5, ["3_5_SINGLE", "5_6_DOUBLE"]),
            (6, ["5_6_DOUBLE", "6_7_SINGLE"]),
            (7, ["6_7_SINGLE", "7_22_AROMATIC", "7_8_AROMATIC"]),
            (8, ["7_8_AROMATIC", "8_9_AROMATIC"]),
            (9, ["8_9_AROMATIC", "9_10_AROMATIC"]),
            (10, ["10_11_SINGLE", "10_21_AROMATIC", "9_10_AROMATIC"]),
            (11, ["10_11_SINGLE", "11_12_SINGLE"]),
            (12, ["11_12_SINGLE", "12_13_SINGLE"]),
            (13, ["12_13_SINGLE", "13_14_AROMATIC", "13_19_AROMATIC"]),
            (14, ["13_14_AROMATIC", "14_15_AROMATIC"]),
            (15, ["14_15_AROMATIC", "15_16_SINGLE", "15_17_AROMATIC"]),
            (16, ["15_16_SINGLE"]),
            (17, ["15_17_AROMATIC", "17_18_AROMATIC"]),
            (18, ["17_18_AROMATIC", "18_19_AROMATIC"]),
            (19, ["13_19_AROMATIC", "18_19_AROMATIC", "19_20_SINGLE"]),
            (20, ["19_20_SINGLE"]),
            (21, ["10_21_AROMATIC", "21_22_AROMATIC"]),
            (22, ["21_22_AROMATIC", "22_23_SINGLE", "7_22_AROMATIC"]),
            (23, ["22_23_SINGLE"]),
        ]
    )

    example_2 = (
        "C1COCCO1.CC(C)(C)[S@@](=O)[NH:1][C@H:2]1[CH2:3][CH2:4][CH2:5][N:6]"
        "([C:7](=[O:8])[O:9][CH2:10][c:11]2[cH:12][cH:13][cH:14][cH:15][cH:16]2)[CH2:17]"
        "[CH2:18]1.CO.Cl"
    )
    answer_2 = OrderedDict(
        [
            (1, ["0_1_SINGLE", "1_2_SINGLE"]),
            (2, ["1_2_SINGLE", "2_18_SINGLE", "2_3_SINGLE"]),
            (3, ["2_3_SINGLE", "3_4_SINGLE"]),
            (4, ["3_4_SINGLE", "4_5_SINGLE"]),
            (5, ["4_5_SINGLE", "5_6_SINGLE"]),
            (6, ["5_6_SINGLE", "6_17_SINGLE", "6_7_SINGLE"]),
            (7, ["6_7_SINGLE", "7_8_DOUBLE", "7_9_SINGLE"]),
            (8, ["7_8_DOUBLE"]),
            (9, ["7_9_SINGLE", "9_10_SINGLE"]),
            (10, ["10_11_SINGLE", "9_10_SINGLE"]),
            (11, ["10_11_SINGLE", "11_12_AROMATIC", "11_16_AROMATIC"]),
            (12, ["11_12_AROMATIC", "12_13_AROMATIC"]),
            (13, ["12_13_AROMATIC", "13_14_AROMATIC"]),
            (14, ["13_14_AROMATIC", "14_15_AROMATIC"]),
            (15, ["14_15_AROMATIC", "15_16_AROMATIC"]),
            (16, ["11_16_AROMATIC", "15_16_AROMATIC"]),
            (17, ["17_18_SINGLE", "6_17_SINGLE"]),
            (18, ["17_18_SINGLE", "2_18_SINGLE"]),
        ]
    )

    example_3 = (
        "Cc1ccccc1.O=C(O)c1ccccc1.O=[CH:8][CH:9]1[CH2:10][CH2:11][CH:12]"
        "([C:13]([CH3:14])([CH3:15])[CH2:16][CH3:17])[CH2:18][CH2:19]1.c1ccc(P(c2ccccc2)"
        "(c2ccccc2)=[C:6]([C:4]([O:3][CH2:2][CH3:1])=[O:5])[CH3:7])cc1"
    )
    answer_3 = OrderedDict(
        [
            (1, ["1_2_SINGLE"]),
            (2, ["1_2_SINGLE", "2_3_SINGLE"]),
            (3, ["2_3_SINGLE", "3_4_SINGLE"]),
            (4, ["3_4_SINGLE", "4_5_DOUBLE", "4_6_SINGLE"]),
            (5, ["4_5_DOUBLE"]),
            (6, ["0_6_DOUBLE", "4_6_SINGLE", "6_7_SINGLE"]),
            (7, ["6_7_SINGLE"]),
            (8, ["0_8_DOUBLE", "8_9_SINGLE"]),
            (9, ["8_9_SINGLE", "9_10_SINGLE", "9_19_SINGLE"]),
            (10, ["10_11_SINGLE", "9_10_SINGLE"]),
            (11, ["10_11_SINGLE", "11_12_SINGLE"]),
            (12, ["11_12_SINGLE", "12_13_SINGLE", "12_18_SINGLE"]),
            (13, ["12_13_SINGLE", "13_14_SINGLE", "13_15_SINGLE", "13_16_SINGLE"]),
            (14, ["13_14_SINGLE"]),
            (15, ["13_15_SINGLE"]),
            (16, ["13_16_SINGLE", "16_17_SINGLE"]),
            (17, ["16_17_SINGLE"]),
            (18, ["12_18_SINGLE", "18_19_SINGLE"]),
            (19, ["18_19_SINGLE", "9_19_SINGLE"]),
        ]
    )

    assert get_atomic_neighbourhoods(example_1) == answer_1
    assert get_atomic_neighbourhoods(example_2) == answer_2
    assert get_atomic_neighbourhoods(example_3) == answer_3


def test_get_all_atom_indices():
    rxn_1 = (
        "CO.O=C(OCc1ccccc1)[NH:9][C@@H:7]([C:5]([NH:4][CH:2]([CH3:1])[CH3:3])="
        "[O:6])[CH3:8].[H][H].[Pd]>>[CH3:1][CH:2]([CH3:3])[NH:4][C:5](=[O:6])[C@@H:7]"
        "([CH3:8])[NH2:9]"
    )
    rxn_1_components = rxn_1.split(">>")

    rxn_2 = (
        "O=C(Cl)c1ccccc1F.O=[C:8](Cl)[c:7]1[c:2]([CH3:1])[cH:3][cH:4][cH:5][cH:6]"
        "1.[NH2:9][c:10]1[c:11]([Cl:12])[cH:13][cH:14][cH:15][c:16]1[C:17](=[O:18])[OH:19]"
        ">>[CH3:1][c:2]1[cH:3][cH:4][cH:5][cH:6][c:7]1-[c:8]1[n:9][c:10]2[c:11]([Cl:12])"
        "[cH:13][cH:14][cH:15][c:16]2[c:17](=[O:18])[o:19]1"
    )
    rxn_2_components = rxn_2.split(">>")

    rxn_3 = (
        "C1CCC2=NCCCN2CC1.C1CCOC1.CO.O.O[CH2:2][c:3]1[cH:4][c:5]([Cl:6])[cH:7]"
        "[cH:8][c:9]1[CH:10]1[CH2:11][CH2:12][N:13]1[C:14](=[O:15])[O:16][CH2:17][C:18]"
        "([Cl:19])([Cl:20])[Cl:21].[N-]=[N+]=[N:1]P(=O)(c1ccccc1)c1ccccc1.[Na+].[OH-]."
        "c1ccc(P(c2ccccc2)c2ccccc2)cc1>>[NH2:1][CH2:2][c:3]1[cH:4][c:5]([Cl:6])[cH:7]"
        "[cH:8][c:9]1[CH:10]1[CH2:11][CH2:12][N:13]1[C:14](=[O:15])[O:16][CH2:17][C:18]"
        "([Cl:19])([Cl:20])[Cl:21]"
    )
    rxn_3_components = rxn_3.split(">>")

    assert get_all_atom_indices(rxn_1_components[0], rxn_1_components[1]) == {
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
    }
    assert get_all_atom_indices(rxn_2_components[0], rxn_2_components[1]) == {
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
    }
    assert get_all_atom_indices(rxn_3_components[0], rxn_3_components[1]) == {
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
    }


def test_get_atom_list():
    rxn_1 = (
        "O=P(Cl)(Cl)[Cl:12].O[c:11]1[c:6]([C:4]([O:3][CH2:2][CH3:1])=[O:5])[cH:7]"
        "[n:8][c:9]2[c:10]1[CH2:13][CH2:14][CH2:15][CH2:16][CH2:17][CH2:18]2>>[CH3:1]"
        "[CH2:2][O:3][C:4](=[O:5])[c:6]1[cH:7][n:8][c:9]2[c:10]([c:11]1[Cl:12])[CH2:13]"
        "[CH2:14][CH2:15][CH2:16][CH2:17][CH2:18]2"
    )
    rxn_1_components = rxn_1.split(">>")

    rxn_2 = (
        "C1CCOC1.O[CH3:18].[CH3:1][CH2:2][N:3]([C:4]([CH3:5])=[O:19])[c:6]1[cH:7]"
        "[c:8]([O:9][CH2:10][c:11]2[cH:12][cH:13][cH:14][cH:15][cH:16]2)[cH:17][cH:20]"
        "[c:21]1[O:22][CH3:23]>>[CH3:1][CH2:2][N:3]([CH2:4][CH3:5])[c:6]1[cH:7][c:8]([O:9]"
        "[CH2:10][c:11]2[cH:12][cH:13][cH:14][cH:15][cH:16]2)[c:17]([CH:18]=[O:19])[cH:20]"
        "[c:21]1[O:22][CH3:23]"
    )
    rxn_2_components = rxn_2.split(">>")

    rxn_3 = (
        "ClC(Cl)Cl.O=S(Cl)[Cl:5].O[CH2:4][C:2]([CH3:1])([CH3:3])[NH:6][C:7]1=[N:8]"
        "[CH2:9][CH2:10][NH:11]1.[IH:12]>>[CH3:1][C:2]([CH3:3])([CH2:4][Cl:5])[NH:6][C:7]"
        "1=[N:8][CH2:9][CH2:10][NH:11]1"
    )
    rxn_3_components = rxn_3.split(">>")

    assert get_atom_list(
        rxn_1_components[0], rxn_1_components[1], AtomEnvironment.SAME
    ) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 15, 16, 17, 18]
    assert get_atom_list(
        rxn_1_components[0], rxn_1_components[1], AtomEnvironment.CHANGED
    ) == [11, 12]
    assert (
        len(
            set(
                get_atom_list(
                    rxn_1_components[0], rxn_1_components[1], AtomEnvironment.SAME
                )
            ).intersection(
                set(
                    get_atom_list(
                        rxn_1_components[0],
                        rxn_1_components[1],
                        AtomEnvironment.CHANGED,
                    )
                )
            )
        )
        == 0
    )

    assert get_atom_list(
        rxn_2_components[0], rxn_2_components[1], AtomEnvironment.SAME
    ) == [1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 20, 21, 22, 23]
    assert get_atom_list(
        rxn_2_components[0], rxn_2_components[1], AtomEnvironment.CHANGED
    ) == [4, 17, 18, 19]
    assert (
        len(
            set(
                get_atom_list(
                    rxn_2_components[0], rxn_2_components[1], AtomEnvironment.SAME
                )
            ).intersection(
                set(
                    get_atom_list(
                        rxn_2_components[0],
                        rxn_2_components[1],
                        AtomEnvironment.CHANGED,
                    )
                )
            )
        )
        == 0
    )

    assert get_atom_list(
        rxn_3_components[0], rxn_3_components[1], AtomEnvironment.SAME
    ) == [1, 2, 3, 6, 7, 8, 9, 10, 11]
    assert get_atom_list(
        rxn_3_components[0], rxn_3_components[1], AtomEnvironment.CHANGED
    ) == [4, 5]
    assert (
        len(
            set(
                get_atom_list(
                    rxn_3_components[0], rxn_3_components[1], AtomEnvironment.SAME
                )
            ).intersection(
                set(
                    get_atom_list(
                        rxn_3_components[0],
                        rxn_3_components[1],
                        AtomEnvironment.CHANGED,
                    )
                )
            )
        )
        == 0
    )
