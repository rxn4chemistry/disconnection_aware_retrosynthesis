from disconnection_aware_retrosynthesis.tagging import (
    find_number_tags,
    get_tagged_products,
    return_tag_combinations,
)


def test_get_tagged_products():
    rxn_1 = (
        "CC[O:23][C:21]([CH2:20][CH:19]1[CH2:18][CH2:17][N:16]([C:14]([CH2:13]"
        "[CH2:12][CH2:11][C@H:10]2[O:9][C@H:8]([c:7]3[cH:6][cH:5][cH:4][c:3]([O:2]"
        "[CH3:1])[c:38]3[O:39][CH3:40])[c:37]3[c:31]([cH:32][cH:33][c:34]([Cl:35])"
        "[cH:36]3)-[n:30]3[c:26]2[cH:27][cH:28][cH:29]3)=[O:15])[CH2:25][CH2:24]1)="
        "[O:22].Cl.O.O=C([O-])[O-].OCC1CCCO1.[K+]>>[CH3:1][O:2][c:3]1[cH:4][cH:5][cH:6]"
        "[c:7]([C@H:8]2[O:9][C@H:10]([CH2:11][CH2:12][CH2:13][C:14](=[O:15])[N:16]3"
        "[CH2:17][CH2:18][CH:19]([CH2:20][C:21](=[O:22])[OH:23])[CH2:24][CH2:25]3)[c:26]"
        "3[cH:27][cH:28][cH:29][n:30]3-[c:31]3[cH:32][cH:33][c:34]([Cl:35])[cH:36][c:37]"
        "32)[c:38]1[O:39][CH3:40]"
    )
    rxn_1_components = rxn_1.split(">>")

    rxn_2 = (
        "CCN(C(C)C)C(C)C.CN(C)C=O.Cc1ccc(S(=O)(=O)O[CH2:10][CH2:11][CH:12]2"
        "[CH2:13][O:14][c:15]3[cH:16][cH:17][cH:18][cH:19][c:20]3[O:21]2)cc1.[CH3:1]"
        "[CH2:2][O:3][C:4](=[O:5])[CH:6]1[CH2:7][CH2:8][NH:9][CH2:22][CH2:23]1>>[CH3:1]"
        "[CH2:2][O:3][C:4](=[O:5])[CH:6]1[CH2:7][CH2:8][N:9]([CH2:10][CH2:11][CH:12]2"
        "[CH2:13][O:14][c:15]3[cH:16][cH:17][cH:18][cH:19][c:20]3[O:21]2)[CH2:22]"
        "[CH2:23]1"
    )
    rxn_2_components = rxn_2.split(">>")

    rxn_3 = (
        "Br[CH2:18][c:19]1[cH:20][cH:21][cH:22][cH:23][c:24]1[C:25]([F:26])"
        "([F:27])[F:28].C1COCCO1.CN(C)C=O.CO.Cl[C:2](=[O:1])[c:29]1[cH:30][cH:31][cH:32]"
        "[cH:33][cH:34]1.O.O=C([O-])[O-].[Cs+].[NH2:3][C@@H:4]([CH2:5][c:6]1[cH:7][nH:8]"
        "[c:9]2[cH:10][cH:11][cH:12][cH:13][c:14]12)[C:15](=[O:16])[OH:17].[Na+]>>[O:1]="
        "[C:2]([NH:3][CH:4]([CH2:5][c:6]1[cH:7][nH:8][c:9]2[cH:10][cH:11][cH:12][cH:13]"
        "[c:14]12)[C:15](=[O:16])[O:17][CH2:18][c:19]1[cH:20][cH:21][cH:22][cH:23][c:24]"
        "1[C:25]([F:26])([F:27])[F:28])[c:29]1[cH:30][cH:31][cH:32][cH:33][cH:34]1"
    )
    rxn_3_components = rxn_3.split(">>")

    assert (
        get_tagged_products(rxn_1_components[0], rxn_1_components[1])
        == "COc1cccc([C@H]2O[C@H](CCCC(=O)N3CCC(CC(=O)[OH:1])CC3)c3cccn3-c3ccc(Cl)"
        "cc32)c1OC"
    )
    assert (
        get_tagged_products(rxn_2_components[0], rxn_2_components[1])
        == "CCOC(=O)C1CC[N:1]([CH2:1]CC2COc3ccccc3O2)CC1"
    )
    assert (
        get_tagged_products(rxn_3_components[0], rxn_3_components[1])
        == "O=C(C(Cc1c[nH]c2ccccc12)[NH:1][C:1](=O)c1ccccc1)[O:1][CH2:1]c1ccccc1C(F)"
        "(F)F"
    )


def test_find_number_tags():
    assert (
        find_number_tags(
            "Cc1nc(C)c(-c2cn(CCC[CH2:1][N:1]3C[C@@H]4C[C@]4(c4ccc(C(F)(F)F)cc4)C3)c"
            "(=O)[nH]c2=O)o1"
        )
        == 2
    )
    assert (
        find_number_tags("Fc1ccc(C(=C(C2[CH2:1][CH2:1]2)[CH2:1][OH:1])c2ccc(F)cc2)cc1")
        == 4
    )
    assert find_number_tags("CC[c:1]1[cH:1][c:1]([NH2:1])[cH:1][cH:1][c:1]1[OH:1]") == 8


def test_return_tag_combinations():
    assert (
        return_tag_combinations(
            "Cc1nc(C)c(-c2cn(CCC[CH2:1][N:1]3C[C@@H]4C[C@]4(c4ccc(C(F)(F)F)cc4)C3)c"
            "(=O)[nH]c2=O)o1"
        )
        == 2
    )
    assert (
        return_tag_combinations(
            "Fc1ccc(C(=C(C2[CH2:1][CH2:1]2)[CH2:1][OH:1])c2ccc(F)cc2)cc1"
        )
        == 14
    )
    assert (
        return_tag_combinations("CC[c:1]1[cH:1][c:1]([NH2:1])[cH:1][cH:1][c:1]1[OH:1]")
        == 92
    )
