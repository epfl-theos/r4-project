import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from ase import Atoms
from mendeleev import element
from pymatgen.core.structure import Structure
from skcosmo.preprocessing import StandardFlexibleScaler
from sklearn.ensemble import (
    ExtraTreesClassifier,
    GradientBoostingClassifier,
    RandomForestClassifier,
)
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier

from constants import PG_LIST_ORDER

sns.set(style="white", palette="bright", color_codes=True)


def _point_group_and_num_operations(sg):
    if sg == 1:  # spacegroup
        pg = "1"  # pointgroup
        o = 1  # number of operations
    elif sg == 2:
        pg = "1_"
        o = 2
    elif 3 <= sg <= 5:
        pg = "2"
        o = 2
    elif 6 <= sg <= 9:
        pg = "m"
        o = 2
    elif 10 <= sg <= 15:
        pg = "2/m"
        o = 4
    elif 16 <= sg <= 24:
        pg = "222"
        o = 4
    elif 25 <= sg <= 46:
        pg = "mm2"
        o = 4
    elif 47 <= sg <= 74:
        pg = "mmm"
        o = 8
    elif 75 <= sg <= 80:
        pg = "4"
        o = 4
    elif 81 <= sg <= 82:
        pg = "4_"
        o = 4
    elif 83 <= sg <= 88:
        pg = "4/m"
        o = 8
    elif 89 <= sg <= 98:
        pg = "422"
        o = 8
    elif 99 <= sg <= 110:
        pg = "4mm"
        o = 8
    elif 111 <= sg <= 122:
        pg = "4_2m"
        o = 8
    elif 123 <= sg <= 142:
        pg = "4/mmm"
        o = 16
    elif 143 <= sg <= 146:
        pg = "3"
        o = 3
    elif 147 <= sg <= 148:
        pg = "3_"
        o = 6
    elif 149 <= sg <= 155:
        pg = "32"
        o = 6
    elif 156 <= sg <= 161:
        pg = "3m"
        o = 6
    elif 162 <= sg <= 167:
        pg = "3_m"
        o = 12
    elif 168 <= sg <= 173:
        pg = "6"
        o = 6
    elif sg == 174:
        pg = "6_"
        o = 6
    elif 175 <= sg <= 176:
        pg = "6/m"
        o = 12
    elif 177 <= sg <= 182:
        pg = "622"
        o = 12
    elif 183 <= sg <= 186:
        pg = "6mm"
        o = 12
    elif 187 <= sg <= 190:
        pg = "6_m2"
        o = 12
    elif 191 <= sg <= 194:
        pg = "6/mmm"
        o = 24
    elif 195 <= sg <= 199:
        pg = "23"
        o = 12
    elif 200 <= sg <= 206:
        pg = "m_3"
        o = 24
    elif 207 <= sg <= 214:
        pg = "432"
        o = 24
    elif 215 <= sg <= 220:
        pg = "4_3m"
        o = 24
    elif 221 <= sg <= 230:
        pg = "m3_m"
        o = 48
    else:
        raise ValueError(f"Unknown space group: {sg}")

    return pg, o


def point_group(space_group):
    return _point_group_and_num_operations(space_group)[0]


def num_operations(space_group):
    return _point_group_and_num_operations(space_group)[1]


def magic_four(num_atoms):
    idx = [i for i in range(len(num_atoms)) if num_atoms[i] % 4 == 0]
    return idx


def non_magic_four(num_atoms):
    idx = [i for i in range(len(num_atoms)) if num_atoms[i] % 4 != 0]
    return idx


def inh_symm(data_frame):
    counts = data_frame.value_counts().to_frame().reindex(PG_LIST_ORDER).fillna(0)
    for i in range(len(data_frame)):
        elem = data_frame.iloc[i]
        if elem == "m" or elem == "2" or elem == "3" or elem == "1_":
            counts.loc["1"] += 1
        elif elem == "222" or elem == "4" or elem == "4_":
            counts.loc["1"] += 1
            counts.loc["2"] += 1
        elif elem == "mm2":
            counts.loc["1"] += 1
            counts.loc["2"] += 1
            counts.loc["m"] += 1
        elif elem == "2/m":
            counts.loc["1"] += 1
            counts.loc["2"] += 1
            counts.loc["m"] += 1
            counts.loc["1_"] += 1
        elif elem == "4_2m":
            counts.loc["4"] += 1
            counts.loc["mm2"] += 1
            counts.loc["222"] += 1
            counts.loc["1"] += 1
            counts.loc["2"] += 1
            counts.loc["m"] += 1
        elif elem == "4mm":
            counts.loc["4"] += 1
            counts.loc["mm2"] += 1
            counts.loc["1"] += 1
            counts.loc["2"] += 1
            counts.loc["m"] += 1
        elif elem == "422":
            counts.loc["4"] += 1
            counts.loc["222"] += 1
            counts.loc["1"] += 1
            counts.loc["2"] += 1
        elif elem == "4/m":
            counts.loc["4"] += 1
            counts.loc["2/m"] += 1
            counts.loc["4_"] += 1
            counts.loc["1"] += 1
            counts.loc["2"] += 1
            counts.loc["m"] += 1
            counts.loc["1_"] += 1
        elif elem == "mmm":
            counts.loc["222"] += 1
            counts.loc["2/m"] += 1
            counts.loc["mm2"] += 1
            counts.loc["1"] += 1
            counts.loc["2"] += 1
            counts.loc["m"] += 1
            counts.loc["1_"] += 1
        elif elem == "4/mmm":
            counts.loc["4_2m"] += 1
            counts.loc["4mm"] += 1
            counts.loc["422"] += 1
            counts.loc["4/m"] += 1
            counts.loc["mmm"] += 1
            counts.loc["222"] += 1
            counts.loc["4"] += 1
            counts.loc["4_"] += 1
            counts.loc["2/m"] += 1
            counts.loc["mm2"] += 1
            counts.loc["1"] += 1
            counts.loc["2"] += 1
            counts.loc["m"] += 1
            counts.loc["1_"] += 1
        elif elem == "3m" or elem == "6_":
            counts.loc["1"] += 1
            counts.loc["3"] += 1
            counts.loc["m"] += 1
        elif elem == "32" or elem == "6":
            counts.loc["1"] += 1
            counts.loc["3"] += 1
            counts.loc["2"] += 1
        elif elem == "3_":
            counts.loc["1"] += 1
            counts.loc["3"] += 1
            counts.loc["1_"] += 1
        elif elem == "3_m":
            counts.loc["3_"] += 1
            counts.loc["32"] += 1
            counts.loc["3m"] += 1
            counts.loc["2/m"] += 1
            counts.loc["1_"] += 1
            counts.loc["1"] += 1
            counts.loc["3"] += 1
            counts.loc["m"] += 1
            counts.loc["2"] += 1
        elif elem == "6/m":
            counts.loc["3_"] += 1
            counts.loc["6_"] += 1
            counts.loc["6"] += 1
            counts.loc["2/m"] += 1
            counts.loc["1_"] += 1
            counts.loc["1"] += 1
            counts.loc["3"] += 1
            counts.loc["m"] += 1
            counts.loc["2"] += 1
        elif elem == "622":
            counts.loc["222"] += 1
            counts.loc["6"] += 1
            counts.loc["32"] += 1
            counts.loc["1"] += 1
            counts.loc["3"] += 1
            counts.loc["2"] += 1
        elif elem == "6_m2":
            counts.loc["32"] += 1
            counts.loc["6_"] += 1
            counts.loc["3m"] += 1
            counts.loc["mm2"] += 1
            counts.loc["1"] += 1
            counts.loc["3"] += 1
            counts.loc["m"] += 1
            counts.loc["2"] += 1
        elif elem == "6mm":
            counts.loc["3_"] += 1
            counts.loc["6"] += 1
            counts.loc["2/m"] += 1
            counts.loc["1_"] += 1
            counts.loc["1"] += 1
            counts.loc["3"] += 1
            counts.loc["m"] += 1
            counts.loc["2"] += 1
        elif elem == "6/mmm":
            counts.loc["6mm"] += 1
            counts.loc["6_m2"] += 1
            counts.loc["3_m"] += 1
            counts.loc["622"] += 1
            counts.loc["6/m"] += 1
            counts.loc["mmm"] += 1
            counts.loc["222"] += 1
            counts.loc["2/m"] += 1
            counts.loc["mm2"] += 1
            counts.loc["1"] += 1
            counts.loc["2"] += 1
            counts.loc["m"] += 1
            counts.loc["1_"] += 1
            counts.loc["3"] += 1
            counts.loc["3m"] += 1
            counts.loc["32"] += 1
            counts.loc["3_"] += 1
            counts.loc["6_"] += 1
            counts.loc["6"] += 1

        elif elem == "23":
            counts.loc["3"] += 1
            counts.loc["222"] += 1
            counts.loc["2"] += 1
            counts.loc["1"] += 1
        elif elem == "4_3m":
            counts.loc["4_2m"] += 1
            counts.loc["23"] += 1
            counts.loc["3m"] += 1
            counts.loc["4"] += 1
            counts.loc["mm2"] += 1
            counts.loc["222"] += 1
            counts.loc["1"] += 1
            counts.loc["2"] += 1
            counts.loc["m"] += 1
            counts.loc["3"] += 1
        elif elem == "432":
            counts.loc["422"] += 1
            counts.loc["23"] += 1
            counts.loc["32"] += 1
            counts.loc["4"] += 1
            counts.loc["222"] += 1
            counts.loc["1"] += 1
            counts.loc["2"] += 1
            counts.loc["3"] += 1
        elif elem == "m_3":
            counts.loc["mmm"] += 1
            counts.loc["23"] += 1
            counts.loc["3_"] += 1
            counts.loc["2/m"] += 1
            counts.loc["mm2"] += 1
            counts.loc["m"] += 1
            counts.loc["222"] += 1
            counts.loc["1"] += 1
            counts.loc["2"] += 1
            counts.loc["3"] += 1
            counts.loc["1_"] += 1

        elif elem == "m3_m":
            counts.loc["4_3m"] += 1
            counts.loc["432"] += 1
            counts.loc["m_3"] += 1
            counts.loc["3_m"] += 1
            counts.loc["4/mmm"] += 1
            counts.loc["23"] += 1
            counts.loc["4_2m"] += 1
            counts.loc["4mm"] += 1
            counts.loc["422"] += 1
            counts.loc["4/m"] += 1
            counts.loc["mmm"] += 1
            counts.loc["222"] += 1
            counts.loc["4"] += 1
            counts.loc["4_"] += 1
            counts.loc["2/m"] += 1
            counts.loc["mm2"] += 1
            counts.loc["1"] += 1
            counts.loc["2"] += 1
            counts.loc["m"] += 1
            counts.loc["1_"] += 1
            counts.loc["3_"] += 1
            counts.loc["32"] += 1
            counts.loc["3m"] += 1
            counts.loc["3"] += 1
    return counts


def get_pymatgen(atoms, cls=None):
    """
    Returns pymatgen structure from ASE Atoms.

    Args:
        atoms: ASE Atoms object
        cls: The Structure class to instantiate (defaults to pymatgen structure)

    Returns:
        Equivalent pymatgen.core.structure.Structure
    """
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()
    lattice = atoms.get_cell()

    cls = Structure if cls is None else cls
    return cls(lattice, symbols, positions, coords_are_cartesian=True)


def group_conv(a):
    if a == 1:
        group = "Alkali metals"
    elif a == 2:
        group = "Alkaline earth metals"
    elif a is None or a >= 3 or a <= 12:
        group = "Transition metals"
    elif a == 13:
        group = "Icosagens"
    elif a == 14:
        group = "Crystallogens"
    elif a == 15:
        group = "Pnictogens"
    elif a == 16:
        group = "Chalcogens"
    elif a == 17:
        group = "Halogens"
    elif a == 18:
        group = "Nobles gases"
    return group


def get_atoms(structure, **kwargs):
    """
    Returns ASE Atoms object from pymatgen structure or molecule.
    Args:
        structure: pymatgen.core.structure.Structure or pymatgen.core.structure.Molecule
        **kwargs: other keyword args to pass into the ASE Atoms constructor
    Returns:
        ASE Atoms object
    """
    symbols = [str(site.specie.symbol) for site in structure]
    positions = [site.coords for site in structure]
    if hasattr(structure, "lattice"):
        cell = structure.lattice.matrix
        pbc = True
    else:
        cell = None
        pbc = None
    return Atoms(symbols=symbols, positions=positions, pbc=pbc, cell=cell, **kwargs)


def unit_normal(a, b, c):
    x = np.linalg.det([[1, a[1], a[2]], [1, b[1], b[2]], [1, c[1], c[2]]])
    y = np.linalg.det([[a[0], 1, a[2]], [b[0], 1, b[2]], [c[0], 1, c[2]]])
    z = np.linalg.det([[a[0], a[1], 1], [b[0], b[1], 1], [c[0], c[1], 1]])
    magnitude = (x**2 + y**2 + z**2) ** 0.5
    return (x / magnitude, y / magnitude, z / magnitude)


# area of polygon poly
def poly_area(poly):
    if len(poly) < 3:  # not a plane - no area
        return 0
    total = [0, 0, 0]
    N = len(poly)
    for i in range(N):
        vi1 = poly[i]
        vi2 = poly[(i + 1) % N]
        prod = np.cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = np.dot(total, unit_normal(poly[0], poly[1], poly[2]))
    return abs(result / 2)


"""Counting number of atoms from string, for the pandas Materials Project
dataset specifically"""


def natoms(string):
    if string[(string.find("Sites")) + 8] == ")":
        n = (string.find("Sites")) + 7
        return int(string[n])
    else:
        n1 = (string.find("Sites")) + 7
        n2 = (string.find("Sites")) + 8
        return int(string[n1] + string[n2])


def system(bl):
    if bl == "cP" or bl == "cI" or bl == "cF":
        return "cubic"
    elif bl == "tP" or bl == "tI":
        return "tetragonal"
    elif bl == "oP" or bl == "oI" or bl == "oF" or bl == "oC" or bl == "oA":
        return "orthorhombic"
    elif bl == "mP" or bl == "mC":
        return "monoclinic"
    elif bl == "hP":
        return "hexagonal"
    elif bl == "hR":
        return "trigonal"
    elif bl == "aP":
        return "triclinic"


# def primitive_structure_from_cif(cif, parse_engine, symprec, site_tolerance):
#     """Attempt to parse the given `CifData` and create a `StructureData` from it.
#     First the raw CIF file is parsed with the given `parse_engine`. The resulting `StructureData` is then passed through
#     SeeKpath to try and get the primitive cell. If that is successful, important structural parameters as determined by
#     SeeKpath will be set as extras on the structure node which is then returned as output.
#     :param cif: the `CifData` node
#     :param parse_engine: the parsing engine, supported libraries 'ase' and 'pymatgen'
#     :param symprec: a `Float` node with symmetry precision for determining primitive cell in SeeKpath
#     :param site_tolerance: a `Float` node with the fractional coordinate distance tolerance for finding overlapping
#         sites. This will only be used if the parse_engine is pymatgen
#     :return: the primitive `StructureData` as determined by SeeKpath
#     """
#     CifCleanWorkChain = WorkflowFactory(
#         "codtools.cif_clean"
#     )  # pylint: disable=invalid-name

#     try:
#         structure = cif.get_structure(
#             converter=parse_engine.value,
#             site_tolerance=site_tolerance.value,
#             store=False,
#         )
#     except exceptions.UnsupportedSpeciesError:
#         return CifCleanWorkChain.exit_codes.ERROR_CIF_HAS_UNKNOWN_SPECIES
#     except InvalidOccupationsError:
#         return CifCleanWorkChain.exit_codes.ERROR_CIF_HAS_INVALID_OCCUPANCIES
#     except Exception:  # pylint: disable=broad-except
#         return CifCleanWorkChain.exit_codes.ERROR_CIF_STRUCTURE_PARSING_FAILED

#     try:
#         seekpath_results = get_kpoints_path(structure, symprec=symprec)
#     except ValueError:
#         return CifCleanWorkChain.exit_codes.ERROR_SEEKPATH_INCONSISTENT_SYMMETRY
#     except SymmetryDetectionError:
#         return CifCleanWorkChain.exit_codes.ERROR_SEEKPATH_SYMMETRY_DETECTION_FAILED

#     # Store important information that should be easily queryable as attributes in the StructureData
#     parameters = seekpath_results["parameters"].get_dict()
#     structure = seekpath_results["primitive_structure"]

#     # Store the formula as a string, in both hill as well as hill-compact notation, so it can be easily queried for
#     extras = {
#         "formula_hill": structure.get_formula(mode="hill"),
#         "formula_hill_compact": structure.get_formula(mode="hill_compact"),
#         "chemical_system": "-{}-".format("-".join(sorted(structure.get_symbols_set()))),
#     }

#     for key in [
#         "spacegroup_international",
#         "spacegroup_number",
#         "bravais_lattice",
#         "bravais_lattice_extended",
#     ]:
#         try:
#             extras[key] = parameters[key]
#         except KeyError:
#             pass

#     structure.set_extra_many(extras)

#     return structure


def get_r(list):
    r = []
    for i in list:
        el_str = element(str(i))
        rad = el_str.atomic_radius
        if i == "Xe":
            rad = float(108)
        if i == "Kr":
            rad = float(88)
        if i == "Rn":
            rad = float(120)
        r.append(rad)
    return r


def comp_classif(X, Y):

    models = [
        "Logistic Regression",
        "Support Vector Machine",
        "Decision Tree",
        "Random Forest",
        "Extra Tree",
        "Gradient Boost",
    ]
    acc_train = []
    acc_test = []

    i_train, i_test, X_train, X_test, y_train, y_test = train_test_split(
        np.arange(X.shape[0]), X, Y, train_size=0.8
    )
    x_scaler = StandardFlexibleScaler(column_wise=False).fit(X)
    X = x_scaler.transform(X)
    X_train = x_scaler.transform(X_train)
    X_test = x_scaler.transform(X_test)

    logr = LogisticRegression(solver="liblinear", C=0.5, max_iter=1000)
    logr.fit(X_train, y_train)
    logr.predict(X_test)
    acc_train.append(logr.score(X_train, y_train))
    acc_test.append(logr.score(X_test, y_test))

    svm = SVC(kernel="rbf")
    svm.fit(X_train, y_train)
    svm.predict(X_test)
    acc_train.append(svm.score(X_train, y_train))
    acc_test.append(svm.score(X_test, y_test))

    tree = DecisionTreeClassifier(max_depth=None, random_state=2)
    tree.fit(X_train, y_train)
    tree.predict(X_test)
    acc_train.append(tree.score(X_train, y_train))
    acc_test.append(tree.score(X_test, y_test))

    RF = RandomForestClassifier(n_estimators=100, max_depth=None, random_state=2)
    RF.fit(X_train, y_train)
    RF.predict(X_test)
    acc_train.append(RF.score(X_train, y_train))
    acc_test.append(RF.score(X_test, y_test))

    extra_tree = ExtraTreesClassifier(n_estimators=100, max_depth=None, random_state=2)
    extra_tree.fit(X_train, y_train)
    extra_tree.predict(X_test)
    acc_train.append(extra_tree.score(X_train, y_train))
    acc_test.append(extra_tree.score(X_test, y_test))

    gr_boost = GradientBoostingClassifier(loss="deviance", validation_fraction=0.2)
    gr_boost.fit(X_train, y_train)
    gr_boost.predict(X_test)
    acc_train.append(gr_boost.score(X_train, y_train))
    acc_test.append(gr_boost.score(X_test, y_test))

    sns.scatterplot(x=acc_train, y=acc_test, hue=models, palette="Dark2", s=100)
    plt.xlabel("Train set accuracy")
    plt.ylabel("Test set accuracy")
    plt.grid(True)
    sns.despine(left=True, bottom=True)
    #     plt.legend(frameon=False)
    plt.show()
    plt.close()

    sns.barplot(x=acc_test, y=models, palette="Dark2")
    plt.xlabel("Test set accuracy")
    plt.ylabel("Model")
    plt.grid(True)
    sns.despine(left=True, bottom=True)
    #     plt.legend(frameon=False)

    f = print(
        "Classification models: "
        + str(models)
        + "\n Train accuracies: "
        + str(acc_train)
        + "\n Test accuracies: "
        + str(acc_test)
    )

    return (sns, f)
