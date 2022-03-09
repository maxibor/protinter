from math import sqrt

amino = {
    "hydrophobic": ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"],
    "disulphide": ["CYS"],
    "ionic": ["ARG", "LYS", "HIS", "ASP", "GLU"],
    "aroaro": ["PHE", "TRP", "TYR"],
    "arosul": ["PHE", "TRP", "TYR", "CYS", "MET"],
    "cationpi": ["LYS", "ARG", "PHE", "TRP", "TYR"],
    "all": [
        "ALA",
        "GLY",
        "ILE",
        "LEU",
        "PRO",
        "VAL",
        "PHE",
        "TRP",
        "TYR",
        "ASP",
        "GLU",
        "ARG",
        "HIS",
        "LYS",
        "SER",
        "THR",
        "CYS",
        "MET",
        "ASN",
        "GLN",
    ],
}


def dist_cal(x1, y1, z1, x2, y2, z2):
    """
    Calculates distance between two points in 3D
    x1, y1, z1, x2, y2, z2: float
    return: distance(float)
    """
    dist = sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)
    return dist


def dist_center_mass_calc(resid1, resid2):
    """
    Tranforms residues center of masses to x,y,z coord
    resid1, resid2: biopython class residues
    return: distance(float)
    """
    dist = dist_cal(
        x1=resid1.center_mass[0],
        y1=resid1.center_mass[1],
        z1=resid1.center_mass[2],
        x2=resid2.center_mass[0],
        y2=resid2.center_mass[1],
        z2=resid2.center_mass[2],
    )
    return dist


def hydrophobicfun(atom1, atom2, dist=5):
    """
    Check for hydrophobic interactions
    atom1, atom2: biopython class atom
    dist = distance(float) cutoff in Angstrom
    return: distance(float) if < dist
    """
    no_inter_atom = ["N", "CA", "C", "O"]
    if (atom1.get_name() not in no_inter_atom and "H" not in atom1.get_name()) and (
        atom2.get_name() not in no_inter_atom and "H" not in atom2.get_name()
    ):
        d = atom1 - atom2
        if d < dist:
            return d


def disulphidefun(atom1, atom2, dist=2.2):
    """
    Check for disulphide bridges
    atom1, atom2: biopython class atom
    dist = distance(float) cutoff in Angstrom
    return: distance(float)  if < dist
    """
    if "S" in atom1.get_name() and "S" in atom2.get_name():
        d = atom1 - atom2
        if d < dist:
            return d


def ionicfun(atom1, atom2, resid1, resid2, dist=6):
    """
    Check for ionic interactions
    atom1, atom2: biopython class atom
    resid1, resid2: biopython class residues
    dist = distance(float) cutoff in Angstrom
    return: distance(float)  if < dist
    """
    if (
        (
            "N" in atom1.get_name()
            and len(atom1.get_name()) > 1
            and resid1.get_resname() in ["LYS", "ARG", "HIS"]
        )
        and (
            "O" in atom2.get_name()
            and len(atom2.get_name()) > 1
            and resid2.get_resname() in ["ASP", "GLU"]
        )
    ) or (
        (
            "O" in atom1.get_name()
            and len(atom1.get_name()) > 1
            and resid1.get_resname() in ["ASP", "GLU"]
        )
        and (
            "N" in atom2.get_name()
            and len(atom2.get_name()) > 1
            and resid2.get_resname() in ["LYS", "ARG", "HIS"]
        )
    ):
        d = atom1 - atom2
        if d < dist:
            return d


def cationpifun(atom1, atom2, resid1, resid2, dist=6):
    """
    Check for cation-pi interactions
    atom1, atom2: biopython class atom
    resid1, resid2: biopython class residues
    dist: distance(float) cutoff in Angstrom
    return: distance(float)  if < dist
    """
    interac = str(atom1.get_name()) + str(atom2.get_name())
    if (
        (
            (
                resid1.get_resname() in ["LYS", "ARG"]
                and resid2.get_resname() in ["PHE", "TYR", "TRP"]
            )
            and len(atom1.get_name()) > 1
            and atom1.get_name() != "CA"
            and len(atom2.get_name()) > 1
            and atom2.get_name() != "CA"
            and "H" not in atom1.get_name()
            and "H" not in atom2.get_name()
        )
        or (
            (
                resid2.get_resname() in ["LYS", "ARG"]
                and resid1.get_resname() in ["PHE", "TYR", "TRP"]
            )
            and len(atom1.get_name()) > 1
            and atom1.get_name() != "CA"
            and len(atom2.get_name()) > 1
            and atom2.get_name() != "CA"
            and "H" not in atom1.get_name()
            and "H" not in atom2.get_name()
        )
    ) and interac.count("C") != 2:
        d = atom1 - atom2
        if d < dist:
            return d


def aroarofun(resid1, resid2, dmin=4.5, dmax=7):
    """
    Check for aromatic-aromatic interactions
    resid1, resid2: biopython class residues
    dmin, dmax: distance(float) in Angstrom, min and max for cutoff
    return: distance(float)  if < dist
    """

    d = dist_center_mass_calc(resid1, resid2)

    if (d > dmin) and (d < dmax):
        return d


def arosulfun(resid1, resid2, dist=5.3):
    """
    Check for aromatic-sulphur interactions
    resid1, resid2: biopython class residues
    dist: distance(float) cutoff in Angstrom
    return: distance(float)  if < dist
    """
    sulres = ["CYS", "MET"]
    if resid1.get_resname() in sulres and resid2.get_resname() not in sulres:
        for atom in resid1:
            if "S" in atom.get_name():
                s_coord = atom.get_coord()
                d = dist_cal(
                    s_coord[0],
                    s_coord[1],
                    s_coord[2],
                    resid2.center_mass[0],
                    resid2.center_mass[1],
                    resid2.center_mass[2],
                )
                if d < dist:
                    return d
    elif resid2.get_resname() in sulres and resid1.get_resname() not in sulres:
        for atom in resid2:
            if "S" in atom.get_name():
                s_coord = atom.get_coord()
                d = dist_cal(
                    s_coord[0],
                    s_coord[1],
                    s_coord[2],
                    resid1.center_mass[0],
                    resid1.center_mass[1],
                    resid1.center_mass[2],
                )
                if d < dist:
                    return d


def center_mass(resid):
    """
    Calculate the center of mass of a Phenyl ring in aromatic amino-acids
    resid: biopython class residues
    return: extented biopython class residues with .center_mass(float)
    """
    xmean = 0
    ymean = 0
    zmean = 0
    if resid.get_resname() in ["PHE", "TYR"]:
        for atom in resid:
            if atom.get_name() in ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]:
                xmean = xmean + atom.get_coord()[0]
                ymean = ymean + atom.get_coord()[1]
                zmean = zmean + atom.get_coord()[2]
        xmean = xmean / 6
        ymean = ymean / 6
        zmean = zmean / 6
    elif resid.get_resname() == "TRP":
        for atom in resid:
            if atom.get_name() in ["CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3"]:
                xmean = xmean + atom.get_coord()[0]
                ymean = ymean + atom.get_coord()[1]
                zmean = zmean + atom.get_coord()[2]
        xmean = xmean / 6
        ymean = ymean / 6
        zmean = zmean / 6
    resid.center_mass = (xmean, ymean, zmean)
    return resid


def hbondfun(acceptor, donnor, distON=3.5, distS=4):
    """
    Check for hydrogen bonds
    acceptor, donnor: biopython class atoms
    distON: distance(float) cutoff in Angstrom if atom in ["O","N"]
    distS: distance(float) cutoff in Angstrom if atom in ["S"]
    return: distance(float)  if < distS or distON
    """
    d = acceptor - donnor
    if "O" in acceptor.get_name() or "N" in acceptor.get_name():
        if d < distON:
            return d
    elif "S" in acceptor.get_name():
        if d < distS:
            return d


def hbond_main_mainfun(atom1, atom2, distON=3.5, distS=4):
    """
    Check for main-chain/main-chain hydrogen bonds using hbondfun()
    atom1, atom2: biopython class atoms
    distON: distance(float) cutoff in Angstrom if atom in ["O","N"]
    distS: distance(float) cutoff in Angstrom if atom in ["S"]
    return: distance(float)  if < distS or distON
    """
    name1 = atom1.get_name()
    name2 = atom2.get_name()
    if len(name1) == 1 and len(name2) == 1:
        if name1 == "O" and name2 == "H":
            d = hbondfun(acceptor=atom1, donnor=atom2, distON=distON, distS=distS)
            if d:
                return d
        elif name2 == "O" and name1 == "H":
            d = hbondfun(acceptor=atom2, donnor=atom1, distON=distON, distS=distS)
            if d:
                return d


def hbond_main_sidefun(atom1, atom2, distON=3.5, distS=4):
    """
    Check for main-chain/side-chain hydrogen bonds using hbondfun()
    atom1, atom2: biopython class atoms
    distON: distance(float) cutoff in Angstrom if atom in ["O","N"]
    distS: distance(float) cutoff in Angstrom if atom in ["S"]
    return: distance(float)  if < distS or distON
    """
    name1 = atom1.get_name()
    name2 = atom2.get_name()
    if len(name1) == 1 and len(name2) == 2:
        if (name1 in ["O", "N", "S"]) and (
            "O" in name2 or "N" in name2 or "S" in name2
        ):
            d = hbondfun(acceptor=atom1, donnor=atom2, distON=distON, distS=distS)
            if d:
                return d
    if len(name1) == 2 and len(name2) == 1:
        if (name2 in ["O", "N", "S"]) and (
            "O" in name1 or "N" in name1 or "S" in name1
        ):
            d = hbondfun(acceptor=atom1, donnor=atom2, distON=distON, distS=distS)
            if d:
                return d


def hbond_side_sidefun(atom1, atom2, distON=3.5, distS=4):
    """
    Check for side-chain/side-chain hydrogen bonds using hbondfun()
    atom1, atom2: biopython class atoms
    distON: distance(float) cutoff in Angstrom if atom in ["O","N"]
    distS: distance(float) cutoff in Angstrom if atom in ["S"]
    return: distance(float)  if < distS or distON
    """
    name1 = atom1.get_name()
    name2 = atom2.get_name()
    if len(name1) == 2 and len(name2) == 2:
        if ("O" in name2 or "N" in name2 or "S" in name2) and (
            "O" in name1 or "N" in name1
        ):
            d = hbondfun(acceptor=atom2, donnor=atom1, distON=distON, distS=distS)
            if d:
                return d
        elif ("O" in name1 or "N" in name1 or "S" in name1) and (
            "O" in name2 or "N" in name2
        ):
            d = hbondfun(acceptor=atom1, donnor=atom2, distON=distON, distS=distS)
            if d:
                return d


def within_radiusfun(atom1, atom2, dist_min=4):
    """
    get residues that have atom pairs that are within a certain distance
    """
    d = atom1 - atom2
    d = abs(d)
    if d < dist_min:
        return d


def get_res(chain, amino_type, amino=amino):
    """
    Using the amino dictionary, sort amino acids in sequence by chemical properties
    chain: biopython class chain
    amino_type: A key(str) of the amino dictionary
    amino: dictionary
    return: A dictionary with amino acids objects binned by properties
    """
    mydic = {}
    for resid in chain:
        if resid.get_resname() in amino[amino_type]:
            resseq = str(resid).split(" ")[4].split("=")[1]
            mydic[resseq] = resid
    return mydic


def calc_inter(
    residue_dict,
    amino_type,
    filename,
    distmin=0,
    distmax=7,
    distON=3.5,
    distS=4,
    inter=0,
    csv=False,
    atommindist=20,
):
    """
    Main function
    residue_dict: output(dict) of get_res() function
    amino_type: A key(str) of the amino dictionary
    distmin: min distance(float) for interaction
    distmax: max distance(float) for interaction
    distON: distance(float) cutoff in Angstrom if atom in ["O","N"]
    distS: distance(float) cutoff in Angstrom if atom in ["S"]
    inter: minimum interval(int) between two AA for interaction
    """
    found = []
    if csv:
        to_csv = ["RES1 , idRES1 , RES2 , idRES2 , dist(Angstrom)\n"]
    keys = sorted([int(x) for x in list(residue_dict.keys())])
    for i in keys:
        for j in keys:
            if i != j and abs(i - j) > inter:
                InterOfResisFound = False
                resid1 = residue_dict[str(i)]
                resid2 = residue_dict[str(j)]
                ResiduesToFarApart = False
                WithinMinimumDist = False
                for atom1 in resid1:
                    if (
                        (amino_type == "hbond_main_side")
                        or (amino_type == "hbond_side_side")
                        or (amino_type == "hbond_main_main")
                    ):
                        if atom1.get_name() == "C":
                            continue
                    if InterOfResisFound:
                        break
                    if ResiduesToFarApart:
                        ResiduesToFarApart = False
                        res = None
                        break
                    for atom2 in resid2:
                        if WithinMinimumDist == False:
                            if abs(atom1 - atom2) > atommindist:
                                ResiduesToFarApart = True
                                res = None
                                break
                            else:
                                WithinMinimumDist = True
                        if amino_type == "hydrophobic":
                            res = hydrophobicfun(atom1, atom2, dist=distmax)
                            if res:
                                InterOfResisFound = True
                                break
                        elif amino_type == "disulphide":
                            res = disulphidefun(atom1, atom2, dist=distmax)
                            if res:
                                InterOfResisFound = True
                                break
                        elif amino_type == "ionic":
                            res = ionicfun(atom1, atom2, resid1, resid2, dist=distmax)
                            if res:
                                InterOfResisFound = True
                                break
                        elif amino_type == "cationpi":
                            res = cationpifun(
                                atom1, atom2, resid1, resid2, dist=distmax
                            )
                            if res:
                                InterOfResisFound = True
                                break
                        elif amino_type == "hbond_main_main":
                            res = hbond_main_mainfun(
                                atom1, atom2, distON=distON, distS=distS
                            )
                            if res:
                                InterOfResisFound = True
                                break
                        elif amino_type == "hbond_main_side":
                            res = hbond_main_sidefun(
                                atom1, atom2, distON=distON, distS=distS
                            )
                            if res:
                                InterOfResisFound = True
                                break
                        elif amino_type == "hbond_side_side":
                            res = hbond_side_sidefun(
                                atom1, atom2, distON=distON, distS=distS
                            )
                            if res:
                                InterOfResisFound = True
                                break
                        elif amino_type == "within_radius":
                            res = within_radiusfun(atom1, atom2, dist_min=distmax)
                            if res:
                                InterOfResisFound = True
                                break

                if amino_type == "aroaro":
                    res = aroarofun(resid1, resid2, dmin=distmin, dmax=distmax)
                elif amino_type == "arosul":
                    res = arosulfun(resid1, resid2, dist=distmax)

                if res:
                    if (resid1, resid2) not in found and (resid2, resid1) not in found:
                        found.append((resid1, resid2))
                        if csv:
                            to_print = (
                                "{:<5}".format(str(resid1).split()[1])
                                + ", "
                                + "{:<7}".format(str(resid1).split()[3].split("=")[1])
                            )
                            to_print = (
                                to_print
                                + ", "
                                + "{:<5}".format(str(resid2).split()[1])
                                + ", "
                                + "{:<7}".format(str(resid2).split()[3].split("=")[1])
                            )
                            to_print = (
                                to_print + ", " + str("{:06.2f}".format(res)) + "\n"
                            )
                            to_csv.append(to_print)

                        to_print = (
                            "| "
                            + "{:<5}".format(str(resid1).split()[1])
                            + "| "
                            + "{:<7}".format(str(resid1).split()[3].split("=")[1])
                        )
                        to_print = (
                            to_print
                            + "| "
                            + "{:<5}".format(str(resid2).split()[1])
                            + "| "
                            + "{:<7}".format(str(resid2).split()[3].split("=")[1])
                        )
                        to_print = (
                            to_print + "| " + str("{:06.2f}".format(res)) + "        |"
                        )
                        print(to_print)

    print("|-----------------------------------------------|")
    print(
        "| Total: "
        + "{:<5}".format(str(len(found)))
        + "interactions                     |"
    )
    print(" -----------------------------------------------\n ")
    if csv:
        with open("result" + "_" + amino_type + ".csv", "w") as fw:
            for line in to_csv:
                fw.write(line)
