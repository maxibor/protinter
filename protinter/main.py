#!/usr/bin/python3

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import protinter.interlib as pi
import textwrap


def run_protinter(
    pdb,
    csv=False,
    hydrophobic=False,
    disulphide=False,
    ionic=False,
    aroaro=False,
    arosul=False,
    catpi=False,
    hb1=False,
    hb2=False,
    hb3=False,
    within_radius=False,
    a=5.0,
    b=2.2,
    c=6.0,
    d=4.5,
    e=7.0,
    f=5.3,
    g=6.0,
    i=3.5,
    j=4.0,
    w=4.0,
    interval=0,
    atommindist=30,
    printsequence=False,
):
    """
    Main function to run protinter
    Args:
        pdb (str): path to the pdb file
        csv (bool): if True, output a csv file
        hydrophobic (bool): if True, compute hydrophobic interactions [a]
        disulphide (bool): if True, compute disulphide bridges [b]
        ionic (bool): if True, compute ionic interactions [c]
        aroaro (bool): if True, compute aromatic-aromatic interactions [d] [e]
        arosul (bool): if True, compute cation-pi interactions
        catpi (bool): if True, compute catio-pi interactions [g]
        hb1 (bool): if True, compute main chain-main chain H-bonds [i] [j]
        hb2 (bool): if True, compute main chain-side chain H-bonds [i] [j]
        hb3 (bool): if True, compute side chain-side chain H-bonds [i] [j]
        within_radius (float): return residues that are within distance in Angstrom of each other
        a (float): hydrophobic interactions max distance (Angstrom)
        b (float): disulphide bridges max distance (Angstrom)
        c (float): ionic interactions max distance (Angstrom)
        d (float): aromatic-aromatic interactions max distance (Angstrom)
        e (float): aromatic-aromatic interactions max distance (Angstrom)
        f (float): aromatic-sulphur interactions max distance (Angstrom)
        g (float): catio-pi interactions max distance (Angstrom)
        i (float): Donor-acceptor distance cutoff (N and O) (Angstrom)
        j (float): Donor-acceptor distance cutoff (S) (Angstrom)
        w (float): Distance threshold for within_radius (any atom) (Angstrom)
        interval (int): minimum interval separation two AA for interaction. Default = 0
        atommindist (float): minimum distance between two atoms (Angstrom)
        printsequence (bool): if True, print the sequence of the protein
    """
    p = PDBParser()
    structure = p.get_structure("X", pdb)

    def get_residues(type_get):
        for model in structure:
            for chain in model:
                for resid in chain:
                    if (type_get == "aroaro") | (type_get == "arosul"):
                        if resid.get_resname() in pi.amino["aroaro"]:
                            pi.center_mass(resid)
                    return_resid = pi.get_res(chain, amino_type=type_get)
                    # disulphide = pi.get_res(chain, amino_type="disulphide")
                    # ionic = pi.get_res(chain, amino_type="ionic")
                    # catpi = pi.get_res(chain, amino_type="cationpi")
                    # aroaro = pi.get_res(chain, amino_type="aroaro")
                    # arosul = pi.get_res(chain, amino_type="arosul")
                    # hbond = pi.get_res(chain, amino_type="all")
        return return_resid

    if printsequence:
        ppb = PPBuilder()

        print(" ----------------------------------------------- ")
        print("|Sequence of the protein in " + "{:20.20}".format(str(pdb)) + "|")
        print(" ----------------------------------------------- ")
        for pp in ppb.build_peptides(structure):
            print(textwrap.fill(str(pp.get_sequence()), 49))

        print("\n")

    if hydrophobic:
        print(" ----------------------------------------------- ")
        print("| Hydrophobic Interactions                      |")
        print("|-----------------------------------------------|")
        print("| RES1 | idRES1 | RES2 | idRES2 | dist(Angstrom)|")
        pi.calc_inter(
            residue_dict=get_residues("hydrophobic"),
            csv=csv,
            filename=pdb,
            distmax=a,
            amino_type="hydrophobic",
            inter=interval,
            atommindist=atommindist,
        )

    if disulphide:
        print(" ----------------------------------------------- ")
        print("| Disulphide Interactions                       |")
        print("|-----------------------------------------------|")
        print("| RES1 | idRES1 | RES2 | idRES2 | dist(Angstrom)|")
        pi.calc_inter(
            residue_dict=get_residues("disulphide"),
            csv=csv,
            filename=pdb,
            distmax=b,
            amino_type="disulphide",
            inter=interval,
            atommindist=atommindist,
        )

    if ionic:
        print(" ----------------------------------------------- ")
        print("| Ionic Interactions                            |")
        print("|-----------------------------------------------|")
        print("| RES1 | idRES1 | RES2 | idRES2 | dist(Angstrom)|")
        pi.calc_inter(
            residue_dict=get_residues("ionic"),
            csv=csv,
            filename=pdb,
            distmax=c,
            amino_type="ionic",
            inter=interval,
            atommindist=atommindist,
        )

    if catpi:
        print(" ----------------------------------------------- ")
        print("| Cation Pi Interactions                        |")
        print("|-----------------------------------------------|")
        print("| RES1 | idRES1 | RES2 | idRES2 | dist(Angstrom)|")
        pi.calc_inter(
            residue_dict=get_residues("cationpi"),
            csv=csv,
            filename=pdb,
            distmax=g,
            amino_type="cationpi",
            inter=interval,
            atommindist=atommindist,
        )

    if aroaro:
        print(" ----------------------------------------------- ")
        print("| Aromatic-aromatic Interactions                |")
        print("|-----------------------------------------------|")
        print("| RES1 | idRES1 | RES2 | idRES2 | dist(Angstrom)|")
        pi.calc_inter(
            residue_dict=get_residues("aroaro"),
            csv=csv,
            filename=pdb,
            distmin=d,
            distmax=e,
            amino_type="aroaro",
            inter=interval,
            atommindist=atommindist,
        )

    if arosul:
        print(" ----------------------------------------------- ")
        print("| Aromatic-sulphur Interactions                 |")
        print("|-----------------------------------------------|")
        print("| RES1 | idRES1 | RES2 | idRES2 | dist(Angstrom)|")
        pi.calc_inter(
            residue_dict=get_residues("arosul"),
            csv=csv,
            filename=pdb,
            distmax=f,
            amino_type="arosul",
            inter=interval,
            atommindist=atommindist,
        )

    if hb1:
        print(" ----------------------------------------------- ")
        print("| Hydrogen bonds MainChain - MainChain          |")
        print("|-----------------------------------------------|")
        print("| RES1 | idRES1 | RES2 | idRES2 | dist(Angstrom)|")
        pi.calc_inter(
            residue_dict=get_residues("all"),
            csv=csv,
            filename=pdb,
            distON=i,
            distS=j,
            amino_type="hbond_main_main",
            inter=interval,
            atommindist=atommindist,
        )

    if hb2:
        print(" ----------------------------------------------- ")
        print("| Hydrogen bonds MainChain - SideChain          |")
        print("|-----------------------------------------------|")
        print("| RES1 | idRES1 | RES2 | idRES2 | dist(Angstrom)|")
        pi.calc_inter(
            residue_dict=get_residues("all"),
            csv=csv,
            filename=pdb,
            distON=i,
            distS=j,
            amino_type="hbond_main_side",
            inter=interval,
            atommindist=atommindist,
        )

    if hb3:
        print(" ----------------------------------------------- ")
        print("| Hydrogen bonds SideChain - SideChain          |")
        print("|-----------------------------------------------|")
        print("| RES1 | idRES1 | RES2 | idRES2 | dist(Angstrom)|")
        pi.calc_inter(
            residue_dict=get_residues("all"),
            csv=csv,
            filename=pdb,
            distON=i,
            distS=j,
            amino_type="hbond_side_side",
            inter=interval,
            atommindist=atommindist,
        )

    if within_radius:
        print(" ----------------------------------------------- ")
        print("| Residues close to each other                  |")
        print("|-----------------------------------------------|")
        print("| RES1 | idRES1 | RES2 | idRES2 | dist(Angstrom)|")
        pi.calc_inter(
            residue_dict=get_residues("all"),
            csv=csv,
            filename=pdb,
            distmax=w,
            amino_type="within_radius",
            inter=interval,
            atommindist=atommindist,
        )
