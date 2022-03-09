<img src="./img/logo.png" width="300">


# ProtInter : Protein interaction calculator

**ProtInter** is a tool designed to calculate non-covalent interactions of a single chain protein in a `.pdb` file.

## How to install

```c
git clone https://github.com/maxibor/protinter.git
pip install .
```

## Example

```bash
protinter --hydrophobic -csv ./data/1bta.pdb  
```

The above example computes the hydrophobic interactions in `1BTA`, displays the results, and writes them in the file `result_hydrophobic.csv`

## List of interactions available

- hydrophobic interactions
- disulphide bridges
- ionic interactions
- aromatic-aromatic interactions
- aromatic-sulphur interactions
- cation-pi interactions
- hydrogen bonds (use results with caution)

Further informations on the interactions can be found [here](./doc/report.pdf)

## Get help

The help menu of **ProtInter** is accessible with the --help flag.

```bash
Usage: protinter [OPTIONS] PDB

  protinter: Protein interaction calculator
  * Homepage: https://github.com/maxibor/protinter
  * Author: Maxime Borry
  PDB: .pdb entry file

Options:
  --version            Show the version and exit.
  -csv, --csv          Output CSV file
  --within_radius      Return residues that are within distance specified
                       distance (in Angstrom) of each other
  --hydrophobic        compute hydrophobic interactions [a]
  --disulphide         compute disulphide bridges [b]
  --ionic              compute ionic interactions [c]
  --aroaro             compute aromatic-aromatic interactions [d] [e]
  --arosul             compute aromatic-sulphur interactions [f]
  --catpi              compute catio-pi interactions [g]
  --hb1                compute main chain-main chain H-bonds [i] [j]
  --hb2                compute main chain-side chain H-bonds [i] [j]
  --hb3                compute side chain-side chain H-bonds [i] [j]
  --interval INTEGER   minimum interval separation two AA for interaction.
                       Default = 0  [default: 0]
  --a FLOAT            hydrophobic interactions max distance (Angstrom)
                       [default: 5.0]
  --b FLOAT            disulphide bridges max distance (Angstrom)  [default:
                       2.2]
  --c FLOAT            ionic interactions max distance (Angstrom)  [default:
                       6.0]
  --d FLOAT            aromatic-aromatic interactions min distance (Angstrom)
                       [default: 4.5]
  --e FLOAT            aromatic-aromatic interactions max distance (Angstrom)
                       [default: 7.0]
  --f FLOAT            aromatic-sulphur interactions max distance (Angstrom)
                       [default: 5.3]
  --g FLOAT            cation-pi interactions max distance (Angstrom)
                       [default: 6.0]
  --i FLOAT            Donor-acceptor distance cutoff (N and O) (Angstrom)
                       [default: 3.5]
  --j FLOAT            Donor-acceptor distance cutoff (S) (Angstrom)
                       [default: 4]
  --w FLOAT            Distance threshold for within_radius (any atom)
                       (Angstrom)  [default: 4]
  --atommindist FLOAT  Two potentially interacting amino acids are disregarded
                       if two randomly selected atoms of these two amino acids
                       have a distance greater than this value (Angstrom)
                       [default: 30]
  --printsequence      Print the sequence of the protein.
  --help               Show this message and exit.
```

## Contributors

Thanks to everyone who also contributed to making ProtInter !

- [Axel Schmidt](https://github.com/Ax-Sch)