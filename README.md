# ProtInter : Protein interaction calculator

**ProtInter** is a tool designed to calculate interactions of a single chain protein in a `.pdb` file.

## Dependancies
- Python 3
- Biopython

## How to install

```c
git clone https://github.com/maxibor/proteininteraction.git
cd proteininteraction
chmod +x protinter
```

## Example

```c
./protinter --hydrophobic ./data/1bta.pdb > out.txt
```
The above example computes the hydrophobic interactions in `1BTA` and writes them in the file `out.txt`

## List of interactions available

- hydrophobic interactions
- disulphide interactions
- hydrogen bonds
- ionic interactions
- aromatic-aromatic interactions
- aromatic-sulphur interactions
- cation-pi interactions

Further informations on the interactions can be found in the `doc` directory

## Get help

The help menu of **ProtInter** is accessible with the -h or --help flag.

`./protinter -h`
