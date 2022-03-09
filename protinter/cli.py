from protinter import __version__
from protinter.main import run_protinter
import click


@click.command()
@click.version_option(__version__)
@click.argument("pdb", type=click.Path(exists=True))
@click.option("-csv", "--csv", is_flag=True, help="Output CSV file")
@click.option(
    "--within_radius",
    is_flag=True,
    help="Return residues that are within distance specified distance (in Angstrom) of each other",
)
@click.option(
    "--hydrophobic", is_flag=True, help="compute hydrophobic interactions [a]"
)
@click.option("--disulphide", is_flag=True, help="compute disulphide bridges [b]")
@click.option("--ionic", is_flag=True, help="compute ionic interactions [c]")
@click.option(
    "--aroaro",
    is_flag=True,
    help="compute aromatic-aromatic interactions [d] [e]",
)
@click.option(
    "--arosul", is_flag=True, help="compute aromatic-sulphur interactions [f]"
)
@click.option("--catpi", is_flag=True, help="compute catio-pi interactions [g]")
@click.option(
    "--hb1",
    is_flag=True,
    help="compute main chain-main chain H-bonds [i] [j]",
)
@click.option(
    "--hb2",
    is_flag=True,
    help="compute main chain-side chain H-bonds [i] [j]",
)
@click.option(
    "--hb3",
    is_flag=True,
    help="compute side chain-side chain H-bonds [i] [j]",
)
@click.option(
    "--interval",
    default=0,
    type=int,
    show_default=True,
    help="minimum interval separation two AA for interaction. Default = 0",
)
@click.option(
    "--a",
    default=5.0,
    type=float,
    show_default=True,
    help="hydrophobic interactions max distance (Angstrom)",
)
@click.option(
    "--b",
    default=2.2,
    type=float,
    show_default=True,
    help="disulphide bridges max distance (Angstrom)",
)
@click.option(
    "--c",
    default=6.0,
    type=float,
    show_default=True,
    help="ionic interactions max distance (Angstrom)",
)
@click.option(
    "--d",
    default=4.5,
    type=float,
    show_default=True,
    help=" aromatic-aromatic interactions min distance (Angstrom)",
)
@click.option(
    "--e",
    default=7.0,
    type=float,
    show_default=True,
    help=" aromatic-aromatic interactions max distance (Angstrom)",
)
@click.option(
    "--f",
    default=5.3,
    type=float,
    show_default=True,
    help="aromatic-sulphur interactions max distance (Angstrom)",
)
@click.option(
    "--g",
    default=6.0,
    type=float,
    show_default=True,
    help="cation-pi interactions max distance (Angstrom)",
)
@click.option(
    "--i",
    default=3.5,
    type=float,
    show_default=True,
    help="Donor-acceptor distance cutoff (N and O) (Angstrom)",
)
@click.option(
    "--j",
    default=4,
    type=float,
    show_default=True,
    help="Donor-acceptor distance cutoff (S) (Angstrom)",
)
@click.option(
    "--w",
    default=4,
    type=float,
    show_default=True,
    help="Distance threshold for within_radius (any atom) (Angstrom)",
)
@click.option(
    "--atommindist",
    default=30,
    type=float,
    show_default=True,
    help="Two potentially interacting amino acids are disregarded if two randomly selected atoms of these two amino acids have a distance greater than this value (Angstrom)",
)
@click.option(
    "--printsequence", is_flag=True, help="Print the sequence of the protein."
)
def cli(**kwargs):
    """\b
    protinter: Protein interaction calculator
    * Homepage: https://github.com/maxibor/protinter
    * Author: Maxime Borry
    PDB: .pdb entry file
    """
    run_protinter(**kwargs)


if __name__ == "__main__":
    cli()
