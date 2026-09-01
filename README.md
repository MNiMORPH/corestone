# corestone

**ONE LINE: what this model evolves, and of what.**

A short paragraph: the physics in words. What is conserved, what drives the
change, what sets the rate, and what question the model was built to answer.

## Sources to cite

### Base / always

**Author, A. (YEAR), [Title](https://doi.org/...), *Journal*, *vol*, pages,
doi:...**

### Software

**Wickert, A. D. (YEAR). [corestone](https://doi.org/...) (Version X.Y.Z).**

## Installation

### Via pip and PyPI

```sh
pip install corestone
```

If your computer shields the core Python install from external packages, either
pass `--break-system-packages` (fine in my experience, but packages can clash)
or build a separate environment.

### Locally with pip, incorporating ongoing code modifications

```sh
git clone https://github.com/MNiMORPH/corestone.git
cd corestone
pip install -e ".[test]"
```

Run the tests with:

```sh
pytest
```

## Learning how to use corestone

See `examples/`. Each example runs standalone and is meant to be read as much
as run.

## Conventions

Units, signs, and grid orientation are stated in [CLAUDE.md](CLAUDE.md) and are
assumed everywhere in the code. Read that before changing anything numerical.
