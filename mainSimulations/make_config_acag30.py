#!/usr/bin/env python3
"""
Generate a LAMMPS data file for (ACAG)30 x 64 at 200 uM
by remapping the atom types in an existing (CAG)40 x 64 file.

Usage (from the repository root or mainSimulations/):

    python mainSimulations/make_config_acag30.py \
        --input mainSimulations/config_cag40.dat \
        --output mainSimulations/config_acag30.dat
"""

import argparse
from pathlib import Path

def parse_num_atoms(lines):
    for line in lines:
        parts = line.split()
        if len(parts) == 2 and parts[1] == "atoms":
            return int(parts[0])
    raise ValueError("Could not find '<N> atoms' line in header.")

def find_atoms_block(lines, natoms):
    atoms_header_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Atoms"):
            atoms_header_idx = i
            break
    if atoms_header_idx is None:
        raise ValueError("Could not find 'Atoms' section in data file.")
    atoms_start = atoms_header_idx + 2
    atoms_end = atoms_start + natoms
    if atoms_end > len(lines):
        raise ValueError(
            f"Atoms block (expected {natoms} lines) extends past EOF "
            f"(start={atoms_start}, end={atoms_end}, total_lines={len(lines)})."
        )
    return atoms_start, atoms_end

def infer_chain_length_and_count(lines, atoms_start, atoms_end):
    n_atoms = atoms_end - atoms_start
    mol_ids = []
    for i in range(atoms_start, atoms_end):
        line = lines[i].strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        mol_ids.append(int(parts[1]))
    if not mol_ids:
        raise ValueError("No atom entries found in Atoms block.")
    n_chains = max(mol_ids)
    if n_atoms % n_chains != 0:
        raise ValueError(
            f"Total atoms ({n_atoms}) is not divisible by number of chains ({n_chains})."
        )
    chain_len = n_atoms // n_chains
    return chain_len, n_chains

def remap_to_acag30(lines, atoms_start, atoms_end, chain_len):
    base_pattern = [1, 2, 1, 3]  # (ACAG), A=1, C=2, G=3, U=4
    pattern_len = len(base_pattern)
    new_lines = list(lines)
    for i in range(atoms_start, atoms_end):
        raw = lines[i]
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = raw.split()
        if len(parts) < 6:
            raise ValueError(
                f"Unexpected atom line format at line {i+1}: {raw!r}"
            )
        atom_id = int(parts[0])
        idx_in_chain = (atom_id - 1) % chain_len
        new_type = base_pattern[idx_in_chain % pattern_len]
        parts[2] = str(new_type)
        new_lines[i] = " ".join(parts) + "\n"
    return new_lines

def main():
    parser = argparse.ArgumentParser(
        description="Generate config_acag30.dat from a (CAG)40-style LAMMPS data file."
    )
    parser.add_argument(
        "--input", "-i", type=Path, required=True,
        help="Input LAMMPS data file (e.g., config_cag40.dat).",
    )
    parser.add_argument(
        "--output", "-o", type=Path, required=True,
        help="Output LAMMPS data file (e.g., config_acag30.dat).",
    )
    args = parser.parse_args()
    lines = args.input.read_text().splitlines(keepends=True)
    natoms = parse_num_atoms(lines)
    atoms_start, atoms_end = find_atoms_block(lines, natoms)
    chain_len, n_chains = infer_chain_length_and_count(
        lines, atoms_start, atoms_end
    )
    print(
        f"Detected {natoms} atoms, {n_chains} chains, "
        f"{chain_len} nucleotides per chain."
    )
    new_lines = remap_to_acag30(lines, atoms_start, atoms_end, chain_len)
    args.output.write_text("".join(new_lines))
    print(f"Wrote remapped data file to: {args.output}")

if __name__ == "__main__":
    main()
