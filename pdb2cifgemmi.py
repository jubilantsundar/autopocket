"""
Convert PDB to mmCIF using Gemmi, with Chimera-like setup
References: Boltz community discussion (slack)"""
#!/usr/bin/env python3
import argparse, re
from pathlib import Path
import gemmi

def _sanitize_chain_names(st: gemmi.Structure):
    for model in st:
        for ch in model:
            name = (ch.name or "").strip()
            if not name or name in {".", "?"}:
                ch.name = f"X{model.index}-{ch.index}"
            else:
                ch.name = re.sub(r"[^A-Za-z0-9_-]", "_", name)

def _populate_entity_sequences_from_model(st: gemmi.Structure):
    """
    For each polymer entity, build a 3-letter-code sequence from the longest
    polymer subchain in model 0 and store it in Entity.full_sequence.
    Gemmi will then write _entity_poly/_entity_poly_seq to mmCIF.

    This implementation is defensive about Gemmi versions where
    Chain.subchains() may yield either Subchain objects (with .get_polymer())
    or ResidueSpan directly.
    """
    if len(st) == 0:
        return
    model0: gemmi.Model = st[0]

    # Build an index: subchain_id -> polymer ResidueSpan
    subchain_by_id: dict[str, gemmi.ResidueSpan] = {}
    for ch in model0:
        for sc in ch.subchains():
            # Normalize: obtain a ResidueSpan for the polymer portion
            if hasattr(sc, "get_polymer"):
                # sc is a Subchain
                poly = sc.get_polymer()
                sid = sc.subchain_id() if hasattr(sc, "subchain_id") else None
            else:
                # sc is already a ResidueSpan (older/newer Gemmi variants)
                poly = sc
                # Try to obtain subchain id from the span itself
                sid = sc.subchain_id() if hasattr(sc, "subchain_id") else None

            if sid is None:
                # Fallback: derive an ID from chain name and first/last positions
                try:
                    chain_name = ch.name
                except Exception:
                    chain_name = "?"
                sid = f"{chain_name}:{poly.first().seqid}" if len(poly) else f"{chain_name}:empty"

            if len(poly) == 0:
                continue
            subchain_by_id[str(sid)] = poly

    for ent in st.entities:
        # only polymers participate in _entity_poly/_entity_poly_seq
        if ent.entity_type != gemmi.EntityType.Polymer:
            continue

        candidate_seqs: list[list[str]] = []
        for sid in ent.subchains:
            poly = subchain_by_id.get(str(sid))
            if poly is None or len(poly) == 0:
                continue
            # 3-letter monomer IDs in coordinate order
            seq = [res.name for res in poly]
            candidate_seqs.append(seq)

        if not candidate_seqs:
            # no coordinates for this entity → leave sequence empty
            continue

        # pick the longest as the representative entity sequence
        best = max(candidate_seqs, key=len)
        ent.full_sequence = best  # list[str] of 3-letter mon_ids

def _chimera_like_setup(st: gemmi.Structure):
    """
    Make mmCIF-compatible like Chimera exports:
      - non-empty chain/subchain ids
      - entities rebuilt and deduplicated
      - label sequence numbers present
      - entity polymer sequences populated so _entity_poly_seq is written
    """
    _sanitize_chain_names(st)
    st.assign_subchains(force=True, fail_if_unknown=False)
    st.setup_entities()
    st.ensure_entities()
    st.add_entity_types(overwrite=True)
    st.add_entity_ids(overwrite=True)
    st.deduplicate_entities()

    # NEW: infer and set polymer sequences for entities
    _populate_entity_sequences_from_model(st)

    # Ensure residues have label_seq_id after the entity sequences are known
    st.assign_label_seq_id(force=True)

def convert_one(inp: Path, out: Path):
    st = gemmi.read_structure(str(inp))
    _chimera_like_setup(st)

    doc = st.make_mmcif_document()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(doc.as_string())
    print(f"Wrote {out}")

def main():
    p = argparse.ArgumentParser(description="PDB → mmCIF (Gemmi, Chimera-like)")
    p.add_argument("inputs", nargs="+")
    p.add_argument("-o", "--outdir", default="cif_out")
    args = p.parse_args()

    outdir = Path(args.outdir)
    files = []
    for item in args.inputs:
        pth = Path(item)
        files += sorted(pth.glob("*.pdb")) if pth.is_dir() else [pth]

    if not files:
        raise SystemExit("No PDB inputs found.")
    for f in files:
        convert_one(f, outdir / (f.stem + ".cif"))

if __name__ == "__main__":
    main()