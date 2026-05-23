#!/usr/bin/env python3
"""
Normalize protein FASTA headers to >{ACC}|{gene_id}|{transcript_id}|{product}
Sources:
  - PlasmoDB:  >PVP01_0100100.1-p1 | transcript=PVP01_0100100.1 | gene=PVP01_0100100 | ... | gene_product=PIR protein | ...
  - NCBI:      >SCA59564.1 VIR protein [Plasmodium vivax]   (needs GFF cross-ref for gene_id)
  - Liftoff:   >PVP01_0006930.1   (gene_id = transcript - .NNN suffix; product unknown)
"""
import sys, re, os
from collections import defaultdict

def parse_plasmodb(line):
    # >PVP01_0100100.1-p1 | transcript=PVP01_0100100.1 | gene=PVP01_0100100 | gene_product=PIR protein
    parts = [p.strip() for p in line[1:].split('|')]
    h = parts[0].split()[0]
    trans, gene, prod = h, None, ""
    for p in parts[1:]:
        k, _, v = p.partition('=')
        if k == 'transcript': trans = v
        elif k == 'gene': gene = v
        elif k == 'gene_product': prod = v.replace(' ', '_')
    if gene is None: gene = trans.rsplit('.', 1)[0]
    return gene, trans, prod

def parse_liftoff(line):
    # >PVP01_0100100.1
    h = line[1:].split()[0]
    gene = h.rsplit('.', 1)[0]
    return gene, h, ""

def parse_ncbi(line, gff_protein2gene, gff_protein2product):
    # >SCA59564.1 VIR protein [Plasmodium vivax]
    parts = line[1:].split(None, 1)
    h = parts[0]
    rest = parts[1] if len(parts) > 1 else ""
    prod = rest.split('[')[0].strip().replace(' ', '_')
    gene = gff_protein2gene.get(h, h)
    if h in gff_protein2product and not prod:
        prod = gff_protein2product[h]
    return gene, h, prod

def load_ncbi_gff(gff_path):
    """Build protein_id → gene_id lookup from NCBI GFF.
    NCBI CDS line has: protein_id=XYZ.1; Parent=gene-LOCUS_TAG; product=...
    """
    p2g, p2prod = {}, {}
    with open(gff_path) as fh:
        for ln in fh:
            if ln.startswith('#') or '\t' not in ln: continue
            f = ln.split('\t')
            if len(f) < 9 or f[2] != 'CDS': continue
            attrs = dict(kv.split('=', 1) for kv in f[8].rstrip(';\n').split(';') if '=' in kv)
            prot = attrs.get('protein_id')
            if not prot: continue
            # Try locus_tag first; fall back to Parent's gene-LOCUS
            gene = attrs.get('locus_tag') or attrs.get('Name') or attrs.get('gene')
            if not gene:
                parent = attrs.get('Parent', '')
                # Parent could be "gene-PVT01_xxx" or "rna-xyz"
                if parent.startswith('gene-'): gene = parent[5:]
            if not gene: gene = prot.rsplit('.', 1)[0]
            p2g[prot] = gene
            if 'product' in attrs: p2prod[prot] = attrs['product'].replace(' ', '_')
    return p2g, p2prod

def main():
    acc, in_fa, out_fa, source, gff = sys.argv[1:6] if len(sys.argv) >= 6 else (*sys.argv[1:5], None)
    if source == 'ncbi':
        p2g, p2prod = load_ncbi_gff(gff)
    else:
        p2g, p2prod = {}, {}
    n_out = 0
    with open(in_fa) as fin, open(out_fa, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                if source == 'plasmodb':
                    g, t, p = parse_plasmodb(line)
                elif source == 'ncbi':
                    g, t, p = parse_ncbi(line, p2g, p2prod)
                else:  # liftoff
                    g, t, p = parse_liftoff(line)
                p = p or "unknown"
                fout.write(f">{acc}|{g}|{t}|{p}\n")
                n_out += 1
            else:
                fout.write(line)
    sys.stderr.write(f"{acc}: wrote {n_out} proteins\n")

if __name__ == '__main__':
    main()
