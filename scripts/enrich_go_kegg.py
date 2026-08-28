#!/usr/bin/env python3
"""Local GO/KEGG enrichment for DEG tables emitted by deseq2_report.R."""
import argparse, csv, math, re
from collections import defaultdict
from pathlib import Path

ONTOLOGIES = {"biological_process": "BP", "molecular_function": "MF", "cellular_component": "CC"}

def normalize_gene_id(identifier):
    """Match locus IDs between count tables and annotation tables when a GFF prefix is present."""
    match = re.search(r"(Glyma\.[^.\s]+)$", identifier.strip())
    return match.group(1) if match else identifier.strip()

def read_ids(path):
    with open(path, newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        if "Geneid" not in (rd.fieldnames or []): raise SystemExit(f"Geneid column missing in {path}")
        return [normalize_gene_id(r["Geneid"]) for r in rd]

def bh(p):
    n=len(p); order=sorted(range(n),key=p.__getitem__); q=[0.0]*n; prev=1.0
    for rank, i in reversed(list(enumerate(order,1))):
        prev=min(prev, p[i]*n/rank); q[i]=min(1.0,prev)
    return q

def logcomb(n,k):
    if k<0 or k>n: return float("-inf")
    return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)

def hypergeom_sf(kminus1,N,K,n):
    # P(X >= k), evaluated stably enough for ordinary RNA-seq gene sets.
    start=max(0,n-(N-K)); end=min(n,K); terms=[]
    denom=logcomb(N,n)
    for x in range(max(kminus1+1,start),end+1): terms.append(math.exp(logcomb(K,x)+logcomb(N-K,n-x)-denom))
    return min(1.0,sum(terms))

def parse_obo(path):
    data={}; term=None; name=ns=""
    for raw in open(path):
        line=raw.rstrip("\n")
        if line=="[Term]":
            if term and name: data[term]=(name,ns)
            term=None; name=ns=""
        elif line.startswith("id: GO:"): term=line[4:].strip()
        elif line.startswith("name: "): name=line[6:].strip()
        elif line.startswith("namespace: "): ns=line[11:].strip()
    if term and name: data[term]=(name,ns)
    return data

def enrich(foreground, background, mapping, labels, ontology=None):
    bg=[g for g in background if mapping.get(g)]; fg=[g for g in foreground if mapping.get(g)]
    N,n=len(bg),len(fg); term_bg=defaultdict(set)
    for g in bg:
        for term in mapping[g]: term_bg[term].add(g)
    rows=[]
    for term, members in term_bg.items():
        hit=sorted(g for g in fg if term in mapping[g]); k=len(hit)
        if k: rows.append([ontology or "",term,labels.get(term,""),f"{k}/{n}",f"{len(members)}/{N}",hypergeom_sf(k-1,N,len(members),n),"",",".join(hit)])
    q=bh([r[5] for r in rows])
    for r,x in zip(rows,q): r[6]=x
    return sorted(rows,key=lambda r:r[5])

def write(path, header, rows, drop_ontology=False):
    with open(path,"w",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n"); w.writerow(header)
        for r in rows:
            if drop_ontology:
                values = r[1:]
                w.writerow([*values[:4], f"{values[4]:.6g}", f"{values[5]:.6g}", values[6]])
            else:
                w.writerow([*r[:5], f"{r[5]:.6g}", f"{r[6]:.6g}", r[7]])

def main():
    p=argparse.ArgumentParser()
    for key in ("annotation","go_obo","ko_pathway","pathway_names","background","up","down","outdir"): p.add_argument("--"+key.replace("_","-"),required=True)
    p.add_argument("--gene-col",type=int,required=True); p.add_argument("--ko-col",type=int,required=True); p.add_argument("--go-col",type=int,required=True); a=p.parse_args()
    gene_go=defaultdict(set); gene_ko=defaultdict(set)
    for row in csv.reader(open(a.annotation),delimiter="\t"):
        if not row or row[0].startswith("#") or max(a.gene_col,a.ko_col,a.go_col)>=len(row): continue
        gene=row[a.gene_col].strip(); gene_ko[gene].update(x.strip().replace("ko:","") for x in row[a.ko_col].split(",") if x.strip()); gene_go[gene].update(x.strip() for x in row[a.go_col].split(",") if x.strip())
    goinfo=parse_obo(a.go_obo); ko_path=defaultdict(set)
    for row in csv.reader(open(a.ko_pathway),delimiter="\t"):
        if len(row)>=2 and row[1].replace("path:","").startswith("ko"): ko_path[row[0].replace("ko:","").strip()].add(row[1].replace("path:","").strip())
    pathnames={r[0].strip():r[1].strip() for r in csv.reader(open(a.pathway_names),delimiter="\t") if len(r)>=2}
    gene_path={g:{path for ko in ks for path in ko_path[ko]} for g,ks in gene_ko.items()}
    bg,up,dn=read_ids(a.background),read_ids(a.up),read_ids(a.down); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    summary=[]
    for foreground,tag in ((sorted(set(up+dn)),"all"),(up,"up"),(dn,"down")):
        allgo=[]
        for ont,short in ONTOLOGIES.items():
            mapping={g:{t for t in ts if goinfo.get(t,("", ""))[1]==ont} for g,ts in gene_go.items()}
            allgo.extend(enrich(foreground,bg,mapping,{t:v[0] for t,v in goinfo.items()},short))
        write(out/f"enrichment_GO_{tag}.tsv",["ontology","GO_ID","Term","GeneRatio","BgRatio","pvalue","padj","genes"],sorted(allgo,key=lambda r:r[5]))
        rows=enrich(foreground,bg,gene_path,pathnames)
        write(out/f"enrichment_KEGG_{tag}.tsv",["pathway_ID","Pathway","GeneRatio","BgRatio","pvalue","padj","genes"],rows,drop_ontology=True)
        summary.append([tag,"GO",len(foreground),sum(1 for g in foreground if gene_go.get(g)),sum(1 for g in bg if gene_go.get(g)),len(allgo),sum(1 for r in allgo if r[6] < 0.05)])
        summary.append([tag,"KEGG",len(foreground),sum(1 for g in foreground if gene_path.get(g)),sum(1 for g in bg if gene_path.get(g)),len(rows),sum(1 for r in rows if r[6] < 0.05)])
    with open(out/"enrichment_summary.tsv","w",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n"); w.writerow(["gene_set","database","input_genes","mapped_input_genes","mapped_background_genes","tested_terms","fdr_significant_terms"]); w.writerows(summary)
    print(f"Enrichment completed for {len(up)} upregulated and {len(dn)} downregulated genes.")
if __name__=="__main__": main()
