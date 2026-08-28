#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create an auditable Chinese Markdown report from RNA-seq pipeline outputs."""
import argparse, csv, json, math, os, re
from collections import Counter, defaultdict
from pathlib import Path

def tab(path):
    if not path.is_file() or not path.stat().st_size: return []
    with path.open(newline="") as h: return list(csv.DictReader(h, delimiter="\t"))
def n(x):
    try: return float(x)
    except (TypeError, ValueError): return float("nan")
def f(x, digits=2): return "—" if not math.isfinite(x) else f"{x:.{digits}f}"
def pct(x, digits=2): return "—" if not math.isfinite(x) else f"{x:.{digits}f}%"
def link(path, report): return Path(os.path.relpath(path, report.parent)).as_posix()
def table(head, body):
    return "" if not body else "\n".join(["| " + " | ".join(head) + " |", "|"+"|".join("---" for _ in head)+"|", *["| " + " | ".join(str(x) for x in r) + " |" for r in body]])
def fig(project, report, name, caption):
    p=project/"figures"/name
    return f"![{caption}]({link(p,report)})\n\n*图注：{caption.rstrip('。')}。*" if p.is_file() else ""
def samples(path):
    out=[]
    with Path(path).open() as h:
        if h.readline().split()[:2] != ["sample","group"]: raise SystemExit("Sample table needs sample and group columns")
        for line in h:
            z=line.split()
            if z: out.append((z[0],z[1]))
    return out
def mean(xs): return sum(xs)/len(xs) if xs else float("nan")
def pearson(a,b):
    x,y=mean(a),mean(b); den=math.sqrt(sum((v-x)**2 for v in a)*sum((v-y)**2 for v in b))
    return sum((u-x)*(v-y) for u,v in zip(a,b))/den if den else float("nan")
def nested(d,*keys):
    for k in keys:
        if not isinstance(d,dict): return float("nan")
        d=d.get(k,{})
    return d if isinstance(d,(int,float)) else float("nan")

def parse_fastp(path):
    d=json.loads(path.read_text()); before=d.get("summary",{}).get("before_filtering",{}); after=d.get("summary",{}).get("after_filtering",{}); fr=d.get("filtering_result",{}); ac=d.get("adapter_cutting",{})
    raw=nested(d,"summary","before_filtering","total_reads"); clean=nested(d,"summary","after_filtering","total_reads")
    return {"raw":raw,"raw_bp":nested(d,"summary","before_filtering","total_bases"),"clean":clean,"clean_bp":nested(d,"summary","after_filtering","total_bases"),"retain":100*clean/raw if raw else float("nan"),"q20_raw":100*nested(d,"summary","before_filtering","q20_rate"),"q30_raw":100*nested(d,"summary","before_filtering","q30_rate"),"q20":100*nested(d,"summary","after_filtering","q20_rate"),"q30":100*nested(d,"summary","after_filtering","q30_rate"),"gc_raw":100*nested(d,"summary","before_filtering","gc_content"),"gc":100*nested(d,"summary","after_filtering","gc_content"),"adapter":100*n(ac.get("adapter_trimmed_reads"))/raw if raw else float("nan"),"lowq":100*n(fr.get("low_quality_reads"))/raw if raw else float("nan"),"n_reads":100*n(fr.get("too_many_N_reads"))/raw if raw else float("nan"),"dup":100*nested(d,"duplication","rate")}
def parse_hisat2(path):
    t=path.read_text(); grab=lambda p: n(re.search(p,t,re.M).group(1)) if re.search(p,t,re.M) else float("nan")
    total=grab(r"^(\d+) reads;"); zero=grab(r"^\s*(\d+) \([^\n]+\) aligned concordantly 0 times"); unique=grab(r"^\s*(\d+) \([^\n]+\) aligned concordantly exactly 1 time"); multi=grab(r"^\s*(\d+) \([^\n]+\) aligned concordantly >1 times")
    return {"input":total,"overall":grab(r"([0-9.]+)% overall alignment rate"),"unique":100*unique/total if total else float("nan"),"multi":100*multi/total if total else float("nan"),"unmapped":100*zero/total if total else float("nan"),"concordant":100*(unique+multi)/total if total else float("nan")}
def count_matrix(path):
    if not path.is_file(): return None
    with path.open() as h:
        rd=csv.reader(h,delimiter="\t"); header=next(rd); total=bad=0
        for row in rd:
            total+=1
            if len(row)!=len(header) or not row[0] or any(not x.isdigit() for x in row[1:]): bad+=1
    return header,total,bad
def get_manifest(path): return {r["key"]:r["value"] for r in tab(path) if "key" in r and "value" in r}

def main():
 p=argparse.ArgumentParser()
 for key in ("project","samples","control","treatment","output","project-name","species","mode","config","reference","annotation","reference-version","strandedness","min-count","min-samples","fdr","lfc"): p.add_argument("--"+key, required=key in {"project","samples","control","treatment","output"})
 for key,default in (("warn-min-replicates",3),("warn-min-retention",70),("warn-min-q30",80),("warn-min-alignment",70),("warn-min-assignment",50),("warn-min-filtered-fraction",.05),("warn-max-depth-ratio",2),("warn-min-within-correlation",.8),("warn-min-mapping-rate",.3)): p.add_argument("--"+key,type=float,default=default)
 a=p.parse_args(); project=Path(a.project).resolve(); report=Path(a.output).resolve(); res=project/"results"; groups=samples(a.samples); manifest=get_manifest(res/"run_manifest.tsv"); warnings=[]; notices=[]
 bygroup=defaultdict(list)
 for s,g in groups: bygroup[g].append(s)
 for g,ss in bygroup.items():
  if len(ss)<a.warn_min_replicates: warnings.append(f"{g} 组仅有 {len(ss)} 个生物学重复，低于设置的预警线 {a.warn_min_replicates:g}；方差估计和差异检验的稳定性可能受限。")
 if a.control not in bygroup or a.treatment not in bygroup: warnings.append("配置的对照组或处理组未在样本表中找到；比较方向需要核查。")

 qc={}; align={}
 for s,g in groups:
  q=res/"clean_fastq"/"reports"/f"{s}_fastp.json"; h=res/"alignment"/f"{s}_alignment_summary.txt"
  if q.is_file(): qc[s]=parse_fastp(q)
  if h.is_file(): align[s]=parse_hisat2(h)
 if not qc: notices.append("本次为计数矩阵模式或未提供 FASTQ 质控文件，因此无法评价原始读段质量。")
 else:
  for s,x in qc.items():
   if x["retain"]<a.warn_min_retention: warnings.append(f"{s} 的 clean reads 保留率为 {pct(x['retain'])}，低于经验预警线 {a.warn_min_retention:g}%。")
   if x["q30"]<a.warn_min_q30: warnings.append(f"{s} 的 clean Q30 为 {pct(x['q30'])}，低于经验预警线 {a.warn_min_q30:g}%。")
  for metric,label in (("q30","Q30"),("gc","GC 含量"),("dup","重复率")):
   vals=[x[metric] for x in qc.values() if math.isfinite(x[metric])]
   if vals and max(vals)-min(vals)>20: warnings.append(f"样本间 {label} 极差为 {max(vals)-min(vals):.2f} 个百分点，存在样本间技术差异的经验性预警。")
 if align:
  for s,x in align.items():
   if x["overall"]<a.warn_min_alignment: warnings.append(f"{s} 的总体比对率为 {pct(x['overall'])}，低于经验预警线 {a.warn_min_alignment:g}%。")
  vals=[x["overall"] for x in align.values() if math.isfinite(x["overall"])]
  if vals and max(vals)-min(vals)>10: warnings.append(f"样本总体比对率极差为 {max(vals)-min(vals):.2f} 个百分点，建议核查参考匹配、污染和文库类型。")

 assign={}; unassigned={}; csum=res/"counts"/"gene_counts.txt.summary"
 if csum.is_file():
  raw=list(csv.reader(csum.open(),delimiter="\t")); names=[Path(x).name.replace(".sorted.bam","") for x in raw[0][1:]]; d={r[0]:[int(v) for v in r[1:]] for r in raw[1:]}
  for i,s in enumerate(names):
   total=sum(v[i] for v in d.values()); assign[s]=100*d.get("Assigned",[0]*len(names))[i]/total if total else float("nan"); reasons=sorted(((k,v[i]) for k,v in d.items() if k.startswith("Unassigned")),key=lambda x:x[1],reverse=True); unassigned[s]=reasons[0][0].replace("Unassigned_","") if reasons else "—"
   if assign[s]<a.warn_min_assignment: warnings.append(f"{s} 的基因分配率为 {pct(assign[s])}，低于经验预警线 {a.warn_min_assignment:g}%。")
 rawm=count_matrix(res/"counts"/"gene_counts_matrix.tsv"); filt=count_matrix(res/"gene_counts_matrix_filtered.tsv")
 rawgenes=rawm[1] if rawm else float("nan"); filtgenes=filt[1] if filt else float("nan"); frac=100*filtgenes/rawgenes if rawgenes else float("nan")
 if rawm and rawm[2]: warnings.append(f"原始计数矩阵有 {rawm[2]} 行格式或整数计数异常；差异分析输入应重新核查。")
 if rawm and [x for x in rawm[0][1:] if x not in {s for s,_ in groups}]: warnings.append("原始计数矩阵存在样本表之外的列名；请核查样本匹配。")
 if math.isfinite(frac) and frac/100<a.warn_min_filtered_fraction: warnings.append(f"低表达过滤后仅保留 {filtgenes}/{rawgenes} 个基因（{pct(frac)}），低于经验预警线 {pct(100*a.warn_min_filtered_fraction)}。已核查矩阵格式后，该结果仍应视为重要限制：请复核测序深度、过滤参数和计数分布。")

 corrows=tab(res/"sample_correlation_matrix.tsv"); corr={r["sample"]:{k:n(v) for k,v in r.items() if k!="sample"} for r in corrows}; within=[]; between=[]
 for s,g in groups:
  for t,h in groups:
   if s<t and s in corr and t in corr[s]: (within if g==h else between).append(corr[s][t])
 pca=tab(res/"pca_coordinates.tsv"); pc1=n(pca[0].get("PC1_variance_percent")) if pca else float("nan"); pc2=n(pca[0].get("PC2_variance_percent")) if pca else float("nan")
 if within and mean(within)<a.warn_min_within_correlation: warnings.append(f"组内 Pearson 相关性平均值为 {mean(within):.3f}，低于经验预警线 {a.warn_min_within_correlation:.2f}。")
 if within and between and mean(within)<=mean(between): warnings.append(f"组内平均相关性（{mean(within):.3f}）未高于组间平均相关性（{mean(between):.3f}），分组一致性不足。")
 sqc=tab(res/"sample_qc_metrics.tsv")
 if sqc:
  mix=defaultdict(set)
  for r in sqc: mix[r["cluster_k2"]].add(r["group"])
  if any(len(x)>1 for x in mix.values()): warnings.append("VST 欧氏距离的 k=2 聚类未与组别完全一致；建议结合 PCA、样本元数据和原始质控复核，不应仅据此删除样本。")

 de=tab(res/"deseq2_results.tsv"); cls=tab(res/"deseq2_genes_classified.tsv"); sig=tab(res/"deseq2_significant_genes.tsv"); up=[r for r in sig if r.get("regulation")=="Up-regulated"]; down=[r for r in sig if r.get("regulation")=="Down-regulated"]; nonmiss_p=[r for r in de if math.isfinite(n(r.get("pvalue")))]; nonmiss_q=[r for r in de if math.isfinite(n(r.get("padj")))]; padjsig=[r for r in de if n(r.get("padj"))<n(a.fdr)]
 if len(sig)<5: warnings.append(f"仅识别到 {len(sig)} 个 DEG；富集检验的统计功效可能不足。")
 if de and len(sig)/len(de)>0.5: warnings.append(f"DEG 占检验基因的 {pct(100*len(sig)/len(de))}，比例较高；建议检查样本分组、批次效应和模型设计。")

 esum=tab(res/"enrichment_summary.tsv"); enrich_by={(r["gene_set"],r["database"]):r for r in esum}; candidates={}
 for kind in ("GO","KEGG"):
  rr=[]
  for tag in ("all","up","down"): rr+=tab(res/f"enrichment_{kind}_{tag}.tsv")
  if rr: candidates[kind]=min(rr,key=lambda r:n(r["padj"]))
  for tag in ("all","up","down"):
   item=enrich_by.get((tag,kind))
   if item and n(item["input_genes"]) and n(item["mapped_input_genes"])/n(item["input_genes"])<a.warn_min_mapping_rate: warnings.append(f"{tag} DEG 的 {kind} 注释映射率为 {pct(100*n(item['mapped_input_genes'])/n(item['input_genes']))}，低于经验预警线 {pct(100*a.warn_min_mapping_rate)}。")
  if kind in candidates and n(candidates[kind]["padj"])>=n(a.fdr) and n(candidates[kind]["pvalue"])<.05: notices.append(f"{kind} 的最小原始 p 值条目未通过 FDR 校正，只能作为探索性线索，不能表述为显著富集。")

 missing=[]; required=[res/"deseq2_results.tsv",res/"deseq2_summary.txt",res/"gene_counts_matrix_filtered.tsv",project/"figures"/"PCA_plot.png",project/"figures"/"volcano_plot.png"]
 for x in required:
  if not x.exists(): missing.append(link(x,report))
 if missing: warnings.append("以下关键报告引用不存在："+"、".join(missing)+"。")
 title=a.project_name or manifest.get("project_name","RNA-seq项目"); species=a.species or manifest.get("species","未说明物种"); mode=a.mode or manifest.get("analysis_mode","未说明")
 qcrows=[]
 for s,g in groups:
  x=qc.get(s,{}); qcrows.append([s,g,f(x.get("raw"),0),f(x.get("raw_bp"),0),f(x.get("clean"),0),f(x.get("clean_bp"),0),pct(x.get("retain",float("nan"))),pct(x.get("q30_raw",float("nan"))),pct(x.get("q30",float("nan"))),pct(x.get("gc",float("nan"))),pct(x.get("adapter",float("nan"))),pct(x.get("lowq",float("nan"))),pct(x.get("n_reads",float("nan"))),pct(x.get("dup",float("nan")),2),link(res/"clean_fastq"/"reports"/f"{s}_fastp.html",report) if s in qc else "—"])
 alignrows=[[s,g,f(align.get(s,{}).get("input",float("nan")),0),pct(align.get(s,{}).get("overall",float("nan"))),pct(align.get(s,{}).get("unique",float("nan"))),pct(align.get(s,{}).get("multi",float("nan"))),pct(align.get(s,{}).get("concordant",float("nan"))),pct(align.get(s,{}).get("unmapped",float("nan"))),link(res/"alignment"/f"{s}.sorted.bam",report)] for s,g in groups]
 qrows=[[s,g,pct(assign.get(s,float("nan"))),unassigned.get(s,"—")] for s,g in groups]
 top=sorted(sig,key=lambda r:n(r["padj"]))[:10]
 files=[("分析配置",Path(a.config) if a.config else None,"完整参数与输入路径"),("运行清单",res/"run_manifest.tsv","时间、软件版本、参考和关键参数"),("原始计数矩阵",res/"counts"/"gene_counts_matrix.tsv","featureCounts gene-level raw integer counts"),("DESeq2完整结果",res/"deseq2_results.tsv","所有检验基因的统计结果"),("DEG清单",res/"deseq2_significant_genes.tsv","同时满足 FDR 与 fold-change 阈值的基因"),("富集汇总",res/"enrichment_summary.tsv","GO/KEGG 映射、背景与显著条目统计")]
 filetable=[[name,link(p.resolve(),report) if p and p.exists() else "未生成",what] for name,p,what in files]
 qinterp=(f"clean Q30 范围 {f(min(x['q30'] for x in qc.values()),2)}%–{f(max(x['q30'] for x in qc.values()),2)}%，clean GC 范围 {f(min(x['gc'] for x in qc.values()),2)}%–{f(max(x['gc'] for x in qc.values()),2)}%，保留率范围 {f(min(x['retain'] for x in qc.values()),2)}%–{f(max(x['retain'] for x in qc.values()),2)}%。自动异常判断见报告开头。" if qc else "不适用。")
 ainterp=(f"总体比对率范围 {f(min(x['overall'] for x in align.values()),2)}%–{f(max(x['overall'] for x in align.values()),2)}%。样本间和组间比较的自动预警见报告开头。" if align else "当前结果不足以判断。")
 versions=[[k.replace('software_',''),v] for k,v in manifest.items() if k.startswith('software_')]
 qinterp=(f"clean Q30 范围 {f(min(x['q30'] for x in qc.values()),2)}%–{f(max(x['q30'] for x in qc.values()),2)}%，clean GC 范围 {f(min(x['gc'] for x in qc.values()),2)}%–{f(max(x['gc'] for x in qc.values()),2)}%，保留率范围 {f(min(x['retain'] for x in qc.values()),2)}%–{f(max(x['retain'] for x in qc.values()),2)}%。" if qc else "不适用。")
 ainterp=(f"总体比对率范围 {f(min(x['overall'] for x in align.values()),2)}%–{f(max(x['overall'] for x in align.values()),2)}%。" if align else "当前结果不足以判断。")
 report_lines=[f"# {title}：RNA-seq 分析报告","",f"生成日期：{manifest.get('run_date','未记录')}。报告全部数值来自结果文件；经验阈值只用于预警，不构成硬性判定。","","## 质量警告与重要发现","",*( [f"> **质量警告：** {x}" for x in warnings] if warnings else ["> 未触发配置的自动质量预警；仍需人工核查实验设计与元数据。"] ),*( [f"> **重要说明：** {x}" for x in notices] if notices else []),"","## 1. 项目概览","",f"项目：{title}；物种：{species}；模式：{mode}。比较为 **{a.treatment} vs {a.control}**；正 log2FoldChange 表示 {a.treatment} 组高于 {a.control} 组。样本表仅提供 `{a.control}` / `{a.treatment}` 标签，若无实验元数据，当前结果不足以推断处理机制。","",f"参考基因组：`{a.reference or manifest.get('reference_fasta','未记录')}`；注释：`{a.annotation or manifest.get('annotation_gtf','未记录')}`；版本/来源：{a.reference_version or manifest.get('reference_version','未记录')}。","",table(["样本","分组"],groups),"",f"核心结果摘要：共 {len(groups)} 个样本；过滤后 {filtgenes if math.isfinite(filtgenes) else '—'} 个基因进入 DESeq2，FDR 显著基因 {len(padjsig)} 个，DEG {len(sig)} 个（上调 {len(up)}、下调 {len(down)}）。"+(f" PCA 的 PC1/PC2 解释 {f(pc1,1)}%/{f(pc2,1)}% 方差。" if math.isfinite(pc1) else "")+" 质量限制和富集结论见后续章节。","","## 2. 原始数据与测序质量控制","",("未执行 FASTQ 质控：本次从计数矩阵开始分析，当前结果不足以判断原始读段质量。" if not qc else table(["样本","组","raw reads","raw bases","clean reads","clean bases","保留率","raw Q30","clean Q30","clean GC","接头修剪率","低质率","N异常率","重复率","fastp报告"],qcrows)),"",f"质量解读：{qinterp} 自动异常判断见报告开头；fastp JSON/HTML 位于 `results/clean_fastq/reports/`。","","## 3. 参考基因组比对","",f"HISAT2 用于{'双端' if mode=='full' else '已完成计数矩阵模式，不含本步'}比对；排序 BAM/BAI 与摘要在 `results/alignment/`。一致配对率为 concordantly aligned 的比例。", "",table(["样本","组","输入pairs","总体比对率","唯一比对率","多重比对率","一致配对率","未一致比对率","BAM"],alignrows) if align else "未提供 HISAT2 摘要，当前结果不足以判断比对质量。","",f"比对解读：{ainterp} 样本/组间异常和可能的参考匹配、污染或文库类型问题在质量警告中提示。","","## 4. 基因表达定量","",f"featureCounts 设置：feature=`exon`，聚合字段=`gene_id`，双端按 fragment 计数（`--countReadPairs`），链特异性 `-s {a.strandedness or manifest.get('strandedness','未记录')}`；默认不计多重比对或多重重叠 reads。原始矩阵 {rawgenes if math.isfinite(rawgenes) else '—'} 基因，过滤规则为 count ≥ {a.min_count or manifest.get('min_count','未记录')} 且至少 {a.min_samples or manifest.get('min_samples','未记录')} 个样本满足；保留 {filtgenes if math.isfinite(filtgenes) else '—'} 个（{pct(frac)}）。低表达过滤可减少低信息基因对方差估计和多重检验的影响。","",table(["样本","组","基因分配率","主要未分配原因"],qrows) if assign else "未提供 featureCounts 汇总，当前结果不足以判断基因分配质量。","","## 5. 表达数据预处理与样本质量评价","",f"PCA、Pearson 相关性和层次聚类使用 DESeq2 的 blind VST 数据，而非 raw counts；VST 用于减弱均值—方差关系。相关性为 Pearson r；聚类为 VST 欧氏距离、complete linkage。"+(f" PC1/PC2 解释 {f(pc1,1)}%/{f(pc2,1)}% 方差。" if math.isfinite(pc1) else ""),"",f"组内 Pearson r：均值 {f(mean(within),3)}，范围 {f(min(within),3)}–{f(max(within),3)}；组间均值 {f(mean(between),3)}。" if within and between else "相关性矩阵不可用，当前结果不足以判断重复一致性。","",table(["样本","分组","PC1","PC2"],[[r['sample'],r['group'],f(n(r['PC1']),2),f(n(r['PC2']),2)] for r in pca]),"",fig(project,report,"PCA_plot.png","PCA 图：VST 表达矩阵；点颜色表示组别；PC1/PC2 解释率见坐标轴。"),"",fig(project,report,"correlation_heatmap.png","样本 Pearson 相关性热图：基于 VST；颜色表示相关系数，树状图按相关性模式聚类。"),"",fig(project,report,"hierarchical_clustering.png","样本层次聚类图：基于 VST 欧氏距离、complete linkage；用于检查聚类与设计分组是否一致。"),"","## 6. 差异表达分析","",f"DESeq2 基于过滤后的原始整数 counts，设计式 `~ group`，比较 {a.treatment} vs {a.control}，Wald 检验，Benjamini–Hochberg 校正；未实施 log2FoldChange shrinkage。阈值：padj < {a.fdr or manifest.get('fdr','未记录')} 且 |log2FoldChange| ≥ {a.lfc or manifest.get('lfc_cutoff','未记录')}。","",f"统计检验基因 {len(de)} 个；pvalue 非缺失 {len(nonmiss_p)} 个；padj 非缺失 {len(nonmiss_q)} 个；仅满足 padj 阈值 {len(padjsig)} 个；DEG {len(sig)} 个（上调 {len(up)}、下调 {len(down)}）。低 baseMean 而高 fold-change 的基因应与其 baseMean 一并审阅。","",table(["Gene ID","方向","baseMean","log2FC","pvalue","padj"],[[r['Geneid'],r.get('regulation','—'),f(n(r.get('baseMean'))),f(n(r.get('log2FoldChange'))),f"{n(r.get('pvalue')):.2e}",f"{n(r.get('padj')):.2e}"] for r in top]),"",fig(project,report,"MA_plot.png","MA 图：所有 DESeq2 检验基因；横轴 log10(baseMean+1)，纵轴 log2FC；颜色为 DEG 方向，虚线为 |log2FC| 阈值。"),"",fig(project,report,"volcano_plot.png","火山图：所有 DESeq2 检验基因；横轴 log2FC，纵轴 −log10(padj)；颜色为 DEG 方向，虚线标识 FDR 与 fold-change 阈值。"),"",fig(project,report,"heatmap_deg.png",f"DEG 热图：{len(sig)} 个 DEG 的 VST 值按基因行 Z-score 标准化；行/列均层次聚类，展示相对表达模式。"),"","## 7. 功能富集分析","","富集背景为实际进入 DESeq2 且有对应 GO/KEGG 注释的基因；结果分别针对全部、上调和下调 DEG。FDR 显著与仅原始 p 值较小的候选严格区分。","",table(["基因集","数据库","输入DEG","已映射","映射背景","检验条目","FDR显著条目"],[[r['gene_set'],r['database'],r['input_genes'],r['mapped_input_genes'],r['mapped_background_genes'],r['tested_terms'],r['fdr_significant_terms']] for r in esum]),"",table(["数据库","最小 padj 条目","padj","命中基因"],[[k,v.get('Term',v.get('Pathway','—')),v['padj'],v.get('genes','—')] for k,v in candidates.items()]),"",fig(project,report,"GO_enrichment_barplot.png","GO 富集柱状图：横轴为 −log10(pvalue)，条形颜色为 BP/MF/CC；星号仅表示 padj < 0.05。"),"",fig(project,report,"GO_enrichment_bubbleplot.png","GO 富集气泡图：横轴为 GeneRatio，点大小为命中 DEG 数，颜色为 −log10(padj)。"),"",fig(project,report,"KEGG_enrichment_barplot.png","KEGG 富集柱状图：横轴为 −log10(pvalue)，颜色表示上/下调方向；FDR 显著性以星号标注。"),"",fig(project,report,"KEGG_enrichment_bubbleplot.png","KEGG 富集气泡图：横轴为 GeneRatio，点大小为命中 DEG 数，颜色为 −log10(padj)。"),"","## 8. 结果整合与生物学解读","",f"可由当前数据直接支持的结果是：在当前参考、计数、过滤和统计阈值下，{a.treatment} 与 {a.control} 之间检测到 {len(sig)} 个 DEG。样本质量与分组表现应结合第 2–5 节量化指标判断。只有 FDR 显著的富集条目可作为统计支持；未通过 FDR 的条目仅为待验证线索。若没有实验背景或可靠功能注释，当前结果不足以推断具体处理机制。","","## 9. 局限性、异常与后续建议","",*( [f"- {x}" for x in warnings] if warnings else ["- 未触发自动预警；仍建议核对批次信息、实验元数据和原始 fastp/HISAT2 报告。"] ),"- 不应仅依据 PCA 或聚类图删除样本；如出现不一致，应先核查样本身份、批次、文库构建和测序深度。","- 建议对优先候选 DEG 进行独立重复或 qRT-PCR 验证。","","## 10. 文件清单与可重复性信息","",table(["文件","相对路径","内容/推荐用途"],filetable),"",table(["软件","版本/路径"],versions),"",f"关键参数：threads={manifest.get('threads','未记录')}；strandedness={a.strandedness or manifest.get('strandedness','未记录')}；MIN_COUNT={a.min_count or manifest.get('min_count','未记录')}；MIN_SAMPLES={a.min_samples or manifest.get('min_samples','未记录')}；FDR={a.fdr or manifest.get('fdr','未记录')}；LFC cutoff={a.lfc or manifest.get('lfc_cutoff','未记录')}。软件版本、运行时间、参考路径和配置文件记录在 `{link(res/'run_manifest.tsv',report)}`。本流程未使用随机抽样或随机种子。"]
 report.write_text("\n".join(x for x in report_lines if x is not None),encoding="utf-8"); print(f"Report written: {report}")
if __name__=="__main__": main()
