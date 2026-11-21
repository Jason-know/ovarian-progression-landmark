"""Create the Chinese peer-review report with verified analyses and source links."""
from pathlib import Path
import csv

from docx import Document
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

HERE = Path(__file__).resolve().parent
doc = Document()


def setup():
    s = doc.sections[0]
    s.page_width, s.page_height = Inches(8.27), Inches(11.69)
    s.top_margin, s.bottom_margin = Inches(.72), Inches(.68)
    s.left_margin = s.right_margin = Inches(.77)
    for name, size, bold in [("Normal", 10.5, False), ("Title", 29, True), ("Subtitle", 12, False),
                              ("Heading 1", 19, True), ("Heading 2", 12, True), ("Table Text", 9.5, False)]:
        styles = doc.styles
        style = styles[name] if name in styles else styles.add_style(name, WD_STYLE_TYPE.PARAGRAPH)
        style.font.name = "Calibri"
        style.font.size, style.font.bold = Pt(size), bold
        style.font.color.rgb = RGBColor.from_string("11151B" if name != "Subtitle" else "4C5870")
        rfonts = style.element.get_or_add_rPr().rFonts
        rfonts.set(qn("w:eastAsia"), "Microsoft YaHei")
        pf = style.paragraph_format
        pf.space_before, pf.space_after = Pt(0), Pt(8)
        pf.line_spacing = 1.16
        pp = style.element.get_or_add_pPr()
        for border in pp.findall(qn("w:pBdr")):
            pp.remove(border)
    doc.styles["Heading 2"].paragraph_format.space_before = Pt(8)
    doc.styles["Table Text"].paragraph_format.space_after = Pt(0)
    doc.styles["Subtitle"].font.italic = False
    footer = s.footer.paragraphs[0]
    footer.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    f = OxmlElement("w:fldSimple")
    f.set(qn("w:instr"), "PAGE")
    footer._p.append(f)
    doc.core_properties.title = "晚期浆液性卵巢癌稿件同行预审与病理补充"
    doc.core_properties.author = ""


def h(text, level=2):
    doc.add_heading(text, level)


def para(text, links=()):
    p = doc.add_paragraph(text)
    for label, url in links:
        p.add_run("  ")
        el = OxmlElement("w:hyperlink")
        rel = "http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink"
        el.set(qn("r:id"), p.part.relate_to(url, rel, is_external=True))
        run = OxmlElement("w:r")
        pr = OxmlElement("w:rPr")
        col = OxmlElement("w:color")
        col.set(qn("w:val"), "254F88")
        pr.append(col)
        run.append(pr)
        t = OxmlElement("w:t")
        t.text = label
        run.append(t)
        el.append(run)
        p._p.append(el)
    return p


def table(headers, rows, widths):
    t = doc.add_table(rows=0, cols=len(headers))
    t.autofit = False
    for c, width in zip(t.columns, widths, strict=True):
        c.width = Inches(width)
    for i, values in enumerate([headers, *rows]):
        row = t.add_row()
        row._tr.get_or_add_trPr().append(OxmlElement("w:cantSplit"))
        for j, (cell, value, width) in enumerate(zip(row.cells, values, widths, strict=True)):
            cell.width = Inches(width)
            p = cell.paragraphs[0]
            p.style = "Table Text"
            p.paragraph_format.space_before = Pt(4)
            p.paragraph_format.space_after = Pt(4)
            p.add_run(str(value)).bold = i == 0
            if i % 2 == 0:
                shade = OxmlElement("w:shd")
                shade.set(qn("w:fill"), ("E4EDF6" if i == 0 else "F3F6F9"))
                cell._tc.get_or_add_tcPr().append(shade)
    doc.add_paragraph().paragraph_format.space_after = Pt(0)


def rows(name):
    with (HERE / "results" / (name + ".csv")).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def ci(row, term, digits=3, signed=False):
    spec = f"{'+' if signed else ''}.{digits}f"
    out = f"{float(row[term]):{spec}} [{float(row[term+'_lo']):.{digits}f}, {float(row[term+'_hi']):.{digits}f}]"
    return out.replace("-", "−")


def review():
    setup()
    doc.add_paragraph("晚期浆液性卵巢癌稿件\n同行预审与病理补充", style="Title")
    doc.add_paragraph("2026 年 11 月 20 日 · 三个独立 AI 审稿视角与补充分析", style="Subtitle")
    h("编辑决定建议大修")
    para("固定模型移至外部队列后，进展状态改善了风险排序，绝对风险估计仍有偏差。故事应围绕两个问题展开。临床医生已经知道患者是否进展以后，精确到进展日期能否再改变排序？能否将一个队列得到的绝对风险用于另一个队列？三个审稿视角均建议大修。")
    para("本次属投稿前模拟预审。第 5 页列有三位真实专家候选及其机构资料，未向候选人发送稿件，也未取得他们的个人审稿意见。")
    h("查新结论")
    para("上一版查新漏引 Su 等于 2026 年 9 月发表的《Cancer》论文。该论文以一线铂类治疗起始为时间零点，在 487 名 BRCA 野生型晚期上皮性卵巢癌患者中使用 9、12、15 和 18 个月地标，检验进展状态与后续总生存的关联。稿件应撤去任何暗示首次证明 12 个月进展状态具有预后价值的表述。", [("Su 等原文", "https://pubmed.ncbi.nlm.nih.gov/42804616/")])
    table(["最接近的既有工作", "与本稿的实质关系"], [
        ["Su 等，Cancer 2026", "直接重合于地标进展状态与后续 OS。本稿仍可比较时间项的增量和锁定绝对风险的外部校准。"],
        ["Rose 等，2019", "复发后列线图已纳入复发时间和外部验证；以复发时为起点，适用于已复发者。"],
        ["Rosanu 等，2025", "单中心会议摘要已用早期复发动态预测长期生存。本稿的窗口和外部验证设计不同。"],
        ["Szamreta 等，2024；GCIG，2025", "既往条件生存及铂耐药人群地标研究均涉及进展状态，临床时点和估计目标有差异。"],
    ], [2.02, 4.68])
    para("截至本次定向检索，尚未检得在同一固定风险集用锁定方程同时比较进展状态与具体时间，并报告多队列固定时窗区分度及绝对风险校准的相同研究。这个范围内的差异可写入引言。检索结果不足以证明研究唯一。", [
        ("Rose 等", "https://pubmed.ncbi.nlm.nih.gov/30633128/"),
        ("Rosanu 等", "https://doi.org/10.1016/j.esmorw.2025.100376"),
        ("Szamreta 等", "https://pubmed.ncbi.nlm.nih.gov/38211405/"),
        ("GCIG 摘要", "https://doi.org/10.1016/j.annonc.2025.10.1041"),
    ])

    doc.add_page_break()
    h("审稿人一 妇科肿瘤临床", 1)
    para("建议大修。该题有临床可读性，尤其适合说明历史术后信息为何需要随病程更新。最需修正的是应用场景、时间口径和进展状态已给临床带来的信息量。")
    h("把研究问题落实到一次随访")
    para("开篇可比较两位年龄、分期和术后残余相近的患者。到研究时间零点后第 12 个月，一位仍无明确进展，另一位已进展。第一个问题是知道这个状态后能增加多少风险区分。第二个问题是进展发生在第 3 个月或第 9 个月会不会进一步改变未来一年风险判断。将这一比较写入摘要目的，主文先报告实际分层风险，再报告预测和校准。")
    h("核对进展事件及治疗背景")
    para("TCGA PFI 允许复发、新肿瘤事件及伴瘤死亡；ICON7 公开子集从随机化计时，原试验在第 9 和第 12 个月安排评估，影像或临床标准与 TCGA 不同。应主表列明队列零点、进展定义、评估频率、删失规则和地标前事件数。将 TCGA 明确复发或进展标签的敏感性移至主结果，增加边界日期窗口分析。地标时间也不能解释为末次铂类用药后的无铂间期。", [
        ("TCGA-CDR", "https://pmc.ncbi.nlm.nih.gov/articles/PMC6066282/"),
        ("ICON7", "https://pmc.ncbi.nlm.nih.gov/articles/PMC4648090/"),
    ])
    para("ICON7 验证样本来自德国 AGO-OVAR11 组织子集。稿件须解释 GEO 380 条记录与发表的分子分析 359 例之间的区别，以及本分析 224 人的筛选。UKE 也是德国队列，缺共享患者键，无法逐人核验跨库重叠。补上分期、治疗年代、贝伐珠单抗、现代 BRCA/HRD 与 PARP 维持资料的可用性，已知缺失明确写入局限。", [
        ("AGO-OVAR11 组织研究", "https://pubmed.ncbi.nlm.nih.gov/28159814/"),
        ("GEO 数据记录", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140082"),
    ])
    h("明确现在可用的临床价值")
    para("当前结果支持把进展状态作为预后研究及临床试验分层的候选变量。三个主要外部队列的绝对风险点估计均提示低估，其中 UKE 实际估计 20.1%、模型预测 9.5%。讨论中应明确，临床决策曲线仍是阈值模拟，不能据此缩短扫描间隔、启动治疗或告知个体死亡概率。若拟作为临床工具，需在现代治疗人群先进行独立校准和前瞻性评价。", [
        ("ASCO 复发指南", "https://doi.org/10.1200/JCO-26-01591"),
        ("早治疗随机试验", "https://doi.org/10.1016/S0140-6736(10)61268-8"),
    ])

    doc.add_page_break()
    h("审稿人二 临床预测与生存统计", 1)
    para("建议大修。锁定 TCGA 方程和带删失评估是优势。值得保留的统计比较是同一队列、同一风险集、同一结局上的四步基准。")
    h("把单独进展状态列入主结果")
    b = rows("clinical_benchmark")
    table(["外部队列", "仅状态 AUC", "临床变量加状态 AUC", "ΔAUC 与 95% CI", "ΔBrier 与 95% CI"],
          [[r["cohort"], f'{float(r["auc_reference"]):.3f}', f'{float(r["auc_added"]):.3f}',
            ci(r, "delta_auc", signed=True), ci(r, "delta_brier", digits=4, signed=True)] for r in b],
          [.77, .95, 1.13, 1.77, 2.08])
    para("此为预审后的补做对照，外部患者配对自助抽样 1,000 次。只有 ICON7 对加入年龄、分期及残余病灶的 AUC 增益有较明确证据。三个队列的 ΔBrier 区间均跨零，UKE 点估计反而变差。主文应完整呈现临床基线、单独进展状态、临床基线加状态、再加时长四种比较。")
    h("把已进展患者单独评估")
    para("进展时长只改变已进展者的排序。在该组中，AOCS 43 人的 AUC 为 0.548→0.520，UKE 23 人为 0.398→0.380，ICON7 52 人为 0.636→0.648，三组增量区间跨零。分别按进展状态估计删失权重后，三个主要外部的 Status AUC 为 0.860、0.761、0.853，与原分析接近。固定方程的区间只反映评价队列抽样，不能据此建立时间项在其他人群或编码上的等效性。")
    h("把绝对风险评估提到与 AUC 同一优先级")
    para("附图给出锁定 Status 模型的 O/E 区间。AOCS 为 1.39 [1.03, 1.73]，UKE 为 2.12 [1.51, 2.80]，ICON7 为 1.26 [0.91, 1.59]。正式修稿应同时给出校准截距、斜率和区间，外部校准更新与性能评价分别进行。继续评估时间项时，加入开发与外部评价两阶段重抽样，预先确定有临床依据的编码，并考察诊断时间、记录取整及人群选择。")
    para("检查 IV 期负系数。TCGA 到地标前已观测早死在 III 期及 IV 期分别为 30/415 与 17/78。地标后的系数是条件于仍被观察者的模型参数。建议报告地标进入人数、失访及各阶段事件数，做不同时间窗与分层敏感性分析。八算法、筛选流程和不稳定的亚组结果移至补充。", [
        ("外部验证指南", "https://www.bmj.com/content/384/bmj-2023-074820"),
        ("PROBAST+AI", "https://www.bmj.com/content/388/bmj-2024-082505"),
    ])

    doc.add_page_break()
    h("审稿人三 妇科病理与转化研究", 1)
    para("建议大修。历史队列的浆液性标签、数字分级和现代 HGSC 的病理诊断口径不完全一致。主稿使用来源记录的「晚期浆液性输卵管卵巢癌」。建议附队列组织类型、分级和中央复核映射。ICON7 有 53/224 条来源标记为低等级，UKE 146 人均无可用级别。只保留 ICON7 来源高等级的 169 人及 AOCS 历史 3 级的 85 人分别做事后分析，状态带来的 AUC 增量为 0.219、0.244，再加时间的增量分别为 0.001、−0.004。历史数字分级的结果应按代理指标解释。", [
        ("ICCR 组织学报告", "https://www.iccr-cancer.org/datasets/published-datasets/female-reproductive/ovary-ft-pp/"),
    ])
    h("TCGA 原发肿瘤切片定量已分析")
    para("按患者汇总 GDC 原发肿瘤质控切片的间质细胞百分比，中位每人 2 张。原地标人群 347 人中 342 人匹配到病理字段，含 49 例一年内死亡。间质比例中位数 10%，四分位数 5%–17.5%。同病例、同折数的五折比较，加入间质后 AUC 从 0.7721 到 0.7713，ΔAUC 为 −0.0008（95% CI −0.0048 至 0.0025），Brier 增加约 0.0008。病理调整后，进展状态的 HR 为 4.30（2.90–6.40），间质每增加 10 个百分点的 HR 为 1.07（0.92–1.24）。这是开发队列的探索性分析。", [
        ("GDC 病理字段", "https://gdc.cancer.gov/about-data/gdc-data-processing/biospecimen-data-standardization"),
    ])
    h("APOLLO-2 独立来源病理队列已分析")
    para("公开资料有 70 名经病理复核的 III–IV 期 HGSC 患者，PFS 和 OS 均自诊断起算。残余病灶仅有 R0 和 >1 cm 两类公开标签。纳入地标的 59 人中仅 5 人在随后一年死亡。应用锁定方程，AUC 从基线 0.528 升至加入状态后 0.626，配对差 0.098（−0.079 至 0.311）。再加时间约为 0.626。O/E 为 0.67（0.15–1.29），点估计方向与三个主要外部队列不同。仅原发部位标本的 44 人中仅 4 例死亡。适合补充呈现异质性和样本量局限，无足够事件拟合独立的多变量病理模型。", [
        ("APOLLO-2 原文", "https://doi.org/10.1038/s41698-024-00519-8"),
    ])
    h("其他病理数据库的取舍")
    para("SIGNPOST 公开 80 名 HGSC 患者的 CRS 与定量免疫组化，适合后续做新辅助治疗专门研究。现表个体年龄为区间，72 人有成对 PFS 和 OS 及状态，其中 4 人的 PFS 时长超过 OS，需先核准定义。PTRC-HGSOC 的 158 位患者有 H&E 和铂响应标签，公开表缺本研究所需时间。台湾贝伐珠单抗影像队列存活病例缺末次随访日期。两者适合各自病理任务，无法据当前字段复现主要死亡终点。", [
        ("SIGNPOST", "https://doi.org/10.1038/s41467-024-54295-x"),
        ("PTRC-HGSOC", "https://www.cancerimagingarchive.net/collection/ptrc-hgsoc/"),
        ("台湾队列", "https://doi.org/10.7937/TCIA.985G-EY35"),
    ])

    doc.add_page_break()
    h("修改顺序与真实专家候选", 1)
    h("投稿前优先修改")
    para("第一，改标题、摘要和引言，将 Su 等列为最接近研究。用「进展信息的增量和外部校准」替代「首次发现」。候选英文题目为 Incremental Prognostic Value of Progression Timing Beyond Progression Status in Advanced Serous Tubo-Ovarian Carcinoma。")
    para("第二，将四模型对照和外部校准置于主图，原八算法和变量筛选移至补充。第三，将 TCGA 间质分析及 APOLLO-2 小样本结果列为预审后探索性分析。主结论维持在已观察的地标、时长编码与验证人群范围内。")
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    pic = p.add_run().add_picture(str(HERE / "Review_Evidence.png"), width=Inches(6.64))
    pic._inline.docPr.set("descr", "Added AUC for the clinical and pathology benchmarks and external calibration of the locked Status model.")
    para("图中括号为一年内死亡数。A 图前三行比较完整 Status 与单独进展状态，TCGA 行比较加间质前后折外 AUC。B 图 O/E 小于 1 表示高估，大于 1 表示低估。区间只反映评价病例重抽样。")
    h("真实专家候选")
    para("以下是可考虑邀请的专家，请按投稿期刊规则及实际作者名单核对同机构、合著和共享队列等利益冲突。其身份已按机构公开页面核实。")
    para("Susana Banerjee，Royal Marsden 与 ICR 妇科肿瘤临床试验。适合评估一线治疗、维持治疗背景和随访场景。", [
        ("机构主页", "https://www.icr.ac.uk/research-and-discoveries/find-a-researcher/detail/dr-susana-banerjee"),
    ])
    para("Richard Riley，伯明翰大学生物统计教授。适合评估预后建模、外部验证精度和再校准。", [
        ("机构主页", "https://www.birmingham.ac.uk/staff/profiles/applied-health/riley-richard"),
    ])
    para("Martin Köbel，卡尔加里大学妇科病理学者。适合评估卵巢癌组织类型映射和组织标志物的增量证据。", [
        ("机构主页", "https://profiles.ucalgary.ca/martin-koebel"),
    ])
    doc.save(HERE / "三位同行预审与病理补充报告.docx")


if __name__ == "__main__":
    review()
