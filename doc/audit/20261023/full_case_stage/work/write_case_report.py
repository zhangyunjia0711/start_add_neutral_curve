"""Refresh traceable stage plots and report. Unknowns remain unknown."""
import csv,datetime,json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from full_driver import ROOT,sha
from export_neutral_results import dumpcsv,FIELDS
O=ROOT/"outputs"

def main():
    source=O/"neutral_variable_beta_candidates.json";d=json.loads(source.read_text())
    grid=json.loads((O/'seed_catalog.json').read_text())['grids']['main_841']['x_nondim']
    rows=d["unique_candidates_data"];by={}
    for r in rows:
        if r["omega_r"]>=0 and r["frequency_boundary"]=="upper_local":
            if r["station"] not in by or r["omega_r"]>by[r["station"]]["omega_r"]:by[r["station"]]=r
    upper=[by[i] for i in sorted(by)]
    families={name:sorted([r for r in rows if r['source']==name and r['omega_r']>=0],key=lambda r:r['station'])
        for name in ['tail841_lobe_lower','tail841_low_lobe_upper']}
    band_rows=[dict(r,family=name) for name,points in families.items() for r in points]
    dumpcsv(O/'neutral_variable_beta_additional_bands.csv',band_rows,FIELDS+['family'])
    lookup={name:{r['station']:r for r in points} for name,points in families.items()}
    comparisons=[]
    for i in range(1,842):
        lo=lookup['tail841_lobe_lower'].get(i);up=lookup['tail841_low_lobe_upper'].get(i)
        gap=lo['omega_r']-up['omega_r'] if lo and up else None
        comparisons.append(dict(station=i,x_mm=35*grid[i-1],
            lower_of_high_lobe=lo['omega_r'] if lo else None,
            upper_of_low_lobe=up['omega_r'] if up else None,
            difference=gap,status='candidate_gap_between_two_followed_peaks' if gap is not None and gap>0
                else 'overlap_of_two_local_neutral_intervals' if gap is not None else 'unknown_missing_family',
            qualification='No all-beta or full-spectrum stability assertion'))
    dumpcsv(O/'neutral_variable_beta_two_family_comparison.csv',comparisons,list(comparisons[0]))
    mirrored=[dict(r,omega_r=-r["omega_r"],beta_r=-r["beta_r"],alpha_r=-r["alpha_r"],
        g_omega=-r["g_omega"],g_beta=-r["g_beta"],
        frequency_boundary="lower_local" if r["g_omega"]<0 else "upper_local",
        kind="conjugate_derived_not_an_extra_solve",original_omega_r=r["omega_r"],
        original_beta_r=r["beta_r"]) for r in rows if r["omega_r"]<0]
    dumpcsv(O/"neutral_conjugate_derived.csv",mirrored,FIELDS+["kind","original_omega_r","original_beta_r"])
    zero={};nominal=O/"axial_wall_local_v10_near_homotopy.json"
    for r in json.loads(nominal.read_text())["records"]:
        if r["growth"]>1e-7:zero[int(r["iloc"])]=dict(point=r,result_file=str(nominal))
    for name in ["zero_frequency_witness_v10_20261023","zero_frequency_expanded_v10_20261023",
                 "zero_frequency_neutral_witness_v10_20261023"]:
        for p in (O/name/"stations").glob("i*.json"):
            z=json.loads(p.read_text())
            if z["status"]=="positive_growth_witness":zero[z["station"]]=dict(point=z["witness"],result_file=str(p))
    grid=json.loads((O/"seed_catalog.json").read_text())["grids"]["main_841"]["x_nondim"]
    zrows=[dict(station=i,x_mm=35*grid[i-1],status="positive_growth_witness" if i in zero else "unknown",
        omega_r=0,beta_r=zero[i]["point"]["beta_r"] if i in zero else None,
        growth=zero[i]["point"]["growth"] if i in zero else None,
        result_file=zero[i]["result_file"] if i in zero else None) for i in range(1,842)]
    dumpcsv(O/"zero_frequency_witness_merged.csv",zrows,list(zrows[0]))
    audit_path=O/"upper_boundary_139_v10_20261023/compensation_20261023_2148/boundary_audit_after_compensation.json"
    audit=json.loads(audit_path.read_text());amap={int(r["station"]):r for r in audit["results"]}
    passed=[i for i,r in amap.items() if r["after_state"]=="sampled_nonpositive_known_seeds"]
    unknown=[i for i,r in amap.items() if r["after_state"]=="unknown"]
    cut_path=O/"native10_beta_cuts_v10_20261023_2233/cut_progress.json"
    cap=None
    if cut_path.exists():
        cuts=json.loads(cut_path.read_text())
        centers=[r for r in cuts.get('results',[]) if r.get('tag')=='w40' and r.get('kind')=='center' and r.get('complete')]
        if centers:
            cap=centers[0].get('peak')
    fig,axs=plt.subplots(4,1,figsize=(10.5,12),constrained_layout=True)
    positive=[r for r in rows if r["omega_r"]>=0]
    axs[0].scatter([r["x_mm"] for r in positive],[r["omega_r"] for r in positive],s=3,c=".72",alpha=.6,
        label="Other verified local neutral candidates")
    axs[0].plot([r["x_mm"] for r in upper],[r["omega_r"] for r in upper],c="#172b4d",lw=1.4,
        label="Highest verified upper candidate")
    if cap:
        axs[0].scatter([35*cap['x']],[cap['omega_r']],s=37,marker='^',facecolors='none',
            edgecolors='#b32a2a',lw=1,label='Positive-g test: upper root unresolved')
    for stations,color,marker,label in [(passed,"#237443","o","Cap slice: known sampled paths nonpositive"),
                                       (unknown,"#b32a2a","x","Cap slice: unresolved path")]:
        q=[by[i] for i in stations if i in by]
        axs[0].scatter([r["x_mm"] for r in q],[r["omega_r"] for r in q],s=17,c=color,marker=marker,lw=.6,label=label)
    axs[0].set(xlim=(-35,1100),ylim=(0,None),ylabel=r"$\omega$",
        title="Variable real-beta neutral candidates: 75-degree swept blunt plate")
    axs[0].legend(ncol=2,fontsize=7.5,frameon=False)
    nose=[r for r in upper if r["x_mm"]<-30]
    axs[1].plot([r["x_mm"]+35 for r in nose],[r["omega_r"] for r in nose],".-",ms=3,lw=1,c="#172b4d",
        label="Highest verified upper candidate")
    q=[r for r in mirrored if r["frequency_boundary"]=="lower_local" and r["x_mm"]<-30]
    axs[1].scatter([r["x_mm"]+35 for r in q],[r["omega_r"] for r in q],s=13,c="#b36b13",marker="s",
        label="Conjugate-derived lower candidates")
    first=min(r["x_mm"]+35 for r in upper)
    axs[1].axvspan(.01,first,facecolor=".9",hatch="///",edgecolor=".7",lw=0,label="Upper candidate unknown")
    axs[1].set(xscale="log",xlim=(.01,5),ylim=(0,None),ylabel=r"$\omega$",
        xlabel="Global wall x + 35 mm (mm; logarithmic scale)",title="First 5 mm of the nose")
    axs[1].legend(fontsize=7.5,frameon=False)
    axs[2].plot([r['x_mm'] for r in upper],[r['omega_r'] for r in upper],c='#172b4d',lw=1,
        label='Highest known upper candidate')
    for name,color,marker,label in [
        ('tail841_lobe_lower','#b66b14','x','High-frequency lobe: local lower candidates'),
        ('tail841_low_lobe_upper','#6c4d95','o','Low-frequency lobe: local upper candidates')]:
        points=families[name]
        axs[2].scatter([r['x_mm'] for r in points],[r['omega_r'] for r in points],
            s=5,c=color,marker=marker,label=label)
    axs[2].set(xlim=(200,1100),ylim=(0,10),ylabel=r'$\omega$',
        title='Additional frequency boundaries; gaps are not certified stable')
    axs[2].legend(fontsize=7.5,frameon=False)
    axs[3].plot([r["x_mm"] for r in upper],[r["beta_r"] for r in upper],c="#216444",lw=1.1)
    axs[3].set(xlim=(-35,1100),xlabel="Global wall x (mm)",ylabel=r"$\beta_*$",
        title="Beta belonging to the highest verified upper candidate")
    for ax in axs:ax.grid(alpha=.2)
    for ax in (axs[0],axs[3]):ax.axvline(0,c=".6",lw=.7)
    fig.suptitle("Stage result: local parallel / zero curvature. Finite paths do not certify the full spectrum.",
        fontsize=9.3,color="#94352d")
    fig.savefig(O/"neutral_variable_beta_stage.png",dpi=200);fig.savefig(O/"neutral_variable_beta_stage.svg");plt.close(fig)
    switches=[dict(stations=[i-1,i],x_mm=[by[i-1]["x_mm"],by[i]["x_mm"]],
        beta=[by[i-1]["beta_r"],by[i]["beta_r"]]) for i in sorted(by) if i-1 in by
        and abs(by[i]["beta_r"]-by[i-1]["beta_r"])>1.5]
    maximum=max(upper,key=lambda r:r["omega_r"]);stamp=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
    status=dict(created=stamp,overall_status="incomplete_nose_and_domain_qualification",
        candidate_input_sha256=sha(source),audit_path=str(audit_path),audit_sha256=sha(audit_path),
        native_upper_count=len(upper),unknown_native_upper=[i for i in range(1,842) if i not in by],
        upper_x_mm=[upper[0]["x_mm"],upper[-1]["x_mm"]],unique_candidates=len(rows),
        beta_cap_audit_sampled_nonpositive=passed,beta_cap_audit_unknown=unknown,
        zero_frequency_positive_count_32_841=sum(i in zero for i in range(32,842)),
        zero_frequency_unknown_32_841=[i for i in range(32,842) if i not in zero],
        largest_observed_upper_candidate={k:maximum[k] for k in ["station","x_mm","omega_r","beta_r"]},
        candidate_winner_switches=switches)
    status['native10_positive_cap_test']=dict(point=cap,source_path=str(cut_path),
        qualification='Actual numerical spatial growth; no neutral point or causal mode-direction assertion.') if cap else None
    status['additional_frequency_families']={name:dict(native_count=len(points),
        station_range=[points[0]['station'],points[-1]['station']] if points else None,
        x_mm_range=[points[0]['x_mm'],points[-1]['x_mm']] if points else None)
        for name,points in families.items()}
    status['two_family_gap_candidate_count']=sum(r['status']=='candidate_gap_between_two_followed_peaks' for r in comparisons)
    status['two_family_interval_crossings']=[dict(stations=[a['station'],b['station']],x_mm=[a['x_mm'],b['x_mm']],
        differences=[a['difference'],b['difference']]) for a,b in zip(comparisons,comparisons[1:])
        if a['difference'] is not None and b['difference'] is not None and a['difference']*b['difference']<0]
    gap_path=O/'tail841_gap_midfrequency_v10_20261024/gap_audit_report.json'
    gap=json.loads(gap_path.read_text()) if gap_path.exists() else None
    if gap and gap.get('finished'):
        status['tail841_midfrequency_audit']={k:gap[k] for k in ['target_omega','beta_domain',
            'known_seed_count','state_counts','total_output_points','vector_bound_points']}
        status['tail841_midfrequency_audit']['path']=str(gap_path)
        status['tail841_midfrequency_audit']['sha256']=sha(gap_path)
    (O/"case_stage_status.json").write_text(json.dumps(status,indent=2,allow_nan=False)+"\n")
    text=f"""# 本工况可变实β中性边界：阶段结果

更新：{stamp}（Asia/Shanghai）。**整工况边界尚未完成验收。**

已有{len(upper)}/841个原生站的减半差分步长验证上中性候选，实际范围x={upper[0]['x_mm']:.9f}～1100mm。未知原生站：{status['unknown_native_upper']}。共有{len(rows)}个去重局部中性候选；其来源标签不是永久物理模态身份。

139代表站的β∈[−40,40]、最高已知上根外侧ω+.05单频切片审计：{len(passed)}站已知路径采样非正，{len(unknown)}站未知。未知站为{unknown}。这不是连续β域、所有ω或全谱穷尽证明。原122/17首轮报告、一次有界补偿和不同返回根/超时均保留。

ω=0：第32～841站中{status['zero_frequency_positive_count_32_841']}/810站已有实际正增长见证。这只证明G(x,0)>0，不证明0到上根之间不存在稳定孔洞；图不填整片不稳定面积。

## 图表与定义

[阶段图](neutral_variable_beta_stage.png)；[全部实际候选](neutral_variable_beta_candidates.csv)；
[最高已知上根](neutral_variable_beta_candidates_highest_known_upper.csv)；
[841站覆盖](neutral_variable_beta_candidates_coverage_841.csv)；
[零频正增长见证](zero_frequency_witness_merged.csv)；
[共同共轭派生下根](neutral_conjugate_derived.csv)。

新增实际频带边界见[附加上下中性候选](neutral_variable_beta_additional_bands.csv)、[两条谱支的频率间隙比较](neutral_variable_beta_two_family_comparison.csv)。第841站已验证低频上根ω≈2.095632566476、β≈2.868800748341，以及高频下根ω≈3.476081395914、β≈3.624489572181。它们满足负β曲率，gω分别为负、正；均以减半步长与原生PEP回代核对。两条已追峰之间的负增长段只是候选间隙，其他β或谱支可能填补，图不涂成全谱稳定区。

附加高频下边界已有{len(families['tail841_lobe_lower'])}个原生站（i508～841，x≈204.035676～1100mm），低频上边界已有{len(families['tail841_low_lobe_upper'])}个原生站（i561～841，x≈285.285445～1100mm）；各自区间无缺站。{len(band_rows)}个点均逐一重新读取h/2实际9点field，核对α、g、导数、曲率、原始运行摘要、冻结v10与原输入SHA；见additional_boundaries_field_audit_20261024.json。两条局部区间交叠/分离的转变括在i620/621，即x≈390.310559～392.270575mm。更上游相邻i507/i560的有限局部峰校正未通过；正曲率/无下降失败均保留，不能将其当作谱支不存在的证明。

实ω、实β、复α，g=−Imα。目标有限域G_B=max g；实际候选满足g=0、gβ=0、gββ<0，上根还要求gω<0。[Nomura 2006 §2.3](https://www.jstage.jst.go.jp/article/tjsass/48/162/48_162_193/_pdf/-char/en)，[Mack 1984 §6.6](https://transition.arizona.edu/sites/default/files/2021-08/mack-bllst.pdf)。

ω和β均无量纲，未换算Hz；x是全局壁面坐标mm。负频实际根保持符号；共同共轭(ω,β,α)→(−ω,−β,−conjα)的派生点单独标记，gω符号反转，不能直接取|ω|当上根。最高候选的β切换相邻站：{[r['stations'] for r in switches]}，精确x/β见case_stage_status.json。最大已观察上候选：i{maximum['station']}，x={maximum['x_mm']:.9f}mm，ω={maximum['omega_r']:.12g}、β={maximum['beta_r']:.12g}；不是连续全谱极值证明。

## 模型、验证与未解决项

实际模型为冻结壁面切/法向、重算法向导数、零显式曲率、局部平行近似。没有因此补齐完整附着线/非平行模型；鼻尖切向速度接近零、空间根趋近和极窄峰导致最前端尚未验收。[Xi, Ren & Fu 2021](https://arxiv.org/abs/2006.05970)说明曲率及外无粘流可能重要，不能仅因模态远离壁面就删除。

名义ω0、β−7谱支已在841站取得实际原生端点；辅助剖面仅作数值延拓并在原H5验证，不计物理站。实际400/799/1597法向检查来自同一CFD场插值，不是CFD网格独立；部分中性点只固定原ωβ比较并线性估频移，未在细网格重新中性/β峰收敛。外层52.66根的增长符号随法向离散改变、向量身份敏感，单独保留未分类。

MAC≥.95是工程连续性门槛，小残差不自动证明前向误差或物理身份。尾端同参数的两实际根可分别正/负增长，负采样不能判全谱稳定。[SLEPc可靠性说明](https://slepc.upv.es/release/documentation/manual/pep.html#reliability-of-the-computed-solution)。逐x优化β也不是一条固定β扰动的N因子。[Mack 1984 §14.4](https://transition.arizona.edu/sites/default/files/2021-08/mack-bllst.pdf)。

## 程序与复现

[用户仓库](https://github.com/zhangyunjia0711/start_add_neutral_curve)。
源码已上传codex/full-beta-20261022提交fbf78b13935e8a6dabb467ebc49f35e81bf0577b；原版backup/original-before-neutral-20261021保留。
原程序目录/home/zaza/start/START_ALL_IN_ONE/START_ALL_IN_ONE。运行新版bin/START_LST_wall_local_v10，SHA634b45c5d30e076f5b86ccb136e0734c1905beedde108477dceb8a75adc3853c；原bin/START_LST基线仍保留。8项回归及安装smoke证据归档。

本轮完整记录位于{ROOT}，runs/保存配置/solver.log/summary，work/保存复现脚本和必要派生数据，outputs/保存结果/日志，commands.log保存实际命令。原始H5未改。无需重算可执行：

~~~bash
cd "{ROOT}"
python3 work/export_neutral_results.py
python3 work/write_case_report.py
~~~

矩阵换元已在原生i10/i16取得并用原空间PEP验证中性驻点，但实际β曲率为正，均是局部最小值，禁止并入最大增长包络；详细代数、原始证据与分类见[矩阵换元诊断](矩阵换元验证与鼻尖驻点分类.md)。

剩余：鼻尖原生中性根、有限域未知及未排除的频率孔洞/其他谱支，均不能用连线代替。

复α换元实际求得第10/11站局部鞍点，Imω分别约−1.02548e−5、−4.63671e−5，原空间PEP邻根呈局部平方根收拢。未建立固定β全部实α相关时间奇点的严格高Imω上界，也未把邻支追溯到相反α半平面的因果线，不能据此将鼻部根定为物理下游/绝对/对流不稳定。原文pinch条件见[Huerre与Monkewitz 1990，pp480–484](https://www.annualreviews.org/content/journals/10.1146/annurev.fl.22.010190.002353)，实际验证与局限见[空间根方向核验](spatial_root_causality_saddle_review_20261024.md)。
"""
    (O/"本工况可变β计算说明.md").write_text(text)
    if cap:
        with (O/"本工况可变β计算说明.md").open('a') as f:
            f.write('\n\n## 第10站未闭合的实际证据\n\n'
                'x=−34.939296356mm，ω=40、β≈12.137901的连续源空间根g≈+3.101228；'
                '同ω在β=+40和−40的实际g约+83.602745/+7.158498。'
                'Ωmin+.05切片的两个β端也为正。该测试框不能作为稳定闭合边界；'
                '图中红色三角是正增长试验点，不是中性根。'
                '数值空间根的因果/方向分类仍在核验，不能仅按大正g或Reωα符号宣称物理下游不稳定。\n')
    if gap and gap.get('finished'):
        with (O/'本工况可变β计算说明.md').open('a') as f:
            f.write('\n\n## 尾端中间频率有限多路径检查\n\n'
                f"i841，ω={gap['target_omega']:.12g}，β∈[−40,40]，从6个独立实际来源检查。"
                f"5条完整已知路径未采到正增长，1条M2上峰源到中频的输运120s超时，保留unknown。"
                f"{gap['vector_bound_points']}/{gap['total_output_points']}个输出点已绑定原生Vec/α/残差。"
                '此结果仍不能证明整段频率、连续β域或全谱稳定。'
                '[审计报告](tail841_gap_midfrequency_v10_20261024/gap_audit_report.json)；'
                '[逐源摘要](tail841_gap_midfrequency_v10_20261024/gap_source_summary.csv)。\n')
    print(json.dumps({k:status[k] for k in ["created","native_upper_count","unknown_native_upper",
        "zero_frequency_positive_count_32_841"]}),flush=True)

if __name__=="__main__":main()

