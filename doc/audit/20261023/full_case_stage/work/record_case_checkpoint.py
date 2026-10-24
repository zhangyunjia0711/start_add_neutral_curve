"""Version a small, explicit stage evidence snapshot; no CFD/run-tree copy."""
import datetime,json,shutil,subprocess
from pathlib import Path
from full_driver import ROOT,sha,snapshot
SOURCE=ROOT/"work/solver_source"
TARGET=SOURCE/"doc/audit/20261023/full_case_stage"
OUTPUTS=["neutral_variable_beta_stage.png","neutral_variable_beta_stage.svg",
    "neutral_variable_beta_candidates.csv","neutral_variable_beta_candidates_highest_known_upper.csv",
    "neutral_variable_beta_candidates_coverage_841.csv","neutral_variable_beta_candidates_coverage_139.csv",
    "neutral_conjugate_derived.csv","zero_frequency_witness_merged.csv",
    "neutral_variable_beta_additional_bands.csv","neutral_variable_beta_two_family_comparison.csv",
    "case_stage_status.json","本工况可变β计算说明.md",
    "矩阵换元验证与鼻尖驻点分类.md","owned_legacy_attempt_cleanup_20261024.json",
    "additional_boundaries_field_audit_20261024.json",
    "tail841_all_boundaries_20261024.json","tail841_finite_ridge_audit_20261024.json",
    "tail841_finite_peak_paths_20261024.csv","tail841_finite_peak_paths_20261024.png",
    "spatial_root_causality_saddle_review_20261024.md","spatial_root_causality_saddle_review_20261024.json",
    "complex_alpha_saddle_native10_20261024.json","complex_alpha_saddle_native11_20261024.json",
    "complex_alpha_saddle_native10_neighbors_20261024.json",
    "工作日志_2026-10-24_鼻部谱支核验与完整边界收尾.md",
    "工作日志_2026-10-24_native10切片归约与复频率接口.md",
    "工作日志_2026-10-24_空间根方向与复alpha鞍点核验.md",
    "工作日志_2026-10-24_尾端双峰低频连通性核验.md",
    "工作日志_2026-10-24_尾端中频孔洞检查.md",
    "original_checkout_integration.json","installed_v10_smoke.json",
    "flat_neutral_normal_validation_20261023.json","local_mode_audit_M3near_neutral_normal_validation.json",
    "local_mode_audit_normal_validation.json","outer_grid_audit_20261023.json",
    "primary_crossover_native_i372.json","tail_crossover_native_i773.json",
    "工作日志_2026-10-22_完整工况可变β中性边界.md",
    "工作日志_2026-10-23_完整工况中性边界续算.md"]
EXTRA=["outputs/upper_boundary_139_v10_20261023/boundary_audit_report.json",
    "outputs/upper_boundary_139_v10_20261023/compensation_20261023_2148/compensation_report.json",
    "outputs/upper_boundary_139_v10_20261023/compensation_20261023_2148/boundary_audit_after_compensation.json",
    "outputs/pep_matrix_temporal_probe_20261023.json",
    "outputs/narrow_ridge_second_pair_recovery_20261023.json"]
for folder,names in {
    "native10_beta_cuts_v10_20261023_2233":["native10_cut_report.json","native10_cut_samples.csv",
        "native10_cut_summary.csv","native10_cut_witnesses.csv","native10_beta_cuts.png",
        "reduction_manifest.json","final_source_integrity.json","README.md"],
    "complex_omega_interface_v10_20261024":["complex_omega_probe.json","complex_omega_points.csv",
        "manifest.json","final_source_integrity.json","README.md"],
    "dense_tail841_lower_20261024":["progress.json"],
    "dense_tail841_low_upper_20261024":["progress.json"],
    "tail841_gap_midfrequency_v10_20261024":["gap_audit_report.json","gap_source_summary.csv",
        "gap_samples.csv","gap_positive_witnesses.csv","gap_beta_cuts.png","README.md","manifest.json"],
}.items():
    EXTRA.extend('outputs/'+folder+'/'+name for name in names)
def main():
    snapshot()
    assert SOURCE.resolve() in TARGET.resolve().parents
    TARGET.mkdir(parents=True,exist_ok=True)
    selected=[ROOT/"outputs"/n for n in OUTPUTS]+[ROOT/p for p in EXTRA]
    selected+=list((ROOT/"work").glob("*.py"))
    for folder in ["local_mode_normal_validation_20261022","flat_neutral_normal_validation_20261023",
        "pep_matrix_temporal_probe_20261023","complex_alpha_saddle_20261024",
        "tail841_boundary_native_20261024"]:
        selected+=list((ROOT/"work"/folder).glob("*.py"))
    for folder in ["upper_boundary_139_v10_20261023","upper_boundary_139_v10_20261023/compensation_20261023_2148",
        "native10_beta_cuts_v10_20261023_2233","complex_omega_interface_v10_20261024",
        "tail841_gap_midfrequency_v10_20261024"]:
        selected+=list((ROOT/"outputs"/folder).glob("*.py"))
    copied=[]
    for source in sorted(set(selected)):
        if not source.exists():continue
        rel=source.relative_to(ROOT);dest=TARGET/rel
        assert TARGET.resolve() in dest.resolve().parents
        dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(source,dest)
        copied.append(dict(source=str(source),archive=str(dest.relative_to(SOURCE)),
            bytes=dest.stat().st_size,sha256=sha(dest)))
    manifest=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        stage="incomplete_nose_and_domain_qualification",calculation_root=str(ROOT),
        core_commit="fbf78b13935e8a6dabb467ebc49f35e81bf0577b",
        files=copied,original_inputs=snapshot(),
        qualification="Small explicit evidence/code snapshot; raw CFD and full run trees stay at calculation_root.")
    (TARGET/"archive_manifest.json").write_text(json.dumps(manifest,indent=2,allow_nan=False)+"\n")
    (TARGET/"README.md").write_text("""# 75°工况阶段快照

这是未完成鼻尖/有限域验收时的版本检查点，不是完整认证结果。查看 outputs/本工况可变β计算说明.md 与 case_stage_status.json。

work/保存实际使用的脚本快照，依赖原工况目录、配置、H5及先前归档结果；不能把本目录当独立可运行示例。原完整计算目录和每文件SHA见archive_manifest.json。需要复现时保持原ROOT/work与ROOT/outputs结构和原CASE/M1/M2/M3输入，使用冻结v10（核心fbf78b1）；每个真实run的config.ini/solver.log/summary在原ROOT/runs。未将大型原始CFD或全部求解日志复制入Git。

不同数值/来源模式明确区分：实际中性点、正增长见证、共同共轭派生点、仅辅助剖面点、失败/unknown。此提交只归档证据与脚本，不修改Fortran核心。
""")
    print(json.dumps(dict(files=len(copied),bytes=sum(r["bytes"] for r in copied),
        target=str(TARGET))),flush=True)
if __name__=="__main__":main()

