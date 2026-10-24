"""Copy explicitly archived output artifacts to a task-owned delivery folder."""
import argparse,datetime,json,shutil
from pathlib import Path
from full_driver import ROOT,sha


def main():
    p=argparse.ArgumentParser();p.add_argument('--destination',type=Path,required=True);a=p.parse_args()
    dest=a.destination.resolve();manifest=ROOT/'work/solver_source/doc/audit/20261023/full_case_stage/archive_manifest.json'
    assert dest.name=='variable_beta_full_20261024'
    owner=dest/'delivery_manifest.json'
    if dest.exists():
        assert owner.exists() and json.loads(owner.read_text())['calculation_root']==str(ROOT)
    else:dest.mkdir(parents=True)
    chosen=[Path(r['source']) for r in json.loads(manifest.read_text())['files']
            if (ROOT/'outputs') in Path(r['source']).parents]
    chosen.append(ROOT/'outputs/neutral_variable_beta_candidates.json')
    files=[]
    for source in sorted(set(chosen)):
        rel=source.relative_to(ROOT/'outputs');target=dest/rel
        assert dest in target.resolve().parents
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        assert sha(source)==sha(target)
        files.append(dict(path=str(rel),bytes=target.stat().st_size,sha256=sha(target),source=str(source)))
    text='''本目录为2026-10-24阶段交付。

请先看本工况可变β计算说明.md与neutral_variable_beta_stage.png。
所有omega/beta无量纲，x为全局壁面坐标mm；模型为局部壁面平行/零显式曲率。
这是已验局部谱支的阶段候选集：鼻部最前10站上根及因果分类未解决，有限路径不能认证完整全谱包络。
新增高频下边界334站、低频上边界281站；候选频率间隙不作全谱稳定宣称。
完整原始CFD/运行日志/向量/复现脚本仍在delivery_manifest.json的calculation_root。
源码与小型证据快照在用户GitHub仓库codex/full-beta-20261022分支。
'''
    (dest/'请先阅读.txt').write_text(text,encoding='utf-8')
    files.append(dict(path='请先阅读.txt',bytes=(dest/'请先阅读.txt').stat().st_size,sha256=sha(dest/'请先阅读.txt')))
    result=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        calculation_root=str(ROOT),archive_manifest_sha256=sha(manifest),files=files,
        status='incomplete_nose_and_domain_qualification',no_original_inputs_modified=True)
    owner.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(destination=str(dest),files=len(files),bytes=sum(v['bytes'] for v in files))))


if __name__=='__main__':main()
