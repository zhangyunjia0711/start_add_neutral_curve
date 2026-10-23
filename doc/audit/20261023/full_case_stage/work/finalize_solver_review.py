"""Finish the v4 review without altering historical output records."""
import datetime,hashlib,json,pathlib,re,subprocess
from full_driver import ROOT
freeze=json.loads((ROOT/'outputs/solver_v4_freeze.json').read_text())
vector=json.loads((ROOT/'outputs/vector_validation_v4.json').read_text())
raw=json.loads((ROOT/'outputs/raw_vector_regression_attempt3.json').read_text())
assert vector['sha256']==freeze['sha256']==raw['solver_sha256']
assert hashlib.sha256(pathlib.Path(freeze['solver']).read_bytes()).hexdigest()==freeze['sha256']
assert raw['passed'] and raw['petsc_errors']==raw['unused_options']==0
log=(ROOT/'outputs/raw_vector_regression_attempt3.log').read_text()
line=next(s for s in log.splitlines() if 'NORMALIZATION MEMORY BEFORE/AFTER/DELTA' in s)
values=[float(v.replace('D','E')) for v in line.split(':',1)[1].split()]
assert len(values)==3 and values[2]==0
original=pathlib.Path('/home/zaza/start/START_ALL_IN_ONE/START_ALL_IN_ONE/bin/START_LST')
original_sha=hashlib.sha256(original.read_bytes()).hexdigest();assert original_sha=='807169bb401e558d9509c11f75304e5694d1d8af6925660fca251e8bc29f469a'
stamp=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
body='''

## v4完成补录（{stamp}）

本轮已统一重新构建并通过两项现有CTest。主代理补丁已核查并编译：使用raw API的selected_eigen/error核验，初始向量在saved tracking下捕获，缺origin参考明确失败，所有新增Vec操作逐项检查ierr并在不可恢复错误时终止；连续尝试上限4000。已释放每次raw临时current，history独立拷贝和正常出口释放逻辑保留。

新增raw回归前两次记录为失败：实际创建前将输出赋为PETSC_NULL_VEC或0会被本地Fortran wrapper视作NULL输出地址。最终改为普通局部created_vector接收VecDuplicate，再赋给返回对象；不再将无效返回sentinel直接传作PETSc输出。`raw_vector_regression.json`及`raw_vector_regression_attempt2.json`保留为失败尝试；最终接受记录是`raw_vector_regression_attempt3.json/.log`。另一次测试编译的VecScale实参由real改为PetscScalar complex，属于测试修正。

最终raw回归通过PEP/EPS选中pair与残差绑定、tracking开关不重选、独立复制及跨下一Solve存活检查。5000次normal_u前后PETSc当前分配量为{before:.0f}→{after:.0f}字节，增量{delta:.0f}。这只证明该重复归一化段没有累计分配增长。退出时-malloc_debug仍有其他分配栈；核查现有test_slepc_adapter源码确实调用了两个Finalize，该既有测试独立运行也出现残留记录（`adapter_finalize_diagnostic.log`）。不能据此宣称整个进程无任何泄漏，残留对象的完整归因未完成；没有为此修改PETSc。

v4冻结二进制：`{solver}`，SHA256 `{sha}`。
冻结清单`solver_v4_freeze.json`记录源码hash、分支、未commit状态及原输入/程序核验；最终源码补丁是`source_v4_snapshot.patch`。较早`raw_vector_resource_fix.patch`仅保留最初尝试，包含后来修复的NULL输出问题，不能用它替代最终补丁。

本代理实际使用冻结v4复算nose i119 omega0及neutral i321 omega20，在0.08/0.4延拓步长下各比较25个共同beta样本。最大复alpha差分别3.244606e-11和5.378317e-12；4次均rejected0、PETScerror0、unused0。nose两个峰growth差约6.4e-14，峰beta差3.07e-6；二者beta不完全相同，不能把峰alpha直接当同参数点比较。正式记录为`vector_validation_v4.json`，旧v3记录未覆盖。

原二进制当前hash仍为`{oldsha}`。主代理及scope代理在v4冻结后开始v5壁面局部坐标/sample_only开发；本次v4检验不能代替新增模型的验证。此处完成的是v4求解器支持工作，完整工况计算仍由主代理继续。
'''.format(stamp=stamp,before=values[0],after=values[1],delta=values[2],solver=freeze['solver'],sha=freeze['sha256'],oldsha=original_sha)
with (ROOT/'outputs/solver_build_review.md').open('a') as f:f.write(body)
result=dict(created=stamp,solver=freeze['solver'],sha256=freeze['sha256'],original_sha256=original_sha,normalization_memory_before=values[0],normalization_memory_after=values[1],normalization_memory_delta=values[2],raw_regression='outputs/raw_vector_regression_attempt3.json',vector_validation='outputs/vector_validation_v4.json',review='outputs/solver_build_review.md',remaining='Other malloc_debug finalization allocations not fully attributed; full-case and v5 model work continues in root task')
out=ROOT/'outputs/solver_review_completion.json';assert not out.exists();out.write_text(json.dumps(result,indent=2)+'\n');print('SOLVER REVIEW COMPLETE '+json.dumps(result),flush=True)