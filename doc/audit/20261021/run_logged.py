#!/usr/bin/env python3
"""Run one bash command and keep dated command/output records (Asia/Shanghai)."""
import datetime, json, pathlib, subprocess, sys
root=pathlib.Path(__file__).resolve().parent
stamp=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
command=sys.argv[1]
log=root/'commands.log'
with log.open('a',encoding='utf-8') as out:
    out.write('\n['+stamp+'] COMMAND '+json.dumps(command,ensure_ascii=False)+'\n');out.flush()
    proc=subprocess.Popen(['bash','-lc',command],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,errors='replace')
    for line in proc.stdout:
        sys.stdout.write(line);sys.stdout.flush();out.write(line);out.flush()
    code=proc.wait();out.write('[EXIT '+str(code)+']\n')
sys.exit(code)
