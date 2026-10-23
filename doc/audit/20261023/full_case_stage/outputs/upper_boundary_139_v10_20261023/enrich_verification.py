"""Recover seed verification only from its explicit origin artifact, with hashes."""
import hashlib,json,pathlib
AREA=pathlib.Path(__file__).resolve().parent;RUN=AREA.parent.parent
d=json.loads((AREA/'root_inventory.json').read_text());extra=[]
tail_path=RUN/'outputs/tail_beta_roots_20261023.json'
tail=json.loads(tail_path.read_text()) if tail_path.exists() else {}
for station in d['stations'].values():
    for root in station['roots']:
        if root['g_omega_values']:continue
        for source in root['provenance']:
            name=(source.get('seed_provenance') or {}).get('origin_job')
            if not name:continue
            candidates=[]
            for path in [RUN/'outputs/refined_candidates'/(name+'.json'),RUN/'outputs'/(name+'.json')]:
                if path.exists():candidates.append((path,json.loads(path.read_text())))
            candidates.extend((tail_path,v) for v in tail.get('results',[]) if v.get('job',{}).get('name')==name)
            for path,record in candidates:
                result=record.get('result',record);point=result.get('point');verify=result.get('verification')
                if result.get('status')!='candidate_converged' or not point or not verify:continue
                original=source['actual_point']
                if abs(point['omega_r']-original['omega_r'])>1e-6 or abs(point['beta_r']-original['beta_r'])>1e-4:continue
                content=path.read_bytes();snapshot=AREA/'inputs'/('artifact_'+path.name)
                if snapshot.exists():assert snapshot.read_bytes()==content
                else:snapshot.write_bytes(content)
                derivative=verify['J'][0][0]
                if source['derived_joint_conjugate']:derivative=-derivative
                root['g_omega_values'].append(derivative);root['g_beta_beta_values'].append(verify['J'][1][1])
                evidence=dict(root_id=root['root_id'],source=str(path),snapshot=str(snapshot),sha256=hashlib.sha256(content).hexdigest(),
                              origin_job=name,g_omega=derivative,derived_joint_conjugate=source['derived_joint_conjugate'],verification=verify)
                root.setdefault('additional_verification',[]).append(evidence);extra.append(evidence)
                break
            if root['g_omega_values']:break
    highest=station['highest_positive']
    if highest:
        actual=next(v for v in station['roots'] if v['root_id']==highest['root_id'])
        station['highest_positive']=actual
        station['upper_derivative_verified']=bool(actual['g_omega_values']) and all(v<0 for v in actual['g_omega_values'])
        station['g_omega_state']='negative' if station['upper_derivative_verified'] else 'missing_or_not_negative'
d['additional_verification']=extra
(AREA/'root_inventory_enriched.json').write_text(json.dumps(d,indent=2)+'\n')
print(json.dumps(dict(enriched=[dict(root_id=v['root_id'],source=v['source'],g_omega=v['g_omega']) for v in extra],
    highest={k:dict(upper_verified=v['upper_derivative_verified'],g_omega=v['highest_positive']['g_omega_values'] if v['highest_positive'] else []) for k,v in d['stations'].items()})),flush=True)
