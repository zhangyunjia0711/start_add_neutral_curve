import pathlib,json,collections,datetime
A=pathlib.Path(__file__).resolve().parent
q=json.loads((A/'boundary_progress.json').read_text())
r=json.loads((A/'full8_resources.json').read_text())
print(json.dumps(dict(now=datetime.datetime.now().astimezone().isoformat(),states=dict(collections.Counter(v['state'] for v in q['results'])),
 unknown_stations=[v['station'] for v in q['results'] if v['state']=='unknown'],running_stations=[v['station'] for v in q['results'] if v['state']=='running'],
 rss_peak_MiB=r['max_process_rss_bytes']/2**20,minimum_available_GiB=r['minimum_available_bytes']/2**30,
 current_available_GiB=r['samples'][-1]['mem_available_bytes']/2**30,load=r['samples'][-1]['loadavg'],
 positive_events=json.loads((A/'positive_events.json').read_text())['events'] if (A/'positive_events.json').exists() else []),ensure_ascii=False))

