"""Read-only package/header audit. Extraction only into this task's directory."""
import struct,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
def main():
    path=ROOT/'baseflow_different_sweptangale.lpk';b=path.read_bytes();assert b[:12]==b'#!PK8000\r\n\x1a\n'
    at=24
    def entry():
        nonlocal at
        n=struct.unpack_from('<i',b,at)[0];at+=4
        assert 0<n<1000
        name=b[at:at+n].rstrip(b'\0').decode('ascii');at+=n
        ext=b[at:at+4].rstrip(b'\0').decode('ascii');at+=4
        offset,size=struct.unpack_from('<II',b,at);at+=8
        assert 0<offset<=len(b) and offset+size<=len(b),(name,offset,size)
        return dict(name=name,ext=ext,offset=offset,size=size)
    layout=entry();count=struct.unpack_from('<i',b,at)[0];at+=4
    entries=[entry() for _ in range(count)];dest=OUT/'package_inspection';dest.mkdir(exist_ok=True)
    for r in [layout]+entries:
        data=b[r['offset']:r['offset']+r['size']];r['sha256']=hashlib.sha256(data).hexdigest()
        if r['ext']=='lay':
            assert data.startswith(b'#!MC');(dest/r['name']).write_bytes(data)
        else:
            assert data.startswith(b'#!TDV');pos=16 if data[5:8] in [b'112',b'191'] else 12
            def string():
                nonlocal pos
                chars=[]
                while True:
                    k=struct.unpack_from('<i',data,pos)[0];pos+=4
                    if k==0:return ''.join(chars)
                    assert 0<=k<=65535;chars.append(chr(k))
            r['version']=data[:8].decode('ascii');r['title']=string();nv=struct.unpack_from('<i',data,pos)[0];pos+=4
            assert 0<nv<200;r['variables']=[string() for _ in range(nv)]
            print(r['name'],r['title'],r['variables'],flush=True)
    result=dict(source=str(path),source_sha256=hashlib.sha256(b).hexdigest(),layout=layout,datasets=entries)
    (OUT/'package_inventory.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
