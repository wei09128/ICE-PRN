"""Check the frozen source snapshot and packaged inputs without rerunning inference."""
from pathlib import Path, PurePosixPath
import hashlib, json, zipfile

ROOT=Path(__file__).resolve().parent

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    frozen=json.loads((ROOT/'RELEASE_HASHES.json').read_text())
    for name,want in frozen.items():
        p=ROOT/name
        if not p.is_file() or digest(p)!=want:raise SystemExit('FAIL source: '+name)
    inputs=json.loads((ROOT/'figure_runner/INPUT_HASHES.json').read_text())
    for name,want in inputs.items():
        if digest(ROOT/'figure_runner'/name)!=want:raise SystemExit('FAIL input: '+name)
    count=0
    for p in ROOT.rglob('*'):
        if '.git' in p.parts or not p.is_file():continue
        if p.stat().st_size>=100*1024**2:raise SystemExit('File too large: '+str(p))
        if p.suffix.lower()=='.zip':
            with zipfile.ZipFile(p) as z:
                for item in z.infolist():
                    n=PurePosixPath(item.filename)
                    if n.is_absolute() or '..' in n.parts or '\\' in item.filename:
                        raise SystemExit('Unsafe ZIP path: '+str(p))
                bad=z.testzip()
                if bad:raise SystemExit('Corrupt ZIP member: '+bad)
            count+=1
    print(f'PASS: {len(frozen)} source hashes, {len(inputs)} figure inputs, {count} ZIP archives; no file ≥100 MiB.')

if __name__=='__main__':main()
