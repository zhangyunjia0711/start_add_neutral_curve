"""Publish verified fixed-beta x-omega contours of both geometries."""
import csv
import datetime
import hashlib
import json
import re
import shutil
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

WIN = Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°')
ROOT = Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_xomega_multiple_beta_20261027')
sys.path.insert(0, str(ROOT / 'work'))
from solver import Pencil
from common import snapshot

DEST = WIN / 'outputs/M3_xomega_multiple_beta_20261027'
BETAS = [-30., -20., -10., -7., 0., 10., 20.]
CASES = ['bluntplate', 'rough']
LABELS = ['Baseline blunt plate', 'Optimized nose: continuous curvature']
COLORS = {-30.: '#762a83', -20.: '#2166ac', -10.: '#159daf', -7.: '#161616',
          0.: '#1b9e77', 10.: '#e48120', 20.: '#d73027'}
W = 6.0727867021443425 / 50
NAME = 'M3_xomega_multiple_beta'
NAMES = [NAME] + [case + '_' + NAME for case in CASES]
VARS = ['x_mm', 'omega_nd', 'f_kHz', 'beta', 'alpha_r', 'alpha_i', 'growth',
        'station_coordinate', 'point_kind', 'case_id', 'polynomial_residual']
SOURCES = [
    'bluntplate_beta-30_c0_final.json', 'bluntplate_beta-20_c0_final.json', 'bluntplate_beta-20_c1_final.json',
    'bluntplate_beta-10_c0_final.json', 'bluntplate_beta+0_c0_final.json', 'bluntplate_beta+10_c0_final.json',
    'bluntplate_beta+20_temporal_c0_final.json', 'rough_beta-30_dense_final.json',
    'rough_beta-20_physical_c0_final.json', 'rough_beta-20_c2_retry3_final.json',
    'rough_beta-10_temporal_c0_final.json', 'rough_beta+0_c0_retry1_final.json', 'rough_beta+10_temporal_c0_final.json']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def js(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def parse(path):
    variables, zones, inside = [], [], False
    for line in path.read_text(encoding='ascii').splitlines():
        s = line.strip()
        if s.startswith('VARIABLES'):
            inside = True
            variables += re.findall(r'"([^"]+)"', s)
            continue
        if inside and s.startswith('"'):
            variables += re.findall(r'"([^"]+)"', s)
            continue
        inside = False
        if s.startswith('ZONE '):
            zones.append(dict(name=re.search(r'T\s*=\s*"([^"]+)"', s).group(1), rows=[], double=False))
        elif s.startswith('DT=') or s.startswith('DT ='):
            assert len(re.findall('DOUBLE', s)) == len(variables)
            zones[-1]['double'] = True
        elif zones and re.match(r'^[+-]?(?:[0-9]|\.[0-9])', s):
            row = [float(v) for v in s.split()]
            assert len(row) == len(variables) and np.isfinite(row).all()
            zones[-1]['rows'].append(row)
    assert all(z['double'] for z in zones)
    return variables, zones


def load():
    curves, audit = [], []
    for name in SOURCES:
        path = ROOT / 'outputs' / name
        d = json.loads(path.read_text())
        assert d['status'] == 'closed_verified'
        case, beta = d['case'], d['beta']
        points = d['points']
        assert points[-1]['kind'] == 2
        assert all(points[0][k] == points[-1][k] for k in ['xi', 'x_mm', 'alpha_r', 'omega_r', 'beta_r'])
        assert d['minimum_adjacent_raw_MAC'] > .9
        a = np.array([[p['x_mm'], p['omega_r'], p['omega_r'] / W, beta, p['alpha_r'], p['alpha_i'],
                       p['growth'], p['xi'], p['kind'], CASES.index(case) + 1, p['right_residual']] for p in points])
        checks = []
        for p in points:
            if p['kind'] != 0:
                continue
            pencil = Pencil(case, int(round(p['xi'])), p['omega_r'])
            alpha = complex(p['alpha_r'], p['alpha_i'])
            K = pencil.spatial(beta)
            q = np.load(p['vector_file'])
            residual = float(np.linalg.norm((K[0] + alpha * K[1] + alpha * alpha * K[2]) @ q))
            assert residual < 1e-7
            checks.append(residual)
        curves.append(dict(case=case, beta=beta, rows=a, source=str(path), source_sha256=sha(path), legacy=False))
        audit.append(dict(source=str(path), sha256=sha(path), case=case, beta=beta, points=len(points),
                          native_points_rechecked=len(checks), maximum_native_residual=max(checks),
                          minimum_adjacent_MAC=d['minimum_adjacent_raw_MAC'], physical_tips=len(d['physical_tips']),
                          independent_START_checks=len(d['original_checks']) + len(d.get('dense_original_checks', []))))
    previous = json.loads((WIN / 'outputs/comparison_M3_beta_minus7_20261025/comparison_summary.json').read_text())
    for case in CASES:
        source = next(r for r in previous['sources'] if r['path'].endswith('_curve.csv') and
                      (('/rough/' in r['path']) if case == 'rough' else ('/bluntplate/' in r['path'])))
        path = Path(source['path'])
        assert sha(path) == source['sha256']
        with path.open() as stream:
            records = list(csv.DictReader(stream))
        a = np.array([[float(r['x_mm']), float(r['omega_r']), float(r['omega_r']) / W, -7., float(r['alpha_r']),
                       float(r['alpha_i']), -float(r['alpha_i']), float(r['xi']), float(r['point_kind']),
                       CASES.index(case) + 1, float(r['polynomial_residual'])] for r in records])
        native = [r for r in records if float(r['point_kind']) == 0]
        checks = []
        for j in sorted(set(np.linspace(0, len(native) - 1, 14, dtype=int).tolist())):
            row = native[j]
            metadata = json.loads(Path(row['source_file']).read_text())
            q = np.load(metadata['vector_file'])
            pencil = Pencil(case, int(float(row['xi'])), float(row['omega_r']))
            alpha = complex(float(row['alpha_r']), float(row['alpha_i']))
            K = pencil.spatial(-7.)
            residual = float(np.linalg.norm((K[0] + alpha * K[1] + alpha * alpha * K[2]) @ q))
            assert residual < 1e-7
            checks.append(dict(station=int(float(row['xi'])), residual=residual,
                               source_metadata=row['source_file'], vector_file=metadata['vector_file']))
        curves.append(dict(case=case, beta=-7., rows=a, source=str(path), source_sha256=sha(path), legacy=True))
        audit.append(dict(source=str(path), sha256=sha(path), case=case, beta=-7., points=len(a), legacy=True,
                          current_joint_matrix_checks=checks))
    for c in curves:
        a = c['rows']
        assert np.isfinite(a).all() and np.max(abs(a[:, 5])) < 1e-8 and np.max(a[:, 10]) < 1e-7
        assert np.array_equal(a[0, :4], a[-1, :4]) and a[-1, 8] == 2
    curves.sort(key=lambda c: (CASES.index(c['case']), c['beta'], float(c['rows'][:, 1].mean())))
    assert len(curves) == 15 and {c['beta'] for c in curves if c['case'] == 'bluntplate'} == set(BETAS)
    assert {c['beta'] for c in curves if c['case'] == 'rough'} == set(BETAS) - {20.}
    seeds = json.loads((ROOT / 'outputs/seeds.json').read_text())['results']
    matched = []
    for j, s in enumerate(seeds):
        if s['status'] != 'neutral_corrected':
            continue
        r = s['point']
        hits = []
        for k, c in enumerate(curves):
            if c['case'] != s['case'] or c['beta'] != s['beta']:
                continue
            a = c['rows']
            if np.any((abs(a[:, 7] - r['xi']) < 1e-7) & (abs(a[:, 4] - r['alpha_r']) < 1e-5) &
                      (abs(a[:, 1] - r['omega_r']) < 1e-5)):
                hits.append(k)
        assert hits, ('Unmatched corrected neutral seed', j)
        matched.append(dict(seed_index=j, curve_indices=hits))
    js(DEST / 'numerical_audit.json', dict(curves=audit, matched_seeds=matched,
       all_corrected_seeds_covered=True, protected_inputs=snapshot(),
       scope='Original matrices and actual-profile tips; current grids, zero explicit curvature; tracked representative and resolved pockets'))
    return curves


def plot(curves):
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.fonttype': 'none', 'axes.formatter.useoffset': False})
    bound = json.loads((ROOT / 'outputs/rough_positive_beta_limit_refined.json').read_text())['critical']['beta_r']
    xmax = np.ceil(max(c['rows'][:, 0].max() for c in curves) / 50) * 50 + 5
    for physical in [False, True]:
        column = 2 if physical else 1
        values = np.concatenate([c['rows'][:, column] for c in curves])
        tick = 50 if physical else 5
        ylim = (tick * np.floor((values.min() - .5) / tick), tick * np.ceil((values.max() + .5) / tick))
        fig, axes = plt.subplots(2, 2, figsize=(14.5, 10.5))
        fig.subplots_adjust(left=.075, right=.985, top=.90, bottom=.14, wspace=.19, hspace=.32)
        plotted = []
        for col, case in enumerate(CASES):
            for row, xlim, title in [(0, (-36, xmax), 'full curve extent'), (1, (-36, 170), 'upstream detail')]:
                ax = axes[row, col]
                for j, c in enumerate(curves):
                    if c['case'] != case:
                        continue
                    a = c['rows']
                    ax.plot(a[:, 0], a[:, column], color=COLORS[c['beta']], lw=1.8 if c['beta'] == -7 else 1.5)
                    if row == 0:
                        assert a[:, 0].min() >= xlim[0] and a[:, 0].max() <= xlim[1]
                    plotted.append(dict(case=case, view=title, curve_index=j))
                ax.set(xlim=xlim, ylim=ylim, xlabel='Global x (mm)',
                       ylabel='Physical f (kHz)' if physical else r'$\omega$ (nondimensional)',
                       title=LABELS[col] + ' / ' + title)
                ax.grid(alpha=.18)
                if case == 'rough' and row == 0:
                    ax.text(.98, .97, f'No verified contour at beta = +20\nResolved branch upper beta ~ {bound:.2f}',
                            ha='right', va='top', transform=ax.transAxes, fontsize=9,
                            bbox=dict(facecolor='white', edgecolor='none', alpha=.85))
        fig.suptitle('M3 x-frequency neutral curves at seven fixed spanwise wavenumbers', fontsize=15)
        fig.text(.5, .935, 'Signed frequencies; separate pockets retained; same beta colors in both geometries', ha='center', fontsize=10)
        fig.legend(handles=[Line2D([0], [0], color=COLORS[b], lw=2, label=rf'$\beta={b:g}$') for b in BETAS],
                   loc='lower center', bbox_to_anchor=(.5, .055), ncol=7, frameon=False)
        fig.text(.5, .026, '15 verified closed contours. No added conjugate copies. Zero explicit curvature in the current local-wall LST model.', ha='center', fontsize=9)
        fig.text(.5, .009, 'Subgrid turning points and dense samples use reassembled profiles; they are not a CFD grid-convergence study.', ha='center', fontsize=9)
        suffix = '_physical_kHz' if physical else '_comparison'
        for ext in ['png', 'svg']:
            fig.savefig(DEST / (NAME + suffix + '.' + ext), dpi=220, facecolor='white')
        plt.close(fig)
        js(DEST / ('plot_validation_physical.json' if physical else 'plot_validation.json'),
           dict(plotted=plotted, full_view_contains_all_points=True, x_max=xmax, ylim=ylim))


def dat(name, curves):
    with (DEST / (name + '.dat')).open('w') as stream:
        stream.write('TITLE="Fixed beta M3 x-omega neutral curves"\nVARIABLES=' + ','.join('"' + v + '"' for v in VARS) + '\n')
        for c in curves:
            a = c['rows']
            stream.write(f'ZONE T="{c["zone"]}", I={len(a)}, DATAPACKING=POINT\nDT=(' + ','.join(['DOUBLE'] * len(VARS)) + ')\n')
            np.savetxt(stream, a, fmt='%.17g')


def macros():
    out = ['#!MC 1410']
    colors = ['PURPLE', 'BLUE', 'CYAN', 'BLACK', 'GREEN', 'YELLOW', 'RED']
    for name in NAMES:
        _, zones = parse(DEST / (name + '.dat'))
        out += [f'$!READDATASET \'"{name}.plt"\'', ' READDATAOPTION=NEW', ' RESETSTYLE=YES',
                ' INITIALPLOTTYPE=XYLINE', '$!PLOTTYPE=XYLINE', '$!DELETELINEMAPS']
        for j, z in enumerate(zones, 1):
            color = colors[BETAS.index(z['rows'][0][3])]
            pattern = 'SOLID' if z['rows'][0][9] == 1 else 'DASHED'
            out += ['$!CREATELINEMAP', f'$!LINEMAP [{j}]', f" NAME='{z['name']}'",
                    f' ASSIGN{{ZONE={j} XAXISVAR=1 YAXISVAR=2}}',
                    f' LINES{{SHOW=YES COLOR={color} LINEPATTERN={pattern}}}', ' SYMBOLS{SHOW=NO}']
        out += ['$!ACTIVELINEMAPS=[' + ','.join(map(str, range(1, len(zones) + 1))) + ']', '$!VIEW FIT',
                f'$!SAVELAYOUT "{name}.lay"', f'$!WRITEDATASET "{name}_roundtrip.dat"',
                ' BINARY=NO', ' USEPOINTFORMAT=YES', ' PRECISION=17', ' INCLUDETEXT=NO',
                ' INCLUDEGEOM=NO', ' INCLUDECUSTOMLABELS=NO']
        for j in range(1, len(zones) + 1):
            out.append(f'$!LINEMAP [{j}] ASSIGN{{YAXISVAR=3}}')
        out += ['$!VIEW FIT', f'$!SAVELAYOUT "{name}_physical_kHz.lay"']
    out.append('$!QUIT')
    (DEST / 'tecplot_roundtrip.mcr').write_text('\n'.join(out) + '\n', encoding='ascii')


def verify():
    assert 'ERROR' not in (DEST / 'batch.log').read_text().upper()
    files = []
    for name in NAMES:
        source, back, binary = [DEST / (name + s) for s in ['.dat', '_roundtrip.dat', '.plt']]
        va, za = parse(source)
        vb, zb = parse(back)
        assert va == vb == VARS and za == zb and binary.read_bytes()[:8] == b'#!TDV112'
        assert (DEST / (name + '.lay')).exists() and (DEST / (name + '_physical_kHz.lay')).exists()
        files.append(dict(file=binary.name, sha256=sha(binary), source_sha256=sha(source), roundtrip_sha256=sha(back),
                          zones=len(za), points=sum(len(z['rows']) for z in za), all_values_exact=True))
    js(DEST / 'tecplot_validation.json', dict(tool='Tecplot 360 EX 2024 R1/Preplot', files=files))
    print('NATIVE_VERIFIED', [(r['file'], r['zones'], r['points']) for r in files], flush=True)


def main():
    if '--verify' in sys.argv:
        verify()
        return
    DEST.mkdir(exist_ok=True)
    curves = load()
    summary = []
    component_count = {}
    for c in curves:
        key = c['case'], c['beta']
        component_count[key] = component_count.get(key, 0) + 1
        c['zone'] = f'{c["case"]}_beta{c["beta"]:+g}_component{component_count[key]}'
        np.savetxt(DEST / (c['zone'] + '.csv'), c['rows'], delimiter=',', header=','.join(VARS), comments='', fmt='%.17g')
        a = c['rows']
        summary.append({**{k: v for k, v in c.items() if k != 'rows'}, 'points': len(a),
                        'x_range_mm': [float(a[:,0].min()), float(a[:,0].max())],
                        'omega_range': [float(a[:,1].min()), float(a[:,1].max())]})
    js(DEST / 'selection_manifest.json', dict(time=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
       betas=BETAS, curves=summary, curve_count=len(curves), points=sum(len(c['rows']) for c in curves),
       missing=[dict(case='rough', beta=20., reason='No verified neutral contour on the tracked representative; resolved beta upper bound supplemented separately')],
       no_added_conjugate=True))
    dat(NAME, curves)
    for case in CASES:
        dat(case + '_' + NAME, [c for c in curves if c['case'] == case])
    plot(curves)
    macros()
    shutil.copy2(__file__, DEST / 'publish.py')
    shutil.copy2(ROOT / 'outputs/rough_positive_beta_limit_refined.json', DEST / 'optimized_positive_beta_bound.json')
    print('PUBLISHED', len(curves), sum(len(c['rows']) for c in curves), flush=True)


if __name__ == '__main__':
    main()
