"""Combine every previously delivered signed-frequency M3 curve, without refitting."""
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

ROOT = Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°')
DEST = ROOT / 'outputs/M3_all_neutral_curves_20261027'
EXT = ROOT / 'outputs/M3_xbeta_extended_frequencies_20261027'
LIM = ROOT / 'outputs/M3_frequency_limits_20261027'
CASES = ['bluntplate', 'rough']
LABELS = ['Baseline blunt plate', 'Optimized nose: continuous curvature']
NAME = 'M3_all_neutral_curves'
NAMES = [NAME] + [case + '_M3_all_neutral_curves' for case in CASES]
VARS = ['x_mm', 'beta', 'f_kHz', 'omega_nd', 'alpha_r', 'alpha_i', 'growth',
        'x_flat_mm', 'station_coordinate', 'point_kind', 'case_id', 'right_residual', 'raw_MAC']
COLORS = {-240.2: '#0a3154', -237.14: '#385a7c', -200.: '#007b8a', -150.: '#498d6d',
          -112.: '#00a6a6', -110.: '#4c63d2', -100.: '#762a83', -75.: '#5954a6',
          -50.: '#2166ac', -25.: '#149fb3', 0.: '#161616', 25.: '#658e28',
          50.: '#1b9e77', 75.: '#e48120', 100.: '#d73027', 125.: '#8c564b',
          150.: '#e377c2', 175.: '#bc8f00', 200.: '#7f3c8d', 230.: '#006d77',
          250.: '#a44a3f', 264.53: '#a52a6a', 275.: '#d29d00', 290.: '#e65b00', 291.6: '#c51b7d'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def js(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')


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


def checked_csv(folder, name, manifest):
    path = folder / name
    record = next(r for r in manifest['files'] if r['path'] == name)
    assert sha(path) == record['sha256'], str(path)
    with path.open() as stream:
        reader = csv.DictReader(stream)
        assert reader.fieldnames == VARS
        rows = np.array([[float(row[k]) for k in VARS] for row in reader])
    return rows, dict(path=str(path), sha256=sha(path), bytes=path.stat().st_size)


def load():
    ext_selection = json.loads((EXT / 'selection_manifest.json').read_text())
    lim_selection = json.loads((LIM / 'selection_manifest.json').read_text())
    ext_delivery = json.loads((EXT / 'delivery_manifest.json').read_text())
    lim_delivery = json.loads((LIM / 'delivery_manifest.json').read_text())
    records, duplicates = {}, []
    for r in ext_selection['source_files']:
        case, f = r['case'], float(r['f_kHz'])
        rows, source = checked_csv(EXT, f'{case}_f{f:+g}kHz.csv', ext_delivery)
        region = 'nose' if f == -112 else ('positive' if f > 0 else 'stationary' if f == 0 else 'tracked')
        assert len(rows) == r['points']
        records[case, f, region] = dict(case=case, f_kHz=f, region=region, rows=rows, source=source)
    for r in lim_selection['curves']:
        case, f, region = r['case'], float(r['f_kHz']), r['region']
        rows, source = checked_csv(LIM, f'{case}_{region}_f{f:+g}kHz.csv', lim_delivery)
        key = case, f, region
        if key in records:
            assert np.array_equal(rows, records[key]['rows']), 'Duplicate references differ'
            duplicates.append(dict(case=case, f_kHz=f, region=region, source=source,
                                   retained=records[key]['source']))
        else:
            records[key] = dict(case=case, f_kHz=f, region=region, rows=rows, source=source)
    curves = sorted(records.values(), key=lambda c: (CASES.index(c['case']), c['f_kHz'], c['region']))
    assert len(curves) == 45 and len(duplicates) == 4
    for j, c in enumerate(curves, 1):
        a, f = c['rows'], c['f_kHz']
        c['id'] = j
        c['zone'] = f'{c["case"]}_{c["region"]}_f{f:+g}kHz'
        assert np.isfinite(a).all() and np.max(abs(a[:, 5])) < 1e-8
        assert np.all(a[:, 2] == f) and np.all(a[:, 10] == CASES.index(c['case']) + 1)
        assert np.allclose(a[:, 3], 6.0727867021443425 * f / 50, rtol=1e-13, atol=1e-13)
        if f:
            assert np.array_equal(a[0, :4], a[-1, :4]) and a[-1, 9] == 2
        else:
            assert abs(a[0, 0] - 1100) < 1e-8 and abs(a[-1, 0] - 1100) < 1e-8
            assert not np.array_equal(a[0, :2], a[-1, :2])
        target = DEST / 'curves' / (c['zone'] + '.csv')
        target.parent.mkdir(exist_ok=True)
        shutil.copy2(c['source']['path'], target)
        assert sha(target) == c['source']['sha256']
        c['output_csv'] = str(target.relative_to(DEST))
    source_dat = LIM / 'M3_frequency_limits_neutral_curves.dat'
    entry = next(r for r in lim_delivery['files'] if r['path'] == source_dat.name)
    assert sha(source_dat) == entry['sha256']
    variables, zones = parse(source_dat)
    assert variables == VARS
    limits = [z for z in zones if z['rows'][0][9] == 9]
    assert len(limits) == 4
    return curves, limits, duplicates, source_dat


def styles(frequencies):
    result = {0.: '-'}
    for sign in [-1, 1]:
        fs = sorted([f for f in frequencies if np.sign(f) == sign])
        result.update({f: ['-', '--', '-.'][i % 3] for i, f in enumerate(fs)})
    return result


def save(fig, name):
    for ext in ['png', 'svg']:
        fig.savefig(DEST / (name + '.' + ext), dpi=220, facecolor='white')
    plt.close(fig)


def plot(curves, limits):
    frequencies = sorted({c['f_kHz'] for c in curves})
    assert frequencies == sorted(COLORS)
    style = styles(frequencies)
    available = {f: {c['case'] for c in curves if c['f_kHz'] == f} for f in frequencies}
    labels = {f: ('0' if f == 0 else f'{f:+g}') +
              (' (B)' if available[f] == {'bluntplate'} else ' (O)' if available[f] == {'rough'} else '')
              for f in frequencies}
    handles = [Line2D([0], [0], color=COLORS[f], ls=style[f], lw=1.8, label=labels[f]) for f in frequencies]
    ymin = min(c['rows'][:, 1].min() for c in curves)
    ymax = max(c['rows'][:, 1].max() for c in curves)
    ylim = (5 * np.floor((ymin - 1) / 5), 5 * np.ceil((ymax + 1) / 5))
    full_xlim = (-36, 1120)
    receipt = []
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.fonttype': 'none', 'axes.formatter.useoffset': False})
    fig, axes = plt.subplots(2, 2, figsize=(16, 12.8))
    fig.subplots_adjust(left=.065, right=.985, top=.91, bottom=.22, wspace=.17, hspace=.32)
    for col, case in enumerate(CASES):
        cs = [c for c in curves if c['case'] == case]
        for row, xlim, detail in [(0, full_xlim, 'full computed extent'), (1, (-36, 200), 'upstream detail')]:
            ax = axes[row, col]
            for c in cs:
                a, f = c['rows'], c['f_kHz']
                ax.plot(a[:, 0], a[:, 1], color=COLORS[f], ls=style[f], lw=1.8 if f == 0 else 1.25,
                        zorder=4 if f == 0 else 2)
                if row == 0:
                    assert a[:, 0].min() >= xlim[0] and a[:, 0].max() <= xlim[1]
                    assert a[:, 1].min() >= ylim[0] and a[:, 1].max() <= ylim[1]
            for limit in limits:
                a = np.array(limit['rows'])
                if a[0, 10] == CASES.index(case) + 1:
                    ax.plot(a[:, 0], a[:, 1], 'k*', ms=8, zorder=6)
            ax.set(xlim=xlim, ylim=ylim, xlabel='Global x (mm)', ylabel=r'$\beta$ (nondimensional)',
                   title=LABELS[col] + ' / ' + detail)
            ax.grid(alpha=.18)
            receipt.append(dict(figure='comparison', case=case, view=detail,
                                plotted_curve_ids=[c['id'] for c in cs], xlim=xlim, ylim=ylim))
    fig.suptitle('All computed M3 x-beta neutral curves: both frequency signs together', fontsize=16)
    fig.text(.5, .945, '45 curves at 25 signed physical frequencies; same frequency colors in both geometries',
             ha='center', fontsize=11)
    fig.legend(handles=handles, title='f (kHz)     B: baseline only; O: optimized only',
               loc='lower center', bbox_to_anchor=(.5, .052), ncol=5, frameon=False, fontsize=10,
               handlelength=3.5, columnspacing=2.8)
    fig.text(.5, .023, 'Every delivered curve is retained; separate pockets remain separate. Black stars: numerical frequency limits.',
             ha='center', fontsize=9)
    fig.text(.5, .006, 'No conjugate copies added. f = 0 is open at x = 1100 mm. Same zero-explicit-curvature LST model.',
             ha='center', fontsize=9)
    save(fig, NAME + '_comparison')
    # A literal superposition of both geometries is supplied as an additional figure.
    fig, axes = plt.subplots(1, 2, figsize=(16, 8.5))
    fig.subplots_adjust(left=.065, right=.985, top=.84, bottom=.32, wspace=.17)
    for ax, xlim, detail in zip(axes, [full_xlim, (-36, 200)], ['full computed extent', 'upstream detail']):
        for c in curves:
            a, f = c['rows'], c['f_kHz']
            ax.plot(a[:, 0], a[:, 1], color=COLORS[f], ls='-' if c['case'] == 'bluntplate' else '--',
                    lw=1.3, alpha=.9)
        ax.set(xlim=xlim, ylim=ylim, xlabel='Global x (mm)', ylabel=r'$\beta$ (nondimensional)', title=detail)
        ax.grid(alpha=.18)
        receipt.append(dict(figure='both_geometries_overlay', view=detail,
                            plotted_curve_ids=[c['id'] for c in curves], xlim=xlim, ylim=ylim))
    fig.suptitle('All 45 computed curves superposed: frequency color, geometry line style', fontsize=16)
    fig.text(.5, .885, 'Solid: baseline blunt plate     Dashed: optimized continuous-curvature nose', ha='center', fontsize=11)
    color_handles = [Line2D([0], [0], color=COLORS[f], lw=2, label=labels[f]) for f in frequencies]
    fig.legend(handles=color_handles, title='f (kHz)     B: baseline only; O: optimized only', loc='lower center',
               bbox_to_anchor=(.5, .075), ncol=5, frameon=False, fontsize=10, handlelength=3.5, columnspacing=2.8)
    fig.text(.5, .027, 'All supplied curves; no added conjugate copies or interpolated frequencies. f = 0 remains open at the domain end.',
             ha='center', fontsize=9)
    save(fig, NAME + '_both_geometries_overlay')
    js(DEST / 'plot_validation.json', dict(frequencies_kHz=frequencies, axes=receipt,
                                         all_45_curves_in_full_views=True, full_view_bounds_cover_all_points=True))


def dat(name, zones):
    with (DEST / (name + '.dat')).open('w') as stream:
        stream.write('TITLE="All previously delivered M3 x-beta curves, both frequency signs"\nVARIABLES=' +
                     ','.join('"' + v + '"' for v in VARS) + '\n')
        for z in zones:
            a = np.asarray(z['rows'])
            stream.write(f'ZONE T="{z["name"]}", I={len(a)}, DATAPACKING=POINT\nDT=(' +
                         ','.join(['DOUBLE'] * len(VARS)) + ')\n')
            np.savetxt(stream, a, fmt='%.17g')


def macros():
    out = ['#!MC 1410']
    palette = ['BLUE', 'CYAN', 'GREEN', 'PURPLE', 'RED', 'BLACK']
    for name in NAMES:
        _, zones = parse(DEST / (name + '.dat'))
        frequencies = sorted({z['rows'][0][2] for z in zones if z['rows'][0][9] != 9})
        out += [f'$!READDATASET \'"{name}.plt"\'', ' READDATAOPTION=NEW', ' RESETSTYLE=YES',
                ' INITIALPLOTTYPE=XYLINE', '$!PLOTTYPE=XYLINE', '$!DELETELINEMAPS']
        for j, z in enumerate(zones, 1):
            critical = z['rows'][0][9] == 9
            f = z['rows'][0][2]
            k = frequencies.index(f) if not critical else 0
            color = palette[k % len(palette)] if not critical else 'BLACK'
            pattern = 'SOLID' if z['rows'][0][10] == 1 else 'DASHED'
            out += ['$!CREATELINEMAP', f'$!LINEMAP [{j}]', f" NAME='{z['name']}'",
                    f' ASSIGN{{ZONE={j} XAXISVAR=1 YAXISVAR=2}}',
                    f' LINES{{SHOW={"NO" if critical else "YES"} COLOR={color} LINEPATTERN={pattern}}}',
                    f' SYMBOLS{{SHOW={"YES" if critical else "NO"}}}']
        out += ['$!ACTIVELINEMAPS=[' + ','.join(map(str, range(1, len(zones) + 1))) + ']', '$!VIEW FIT',
                f'$!SAVELAYOUT "{name}.lay"', f'$!WRITEDATASET "{name}_roundtrip.dat"',
                ' BINARY=NO', ' USEPOINTFORMAT=YES', ' PRECISION=17', ' INCLUDETEXT=NO',
                ' INCLUDEGEOM=NO', ' INCLUDECUSTOMLABELS=NO']
    out.append('$!QUIT')
    (DEST / 'tecplot_roundtrip.mcr').write_text('\n'.join(out) + '\n', encoding='ascii')


def verify():
    assert 'ERROR' not in (DEST / 'batch.log').read_text().upper()
    receipts = []
    for name in NAMES:
        source, back, binary = [DEST / (name + suffix) for suffix in ['.dat', '_roundtrip.dat', '.plt']]
        va, za = parse(source)
        vb, zb = parse(back)
        assert va == vb == VARS and za == zb
        assert binary.read_bytes()[:8] == b'#!TDV112' and (DEST / (name + '.lay')).exists()
        receipts.append(dict(file=binary.name, sha256=sha(binary), source_sha256=sha(source),
                             roundtrip_sha256=sha(back), zones=len(za), points=sum(len(z['rows']) for z in za),
                             all_values_exact=True))
    js(DEST / 'tecplot_validation.json', dict(tool='Tecplot 360 EX 2024 R1 / Preplot', files=receipts))
    print('NATIVE_VERIFIED', [(r['file'], r['zones'], r['points']) for r in receipts], flush=True)


def main():
    if '--verify' in sys.argv:
        verify()
        return
    DEST.mkdir(exist_ok=True)
    curves, limits, duplicates, limit_dat = load()
    records = []
    for c in curves:
        a = c['rows']
        records.append({**{k: v for k, v in c.items() if k != 'rows'}, 'points': len(a),
                        'closed': bool(c['f_kHz']), 'x_range_mm': [float(a[:, 0].min()), float(a[:, 0].max())],
                        'beta_range': [float(a[:, 1].min()), float(a[:, 1].max())]})
    js(DEST / 'selection_manifest.json', dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
       curves=records, curve_count=len(curves), curve_points=sum(c['points'] for c in records),
       critical_points=len(limits), frequency_count=len({c['f_kHz'] for c in records}),
       duplicates_removed=duplicates, critical_source=dict(path=str(limit_dat), sha256=sha(limit_dat)),
       scope='All previously delivered curves of the tracked M3 representative and resolved pockets; no conjugate copy added',
       operation='Exact data combination only; no new eigensolutions, curve fitting or frequency interpolation'))
    zones = [dict(name=c['zone'], rows=c['rows']) for c in curves] + limits
    dat(NAME, zones)
    for i, case in enumerate(CASES, 1):
        dat(case + '_M3_all_neutral_curves', [z for z in zones if z['rows'][0][10] == i])
    plot(curves, limits)
    macros()
    shutil.copy2(__file__, DEST / 'publish_all.py')
    shutil.copy2(LIM / 'frequency_limits.json', DEST / 'frequency_limits.json')
    print('PUBLISHED', len(curves), 'curves', sum(len(c['rows']) for c in curves), 'curve points; 4 critical points', flush=True)


if __name__ == '__main__':
    main()
