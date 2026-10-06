#!/usr/bin/env python3
"""
check_gxtb_logs.py - summarize Gaussian (+ xtb-gaussian / g-xTB) jobs in a folder.

For every .log it reports: termination, optimization convergence, number of steps,
imaginary frequencies, energies (E, H, G), charge/multiplicity, xTB flags, and whether
the metal lost/gained ligands or any other bond broke/formed during the optimization.
It groups isomers/spin states (names like <system>_<kind>_<ArX>_o3_m2) and gives
relative free energies in kcal/mol, and can write restart / displaced inputs.

Usage (run in the folder with the .log/.com files):
    python check_gxtb_logs.py                      # summary table + check_summary.csv
    python check_gxtb_logs.py --restart            # new inputs for unfinished/unconverged jobs
    python check_gxtb_logs.py --displace 0.15      # new inputs pushed off imaginary modes (minima)
    python check_gxtb_logs.py path/to/folder --restart --displace 0.15

New inputs go to ./resubmit/ with filenames_resubmit.txt (Unix line endings).
Only numpy is required.
"""
import argparse, csv, glob, os, re, sys
import numpy as np

HARTREE2KCAL = 627.5095
SYMBOLS = ("X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As "
           "Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb "
           "Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn").split()
# covalent radii (Angstrom, Alvarez 2008; low-spin values for transition metals)
RCOV = dict(H=0.31, B=0.84, C=0.76, N=0.71, O=0.66, F=0.57, Si=1.11, P=1.07, S=1.05, Cl=1.02, Se=1.20,
            Br=1.20, I=1.39, Sc=1.70, Ti=1.60, V=1.53, Cr=1.39, Mn=1.39, Fe=1.32, Co=1.26, Ni=1.24,
            Cu=1.32, Zn=1.22, Y=1.90, Zr=1.75, Nb=1.64, Mo=1.54, Tc=1.47, Ru=1.46, Rh=1.42, Pd=1.39,
            Ag=1.45, Cd=1.44, Hf=1.75, Ta=1.70, W=1.62, Re=1.51, Os=1.44, Ir=1.41, Pt=1.36, Au=1.36, Hg=1.32)
METALS = {"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
          "Rh", "Pd", "Ag", "Cd", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"}
NAME_RE = re.compile(r"^(?P<group>.+)_o(?P<iso>\d+)_m(?P<mult>\d+)(?:_(?:r|d)\d+)*$")

# ----------------------------------------------------------------- parsing

def route_lines(text):
    """All route sections ('#' lines up to the dashed line), joined across Gaussian's line wrapping."""
    routes, lines = [], text.splitlines()
    for i, l in enumerate(lines):
        if l.startswith(" #"):
            block = [l[1:]]
            for nxt in lines[i + 1:i + 8]:
                if nxt.startswith(" ---"):
                    break
                block.append(nxt[1:])
            routes.append("".join(block))
    return routes

def geometries(text, frame=None):
    """List of (Z array, coords) for every printed orientation (Input preferred, else Standard),
    or only the given frame ("Input"/"Standard")."""
    pat = r"{} orientation:.*?-{{60,}}\n.*?-{{60,}}\n(.*?)\n -{{60,}}"
    if frame:
        blocks = re.findall(pat.format(frame), text, re.S)
    else:
        blocks = re.findall(pat.format("Input"), text, re.S) or re.findall(pat.format("Standard"), text, re.S)
    out = []
    for b in blocks:
        rows = [r.split() for r in b.splitlines()]
        out.append((np.array([int(r[1]) for r in rows]), np.array([[float(v) for v in r[3:6]] for r in rows])))
    return out

def last_freq_job(text):
    """Text of the last frequency calculation (after the last 'Harmonic frequencies' header)."""
    i = text.rfind("Harmonic frequencies")
    return text[i:] if i >= 0 else ""

def frequencies(ftext):
    vals = []
    for m in re.finditer(r"Frequencies --\s+(.*)", ftext):
        vals += [float(v) for v in m.group(1).split()]
    return vals

def first_mode_vector(ftext, natoms):
    """Displacement vector of the first (lowest) normal mode in the last frequency job."""
    i = ftext.find("Frequencies --")
    if i < 0:
        return None
    lines = ftext[i:].splitlines()
    k = next((n for n, l in enumerate(lines) if l.split()[:2] == ["Atom", "AN"]), None)
    if k is None:
        return None
    vec = [[float(v) for v in l.split()[2:5]] for l in lines[k + 1:k + 1 + natoms]]
    return np.array(vec) if len(vec) == natoms else None

def last_float(pattern, text):
    m = re.findall(pattern, text)
    return float(m[-1]) if m else None

def bonds(Z, x, scale=1.15):
    el = [SYMBOLS[z] for z in Z]
    D = np.linalg.norm(x[:, None] - x[None], axis=-1)
    out = set()
    for i in range(len(el)):
        for j in range(i + 1, len(el)):
            if el[i] == "H" and el[j] == "H":
                continue
            ri, rj = RCOV.get(el[i], 1.5), RCOV.get(el[j], 1.5)
            s = 1.3 if (el[i] in METALS or el[j] in METALS) else scale
            if D[i, j] < s * (ri + rj):
                out.add((i, j))
    return out

def label(Z, i):
    return f"{SYMBOLS[Z[i]]}{i + 1}"

# ----------------------------------------------------------------- one log

def analyze(path):
    text = open(path, errors="ignore").read()
    name = os.path.splitext(os.path.basename(path))[0]
    r = {"name": name, "file": path}
    routes = route_lines(text)
    route = " | ".join(routes).lower()
    r["route"] = routes[0] if routes else ""
    r["jobtype"] = ("irc" if "irc" in route else "ts" if re.search(r"opt\s*=?\s*\(?[^|]*\bts\b", route)
                    else "scan" if "modredundant" in route and re.search(r"\n\s*[BADL]\s.*\bS\b", text) else
                    "opt" if "opt" in route else "freq" if "freq" in route else "sp")
    ext = re.search(r'Running external command "([^"]+)"', text)
    r["xtb_cmd"] = ext.group(1) if ext else ""
    flags = re.sub(r"--charge\s+\S+|--chrg\s+\S+|--uhf\s+\S+|-P\s+\S+|\sR$|xtb-gaussian", "", r["xtb_cmd"]).split()
    r["xtb_flags"] = " ".join(flags)
    cm = re.search(r"Charge =\s*(-?\d+)\s+Multiplicity =\s*(\d+)", text)
    r["charge"], r["mult"] = (int(cm.group(1)), int(cm.group(2))) if cm else (None, None)
    uhf = re.search(r"--uhf\s+(\d+)", r["xtb_cmd"])
    chg = re.search(r"--(?:charge|chrg)\s+(-?\d+)", r["xtb_cmd"])
    r["spin_mismatch"] = bool(uhf and r["mult"] and int(uhf.group(1)) != r["mult"] - 1) or \
                         bool(chg and r["charge"] is not None and int(chg.group(1)) != r["charge"])

    r["normal"] = text.count("Normal termination")
    r["error"] = "Error termination" in text
    r["xtb_failed"] = ("abnormal termination of xtb" in text or "FIO-F-217" in text
                       or bool(re.search(r"Recovered energy=\s+0\.0+\s", text)))
    r["steps"] = text.count("Step number")
    r["stationary"] = "Stationary point found" in text
    r["opt_stopped"] = "Optimization stopped" in text or "Number of steps exceeded" in text
    conv = re.findall(r"(Maximum Force|RMS     Force|Maximum Displacement|RMS     Displacement)\s+\S+\s+\S+\s+(YES|NO)", text)
    r["last_criteria"] = "".join("Y" if c[1] == "YES" else "n" for c in conv[-4:]) if conv else ""
    r["irc_minima"] = text.count("PES minimum detected")
    r["scan_points"] = text.count("Optimization completed") if r["jobtype"] == "scan" else None

    ftext = last_freq_job(text)
    fr = frequencies(ftext)
    r["n_freq"] = len(fr)
    r["n_imag"] = sum(1 for f in fr if f < 0)
    r["imag"] = ", ".join(f"{f:.0f}" for f in fr if f < 0)
    r["lowest_freq"] = min(fr) if fr else None
    r["E"] = last_float(r"Recovered energy=\s+(-?\d+\.\d+)", text) or last_float(r"SCF Done:.*?=\s+(-?\d+\.\d+)", text)
    r["ZPE"] = last_float(r"Zero-point correction=\s+(-?\d+\.\d+)", text)
    r["H"] = last_float(r"Sum of electronic and thermal Enthalpies=\s+(-?\d+\.\d+)", text)
    r["G"] = last_float(r"Sum of electronic and thermal Free Energies=\s+(-?\d+\.\d+)", text)
    el_t = re.findall(r"Elapsed time:\s+(\d+) days\s+(\d+) hours\s+(\d+) minutes\s+([\d.]+)", text)
    r["minutes"] = round(sum(int(d) * 1440 + int(h) * 60 + int(m) + float(s) / 60 for d, h, m, s in el_t), 1) if el_t else None

    geos = geometries(text)
    r["_geos"] = geos
    r["_freqtext"] = ftext
    # normal modes are printed in the standard orientation (input orientation if nosymm): keep that frame
    before = text[:text.rfind("Harmonic frequencies")] if ftext else ""
    std = geometries(before, "Standard") or geometries(before, "Input")
    r["_mode_geom"] = std[-1][1] if std else None
    nm = NAME_RE.match(name)
    r["name_mult_mismatch"] = bool(nm and r["mult"] and int(nm.group("mult")) != r["mult"])
    r["metal_change"] = r["bond_change"] = ""
    if len(geos) >= 2:
        (Z, x0), (_, x1) = geos[0], geos[-1]
        b0, b1 = bonds(Z, x0), bonds(Z, x1)
        lost, formed = b0 - b1, b1 - b0
        is_m = lambda b: SYMBOLS[Z[b[0]]] in METALS or SYMBOLS[Z[b[1]]] in METALS
        fmt = lambda bs: ", ".join(f"{label(Z, i)}-{label(Z, j)}" for i, j in sorted(bs))
        ml = [b for b in lost if is_m(b)]; mf = [b for b in formed if is_m(b)]
        ol = [b for b in lost if not is_m(b)]; of = [b for b in formed if not is_m(b)]
        r["metal_change"] = "; ".join(s for s in (f"lost {fmt(ml)}" if ml else "", f"new {fmt(mf)}" if mf else "") if s)
        r["bond_change"] = "; ".join(s for s in (f"broken {fmt(ol)}" if ol else "", f"formed {fmt(of)}" if of else "") if s)
        metal = next((i for i, z in enumerate(Z) if SYMBOLS[z] in METALS), None)
        if metal is not None:
            r["metal_CN"] = sum(1 for b in b1 if metal in b)
            # notebook numbering: atom 2 = ipso C, atom 3 = leaving-group atom (if bound to the metal)
            if (len(Z) > 3 and SYMBOLS[Z[1]] == "C" and SYMBOLS[Z[2]] in ("Cl", "Br", "I", "O")
                    and any(metal in b and 2 in b for b in b0 | b1)):
                r["C2_X3"] = round(float(np.linalg.norm(x1[1] - x1[2])), 3)
    r["status"] = classify(r)
    return r

def classify(r):
    if r["xtb_failed"]:
        return "XTB_FAILED"
    if r["spin_mismatch"]:
        return "CHARGE/SPIN_MISMATCH"
    if r["error"]:
        if r["opt_stopped"]:
            return "NOT_CONVERGED"
        return "ERROR"
    if r["normal"] == 0:
        return "INCOMPLETE"          # still running, or killed (time limit / preemption)
    if r["jobtype"] in ("irc", "scan", "sp"):
        status = "OK"
    elif r["jobtype"] in ("opt", "ts", "freq"):
        if r["jobtype"] != "freq" and not r["stationary"]:
            return "NOT_CONVERGED"
        if r["n_freq"] == 0:
            status = "OK_NO_FREQ"
        elif r["jobtype"] == "ts":
            status = "OK" if r["n_imag"] == 1 else f"TS_NIMAG={r['n_imag']}"
        else:
            status = "OK" if r["n_imag"] == 0 else f"IMAGINARY({r['n_imag']})"
    else:
        status = "OK"
    if status.startswith("OK") and (r["metal_change"] or r["bond_change"]) and r["jobtype"] != "irc":
        status += "+CHECK_BONDS"
    return status

# ----------------------------------------------------------------- new inputs

def read_com_template(name, folder):
    """Header (through charge/mult line) and tail (after coordinates) of the original .com, if present."""
    path = os.path.join(folder, name + ".com")
    if not os.path.isfile(path):
        return None
    lines = open(path).read().replace("\r\n", "\n").split("\n")
    cm = next(i for i, l in enumerate(lines) if re.match(r"^\s*-?\d+\s+\d+\s*$", l))
    end = next((i for i in range(cm + 1, len(lines)) if not lines[i].strip()), len(lines))
    return lines[:cm + 1], lines[end:]

def write_input(r, folder, outdir, newname, x_new, note):
    tpl = read_com_template(r["name"], folder)
    if tpl is None:
        print(f"  {r['name']}: no {r['name']}.com found, cannot write {newname}.com")
        return None
    head, tail = tpl
    head = [re.sub(r"%chk=.*", f"%chk={newname}.chk", l, flags=re.I) for l in head]
    tail = [re.sub(r"%chk=.*", f"%chk={newname}.chk", l, flags=re.I) for l in tail]   # --Link1-- steps too
    # every xtb-gaussian call in the file must use the same charge and --uhf
    whole = "\n".join(head + tail)
    chg = set(re.findall(r"--(?:charge|chrg)\s+(-?\d+)", whole)); uhf = set(re.findall(r"--uhf\s+(\d+)", whole))
    if len(chg) > 1 or len(uhf) > 1:
        print(f"  WARNING {newname}: the input's steps use different --charge {sorted(chg)} / --uhf {sorted(uhf)}; "
              f"fix this in {newname}.com before submitting (e.g. a --Link1-- section with the wrong charge).")
    Z = r["_geos"][-1][0]
    coords = [f"{SYMBOLS[z]:<3s}{a:14.6f}{b:14.6f}{c:14.6f}" for z, (a, b, c) in zip(Z, x_new)]
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, newname + ".com"), "w", newline="\n") as f:
        f.write("\n".join(head + coords + tail).rstrip("\n") + "\n\n")
    print(f"  wrote resubmit/{newname}.com  ({note})")
    return newname + ".com"

def next_name(name, tag):
    m = re.match(rf"^(.*)_{tag}(\d+)$", name)
    return f"{m.group(1)}_{tag}{int(m.group(2)) + 1}" if m else f"{name}_{tag}1"

# ----------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", nargs="?", default=".")
    ap.add_argument("--csv", default="check_summary.csv")
    ap.add_argument("--restart", action="store_true", help="restart NOT_CONVERGED/INCOMPLETE jobs from their last geometry")
    ap.add_argument("--displace", type=float, default=0.0, metavar="A",
                    help="for minima with imaginary modes: push along the lowest mode (max atom shift A, e.g. 0.15)")
    args = ap.parse_args()

    logs = sorted(glob.glob(os.path.join(args.folder, "*.log")))
    logs = [l for l in logs if not os.path.basename(l).startswith("slurm-")]
    if not logs:
        sys.exit(f"No .log files in {args.folder}")
    res = [analyze(l) for l in logs]

    # relative free energies within each isomer/spin group
    groups = {}
    for r in res:
        m = NAME_RE.match(r["name"])
        r["group"] = m.group("group") if m else r["name"]
        if r["status"].startswith("OK") and r["G"] is not None:
            groups.setdefault(r["group"], []).append(r)
    for g, members in groups.items():
        gmin = min(m["G"] for m in members)
        for m in members:
            m["dG_kcal"] = round((m["G"] - gmin) * HARTREE2KCAL, 2)

    # print table
    cols = [("name", 38), ("status", 22), ("job", 5), ("chg", 3), ("mult", 4), ("steps", 5), ("imag", 12),
            ("G_hartree", 15), ("dG", 7), ("notes", 0)]
    print("".join(f"{c:<{w}s} " if w else c for c, w in cols))
    for r in sorted(res, key=lambda r: (r["group"], r.get("dG_kcal", 1e9), r["name"])):
        notes = "; ".join(s for s in (
            f"metal: {r['metal_change']}" if r["metal_change"] else "",
            f"bonds: {r['bond_change']}" if r["bond_change"] else "",
            f"C2-X3 {r['C2_X3']}" if r.get("C2_X3") else "",
            f"IRC minima {r['irc_minima']}" if r["jobtype"] == "irc" else "",
            f"criteria {r['last_criteria']}" if r["status"] == "NOT_CONVERGED" else "",
            f"NAME SAYS m{NAME_RE.match(r['name']).group('mult')} BUT LOG HAS mult {r['mult']}" if r["name_mult_mismatch"] else "") if s)
        vals = [r["name"], r["status"], r["jobtype"], str(r["charge"]), str(r["mult"]), str(r["steps"]),
                r["imag"] or "-", f"{r['G']:.6f}" if r["G"] is not None else "-",
                f"{r['dG_kcal']:.1f}" if "dG_kcal" in r else "-", notes]
        print("".join(f"{v:<{w}s} " if w else v for v, (c, w) in zip(vals, cols)))

    # warnings
    flagsets = {r["xtb_flags"] for r in res if r["xtb_flags"]}
    if len(flagsets) > 1:
        print("\nWARNING: logs use different xTB settings (energies are NOT comparable across them):")
        for f in sorted(flagsets):
            print(f"   [{f}]  " + ", ".join(r["name"] for r in res if r["xtb_flags"] == f)[:200])
    counts = {}
    for r in res:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    print("\nSummary: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
    print("dG = free energy relative to the lowest OK structure with the same name stem "
          "(same species, different isomers o#/spin states m#), kcal/mol.")

    # csv
    keys = ["name", "group", "status", "jobtype", "charge", "mult", "xtb_flags", "normal", "steps", "stationary",
            "last_criteria", "n_imag", "imag", "lowest_freq", "E", "ZPE", "H", "G", "dG_kcal", "metal_CN",
            "metal_change", "bond_change", "C2_X3", "irc_minima", "scan_points", "minutes", "file"]
    with open(os.path.join(args.folder, args.csv), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in res:
            w.writerow(r)
    print(f"Wrote {os.path.join(args.folder, args.csv)}")

    # new inputs
    if args.restart or args.displace:
        outdir = os.path.join(args.folder, "resubmit")
        new = []
        print()
        for r in res:
            if not r["_geos"]:
                continue
            x_last = r["_geos"][-1][1]
            if args.restart and r["status"] in ("NOT_CONVERGED", "INCOMPLETE"):
                if r["jobtype"] == "ts" and r["metal_change"]:
                    print(f"  {r['name']}: NOT restarted - the TS search changed metal-ligand bonding "
                          f"({r['metal_change']}); it has left the saddle region. Build a better guess (scan).")
                    continue
                n = write_input(r, args.folder, outdir, next_name(r["name"], "r"), x_last,
                                f"restart from last geometry, step {r['steps']}")
                if n: new.append(n)
            if args.displace and r["status"].startswith("IMAGINARY") and r["jobtype"] == "opt":
                v = first_mode_vector(r["_freqtext"], len(x_last))
                if v is None:
                    print(f"  {r['name']}: could not read the imaginary mode"); continue
                base = r["_mode_geom"] if r["_mode_geom"] is not None else x_last
                v = v / np.abs(v).max() * args.displace
                n = write_input(r, args.folder, outdir, next_name(r["name"], "d"), base + v,
                                f"displaced {args.displace} A along {r['imag'].split(',')[0]} cm-1 mode")
                if n: new.append(n)
        if new:
            with open(os.path.join(outdir, "filenames_resubmit.txt"), "w", newline="\n") as f:
                f.write("\n".join(new) + "\n")
            print(f"\n{len(new)} new inputs in {outdir}. Submit from there with your SLURM script, e.g.\n"
                  f"   cd {outdir}; sbatch --array=0-{len(new)-1}%20 ../g16_xtb.slurm filenames_resubmit.txt")
        elif args.restart or args.displace:
            print("Nothing to resubmit.")

if __name__ == "__main__":
    main()
