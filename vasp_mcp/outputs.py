"""Summaries of finished, crashed or running VASP calculations.

Reads OUTCAR, OSZICAR, vasprun.xml and the job's standard output. Like the
input readers in vasp-lint, nothing here needs more than the standard library,
and a file that cannot be parsed is skipped rather than raised on.
"""

import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

from vasp_lint.inputs import read_incar, read_poscar

# Messages VASP writes to standard output when a run fails. The selection
# follows the error table of custodian's VaspErrorHandler (MIT licensed).
STDOUT_ERRORS = {
    "BRMIX: very serious problems": "charge-density mixing failed",
    "ZBRENT: fatal": "line search in the ionic relaxation failed",
    "Error EDDDAV: Call to ZHEGV failed": "Davidson diagonalisation failed",
    "WARNING in EDDRMM: call to ZHEGV failed": "RMM-DIIS diagonalisation failed",
    "Sub-Space-Matrix is not hermitian": "Davidson subspace matrix not hermitian",
    "ERROR in subspace rotation PSSYEVX": "ScaLAPACK subspace rotation failed",
    "ERROR in subspace rotation PDSYEVX": "ScaLAPACK subspace rotation failed",
    "LAPACK: Routine ZPOTRF failed": "Cholesky decomposition failed",
    "EDWAV: internal error, the gradient is not orthogonal": "gradient not orthogonal",
    "TOO FEW BANDS": "too few bands",
    "number of bands is not sufficient": "too few bands",
    "Tetrahedron method fails": "tetrahedron method failed",
    "Routine TETIRR needs special values": "tetrahedron method failed",
    "Fatal error detecting k-mesh": "k-mesh symmetry error",
    "Fatal error: unable to match k-point": "k-mesh symmetry error",
    "internal error in subroutine PRICEL": "symmetry error (PRICEL)",
    "rotation matrix was not found": "symmetry error (increase SYMPREC)",
    "Found some non-integer element in rotation matrix": "symmetry error (SGRCON)",
    "Inconsistent Bravais lattice": "symmetry error (Bravais lattice)",
    "POSMAP": "symmetry error (POSMAP)",
    "RHOSYG": "symmetry error (RHOSYG)",
    "BRIONS problems: POTIM should be increased": "BRIONS: POTIM should be increased",
    "ERROR RSPHER": "real-space projection error (RSPHER)",
    "REAL_OPTLAY: internal error": "real-space projection error (REAL_OPTLAY)",
    "ERROR: the triple product of the basis vectors": "lattice vectors are left-handed or degenerate",
    "while reading WAVECAR": "WAVECAR could not be read",
    "Error reading item": "INCAR value could not be read",
    "segmentation fault": "segmentation fault",
    "Killed": "process killed",
    "DUE TO TIME LIMIT": "job hit the wall-time limit",
    "oom-kill": "out of memory",
    "Out Of Memory": "out of memory",
}
STDOUT_NAMES = ("vasp.out", "stdout", "OUTPUT", "out", "log")
STDOUT_GLOBS = ("slurm-*.out", "*.o[0-9]*", "*.out", "*.log")

PARTIAL = (0.01, 0.99)  # occupations that mean a state sits at the Fermi level
BANNER = re.compile(r"^[WARNIG! ]*$")  # the ASCII-art "WARNING" banner inside a box
FLOAT = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?"


def _floats(text):
    return [float(x) for x in re.findall(FLOAT, text)]


def _tag(text, name):
    """An integer or float parameter from the OUTCAR parameter section."""
    m = re.search(rf"^\s*{name}\s*=\s*({FLOAT})", text, re.MULTILINE)
    return float(m.group(1)) if m else None


def _selective_free(path):
    """Per-atom True if the atom may move in any direction, None without selective dynamics."""
    try:
        lines = Path(path).read_text(errors="replace").splitlines()
    except OSError:
        return None
    start = 7 if not lines[5].split()[0].isdigit() else 6
    if lines[start].strip()[:1].upper() != "S":
        return None
    poscar = read_poscar(path)
    if poscar is None:
        return None
    free = []
    for line in lines[start + 2 : start + 2 + poscar.nions]:
        flags = line.split()[3:6]
        free.append(any(f.upper().startswith("T") for f in flags))
    return free


def _box_warnings(text):
    """The text of each boxed warning or advice block in an OUTCAR, without duplicates."""
    found, block = [], []
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("|") and stripped.endswith("|"):
            inner = stripped[1:-1].strip()
            if inner and not BANNER.match(inner):
                block.append(inner)
        elif block:
            message = " ".join(" ".join(block).split())
            message = re.sub(r"^ADVICE TO THIS USER RUNNING .*?VOICE \.\.\.\):\s*", "", message)
            if message not in found:
                found.append(message)
            block = []
    return found


def read_outcar(path):
    text = Path(path).read_text(errors="replace")
    out = {"file": Path(path).name}
    first = text.split("\n", 1)[0].split()
    out["vasp_version"] = first[0] if first and first[0].startswith("vasp") else None
    out["finished"] = "General timing and accounting informations for this job" in text
    m = re.search(rf"Elapsed time \(sec\):\s*({FLOAT})", text)
    out["elapsed_s"] = float(m.group(1)) if m else None

    params = {}
    for name in (
        "NIONS",
        "ISPIN",
        "NELM",
        "EDIFF",
        "EDIFFG",
        "NSW",
        "IBRION",
        "ISIF",
        "POTIM",
        "NBANDS",
    ):
        params[name] = _tag(text, name)
    m = re.search(r"number of ions\s+NIONS\s*=\s*(\d+)", text)
    if m:
        params["NIONS"] = float(m.group(1))
    m = re.search(r"NBANDS\s*=\s*(\d+)", text)
    if m:
        params["NBANDS"] = float(m.group(1))
    m = re.search(r"LNONCOLLINEAR\s*=\s*([TF])", text)
    params["LNONCOLLINEAR"] = m.group(1) == "T" if m else None
    m = re.search(r"ML_LMLFF\s*=\s*([TF])", text)
    params["ML_LMLFF"] = m.group(1) == "T" if m else None
    out["parameters"] = {
        k: (int(v) if isinstance(v, float) and v.is_integer() and k != "EDIFF" else v) for k, v in params.items()
    }

    # electronic steps of each ionic step: "Iteration  ionic( electronic)"
    electronic = {}
    for ionic, elec in re.findall(r"Iteration\s+(\d+)\(\s*(\d+)\)", text):
        electronic[int(ionic)] = int(elec)
    out["electronic_steps"] = dict(sorted(electronic.items()))  # ionic step -> electronic steps

    # one block per completed DFT ionic step; machine-learned steps print "ML energy"
    dft = re.findall(
        rf"^\s+energy\s+without entropy=\s*({FLOAT})\s+energy\(sigma->0\)\s*=\s*({FLOAT})",
        text,
        re.MULTILINE,
    )
    ml = re.findall(
        rf"ML energy\s+without entropy=\s*({FLOAT})\s+ML energy\(sigma->0\)\s*=\s*({FLOAT})",
        text,
    )
    toten = re.findall(rf"free\s+energy\s+TOTEN\s*=\s*({FLOAT})", text)
    out["ionic_steps"] = len(dft)
    out["ml_steps"] = len(ml)
    if dft:
        out["energy"] = {
            "free_energy_TOTEN": float(toten[-1]) if toten else None,
            "energy_without_entropy": float(dft[-1][0]),
            "energy_sigma_to_0": float(dft[-1][1]),
            "unit": "eV",
        }
    out["reached_required_accuracy"] = "reached required accuracy" in text

    # forces of the last DFT ionic step
    blocks = re.findall(r"TOTAL-FORCE \(eV/Angst\)\s*\n\s*-+\n(.*?)\n\s*-+", text, re.DOTALL)
    if blocks:
        out["forces"] = [_floats(line)[3:6] for line in blocks[-1].strip().splitlines()]
    pressures = re.findall(rf"external pressure =\s*({FLOAT}) kB", text)
    out["external_pressure_kB"] = float(pressures[-1]) if pressures else None
    fermi = re.findall(rf"E-fermi :\s*({FLOAT})", text)
    out["efermi_eV"] = float(fermi[-1]) if fermi else None
    mags = re.findall(r"number of electron[ \t]+\S+[ \t]+magnetization[ \t]*([^\n]*)", text)
    if mags:
        values = _floats(mags[-1])
        if values:
            out["magnetization_muB"] = values[0] if len(values) == 1 else values
    temps = [float(t) for t in re.findall(rf"\(temperature\s+({FLOAT})\s*K\)", text)]
    if temps:
        out["temperature_K"] = {
            "first": temps[0],
            "last": temps[-1],
            "mean": sum(temps) / len(temps),
            "steps": len(temps),
        }
    out["warnings"] = _box_warnings(text) + sorted(
        {" ".join(line.split()) for line in re.findall(r"^\s*WARNING[^\n]*", text, re.MULTILINE)}
    )
    return out


def read_oszicar(path):
    """Energies of each ionic step when there is no usable OUTCAR."""
    steps, electronic, count = [], [], 0
    for line in Path(path).read_text(errors="replace").splitlines():
        if re.match(r"^\s*(DAV|RMM|CG|SDA|ALG|DMP|ORT|EDD):", line):
            count += 1
        m = re.search(rf"F=\s*({FLOAT})\s+E0=\s*({FLOAT})", line)
        if m:
            steps.append({"F": float(m.group(1)), "E0": float(m.group(2))})
            electronic.append(count)
            count = 0
    return {
        "file": Path(path).name,
        "ionic_steps": len(steps),
        "electronic_steps": electronic,
        "last": steps[-1] if steps else None,
    }


def _band_edges(states):
    """VBM and CBM of (energy, occupation, k index) states; occupied means occupation > 0.5."""
    occupied = [s for s in states if s[1] > 0.5]
    empty = [s for s in states if s[1] <= 0.5]
    if not occupied or not empty:
        return None
    vbm, cbm = max(occupied), min(empty)
    partial = sum(1 for s in states if PARTIAL[0] < s[1] < PARTIAL[1])
    return vbm, cbm, partial


def _edges_dict(edges, kpoints):
    vbm, cbm, partial = edges
    out = {
        "vbm_eV": round(vbm[0], 4),
        "cbm_eV": round(cbm[0], 4),
        "vbm_kpoint": kpoints[vbm[2]] if vbm[2] < len(kpoints) else vbm[2],
        "cbm_kpoint": kpoints[cbm[2]] if cbm[2] < len(kpoints) else cbm[2],
        "partially_occupied_states": partial,
    }
    if partial:
        out["gap_eV"] = 0.0
        out["metallic"] = True
    else:
        out["gap_eV"] = round(max(0.0, cbm[0] - vbm[0]), 4)
        out["metallic"] = out["gap_eV"] == 0.0
        out["direct"] = vbm[2] == cbm[2] if out["gap_eV"] > 0 else None
    return out


def read_vasprun_gap(path):
    """Band gap, per spin and overall, from the last eigenvalues and occupations in vasprun.xml.

    A state counts as occupied when its occupation is above 0.5. States with an
    occupation between 0.01 and 0.99 lie at the Fermi level: the channel is then
    reported as metallic instead of giving the spurious gap between two partly
    filled states. Returns None if the file is incomplete (running or crashed).
    """
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError:
        return None
    calcs = root.findall("calculation")
    eig = calcs[-1].find("eigenvalues") if calcs else None
    if eig is None:
        return None
    kpoints = [_floats(v.text) for v in root.findall("kpoints/varray[@name='kpointlist']/v")]
    spins = {}
    for spin_set in eig.findall("array/set/set"):
        states = []
        for k, k_set in enumerate(spin_set.findall("set")):
            for r in k_set.findall("r"):
                energy, occ = (float(x) for x in r.text.split()[:2])
                states.append((energy, occ, k))
        spins[spin_set.get("comment")] = states
    overall = _band_edges([s for states in spins.values() for s in states])
    if overall is None:
        return None
    result = _edges_dict(overall, kpoints)
    if len(spins) > 1:
        result["spins"] = {
            name: _edges_dict(e, kpoints) for name, states in spins.items() if (e := _band_edges(states))
        }
    efermi = root.find("calculation/dos/i[@name='efermi']")
    result["efermi_eV"] = float(efermi.text) if efermi is not None else None
    result["occupied_if"] = "occupation > 0.5"
    if result["metallic"] and result["partially_occupied_states"]:
        result["note"] = (
            "Partly occupied states at the Fermi level: the system is metallic in this calculation, "
            "or SIGMA is larger than the gap."
        )
    return result


def scan_stdout(directory):
    """Known VASP error messages in the job's standard output files."""
    directory = Path(directory)
    files = [directory / n for n in STDOUT_NAMES if (directory / n).is_file()]
    for pattern in STDOUT_GLOBS:
        files += [p for p in sorted(directory.glob(pattern)) if p not in files and p.name != "OUTCAR"]
    hits = []
    for path in files:
        try:
            text = path.read_text(errors="replace")
        except OSError:
            continue
        for message, meaning in STDOUT_ERRORS.items():
            count = text.count(message)
            if count:
                hits.append(
                    {
                        "file": path.name,
                        "message": message,
                        "meaning": meaning,
                        "count": count,
                    }
                )
    return hits


def calculation_kind(ibrion, nsw, images=0):
    if images:
        return "nudged elastic band"
    if nsw in (None, 0) or ibrion == -1:
        return "single point"
    if ibrion == 0:
        return "molecular dynamics"
    if ibrion in (5, 6):
        return "phonons (finite differences)"
    if ibrion in (7, 8):
        return "phonons (DFPT)"
    if ibrion in (1, 2, 3):
        return "relaxation"
    if ibrion == 44:
        return "transition state (improved dimer)"
    return f"IBRION = {ibrion}"


def summarize(directory):
    """What happened in one calculation directory, with a list of problems found."""
    directory = Path(directory)
    summary = {"path": str(directory), "problems": []}
    problems = summary["problems"]
    incar = read_incar(directory / "INCAR") if (directory / "INCAR").is_file() else None

    outcar = directory / "OUTCAR"
    if outcar.is_file() and outcar.stat().st_size > 0:
        o = read_outcar(outcar)
    elif (directory / "OSZICAR").is_file():
        o = None
        summary["oszicar"] = read_oszicar(directory / "OSZICAR")
        problems.append("No OUTCAR: only OSZICAR energies are available.")
    else:
        summary["status"] = "not run" if incar else "no VASP files"
        summary["stdout_errors"] = scan_stdout(directory)
        return summary

    if o is not None:
        p = o["parameters"]
        images = incar.int("IMAGES", 0) if incar else 0
        kind = calculation_kind(p.get("IBRION"), p.get("NSW"), images)
        if p.get("ML_LMLFF"):
            kind += " with machine-learned force field"
        summary.update(
            {
                "kind": kind,
                "vasp_version": o["vasp_version"],
                "finished": o["finished"],
                "elapsed_s": o["elapsed_s"],
                "nions": p.get("NIONS"),
                "ionic_steps": o["ionic_steps"],
            }
        )
        if o["ml_steps"]:
            summary["ml_steps"] = o["ml_steps"]
        if "energy" in o:
            summary["energy"] = o["energy"]
            if o["ml_steps"]:
                summary["energy"]["note"] = "from the last DFT step; machine-learned steps are not included"
            if p.get("NIONS"):
                summary["energy"]["energy_sigma_to_0_per_atom"] = o["energy"]["energy_sigma_to_0"] / p["NIONS"]
        for key in (
            "external_pressure_kB",
            "efermi_eV",
            "magnetization_muB",
            "temperature_K",
        ):
            if o.get(key) is not None:
                summary[key] = o[key]

        if not o["finished"]:
            problems.append(
                "OUTCAR has no timing section: the job is still running, crashed, or was killed "
                "(wall time or memory). Check the job's standard output."
            )

        nelm, steps = p.get("NELM"), o["electronic_steps"]
        if nelm and steps:
            hit = [i for i, n in steps.items() if n >= nelm]
            summary["electronic"] = {
                "NELM": nelm,
                "max_steps": max(steps.values()),
                "ionic_steps_at_NELM": hit,
            }
            if hit:
                where = "the only ionic step" if len(steps) == 1 else f"ionic step(s) {_ranges(hit)}"
                problems.append(
                    f"Electronic loop reached NELM = {nelm} in {where}: that SCF is not converged to EDIFF "
                    f"and its energy and forces are unreliable."
                )

        if kind.startswith("relaxation") and o["ionic_steps"]:
            converged = o["reached_required_accuracy"]
            summary["ionic_converged"] = converged
            nsw = p.get("NSW")
            if not converged and o["finished"]:
                problems.append(
                    f"Relaxation stopped after {o['ionic_steps']} ionic steps (NSW = {nsw}) without reaching "
                    f"EDIFFG = {p.get('EDIFFG')}. Restart from CONTCAR."
                )

        forces = o.get("forces") if not kind.startswith("molecular dynamics") else None
        if forces:
            free = _selective_free(directory / "CONTCAR") or _selective_free(directory / "POSCAR")
            norms = [math.sqrt(sum(c * c for c in f)) for f in forces]
            if free and len(free) == len(norms):
                norms = [n for n, f in zip(norms, free) if f]
                summary["selective_dynamics_free_atoms"] = len(norms)
            if norms:
                summary["max_force_eV_per_A"] = round(max(norms), 5)
            ediffg = p.get("EDIFFG")
            if norms and ediffg is not None and ediffg < 0 and max(norms) > -ediffg and kind.startswith("relaxation"):
                problems.append(
                    f"Largest force on a free atom is {max(norms):.3f} eV/A, above |EDIFFG| = {-ediffg} eV/A."
                )

        if o["warnings"]:
            summary["outcar_warnings"] = o["warnings"]

    vasprun = directory / "vasprun.xml"
    if vasprun.is_file():
        with open(vasprun, "rb") as f:
            f.seek(max(0, vasprun.stat().st_size - 200))
            complete = b"</modeling>" in f.read()
        gap = read_vasprun_gap(vasprun) if complete else None
        if gap is not None:
            summary["band_gap"] = gap
        elif not complete and o is not None and o["finished"]:
            problems.append("vasprun.xml is incomplete although OUTCAR finished (disk full or killed while writing).")

    errors = scan_stdout(directory)
    if errors:
        summary["stdout_errors"] = errors
        problems.append("Error messages in standard output: " + "; ".join(sorted({e["meaning"] for e in errors})) + ".")

    if o is not None:
        summary["status"] = (
            "finished"
            if o["finished"] and not problems
            else ("finished with problems" if o["finished"] else "incomplete")
        )
    return summary


def _ranges(numbers):
    """[1, 2, 3, 7] -> '1-3, 7'"""
    parts, start = [], numbers[0]
    for a, b in zip(numbers, numbers[1:] + [None]):
        if b != a + 1:
            parts.append(f"{start}-{a}" if a != start else str(a))
            start = b
    return ", ".join(parts)
