"""Fixtures copy the line formats of real VASP 5.4.4 and 6.5.1 OUTCARs."""

import asyncio
import json

import pytest
from mcp import Client

from vasp_mcp import server as srv
from vasp_mcp.outputs import read_outcar, read_vasprun_gap, summarize

HEADER = """ vasp.6.4.2 20Jul23 (build Sep 01 2023 10:00:00) complex
   number of dos      NEDOS =    301   number of ions     NIONS =      2
   ISPIN  =      2    spin polarized calculation?
   LNONCOLLINEAR =      F non collinear calculations
   NELM   =      5;   NELMIN=  2; NELMDL= -5     # of ELM steps
   EDIFF  = 0.1E-05   stopping-criterion for ELM
   EDIFFG = -.1E-01   stopping-criterion for IOM
   NSW    =     50    number of steps for IOM
   IBRION =      2    ionic relax: 0-MD 1-quasi-New 2-CG
"""
BOX = """ -----------------------------------------------------------------------------
|                                                                             |
|           W    W    AA    RRRRR   N    N  II  N    N   GGGG   !!!           |
|                                                                             |
|      Your generating k-point grid is not commensurate to the symmetry      |
|      of the lattice.                                                        |
|                                                                             |
 -----------------------------------------------------------------------------
"""
TIMING = """ General timing and accounting informations for this job:
                         Elapsed time (sec):      123.456
"""


def ionic_step(n, electronic, energy, fz):
    lines = [f"--------------------- Iteration {n:6d}({e:4d})  ---------------------" for e in range(1, electronic + 1)]
    lines += [
        "  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)",
        f"  free  energy   TOTEN  =       {energy - 0.001:.8f} eV",
        f"  energy  without entropy=      {energy + 0.001:.8f}  energy(sigma->0) =      {energy:.8f}",
        " POSITION                                       TOTAL-FORCE (eV/Angst)",
        " -----------------------------------------------------------------------------------",
        f"      0.00000      0.00000      0.00000         0.000000      0.000000     {fz:.6f}",
        f"      0.00000      0.00000      2.00000         0.000000      0.000000     {-fz:.6f}",
        " -----------------------------------------------------------------------------------",
        "  external pressure =       -0.38 kB  Pullay stress =        0.00 kB",
        " number of electron      12.0000000 magnetization       2.0000000",
        " augmentation part        0.9218444 magnetization       0.6538737",
        " E-fermi :  -1.2345     XC(G=0):  -4.3872     alpha+bet : -4.2350",
    ]
    return "\n".join(lines) + "\n"


def write(folder, outcar, **files):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "OUTCAR").write_text(outcar)
    for name, text in files.items():
        (folder / name).write_text(text)
    return folder


def test_converged_relaxation(tmp_path):
    outcar = HEADER + BOX + ionic_step(1, 3, -10.0, 0.5) + ionic_step(2, 2, -10.2, 0.004)
    outcar += " reached required accuracy - stopping structural energy minimisation\n" + TIMING
    s = summarize(write(tmp_path, outcar))
    assert s["status"] == "finished" and s["problems"] == []
    assert s["kind"] == "relaxation" and s["ionic_converged"] is True
    assert s["ionic_steps"] == 2 and s["energy"]["energy_sigma_to_0"] == -10.2
    assert s["energy"]["energy_sigma_to_0_per_atom"] == -5.1
    assert s["max_force_eV_per_A"] == 0.004
    assert s["magnetization_muB"] == 2.0 and s["efermi_eV"] == -1.2345
    assert s["outcar_warnings"] == ["Your generating k-point grid is not commensurate to the symmetry of the lattice."]


def test_unconverged_relaxation(tmp_path):
    outcar = HEADER + ionic_step(1, 5, -10.0, 0.5) + ionic_step(2, 3, -10.1, 0.2) + TIMING
    s = summarize(write(tmp_path, outcar))
    assert s["status"] == "finished with problems"
    assert s["electronic"]["ionic_steps_at_NELM"] == [1]
    text = " ".join(s["problems"])
    assert "reached NELM = 5 in ionic step(s) 1" in text
    assert "without reaching EDIFFG" in text
    assert "0.200 eV/A, above |EDIFFG| = 0.01" in text


def test_incomplete_run_and_stdout_error(tmp_path):
    outcar = HEADER + ionic_step(1, 2, -10.0, 0.5)
    s = summarize(write(tmp_path, outcar, **{"slurm-42.out": "ZBRENT: fatal error in bracketing\n"}))
    assert s["status"] == "incomplete" and s["finished"] is False
    assert s["stdout_errors"][0]["meaning"] == "line search in the ionic relaxation failed"


def test_selective_dynamics_ignores_fixed_atoms(tmp_path):
    poscar = "x\n1.0\n5 0 0\n0 5 0\n0 0 5\nH\n2\nSelective dynamics\nCartesian\n0 0 0 F F F\n0 0 2 T T T\n"
    outcar = HEADER + ionic_step(1, 2, -10.0, 0.5) + TIMING
    s = summarize(write(tmp_path, outcar, POSCAR=poscar))
    assert s["selective_dynamics_free_atoms"] == 1 and s["max_force_eV_per_A"] == 0.5


def test_ml_md_counts_dft_and_ml_steps(tmp_path):
    md = HEADER.replace("IBRION =      2", "IBRION =      0").replace("ISPIN  =      2", "ISPIN  =      1")
    md += "   ML_LMLFF  =      T\n"
    md += ionic_step(1, 2, -10.0, 0.5).replace("magnetization       2.0000000", "magnetization ")
    for _ in range(3):
        md += "  ML energy  without entropy=     -10.1  ML energy(sigma->0) =     -10.1\n"
        md += "  kin. lattice  EKIN_LAT=         0.000000  (temperature  300.00 K)\n"
    s = summarize(write(tmp_path, md + TIMING))
    assert s["kind"] == "molecular dynamics with machine-learned force field"
    assert s["ionic_steps"] == 1 and s["ml_steps"] == 3
    assert s["temperature_K"]["mean"] == 300.0
    assert "magnetization_muB" not in s and "max_force_eV_per_A" not in s


def vasprun(states):
    """Minimal vasprun.xml: states[spin][k] = [(energy, occupation), ...]."""
    kpts = "".join(f"<v>{k * 0.5} 0 0</v>" for k in range(len(states[0])))
    spins = ""
    for s, spin in enumerate(states, 1):
        ks = "".join(
            f'<set comment="kpoint {k}">' + "".join(f"<r>{e} {o}</r>" for e, o in band) + "</set>"
            for k, band in enumerate(spin, 1)
        )
        spins += f'<set comment="spin {s}">{ks}</set>'
    return (
        f'<modeling><kpoints><varray name="kpointlist">{kpts}</varray></kpoints><calculation>'
        f"<eigenvalues><array><set>{spins}</set></array></eigenvalues>"
        '<dos><i name="efermi"> 0.5 </i></dos></calculation></modeling>'
    )


def test_gap_indirect(tmp_path):
    (tmp_path / "vasprun.xml").write_text(vasprun([[[(0.0, 1), (1.4, 0)], [(0.2, 1), (1.5, 0)]]]))
    gap = read_vasprun_gap(tmp_path / "vasprun.xml")
    assert gap["gap_eV"] == 1.2 and gap["direct"] is False
    assert gap["vbm_kpoint"] == [0.5, 0.0, 0.0] and gap["cbm_kpoint"] == [0.0, 0.0, 0.0]


def test_gap_half_metal(tmp_path):
    up = [[(0.0, 1), (0.45, 0.6), (0.6, 0.3)]]
    down = [[(-1.0, 1), (2.0, 0)]]
    (tmp_path / "vasprun.xml").write_text(vasprun([up, down]))
    gap = read_vasprun_gap(tmp_path / "vasprun.xml")
    assert gap["metallic"] is True and gap["gap_eV"] == 0.0
    assert gap["spins"]["spin 1"]["metallic"] is True
    assert gap["spins"]["spin 2"]["gap_eV"] == 3.0


def test_truncated_vasprun_is_ignored(tmp_path):
    (tmp_path / "vasprun.xml").write_text("<modeling><calculation>")
    assert read_vasprun_gap(tmp_path / "vasprun.xml") is None


def test_outcar_vasp5_format(tmp_path):
    text = HEADER.replace(" vasp.6.4.2", " vasp.5.4.4.18Apr17") + ionic_step(1, 2, -1.0, 0.0)
    assert read_outcar(write(tmp_path, text) / "OUTCAR")["vasp_version"] == "vasp.5.4.4.18Apr17"


def call(name, args):
    async def go():
        async with Client(srv.server) as client:
            result = await client.call_tool(name, args)
            return result.is_error, result.content[0].text

    return asyncio.run(go())


def test_server_tools(tmp_path):
    calc = tmp_path / "relax"
    done = " reached required accuracy - stopping structural energy minimisation\n" + TIMING
    write(calc, HEADER + ionic_step(1, 2, -1.0, 0.0) + done, INCAR="ENCUT = 400\nISMEAR = 0\nLREAL = .TRUE.\n")
    error, text = call("check_vasp_inputs", {"path": str(tmp_path)})
    assert not error and {f["rule"] for f in json.loads(text)["calculations"][0]["findings"]} == {"VL008"}
    error, text = call("summarize_vasp_run", {"path": str(calc)})
    assert not error and json.loads(text)["status"] == "finished"
    error, text = call("scan_vasp_project", {"path": str(tmp_path)})
    assert not error and json.loads(text)["counts"] == {"finished": 1}
    error, text = call("read_vasp_inputs", {"path": str(calc)})
    assert not error and json.loads(text)["incar"]["ENCUT"] == "400"
    error, text = call("list_vasp_lint_rules", {})
    assert not error and len(json.loads(text)["rules"]) >= 30


def test_server_root_restriction(tmp_path, monkeypatch):
    monkeypatch.setattr(srv, "ROOTS", [tmp_path / "allowed"])
    error, text = call("summarize_vasp_run", {"path": str(tmp_path)})
    assert error and "outside the folders" in text


@pytest.mark.parametrize("name", ["summarize_vasp_run", "check_vasp_inputs"])
def test_server_missing_folder(name, tmp_path):
    error, text = call(name, {"path": str(tmp_path / "nope")})
    assert error and "not a folder" in text
