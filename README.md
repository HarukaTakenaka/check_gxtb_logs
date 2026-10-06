# check_gxtb_logs
A command-line tool that checks a folder of Gaussian 16 jobs run with xTB as the external energy engine (via xtb-gaussian, e.g. g-xTB or GFN2-xTB) and tells you which ones finished correctly, which failed and why, and how the successful structures compare in energy. It can also write restart inputs and tidy failed jobs into a separate folder.

It was written for screening transition-metal complexes (e.g. Ni/Pd oxidative addition products, isomers, and spin states), but works for any Gaussian + xtb-gaussian optimization, frequency, TS, scan, or IRC job.

Requirements
Python 3 with numpy (nothing else)
Gaussian 16 output files (.log). Tested with Gaussian 16 C.01, xtb-gaussian, and xtb 6.7.1 (--gxtb). Ordinary Gaussian DFT logs should mostly work too (energies are read from SCF Done when no external energy is present), but this has not been tested.
Installation
bash
git clone <this-repo-url>
cp check_gxtb_logs.py ~/bin/          # or anywhere on your PATH
python ~/bin/check_gxtb_logs.py -h    # prints usage if it works

If the file passed through Windows, convert line endings first: dos2unix check_gxtb_logs.py.

Quick start

Run it in the folder that contains the .log (and matching .com) files:

bash
python check_gxtb_logs.py

Example output:

name                          status          job  chg mult steps imag  G_hartree     dG    notes
terpyNiBr_OA_PhI_o5_m2        OK              opt  0   2    37    -     -5353.945896  0.0   C2-X3 2.824
terpyNiBr_OA_PhI_o3_m2        OK              opt  0   2    26    -     -5353.940363  3.5   C2-X3 4.459
terpyNiBr_OA_PhI_o5_m4        OK              opt  0   4    57    -     -5353.915752  18.9  C2-X3 3.558
terpyNiBr_OA_PhI_o1_m2        XTB_FAILED      opt  0   2    5     -     -             -
TS_1_opt                      NOT_CONVERGED   ts   0   2    200   -     -             -     metal: lost C2-Ni7, Ni7-Cl22; ...

A full table is also written to check_summary.csv.

What it checks

For every .log file in the folder (SLURM slurm-*.log files are ignored):

Termination: number of Normal termination lines, Error termination, and xTB crashes (abnormal termination of xtb, the FIO-F-217 read error, zero energies).
Optimization: number of steps, Stationary point found, and the last four convergence criteria.
Frequencies (last frequency job): number and values of imaginary frequencies, lowest frequency.
Energies: final electronic energy, ZPE, enthalpy, and Gibbs free energy (Hartree).
Charge and spin: Gaussian's charge/multiplicity vs. the --charge/--uhf passed to xTB.
Bonding changes: compares the first and last geometries and reports metal–ligand bonds lost or formed and other bonds broken or formed (distance-based, using covalent radii).
Job type: opt, ts, scan, irc, freq, or sp, detected from the route section.
Consistency: warns if logs in the same folder used different xTB settings (energies from different settings are not comparable).
Status codes
Status	Meaning
OK	Normal termination, converged, and the expected number of imaginary frequencies (0 for minima, 1 for TS)
OK_NO_FREQ	Converged, but no frequency calculation was found
OK+CHECK_BONDS	Finished, but bonding changed during the run; look at the structure
NOT_CONVERGED	Hit the step limit or stopped before converging
IMAGINARY(n)	A minimum optimization that ended with n imaginary frequencies
TS_NIMAG=n	A TS optimization with other than one imaginary frequency
XTB_FAILED	xTB crashed inside the Gaussian job
ERROR	Gaussian error termination for another reason
INCOMPLETE	No termination message: still running, or killed (time limit, preemption)
CHARGE/SPIN_MISMATCH	xTB's --charge/--uhf disagree with Gaussian's charge/multiplicity
Relative energies of isomers and spin states

Files named with the convention <species>_o<isomer>_m<multiplicity> (for example terpyNiBr_OA_PhI_o5_m2) are grouped by <species>. Within each group, dG is the free energy of each OK structure relative to the lowest one, in kcal/mol (1 Hartree = 627.5095 kcal/mol). This lets you read off the most stable isomer and spin state directly. If the multiplicity in the file name differs from the one in the log, the script says so in the notes.

Files that don't follow the convention get their own one-member group (dG = 0.0).

Two additional notes are printed when they apply:

C2-X3: distance between atoms 2 and 3 when atom 2 is a carbon and atom 3 is a halogen or O bound to the metal (the numbering used for oxidative addition products: metal = 1, ipso C = 2, leaving-group atom = 3). Roughly 3.5 Å or less means X is cis to the aryl; 4 Å or more means trans.
IRC minima: how many directions of an IRC reached a minimum.
Options
python check_gxtb_logs.py [folder] [--csv FILE] [--restart] [--displace A]
                          [--move-failed] [--fail-statuses LIST] [--dry-run]
Option	What it does
folder	Folder to check (default: current folder)
--csv FILE	Name of the summary CSV (default check_summary.csv)
--restart	Write new inputs for NOT_CONVERGED and INCOMPLETE jobs, starting from their last geometry
--displace A	For minima with imaginary frequencies, write new inputs displaced along the lowest mode (largest atomic shift A Å, e.g. 0.15)
--move-failed	Move all files of failed jobs (name.log, name.com, name.chk, name.xyz, ...) to ./failed/
--fail-statuses LIST	Comma-separated statuses counted as failed (default XTB_FAILED,ERROR,NOT_CONVERGED,CHARGE/SPIN_MISMATCH)
--dry-run	With --move-failed, only show what would be moved
Restarting and displacing
bash
python check_gxtb_logs.py --restart --displace 0.15

New inputs are written to ./resubmit/ together with filenames_resubmit.txt, and the script prints the sbatch command to submit them. Each new input copies the original .com file (route, AddRedun/ModRedundant lines, and --Link1-- steps) and replaces only the coordinates and the %chk names. New names get _r1, _r2, ... (restart) or _d1, _d2, ... (displaced).

Safety checks:

A TS optimization whose metal–ligand bonding changed is not restarted: it has left the saddle region, and restarting it will not help. Build a better guess (e.g. a relaxed scan) instead.
If the steps of an input pass different --charge or --uhf values to xTB (for example a --Link1-- frequency step with the wrong charge), a warning is printed.
Displacements use the normal-mode vectors in the frame Gaussian printed them in (standard orientation, or input orientation with nosymm).

The original .com file must be in the same folder as its .log.

Moving failed jobs
bash
python check_gxtb_logs.py --move-failed --dry-run     # preview
python check_gxtb_logs.py --move-failed               # move

Files are moved to ./failed/, and the reason for each job is appended to failed/failed_reasons.txt. Nothing is overwritten: if a file of the same name is already in failed/, the new one gets an .old1, .old2, ... suffix.

INCOMPLETE jobs are not moved by default, because they may still be running. Check squeue -u $USER before adding INCOMPLETE to --fail-statuses.
IMAGINARY jobs are not moved by default (they usually only need --displace). Add them with --fail-statuses XTB_FAILED,ERROR,NOT_CONVERGED,IMAGINARY.
When combined with --restart, restart inputs are written before anything is moved.
SLURM .out/.err files are not moved, since their names don't contain the structure name.
Typical workflow
bash
cd my_results_folder
python ~/bin/check_gxtb_logs.py                          # 1. look at the table
python ~/bin/check_gxtb_logs.py --restart --displace 0.15 --move-failed
cd resubmit
sbatch --array=0-<N-1>%20 ../g16_xtb.slurm filenames_resubmit.txt   # 2. resubmit

Then open anything marked +CHECK_BONDS (or any surprising result) in GaussView before using it.

Limitations
Bonding changes are detected from distances (1.15 × sum of covalent radii for non-metal pairs, 1.3 × for metal–ligand pairs). Unusually long but genuine bonds (e.g. weakly bound triflate) or short non-bonded contacts can be misreported. Always confirm flagged structures by eye.
Only the first metal atom in a structure is used for the coordination-number and C2-X3 checks.
XTB_FAILED means xTB crashed, but xTB's own error message is not in the Gaussian log. To see it, run xTB directly on the structure, e.g. xtb mol.xyz --gxtb --chrg 0 --uhf 1 --sp.

Acknowledgement: Claude
Energies are compared only within a name group. Comparing across species (e.g. a reaction energy) is up to you, and only meaningful when all structures used identical xTB settings.
The dG values are xTB-level free energies, useful for screening, not for reporting.
