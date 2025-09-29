#!/bin/bash
#SBATCH -p gpu-vladimir
#SBATCH --gres=gpu:1
#SBATCH -t 2-12:00:00
#SBATCH --job-name=seqDesignA_kpt_allcombos
#SBATCH --mem=50G
#SBATCH --mail-user=ahgonzalez@ucdavis.edu
#SBATCH --mail-type=END

# =============================================================
#  A-chain k-point mutational design over ALL site combinations
#  - Iterates all nCk position-combos from a provided A-chain site list
#  - Warns if combos exceed a threshold (default 1000)
#  - Accepts sequences only if ALL k sites differ from WT (A-chain only)
#  - Deduplicates globally (SHA1 on A-chain sequence)
#  - Excludes cysteine via LigandMPNN's --omit_AA "C"
#  - Keeps chain B fully fixed via --chains_to_design A
# =============================================================

# ----------------------------- Usage -----------------------------
usage() {
  cat <<USAGE
Usage: $0 --pdb <path> --outdir <path> --positions "A14 A15 ..." --num_mutation_pos <K> [options]

Required:
  --pdb <path>                 Path to input PDB (chain A mutates, chain B fixed)
  --outdir <path>              Output directory
  --positions "A14 A15 ..."   Space/comma separated A-chain sites (or use --positions-file)
  --num_mutation_pos <K>       Number of positions to mutate simultaneously (K >= 2)

Optional:
  --positions-file <file>      File with A-chain positions (one per line)
  --seqs-per-combo <int>       Target sequences per combination (default 400)
  --batch-size <int>           LigandMPNN batch size (default 100)
  --temperature <float>        Sampling temperature (default 0.15)
  --model-type <name>          ligand_mpnn | soluble_mpnn | protein_mpnn (default soluble_mpnn)
  --checkpoint <path>          Explicit model .pt path (overrides defaults)
  --seed <int>                 Base RNG seed (default 57)
  --warn-combos-threshold <N>  Print BIG warning if total combos > N (default 1000)
  --abort-if-combos-over <N>   Abort run if total combos > N (no default)
  --max-accepted <int>         Stop early once this many unique accepted sequences gathered (0 = no cap; default 0)
  -h | --help                  Show this help and exit
USAGE
}

# -------------------------- Parse args ---------------------------
PDB=""; OUTDIR=""; POS_RAW=""; POSITIONS_FILE=""; NUM_MUT=0
SEQS_PER_COMBO=400; BATCH_SIZE=100; TEMP=0.15
MODEL_TYPE="soluble_mpnn"; CHECKPOINT=""; BASE_SEED=57
WARN_COMBOS_THRESH=1000; ABORT_IF_COMBOS_OVER=""; MAX_ACCEPTED=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pdb) PDB="$2"; shift ;;
    --outdir) OUTDIR="$2"; shift ;;
    --positions) POS_RAW="$2"; shift ;;
    --positions-file) POSITIONS_FILE="$2"; shift ;;
    --num_mutation_pos) NUM_MUT="$2"; shift ;;
    --seqs-per-combo) SEQS_PER_COMBO="$2"; shift ;;
    --batch-size) BATCH_SIZE="$2"; shift ;;
    --temperature) TEMP="$2"; shift ;;
    --model-type) MODEL_TYPE="$2"; shift ;;
    --checkpoint) CHECKPOINT="$2"; shift ;;
    --seed) BASE_SEED="$2"; shift ;;
    --warn-combos-threshold) WARN_COMBOS_THRESH="$2"; shift ;;
    --abort-if-combos-over) ABORT_IF_COMBOS_OVER="$2"; shift ;;
    --max-accepted) MAX_ACCEPTED="$2"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown parameter: $1" >&2; usage; exit 2 ;;
  esac
  shift
done

[[ -z "$PDB" || -z "$OUTDIR" ]] && { echo "Error: --pdb and --outdir are required." >&2; usage; exit 2; }
[[ $NUM_MUT -lt 2 ]] && { echo "Error: --num_mutation_pos K must be >= 2" >&2; exit 2; }

# ------------------------ Input positions ------------------------
if [[ -n "$POSITIONS_FILE" ]]; then
  [[ -f "$POSITIONS_FILE" ]] || { echo "Positions file not found: $POSITIONS_FILE" >&2; exit 2; }
  POS_RAW="$(tr '\n' ' ' < "$POSITIONS_FILE")"
fi
IFS=', ' read -r -a POS_ARR <<< "$POS_RAW"

# normalize & dedupe (order-preserving), keep only A-chain numeric entries
declare -A seen_pos
POS=()
for tok in "${POS_ARR[@]}"; do
  [[ "$tok" =~ ^A[0-9]+$ ]] || continue
  if [[ -z "${seen_pos[$tok]}" ]]; then
    seen_pos[$tok]=1
    POS+=("$tok")
  fi
done

if ((${#POS[@]} < NUM_MUT)); then
  echo "Error: you provided ${#POS[@]} usable positions, which is < K=$NUM_MUT" >&2; exit 2
fi

mkdir -p "$OUTDIR"
LOG="$OUTDIR/run.log"

# -------------------- WT map for chain A only --------------------
WTMAP_JSON="$OUTDIR/wt_map.json"
python - "$PDB" "$WTMAP_JSON" <<'PY'
import sys,json
pdb_path, out_json = sys.argv[1], sys.argv[2]
three={"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
       "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
       "THR":"T","TRP":"W","TYR":"Y","VAL":"V","SEC":"U","PYL":"O"}
aseq=[]; resnums=[]; seen=set()
with open(pdb_path) as f:
  for line in f:
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
      continue
    if line[21] != 'A':
      continue
    resnum=line[22:26].strip()
    resname=line[17:20].strip()
    key=(resnum,resname)
    if key in seen: continue
    seen.add(key)
    aseq.append(three.get(resname,'X'))
    try:
      resnums.append(int(resnum))
    except ValueError:
      pass
resmap={r:i for i,r in enumerate(resnums)}
json.dump({"aseq":"".join(aseq),"resnums":resnums,"resnum_to_idx":resmap},open(out_json,"w"))
PY

# verify requested positions exist in chain A
python - "$WTMAP_JSON" ${POS[@]} <<'PY'
import sys,json
m=json.load(open(sys.argv[1]))
idx={int(k):v for k,v in m["resnum_to_idx"].items()}
missing=[p for p in sys.argv[2:] if int(p[1:]) not in idx]
if missing:
  print("Error: missing positions", missing)
  sys.exit(2)
print("All requested positions found.")
PY

# ---------------------- Combos math & warn -----------------------
COMBO_COUNT=$(python - <<PY
import math
from math import comb
n=${#POS[@]}
k=${NUM_MUT}
print(comb(n,k))
PY
)

EST_CANDIDATES=$(( COMBO_COUNT * SEQS_PER_COMBO ))
{
  echo "PDB: $PDB"
  echo "Outdir: $OUTDIR"
  echo "A-chain positions (N=${#POS[@]}): ${POS[*]}"
  echo "K (num_mutation_pos): $NUM_MUT"
  echo "Total combinations (nCk): $COMBO_COUNT"
  echo "Est. raw candidates to parse (combos * seqs-per-combo): $COMBO_COUNT * $SEQS_PER_COMBO = $EST_CANDIDATES"
} | tee "$LOG"

if [[ -n "$ABORT_IF_COMBOS_OVER" && "$COMBO_COUNT" -gt "$ABORT_IF_COMBOS_OVER" ]]; then
  echo "ABORT: total combos ($COMBO_COUNT) > --abort-if-combos-over=$ABORT_IF_COMBOS_OVER" | tee -a "$LOG"
  exit 3
fi

if [[ "$COMBO_COUNT" -gt "$WARN_COMBOS_THRESH" ]]; then
  cat <<WARN | tee -a "$LOG"
=====================================================================
  WARNING: total site-combinations = $COMBO_COUNT  (> $WARN_COMBOS_THRESH)
  LigandMPNN is fast, but iterating many combinations is the bottleneck.
  You indicated >1000 combos may take >1 day on your setup.
  Consider:
    • Lower K or reduce site list size
    • Increase --seqs-per-combo only if necessary
    • Use --abort-if-combos-over N to stop automatically
    • Or set a soft cap with --max-accepted <N> (e.g., 30000)
=====================================================================
WARN
fi

# ------------------------ Model checkpoints ----------------------
DEFAULT_CKPT_LIGAND="/share/yarovlab/ahgz/apps/LigandMPNN/model_params/ligandmpnn_v_48_020.pt"
DEFAULT_CKPT_SOLUBLE="/share/yarovlab/ahgz/apps/LigandMPNN/model_params/solublempnn_v_48_020.pt"
DEFAULT_CKPT_PROTEIN="/share/yarovlab/ahgz/apps/LigandMPNN/model_params/proteinmpnn_v_48_020.pt"

MODEL_FLAGS=( --model_type "$MODEL_TYPE" )
if [[ -z "$CHECKPOINT" ]]; then
  case "$MODEL_TYPE" in
    ligand_mpnn)  CHECKPOINT="$DEFAULT_CKPT_LIGAND" ;;
    soluble_mpnn) CHECKPOINT="$DEFAULT_CKPT_SOLUBLE" ;;
    protein_mpnn) CHECKPOINT="$DEFAULT_CKPT_PROTEIN" ;;
    *) echo "Error: unknown --model-type '$MODEL_TYPE'" >&2; exit 2 ;;
  esac
fi
[[ -f "$CHECKPOINT" ]] || { echo "Error: checkpoint not found: $CHECKPOINT" >&2; exit 2; }
case "$MODEL_TYPE" in
  ligand_mpnn)  MODEL_FLAGS+=( --checkpoint_ligand_mpnn "$CHECKPOINT" ) ;;
  soluble_mpnn) MODEL_FLAGS+=( --checkpoint_soluble_mpnn "$CHECKPOINT" ) ;;
  protein_mpnn) MODEL_FLAGS+=( --checkpoint_protein_mpnn "$CHECKPOINT" ) ;;
  *) echo "Error: unknown --model-type '$MODEL_TYPE'" >&2; exit 2 ;;
 esac

# --------------------------- Enumerate ---------------------------
COMBOS_TXT="$OUTDIR/combos.txt"
python - "$COMBOS_TXT" $NUM_MUT ${POS[@]} <<'PY'
import sys,itertools
out=sys.argv[1]; k=int(sys.argv[2]); pos=sys.argv[3:]
with open(out,'w') as f:
  for comb in itertools.combinations(pos,k):
    f.write(":".join(comb)+"\n")
PY
TOTAL_COMBOS=$(wc -l < "$COMBOS_TXT")
echo "Combo list saved: $COMBOS_TXT ($TOTAL_COMBOS lines)" | tee -a "$LOG"

# --------------------------- Collectors --------------------------
ACCEPTED_FASTA="$OUTDIR/accepted.fa"; > "$ACCEPTED_FASTA"
ACCEPTED_CSV="$OUTDIR/accepted.csv"; echo "seq_id,positions,muts,overall_confidence" > "$ACCEPTED_CSV"
SEEN_HASH="$OUTDIR/seen.sha1"; > "$SEEN_HASH"

accepted=0; combo_idx=0

# ---------------------------- Main loop --------------------------
while IFS= read -r combo; do
  ((combo_idx++))
  IFS=':' read -r -a sites <<< "$combo"
  subdir="$OUTDIR/run_${combo_idx}_$(echo "$combo" | tr ':' '_')"
  mkdir -p "$subdir"

  echo "[$(date +'%F %T')] Combo $combo_idx/$TOTAL_COMBOS :: ${sites[*]} | accepted=$accepted" | tee -a "$LOG"

  # Launch LigandMPNN for this combo
  python /share/yarovlab/ahgz/apps/LigandMPNN/run.py \
    "${MODEL_FLAGS[@]}" \
    --pdb_path "$PDB" \
    --out_folder "$subdir" \
    --redesigned_residues "${sites[*]}" \
    --chains_to_design "A" \
    --omit_AA "C" \
    --batch_size "$BATCH_SIZE" \
    --number_of_batches $(((SEQS_PER_COMBO + BATCH_SIZE - 1) / BATCH_SIZE)) \
    --temperature "$TEMP" \
    --seed $((BASE_SEED + combo_idx)) \
    --save_stats 1

  # Parse FASTAs: chain A only, require ALL k sites != WT; store unique
  python - "$WTMAP_JSON" "$subdir" "$ACCEPTED_FASTA" "$ACCEPTED_CSV" "$SEEN_HASH" "$combo" <<'PY'
import sys,os,glob,re,json,hashlib
wt_json, subdir, outfa, outcsv, seenf, combo = sys.argv[1:7]
WT=json.load(open(wt_json))
aseq_wt=WT["aseq"]
idx={int(k):v for k,v in WT["resnum_to_idx"].items()}
sites=combo.split(':')
Aidx=[ idx[int(s[1:])] for s in sites ]
wt_at=[ aseq_wt[i] for i in Aidx ]

seen=set(open(seenf).read().split()) if os.path.exists(seenf) else set()
sha1=lambda s: hashlib.sha1(s.encode()).hexdigest()

# collect all fasta-like files (including seqs/*.fa*)
paths=sorted([p for p in glob.glob(f"{subdir}/**/*.fa*", recursive=True) if os.path.isfile(p)])

def iter_fa(path):
  with open(path) as f:
    hdr=None; seq=[]
    for line in f:
      line=line.strip()
      if not line: continue
      if line.startswith('>'):
        if hdr: yield hdr, ''.join(seq)
        hdr=line[1:]; seq=[]
      else:
        seq.append(line)
    if hdr: yield hdr, ''.join(seq)

# helper to parse confidence field from header if present
conf_re=re.compile(r"overall_confidence=([0-9.]+)")
get_conf=lambda h: (conf_re.search(h).group(1) if conf_re.search(h) else "")

new=0
with open(outfa,'a') as fa, open(outcsv,'a') as csv, open(seenf,'a') as sf:
  for path in paths:
    for h,seq_full in iter_fa(path):
      # split A vs others by the first ':' if present
      aseq = seq_full.split(':',1)[0]
      if len(aseq) != len(aseq_wt):
        continue
      aas=[ aseq[i] for i in Aidx ]
      # require ALL sites mutated vs WT
      if any(a==w for a,w in zip(aas, wt_at)):
        continue
      key=sha1(aseq)
      if key in seen:
        continue
      seen.add(key); sf.write(key+'\n')
      muts=';'.join(f"{w}{sites[i][1:]}{aas[i]}" for i,w in enumerate(wt_at))
      conf=get_conf(h)
      fa.write(f">{h}\n{aseq}\n")
      csv.write(f"{h},{','.join(sites)},{muts},{conf}\n")
      new+=1
print(new)
PY

  # update accepted count from FASTA headers
  if [[ -s "$ACCEPTED_FASTA" ]]; then
    accepted=$(grep -c '^>' "$ACCEPTED_FASTA")
  else
    accepted=0
  fi
  echo "Accepted so far: $accepted" | tee -a "$LOG"

  if [[ "$MAX_ACCEPTED" -gt 0 && "$accepted" -ge "$MAX_ACCEPTED" ]]; then
    echo "Reached --max-accepted=$MAX_ACCEPTED; stopping early." | tee -a "$LOG"
    break
  fi

done < "$COMBOS_TXT"

# ---------------------------- Summary ----------------------------
echo "Done. Accepted $accepted unique sequences." | tee -a "$LOG"
echo "FASTA: $ACCEPTED_FASTA"
echo "CSV:   $ACCEPTED_CSV"