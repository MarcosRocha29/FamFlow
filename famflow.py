#!/usr/bin/env python3
"""
famflow.py -- fluxo generico de exploracao de familias de proteinas

Da semente ao modelo, passando por rede de similaridade, comunidades,
taxonomia, arquitetura de dominio e vizinhanca genomica.

NAO e um pipeline com um objetivo fixo. E um fluxo com CONTRATOS entre
estagios: cada estagio le um artefato definido, escreve um artefato
definido, e nao sabe nada sobre o que vem depois. Voce para onde quiser.

=============================================================================
OS QUATRO PRINCIPIOS
=============================================================================
1. TABELA CENTRAL
   Existe uma unica tabela indexada por sequencia (00_table/sequences.tsv)
   que ACUMULA colunas conforme os estagios rodam. Ela e o entregavel,
   nao um subproduto. Qualquer estagio pode ser pulado ou substituido por
   uma analise externa, desde que a coluna correspondente apareca nela.

2. TODO ESTAGIO E TERMINAL
   Cada estagio escreve, alem do artefato, um REPORT.md legivel. Parar na
   SSN nao e parar no meio -- e parar com um relatorio de rede na mao.

3. BIOLOGIA EM CONFIG, NAO EM CODIGO
   O vocabulario de contexto (quais genes vizinhos importam) vem de um
   JSON de marcadores, nao de regex hardcoded. SEM config, o 'profile'
   faz descoberta NAO-SUPERVISIONADA: compara a frequencia de termos de
   cada comunidade contra o fundo e reporta o que esta enriquecido.
   Sem vocabulario, ele ACHA o vocabulario.

4. ANCORA OPCIONAL
   Nao ha "a query" privilegiada. Se voce der uma ancora (--anchor), a
   comunidade dela e destacada nos relatorios. Se nao der, todas as
   comunidades sao cidadas iguais.

=============================================================================
ESTAGIOS
=============================================================================
  seed      Normaliza a entrada (fasta, multi-fasta, HMM pronto, chopping)
  search    hmmbuild (se preciso) + hmmsearch contra o banco
  inspect   Score x cobertura por hit, ANTES de qualquer corte [grafico]
  collect   Recupera os hits (envelope OU proteina inteira) -> inicia a tabela
  cluster   Reduz redundancia PRESERVANDO o mapa de membros
  matrix    all-vs-all (diamond/blast/mmseqs)
  ssn       Rede + comunidades + varredura de corte salva
  annotate  Taxonomia e arquitetura de dominio (best effort, opcional)
  context   Vizinhanca genomica (escopo configuravel)
  profile   Descreve CADA comunidade -> tabela comparativa   [muito util como fim]
  select    Filtro OPCIONAL guiado por config (nao por codigo)
  build     Alinha + hmmbuild por grupo
  scan      Modelos vs alvo, + auto-scan de QC (modelos x comunidades)
  table     Exporta/inspeciona a tabela central
  report    Monta Metodos + Resultados em Markdown a partir dos checkpoints
  status    Progresso de todos os estagios

=============================================================================
USO
=============================================================================
    F="python famflow.py"

    # minimo: ja tenho o fasta dos hits, quero so a rede e as comunidades
    $F collect --from-fasta meus_hits.fasta --unit calicina
    $F cluster
    $F matrix  --cpu 24
    $F ssn     --cutoff-steps 30
    $F profile                       # <- pode parar aqui

    # fluxo completo, do zero
    $F seed    --input sementes.fasta --unit calicina
    $F search  --workers 4 --cpus-per-worker 6 --evalue 10
    $F inspect --min-coverage 0.6    # olha score x cobertura ANTES de cortar
    $F collect --mode envelope --min-coverage 0.6
    $F cluster
    $F matrix
    $F ssn
    $F annotate --taxonomy --arch
    $F context  --scope all
    $F profile  --markers perfis/t2ss.json     # ou sem --markers: descoberta
    $F select   --markers perfis/t2ss.json --min-markers 1
    $F build
    $F scan     --target genomas.faa
    $F report   -o relatorio.md
"""

import argparse
import concurrent.futures
import glob
import json
import math
import os
import re
import subprocess
import sys
import time
import traceback
from collections import Counter, defaultdict
from shutil import which

# banco default -- ajuste para o seu servidor
NR = "/scratch/global/databases/fadb/nr/nr"

DIRS = {
    "table":     "00_table",
    "seed":      "01_seed",
    "search":    "02_search",
    "inspect":   "02_inspect",
    "hits":      "03_hits",
    "clustered": "04_clustered",
    "matrix":    "05_matrix",
    "ssn":       "06_ssn",
    "annot":     "07_annot",
    "context":   "08_context",
    "profile":   "09_profile",
    "curated":   "10_curated",
    "models":    "11_models",
    "scan":      "12_scan",
}

TABLE = os.path.join(DIRS["table"], "sequences.tsv")


# =============================================================================
# tabela central -- a espinha do fluxo
#
# Uma linha por sequencia, indexada por 'id'. Cada estagio ACRESCENTA
# colunas. Nenhum estagio reescreve a tabela inteira do zero: o merge e
# por 'id', preservando o que ja estava la.
#
# Consequencia pratica: voce pode pular um estagio e preencher a coluna
# dele por fora (um script seu, um resultado de outra ferramenta), que o
# resto do fluxo continua funcionando. O contrato e a coluna, nao o codigo.
# =============================================================================

def load_table(path=TABLE):
    import pandas as pd
    if not os.path.exists(path):
        return pd.DataFrame(columns=["id"])
    return pd.read_csv(path, sep="\t", low_memory=False)


def save_table(df, path=TABLE):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    return path


def update_table(new_df, key="id", path=TABLE, overwrite=True):
    """
    Funde colunas novas na tabela central, casando por 'id'.

    overwrite=True   colunas homonimas sao substituidas pelos valores novos
    overwrite=False  colunas ja existentes sao preservadas (nao sobrescreve)

    Linhas de 'id' que ainda nao estavam na tabela sao ACRESCENTADAS -- assim
    'collect' cria a tabela e os demais estagios so a enriquecem.
    """
    import pandas as pd

    old = load_table(path)
    new_df = new_df.copy()
    new_df[key] = new_df[key].astype(str)

    if len(old) == 0:
        return save_table(new_df, path), new_df

    old[key] = old[key].astype(str)
    shared = [c for c in new_df.columns if c != key and c in old.columns]
    if shared:
        if overwrite:
            old = old.drop(columns=shared)
        else:
            new_df = new_df.drop(columns=shared)

    merged = old.merge(new_df, on=key, how="outer")
    return save_table(merged, path), merged


def table_columns_report(df):
    """Quantas linhas tem valor em cada coluna -- o mapa de cobertura do fluxo."""
    lines = []
    n = len(df)
    for c in df.columns:
        filled = int(df[c].notna().sum())
        pct = filled / n if n else 0
        lines.append(f"  {c:24} {filled:>8}/{n} ({pct:5.1%})")
    return "\n".join(lines)


# =============================================================================
# relatorio por estagio -- e o que torna cada estagio um FIM possivel
# =============================================================================

class Report:
    """
    Acumula linhas e grava REPORT.md no diretorio do estagio.

    Existe porque o criterio de 'parar aqui' e ter algo legivel na mao.
    Um .graphml nao e um resultado; um .graphml + um resumo do que ha nele e.
    """

    def __init__(self, stage, out_dir, title=None):
        self.stage = stage
        self.out_dir = out_dir
        self.lines = []
        self.title = title or stage
        os.makedirs(out_dir, exist_ok=True)

    def h(self, text, level=2):
        self.lines.append(f"\n{'#' * level} {text}\n")

    def p(self, text):
        self.lines.append(str(text))

    def kv(self, key, value):
        self.lines.append(f"- **{key}:** {value}")

    def table(self, rows, headers):
        if not rows:
            return
        self.lines.append("")
        self.lines.append("| " + " | ".join(str(h) for h in headers) + " |")
        self.lines.append("|" + "|".join("---" for _ in headers) + "|")
        for r in rows:
            self.lines.append("| " + " | ".join("" if v is None else str(v) for v in r) + " |")
        self.lines.append("")

    def code(self, text):
        self.lines.append(f"\n```\n{text}\n```\n")

    def save(self):
        path = os.path.join(self.out_dir, "REPORT.md")
        header = (f"# {self.title}\n\n_estagio `{self.stage}` -- "
                  f"gerado em {time.strftime('%Y-%m-%d %H:%M:%S')}_\n")
        with open(path, "w") as fh:
            fh.write(header + "\n".join(self.lines) + "\n")
        print(f"\n[relatorio] {path}")
        return path


# =============================================================================
# checkpoint / retomada / parametros
# =============================================================================

def save_run_params(stage, args, out_dir):
    """
    Registra os parametros EXATOS desta execucao. Mantem historico (lista),
    nao sobrescreve -- o 'report' final consegue mostrar TODAS as execucoes
    de um estagio, nao so a ultima. Sem isso, um relatorio de metodologia
    teria que assumir os defaults do codigo, que podem nao bater com o que
    foi de fato rodado.
    """
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f".{stage}_params.json")
    params = {k: v for k, v in vars(args).items() if k not in ("func", "stage")}

    history = []
    if os.path.exists(path):
        try:
            with open(path) as fh:
                history = json.load(fh)
            if not isinstance(history, list):
                history = [history]
        except Exception:
            history = []

    history.append({"when": time.strftime("%Y-%m-%d %H:%M:%S"), "params": params})
    with open(path, "w") as fh:
        json.dump(history, fh, indent=2, default=str)


class Progress:
    """Registra o que cada estagio concluiu, com traceback do que quebrou."""

    def __init__(self, stage, out_dir):
        self.path = os.path.join(out_dir, f".{stage}_checkpoint.json")
        self.stage = stage
        os.makedirs(out_dir, exist_ok=True)
        self.data = {"stage": stage, "done": [], "failed": {}, "started": None}
        if os.path.exists(self.path):
            try:
                with open(self.path) as fh:
                    self.data = json.load(fh)
            except Exception:
                pass
        self.data["started"] = time.strftime("%Y-%m-%d %H:%M:%S")

    def is_done(self, item):
        return item in self.data["done"]

    def mark_done(self, item):
        if item not in self.data["done"]:
            self.data["done"].append(item)
        self.data["failed"].pop(item, None)
        self._save()

    def mark_failed(self, item, err):
        self.data["failed"][item] = {
            "error": str(err),
            "traceback": traceback.format_exc()[-2000:],
            "when": time.strftime("%Y-%m-%d %H:%M:%S"),
        }
        self._save()

    def _save(self):
        with open(self.path, "w") as fh:
            json.dump(self.data, fh, indent=2)

    def report(self):
        d, f = len(self.data["done"]), len(self.data["failed"])
        print(f"\n{'=' * 70}")
        print(f"[{self.stage}] {d} concluido(s), {f} com erro")
        if f:
            print(f"\nFALHARAM (detalhes em {self.path}):")
            for item, info in self.data["failed"].items():
                print(f"  {item}: {info['error']}")
            print("\nRode o mesmo comando de novo -- os concluidos sao pulados.")
        return f == 0


# =============================================================================
# utilidades de sequencia
# =============================================================================

def read_fasta(path, keep_desc=False):
    out, header, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    out[header] = "".join(buf)
                header = line[1:] if keep_desc else line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if header is not None:
        out[header] = "".join(buf)
    return out


def write_fasta(path, records, width=60):
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w") as fh:
        for h, s in records:
            fh.write(f">{h}\n")
            for i in range(0, len(s), width):
                fh.write(s[i:i + width] + "\n")


def parse_chopping(s):
    """'1-169_323-361,171-318' -> [[(1,169),(323,361)], [(171,318)]]"""
    out = []
    for dom in s.strip().split(","):
        segs = []
        for seg in dom.split("_"):
            a, b = seg.split("-")
            segs.append((int(a), int(b)))
        out.append(segs)
    return out


def unit_list(d, ext=".fasta"):
    return sorted(os.path.splitext(os.path.basename(p))[0]
                  for p in glob.glob(os.path.join(d, f"*{ext}")))


def real_identity(a, b):
    """Identidade real por alinhamento par a par (nao posicional)."""
    try:
        from Bio import Align
        al = Align.PairwiseAligner()
        al.mode = "local"
        al.open_gap_score, al.extend_gap_score = -10, -0.5
        try:
            from Bio.Align import substitution_matrices
            al.substitution_matrix = substitution_matrices.load("BLOSUM62")
        except Exception:
            al.match_score, al.mismatch_score = 2, -1
        aln = al.align(a, b)[0]
        x, y = str(aln[0]), str(aln[1])
        matches = sum(1 for i, j in zip(x, y) if i == j and i != "-")
        return matches / max(len(a.replace("-", "")), 1)
    except ImportError:
        n = min(len(a), len(b))
        return (sum(1 for i in range(n) if a[i] == b[i]) / max(len(a), len(b), 1)
                if n else 0.0)


# organismo no padrao do nr: ">ACC descricao [Genero especie]"
RE_ORGANISM = re.compile(r"\[([^\[\]]+)\]\s*$")


def organism_from_desc(desc):
    m = RE_ORGANISM.search(desc.strip())
    return m.group(1) if m else None


# =============================================================================
# vocabulario de contexto -- CONFIG, nao codigo
#
# O pipeline original tinha gsp/pul/xcp/pilT escritos no fonte. Isso amarra
# o fluxo a um sistema biologico so. Aqui o vocabulario e um JSON:
#
# {
#   "name": "t2ss",
#   "markers": [
#     {"label": "gspD", "class": "mandatory",
#      "patterns": ["\\b(gsp|pul|xcp|out|eps|xps|hxc)\\s*d\\b"]},
#     {"label": "gspO", "class": "accessory",
#      "patterns": ["\\bgsp\\s*o\\b", "\\bpil\\s*d\\b", "\\bcom\\s*c\\b"]}
#   ],
#   "veto": [
#     {"label": "pilT", "patterns": ["\\bpil\\s*t\\b"],
#      "why": "marcador de T4P, nao T2SS"}
#   ]
# }
#
# Sem arquivo de marcadores, o fluxo NAO fica cego: o 'profile' cai em modo
# de DESCOBERTA e reporta os termos enriquecidos em cada comunidade contra
# o fundo. E o modo certo quando voce ainda nao sabe o que procurar.
# =============================================================================

def load_markers(path):
    if not path:
        return None
    with open(path) as fh:
        cfg = json.load(fh)
    for m in cfg.get("markers", []) + cfg.get("veto", []):
        m["_re"] = [re.compile(p, re.I) for p in m["patterns"]]
    cfg.setdefault("name", os.path.splitext(os.path.basename(path))[0])
    return cfg


def evaluate_context(texts, cfg):
    """
    Descreve a vizinhanca segundo o vocabulario configurado.

    Retorna CONTAGENS e rotulos, nunca um veredito -- a decisao de quantos
    marcadores bastam fica com quem analisa, depois de ver a distribuicao
    real. Quem transforma isso em aceito/rejeitado e o 'select', e so se
    voce pedir.
    """
    if cfg is None:
        return None

    found, classes, vetoed = set(), defaultdict(set), []
    for raw in texts:
        if not isinstance(raw, str) or not raw:
            continue
        for m in cfg.get("veto", []):
            if any(r.search(raw) for r in m["_re"]):
                vetoed.append(m["label"])
        for m in cfg.get("markers", []):
            if any(r.search(raw) for r in m["_re"]):
                found.add(m["label"])
                classes[m.get("class", "marker")].add(m["label"])

    return {
        "markers": sorted(found),
        "n_markers": len(found),
        "by_class": {k: sorted(v) for k, v in classes.items()},
        "n_by_class": {k: len(v) for k, v in classes.items()},
        "veto": sorted(set(vetoed)),
        "has_veto": bool(vetoed),
    }


# =============================================================================
# descoberta de vocabulario -- o modo sem config
#
# Para cada grupo (comunidade), compara a frequencia de cada termo dentro do
# grupo com a frequencia dele no conjunto todo. Termos muito acima do fundo
# sao o que DISTINGUE aquele grupo.
#
# Usa log-odds com suavizacao (add-1). Nao e teste estatistico -- e ranking
# exploratorio, para voce olhar e decidir o que virou hipotese.
# =============================================================================

def discover_vocabulary(items_by_group, top=15, min_count=3):
    """
    items_by_group: {grupo: [termo, termo, ...]}  (termos repetidos contam)
    -> {grupo: [(termo, n_no_grupo, freq_grupo, freq_fundo, log_odds), ...]}
    """
    bg = Counter()
    per_group = {}
    for g, items in items_by_group.items():
        c = Counter(t for t in items if t)
        per_group[g] = c
        bg.update(c)

    total_bg = sum(bg.values()) or 1
    out = {}
    for g, c in per_group.items():
        total_g = sum(c.values()) or 1
        rows = []
        for term, n in c.items():
            if n < min_count:
                continue
            f_g = (n + 1) / (total_g + 1)
            f_bg = (bg[term] + 1) / (total_bg + 1)
            rows.append((term, n, f_g, f_bg, math.log2(f_g / f_bg)))
        rows.sort(key=lambda r: (-r[4], -r[1]))
        out[g] = rows[:top]
    return out


def split_terms(value):
    """Quebra celulas de anotacao ('pfam+pfam', 'a; b') em termos."""
    if not isinstance(value, str) or not value.strip():
        return []
    parts = re.split(r"[+;,|]", value)
    return [p.strip() for p in parts if p.strip()]


# aliases de coluna tolerados na tabela de busca -- em vez de exigir que a
# tabela chegue com os nomes exatos, o famflow reconhece as variantes mais
# comuns. 'cov' -> 'coverage' porque e o nome usado no projeto de calicinas.
COLUMN_ALIASES = {"cov": "coverage"}


def normalize_search_columns(df):
    ren = {k: v for k, v in COLUMN_ALIASES.items()
           if k in df.columns and v not in df.columns}
    return df.rename(columns=ren) if ren else df


def add_profile_coverage(df, seed_dir=None, unit_name=None):
    """
    Preenche df['qcov'] (cobertura do PERFIL, nao do alvo) e devolve
    (df, fonte). Tenta, em ordem de preferencia:

    1. 'coverage' ja presente na tabela (inclui o alias 'cov') -- confia no
       numero que a pessoa trouxe, e o mesmo numero que o 'inspect' plotou.
    2. 'aln_hmm_length' + 'qlen', ambos JA NA TABELA -- nao precisa de
       nenhum arquivo externo.
    3. 'qstart'/'qend' + 'qlen', JA NA TABELA -- idem, calcula o tamanho do
       trecho alinhado pela diferenca das coordenadas.
    4. 'aln_hmm_length' + fasta semente em seed_dir/{unit_name}.fasta --
       unico caso que ainda depende de um arquivo externo, mantido para
       quando a tabela nao trouxer 'qlen'.

    Se nada disso existir, devolve o df sem 'qcov' e fonte=None -- quem
    chamou decide se isso e erro ou so um filtro que fica sem efeito.
    """
    if "coverage" in df.columns:
        df = df.copy()
        df["qcov"] = df["coverage"]
        return df, "coluna 'coverage' da tabela"

    if "aln_hmm_length" in df.columns and "qlen" in df.columns:
        df = df.copy()
        df["qcov"] = df["aln_hmm_length"] / df["qlen"]
        return df, "aln_hmm_length / qlen (tabela)"

    if {"qstart", "qend", "qlen"}.issubset(df.columns):
        df = df.copy()
        df["qcov"] = (df["qend"] - df["qstart"]).abs() / df["qlen"]
        return df, "(qend - qstart) / qlen (tabela)"

    if "aln_hmm_length" in df.columns and seed_dir and unit_name:
        qpath = os.path.join(seed_dir, f"{unit_name}.fasta")
        if os.path.exists(qpath):
            qlen = len(next(iter(read_fasta(qpath).values())))
            if qlen:
                df = df.copy()
                df["qcov"] = df["aln_hmm_length"] / qlen
                return df, f"aln_hmm_length / comprimento de {qpath}"

    return df, None


# =============================================================================
# seed -- normaliza a entrada
#
# Aceita qualquer ponto de partida: uma sequencia, um conjunto de sementes,
# um HMM ja construido, ou um fasta que voce quer fatiar em dominios.
# Tudo vira o mesmo contrato: um arquivo por UNIDADE em 01_seed/.
#
# "Unidade" e so o nome do que voce esta perseguindo. Pode ser um dominio,
# uma proteina, uma familia inteira. O fluxo nao se importa.
# =============================================================================

def stage_seed(args):
    out_dir = args.out_dir
    save_run_params("seed", args, out_dir)
    prog = Progress("seed", out_dir)
    rep = Report("seed", out_dir, "Sementes")
    made = []

    if args.hmm:
        # HMM pronto: so registra, nada a construir
        for h in args.hmm:
            name = args.unit or os.path.splitext(os.path.basename(h))[0]
            dst = os.path.join(out_dir, f"{name}.hmm")
            os.makedirs(out_dir, exist_ok=True)
            if os.path.abspath(h) != os.path.abspath(dst):
                with open(h) as s, open(dst, "w") as d:
                    d.write(s.read())
            print(f"[ok] {name}: HMM pronto registrado ({dst})")
            made.append((name, "hmm-pronto", "-"))
            prog.mark_done(name)

    if args.input:
        recs = read_fasta(args.input, keep_desc=True)
        if not recs:
            sys.exit(f"[erro] fasta vazio: {args.input}")

        if args.chopping:
            # modo dominio: fatia UMA sequencia nos intervalos dados
            if len(recs) != 1:
                sys.exit("[erro] --chopping exige um fasta de sequencia unica")
            desc, seq = next(iter(recs.items()))
            base = args.unit or desc.split()[0]
            for i, segs in enumerate(parse_chopping(args.chopping), 1):
                name = f"{base}_dominio_{i}"
                if prog.is_done(name) and not args.force:
                    continue
                try:
                    if max(b for _, b in segs) > len(seq):
                        print(f"  [ATENCAO] {name}: chopping excede o comprimento")
                    dseq = "".join(seq.upper()[a - 1:b] for a, b in segs)
                    if not dseq:
                        raise ValueError("dominio vazio apos o corte")
                    coords = "_".join(f"{a}-{b}" for a, b in segs)
                    write_fasta(os.path.join(out_dir, f"{name}.fasta"),
                                [(f"{name} residuos={coords}", dseq)])
                    print(f"[ok] {name}: {len(dseq)} aa ({coords})")
                    made.append((name, "dominio", f"{len(dseq)} aa"))
                    prog.mark_done(name)
                except Exception as e:
                    print(f"[ERRO] {name}: {e}")
                    prog.mark_failed(name, e)
        else:
            # modo conjunto: o fasta inteiro e UMA unidade
            name = args.unit or os.path.splitext(os.path.basename(args.input))[0]
            write_fasta(os.path.join(out_dir, f"{name}.fasta"), list(recs.items()))
            print(f"[ok] {name}: {len(recs)} sequencia(s) semente")
            made.append((name, "conjunto", f"{len(recs)} seq"))
            prog.mark_done(name)

    if not made:
        sys.exit("[erro] nada a fazer -- use --input e/ou --hmm")

    rep.h("Unidades criadas")
    rep.table(made, ["unidade", "tipo", "tamanho"])
    rep.save()
    prog.report()


# =============================================================================
# search -- hmmbuild (se preciso) + hmmsearch
# =============================================================================

def _hmmsearch_unit(name, seed_dir, out_dir, db, cpus, evalue, incE):
    """
    hmmbuild + hmmsearch via pyhmmer.

    Por que pyhmmer e nao mmseqs: o prefilter do mmseqs carrega um indice de
    k-mer do banco inteiro na memoria antes de comecar (~180GB+ pro nr) e
    morre no cgroup. HMMER processa o banco em fluxo. Mais lento por residuo,
    sem o modo de falha catastrofico.

    Se ja existe {name}.hmm em 01_seed/, ele e usado como esta -- assim voce
    pode trazer um modelo construido fora do fluxo sem gambiarra.
    """
    import pyhmmer as ph
    from rotifer.devel.alpha import epsoares as rdae
    import rotifer.devel.beta.sequence as rdbs

    out_tsv = os.path.join(out_dir, f"{name}.tsv")
    hmm_dir = os.path.join(out_dir, "_hmm")
    os.makedirs(hmm_dir, exist_ok=True)
    hmm_path = os.path.join(hmm_dir, f"{name}.hmm")

    t0 = time.time()
    given = os.path.join(seed_dir, f"{name}.hmm")
    if os.path.exists(given):
        hmm_path, origem = given, "HMM fornecido"
    else:
        query = os.path.join(seed_dir, f"{name}.fasta")
        seqobj = rdbs.sequence(query)
        rdae.hmmbuild(seqobj, hmm_name=name, save=hmm_path)
        origem = f"hmmbuild de {len(seqobj.df)} seq"

    # NAO usar rdae.hmmsearch(): ela checa hmm_file.is_pressed antes de ler o
    # modelo, e a checagem da falso positivo em HMM recem-salvo (nunca passou
    # por hmmpress) -- cai no ramo optimized_profiles() e quebra. Para 1 HMM
    # por busca, ler direto resolve.
    with ph.plan7.HMMFile(hmm_path) as hf:
        hmms = list(hf)
    kw = {}
    if evalue is not None:
        kw.update(E=evalue, domE=evalue)
    if incE is not None:
        kw.update(incE=incE, incdomE=incE)
    with ph.easel.SequenceFile(db, digital=True,
                               alphabet=ph.easel.Alphabet.amino()) as sf:
        hits = list(ph.hmmer.hmmsearch(hmms, sf, cpus=cpus, **kw))

    df = rdae.pyhmmer_to_df(hits)
    df.to_csv(out_tsv, sep="\t", index=False)
    return {"name": name, "n_hits": len(df), "seconds": time.time() - t0,
            "origem": origem}


def stage_search(args):
    names = sorted(set(unit_list(args.seed_dir, ".fasta")) |
                   set(unit_list(args.seed_dir, ".hmm")))
    if not names:
        sys.exit(f"[erro] nenhuma unidade em {args.seed_dir} (rode 'seed')")
    if not os.path.exists(args.db) and not glob.glob(args.db + "*"):
        sys.exit(f"[erro] banco nao encontrado: {args.db}")

    os.makedirs(args.out_dir, exist_ok=True)
    save_run_params("search", args, args.out_dir)
    prog = Progress("search", args.out_dir)
    rep = Report("search", args.out_dir, "Busca no banco")

    pending = [n for n in names if not prog.is_done(n) or args.force]
    if not pending:
        print("[fim] todas as buscas ja concluidas")
        prog.report()
        return

    ev = "default do HMMER (~10)" if args.evalue is None else args.evalue
    print(f"[info] banco: {args.db}")
    print(f"[info] e-value de reporte: {ev}")
    print(f"[info] {len(pending)}/{len(names)} pendente(s), {args.workers} "
          f"em paralelo, {args.cpus_per_worker} cpu(s) cada")

    rows, done = [], 0
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(_hmmsearch_unit, n, args.seed_dir, args.out_dir,
                            args.db, args.cpus_per_worker, args.evalue,
                            args.inc_evalue): n for n in pending}
        for f in concurrent.futures.as_completed(futs):
            name = futs[f]
            done += 1
            try:
                r = f.result()
                print(f"[{done}/{len(pending)}] {name}: {r['n_hits']} hit(s) "
                      f"({r['seconds']:.0f}s, {r['origem']})")
                rows.append((name, r["n_hits"], f"{r['seconds']:.0f}s", r["origem"]))
                prog.mark_done(name)
            except Exception as e:
                print(f"[{done}/{len(pending)}] {name}: ERRO {e}")
                prog.mark_failed(name, e)

    rep.h("Parametros")
    rep.kv("banco", args.db)
    rep.kv("e-value de reporte", ev)
    rep.kv("e-value de inclusao", args.inc_evalue or "default")
    rep.h("Resultados por unidade")
    rep.table(rows, ["unidade", "hits", "tempo", "modelo"])
    rep.p("\n> O numero de hits aqui e do e-value de REPORTE. Filtros mais "
          "rigorosos (e-value, cobertura do perfil) podem ser aplicados no "
          "estagio `collect` sem refazer a busca.")
    rep.save()
    prog.report()


# =============================================================================
# inspect -- score x cobertura por hit, ANTES de qualquer corte
#
# Visualiza o funil de hmmsearch cru: cada ponto e um hit, coverage do
# PERFIL no eixo x, score no eixo y. E o mesmo grafico usado para decidir
# --max-evalue/--min-coverage no 'collect' -- rodar aqui primeiro evita
# escolher esses parametros as cegas.
#
# Ordenacao: os pontos sao desenhados na ordem dada por --sort-by (default
# 'evalue'), do PIOR para o MELHOR -- quem e desenhado por ultimo fica por
# cima. Isso importa de verdade quando ha coluna 'iteration' (busca
# incremental tipo jackhmmer, onde um hit que ja apareceu numa iteracao
# nao aparece mais nas seguintes): sem essa ordem, a massa de pontos da
# iteracao 1 pode enterrar visualmente os poucos pontos novos das
# iteracoes seguintes -- que sao justamente os mais interessantes, porque
# sao os hits que so o perfil relaxado das rodadas seguintes capturou.
# =============================================================================

def _fmt_cell(v):
    if isinstance(v, float):
        return f"{v:.2f}" if not v.is_integer() else f"{int(v)}"
    return str(v)


def _save_table_png(plt, df, path, title, dpi=150, highlight_row=None):
    """
    Renderiza um DataFrame como imagem -- pronta pra colar num slide junto
    do scatter. A linha do corte escolhido (se houver) fica destacada, no
    mesmo espirito da linha solida preta no scatter: a mesma decisao
    marcada nos dois lugares.
    """
    fig, ax = plt.subplots(figsize=(1.3 * len(df.columns) + 1,
                                    0.38 * (len(df) + 1) + 0.6))
    ax.axis("off")
    cell_text = [[_fmt_cell(v) for v in row] for row in df.values]
    # bbox=[0,0,1,1] preenche a area inteira do eixo -- sem isso, loc='center'
    # centraliza a tabela deixando faixas de vazio em cima e embaixo
    tbl = ax.table(cellText=cell_text, colLabels=list(df.columns),
                   cellLoc="center", bbox=[0, 0, 1, 1])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#dddddd")
        if r == 0:
            cell.set_facecolor("#2ca02c")
            cell.set_text_props(color="white", weight="bold")
        elif highlight_row is not None and r - 1 == highlight_row:
            cell.set_facecolor("#fff3b0")
        else:
            cell.set_facecolor("#f7f7f7" if r % 2 == 0 else "white")
    fig.subplots_adjust(top=0.85)
    fig.suptitle(title, fontsize=11, y=0.97)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def stage_inspect(args):
    import pandas as pd
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("[erro] matplotlib nao instalado (pip install matplotlib "
                  "--break-system-packages)")

    names = unit_list(args.search_dir, ".tsv")
    if not names:
        sys.exit(f"[erro] nenhum .tsv em {args.search_dir} (rode 'search')")

    save_run_params("inspect", args, args.out_dir)
    rep = Report("inspect", args.out_dir, "Score x cobertura (pre-corte)")
    os.makedirs(args.out_dir, exist_ok=True)

    sort_col = args.sort_by  # 'evalue' (default) ou 'score'
    ascending = sort_col == "evalue"  # pior evalue = maior; pior score = menor

    # ciclo de estilo estavel: mesma iteracao sempre com a mesma cor/marcador,
    # independente de quantas iteracoes existirem
    palette = [("#2ca02c", "o"), ("#1f77b4", "s"), ("#d62728", "^"),
              ("#ff7f0e", "D"), ("#9467bd", "v"), ("#8c564b", "P")]

    cov_grid = [0.0, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9]
    ev_grid = [None, 1e-3, 1e-5, 1e-10]

    rows = []
    for name in names:
        tsv = os.path.join(args.search_dir, f"{name}.tsv")
        df = pd.read_csv(tsv, sep="\t")
        df = normalize_search_columns(df)
        if len(df) == 0:
            print(f"[skip] {name}: tsv vazio")
            continue

        if sort_col not in df.columns:
            print(f"[aviso] {name}: coluna '{sort_col}' ausente, usando a "
                  f"outra disponivel")
            sort_col_eff = "score" if sort_col == "evalue" else "evalue"
        else:
            sort_col_eff = sort_col

        # cobertura do PERFIL (nao do alvo) -- mesma conta do 'collect', para
        # que a linha de corte desenhada aqui seja a mesma que --min-coverage
        # vai de fato aplicar
        df, cov_src = add_profile_coverage(df, args.seed_dir, name)
        if cov_src is None:
            print(f"[aviso] {name}: sem como calcular cobertura do perfil "
                  f"(faltam qlen+aln_hmm_length/qstart+qend na tabela, e "
                  f"tambem 'coverage' e fasta semente) -- pulando")
            continue
        df["coverage"] = df["qcov"]

        # ordena do PIOR para o MELHOR: o melhor fica desenhado por cima
        asc = sort_col_eff == "evalue"
        df = df.sort_values(sort_col_eff, ascending=asc)

        # sensibilidade do corte: quantas sequencias sobram em cada
        # combinacao de --min-coverage x --max-evalue -- calculada ANTES do
        # scatter para poder desenhar as mesmas linhas nos dois lugares
        sens = []
        for c in cov_grid:
            row = {"min_coverage": c}
            for e in ev_grid:
                sub = df[df["coverage"] >= c]
                if e is not None and "evalue" in sub.columns:
                    sub = sub[sub["evalue"] <= e]
                label = "sem_filtro" if e is None else f"evalue<={e:g}"
                row[label] = len(sub)
            sens.append(row)
        sens_df = pd.DataFrame(sens)
        sens_path = os.path.join(args.out_dir, f"{name}_sensibilidade.tsv")
        sens_df.to_csv(sens_path, sep="\t", index=False)

        # equivalente em SCORE de cada corte de e-value testado -- nao existe
        # uma formula fechada sem saber o tamanho exato do banco, entao pega
        # o menor score observado nos dados que de fato passam naquele
        # e-value. E o mesmo numero, so que emprestado dos proprios dados.
        ev_score_lines = []
        for e in ev_grid:
            if e is None or "evalue" not in df.columns:
                continue
            sub = df[df["evalue"] <= e]
            if len(sub):
                ev_score_lines.append((e, float(sub["score"].min())))

        if "iteration" in df.columns:
            # desenha as iteracoes na ordem inversa da numeracao (as
            # ultimas primeiro), para que os hits NOVOS de cada iteracao
            # (que sao poucos, por definicao de busca incremental) fiquem
            # por cima da massa da iteracao 1, em vez de soterrados por ela
            iters = sorted(df["iteration"].unique(), reverse=True)
        else:
            df["iteration"] = 1
            iters = [1]

        # cor/marcador fixos por NUMERO da iteracao, nao pela ordem de
        # desenho -- senao a Iter 1 muda de cor conforme quais outras
        # iteracoes existem no tsv, o que so confunde a leitura
        def style_for(it):
            try:
                idx = int(it) - 1
            except (TypeError, ValueError):
                idx = sorted(df["iteration"].unique()).index(it)
            return palette[idx % len(palette)]

        fig, ax = plt.subplots(figsize=(7, 6))
        for it in iters:
            color, marker = style_for(it)
            sub = df[df["iteration"] == it]
            if len(sub) == 0:
                continue
            ax.scatter(sub["coverage"], sub["score"], s=25, alpha=0.6,
                      color=color, marker=marker, label=f"Iter {it}",
                      edgecolors="none")

        # linhas FINAS de referencia -- todos os cortes que a tabela de
        # sensibilidade testou, para servir de contexto visual. Ficam
        # discretas (pontilhadas, baixo alpha) porque o corte de fato
        # escolhido (abaixo) e o que precisa saltar aos olhos.
        for c in cov_grid:
            if c == 0.0:
                continue
            ax.axvline(c, linestyle=":", color="0.6", alpha=0.5, linewidth=0.8)
            ax.text(c, 0.98, f"{c:g}", transform=ax.get_xaxis_transform(),
                   ha="center", va="top", fontsize=7, color="0.4")
        for e, sc in ev_score_lines:
            ax.axhline(sc, linestyle=":", color="steelblue", alpha=0.4, linewidth=0.8)
            ax.text(1.005, sc, f"e\u2264{e:g}", transform=ax.get_yaxis_transform(),
                   ha="left", va="center", fontsize=7, color="steelblue",
                   clip_on=False)

        # o corte que de fato foi escolhido -- solido, preto, por cima das
        # linhas de referencia pontilhadas
        if args.min_coverage is not None:
            ax.axvline(args.min_coverage, linestyle="--", color="black", alpha=0.8,
                      linewidth=1.3)
        if args.score_cutoff is not None:
            ax.axhline(args.score_cutoff, linestyle="--", color="black", alpha=0.8,
                      linewidth=1.3)
        ax.axhline(0, linestyle="--", color="gray", alpha=0.4)
        ax.axvline(0, linestyle="--", color="gray", alpha=0.4)
        ax.set_xlabel("Coverage (perfil)")
        ax.set_ylabel("Score")
        ax.set_title(f"{name} -- {len(df)} hit(s), ordenado por {sort_col_eff}")
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        if len(by_label) > 1 or "iteration" in df.columns:
            ax.legend(by_label.values(), by_label.keys(), frameon=False)
        fig.tight_layout()
        png = os.path.join(args.out_dir, f"{name}_score_coverage.png")
        fig.savefig(png, dpi=args.dpi)
        plt.close(fig)

        # densidade -- o scatter satura com muitos pontos e some a
        # informacao de QUANTOS hits ha em cada regiao (ex: a listra em
        # coverage~1.0 pode ser 500 ou 8000 sequencias, indistinguivel no
        # scatter). Hexbin em log de contagem resolve.
        fig2, ax2 = plt.subplots(figsize=(7, 6))
        hb = ax2.hexbin(df["coverage"], df["score"], gridsize=40,
                        bins="log", mincnt=1, cmap="viridis")
        fig2.colorbar(hb, ax=ax2, label="log10(n hits)")
        if args.min_coverage is not None:
            ax2.axvline(args.min_coverage, linestyle="--", color="white", alpha=0.7)
        ax2.set_xlabel("Coverage (perfil)")
        ax2.set_ylabel("Score")
        ax2.set_title(f"{name} -- densidade ({len(df)} hit(s))")
        fig2.tight_layout()
        png2 = os.path.join(args.out_dir, f"{name}_score_coverage_density.png")
        fig2.savefig(png2, dpi=args.dpi)
        plt.close(fig2)

        # figura da tabela de sensibilidade -- pronta pra slide, com a
        # linha do corte escolhido destacada (mesma decisao marcada tanto
        # aqui quanto na linha solida do scatter acima)
        highlight_row = None
        if args.min_coverage is not None:
            highlight_row = min(range(len(cov_grid)),
                                key=lambda i: abs(cov_grid[i] - args.min_coverage))
        table_png = os.path.join(args.out_dir, f"{name}_sensibilidade.png")
        _save_table_png(plt, sens_df, table_png,
                        f"{name} -- sequencias restantes por corte",
                        dpi=args.dpi, highlight_row=highlight_row)

        n_low_cov = int((df["coverage"] < 0.5).sum())
        rows.append({"unidade": name, "hits": len(df),
                     "coverage_fonte": cov_src, "ordenado_por": sort_col_eff,
                     "cov_mediana": float(df["coverage"].median()),
                     "score_mediano": float(df["score"].median()),
                     "hits_cov<0.5": n_low_cov,
                     "png": os.path.basename(png),
                     "table_png": os.path.basename(table_png)})
        print(f"[ok] {name}: {len(df)} hit(s), cov mediana="
              f"{df['coverage'].median():.2f}, {n_low_cov} com cov<0.5 "
              f"-> {png}, {table_png}")

    if not rows:
        sys.exit("[erro] nenhuma unidade produziu grafico")

    rep.h("Parametros")
    rep.kv("ordenado por", sort_col)
    rep.kv("linha de corte (cobertura)", args.min_coverage or "nao desenhada")
    rep.kv("linha de corte (score)", args.score_cutoff or "nao desenhada")
    rep.h("Resumo por unidade")
    rep.table([(r["unidade"], r["hits"], f"{r['cov_mediana']:.2f}",
               f"{r['score_mediano']:.1f}", r["hits_cov<0.5"], r["png"],
               r["table_png"])
              for r in rows],
              ["unidade", "hits", "cov mediana", "score mediano",
               "hits com cov<0.5", "scatter", "tabela (png)"])
    rep.p("\n> `hits com cov<0.5` e a cauda de motivo isolado -- o que "
          "`--min-coverage` no `collect` remove. As linhas pontilhadas no "
          "scatter marcam TODOS os cortes testados em `*_sensibilidade.tsv`; "
          "a linha solida preta marca o corte que voce de fato passou em "
          "`--min-coverage`/`--score-cutoff`. `*_sensibilidade.png` e a "
          "mesma tabela pronta pra colar num slide, com a linha escolhida "
          "destacada.")
    rep.p("\n> Confirme se a cobertura acima e do PERFIL, nao do alvo: "
          "cobertura do alvo perto de 1.0 pode ser so proteina fragmentada "
          "no banco, nao hit completo.")
    rep.save()


# =============================================================================
# collect -- recupera as sequencias e INICIA a tabela central
# =============================================================================

def _fetch_full_lengths(ids, db, out_dir, tag):
    """
    Busca as sequencias INTEIRAS so para extrair o comprimento delas.

    Usado quando o que foi coletado e um recorte (envelope) e o comprimento
    da proteina de origem nao e conhecido de outra forma. Mesmo custo de
    E/S que recuperar em --mode full -- e o preco de saber se um envelope
    de 80 residuos e 90% de uma proteina pequena ou 15% de uma grande.

    NAO existe atalho pela tabela de busca para isso: colunas como
    'talilen' (total alignment length) sao o comprimento do trecho
    ALINHADO no alvo -- primo de 'aln_hmm_length' mas do lado do alvo,
    redundante com eend-estart -- nao o tamanho da proteina inteira. Usar
    'talilen' aqui daria um target_coverage sempre proximo de 1.0, mascarando
    exatamente a diferenca que este numero deveria mostrar.
    """
    import rotifer.devel.beta.sequence as rdbs
    ids = list(dict.fromkeys(ids))  # unicos, preservando ordem
    if not ids:
        return {}
    tmp = os.path.join(out_dir, f".{tag}_fulllen.tmp.fasta")
    rdbs.sequence(ids, local_database_path=db).to_file(tmp)
    lens = {h: len(s) for h, s in read_fasta(tmp).items()}
    os.remove(tmp)
    return lens


def stage_collect(args):
    import pandas as pd
    # o rotifer so e importado no ramo que precisa dele: o atalho
    # --from-fasta tem que funcionar em qualquer maquina, sem o stack todo

    save_run_params("collect", args, args.out_dir)
    prog = Progress("collect", args.out_dir)
    rep = Report("collect", args.out_dir, "Recuperacao dos hits")
    os.makedirs(args.out_dir, exist_ok=True)
    rows, all_meta = [], []

    # ATALHO: entrar no fluxo direto por um fasta pronto, sem passar por
    # seed/search. E o caso de quem ja fez a busca por fora.
    if args.from_fasta:
        name = args.unit or os.path.splitext(os.path.basename(args.from_fasta))[0]
        recs = read_fasta(args.from_fasta, keep_desc=True)
        out = os.path.join(args.out_dir, f"{name}.fasta")
        clean = {h.split()[0]: s for h, s in recs.items()}
        write_fasta(out, list(clean.items()))

        # mesmo raciocinio do ramo normal: em 'full' o que foi importado JA
        # e a proteina inteira (target_length trivial); em 'envelope' o
        # fasta e um recorte e o tamanho total precisa vir de um fetch no
        # banco pelos proprios IDs do fasta (assume que sao accessions
        # validas). --target-length skip evita esse fetch.
        target_len = {}
        if args.mode == "full":
            target_len = {h: len(s) for h, s in clean.items()}
        elif args.target_length == "auto":
            try:
                target_len = _fetch_full_lengths(clean.keys(), args.db,
                                                 args.out_dir, name)
            except Exception as e:
                print(f"  [aviso] {name}: nao consegui recuperar tamanho "
                      f"total das proteinas ({type(e).__name__}: {e})")

        meta = pd.DataFrame([{
            "id": h.split()[0], "unit": name, "length": len(s),
            "target_length": target_len.get(h.split()[0]),
            "target_coverage": (len(s) / target_len[h.split()[0]]
                                if target_len.get(h.split()[0]) else None),
            "organism": organism_from_desc(h),
            "description": h.split(" ", 1)[1] if " " in h else None,
        } for h, s in recs.items()])
        all_meta.append(meta)
        rows.append((name, len(recs), "-", "importado", len(recs)))
        n_no_tlen = sum(1 for h in clean if target_len.get(h) is None)
        print(f"[ok] {name}: {len(recs)} sequencia(s) importada(s) de "
              f"{args.from_fasta} [{args.mode}]")
        if args.mode == "envelope" and n_no_tlen:
            print(f"       [aviso] {n_no_tlen}/{len(clean)} sem target_length")
        prog.mark_done(name)

    else:
        from rotifer.devel.alpha import epsoares as rdae

        # unidades vem do que EXISTE em search-dir -- e o insumo que o
        # 'collect' realmente precisa. O seed-dir e so enriquecimento
        # opcional (a ancora, se a semente estiver la); nao deveria ser
        # pre-requisito para processar uma busca que ja foi feita.
        names = unit_list(args.search_dir, ".tsv")
        if not names:
            sys.exit(f"[erro] nenhum .tsv em {args.search_dir} (rode 'search', "
                      f"ou copie sua tabela para {args.search_dir}/<unidade>.tsv)")

        for name in names:
            if prog.is_done(name) and not args.force:
                continue
            tsv = os.path.join(args.search_dir, f"{name}.tsv")
            if not os.path.exists(tsv):
                print(f"[skip] {name}: sem resultado de busca")
                continue

            try:
                df = pd.read_csv(tsv, sep="\t")
                df = normalize_search_columns(df)
                n_raw = len(df)

                if args.max_evalue is not None:
                    df = df[df["evalue"] <= args.max_evalue]
                n_ev = len(df)

                # COBERTURA DO PERFIL (nao do alvo): quanto do HMM foi de fato
                # alinhado. Sem isso, um hit que bate so um motivo curto de 15
                # residuos de um perfil de 200 passa no e-value mas e ruido
                # para o hmmbuild final -- entra so gap no alinhamento.
                #
                # Calculada sempre que a propria tabela tiver o necessario
                # (qlen + aln_hmm_length ou qstart/qend, ou uma coluna
                # 'coverage'/'cov' ja pronta) -- sem precisar de nenhum
                # arquivo externo. So cai para a fasta semente se a tabela
                # nao trouxer 'qlen'.
                df, cov_src = add_profile_coverage(df, args.seed_dir, name)
                if args.min_coverage is not None:
                    if cov_src is not None:
                        df = df[df["qcov"] >= args.min_coverage]
                    else:
                        print(f"  [aviso] {name}: sem como calcular cobertura "
                              f"do perfil (faltam qlen+aln_hmm_length/"
                              f"qstart+qend, coverage/cov, ou fasta semente); "
                              f"--min-coverage ignorado")
                n_cov = len(df)

                # um hit por alvo (melhor score), senao o mesmo alvo e
                # recortado mais de uma vez
                df = df.sort_values("score", ascending=False).drop_duplicates("sequence")
                if args.max_seqs and len(df) > args.max_seqs:
                    print(f"  [aviso] {name}: {len(df)} hits, truncando para "
                          f"{args.max_seqs}")
                    df = df.head(args.max_seqs)
                if len(df) == 0:
                    raise ValueError("nenhum hit apos os filtros")
                n_targets = len(df)

                targets = df["sequence"].astype(str).tolist()

                if args.mode == "envelope":
                    if not {"estart", "eend"}.issubset(df.columns):
                        raise ValueError(
                            "--mode envelope exige colunas 'estart'/'eend' "
                            "(limites do alinhamento), ausentes nesta tabela "
                            "-- use --mode full para recuperar a proteina "
                            "inteira por accession")
                    # so a regiao que alinhou (env_from/env_to) + folga
                    seqobj = rdae.extract_envelope(
                        df, start="estart", end="eend",
                        expand=args.expand, local_database_path=args.db)
                else:
                    # proteina inteira -- necessario se voce quer ver
                    # diferencas ARQUITETURAIS entre familias, nao so o fold
                    import rotifer.devel.beta.sequence as rdbs
                    seqobj = rdbs.sequence(targets, local_database_path=args.db)

                out = os.path.join(args.out_dir, f"{name}.fasta")
                seqobj.to_file(out)
                recs = read_fasta(out)

                # tamanho da proteina INTEIRA, nao so do pedaco recortado.
                # No 'full' e trivial (o que foi recuperado JA e a proteina
                # inteira). No 'envelope' o recorte perde essa informacao --
                # sem ela, dado o tamanho do envelope isolado, nao da pra
                # saber se ele e 90% de uma proteina pequena ou 20% de uma
                # grande, e isso importa para diferenciar arquitetura entre
                # familias e para flagrar fragmentacao de assembly. Custa um
                # fetch a mais no banco (mesmo custo do 'full'); desative com
                # --target-length skip se a escala nao compensar.
                target_len = {}
                if args.mode == "full":
                    target_len = {h: len(s) for h, s in recs.items()}
                elif args.target_length == "auto":
                    try:
                        target_len = _fetch_full_lengths(targets, args.db,
                                                         args.out_dir, name)
                    except Exception as e:
                        print(f"  [aviso] {name}: nao consegui recuperar o "
                              f"tamanho total das proteinas ({type(e).__name__}: "
                              f"{e}) -- target_length ficara vazio")

                # rastreia, por chave final em 'recs', qual accession original
                # ela veio de -- necessario porque a ancora pode RENOMEAR uma
                # chave (self-hit exato/assumido) antes de eu conseguir casar
                # com target_len
                source = {h: h for h in recs}

                # ancora: garante que a semente esta no conjunto, identificavel.
                # Sem isso nao ha como saber depois qual comunidade e a dela.
                qmode = "sem ancora"
                qpath = os.path.join(args.seed_dir, f"{name}.fasta")
                if args.anchor and os.path.exists(qpath):
                    qseq = next(iter(read_fasta(qpath).values())).upper()
                    exact = [h for h, s in recs.items() if s.upper() == qseq]
                    if exact:
                        old = exact[0]
                        recs = {(name if h == old else h): s
                                for h, s in recs.items()}
                        source[name] = source.pop(old, old)
                        qmode = "self-hit exato"
                    elif recs:
                        top_h, top_s = max(recs.items(), key=lambda kv: len(kv[1]))
                        ident = real_identity(qseq, top_s)
                        if ident >= args.anchor_identity:
                            recs = {(name if h == top_h else h): s
                                    for h, s in recs.items()}
                            source[name] = source.pop(top_h, top_h)
                            qmode = f"assumida = melhor hit ({ident:.1%})"
                        else:
                            recs[name] = qseq
                            # ancora inserida do zero (nao veio de um alvo
                            # recuperado) -- e a propria semente, seu
                            # "tamanho total" e o dela mesma
                            target_len[name] = len(qseq)
                            source[name] = name
                            qmode = f"ancora inserida (melhor hit so {ident:.1%})"
                    else:
                        recs[name] = qseq
                        target_len[name] = len(qseq)
                        source[name] = name
                        qmode = "ancora inserida (sem hits)"
                    write_fasta(out, list(recs.items()))

                n = len(recs)
                n_missing = n_targets - n

                sub = df.set_index("sequence")
                meta = pd.DataFrame([{
                    "id": h, "unit": name, "length": len(s),
                    "target_length": target_len.get(source.get(h, h)),
                    "target_coverage": (
                        len(s) / target_len[source[h]]
                        if source.get(h) in target_len and target_len[source[h]]
                        else None),
                    "search_evalue": sub["evalue"].get(h),
                    "search_score": sub["score"].get(h),
                    "search_qcov": (sub["qcov"].get(h) if "qcov" in sub.columns
                                    else None),
                    "is_anchor": bool(args.anchor and h == name),
                } for h, s in recs.items()])
                all_meta.append(meta)

                n_no_tlen = sum(1 for h in recs if target_len.get(source.get(h, h)) is None)

                print(f"[ok] {name}: {n_raw} hit(s) -> {n} sequencia(s) "
                      f"[{args.mode}] -- {qmode}")
                if args.mode == "envelope" and n_no_tlen:
                    tag = " (--target-length skip)" if args.target_length == "skip" else ""
                    print(f"       [aviso] {n_no_tlen}/{n} sem target_length{tag}")
                if n_missing > 0:
                    pct = n_missing / n_targets
                    # o snapshot do nr usado na busca e o usado aqui para
                    # recuperar sequencia podem ser DIFERENTES -- accessions
                    # novos nao existem no mais antigo
                    flag = " <-- ATENCAO" if pct > 0.1 else ""
                    print(f"       [aviso] {n_missing} ({pct:.1%}) nao "
                          f"recuperados: provavel diferenca de snapshot do "
                          f"banco{flag}")
                rows.append((name, n_raw, n_ev, n_cov, n))
                prog.mark_done(name)

            except Exception as e:
                print(f"[ERRO] {name}: {e}")
                prog.mark_failed(name, e)

    if all_meta:
        meta = pd.concat(all_meta, ignore_index=True)
        _, full = update_table(meta)
        print(f"\n[tabela] {len(full)} sequencia(s) em {TABLE}")

    rep.h("Filtros aplicados")
    rep.kv("modo", args.mode)
    rep.kv("e-value maximo", args.max_evalue or "sem filtro")
    rep.kv("cobertura minima do perfil", args.min_coverage or "sem filtro")
    rep.kv("folga do envelope", f"{args.expand} residuos"
           if args.mode == "envelope" else "n/a")
    rep.kv("tamanho da proteina de origem (target_length)",
           "trivial (= length, modo full)" if args.mode == "full"
           else f"--target-length {args.target_length}")
    rep.h("Funil por unidade")
    rep.table(rows, ["unidade", "hits brutos", "apos e-value",
                     "apos cobertura", "recuperadas"])
    if all_meta and "target_coverage" in full.columns:
        n_tc = int(full["target_coverage"].notna().sum())
        rep.p(f"\n> `target_length`/`target_coverage` na tabela central: "
              f"{n_tc}/{len(full)} sequencia(s). `target_coverage` baixo "
              f"significa que o envelope recortado e uma fracao pequena da "
              f"proteina original -- util para separar domino isolado de "
              f"proteina multi-dominio na hora de perfilar as comunidades.")
    rep.save()
    prog.report()


# =============================================================================
# cluster -- reduz redundancia PRESERVANDO o mapa de membros
#
# O pipeline original escrevia so os representantes e jogava fora quem estava
# em cada cluster. Para um estudo de DIVERSIDADE isso e fatal: o tamanho de
# uma comunidade deixa de significar qualquer coisa, porque uma familia com
# 5000 sequencias quase identicas e uma com 5 divergentes viram ambas "1
# representante". Aqui o mapa e salvo e vira coluna na tabela central.
#
# Dois diagnosticos, no mesmo espirito do 'inspect' antes do 'collect':
#
# --sweep    ANTES de decidir --identity/--coverage: roda o mmseqs numa
#            grade de combinacoes e planta quantos clusters (comunidades)
#            cada uma produz -- um heatmap identity x coverage. Nao escreve
#            fasta nenhum, e so diagnostico. Cada combinacao e um comando
#            de mmseqs inteiro; a grade default e deliberadamente pequena
#            (4x3=12) por causa disso -- aumente com --sweep-identity/
#            --sweep-coverage se a escala permitir.
#
# (padrao)   Depois de decidido o limiar, o clustering de verdade roda e
#            SEMPRE sai com a distribuicao de tamanho dos clusters (quantos
#            tem 1 membro, quantos tem dezenas, etc) -- e o "quantas
#            comunidades em cada agrupamento" para o limiar escolhido.
# =============================================================================

def _cluster_grid_stats(seqobj, ident, cov):
    """Roda add_cluster pra uma combinacao e devolve (n_clusters, n_singletons, biggest)."""
    clustered = seqobj.add_cluster(coverage=cov, identity=ident)
    col = f"c{int(cov * 100)}i{int(ident * 100)}"
    if col not in clustered.df.columns:
        new_cols = [c for c in clustered.df.columns if c not in seqobj.df.columns]
        if not new_cols:
            raise RuntimeError("add_cluster nao criou coluna de cluster")
        col = new_cols[-1]
    sizes = clustered.df.groupby(col).size()
    return len(sizes), int((sizes == 1).sum()), int(sizes.max())


def _size_bins(sizes):
    """Agrupa tamanhos de cluster em faixas fixas, comparaveis entre unidades."""
    edges = [(1, 1, "1"), (2, 2, "2"), (3, 5, "3-5"), (6, 10, "6-10"),
             (11, 20, "11-20"), (21, 50, "21-50"), (51, 100, "51-100"),
             (101, 500, "101-500"), (501, float("inf"), "501+")]
    counts = []
    for lo, hi, label in edges:
        n = int(((sizes >= lo) & (sizes <= hi)).sum())
        if n or label in ("1", "2"):
            counts.append((label, n))
    return counts


def stage_cluster(args):
    import pandas as pd
    import rotifer.devel.beta.sequence as rdbs

    fastas = sorted(glob.glob(os.path.join(args.hits_dir, "*.fasta")))
    if not fastas:
        sys.exit(f"[erro] nenhum .fasta em {args.hits_dir} (rode 'collect')")

    save_run_params("cluster", args, args.out_dir)

    # ---- modo diagnostico: so a grade, nao clusteriza de verdade ----
    if args.sweep:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            sys.exit("[erro] matplotlib nao instalado (pip install matplotlib "
                      "--break-system-packages)")

        sweep_dir = os.path.join(args.out_dir, "_sweep")
        os.makedirs(sweep_dir, exist_ok=True)
        prog = Progress("cluster_sweep", sweep_dir)
        rep = Report("cluster_sweep", sweep_dir, "Sweep de limiares de clustering")

        idents = sorted(args.sweep_identity)
        covs = sorted(args.sweep_coverage)
        total = len(fastas) * len(idents) * len(covs)
        print(f"[info] {len(idents)} identidade(s) x {len(covs)} cobertura(s) "
              f"x {len(fastas)} unidade(s) = {total} combinacao(oes) de mmseqs")

        rows = []
        for fa in fastas:
            name = os.path.splitext(os.path.basename(fa))[0]
            seqobj = rdbs.sequence(fa)
            n_seqs = len(seqobj.df)
            for ident in idents:
                for cov in covs:
                    item = f"{name}_i{int(ident*100)}c{int(cov*100)}"
                    if prog.is_done(item) and not args.force:
                        continue
                    try:
                        n_cl, n_single, biggest = _cluster_grid_stats(seqobj, ident, cov)
                        print(f"[ok] {name} id={ident:.0%} cov={cov:.0%}: "
                              f"{n_cl} cluster(s), {n_single} singleton(s), "
                              f"maior={biggest}")
                        rows.append({"unit": name, "n_seqs": n_seqs,
                                    "identity": ident, "coverage": cov,
                                    "n_clusters": n_cl, "n_singletons": n_single,
                                    "biggest_cluster": biggest})
                        prog.mark_done(item)
                    except Exception as e:
                        print(f"[ERRO] {name} id={ident:.0%} cov={cov:.0%}: {e}")
                        prog.mark_failed(item, e)

        if not rows and not os.path.exists(os.path.join(sweep_dir, "sweep.tsv")):
            sys.exit("[erro] sweep nao produziu nenhum resultado")

        # funde com o que ja estava salvo -- rows so tem as combinacoes
        # NOVAS desta execucao (as ja feitas foram puladas pelo checkpoint
        # e nao devem sumir do tsv)
        prior_path = os.path.join(sweep_dir, "sweep.tsv")
        if not rows:
            sdf = pd.read_csv(prior_path, sep="\t")
            print(f"[info] nada novo -- reaproveitando {prior_path} "
                  f"({len(sdf)} combinacao(oes) ja computada(s))")
        else:
            new_df = pd.DataFrame(rows)
            if os.path.exists(prior_path):
                prior = pd.read_csv(prior_path, sep="\t")
                key = ["unit", "identity", "coverage"]
                prior = prior[~prior.set_index(key).index.isin(
                    new_df.set_index(key).index)]
                sdf = pd.concat([prior, new_df], ignore_index=True)
            else:
                sdf = new_df
            sdf.to_csv(prior_path, sep="\t", index=False)

        rep.h("Grade testada")
        rep.kv("identidades", ", ".join(f"{i:.0%}" for i in idents))
        rep.kv("coberturas", ", ".join(f"{c:.0%}" for c in covs))
        rep.kv("combinacoes por unidade", len(idents) * len(covs))

        for name, sub in sdf.groupby("unit"):
            piv = sub.pivot(index="identity", columns="coverage", values="n_clusters")
            fig, ax = plt.subplots(figsize=(1.2 * len(covs) + 2, 1.0 * len(idents) + 2))
            im = ax.imshow(piv.values, cmap="viridis", aspect="auto")
            ax.set_xticks(range(len(piv.columns)))
            ax.set_xticklabels([f"{c:.0%}" for c in piv.columns])
            ax.set_yticks(range(len(piv.index)))
            ax.set_yticklabels([f"{i:.0%}" for i in piv.index])
            ax.set_xlabel("coverage")
            ax.set_ylabel("identity")
            ax.set_title(f"{name} -- n. de clusters ({sub['n_seqs'].iloc[0]} seq)")
            for r in range(piv.shape[0]):
                for c in range(piv.shape[1]):
                    v = piv.values[r, c]
                    if pd.notna(v):
                        ax.text(c, r, int(v), ha="center", va="center",
                               color="white", fontsize=8)
            fig.colorbar(im, ax=ax, label="n. clusters")
            fig.tight_layout()
            png = os.path.join(sweep_dir, f"{name}_sweep_heatmap.png")
            fig.savefig(png, dpi=150)
            plt.close(fig)
            print(f"[ok] {name}: heatmap -> {png}")

            rep.h(f"{name}", level=3)
            rep.table([(f"{r['identity']:.0%}", f"{r['coverage']:.0%}",
                       r["n_clusters"], r["n_singletons"], r["biggest_cluster"])
                      for _, r in sub.sort_values(["identity", "coverage"]).iterrows()],
                      ["identity", "coverage", "n_clusters", "n_singletons",
                       "maior_cluster"])

        rep.p("\n> Cada celula do heatmap e uma chamada de mmseqs inteira. "
              "Quanto mais estavel o n. de clusters entre celulas vizinhas, "
              "menos a escolha exata do limiar importa; regioes onde o "
              "numero salta bastante de uma celula pra outra sao onde a "
              "decisao de fato muda o resultado.")
        rep.save()
        prog.report()
        print("\n[proximo passo] escolha --identity/--coverage pelo heatmap "
              "e rode 'cluster' sem --sweep para o clustering de verdade.")
        return

    # ---- modo normal: clusteriza de verdade ----
    prog = Progress("cluster", args.out_dir)
    rep = Report("cluster", args.out_dir, "Reducao de redundancia")

    overrides = {}
    for it in args.override or []:
        parts = it.split(":")
        if len(parts) != 3:
            sys.exit(f"[erro] --override mal formado: '{it}' (nome:identity:coverage)")
        overrides[parts[0]] = (float(parts[1]), float(parts[2]))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        have_plt = True
    except ImportError:
        have_plt = False
        print("[aviso] matplotlib nao instalado -- pulando os graficos de "
              "distribuicao de tamanho (tsv ainda e gerado)")

    rows, maps = [], []
    for fa in fastas:
        name = os.path.splitext(os.path.basename(fa))[0]
        if prog.is_done(name) and not args.force:
            continue
        try:
            ident, cov = overrides.get(name, (args.identity, args.coverage))
            seqobj = rdbs.sequence(fa)
            n_before = len(seqobj.df)

            clustered = seqobj.add_cluster(coverage=cov, identity=ident)
            col = f"c{int(cov * 100)}i{int(ident * 100)}"
            if col not in clustered.df.columns:
                new_cols = [c for c in clustered.df.columns
                            if c not in seqobj.df.columns]
                if not new_cols:
                    raise RuntimeError("add_cluster nao criou coluna de cluster")
                col = new_cols[-1]

            # >>> o mapa que o original descartava <<<
            m = clustered.df[["id", col]].rename(columns={col: "cluster_rep"})
            m["id"] = m["id"].astype(str)
            m["cluster_rep"] = m["cluster_rep"].astype(str)
            sizes = m.groupby("cluster_rep").size().rename("cluster_size")
            m = m.merge(sizes, left_on="cluster_rep", right_index=True)
            m["is_cluster_rep"] = m["id"] == m["cluster_rep"]
            m.to_csv(os.path.join(args.out_dir, f"{name}_clusters.tsv"),
                     sep="\t", index=False)
            maps.append(m)

            reps = clustered.df[clustered.df["id"] == clustered.df[col]]
            write_fasta(os.path.join(args.out_dir, f"{name}.fasta"),
                        [(str(r["id"]), str(r["sequence"]).replace("-", ""))
                         for _, r in reps.iterrows()])

            n_after = len(reps)
            biggest = int(sizes.max()) if len(sizes) else 0
            tag = " [OVERRIDE]" if name in overrides else ""
            print(f"[ok] {name}: {n_before} -> {n_after} "
                  f"(id>={ident:.0%} cov>={cov:.0%}) maior cluster={biggest}{tag}")
            rows.append((name, n_before, n_after, f"{ident:.0%}", f"{cov:.0%}",
                         biggest))

            # distribuicao de tamanho -- "quantas comunidades em cada
            # agrupamento" pro limiar que de fato foi usado
            dist = _size_bins(sizes)
            pd.DataFrame(dist, columns=["tamanho", "n_clusters"]).to_csv(
                os.path.join(args.out_dir, f"{name}_size_dist.tsv"),
                sep="\t", index=False)

            if have_plt:
                fig, ax = plt.subplots(figsize=(7, 4))
                labels = [d[0] for d in dist]
                vals = [d[1] for d in dist]
                ax.bar(labels, vals, color="#2ca02c")
                ax.set_xlabel("tamanho do cluster (n. de sequencias)")
                ax.set_ylabel("n. de clusters")
                ax.set_yscale("log")
                ax.set_title(f"{name} -- distribuicao de tamanho "
                            f"({len(sizes)} cluster(s))")
                fig.tight_layout()
                dist_png = os.path.join(args.out_dir, f"{name}_size_dist.png")
                fig.savefig(dist_png, dpi=150)
                plt.close(fig)

                # rank-size log-log: mesma logica da curva de abundancia --
                # mostra se a distribuicao segue algo tipo lei de potencia
                # (poucos clusters grandes, cauda longa de pequenos)
                ranked = sizes.sort_values(ascending=False).reset_index(drop=True)
                fig2, ax2 = plt.subplots(figsize=(6, 5))
                ax2.plot(range(1, len(ranked) + 1), ranked.values, "o",
                        markersize=3, alpha=0.6, color="#1f77b4")
                ax2.set_xscale("log")
                ax2.set_yscale("log")
                ax2.set_xlabel("rank do cluster")
                ax2.set_ylabel("tamanho do cluster")
                ax2.set_title(f"{name} -- rank x tamanho")
                fig2.tight_layout()
                rank_png = os.path.join(args.out_dir, f"{name}_rank_size.png")
                fig2.savefig(rank_png, dpi=150)
                plt.close(fig2)

            prog.mark_done(name)

        except Exception as e:
            print(f"[ERRO] {name}: {e}")
            prog.mark_failed(name, e)

    if maps:
        update_table(pd.concat(maps, ignore_index=True))
        print(f"[tabela] colunas cluster_rep/cluster_size/is_cluster_rep gravadas")

    rep.h("Reducao por unidade")
    rep.table(rows, ["unidade", "antes", "representantes", "identidade",
                     "cobertura", "maior cluster"])
    rep.p("\n> `cluster_size` esta na tabela central. Use-o em qualquer "
          "contagem por comunidade -- numero de representantes NAO e "
          "estimativa de abundancia.")
    if have_plt and maps:
        rep.p("\n> Distribuicao de tamanho e rank x tamanho por unidade em "
              "`*_size_dist.png` e `*_rank_size.png`. Rode `cluster --sweep` "
              "antes deste passo se ainda nao decidiu --identity/--coverage.")
    rep.save()
    prog.report()


# =============================================================================
# matrix -- all-vs-all
# =============================================================================

def stage_matrix(args):
    import rotifer.devel.beta.sequence as rdbs
    from rotifer.devel.alpha import mvroliveira as rdao

    fastas = sorted(glob.glob(os.path.join(args.hits_dir, "*.fasta")))
    if not fastas:
        sys.exit(f"[erro] nenhum .fasta em {args.hits_dir}")

    save_run_params("matrix", args, args.out_dir)
    prog = Progress("matrix", args.out_dir)
    rep = Report("matrix", args.out_dir, "Comparacao todos-contra-todos")
    rows = []

    for fa in fastas:
        name = os.path.splitext(os.path.basename(fa))[0]
        out = os.path.join(args.out_dir, f"{name}_hits.tsv")
        if prog.is_done(name) and os.path.exists(out) and not args.force:
            print(f"[skip] {name}")
            continue
        try:
            seqobj = rdbs.sequence(fa)
            n = len(seqobj.df)
            print(f"\n=== {name}: all-vs-all de {n} sequencia(s) ===")
            t0 = time.time()
            hits = rdao.all_vs_all_search(
                seqobj, tool=args.tool, cpu=args.cpu, evalue=args.evalue,
                keep_self_hits=False, progress=args.progress)
            hits.to_csv(out, sep="\t", index=False)
            print(f"[ok] {name}: {len(hits)} par(es) ({time.time() - t0:.0f}s)")
            rows.append((name, n, len(hits), args.tool, args.evalue))
            prog.mark_done(name)
        except Exception as e:
            print(f"[ERRO] {name}: {e}")
            prog.mark_failed(name, e)

    rep.h("Matrizes")
    rep.table(rows, ["unidade", "sequencias", "pares", "ferramenta", "e-value"])
    rep.p("\n> O e-value aqui e o TETO do que pode virar aresta. Quem escolhe "
          "o corte real e o estagio `ssn`. Se estiver estudando uma "
          "superfamilia divergente, use um valor permissivo aqui e deixe a "
          "varredura de corte trabalhar.")
    rep.save()
    prog.report()


# =============================================================================
# ssn -- rede + comunidades
#
# Duas correcoes importantes em relacao ao pipeline original:
#
# (a) ORFAS. O all-vs-all roda com keep_self_hits=False e o grafo nasce da
#     tabela de pares, entao sequencia sem nenhuma aresta simplesmente NAO
#     vira no -- some em silencio. Num fluxo query-centrico isso e inofensivo;
#     num estudo de diversidade as orfas sao potencialmente as familias mais
#     divergentes, exatamente o que se procura. Aqui o conjunto de nos e
#     reconciliado contra o fasta e as ausentes sao reportadas e marcadas na
#     tabela como 'orphan'.
#
# (b) FALLBACK. Quando o auto_cutoff degenera, o original criava nos sem
#     nenhuma aresta -- um graphml inutil. Aqui o fallback constroi o grafo
#     de verdade a partir dos pares, com corte fixo.
# =============================================================================

def _community_attr(G):
    for _, d in G.nodes(data=True):
        for a in ("Louvain", "Leiden", "leidein", "community"):
            if a in d:
                return a
    return None


def stage_ssn(args):
    import pandas as pd
    import networkx as nx
    from rotifer.devel.alpha import mvroliveira as rdao

    matrices = sorted(glob.glob(os.path.join(args.matrix_dir, "*_hits.tsv")))
    if not matrices:
        sys.exit(f"[erro] nenhuma matriz em {args.matrix_dir} (rode 'matrix')")

    save_run_params("ssn", args, args.out_dir)
    prog = Progress("ssn", args.out_dir)
    rep = Report("ssn", args.out_dir, "Rede de similaridade e comunidades")
    rows, assign = [], []

    for path in matrices:
        name = os.path.basename(path).replace("_hits.tsv", "")
        out = os.path.join(args.out_dir, f"{name}.graphml")
        if prog.is_done(name) and os.path.exists(out) and not args.force:
            print(f"[skip] {name}")
            continue

        try:
            hits = pd.read_csv(path, sep="\t", low_memory=False)
            n_before_range = len(hits)

            # restringe o INTERVALO de bitscore testado no closeness_scan,
            # por PERCENTIL (nao valor absoluto -- bitscore escala com
            # comprimento de alinhamento, e unidades diferentes tem tamanhos
            # bem diferentes entre si). Evita gastar passos de amostragem
            # numa faixa que nunca seria escolhida como pico: cortes frouxos
            # demais no fim de baixo, cortes que so isolam singletons no fim
            # de cima.
            if args.cutoff_min_pct is not None or args.cutoff_max_pct is not None:
                col_bit = next((c for c in hits.columns
                               if "bitscore" in c.lower() or c.lower() == "bits"), None)
                if col_bit:
                    lo = (hits[col_bit].quantile(args.cutoff_min_pct / 100)
                          if args.cutoff_min_pct is not None else hits[col_bit].min())
                    hi = (hits[col_bit].quantile(args.cutoff_max_pct / 100)
                          if args.cutoff_max_pct is not None else hits[col_bit].max())
                    hits = hits[(hits[col_bit] >= lo) & (hits[col_bit] <= hi)]
                    print(f"  [range] bitscore restrito a [{lo:.1f}, {hi:.1f}] "
                          f"(percentil {args.cutoff_min_pct}-{args.cutoff_max_pct}), "
                          f"{n_before_range} -> {len(hits)} pares")
                else:
                    print(f"  [aviso] coluna de bitscore nao encontrada, "
                          f"--cutoff-min-pct/--cutoff-max-pct ignorados")

            print(f"\n=== {name}: {len(hits)} par(es) na matriz ===")
            t0 = time.time()
            fallback = False

            try:
                G = rdao.build_ssn(hits, auto_cutoff=True,
                                   cutoff_steps=args.cutoff_steps,
                                   auto_start=getattr(args, "auto_start", True),
                                   add_community=True, progress=args.progress)
            except ValueError as e:
                if "auto_cutoff" not in str(e) and "closeness" not in str(e):
                    raise
                # N pequeno: a curva de closeness degenera e nao ha pico
                # confiavel para escolher automaticamente. Nesse regime,
                # subdividir em comunidades nao e significativo -- inventar
                # um corte fixo so pra ter alguma rede criaria estrutura que
                # os dados nao sustentam. Uma comunidade so, sem arestas
                # fabricadas.
                #
                # A curva em si (mesmo sem pico confiavel) ainda e calculada
                # e guardada -- serve de ponto de partida para inspecao
                # manual no dashboard (slider de corte), mesmo quando o
                # algoritmo nao confia em nenhum ponto dela como escolha
                # automatica.
                print(f"  [fallback] auto_cutoff degenerou -- comunidade unica "
                      f"(sem pico real na curva de closeness, N pequeno demais)")
                ids = set(hits.iloc[:, 0]).union(set(hits.iloc[:, 1]))
                G = nx.Graph()
                for i in ids:
                    G.add_node(str(i), Louvain="Louvain_0",
                              ssn_fallback="comunidade_unica")
                fallback = True
                try:
                    G.graph["closeness_scan"] = rdao.closeness_scan(
                        hits, n_steps=args.cutoff_steps,
                        auto_start=getattr(args, "auto_start", True),
                        progress=False)
                except Exception as e_scan:
                    print(f"  [aviso] nao consegui calcular a curva mesmo assim: {e_scan}")

            # --- reconciliacao com o fasta: quem ficou de fora ---
            fa = os.path.join(args.hits_dir, f"{name}.fasta")
            orphans = []
            if os.path.exists(fa):
                todos = set(read_fasta(fa).keys())
                nos = {str(n) for n in G.nodes}
                orphans = sorted(todos - nos)
                if orphans:
                    op = os.path.join(args.out_dir, f"{name}_orphans.txt")
                    with open(op, "w") as fh:
                        fh.write("\n".join(orphans) + "\n")
                    print(f"  [orfas] {len(orphans)}/{len(todos)} sem nenhuma "
                          f"aresta acima do corte -> {op}")
                    if args.keep_orphans:
                        for o in orphans:
                            G.add_node(o, Louvain="orphan", is_orphan=1)
                        print(f"  [orfas] adicionadas ao grafo como 'orphan'")

            if hasattr(rdao, "export_ssn"):
                rdao.export_ssn(G, out)
            else:
                nx.write_graphml(G, out)

            # PNG do closeness_scan -- pra inspecao visual de ONDE o
            # auto_cutoff escolheu o corte (ou, no caso de fallback, para
            # decidir manualmente pelo dashboard mesmo sem pico confiavel).
            plot_path = None
            scan = G.graph.get("closeness_scan")
            if scan is not None and hasattr(rdao, "plot_closeness_scan"):
                try:
                    plot_path = os.path.join(args.out_dir, f"{name}_closeness_scan.png")
                    rdao.plot_closeness_scan(scan, path=plot_path, title=name)
                    print(f"  [scan] curva de closeness -> {plot_path}")
                except Exception as e_plot:
                    print(f"  [aviso] falha ao gerar PNG do closeness_scan: {e_plot}")
                    plot_path = None

            # veredito automatico sobre a curva. Pega os dois modos de falha
            # que passam despercebidos ao olhar uma figura isolada: curva
            # ainda subindo no fim (otimo esta ALEM do que foi varrido) e
            # curva plana (nao ha pico, o corte sai quase arbitrario).
            diag = None
            if scan is not None and hasattr(rdao, "scan_diagnostics"):
                try:
                    diag = rdao.scan_diagnostics(scan)
                    flag = "  <-- ATENCAO" if (diag.get("truncated") or diag.get("flat")) else ""
                    print(f"  [scan] {diag['verdict']}{flag}")
                except Exception as e_diag:
                    print(f"  [aviso] scan_diagnostics falhou: {e_diag}")

            attr = _community_attr(G) or "Louvain"
            comms = Counter(str(d.get(attr, "NA")) for _, d in G.nodes(data=True))
            deg = dict(G.degree())
            assign.append(pd.DataFrame([{
                "id": str(nd),
                "community": str(d.get(attr, "NA")),
                "ssn_degree": deg.get(nd, 0),
                "is_orphan": bool(d.get("is_orphan", 0)) or str(nd) in set(orphans),
            } for nd, d in G.nodes(data=True)]))

            print(f"[ok] {name}: {G.number_of_nodes()} nos, "
                  f"{G.number_of_edges()} arestas, {len(comms)} comunidade(s) "
                  f"({time.time() - t0:.0f}s)")
            rows.append((name, G.number_of_nodes(), G.number_of_edges(),
                         len(comms), len(orphans), "sim" if fallback else "nao",
                         os.path.basename(plot_path) if plot_path else "-"))
            prog.mark_done(name)

        except Exception as e:
            print(f"[ERRO] {name}: {e}")
            prog.mark_failed(name, e)

    if assign:
        merged = pd.concat(assign, ignore_index=True)
        update_table(merged)
        top = (merged.groupby("community").size()
               .sort_values(ascending=False).head(20))
        rep.h("Comunidades (20 maiores, por numero de representantes)")
        rep.table([(c, int(n)) for c, n in top.items()],
                  ["comunidade", "representantes"])

    rep.h("Redes")
    rep.table(rows, ["unidade", "nos", "arestas", "comunidades", "orfas",
                     "fallback", "scan (png)"])
    rep.p("\n> Orfas sao sequencias sem nenhuma aresta acima do corte. Num "
          "estudo de diversidade elas nao sao lixo: sao as candidatas mais "
          "divergentes. Estao listadas em `*_orphans.txt` e marcadas na "
          "tabela central (`is_orphan`).")
    rep.p("\n> `fallback=sim` significa que a curva de closeness nao teve "
          "pico (N pequeno demais para o corte automatico fazer sentido) -- "
          "a unidade vira uma comunidade so, sem arestas fabricadas. Nao "
          "tenta simular uma subdivisao que os dados nao sustentam.")
    rep.p("\n> `*_closeness_scan.png` mostra a curva que o auto_cutoff "
          "escaneou e onde escolheu o corte -- use pra conferir se a "
          "escolha faz sentido, principalmente se `--cutoff-min-pct`/"
          "`--cutoff-max-pct` restringiram a faixa testada.")
    rep.save()
    prog.report()


# =============================================================================
# annotate -- taxonomia e arquitetura de dominio
#
# Best effort e degradacao explicita: tenta as fontes na ordem e reporta a
# cobertura obtida. Nunca falha o fluxo inteiro por causa de anotacao.
# =============================================================================

def stage_annotate(args):
    import pandas as pd

    save_run_params("annotate", args, args.out_dir)
    rep = Report("annotate", args.out_dir, "Anotacao (taxonomia e arquitetura)")
    df = load_table()
    if len(df) == 0:
        sys.exit("[erro] tabela central vazia (rode 'collect')")
    ids = df["id"].astype(str).tolist()

    if args.taxonomy:
        got = None
        # 1) organismo ja capturado do header do nr no collect
        if "organism" in df.columns and df["organism"].notna().any():
            print(f"[taxonomia] organismo do header: "
                  f"{df['organism'].notna().sum()}/{len(df)}")
        # 2) rotifer
        try:
            import rotifer.db.ncbi as ncbi
            tax = ncbi.taxonomy() if hasattr(ncbi, "taxonomy") else None
            if tax is not None and hasattr(tax, "fetchall"):
                got = tax.fetchall(ids)
        except Exception as e:
            print(f"[taxonomia] rotifer indisponivel ({type(e).__name__}: {e})")
        if got is not None and len(got):
            cols = [c for c in got.columns
                    if c.lower() in ("id", "taxid", "organism", "lineage",
                                     "superkingdom", "phylum", "class",
                                     "order", "family", "genus")]
            got = got[cols].rename(columns={cols[0]: "id"})
            update_table(got, overwrite=False)
            print(f"[taxonomia] {len(got)} registro(s) anotado(s)")
        else:
            print("[taxonomia] sem fonte externa -- usando so o header do banco")
            if "organism" in df.columns:
                # genero como proxy grosseiro de linhagem
                gen = df[["id", "organism"]].dropna().copy()
                gen["genus"] = gen["organism"].astype(str).str.split().str[0]
                update_table(gen[["id", "genus"]], overwrite=False)

    if args.arch:
        try:
            from rotifer.devel.alpha import rodolfo as rdar
            sub = df[["id"]].copy()
            sub["pid"] = sub["id"]
            out = rdar.add_arch_to_df(sub, db=args.db)
            keep = [c for c in ("id", "pfam", "aravind", "arch")
                    if c in out.columns]
            if len(keep) > 1:
                update_table(out[keep])
                n = out[keep[1]].notna().sum()
                print(f"[arquitetura] {n}/{len(out)} com anotacao")
        except Exception as e:
            print(f"[arquitetura] falhou ({type(e).__name__}: {e})")

    df = load_table()
    rep.h("Cobertura da tabela central apos anotacao")
    rep.code(table_columns_report(df))
    rep.save()


# =============================================================================
# context -- vizinhanca genomica
#
# Diferenca central em relacao ao original: o ESCOPO e argumento, nao regra.
#   --scope all          todas as sequencias  (estudo comparativo)
#   --scope community:X  so uma comunidade
#   --scope anchor       so a comunidade da ancora  (comportamento antigo)
#   --scope top:N        so as N maiores comunidades  (meio-termo pratico)
# =============================================================================

def _scope_ids(df, scope):
    if scope == "all" or "community" not in df.columns:
        return df["id"].astype(str).tolist(), "todas as sequencias"
    if scope == "anchor":
        if "is_anchor" not in df.columns or not df["is_anchor"].any():
            return df["id"].astype(str).tolist(), "sem ancora -- caiu para todas"
        comm = df.loc[df["is_anchor"] == True, "community"].iloc[0]
        sel = df[df["community"] == comm]
        return sel["id"].astype(str).tolist(), f"comunidade da ancora ({comm})"
    if scope.startswith("community:"):
        comm = scope.split(":", 1)[1]
        sel = df[df["community"].astype(str) == comm]
        return sel["id"].astype(str).tolist(), f"comunidade {comm}"
    if scope.startswith("top:"):
        n = int(scope.split(":", 1)[1])
        top = df["community"].value_counts().head(n).index
        sel = df[df["community"].isin(top)]
        return sel["id"].astype(str).tolist(), f"{n} maiores comunidades"
    raise SystemExit(f"[erro] --scope invalido: {scope}")


def stage_context(args):
    import pandas as pd
    import rotifer.db.ncbi as ncbi
    from rotifer.genome.data import NeighborhoodDF

    save_run_params("context", args, args.out_dir)
    prog = Progress("context", args.out_dir)
    rep = Report("context", args.out_dir, "Vizinhanca genomica")

    df = load_table()
    if len(df) == 0:
        sys.exit("[erro] tabela central vazia")
    ids, desc = _scope_ids(df, args.scope)
    ids = sorted(set(ids))
    print(f"[escopo] {desc}: {len(ids)} sequencia(s)")
    if len(ids) > args.warn_above:
        print(f"[aviso] escopo grande -- o cursor do NCBI vai levar tempo. "
              f"Considere --scope top:N para reduzir.")

    try:
        gnc = ncbi.GeneNeighborhoodCursor(before=args.window, after=args.window,
                                          threads=args.threads,
                                          progress=args.progress)
        ndf = gnc.fetchall(ids)
        if ndf is None or len(ndf) == 0:
            raise ValueError("nenhuma vizinhanca retornada")

        nb = ndf["block_id"].nunique() if "block_id" in ndf.columns else "?"
        print(f"[ok] {len(ndf)} feature(s), {nb} vizinhanca(s)")

        if args.arch:
            try:
                from rotifer.devel.alpha import rodolfo as rdar
                ndf = rdar.add_arch_to_df(ndf, db=args.db)
                n_pfam = ndf["pfam"].notna().sum() if "pfam" in ndf.columns else 0
                print(f"      [arquitetura] pfam em {n_pfam}/{len(ndf)} feature(s)")
            except Exception as e:
                print(f"      [arquitetura] falhou ({type(e).__name__}: {e}) "
                      f"-- seguindo so com 'product'")

        out = os.path.join(args.out_dir, "neighborhood.tsv")
        ndf.to_csv(out, sep="\t", index=False)

        ann = args.annotation if args.annotation in ndf.columns else "product"
        try:
            compact = NeighborhoodDF(ndf).compact_neighborhood(annotation=ann)
            if compact is not None:
                cpath = os.path.join(args.out_dir, "neighborhood_compact.tsv")
                compact.to_csv(cpath, sep="\t")
                print(f"      [compact] {len(compact)} bloco(s) -> {cpath}")
        except Exception as e:
            print(f"      [compact] falhou: {type(e).__name__}: {e}")

        if getattr(gnc, "missing", None) is not None and len(gnc.missing):
            gnc.missing.to_csv(os.path.join(args.out_dir, "missing.tsv"),
                               sep="\t", index=False)
            print(f"      {len(gnc.missing)}/{len(ids)} sem vizinhanca")

        rep.h("Escopo")
        rep.kv("selecao", desc)
        rep.kv("sequencias consultadas", len(ids))
        rep.kv("janela", f"{args.window} genes de cada lado")
        rep.kv("features recuperadas", len(ndf))
        rep.kv("vizinhancas", nb)
        prog.mark_done(args.scope)

    except Exception as e:
        print(f"[ERRO] {e}")
        prog.mark_failed(args.scope, e)

    rep.save()
    prog.report()


# =============================================================================
# profile -- descreve CADA comunidade
#
# Este e o estagio que substitui o 'select' do pipeline original na funcao de
# "olhar para as comunidades". A diferenca conceitual:
#
#   select (original) = FILTRO. Escolhe a comunidade da query, aceita/rejeita
#                       membros, descarta o resto. As outras comunidades sao
#                       ruido.
#   profile (aqui)    = DESCRICAO. Nao descarta nada. Perfila cada comunidade
#                       e produz a tabela comparativa. As outras comunidades
#                       sao o resultado.
#
# Dois modos:
#   COM --markers  pontua cada comunidade contra um vocabulario conhecido
#   SEM --markers  DESCOBERTA: reporta os termos enriquecidos em cada
#                  comunidade contra o fundo. Sem vocabulario, acha o
#                  vocabulario.
#
# E um bom lugar para parar.
# =============================================================================

def stage_profile(args):
    import pandas as pd

    save_run_params("profile", args, args.out_dir)
    rep = Report("profile", args.out_dir, "Perfil comparativo das comunidades")
    cfg = load_markers(args.markers)

    df = load_table()
    if len(df) == 0:
        sys.exit("[erro] tabela central vazia")
    if "community" not in df.columns:
        sys.exit("[erro] sem coluna 'community' (rode 'ssn')")

    # ---- contexto por sequencia, se houver ----
    ctx_terms = defaultdict(list)     # id -> termos da vizinhanca
    ctx_eval = {}                     # id -> evaluate_context(...)
    npath = os.path.join(args.context_dir, "neighborhood.tsv")
    if os.path.exists(npath):
        nb = pd.read_csv(npath, sep="\t", low_memory=False)
        # colunas de anotacao. Por padrao usa UMA fonte de termos, nao todas:
        # 'pfam' e 'product' costumam dizer a mesma coisa com palavras
        # diferentes ("PilT" e "PilT family protein"), e somar as duas conta o
        # mesmo evento duas vezes, inflando o ranking de enriquecimento.
        # A deteccao de MARCADORES continua olhando todas as colunas, porque
        # ali redundancia so ajuda.
        all_cols = [c for c in ("pfam", "product", "gene") if c in nb.columns]
        if args.term_cols:
            cols = [c for c in args.term_cols if c in nb.columns]
        else:
            cols = ([c for c in ("pfam", "arch") if c in nb.columns]
                    or [c for c in ("product", "gene") if c in nb.columns])[:1]
        if {"block_id", "pid"}.issubset(nb.columns) and all_cols:
            print(f"[contexto] {npath} -- termos de: {', '.join(cols) or 'nenhuma'} "
                  f"| marcadores de: {', '.join(all_cols)}")
            wanted = set(df["id"].astype(str))
            for _, block in nb.groupby("block_id"):
                marker_texts = []
                for c in all_cols:
                    marker_texts.extend(block[c].dropna().astype(str).tolist())
                term_texts = []
                for c in cols:
                    term_texts.extend(block[c].dropna().astype(str).tolist())
                terms = [t for v in term_texts for t in split_terms(v)]
                ev = evaluate_context(marker_texts, cfg)
                for pid in set(block["pid"].dropna().astype(str)) & wanted:
                    ctx_terms[pid].extend(terms)
                    if ev is not None:
                        prev = ctx_eval.get(pid)
                        if prev is None or ev["n_markers"] > prev["n_markers"]:
                            ctx_eval[pid] = ev
    else:
        print(f"[contexto] {npath} ausente -- perfil sem vizinhanca genomica")

    # genero derivado do organismo, se a taxonomia nao foi anotada. Organismo
    # cru quase nunca agrega ("Genus sp. XYZ" e unico por sequencia) e o
    # resumo taxonomico sai inutil; o genero ja separa.
    if "genus" not in df.columns and "organism" in df.columns:
        df["genus"] = df["organism"].dropna().astype(str).str.split().str[0]

    # ---- perfil por comunidade ----
    has_size = "cluster_size" in df.columns
    rows, vocab_input = [], {}
    for comm, sub in df.groupby("community"):
        ids = sub["id"].astype(str).tolist()
        n_rep = len(sub)
        n_seq = int(sub["cluster_size"].fillna(1).sum()) if has_size else n_rep

        lens = sub["length"].dropna() if "length" in sub.columns else pd.Series(dtype=float)
        tlens = (sub["target_length"].dropna()
                if "target_length" in sub.columns else pd.Series(dtype=float))
        tcovs = (sub["target_coverage"].dropna()
                if "target_coverage" in sub.columns else pd.Series(dtype=float))
        taxcol = next((c for c in ("phylum", "class", "order", "genus",
                                   "organism") if c in sub.columns), None)
        taxtop = ""
        if taxcol:
            vc = sub[taxcol].dropna().value_counts()
            taxtop = "; ".join(f"{k} ({v})" for k, v in vc.head(3).items())

        archcol = next((c for c in ("arch", "pfam") if c in sub.columns), None)
        archtop = ""
        if archcol:
            vc = sub[archcol].dropna().value_counts()
            archtop = "; ".join(f"{k} ({v})" for k, v in vc.head(3).items())

        with_ctx = sum(1 for i in ids if i in ctx_terms)
        vocab_input[comm] = [t for i in ids for t in ctx_terms.get(i, [])]

        row = {
            "community": comm,
            "n_representantes": n_rep,
            "n_sequencias": n_seq,
            "len_mediana": float(lens.median()) if len(lens) else None,
            "len_min": float(lens.min()) if len(lens) else None,
            "len_max": float(lens.max()) if len(lens) else None,
            "target_len_mediana": float(tlens.median()) if len(tlens) else None,
            "target_coverage_mediana": float(tcovs.median()) if len(tcovs) else None,
            "grau_medio": (float(sub["ssn_degree"].mean())
                           if "ssn_degree" in sub.columns else None),
            "com_contexto": with_ctx,
            "taxonomia_top3": taxtop,
            "arquitetura_top3": archtop,
            "tem_ancora": (bool(sub["is_anchor"].any())
                           if "is_anchor" in sub.columns else False),
        }

        if cfg is not None:
            evs = [ctx_eval[i] for i in ids if i in ctx_eval]
            if evs:
                mark = Counter(m for e in evs for m in e["markers"])
                row["marcadores_top"] = "; ".join(f"{k} ({v})"
                                                  for k, v in mark.most_common(5))
                row["frac_com_marcador"] = sum(1 for e in evs
                                               if e["n_markers"] > 0) / len(evs)
                row["frac_com_veto"] = sum(1 for e in evs
                                           if e["has_veto"]) / len(evs)
                row["media_marcadores"] = sum(e["n_markers"]
                                              for e in evs) / len(evs)
        rows.append(row)

    prof = pd.DataFrame(rows).sort_values("n_sequencias", ascending=False)
    ppath = os.path.join(args.out_dir, "communities.tsv")
    os.makedirs(args.out_dir, exist_ok=True)
    prof.to_csv(ppath, sep="\t", index=False)
    print(f"\n[ok] perfil de {len(prof)} comunidade(s) -> {ppath}")

    # ---- descoberta de vocabulario ----
    disc = {}
    if any(vocab_input.values()):
        disc = discover_vocabulary(vocab_input, top=args.top_terms,
                                   min_count=args.min_term_count)
        dpath = os.path.join(args.out_dir, "enriched_terms.tsv")
        with open(dpath, "w") as fh:
            fh.write("community\tterm\tn\tfreq_comunidade\tfreq_fundo\tlog2_odds\n")
            for g, terms in disc.items():
                for t, n, fg, fb, lo in terms:
                    fh.write(f"{g}\t{t}\t{n}\t{fg:.5f}\t{fb:.5f}\t{lo:.3f}\n")
        print(f"[ok] termos enriquecidos -> {dpath}")

    # ---- relatorio ----
    rep.h("Modo")
    rep.kv("vocabulario", cfg["name"] if cfg else
           "DESCOBERTA (sem --markers): termos enriquecidos contra o fundo")
    rep.kv("comunidades", len(prof))
    rep.kv("sequencias (representantes)", int(prof["n_representantes"].sum()))
    if has_size:
        rep.kv("sequencias (antes do cluster)", int(prof["n_sequencias"].sum()))

    show = prof.head(args.top_communities)
    rep.h(f"Comunidades ({len(show)} maiores)")
    headers = ["comunidade", "repr.", "seqs", "len med.", "grau med.",
               "ctx", "taxonomia (top 3)"]
    rep.table([(r["community"], r["n_representantes"], r["n_sequencias"],
                None if r["len_mediana"] is None else f"{r['len_mediana']:.0f}",
                None if r["grau_medio"] is None else f"{r['grau_medio']:.1f}",
                r["com_contexto"], r["taxonomia_top3"] or "-")
               for _, r in show.iterrows()], headers)

    if cfg is not None and "frac_com_marcador" in prof.columns:
        rep.h(f"Marcadores do vocabulario `{cfg['name']}`")
        rep.table([(r["community"],
                    f"{r.get('frac_com_marcador', 0):.0%}",
                    f"{r.get('media_marcadores', 0):.1f}",
                    f"{r.get('frac_com_veto', 0):.0%}",
                    r.get("marcadores_top", "") or "-")
                   for _, r in show.iterrows()],
                  ["comunidade", "% com marcador", "media de marcadores",
                   "% com veto", "marcadores mais frequentes"])

    if disc:
        rep.h("Termos enriquecidos por comunidade (descoberta)")
        rep.p("Log2 da razao entre a frequencia do termo na comunidade e no "
              "conjunto todo. Ranking exploratorio, nao teste estatistico.")
        for comm in show["community"].tolist():
            terms = disc.get(comm) or []
            if not terms:
                continue
            rep.h(f"{comm}", level=3)
            rep.table([(t, n, f"{lo:+.2f}") for t, n, fg, fb, lo in terms],
                      ["termo", "n", "log2 odds"])

    rep.h("Como ler isto")
    rep.p("- `n_sequencias` usa `cluster_size`: e o tamanho real da familia, "
          "nao o numero de representantes.")
    rep.p("- Comunidade com grau medio alto e len mediana estreita = familia "
          "coesa. Grau baixo e len dispersa = provavelmente heterogenea, "
          "candidata a re-analise com outro corte na SSN.")
    rep.p("- `com_contexto` baixo significa que a comparacao de vizinhanca "
          "daquela comunidade esta apoiada em pouca evidencia.")
    rep.p("- `target_coverage_mediana` baixa (longe de 1.0) indica que a "
          "sequencia coletada e uma fracao pequena da proteina de origem -- "
          "comunidade provavelmente formada por um dominio isolado dentro "
          "de proteinas multi-dominio maiores, nao pela proteina inteira. "
          "Compare com `target_len_mediana`: mesma `len_mediana` do "
          "envelope com `target_len_mediana` muito maior e o sinal disso.")
    rep.save()


# =============================================================================
# select -- filtro OPCIONAL
#
# Existe para quem quer o comportamento do pipeline antigo: reduzir a um
# subconjunto curado antes de construir modelo. Diferente do original, as
# regras vem do JSON de marcadores, nao do codigo, e nao ha exigencia de
# ancora.
# =============================================================================

def stage_select(args):
    import pandas as pd

    save_run_params("select", args, args.out_dir)
    rep = Report("select", args.out_dir, "Selecao curada")
    cfg = load_markers(args.markers)

    df = load_table()
    if len(df) == 0:
        sys.exit("[erro] tabela central vazia")
    if "community" not in df.columns:
        sys.exit("[erro] sem coluna 'community' (rode 'ssn')")

    ctx_eval = {}
    npath = os.path.join(args.context_dir, "neighborhood.tsv")
    if os.path.exists(npath) and cfg is not None:
        nb = pd.read_csv(npath, sep="\t", low_memory=False)
        cols = [c for c in ("pfam", "product", "gene") if c in nb.columns]
        if {"block_id", "pid"}.issubset(nb.columns) and cols:
            wanted = set(df["id"].astype(str))
            for _, block in nb.groupby("block_id"):
                texts = []
                for c in cols:
                    texts.extend(block[c].dropna().astype(str).tolist())
                ev = evaluate_context(texts, cfg)
                for pid in set(block["pid"].dropna().astype(str)) & wanted:
                    prev = ctx_eval.get(pid)
                    if prev is None or ev["n_markers"] > prev["n_markers"]:
                        ctx_eval[pid] = ev

    # quais comunidades entram
    if args.communities:
        keep = set(args.communities)
    elif args.anchor_community and "is_anchor" in df.columns and df["is_anchor"].any():
        keep = set(df.loc[df["is_anchor"] == True, "community"].astype(str))
    else:
        keep = set(df["community"].astype(str))
    print(f"[selecao] {len(keep)} comunidade(s) consideradas")

    seqs = {}
    for fa in glob.glob(os.path.join(args.hits_dir, "*.fasta")):
        seqs.update(read_fasta(fa))

    decisoes, por_comm = [], defaultdict(list)
    dist = Counter()
    for _, r in df.iterrows():
        sid, comm = str(r["id"]), str(r["community"])
        if comm not in keep:
            continue
        ev = ctx_eval.get(sid)

        if ev is None:
            verd = "ACEITO" if args.no_context == "accept" else "REJEITADO"
            motivo = f"sem contexto (--no-context={args.no_context})"
            marc, nm = "", ""
        elif ev["has_veto"]:
            verd, motivo = "REJEITADO", f"veto: {','.join(ev['veto'])}"
            marc, nm = ",".join(ev["markers"]), ev["n_markers"]
        elif ev["n_markers"] >= args.min_markers:
            verd = "ACEITO"
            motivo = f"{ev['n_markers']} marcador(es): {','.join(ev['markers']) or '-'}"
            marc, nm = ",".join(ev["markers"]), ev["n_markers"]
            dist[ev["n_markers"]] += 1
        else:
            verd = "REJEITADO"
            motivo = f"so {ev['n_markers']} marcador(es) (--min-markers {args.min_markers})"
            marc, nm = ",".join(ev["markers"]), ev["n_markers"]
            dist[ev["n_markers"]] += 1

        decisoes.append({"id": sid, "community": comm, "veredito": verd,
                         "motivo": motivo, "marcadores": marc, "n_marcadores": nm})
        if verd == "ACEITO" and sid in seqs:
            por_comm[comm].append((sid, seqs[sid]))

    os.makedirs(args.out_dir, exist_ok=True)
    pd.DataFrame(decisoes).to_csv(
        os.path.join(args.out_dir, "decisoes.tsv"), sep="\t", index=False)

    rows = []
    for comm, recs in sorted(por_comm.items(), key=lambda kv: -len(kv[1])):
        if len(recs) < args.min_seqs:
            print(f"[skip] {comm}: {len(recs)} seq(s), abaixo de --min-seqs")
            rows.append((comm, len(recs), "descartada (poucas seqs)"))
            continue
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", comm)
        p = os.path.join(args.out_dir, f"{safe}_curado.fasta")
        write_fasta(p, recs)
        print(f"[ok] {comm}: {len(recs)} seq(s) -> {p}")
        rows.append((comm, len(recs), os.path.basename(p)))

    rep.h("Criterios")
    rep.kv("vocabulario", cfg["name"] if cfg else "nenhum (sem filtro de contexto)")
    rep.kv("minimo de marcadores", args.min_markers)
    rep.kv("sem contexto", args.no_context)
    rep.kv("comunidades consideradas", len(keep))
    if dist:
        rep.h("Distribuicao de marcadores na vizinhanca")
        rep.table(sorted(dist.items()), ["n marcadores", "sequencias"])
    rep.h("Fastas curados")
    rep.table(rows, ["comunidade", "sequencias", "arquivo"])
    rep.save()


# =============================================================================
# build -- alinha + hmmbuild, um modelo por grupo
# =============================================================================

def stage_build(args):
    import rotifer.devel.beta.sequence as rdbs

    if not which("hmmbuild"):
        sys.exit("[erro] hmmbuild nao encontrado no PATH")

    fastas = sorted(glob.glob(os.path.join(args.curated_dir, "*_curado.fasta")))
    if not fastas:
        sys.exit(f"[erro] nenhum *_curado.fasta em {args.curated_dir} (rode 'select')")

    save_run_params("build", args, args.out_dir)
    prog = Progress("build", args.out_dir)
    rep = Report("build", args.out_dir, "Construcao dos modelos")
    built, rows = [], []

    for fa in fastas:
        name = os.path.basename(fa).replace("_curado.fasta", "")
        hmm = os.path.join(args.out_dir, f"{name}.hmm")
        if prog.is_done(name) and os.path.exists(hmm) and not args.force:
            built.append(hmm)
            continue
        try:
            n = len(read_fasta(fa))
            if n < args.min_seqs:
                print(f"[skip] {name}: {n} seq(s), abaixo de --min-seqs "
                      f"(perfil de poucas sequencias nao estima bem as "
                      f"probabilidades)")
                rows.append((name, n, "descartado"))
                continue
            seqobj = rdbs.sequence(fa)
            aln = seqobj.align(method=args.aligner, cpu=args.cpu)
            aln_path = os.path.join(args.out_dir, f"{name}_aln.fasta")
            aln.to_file(aln_path)
            r = subprocess.run(["hmmbuild", "-n", name, hmm, aln_path],
                               capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError(f"hmmbuild: {r.stderr[-400:]}")
            print(f"[ok] {name}: {n} seq(s) -> {hmm}")
            rows.append((name, n, os.path.basename(hmm)))
            built.append(hmm)
            prog.mark_done(name)
        except Exception as e:
            print(f"[ERRO] {name}: {e}")
            prog.mark_failed(name, e)

    if built and args.concat:
        combined = os.path.join(args.out_dir, f"{args.concat_name}.hmm")
        with open(combined, "w") as dst:
            for h in built:
                with open(h) as src:
                    dst.write(src.read())
        subprocess.run(["hmmpress", "-f", combined], capture_output=True)
        print(f"\n[ok] banco combinado: {combined} ({len(built)} modelos)")

    rep.h("Modelos")
    rep.kv("alinhador", args.aligner)
    rep.table(rows, ["grupo", "sequencias", "modelo"])
    rep.save()
    prog.report()


# =============================================================================
# scan -- modelos contra alvo, + QC de separabilidade
#
# O auto-scan e o controle que o pipeline original nao tinha: escaneia cada
# modelo contra o proprio conjunto de partida e cruza o melhor hit com a
# comunidade de origem. Se os modelos nao separam as comunidades, as
# comunidades nao eram familias -- e melhor descobrir isso aqui do que na
# banca.
# =============================================================================

def stage_scan(args):
    import pandas as pd

    if not which("hmmsearch"):
        sys.exit("[erro] hmmsearch nao encontrado no PATH")

    hmms = [h for h in sorted(glob.glob(os.path.join(args.hmm_dir, "*.hmm")))
            if not os.path.basename(h).startswith(args.concat_name)]
    if not hmms:
        sys.exit(f"[erro] nenhum .hmm em {args.hmm_dir} (rode 'build')")

    save_run_params("scan", args, args.out_dir)
    prog = Progress("scan", args.out_dir)
    rep = Report("scan", args.out_dir, "Busca dos modelos")
    os.makedirs(args.out_dir, exist_ok=True)

    targets = []
    if args.target:
        if not os.path.exists(args.target):
            sys.exit(f"[erro] alvo nao encontrado: {args.target}")
        targets.append(("alvo", args.target))
    if args.self_scan:
        pool = os.path.join(args.out_dir, "_selfpool.fasta")
        recs = {}
        for fa in sorted(glob.glob(os.path.join(args.hits_dir, "*.fasta"))):
            recs.update(read_fasta(fa))
        write_fasta(pool, list(recs.items()))
        targets.append(("self", pool))
    if not targets:
        sys.exit("[erro] use --target e/ou --self-scan")

    print(f"[info] {len(hmms)} modelo(s) contra {len(targets)} alvo(s)")
    resumo, best = [], {}

    for tag, tgt in targets:
        for hmm in hmms:
            name = os.path.splitext(os.path.basename(hmm))[0]
            item = f"{tag}:{name}"
            tbl = os.path.join(args.out_dir, f"{tag}__{name}.tbl")
            if prog.is_done(item) and os.path.exists(tbl) and not args.force:
                continue
            try:
                r = subprocess.run(
                    ["hmmsearch", "--cpu", str(args.cpu), "-E", str(args.evalue),
                     "--tblout", tbl, "-o", "/dev/null", hmm, tgt],
                    capture_output=True, text=True)
                if r.returncode != 0:
                    raise RuntimeError(r.stderr[-400:])
                n_hits = 0
                with open(tbl) as fh:
                    for line in fh:
                        if line.startswith("#") or not line.strip():
                            continue
                        n_hits += 1
                        if tag == "self":
                            f = line.split()
                            sid, sc = f[0], float(f[5])
                            if sid not in best or sc > best[sid][1]:
                                best[sid] = (name, sc)
                print(f"[ok] {tag}/{name}: {n_hits} hit(s)")
                resumo.append({"alvo": tag, "modelo": name, "n_hits": n_hits})
                prog.mark_done(item)
            except Exception as e:
                print(f"[ERRO] {item}: {e}")
                prog.mark_failed(item, e)

    if resumo:
        d = pd.DataFrame(resumo).sort_values(["alvo", "n_hits"], ascending=[True, False])
        d.to_csv(os.path.join(args.out_dir, "resumo_scan.tsv"), sep="\t", index=False)
        rep.h("Hits por modelo")
        rep.table(d.values.tolist(), list(d.columns))

    if best:
        df = load_table()
        if "community" in df.columns:
            m = df[["id", "community"]].copy()
            m["id"] = m["id"].astype(str)
            m["melhor_modelo"] = m["id"].map(lambda i: best.get(i, (None, None))[0])
            m = m.dropna(subset=["melhor_modelo"])
            cm = pd.crosstab(m["community"], m["melhor_modelo"])
            cpath = os.path.join(args.out_dir, "confusao_modelo_x_comunidade.tsv")
            cm.to_csv(cpath, sep="\t")
            diag = sum(cm.loc[c, c] for c in cm.index if c in cm.columns)
            acc = diag / cm.values.sum() if cm.values.sum() else 0
            print(f"\n[QC] separabilidade: {acc:.1%} das sequencias tem como "
                  f"melhor modelo o da propria comunidade -> {cpath}")
            rep.h("QC de separabilidade")
            rep.kv("concordancia modelo x comunidade", f"{acc:.1%}")
            rep.p("\n> Concordancia baixa significa que os modelos nao "
                  "distinguem as comunidades entre si. Nesse caso as "
                  "comunidades provavelmente nao correspondem a familias "
                  "separaveis por perfil, e vale revisitar o corte da SSN.")

    rep.save()
    prog.report()


# =============================================================================
# table / report / status
# =============================================================================

def stage_table(args):
    df = load_table()
    if len(df) == 0:
        sys.exit("[erro] tabela central vazia")
    print(f"{len(df)} linha(s), {len(df.columns)} coluna(s)\n")
    print(table_columns_report(df))
    if args.out:
        df.to_csv(args.out, sep="\t", index=False)
        print(f"\n[ok] copia -> {args.out}")


def stage_report(args):
    """Monta Metodos + Resultados a partir dos checkpoints e dos REPORT.md."""
    out = []
    out.append(f"# Relatorio do fluxo\n\n_gerado em "
               f"{time.strftime('%Y-%m-%d %H:%M:%S')}_\n")

    out.append("\n## Metodos (parametros efetivamente usados)\n")
    for stage, d in _STAGE_DIRS:
        p = os.path.join(d, f".{stage}_params.json")
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            hist = json.load(fh)
        out.append(f"\n### `{stage}`\n")
        for h in hist:
            args_str = " ".join(f"--{k.replace('_', '-')} {v}"
                                for k, v in sorted(h["params"].items())
                                if v not in (None, False, "")
                                and not k.startswith("_"))
            out.append(f"- {h['when']}: `{args_str}`")

    out.append("\n## Resultados por estagio\n")
    for stage, d in _STAGE_DIRS:
        r = os.path.join(d, "REPORT.md")
        if not os.path.exists(r):
            continue
        with open(r) as fh:
            body = fh.read()
        # rebaixa os titulos: o REPORT.md de cada estagio e um documento
        # autonomo (comeca em '#'), mas aqui ele vira uma SUBSECAO. Sem isso
        # os '##' do estagio ficam no mesmo nivel dos '##' do consolidado e a
        # hierarquia quebra.
        demoted = []
        for line in body.split("\n"):
            if line.startswith("# "):
                continue
            if line.startswith("#"):
                line = "##" + line
            demoted.append(line)
        out.append(f"\n### Estagio `{stage}`\n")
        out.append("\n".join(demoted))

    try:
        df = load_table()
        if len(df):
            out.append("\n## Tabela central\n")
            out.append(f"{len(df)} sequencia(s), {len(df.columns)} coluna(s).\n")
            out.append("```\n" + table_columns_report(df) + "\n```")
    except Exception:
        pass

    with open(args.out, "w") as fh:
        fh.write("\n".join(out) + "\n")
    print(f"[ok] {args.out}")


_STAGE_DIRS = [
    ("seed", DIRS["seed"]), ("search", DIRS["search"]),
    ("inspect", DIRS["inspect"]),
    ("collect", DIRS["hits"]),
    ("cluster_sweep", os.path.join(DIRS["clustered"], "_sweep")),
    ("cluster", DIRS["clustered"]),
    ("matrix", DIRS["matrix"]), ("ssn", DIRS["ssn"]),
    ("annotate", DIRS["annot"]), ("context", DIRS["context"]),
    ("profile", DIRS["profile"]), ("select", DIRS["curated"]),
    ("build", DIRS["models"]), ("scan", DIRS["scan"]),
]



# =============================================================================
# dashboard -- painel HTML interativo, autocontido
#
# A geracao vive no modulo (rdao.ssn_dashboard), nao aqui: qualquer pipeline
# que produza .graphml no mesmo formato ganha o painel de graca, e a evolucao
# do painel acontece num lugar so.
# =============================================================================

def stage_dashboard(args):
    from rotifer.devel.alpha import mvroliveira as rdao

    save_run_params("dashboard", args, os.path.dirname(os.path.abspath(args.out)) or ".")

    # aliases: quando a sequencia de referencia entrou na rede sob outro
    # accession (o 'seed' registra isso na tabela central quando aplicavel),
    # o painel precisa saber para achar o no de referencia.
    aliases = {}
    if os.path.exists(TABLE):
        try:
            df = load_table()
            if {"unit", "id"}.issubset(df.columns) and "is_reference" in df.columns:
                ref = df[df["is_reference"].astype(str).str.lower().isin(("true", "1", "yes"))]
                for u, sub in ref.groupby("unit"):
                    aliases[str(u)] = [str(x) for x in sub["id"].dropna().unique()]
        except Exception as e:
            print(f"[aviso] nao consegui ler aliases da tabela central: {e}")

    out = rdao.ssn_dashboard(
        ssn_dir=args.ssn_dir,
        out=args.out,
        matrix_dir=args.matrix_dir,
        aliases=aliases or None,
        max_nodes=args.max_nodes,
        layout_iter=args.layout_iter,
    )
    print(f"[fim] {out} -- abra direto no navegador, sem servidor nem internet")
    if aliases:
        print(f"[info] {len(aliases)} unidade(s) com accession de referencia vindo da tabela")


def stage_status(args):
    print("=" * 70)
    print("PROGRESSO DO FLUXO")
    print("=" * 70)
    for stage, d in _STAGE_DIRS:
        cp = os.path.join(d, f".{stage}_checkpoint.json")
        rp = os.path.join(d, "REPORT.md")
        if not os.path.exists(cp) and not os.path.exists(rp):
            print(f"  {stage:10} nao iniciado")
            continue
        if not os.path.exists(cp):
            print(f"  {stage:10} concluido (sem checkpoint por item)")
            continue
        with open(cp) as fh:
            data = json.load(fh)
        n_done, n_fail = len(data["done"]), len(data["failed"])
        s = f"{n_done} ok" + (f", {n_fail} COM ERRO" if n_fail else "")
        print(f"  {stage:10} {s:30} (ultimo: {data.get('started', '?')})")
        if n_fail and args.verbose:
            for item, info in data["failed"].items():
                print(f"      {item}: {info['error']}")
    try:
        df = load_table()
        print(f"\n  tabela central: {len(df)} linha(s), {len(df.columns)} coluna(s)")
    except Exception:
        pass
    if not args.verbose:
        print("\n(use --verbose para ver os erros de cada estagio)")


# =============================================================================
# CLI
# =============================================================================

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="stage", required=True)

    p = sub.add_parser("seed", help="normaliza a entrada")
    p.add_argument("--input", help="fasta semente (uma ou varias sequencias)")
    p.add_argument("--hmm", action="append", help="HMM ja construido (repetivel)")
    p.add_argument("--unit", help="nome da unidade (default: nome do arquivo)")
    p.add_argument("--chopping", help="'1-169_323-361,171-318' -- fatia em dominios")
    p.add_argument("--out-dir", default=DIRS["seed"])
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=stage_seed)

    p = sub.add_parser("search", help="hmmbuild + hmmsearch contra o banco")
    p.add_argument("--seed-dir", default=DIRS["seed"])
    p.add_argument("--db", default=NR)
    p.add_argument("--out-dir", default=DIRS["search"])
    p.add_argument("--workers", type=int, default=4,
                   help="unidades processadas SIMULTANEAMENTE (processos)")
    p.add_argument("--cpus-per-worker", type=int, default=6,
                   help="threads do pyhmmer DENTRO de cada busca")
    p.add_argument("--evalue", type=float, default=None,
                   help="e-value de REPORTE. Default: o do HMMER (~10), "
                        "permissivo de proposito -- filtre depois no 'collect', "
                        "que nao exige refazer a busca")
    p.add_argument("--inc-evalue", type=float, default=None,
                   help="e-value de INCLUSAO (incE/incdomE)")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=stage_search)

    p = sub.add_parser("inspect", help="score x cobertura por hit, antes de cortar")
    p.add_argument("--search-dir", default=DIRS["search"])
    p.add_argument("--seed-dir", default=DIRS["seed"])
    p.add_argument("--out-dir", default=DIRS["inspect"])
    p.add_argument("--sort-by", choices=["evalue", "score"], default="evalue",
                   help="ordem de desenho, do pior para o melhor -- o melhor "
                        "fica por cima. Default: evalue")
    p.add_argument("--min-coverage", type=float, default=None,
                   help="so desenha a linha de referencia (nao filtra nada)")
    p.add_argument("--score-cutoff", type=float, default=None,
                   help="so desenha a linha de referencia (nao filtra nada)")
    p.add_argument("--dpi", type=int, default=200)
    p.set_defaults(func=stage_inspect)

    p = sub.add_parser("collect", help="recupera os hits e inicia a tabela")
    p.add_argument("--from-fasta", help="atalho: entra direto por um fasta pronto")
    p.add_argument("--unit", help="nome da unidade quando usar --from-fasta")
    p.add_argument("--seed-dir", default=DIRS["seed"])
    p.add_argument("--search-dir", default=DIRS["search"])
    p.add_argument("--db", default=NR)
    p.add_argument("--out-dir", default=DIRS["hits"])
    p.add_argument("--mode", choices=["envelope", "full"], default="envelope",
                   help="'envelope' recorta so a regiao alinhada (a rede compara "
                        "o mesmo fold); 'full' guarda a proteina inteira "
                        "(necessario se as familias diferirem em ARQUITETURA)")
    p.add_argument("--expand", type=int, default=10,
                   help="folga em cada borda do envelope")
    p.add_argument("--target-length", choices=["auto", "skip"], default="auto",
                   help="em --mode envelope, busca o tamanho da proteina "
                        "INTEIRA (nao so do envelope) para as colunas "
                        "target_length/target_coverage da tabela central. "
                        "Custa um fetch a mais no banco por unidade -- use "
                        "'skip' se a escala nao compensar. Sem efeito em "
                        "--mode full, onde e trivial")
    p.add_argument("--max-evalue", type=float, default=None)
    p.add_argument("--min-coverage", type=float, default=None,
                   help="cobertura minima DO PERFIL (nao do alvo). 0.6-0.7 "
                        "remove matches de motivo isolado. Default: sem filtro")
    p.add_argument("--max-seqs", type=int, default=None,
                   help="teto por unidade (o all-vs-all cresce com o quadrado)")
    p.add_argument("--anchor", action="store_true", default=True,
                   help="garante a semente no conjunto, identificavel")
    p.add_argument("--no-anchor", dest="anchor", action="store_false")
    p.add_argument("--anchor-identity", type=float, default=0.95)
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=stage_collect)

    p = sub.add_parser("cluster", help="reduz redundancia, preservando o mapa")
    p.add_argument("--hits-dir", default=DIRS["hits"])
    p.add_argument("--out-dir", default=DIRS["clustered"])
    p.add_argument("--identity", type=float, default=0.7)
    p.add_argument("--coverage", type=float, default=0.8)
    p.add_argument("--override", action="append", default=None,
                   help="nome:identity:coverage (repetivel)")
    p.add_argument("--sweep", action="store_true",
                   help="NAO clusteriza -- roda uma grade de identity x "
                        "coverage e planta quantos clusters cada combinacao "
                        "produz (heatmap). Cada celula e uma chamada de "
                        "mmseqs inteira; grade default pequena de proposito")
    p.add_argument("--sweep-identity", type=float, nargs="+",
                   default=[0.3, 0.5, 0.7, 0.9])
    p.add_argument("--sweep-coverage", type=float, nargs="+",
                   default=[0.5, 0.7, 0.9])
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=stage_cluster)

    p = sub.add_parser("matrix", help="all-vs-all")
    p.add_argument("--hits-dir", default=DIRS["clustered"])
    p.add_argument("--out-dir", default=DIRS["matrix"])
    p.add_argument("--tool", default="diamond", choices=["diamond", "blast", "mmseqs"])
    p.add_argument("--cpu", type=int, default=24)
    p.add_argument("--evalue", type=float, default=1e-3,
                   help="TETO do que pode virar aresta. Para superfamilia "
                        "divergente, use permissivo e deixe a SSN cortar")
    p.add_argument("--progress", action="store_true", default=True)
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=stage_matrix)

    p = sub.add_parser("ssn", help="rede + comunidades")
    p.add_argument("--matrix-dir", default=DIRS["matrix"])
    p.add_argument("--hits-dir", default=DIRS["clustered"],
                   help="usado para reconciliar orfas contra o grafo")
    p.add_argument("--out-dir", default=DIRS["ssn"])
    p.add_argument("--cutoff-steps", type=int, default=30)
    p.add_argument("--no-auto-start", dest="auto_start", action="store_false",
                   help="varre o intervalo bruto inteiro de bitscore. Por padrao "
                        "a varredura comeca no menor corte que de fato fragmenta "
                        "a rede (achado por bisseccao): abaixo disso a rede ainda "
                        "e um componente unico, entao nenhum ponto ali pode ser o "
                        "pico -- gastar amostragem la tira resolucao da regiao que "
                        "decide o resultado. Use esta flag para plotar a curva "
                        "completa, incluindo a cabeca plana.")
    p.add_argument("--cutoff-min-pct", type=float, default=None,
                   help="[EVITE] restringe o closeness_scan a partir deste "
                        "PERCENTIL do bitscore. Diferente de --no-auto-start, "
                        "isto REMOVE pares da matriz antes da varredura, entao "
                        "pode descartar o pico bom junto. Prefira o auto_start, "
                        "que apenas escolhe onde comecar a amostrar.")
    p.add_argument("--cutoff-max-pct", type=float, default=None,
                   help="restringe o closeness_scan ate este PERCENTIL do "
                        "bitscore (0-100). Evita amostrar cortes rigidos "
                        "demais, que so isolam singletons. Ex: 95")
    p.add_argument("--keep-orphans", action="store_true",
                   help="adiciona as orfas ao grafo como comunidade 'orphan'")
    p.add_argument("--progress", action="store_true", default=True)
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=stage_ssn, auto_start=True)

    p = sub.add_parser("annotate", help="taxonomia e arquitetura")
    p.add_argument("--out-dir", default=DIRS["annot"])
    p.add_argument("--db", default=NR)
    p.add_argument("--taxonomy", action="store_true", default=True)
    p.add_argument("--no-taxonomy", dest="taxonomy", action="store_false")
    p.add_argument("--arch", action="store_true", default=True)
    p.add_argument("--no-arch", dest="arch", action="store_false")
    p.set_defaults(func=stage_annotate)

    p = sub.add_parser("context", help="vizinhanca genomica")
    p.add_argument("--out-dir", default=DIRS["context"])
    p.add_argument("--db", default=NR)
    p.add_argument("--scope", default="all",
                   help="all | anchor | community:X | top:N")
    p.add_argument("--window", type=int, default=10)
    p.add_argument("--threads", type=int, default=15)
    p.add_argument("--annotation", default="pfam")
    p.add_argument("--arch", action="store_true", default=True)
    p.add_argument("--no-arch", dest="arch", action="store_false")
    p.add_argument("--warn-above", type=int, default=3000)
    p.add_argument("--progress", action="store_true", default=True)
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=stage_context)

    p = sub.add_parser("profile", help="descreve cada comunidade [bom ponto de parada]")
    p.add_argument("--context-dir", default=DIRS["context"])
    p.add_argument("--out-dir", default=DIRS["profile"])
    p.add_argument("--markers", default=None,
                   help="JSON de marcadores. SEM ele, entra em modo "
                        "DESCOBERTA (termos enriquecidos contra o fundo)")
    p.add_argument("--term-cols", nargs="+", default=None,
                   help="colunas da vizinhanca usadas na DESCOBERTA de termos. "
                        "Default: so 'pfam' (ou so 'product' se nao houver). "
                        "Somar pfam+product conta o mesmo evento duas vezes")
    p.add_argument("--top-communities", type=int, default=25)
    p.add_argument("--top-terms", type=int, default=15)
    p.add_argument("--min-term-count", type=int, default=3)
    p.set_defaults(func=stage_profile)

    p = sub.add_parser("select", help="filtro opcional guiado por config")
    p.add_argument("--hits-dir", default=DIRS["clustered"])
    p.add_argument("--context-dir", default=DIRS["context"])
    p.add_argument("--out-dir", default=DIRS["curated"])
    p.add_argument("--markers", default=None)
    p.add_argument("--communities", nargs="+", default=None,
                   help="restringe a estas comunidades")
    p.add_argument("--anchor-community", action="store_true",
                   help="so a comunidade da ancora (comportamento query-centrico)")
    p.add_argument("--min-markers", type=int, default=1)
    p.add_argument("--no-context", choices=["accept", "reject"], default="accept",
                   help="'accept' evita penalizar genoma mal anotado junto "
                        "com o errado")
    p.add_argument("--min-seqs", type=int, default=3)
    p.set_defaults(func=stage_select)

    p = sub.add_parser("build", help="alinha + hmmbuild por grupo")
    p.add_argument("--curated-dir", default=DIRS["curated"])
    p.add_argument("--out-dir", default=DIRS["models"])
    p.add_argument("--aligner", default="famsa")
    p.add_argument("--cpu", type=int, default=12)
    p.add_argument("--min-seqs", type=int, default=3)
    p.add_argument("--concat", action="store_true", default=True)
    p.add_argument("--concat-name", default="familias")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=stage_build)

    p = sub.add_parser("scan", help="modelos vs alvo + QC de separabilidade")
    p.add_argument("--hmm-dir", default=DIRS["models"])
    p.add_argument("--hits-dir", default=DIRS["clustered"])
    p.add_argument("--target", default=None)
    p.add_argument("--self-scan", action="store_true", default=True,
                   help="escaneia o proprio conjunto e cruza com a comunidade")
    p.add_argument("--no-self-scan", dest="self_scan", action="store_false")
    p.add_argument("--out-dir", default=DIRS["scan"])
    p.add_argument("--cpu", type=int, default=16)
    p.add_argument("--evalue", default="1e-5")
    p.add_argument("--concat-name", default="familias")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=stage_scan)

    p = sub.add_parser("table", help="inspeciona/exporta a tabela central")
    p.add_argument("-o", "--out", default=None)
    p.set_defaults(func=stage_table)

    p = sub.add_parser("report", help="monta Metodos + Resultados em Markdown")
    p.add_argument("-o", "--out", default="relatorio.md")
    p.set_defaults(func=stage_report)

    p = sub.add_parser("dashboard", help="painel HTML interativo das redes")
    p.add_argument("--ssn-dir", default=DIRS["ssn"])
    p.add_argument("--matrix-dir", default=DIRS["matrix"])
    p.add_argument("--out", default="painel.html")
    p.add_argument("--max-nodes", type=int, default=None,
                   help="teto OPCIONAL de nos por rede. Por padrao desenha a rede "
                        "inteira; use so se a simulacao ficar lenta. O grupo da "
                        "referencia e sempre preservado inteiro.")
    p.add_argument("--layout-iter", type=int, default=50)
    p.set_defaults(func=stage_dashboard)

    p = sub.add_parser("status", help="progresso de todos os estagios")
    p.add_argument("--verbose", action="store_true")
    p.set_defaults(func=stage_status)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
