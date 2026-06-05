import pandas as pd
import sys
import pytest
from pathlib import Path
ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from deseq_utils import classify_gene, find_extremes, read_tsv

# ── Uso de IA ───────────────────────────────────────────────
# Herramienta : ChatGPT
# Prompt usado: "Ayúdame a escribir pruebas con pytest para algunas funciones
#                de mi archivo deseq_utils.py."
# Qué generó  : Sugerencias de pruebas para validar la clasificación de genes
#               como upregulated, downregulated y no_change, además de una
#               prueba para verificar que find_extremes() regrese las claves
#               esperadas.
# Qué modifiqué: Adapté las pruebas a los nombres reales de mis funciones,
#                usé ejemplos pequeños creados por mí en DataFrames de pandas,
#                revisé los errores de importación y ajusté el archivo para que
#                pytest pudiera encontrar mi módulo deseq_utils.py.
# ────────────────────────────────────────────────────────────

def test_read_tsv_file_not_found():
    with pytest.raises(FileNotFoundError):
        read_tsv("archivo_que_no_existe.tsv")

def test_clasificar_upregulated():
    """Un gen con padj=0.01 y lfc=3.5 debe ser upregulated."""
    df = pd.DataFrame([
        {
            "gene_id": "IFIT1",
            "log2FoldChange": 3.5,
            "padj": 0.01
        }
    ])

    resultado = classify_gene(df, umbral_padj=0.05, umbral_lfc=1.0)

    assert resultado.loc[0, "cambio"] == "upregulated"


def test_clasificar_downregulated():
    """Un gen con padj=0.01 y lfc=-2.0 debe ser downregulated."""
    df = pd.DataFrame([
        {
            "gene_id": "MYC",
            "log2FoldChange": -2.0,
            "padj": 0.01
        }
    ])

    resultado = classify_gene(df, umbral_padj=0.05, umbral_lfc=1.0)

    assert resultado.loc[0, "cambio"] == "downregulated"

def test_clasificar_no_change_by_lfc():
    """Un gen significativo pero con lfc bajo debe ser no_change."""
    df = pd.DataFrame([
        {
            "gene_id": "GAPDH",
            "log2FoldChange": 0.3,
            "padj": 0.01
        }
    ])

    resultado = classify_gene(df, umbral_padj=0.05, umbral_lfc=1.0)

    assert resultado.loc[0, "cambio"] == "no_change"

def test_clasificar_no_change_by_padj():
    """Un gen con padj alto debe ser no_change aunque tenga lfc alto."""
    df = pd.DataFrame([
        {
            "gene_id": "GAPDH",
            "log2FoldChange": 4.0,
            "padj": 0.8
        }
    ])

    resultado = classify_gene(df, umbral_padj=0.05, umbral_lfc=1.0)

    assert resultado.loc[0, "cambio"] == "no_change"


def test_find_extremes_returns_correct_keys():
    """find_extremes debe retornar las claves esperadas."""
    df = pd.DataFrame([
        {
            "gene_id": "IFIT1",
            "cambio": "upregulated",
            "description": "interferon induced protein",
            "log2FoldChange": 5.0,
            "padj": 0.01
        },
        {
            "gene_id": "MYC",
            "cambio": "downregulated",
            "description": "MYC proto-oncogene",
            "log2FoldChange": -3.0,
            "padj": 0.02
        },
        {
            "gene_id": "RPS6",
            "cambio": "downregulated",
            "description": "ribosomal protein S6",
            "log2FoldChange": -1.5,
            "padj": 0.001
        }
    ])

    extremos = find_extremes(df)

    assert "mas_inducido" in extremos
    assert "mas_reprimido" in extremos
    assert "mas_significativo" in extremos