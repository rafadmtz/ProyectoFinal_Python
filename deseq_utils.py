# ── Uso de IA ───────────────────────────────────────────────

# Herramienta : ChatGPT

# Prompt usado: "Ayúdame a revisar, corregir y completar mis bosquejos de código

# para leer un GFF3, parsear la columna attributes, clasificar genes

# de DESeq2, agregar anotaciones funcionales y guardar archivos

# de salida."

# Qué generó  : Explicaciones paso a paso, correcciones de sintaxis y lógica,

# y sugerencias para estructurar funciones como leer_gff(),

# parse_attributes(), clasificar_genes(),

# agregar_anotacion_funcional(), find_extremes()

# y guardar_archivos_salida().

# Qué modifiqué: Adapté los nombres de variables y funciones a mi proyecto,

# integré las partes que entendí en mi propio flujo de trabajo,

# probé el código con mis archivos, ajusté rutas, umbrales,

# columnas de salida y mensajes para cumplir los requisitos

# del examen.

# ────────────────────────────────────────────────────────────

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def read_tsv(filepath: str) -> pd.DataFrame:
    """
    Lee un archivo TSV y lo convierte en un DataFrame de pandas.

    Parameters
    ----------
    filepath : str
        Ruta del archivo TSV que se desea leer.

    Returns
    -------
    pd.DataFrame
        DataFrame con el contenido del archivo TSV.
        Regresa None si el archivo no se encuentra o si hay valores inválidos.
    """
    try:
        df = pd.read_csv(filepath, sep="\t")

        columnas_numericas = [
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj"
        ]

        for columna in columnas_numericas:
            df[columna] = pd.to_numeric(df[columna], errors="raise")

        print(f"Genes cargados del TSV: {len(df)}")

        return df

    except FileNotFoundError as error:
        raise FileNotFoundError(
        f"No se encontró el archivo: {filepath}"
        ) from error

    except ValueError as error:
        raise ValueError(
            f"Error al leer valores numericos del TSV: {error}"
        ) from error

def read_gff(filepath: str) -> pd.DataFrame:
    """
    Lee un archivo GFF3 ignorando las líneas de comentario.

    Parameters
    ----------
    filepath : str
        Ruta del archivo GFF3 que se desea leer.

    Returns
    -------
    pd.DataFrame
        DataFrame con las 9 columnas del archivo GFF3.
        Regresa None si el archivo no se encuentra.
    """
    columnas = [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes"
    ]

    try:
        gff = pd.read_csv(
            filepath,
            sep="\t",
            comment="#",
            header=None,
            names=columnas
        )

    except FileNotFoundError as error:
        raise FileNotFoundError(
            f"No se encontró el archivo: {filepath}"
        ) from error

    except ValueError as error:
        raise ValueError(
            f"Error al leer el archivo GFF: {error}"
        ) from error

    return gff

def parse_attributes(attributes: str) -> dict:
    """
    Convierte la columna attributes de un archivo GFF3 en un diccionario.

    Parameters
    ----------
    attributes : str
        Texto de la columna attributes del GFF3.
        Ejemplo:
        ID=...;Name=IFIT1;description=...;gene_type=...

    Returns
    -------
    dict
        Diccionario con los campos separados.
        Ejemplo:
        {
            "ID": "...",
            "Name": "IFIT1",
            "description": "...",
            "gene_type": "protein_coding"
        }
    """
    info = {}

    fields = attributes.split(";")

    for field in fields:
        if "=" in field:
            key, value = field.split("=", 1)
            info[key] = value

    return info

def classify_gene(
    deseq_df: pd.DataFrame,
    umbral_padj: float,
    umbral_lfc: float
) -> pd.DataFrame:
    """
    Clasifica genes según padj y log2FoldChange.

    Parameters
    ----------
    deseq_df : pd.DataFrame
        DataFrame con los resultados de DESeq2.
    umbral_padj : float
        Umbral máximo de padj para considerar significancia.
    umbral_lfc : float
        Umbral mínimo de log2FoldChange.
    
    Returns
    -------
    pd.DataFrame
        DataFrame con gene_id y cambio.
    """
    genes_clasificados = []

    for index, gene in deseq_df.iterrows():
        if gene["padj"] < umbral_padj and gene["log2FoldChange"] >= umbral_lfc:
            cambio = "upregulated"

        elif gene["padj"] < umbral_padj and gene["log2FoldChange"] <= -umbral_lfc:
            cambio = "downregulated"

        else:
            cambio = "no_change"

        genes_clasificados.append({
            "gene_id": gene["gene_id"],
            "cambio": cambio,
            "log2FoldChange": gene["log2FoldChange"],
            "padj": gene["padj"]
        })

    return pd.DataFrame(genes_clasificados)

def imprimir_resumen(
    clasificacion_df: pd.DataFrame,
    archivo_analizado: str = None,
    umbral_padj: float = None,
    umbral_lfc: float = None
) -> None:
    """
    Imprime el resumen de genes clasificados por categoría.

    Parameters
    ----------
    clasificacion_df : pd.DataFrame
        DataFrame con una columna llamada 'cambio'.
    archivo_analizado : str
        Ruta del archivo analizado.
    umbral_padj : float
        Umbral usado para padj.
    umbral_lfc : float
        Umbral usado para log2FoldChange.

    Returns
    -------
    None
        La función solo imprime información en pantalla.
    """
    total_genes = len(clasificacion_df)

    conteos = clasificacion_df["cambio"].value_counts()

    print("Resumen del analisis")
    print("--------------------")

    if archivo_analizado is not None:
        print(f"Archivo analizado: {archivo_analizado}")

    if umbral_padj is not None:
        print(f"Umbral padj: {umbral_padj}")

    if umbral_lfc is not None:
        print(f"Umbral log2FoldChange: {umbral_lfc}")

    print(f"Total de genes analizados: {total_genes}")
    print()

    categorias = ["upregulated", "downregulated", "no_change"]

    for categoria in categorias:
        cantidad = conteos.get(categoria, 0)

        if total_genes > 0:
            porcentaje = (cantidad / total_genes) * 100
        else:
            porcentaje = 0

        print(f"{categoria}: {cantidad} genes ({porcentaje:.2f}%)")

def agregar_anotacion_funcional(
    clasificacion_df: pd.DataFrame,
    gene_dict: dict
) -> pd.DataFrame:
    """
    Agrega la descripción funcional a los genes upregulated y downregulated.

    Parameters
    ----------
    clasificacion_df : pd.DataFrame
        DataFrame con los genes clasificados. Debe contener las columnas
        gene_id, log2FoldChange, padj y cambio.
    gene_dict : dict
        Diccionario con anotaciones del GFF en formato {gene_id: description}.

    Returns
    -------
    pd.DataFrame
        DataFrame con genes significativos y su descripción funcional.
    """
    anotacion_funcional = []

    for i, gene in clasificacion_df.iterrows():
        if gene["cambio"] == "upregulated" or gene["cambio"] == "downregulated":
            description = gene_dict.get(gene["gene_id"], "sin anotación")

            anotacion_funcional.append({
                "gene_id": gene["gene_id"],
                "cambio": gene["cambio"],
                "description": description,
                "log2FoldChange": gene["log2FoldChange"],
                "padj": gene["padj"]
            })

    return pd.DataFrame(anotacion_funcional)



def find_extremes(anotacion_df):
    """
    Identifica genes extremos entre los genes significativos.

    Parameters
    ----------
    anotacion_df : pd.DataFrame
        DataFrame con genes upregulated y downregulated.

    Returns
    -------
    dict
        Diccionario con el gen más inducido, más reprimido y más significativo.
    """
    if anotacion_df.empty:
        return {
            "mas_inducido": None,
            "mas_reprimido": None,
            "mas_significativo": None
        }

    mas_inducido = anotacion_df.loc[anotacion_df["log2FoldChange"].idxmax()]
    mas_reprimido = anotacion_df.loc[anotacion_df["log2FoldChange"].idxmin()]
    mas_significativo = anotacion_df.loc[anotacion_df["padj"].idxmin()]

    return {
        "mas_inducido": mas_inducido,
        "mas_reprimido": mas_reprimido,
        "mas_significativo": mas_significativo
    }

from pathlib import Path
import pandas as pd


def guardar_archivos_salida(
    anotacion_df: pd.DataFrame,
    clasificacion_df: pd.DataFrame,
    extremos: dict,
    output_dir: str,
    archivo_analizado: str,
    umbral_padj: float,
    umbral_lfc: float
) -> None:
    """
    Guarda los archivos de salida del análisis.

    Parameters
    ----------
    anotacion_df : pd.DataFrame
        DataFrame con genes upregulated y downregulated, incluyendo description.
    clasificacion_df : pd.DataFrame
        DataFrame con todos los genes clasificados.
    extremos : dict
        Diccionario con los genes extremos.
    output_dir : str
        Directorio donde se guardarán los archivos.
    archivo_analizado : str
        Ruta del archivo DESeq2 analizado.
    umbral_padj : float
        Umbral de padj usado.
    umbral_lfc : float
        Umbral de log2FoldChange usado.

    Returns
    -------
    None
        Guarda archivos en el directorio indicado.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    columnas_tsv = ["gene_id", "log2FoldChange", "padj", "description"]

    upregulated_df = anotacion_df[anotacion_df["cambio"] == "upregulated"]
    downregulated_df = anotacion_df[anotacion_df["cambio"] == "downregulated"]

    upregulated_df[columnas_tsv].to_csv(
        output_path / "upregulated_genes.tsv",
        sep="\t",
        index=False
    )

    downregulated_df[columnas_tsv].to_csv(
        output_path / "downregulated_genes.tsv",
        sep="\t",
        index=False
    )

    total_genes = len(clasificacion_df)
    conteos = clasificacion_df["cambio"].value_counts()

    with open(output_path / "summary_report.txt", "w", encoding="utf-8") as file:
        file.write("Resumen del analisis de genes diferencialmente expresados\n")
        file.write("========================================================\n\n")

        file.write(f"Archivo analizado: {archivo_analizado}\n")
        file.write(f"Umbral padj: {umbral_padj}\n")
        file.write(f"Umbral log2FoldChange: {umbral_lfc}\n")
        file.write(f"Total de genes analizados: {total_genes}\n\n")

        file.write("Conteo por categoria\n")
        file.write("--------------------\n")

        for categoria in ["upregulated", "downregulated", "no_change"]:
            cantidad = conteos.get(categoria, 0)
            porcentaje = (cantidad / total_genes) * 100 if total_genes > 0 else 0

            file.write(f"{categoria}: {cantidad} genes ({porcentaje:.2f}%)\n")

        file.write("\nGenes extremos\n")
        file.write("--------------\n")

        for nombre, gen in extremos.items():
            if gen is not None:
                file.write(
                    f"{nombre}: {gen['gene_id']} "
                    f"log2FoldChange={gen['log2FoldChange']} "
                    f"padj={gen['padj']}\n"
                )
            else:
                file.write(f"{nombre}: no disponible\n")
        file.write("\nGenes diferencialmente expresados con descripcion\n")
        file.write("-------------------------------------------------\n")

        anotacion_df.to_csv(
            file,
            sep="\t",
            index=False,
            columns=["gene_id", "cambio", "log2FoldChange", "padj", "description"]
        )

def guardar_volcano_plot(
    clasificacion_df: pd.DataFrame,
    output_dir: str,
    umbral_padj: float,
    umbral_lfc: float
) -> None:
    """
    Genera y guarda un volcano plot del análisis diferencial.

    Parameters
    ----------
    clasificacion_df : pd.DataFrame
        DataFrame con todos los genes clasificados. Debe contener las columnas
        gene_id, log2FoldChange, padj y cambio.
    output_dir : str
        Directorio donde se guardará la imagen.
    umbral_padj : float
        Umbral de padj usado para significancia.
    umbral_lfc : float
        Umbral de log2FoldChange usado para clasificar genes.

    Returns
    -------
    None
        Guarda el gráfico como volcano_plot.png.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    df = clasificacion_df.copy()

    df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
    df["log2FoldChange"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")

    df = df.dropna(subset=["padj", "log2FoldChange"])

    df["padj_plot"] = df["padj"].replace(0, np.nextafter(0, 1))
    df["minus_log10_padj"] = -np.log10(df["padj_plot"])

    colores = {
        "upregulated": "red",
        "downregulated": "blue",
        "no_change": "gray"
    }

    plt.figure(figsize=(8, 6))

    for categoria, color in colores.items():
        subset = df[df["cambio"] == categoria]

        plt.scatter(
            subset["log2FoldChange"],
            subset["minus_log10_padj"],
            c=color,
            label=categoria,
            alpha=0.7,
            s=25
        )

    plt.axvline(
        x=umbral_lfc,
        linestyle="--",
        color="black",
        linewidth=1
    )

    plt.axvline(
        x=-umbral_lfc,
        linestyle="--",
        color="black",
        linewidth=1
    )

    plt.axhline(
        y=-np.log10(umbral_padj),
        linestyle="--",
        color="black",
        linewidth=1
    )

    top_genes = df.sort_values("padj").head(5)

    for _, gene in top_genes.iterrows():
        plt.text(
            gene["log2FoldChange"],
            gene["minus_log10_padj"],
            gene["gene_id"],
            fontsize=8
        )

    plt.title("Volcano plot: IAV vs Mock")
    plt.xlabel("log2FoldChange")
    plt.ylabel("-log10(padj)")
    plt.legend()
    plt.tight_layout()

    plt.savefig(output_path / "volcano_plot.png", dpi=300)
    plt.close()