# ── Uso de IA ───────────────────────────────────────────────

# Herramienta : ChatGPT

# Prompt usado: "Ayúdame a organizar el script principal para ejecutar el análisis

# usando funciones importadas desde deseq_utils.py."

# Qué generó  : Orientación para importar el módulo deseq_utils, llamar funciones

# de lectura, clasificación, anotación, resumen y escritura de

# archivos de salida.

# Qué modifiqué: Organicé el flujo principal del programa, conecté las funciones

# con mis variables, ajusté los argumentos, rutas y parámetros

# para que el script funcionara en mi proyecto.

# ────────────────────────────────────────────────────────────

import argparse
import deseq_utils as du

def parse_args() -> argparse.Namespace:
    """
    Lee los argumentos de línea de comandos para el análisis DESeq2.

    Returns
    -------
    argparse.Namespace
        Objeto con las rutas de entrada, umbrales y directorio de salida.
    """
    parser = argparse.ArgumentParser(
        description="Análisis de genes diferencialmente expresados en respuesta a IAV."
    )

    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Ruta al archivo TSV de resultados de DESeq2."
    )

    parser.add_argument(
        "--gff",
        type=str,
        required=True,
        help="Ruta al archivo GFF3 de anotaciones."
    )

    parser.add_argument(
        "--padj-threshold",
        type=float,
        default=0.05,
        help="Umbral de significancia estadística para padj. Default: 0.05."
    )

    parser.add_argument(
        "--lfc-threshold",
        type=float,
        default=1.0,
        help="Umbral mínimo de cambio en log2FoldChange. Default: 1.0."
    )

    parser.add_argument(
        "--output-dir",
        type=str,
        default="results/",
        help="Directorio donde se guardarán los archivos de salida. Default: results/."
    )

    return parser.parse_args()

def main():
    #Recibe los argumentos de linea de comandos
    args = parse_args()

    # Lee el archivo TSV de DESeq2
    deseq_df = du.read_tsv(args.input)

    # Lee el archivo GFF3 para obtener anotaciones de genes
    gff = du.read_gff(args.gff)
    
    #Diccionario para almacenar las anotaciones de los genes, con gene_id como clave y description como valor
    gene_dict = {}

    for attributes in gff["attributes"]:
        info = du.parse_attributes(attributes)

        gene_id = info.get("Name")
        description = info.get("description", "sin anotacion")

        if gene_id is not None:
            gene_dict[gene_id] = description

    # Define los umbrales para clasificar los genes
    umbral_padj = args.padj_threshold
    umbral_lfc = args.lfc_threshold

    # Clasifica los genes segun los umbrales de padj y log2FoldChange
    clasificacion_df = du.classify_gene(deseq_df, umbral_padj, umbral_lfc)

    du.imprimir_resumen(
        clasificacion_df,
        args.input,
        umbral_padj,
        umbral_lfc
    )

    # Agrega anotacion funcional a los genes clasificados usando el diccionario de anotaciones
    anotacio_df=du.agregar_anotacion_funcional(clasificacion_df, gene_dict)

    # Identifica los genes con cambios extremos (mayor log2FoldChange positivo y negativo entre los significativos)
    extremos = du.find_extremes(anotacio_df)
    print("Genes extremos:")
    
    for categoria, gene in extremos.items():
        if gene is not None:
            print(f"  {categoria}: {gene['gene_id']} ({gene['description']})")
        else:
            print(f"  {categoria}: No se encontro un gen que cumpla los criterios.")
            
    du.guardar_archivos_salida(
        anotacio_df,
        clasificacion_df,
        extremos,
        args.output_dir,
        args.input,
        umbral_padj,
        umbral_lfc
    )

    du.guardar_volcano_plot(
        clasificacion_df,
        args.output_dir,
        umbral_padj,
        umbral_lfc
    )

if __name__ == "__main__":
    main()