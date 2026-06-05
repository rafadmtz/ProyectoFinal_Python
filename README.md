<!--
── Uso de IA ───────────────────────────────────────────────
Herramienta : ChatGPT
Prompt usado: "Ayúdame a hacer un README.md en base a la información y flujo
de mi proyecto."
Qué generó  : Una propuesta de estructura para documentar el objetivo del proyecto,
los archivos de entrada, la instalación con uv, la ejecución del programa, las
pruebas con pytest y los archivos de salida.
Qué modifiqué: Revisé que los comandos, nombres de archivos, funciones y rutas
coincidieran con mi proyecto final. Realice los cambios que el README.md necesitaba para que fuera representativo del programa.
────────────────────────────────────────────────────────────
-->

# Proyecto Final — Análisis de Genes Diferencialmente Expresados

Este proyecto analiza resultados de expresión diferencial obtenidos con DESeq2 para comparar células A549 infectadas con influenza A contra células control. El programa clasifica genes como inducidos, reprimidos o sin cambio, agrega anotaciones funcionales desde un archivo GFF3 y genera archivos de salida interpretables.

## Objetivo

Procesar un archivo TSV de resultados de DESeq2 y un archivo GFF3 de anotaciones para:

* Clasificar genes según `padj` y `log2FoldChange`.
* Identificar genes `upregulated`, `downregulated` y `no_change`.
* Agregar descripciones funcionales a genes significativos.
* Reportar genes extremos.
* Guardar archivos de salida en formato TSV y TXT.
* Generar opcionalmente un volcano plot.

## Estructura del proyecto

```text
ProyectoFinal_python/
├── analyze_degs.py
├── deseq_utils.py
├── pyproject.toml
├── uv.lock
├── .gitignore
├── README.md
├── diagrama_flujo.md
├── pytest.ini
├── datos/
│   ├── iav_deseq2_results.tsv
│   └── human_genes.gff
├── tests/
│   └── test_deseq_utils.py
└── results/
    ├── upregulated_genes.tsv
    ├── downregulated_genes.tsv
    ├── summary_report.txt
    └── volcano_plot.png
```

## Archivos de entrada

Los archivos de entrada se encuentran en la carpeta `datos/`.

### `datos/iav_deseq2_results.tsv`

Archivo TSV con los resultados de expresión diferencial. Contiene columnas como:

* `gene_id`
* `baseMean`
* `log2FoldChange`
* `lfcSE`
* `stat`
* `pvalue`
* `padj`

### `datos/human_genes.gff`

Archivo GFF3 con anotaciones funcionales. El programa ignora las líneas de comentario que empiezan con `#` y utiliza la columna `attributes` para extraer:

* `Name`
* `description`

Con esta información construye un diccionario en formato:

```python
{gene_id: description}
```

## Instalación del entorno

Este proyecto usa `uv` para manejar el entorno y las dependencias.

Para instalar las dependencias desde el proyecto:

```bash
uv sync
```

Si es necesario agregar dependencias manualmente:

```bash
uv add pandas matplotlib numpy
uv add --dev pytest
```

## Ejecución del programa

Para correr el análisis con los valores por defecto:

```bash
uv run python analyze_degs.py --input datos/iav_deseq2_results.tsv --gff datos/human_genes.gff
```

Por defecto se usan:

* `padj-threshold = 0.05`
* `lfc-threshold = 1.0`
* `output-dir = results/`

También se puede ejecutar especificando todos los argumentos:

```bash
uv run python analyze_degs.py \
    --input datos/iav_deseq2_results.tsv \
    --gff datos/human_genes.gff \
    --padj-threshold 0.05 \
    --lfc-threshold 1.0 \
    --output-dir results/
```

En PowerShell, el comando puede escribirse en una sola línea:

```powershell
uv run python analyze_degs.py --input datos/iav_deseq2_results.tsv --gff datos/human_genes.gff --padj-threshold 0.05 --lfc-threshold 1.0 --output-dir results/
```

## Argumentos del programa

| Argumento          | Tipo  | Requerido | Default    | Descripción                                        |
| ------------------ | ----- | --------- | ---------- | -------------------------------------------------- |
| `--input`          | str   | Sí        | —          | Ruta al archivo TSV de DESeq2                      |
| `--gff`            | str   | Sí        | —          | Ruta al archivo GFF3 de anotaciones                |
| `--padj-threshold` | float | No        | `0.05`     | Umbral de significancia estadística                |
| `--lfc-threshold`  | float | No        | `1.0`      | Umbral mínimo de cambio en log2FoldChange          |
| `--output-dir`     | str   | No        | `results/` | Directorio donde se guardan los archivos de salida |

## Flujo general del programa

El programa inicia leyendo los argumentos de línea de comandos, incluyendo las rutas del archivo TSV de DESeq2, el archivo GFF3, los umbrales de significancia y el directorio de salida. Después, carga el TSV con los resultados de expresión diferencial y el GFF3 con las anotaciones funcionales.

A partir del GFF3 se extraen los campos `Name` y `description` de la columna `attributes`, construyendo un diccionario que relaciona cada gen con su descripción. Luego, cada gen del archivo DESeq2 se clasifica como `upregulated`, `downregulated` o `no_change` según los umbrales de `padj` y `log2FoldChange`.

Finalmente, el programa agrega anotaciones funcionales a los genes significativos, identifica los genes extremos, guarda los archivos `upregulated_genes.tsv`, `downregulated_genes.tsv` y `summary_report.txt`, y opcionalmente genera un volcano plot.

## Criterios de clasificación

| Categoría       | Condición                                              |
| --------------- | ------------------------------------------------------ |
| `upregulated`   | `padj < umbral_padj` y `log2FoldChange >= umbral_lfc`  |
| `downregulated` | `padj < umbral_padj` y `log2FoldChange <= -umbral_lfc` |
| `no_change`     | Todo lo demás                                          |

## Archivos de salida

El programa genera los siguientes archivos dentro del directorio indicado por `--output-dir`.

### `upregulated_genes.tsv`

Contiene los genes inducidos con las columnas:

```text
gene_id    log2FoldChange    padj    description
```

### `downregulated_genes.tsv`

Contiene los genes reprimidos con las columnas:

```text
gene_id    log2FoldChange    padj    description
```

### `summary_report.txt`

Contiene un resumen completo del análisis:

* Archivo analizado.
* Umbrales utilizados.
* Total de genes analizados.
* Conteo y porcentaje por categoría.
* Genes extremos.
* Lista de genes diferencialmente expresados con descripción.

### `volcano_plot.png`

* Eje X: `log2FoldChange`
* Eje Y: `-log10(padj)`
* Rojo: genes `upregulated`
* Azul: genes `downregulated`
* Gris: genes `no_change`

## Genes extremos

El programa identifica entre los genes significativos:

* Gen con mayor `log2FoldChange`.
* Gen con menor `log2FoldChange`.
* Gen con menor `padj`.

Estos genes se imprimen en pantalla y también se reportan en `summary_report.txt`.

## Pruebas con pytest

Las pruebas se encuentran en:

```text
tests/test_deseq_utils.py
```

Para ejecutarlas:

```bash
uv run pytest tests/ -v
```

O también:

```bash
uv run python -m pytest tests/ -v
```

Las pruebas verifican funciones de clasificación de genes y la identificación de genes extremos.

## Archivos principales

### `analyze_degs.py`

Script principal del proyecto. Se encarga de:

* Leer argumentos con `argparse`.
* Llamar las funciones de `deseq_utils.py`.
* Ejecutar el flujo completo del análisis.
* Guardar los archivos de salida.

### `deseq_utils.py`

Módulo auxiliar con funciones reutilizables para:

* Leer archivos TSV.
* Leer y procesar archivos GFF3.
* Parsear atributos.
* Clasificar genes.
* Agregar anotaciones funcionales.
* Encontrar genes extremos.
* Guardar archivos de salida.
* Generar volcano plot.
