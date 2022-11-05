# Hemofilia A

## Instalar

Crear el entorno virtual con `virtualenv`:
```
virtualenv .env
o
python3 -m venv "virtualenv"
```

Activar el entorno:
```
source .env/bin/activate
```

Instalar dependencias:
```
pip install -r requirements.txt
```

### Ejemplos de Llamados

Ejemplos de los llamados:
```
python main.py -e 1 -gb 1A
```
Corre el ejercicio 1 con el archivo GenBank del mRNA, factor coagulador VIII isoform A prepotein.
```
python3 main.py -e 2 -q 2B -r blast
```
Corre el ejercicio 2 con el FASTA file de la proteina B generada al correr el ejercicio 1, creando reportes para cada ORF al realizar consultas online y locales a BLAST. 
```
python main.py -e 3 -ss species/species_file.fasta -out ex3.out
```
Corre el ejercicio 3 con el FASTA file generado en el ejercicio 1, y los compara con las especies resultantes del BLAST. 
```
python main.py -e 4 -b sp_ORF1.report -p "mus musculus"
```
Correra el ejercicio 4 con el reporte sp_ORF1.report viendo si encuentra el patron 
musculus en las descripciones de los alineamientos.
```
python main.py -e 5 -seq genbank/NM_000132.4_sequence.gb
```
Corre el ejercicio 5. Utiliza la(s) secuencia(s) en el archivo genbank/NM_000132.4_sequence.gb para hallar los posibles orfs y realizar un análisis de los dominios de los orfs. Se generaran dos archivos de salida con el nombre por default: protein_orf.orf y protein_orf.patmatmotifs.


### Setup de Blast Local - Linux
Descargarse el tar adecuado desde:
```
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
```
Descomprimirlo en la carpeta actual con el comando:
```
tar -xvf ncbi-blast-2.13.0+-x64-linux.tar.gz
```
Descargarse la base de datos Swissport desde ```https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz``` y descomprimirla en ```ncbi-blast-2.13.0+/data```:
```
tar -xvf swissprot.tar.gz -C ncbi-blast-2.13.0+/data
```
Moverse al root del repositorio y correr una consulta en Blast desde la linea de comandos para verificar su funcionamiento correcto:
```
./ncbi-blast-2.13.0+/bin/blastp -db ncbi-blast-2.13.0+/data/swissprot -query <archivo_fasta>  -out reports/<archivo_salida>
```