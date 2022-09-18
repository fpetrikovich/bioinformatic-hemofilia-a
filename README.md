# Hemofilia A

## Instalar

Crear el entorno virtual con `virtualenv`:
```
virtualenv .env
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
python main.py -e 1 -i 1A
```
Corre el ejercicio 1 con el archivo GenBank del mRNA, factor coagulador VIII isoform A prepotein.

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