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

### Setup de la base de datos PROSITE
Instalarse el paquete Emboss:
```
sudo apt-get install emboss
```
Dentro del repositorio, moverse hacia el directorio emboss:
```
cd emboss
```
Bajarse la base de datos PROSITE mediante wget:
```
wget https://ftp.expasy.org/databases/prosite/prosite.dat
wget https://ftp.expasy.org/databases/prosite/prosite.doc
```
Procesar la base de datos PROSITE para que pueda ser usada por patmatmotifs: 
```
sudo prosextract -prositedir . 
```

### Setup de clustalW
Instalación del módulo:
```
sudo apt install clustalw
```

## Ejecucion de cada ejercicio

Para ejecutar los ejercicios se puede usar un archivo de configuracion o poner los parametros de manera explicita en la linea de comandos.

### Usando linea de comandos de manera explicita

El programa se corre como sigue:
```
python main.py --exercise <num> <PARAMS>
```
El ```num``` es el ejercicio que se desea ejecutar, del 1 al 5. Los PARAMS que se deben usar en cada caso se pueden ver corriendo la ayuda:
```
python main.py --help
```
Algunos parametros tienen valores por default si no se agregan. Estos valores default se pueden ver en el help y en el informe. Ver el informe tambien para una explicacion mas detallada de que hace cada parametro y que valores puede tener.

#### Ejemplos de Llamados

Ejemplos de los llamados:
```
python main.py -e 1 -gb 1B
```
Corre el ejercicio 1 con el archivo GenBank del mRNA, factor coagulador VIII isoform A prepotein.
```
python3 main.py -e 2 -q 2A -r blast
```
Corre el ejercicio 2 con el FASTA file de la proteina isoform A generada al correr el ejercicio 1, creando reportes para cada ORF al realizar consultas online a BLAST usando la base de datos swissprot. 
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
Corre el ejercicio 5. Utiliza la(s) secuencia(s) en el archivo genbank/NM_000132.4_sequence.gb para hallar los posibles orfs (con valor minimo de 2700 por default) y realizar un análisis de los dominios de los orfs. Se generaran dos archivos de salida con el nombre por default: protein_orf.orf y protein_orf.patmatmotifs.

### Usando archivos de configuracion

El archivo de configuracion tiene que estar en la carpeta configuration_files/ y debe ser de extension .ini.

Ahi, se pueden especificar valores default y valores por ejercicios. Los valores de cada ejercicio deben ir debajo de su header adecuado ([ex1], [ex2], etc.).

No olvidarse de agregar todos los parametros, ya sea en el default o debajo del header del ejercicio. 

Para mas informacion de que parametros requiere cada ejercicio, puede ver el informe o correr el siguiente comando de ayuda:
```
python main.py --help
```
Los keys de los parametros en el archivo de configuracion deben tener el mismo nombre que los parametros al correrlo en linea de comandos (ej: --PARAM entonces el key es PARAM).

El comando a correr una vez que este todo seteado:
```
python main.py --exercise <num> --use_config --config_file <filename>
``` 
Donde ``` num``` es un numero de 1 a 5, ```use_config``` indica que se quiere usar un archivo de configuracion y ```config_file``` es que archivo se quiere usar, el cual es ```configuration.ini``` por default.
