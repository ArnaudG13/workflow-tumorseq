# Titre du projet

Dépot de code pour les projets NGS du laboratoire d'oncologie moléculaire (CRCM)

### Prérequis

* conda >= 4.3.22 https://conda.io/docs/index.html
* snakemake >= 3.13.3 https://snakemake.readthedocs.io/en/stable/

### Installation

Installation de conda :

```bash Miniconda3-latest-Linux-x86_64.sh```


Installation de snakemake :

```conda install -c bioconda snakemake```

###Fichier de log

Par défaut les fichiers de log du cluster de calcul sont générés dans

 ```
/home/guille/logs/cluster/snakemake/
```

Vous pouvez modifier le chemin de destination de ces logs dans le fichier cluster.yml et créer le dossier s'il n'existe pas.


### Utilisation

Ajouter ces 3 lignes à votre .bashrc afin d'effacer vos variables d'environnement Python et Java

```
export PYTHONNOUSERSITE=True
unset PYTHONPATH
unset JAVA_HOME
```

Attention !!! Il est préférable de lancer cette commande en étant à la racine de son homedir, créant ainsi un dossier .snakemake pour toutes les analyses. 

Exemple de ligne de commande pour le pipeline somatique avec utilisation du cluster de calcul (oar) :

```
snakemake -j 9999 --snakefile /home/data/ngs_om/travail/scripts/snakefiles/Snakefile_tumor_normal_variant_analysis \
--configfile /home/data/ngs_om/travail/projects/test_Quentin/config_somatic.yml \
--use-conda \
--cluster-config /home/dacosta/cluster.yml \
--cluster "oarsub -l host={cluster.host}/core={cluster.core},walltime={cluster.time} -O {cluster.stdout} -E {cluster.stderr}"
```

Exemple de ligne de commande pour le pipeline Tumour-only avec utilisation du cluster de calcul (oar)

```
snakemake -j 9999 --snakefile /home/data/ngs_om/travail/scripts/snakefiles/Snakefile_tumor_only_variant_analysis \
--configfile /home/data/ngs_om/travail/scripts/snakefiles/config/config_tumor_only.yml \
--use-conda \
--cluster-config /home/data/ngs_om/travail/scripts/snakefiles/config/cluster.yml \
--cluster "oarsub -l host={cluster.host}/core={cluster.core},walltime={cluster.time} -O {cluster.stdout} -E {cluster.stderr}"
```

## Fichiers de configurations

* --configfile : fichier de configuration pour l'analyse (Panel, chemins, etc...)
* --cluster-config : fichier de configuration du cluster (chemins des logs, etc...)

Attention !!!  Il faut obligatoirement modifier le fichier de configuration du cluster pour les fichiers de log. 

## Développé Avec

* [Python](https://www.python.org/) - Language de programmation Python
* [bash](http://git.savannah.gnu.org/cgit/bash.git) - Language de programmation Bash
* [snakemake](https://bitbucket.org/johanneskoester/snakemake/wiki/Home) - Gestionnaire de workflow snakemake

## Auteurs

* **Arnaud Guille**
* **Quentin Da Costa**

