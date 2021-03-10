.PHONY: all

all: data/covid_cantab.csv

data/cantab_soc.csv data/cantab_mot.csv: 00_pull_cantab.R
	./00.0_pull_cantab.R

data/cantab_pet7T_motsoc.csv: data/cantab_soc.csv data/cantab_mot.csv
	./00.1_merge_cantab.R

data/qualtrics_adult-kid.csv: ./qualtrics_covid.R
	./qualtrics_covid.R

data/covid_cantab.csv: data/qualtrics_adult-kid.csv data/cantab_pet7T_motsoc.csv
	./00.2_cantab_covid.R
