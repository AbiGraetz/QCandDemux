#!/usr/bin/env bash
set -euo pipefail

# Create conda environment with required packages
if ! conda info --envs | grep -q "QCandDemux"; then
	conda create -n QCandDemux -y
	conda install -c bioconda -y py nanopack seqkit vsearch edlib
fi

# Create output directories
mkdir -p ./DoradoDemux/16S ./DoradoDemux/ITS
mkdir -p ./MinibarDemux/16S ./MinibarDemux/ITS

echo "Please ensure that path variables are set. If you have not edited these, this script may not run correctly."

# Set paths for required reference files 
VSEARCH_DB_16S="./VSEARCH_DB_16S.fasta" # VSEARCH database for 16S sequences
VSEARCH_DB_ITS="./VSEARCH_DB_ITS.fasta"
DEMUX_16S="./MinibarDemux/16S_demux_file.txt" # Minibar demultiplexing file for 16S data
DEMUX_ITS="./MinibarDemux/ITS_demux_file.txt" # Minibar demultiplexing file for ITS data
BARCODE_FILE="./Dorado_BC_file.txt" # Barcode file for renaming Dorado demuxed files
MINIBAR_PATH="./minibar.py" # Path to minibar script
export PATH=/path/to/dorado/install:$PATH # Path to Dorado installation

# Check input files
for f in calls.fastq "$VSEARCH_DB_16S" "$VSEARCH_DB_ITS" "$DEMUX_16S" "$DEMUX_ITS" "$DORADO_BC"; do
    if [[ ! -f "$f" ]]; then
        echo "Error: Required file $f not found!"
        exit 1
    fi
done

# Initial quality visualisation and filtering
nanoplot -t 8 --fastq calls.fastq --plots kde --maxlength 7500
chopper --minlength 800 --maxlength 7500 -q 15 -t 8 calls.fastq | gzip > calls.Qmin15.fq.gz

# Dorado demux by ONT native barcode
dorado demux --kit-name SQK-NBD114-96 --output-dir ./DoradoDemux --emit-fastq calls.Qmin15.fq.gz

# Check sequences per barcode
for file in ./DoradoDemux/*.fastq; do
	if [ -f "$file" ]; then
	lines=$(wc -l < "$file")
	sequences=$((lines / 4))
	echo "File: $file - Number of sequences: $sequences"
	fi
done

# Rename and move Dorado output files
for fq in ./DoradoDemux/*.fastq; do 
	# Extract barcode ID from filename
	barcode=$(basename "$fq" .fastq)
	# Look up marker type
	marker=$(awk -v bc="$barcode" '$1==bc {print $2}' "$BARCODE_FILE")

	if [[ -z "$marker" ]]; then 
		echo "Warning: barcode $barcode not found in barcode map. Skipping..."
		continue
	fi
	newname="${barcode}_${marker}.QC.fastq"
	mv "$fq" "./DoradoDemux/$marker/$newname"
done

# VSEARCH orienting
# Assumes that orienting 'database' is existing and set in "VSEARCH_DB"
# Check that VSEARCH version is correct 
echo "VSEARCH version: $(vsearch --version)"
echo "If VSEARCH orienting throws an error, please check that your version is >2.16"
for marker in 16S ITS; do
	for w in ./DoradoDemux/$marker/*.fastq; do
		base=$(basename "$w" .fastq)
		vsearch --orient "$w" --db "$VSEARCH_DB_${marker}" --fastqout ./MinibarDemux/"${base}.oriented.fastq"
	done
done

# Minibar demultiplexing by barcoded primer
# Assumes that demultiplexing files are existing and set in DEMUX_16S and DEMUX_ITS respectively.
# Check minibar is installed, wget if not
if [[ ! -x "$MINIBAR_PATH" ]]; then
	echo "Minibar script not found at $MINIBAR_PATH, downloading..."
	wget -O "$MINIBAR_PATH" https://raw.githubusercontent.com/calacademy-research/minibar/master/minibar.py
	chmod +x "$MINIBAR_PATH"
else 
echo "Found minibar script at $MINIBAR_PATH"
fi

# Check input files
echo "Checking demux file $DEMUX_16S"
"$MINIBAR_PATH" -info cols "$DEMUX_16S"
echo "Checking demux file $DEMUX_ITS"
"$MINIBAR_PATH" -info cols "$DEMUX_ITS"

echo "This script assumes your primers can be demultiplexed with edit distance of 4. Please exit if this is not correct and edit the script."

for marker in 16S ITS; do
	for fq in ./DoradoDemux/$marker/*.fastq; do
		echo "Demultiplexing $fq with Minibar for $marker..."
		"$MINIBAR_PATH" "$DEMUX_${marker}" -e 4 -E 8 -l 120 -F "$fq"
		mv sample_*.fastq ./MinibarDemux/"$marker"/
	done
done  

echo "QC and demux finished. Please check output files for length and sequence/sample distributions."
