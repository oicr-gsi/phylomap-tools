### need a directory with consensus sequences generated from the pipeline
module load ncov-tools
module load bis-rlibs
module load phylo-tools
module load sars-covid-2/mn908947.3 


refid=MN908947.3
dir_results=results
dir_fasta=consensus
prefix=OICR
rename=rename.txt
REF=$SARS_COVID_2_ROOT/MN908947.3.fasta
plot_title=SHSC_2021-01-11
plot_prefix=SHSC_2021-01-11

echo "processing all"
mkdir -p $dir_results
echo "collecting and filtering consensus sequences from $dir_results"
python preprocess_consensus.py --directory $dir_fasta -o $dir_results -t $rename -x exclude.txt

echo "augur alignments and tree generation"
augur align --sequences $dir_results/consensus.fasta --reference-sequence $REF --output $dir_results/aligned.fasta --fill-gaps
augur tree --alignment $dir_results/aligned.fasta --output $dir_results/tree_raw.nwk
nw_reroot $dir_results/tree_raw.nwk $refid > $dir_results/tree.nwk

echo "get variants/alleles"
python align2alleles.py --reference-name $refid  $dir_results/aligned.fasta > $dir_results/alleles.tsv

Rscript --vanilla variant_maps.R --directory $dir --prefix $plot_prefix --tree tree.nwk  -o $dir --title $plot_title --metadata ../metadata.txt --covariates RunDate  






