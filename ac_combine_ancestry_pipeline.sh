#!/bin/bash
#
# Time needs to be MUCH LONGER

########## INPUTS ##########
# bamdir="/moto/ziab/yk2840/Frog_workspace/"
# parallelc = 10
# varlist=""
# outdir=""
# F2="TRUE"



########### STEP 1: OBTAIN GENOTYPE #############
module load anaconda
source activate /rigel/zi/users/yk2840/software/conda/platypus_specific


# Here, write a function that splits the input bamlist into 12 chunks and parallel submit.

# Make Directories for the allelecount
	mkdir $outdir/ac/
	ls -l $bamdir > $outdir/ac/bam_list.txt

# Split Bam list into 6 chunks
	split --additional-suffix split -n l/6 $outdir/ac/bam_list.txt

# Genotype Samples
	i=0
	for bams in $outdir/ac/bam_list*split*.txt; do
		sbatch ac_sample.sh $outdir/ac/$bams $outdir/ac/ $varlist  
		i=i+1
	done

# Add function to wait until this job finishes using i


############ Step 2: Combine Sample Genotypes for Ancestry HMM #########

module load R
Rscript ./ancestry_combine.R $outdir/ac/ ./ancestry_panel.RData ./fixed_pet_variants_allelecounter.RData



############ Step 3: Remove All "Missing F2" Data if you're getting memory problems ###############

# awk '{s=0; for (i=7;i<=NF;i++) s+=$i; if (s!=0)print}' your sample panel  # Dealing with 0s in the sample regions
#awk '{ for(i=1; i<=NF;i++) j+=$i; print j; j=0 }' your sample panel for indexing and recalculating the cM distance in the ancestral panel output


########### Step 4: Use Interactive script for now ######################
# Keeps crashing, so I don't think in the long-run an interactive script is a good idea
# need 5060000 of mem to get to forward-backward posterior decoding




###########


############ Step 3: Run Ancestry HMM ########
module load ancestry_hmm/0.94

	if [ F2 == "TRUE" ]; then
		ancestry_hmm \
			-i ./2019-07-17_192444_ac_panel.txt \
			-s ./2019-07-17_192530_ac_ploidy.txt \
			-a 2 .5 .5 \
			-p 0 2 \
			-p 1 2 \
			--fix
	else
		echo "Rerun Ancestry for particular ancestry mix"
	fi

############ Step 3: Run R/QTL? ########


date
