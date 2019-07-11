# ===============================
# = simulate genome, phenotypes =
# ===============================

# 1. macs simulation
# ============================================================

for i in `seq 1 1000`
do
  sbatch --export=i="$i" --output log/macs."$i".out --error log/macs."$i".err macs.sbatch
done

# 2. make tfam file
# ============================================================

seq 1 75000 | awk '{print $0" "$0" 0 0 2 -9"}' > plink/plink.tfam

# 3. parse
# ============================================================

for i in `seq 1 1000`
do
  sbatch --export=i="$i" --output log/macs2plink.chr."$i".out --error log/macs2plink.chr."$i".err macs2plink.sbatch
done

# some could not be completed

for i in `seq 1 1000`
do
	count=`grep ^Error plink/chr"$i".log | wc -l`
	if [ $count -gt 0 ]
	then
    sbatch --export=i="$i" --output log/macs2plink.chr."$i".out --error log/macs2plink.chr."$i".err macs2plink.sbatch
	fi
done

# 4. merge
# ============================================================

# create merge list
seq 1 1000 | awk '{print "plink/chr"$1}' > plink/bed.file.list
~/qgg/software/plink-v1.90b5.3/plink --merge-list plink/bed.file.list --make-bed --out plink/all > plink/merge.all.log 2>&1 &

# 5. look at allele frequency
# ============================================================

# the macs set up was such that pop 2, 3 more closely related to each other
# use label C for the one that split from the rest 2,000 gen ago
# then A and B are interchangable

seq 1 25000 | awk '{print $0" "$0}' > plink/popC.id.list
seq 25001 50000 | awk '{print $0" "$0}' > plink/popA.id.list
seq 50001 75000 | awk '{print $0" "$0}' > plink/popB.id.list

~/qgg/software/plink-v1.90b5.3/plink --bfile plink/all --keep plink/popA.id.list --freq --out plink/popA > plink/popA.freq.log 2>&1 &
~/qgg/software/plink-v1.90b5.3/plink --bfile plink/all --keep plink/popB.id.list --freq --out plink/popB > plink/popB.freq.log 2>&1 &
~/qgg/software/plink-v1.90b5.3/plink --bfile plink/all --keep plink/popC.id.list --freq --out plink/popC > plink/popC.freq.log 2>&1 &

# only use common variants for all fitting etc..
# ============================================================

paste plink/popA.frq plink/popB.frq plink/popC.frq | awk '$5 >= 0.01 || $11 >= 0.01 || $17 >= 0.01 {print $2}' | tail -n+2 > plink/common.snp
~/qgg/software/plink-v1.90b5.3/plink --bfile plink/all --extract plink/common.snp --make-bed --out plink/all.common > plink/all.common.log 2>&1 &

# 6. PCA
# ============================================================

# sample variants for pca
shuf plink/all.bim | head -n 100000 | cut -f 2 > plink/random.snp.list

srun --cpus-per-task=16 --time=72:00:00 --mem=128G  ~/qgg/software/plink-v1.90b5.3/plink --bfile plink/all --extract plink/random.snp.list --pca --out plink/all.pca > plink/all.pca.log 2>&1 &

# 7. simulation
# ============================================================

for i in `seq 1 20`
do
  sbatch --export=i="$i",H2=0.8,all=plink/all.common,ncpu=16 --output=log/sim"$i".out --error=log/sim"$i".err h2.sbatch
done

# one of the replicate could not finish
i = 21
sbatch --export=i="$i",H2=0.8,all=plink/all.common,ncpu=16 --output=log/sim"$i".out --error=log/sim"$i".err h2.sbatch
# after it's done
rm -r /mnt/scratch/huangw53/sim/reps/rep3
mv /mnt/scratch/huangw53/sim/reps/rep21 /mnt/scratch/huangw53/sim/reps/rep3


# 8. summarize data
# ============================================================

/mnt/research/qgg/software/R-3.6.0/bin/Rscript summarizeData.R 20 /mnt/scratch/huangw53/sim/reps/rep reportData/summary.RData


