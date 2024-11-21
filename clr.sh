#用SWEED计算CLR
#/data01/wangyf/software/sweed/SweeD -h

#拆分成50条染色体的vcf
source activate vcftools

for i in {572..621}; do
vcftools --gzvcf /data01/wangyf/project2/CyprinusCarpio/15.pop/15.south-irtysh/0.select-sample/nor-irty-sou/irtysh.vcf.gz  --recode --recode-INFO-all --stdout --chr NC_056${i}.1 > NC_056${i}.1.vcf;
done


#将每个chr拆成100kb的窗口，也就是length/100000，计算grid参数大小，保存到文件chr_gird.txt中
NC_056572.1,396
NC_056573.1,285
NC_056574.1,394
NC_056575.1,363
NC_056576.1,400
NC_056577.1,342
NC_056578.1,473
......

#运行SWEED
for line in `cat chr_grid.txt`
do
  chr=${line%%,*}
  grid=${line##*,}
  /data01/wangyf/software/sweed/SweeD -name chr${chr}.100kb \
        -input ${chr}.vcf \
        -grid ${grid} \
        -minsnps 40 \
        -maf 0.05 \
        -missing 0.1
sed -i '1,3d' SweeD_Report.chr${chr}.100kb #删除每个report文件的前3行,"Position Likelihood Alpha StartPos EndPos"
done

#添加染色体号
for line in `cat chr_grid.txt`
do
  chr=${line%%,*}
  awk -v chr="$chr" '{print chr,$4,$5,$3}' SweeD_Report.chr${chr}.100kb > SweeD_Report.chr${chr}.100kb-2
  sed -i "s@ @\t@g" SweeD_Report.chr${chr}.100kb-2
done

#合并
cat *-2 >total_chr.txt
awk '{sub(/\.0000/, "", $2); sub(/\.0000/, "", $3); print}' total_chr.txt  > total_chr.txt-1
sed -i "s@ @\t@g" total_chr.txt-1

#R可视化
conda activate r4.3

chr1 <- read.table(file="SweeD_Report.NC_056572.1.100kb",sep = "\t", header = T)
plot(chr1[,1]/10000,chr1[,2])
