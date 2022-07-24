thread=1
dir="."
name="sample"
readsidlist=("idlist1" "idlist2")
suffix=".fastq.gz"
script="" #script_dir

str=""
echo -n ""> ${dir}/coverage_cutoff_file
echo -n ""> ${dir}/kmc_${name}
for i in `seq 0 $[${#readsidlist[*]}-1]`
do
	sh $script/1.trim $thread $dir ${readsidlist[$i]} $suffix
	sh $script/2.kmc_db $thread $dir ${name}${i} ${readsidlist[$i]} $suffix
	coverageL=$(PloidyFrost cutoffL ${dir}/hist_${name}${i})
	coverageU=$(PloidyFrost cutoffU ${dir}/hist_${name}${i})
	echo -e "$coverageL\t$coverageU" >>${dir}/coverage_cutoff_file
	echo "${dir}/kmc_${name}${i}">>${dir}/kmc_${name}
	sh $script/3.filter $thread $dir ${name}${i} $coverageL
	str=${str}" -r ${name}${i}_filtered.fq "
	
done

Bifrost build -i -d -k 25 -v ${str} -c -o ${name} -t ${thread} 1>${dir}/log_bifrost 2>${dir}/err_bifrost


PloidyFrost -g ${dir}/${name}.gfa -f ${dir}/${name}.bfg_colors -d ${dir}/kmc_${name} -v -o ${name} -C ${dir}/coverage_cutoff_file 1>${dir}/log_ploidy 2>${dir}/err_ploidy

PloidyFrost model -g ${dir}/PloidyFrost_output/${name}_allele_frequency.txt -o ${dir}/${name}

