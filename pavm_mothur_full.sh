#! /bin/bash
# antoine.bridier-nahmias@inserm.fr
# 27/11/2018
# Mothur
echo "ARGS" 
echo "$@"
mothur="$HOME/data/tools/mothur_1.39.5/mothur"
cd $(dirname $0)
echo "path is" $(pwd)
# Step counter init
stepus=1
n_threadus=4

if [ $# == 0 ]; then 
	echo "--groups : Make the reads group file"
	echo "--contigs : Make contigs"
	echo "--screen_1 : Make screening for size and ambiguous positions"
	echo "--unique : Make the sequences unique and generate count file"
	echo "--silva_pcr : Make the DB processing (virtual PCR or trim)"
	echo "--align : Making the alignments to the db"
	echo "--clean-align : Make the alignments clean"
	echo "--pre-clust : Making the pre-clustering and chimera slaying"
	echo "--classify : Make classification"
	echo "--matrix-clust : Make matrix OTU clustering and count"
  echo "--tree : Make a tree"
fi

prefixus="N_pavm_" # Defined for the whole script
################################################################################
################################# Make groups ##################################
################################################################################
read_dir="../data/reads/${prefixus}"
echo "STEP : $stepus"
if [[ "$@" =~ "--groups" ]]; then
	startus=$(date +%s)
	echo "--groups : Make the reads group file"
	$mothur "#make.file(
	            inputdir=${read_dir}, 
	            type=gz, 
		    	numcols = 2,
	            prefix=${prefixus})"
	# Adding the sample name as the first col in predires.files 
	# and changing - to _ in the prefix in order to please the mothur gods
	sed -r -i 's:(^.*/(.*)_S.*_L0.*fastq.gz):\2 \1:g' \
	    "${read_dir}/${prefixus}.files"
	sed -r -i 's/^([0-9]+)-/\1_/' "${read_dir}/${prefixus}.files"
	endus=$(date +%s)
	elapsus=$((endus - startus))
	echo "Time elapsed in s"
	echo "$elapsus"
else
	echo "NO --groups : Make the reads group file"
fi
((stepus++))


################################################################################
################################ Make contigs ##################################
################################################################################
contigs_dir="../data/${prefixus}contigs"
echo "STEP : $stepus"
if [[ "$@" =~ "--contigs" ]]; then
	startus=$(date +%s)
	echo "--contigs : Make contigs"
	if [ ! -d $contigs_dir ]; then mkdir $contigs_dir ; fi
	# V4 region amplified (<300pb) so trimoverlap = T 
	$mothur "#set.dir(output=${contigs_dir});
	        make.contigs(
	            file=${read_dir}/${prefixus}.files,
	            trimoverlap=T,
	            processors=${n_threadus});
	        summary.seqs()"
	endus=$(date +%s)
	elapsus=$((endus - startus))
	echo "Time elapsed in s"
	echo "$elapsus"
else
	echo "NO --contigs : Make contigs"
fi
((stepus++))


################################################################################
############################# Make first screen ################################
################################################################################
echo "STEP : $stepus"
if [[ "$@" =~ "--screen_1" ]]; then
	startus=$(date +%s)
	echo "--screen_1 : Make screening for size and ambiguous positions"
	$mothur "#screen.seqs(
	            fasta=${contigs_dir}/${prefixus}.trim.contigs.fasta,
	            group=${contigs_dir}/${prefixus}.contigs.groups,
	            summary=${contigs_dir}/${prefixus}.trim.contigs.summary,
	            maxambig=0,
	            minlength=170,
	            maxlength=265,
		    maxhomop=5,
		    processors=${n_threadus});
	            summary.seqs()"
	endus=$(date +%s)
	elapsus=$((endus - startus))
	echo "Time elapsed in s"
	echo "$elapsus"
else
	echo "NO --screen_1 : Make screening for size and ambiguous positions"
fi
((stepus++))

################################################################################
############################## Make unique #####################################
################################################################################
	echo "STEP : $stepus"
if [[ "$@" =~ "--unique" ]]; then
	startus=$(date +%s)
	echo "--unique : Make the sequences unique and generate count file"
	$mothur "#unique.seqs(
	            fasta=${contigs_dir}/${prefixus}.trim.contigs.good.fasta);
	        count.seqs(
	            name=${contigs_dir}/${prefixus}.trim.contigs.good.names,
	            group=${contigs_dir}/${prefixus}.contigs.good.groups,
	            processors=${n_threadus});
	        summary.seqs(count=${contigs_dir}/${prefixus}.trim.contigs.good.count_table)"
	endus=$(date +%s)   
	elapsus=$((endus - startus))
	echo "Time elapsed in s"
	echo "$elapsus"
else	
	echo "NO --unique : Make the sequences unique and generate count file"
fi
((stepus++))

################################################################################
######################### Make Silva customization #############################
################################################################################
silva_dir="../data/silva_db"
echo "STEP : $stepus"
if [[ "$@" =~ "--silva_pcr" ]]; then
	silva_db="silva.nr_v132.align"
	echo "--silva_pcr : Make the DB processing (virtual PCR or trim)"
	# Start and end can be determind by running the alignement first on the full db
	# long and painful
	$mothur "#set.dir(output=${silva_dir});
		summary.seqs(fasta=${silva_dir}/${silva_db}, processors=${n_threadus});
	        pcr.seqs(
	            fasta=${silva_dir}/${silva_db}, 
	            start=11895,
	            end=25318,
	            keepdots=F, 
	            processors=${n_threadus});
	        summary.seqs(processors=${n_threadus})"
else
	echo "NO --silva_pcr : Make the DB processing (virtual PCR)"
fi
((stepus++))

################################################################################
######################### Make alignment to Silva ##############################
################################################################################
align_dir="../data/${prefixus}align"
echo "STEP : $stepus"
if [[ "$@" =~ "--align" ]]; then
	startus=$(date +%s)
	silva_db="silva.nr_v132.align"
	echo "--align : Making the alignments to the db"
	if [ ! -d $align_dir ]; then mkdir $align_dir ; fi
	$mothur "#set.dir(output=${align_dir});
	        align.seqs(
	            fasta=${contigs_dir}/${prefixus}.trim.contigs.good.unique.fasta,
	            reference=${silva_dir}/${silva_db},
 		    flip=T,	
	            processors=${n_threadus});
		summary.seqs(fasta=${align_dir}/${prefixus}.trim.contigs.good.unique.align,
			count=${contigs_dir}/${prefixus}.trim.contigs.good.count_table,
			processors=${n_threadus})"	
	endus=$(date +%s)
	elapsus=$((endus - startus))
	echo "Time elapsed in s"
	echo "$elapsus"
else
	echo "NO --align : Making the alignments to the db"
fi
((stepus++))


################################################################################
########################## Make alignment clean ################################
################################################################################
if [[ "$@" =~ "--clean-align" ]]; then
	echo "STEP : $stepus"
	echo "--clean-align : Make the alignments clean"
	$mothur "#screen.seqs(
	            fasta=${align_dir}/${prefixus}.trim.contigs.good.unique.align,
	            count=${contigs_dir}/${prefixus}.trim.contigs.good.count_table,
	            summary=${align_dir}/${prefixus}.trim.contigs.good.unique.summary,
		    optimize=start-end,
		    criteria=90,	
		    maxambig=0,
		    maxhomop=6,
	            processors=${n_threadus});
	        filter.seqs(fasta=${align_dir}/${prefixus}.trim.contigs.good.unique.good.align,
	            vertical=T, 
	            trump=.,
	            processors = ${n_threadus});
	        summary.seqs(
	            fasta=${align_dir}/${prefixus}.trim.contigs.good.unique.good.filter.fasta,
	            count=${align_dir}/${prefixus}.trim.contigs.good.good.count_table,
	            processors = ${n_threadus})"
else
	echo "STEP : $stepus"
	echo "NO --clean-align : Make the alignments clean"
fi
((stepus++))

################################################################################
#################### Make pre-cluster and chimera slaying ######################
################################################################################
pre_clust_dir="../data/${prefixus}pre_clust"
echo "STEP : $stepus"
if [[ "$@" =~ "--pre-clust" ]]; then
	startus=$(date +%s)
	if [ ! -d $pre_clust_dir ]; then mkdir $pre_clust_dir ; fi
	echo "--pre-clust : Making the pre-clustering and chimera slaying"
	# diffs=3 means 6 dist max within the pre-cluster meaning 6*100/114(align size)=2.80% diff within the pre-cluster 
	#this value is < to our 3% threshold for OTU clustering.
	$mothur "#set.dir(output=${pre_clust_dir});
	        pre.cluster(
	            fasta=${align_dir}/${prefixus}.trim.contigs.good.unique.good.filter.fasta,
	            count=${align_dir}/${prefixus}.trim.contigs.good.good.count_table,
	            diffs=3,
	            processors=${n_threadus});
	        summary.seqs(count=${pre_clust_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.count_table,
	        	processors=${n_threadus});
	        chimera.vsearch(
	            fasta=${pre_clust_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.fasta,
	            count=${pre_clust_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.count_table,
	            dereplicate=t,
	            processors=${n_threadus});
	        remove.seqs(
	            fasta=${pre_clust_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.fasta,
	            count=${pre_clust_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.count_table,
	            accnos=${pre_clust_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.accnos);
	        summary.seqs(count=${pre_clust_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.pick.count_table,
	        	processors=${n_threadus})"
	endus=$(date +%s)   
	elapsus=$((endus - startus))
	echo "Time elapsed in s"
	echo "$elapsus"
else
	echo "NO --pre-clust : Making the pre-clustering and chimera slaying"
fi
((stepus++))

################################################################################
########################### Make classification ################################
################################################################################
class_dir="../data/${prefixus}classification"
silva_fasta="../data/silva_db/silva.nr_v132.align"
silva_tax="../data/silva_db/silva.nr_v132.tax"
echo "STEP : $stepus"
if [[ "$@" =~ "--classify" ]]; then
	startus=$(date +%s)
	if [ ! -d $class_dir ]; then mkdir $class_dir ; fi
	echo "--classify : Make classification"
	$mothur "#set.dir(output=${class_dir});
			classify.seqs(
				fasta=${pre_clust_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.fasta,
				count=${pre_clust_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.count_table,
				reference=${silva_fasta},
				taxonomy=${silva_tax},
				cutoff=80,
				processors=${n_threadus});
			remove.lineage(
				fasta=${pre_clust_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.fasta,
				count=${pre_clust_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.count_table,
				taxonomy=${class_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.nr_v132.wang.taxonomy,
				taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota')"
	endus=$(date +%s)   
	elapsus=$((endus - startus))
	echo "Time elapsed in s"
	echo "$elapsus"
else	
	echo "NO --classify : Make classification"
fi
((stepus++))

################################################################################
###################### Make Matrix OTU cluster and count #######################
################################################################################
echo "STEP : $stepus"
mat_otu_dir="../data/${prefixus}mat_cluster_otu"
if [[ "$@" =~ "--matrix-clust" ]]; then
    startus=$(date +%s)
	if [ ! -d $mat_otu_dir ]; then mkdir $mat_otu_dir ; fi
	echo "--matrix-clust : Make matrix OTU clustering and count"
	# cutoff for dist seq allows not to store too much info
	# the dist file can be big and has to fit in ram once per processor used... careful!
	$mothur "#set.dir(output=${mat_otu_dir});
	        dist.seqs(
	            fasta=${class_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.pick.fasta,
	            cutoff=0.03,
	            processors=${n_threadus});        
	        cluster(
	            column=${mat_otu_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.pick.dist,
	            count=${class_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.pick.count_table,
	            cutoff=0.03);
	        make.shared(
	            list=${mat_otu_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.list,
	            count=${class_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.pick.count_table,
	            label=0.03)"
    endus=$(date +%s)
    elapsus=$((endus - startus))
    echo "Time elapsed in s"
    echo "$elapsus"
else
	echo "NO --matrix-clust : Make matrix OTU clustering and count"	
fi
((stepus++))


################################################################################
###################### Make taxonomic OTU cluster and count ####################
################################################################################
tax_otu_dir="../data/${prefixus}tax_cluster_otu"
echo "STEP : $stepus"
if [ "$1" == 1 ]; then
	startus=$(date +%s)
	if [ ! -d $tax_otu_dir ]; then mkdir $tax_otu_dir ; fi
	echo "Make taxonomic OTU clustering and count"
	$mothur "#set.dir(output=${tax_otu_dir});
			cluster.split(
				fasta=${pre_clust_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.fasta,
				count=${pre_clust_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.count_table,
				splitmethod=classify,
				taxonomy=${class_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.nr_v132.wang.taxonomy,
				taxlevel=6,
				cutoff=0.03,
				processors=${n_threadus});
			make.shared(
		            list=.${tax_otu_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.unique_list.list,
		            count=${pre_clust_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.count_table,
		            label=0.03)"
	endus=$(date +%s)
	elapsus=$((endus - startus))
	echo "Time elapsed in s"
	echo "$elapsus"
else	
	echo "NO Make taxonomic OTU clustering and count"
fi
((stepus++))



################################################################################
########################## Make shared file clean ##################################
################################################################################
# In this step, we are going to clean the shared file by removing some samples and some OTU
# using the filter.shared() and remove.groups() commands.
clean_dir="../data/${prefixus}clean_shared"
echo "STEP : $stepus"
if [ "$1" == 1 ]; then
    startus=$(date +%s)
    if [ ! -d $clean_dir ]; then mkdir $clean_dir ; fi
    echo "Make shared file clean"
    $mothur "#set.dir(output=${clean_dir});
            remove.groups(
            	shared=${tax_otu_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.opti_mcc.unique_list.shared,
            	count=${pre_clust_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.count_table,
            	column=${mat_otu_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.dist,
            	list=${mat_otu_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.opti_mcc.list,
            	accnos=${clean_dir}/bad_groups.accnosgroups);
            filter.shared(
            	shared=${clean_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.opti_mcc.unique_list.0.03.pick.shared,
            	mintotal=3)"
    endus=$(date +%s)
    elapsus=$((endus - startus))
    echo "Time elapsed in s"
    echo "$elapsus"
else    
    echo "NO Make shared file clean"
fi
((stepus++))

################################################################################
############################ Make alpha div ####################################
################################################################################
alpha_dir="../data/${prefixus}alpha"
echo "STEP : $stepus"
if [[ "$@" =~ "--alpha" ]]; then
    startus=$(date +%s)
    if [ ! -d $alpha_dir ]; then mkdir $alpha_dir ; fi
    echo "Make alpha diversity calculations"
    $mothur "#set.dir(output=${alpha_dir});
            rarefaction.single(
                shared=${mat_otu_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.shared,
                calc=sobs-chao-ace-shannon-coverage-nseqs,
                iters=1000,
                freq=1000,
                processors=${n_threadus}
            )"
    endus=$(date +%s)
    elapsus=$((endus - startus))
    echo "Time elapsed in s"
    echo "$elapsus"
else    
    echo "NO Make alpha diversity calculations"
fi
((stepus++))



################################################################################
#########################     Make a tree     ##################################
################################################################################
# The command tree.shared() makes a tree of samples with distances being calculated by BC or TYC
# in our case we need a tree with OTUs as tips to feed to the Unifrac algorithm, we'll try the clearcut command.
# We'll start with get.oturep()
tree_dir="../data/${prefixus}tree"
echo "STEP : $stepus"
if [[ "$@" =~ "--tree" ]]; then
    startus=$(date +%s)
    if [ ! -d $clean_dir ]; then mkdir $clean_dir ; fi
    echo "--tree : Make a tree"
    $mothur "#set.dir(output=${tree_dir});
  			get.oturep(
	            		fasta=${class_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.pick.fasta,
				column=${mat_otu_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.pick.dist,
  				list=${mat_otu_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.list,
  				count=${class_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.pick.count_table);
			clearcut(
				fasta=${tree_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.0.03.rep.fasta,
				DNA=T,
				verbose=T)"
    endus=$(date +%s)
    elapsus=$((endus - startus))
    echo "Time elapsed in s"
    echo "$elapsus"
else    
    echo "NO --tree : Make a tree"
fi
((stepus++))

################################################################################
############################ Make beta div #####################################
################################################################################
echo $@
beta_dir="../data/${prefixus}beta_1000"
echo "STEP : $stepus"
if [[ "$@" =~ "--beta" ]]; then
    startus=$(date +%s)
    if [ ! -d $beta_dir ]; then mkdir $beta_dir ; fi
    echo "Make beta diversity calculations"
    $mothur "#set.dir(output=${beta_dir});
            unifrac.weighted(
		tree=${tree_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.0.03.rep.tre,
		count=${tree_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.0.03.rep.count_table,
		subsample=50000,
		iters=100,
		distance=square,
		processors=${n_threadus});
            unifrac.unweighted(
		tree=${tree_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.0.03.rep.tre,
		count=${tree_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.0.03.rep.count_table,
		subsample=50000,
		iters=100,
		distance=square,
		processors=${n_threadus});
	    dist.shared(
		shared=${mat_otu_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.shared, 
		calc=thetayc-braycurtis,
		output=square,
		subsample=50000,
		iters=100,
		processors=${n_threadus})"
    endus=$(date +%s)
    elapsus=$((endus - startus))
    echo "Time elapsed in s"
    echo "$elapsus"
else    
    echo "NO Make beta diversity calculations"
fi
((stepus++))

################################################################################
###################     Make Taxonomic Assignation    ##########################
################################################################################
assign_dir="../data/${prefixus}assign"
echo "STEP : $stepus"
if [[ "$@" =~ "--assign" ]]; then
    startus=$(date +%s)
    if [ ! -d $assign_dir ]; then mkdir $assign_dir ; fi
    echo "Make Taxonomic Assignation"
    $mothur "#set.dir(output=${assign_dir});
	classify.otu(
	taxonomy=${class_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.nr_v132.wang.taxonomy,
 	list=${mat_otu_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.list,
  	count=${class_dir}/${prefixus}.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.pick.count_table,
	cutoff=51)"
    endus=$(date +%s)
    elapsus=$((endus - startus))
    echo "Time elapsed in s"
    echo "$elapsus"
else    
    echo "NO Make Taxonomic Assignation"
fi
((stepus++))