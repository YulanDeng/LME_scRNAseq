## Script for the WES analysis, neoantigen recognition and clonal evolution
##Yulan Deng, last updated 2023-4-20
##my e-mail:kndeajs@163.com

####################################
#identification of somatic mutation#
####################################
#${hg38RefFile} is the reference file for hg38
#${fastqFile1} is the fastq file for read 1
#${fastqFile2} is the fastq file for read 2
#${smp} is the sampleID
#${outputDir} is the output directory 
#${targetFile} is the bed file of target region for WES
#${cancerBamFile} is the process bam file for cancer
#${normalBamFile} is the process bam file for normal
#${humandbDir} is the directory of humandb of ANNOVAR

#align by BWA
bwa mem -o ${smp}.sam -t 32 -M ${hg38RefFile} ${fastqFile1} ${fastqFile2}
samtools view -@ 16 -bS ${smp}.sam > ${smp}.bam
samtools sort -@ 16 -o ${smp}.sort.bam ${smp}.bam
samtools index ${smp}.sort.bam
rm ${smp}.bam ${smp}.sam

#quality control by qualimap
qualimap --java-mem-size=35G bamqc -bam ${smp}.sort.bam -outdir ${outputDir} -nt 60 \
-outformat PDF -gff ${targetFile}

#markduplicated by picard
java -jar picard.jar MarkDuplicates I=${smp}.sort.bam O=${smp}.sort.dedup.bam \
CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=. \
METRICS_FILE=${smp}.sort.metrics ASSUME_SORTED=true MAX_RECORDS_IN_RAM=8000000

#BaseRecalibrator by GATK
gatk BaseRecalibrator -I ${smp}.sort.dedup.bam -R ${hg38RefFile} \
--known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
--known-sites dbsnp_146.hg38.vcf.gz \
--known-sites 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
--output ${smp}_recalibration_report.grp;  
gatk ApplyBQSR -R ${hg38RefFile} -I ${smp}.sort.dedup.bam \
--bqsr-recal-file ${smp}_recalibration_report.grp \
--output ${smp}.sort.dedup.baserecal.bam

#call somatic mutation by mutect2
gatk Mutect2 -R ${hg38RefFile} -I ${cancerBamFile} -I ${normalBamFile} \
--tumor-sample ${smp}_C -normal ${smp}_A \
--TMP_DIR . --germline-resource af-only-gnomad.hg38.vcf.gz \
--output ${smp}.vcf

#filtation the mutation results by FilterMutectCalls and VariantFiltration
gatk FilterMutectCalls -V ${smp}.vcf -O ${smp}.filter.vcf; 
gatk VariantFiltration -R ${hg38RefFile} \
-V ${smp}.filter.vcf \
-O ${smp}.filter1.vcf \
--filter-name "alternativeDepth" \
--filter-expression "vc.getGenotype(${smp}_C).getAD().1 < 6" \
--filter-name "depth" \
--filter-expression "DP < 30" 

#annoation the somatic mutation by ANNOVAR
grep -v "^#" ${smp}.filter1.vcf | \
awk '($7=="PASS" ){print $1 "\t" $2 "\t" $2 "\t" $4 "\t" $5}' > ${smp}.for.anno1;
table_annovar.pl ${smp}.for.anno1 ${humandbDir} \
-buildver hg38 -out  ${smp} -remove -protocol refGene,cosmic70 -operation g,f -nastring .

###################################
#generate peptide from somatic SNP#
###################################
#${smp} is sample ID
#${inputDir} is the directory of ANNOVAR results
#${chrDir} is the directory of output
#${smp2} is sample ID for matched normal sample
#${refseqFile} is the annotation of amino acids from refSeq
#${refseqID} is the file of gene ID of refSeq

#Load required packages
library(Biostrings)

#Extract the information of amino acids
AAstring <- readAAStringSet(filepath=refseqFile)
nm2np <- read.table(file =refseqID,sep="\t",quote ="",stringsAsFactors=F,fill=T,header=F)
nm2np <- nm2np[nm2np[,2]=="Reviewed",]
nm2np[,1] <- sapply(strsplit(nm2np[,1],".",fixed=T),function(x) x[1])
nm2np <- nm2np[nm2np[,6]%in%names(AAstring),]

#Load the result of ANNOVAR
flnm <- paste0(inputDir,smp,".hg38_multianno.txt")
fl <- read.table(file=flnm,sep="\t",quote = "",stringsAsFactors=F,header=T)
fl <- cbind(fl,rep(smp,nrow(fl)))
fl <- fl[fl[,6]%in%c(  "exonic","exonic;splicing"),]
fl <- fl[fl[,9]=="nonsynonymous SNV",]
colnames(fl)[12] <- "Tumor_Sample"

#make the amino acids information for each somatic mutation
nm_list <- strsplit(unlist(strsplit(fl[,"AAChange.refGene"],",")),":")
aapos <- sapply(nm_list,function(x) x[5])
gene <- sapply(nm_list,function(x) x[1])
nm_df <- data.frame(nm=sapply(nm_list,function(x) x[2]),
pos=sapply(aapos,function(x) as.numeric(sub(".","",gsub("p|[A-Z]","",x)),fixed=T)),
mut=sapply(aapos,function(x) substr(x,nchar(x),nchar(x))),
gene=gene,aapos=aapos,
stringsAsFactors=F)

#generate peptide from somatic SNP
nm_df <- nm_df[nm_df[,1]%in%nm2np[,1],]
np <- nm2np[match(nm_df[,1],nm2np[,1]),6]
nm_df <- cbind(nm_df,np)
seq_list <- lapply(1:nrow(nm_df),function(x){
		aa=nm_df[x,"np"]
		pos=nm_df[x,"pos"]
		mut=nm_df[x,"mut"]
		seq_fun2 <- function(AAstringt,aat,post,lent,mutt){			
			AAstringlist=unlist(strsplit(AAstringt,""))
			AAstringlist[post]=mutt
			AAstringt = paste(AAstringlist,collapse="")			
			if(lent > nchar(AAstringt))
			{
				return(c())
			}else{
				if(post < lent)
				{
					aaMin=1
				}else{
					aaMin= post-lent
				}
				
				if(post + lent -1 > nchar(AAstringt))
				{
					aaMax=nchar(AAstringt)
				}else{
					aaMax= post+lent-1
				}				
				seqRes <- substr(AAstringt,aaMin,aaMax)
				return(seqRes)
			}
			
		}		
		cut_list <- seq_fun2(AAstringt=as.character(AAstring[aa][[1]][1:(length(AAstring[aa][[1]]))]),
			aat=aa,post=pos,lent=14,mutt=mut)
	return(cut_list)
})
seq_v <- unlist(seq_list)
seq_label <- paste0(">",paste(nm_df[,1],nm_df[,4],nm_df[,5],sep="_")[sapply(seq_list,length)==1]) 
seq_df <- as.character(t(matrix(c(seq_label,seq_v),ncol=2)))

#write the result
setwd(paste0(chrDir,smp2))

#######################################################################
#candidate peptides were filtered by proteasomal cleavage from NetChop#
#######################################################################
#${smp} is sample ID
#${chrDir} is the directory of output
#${smp2} is sample ID for matched normal sample
#run NetChop
netchop ${smp}.for.NetChop.pep > ${smp}.NetChop.myout

#filter candidate peptides by NetChop
setwd(paste0(chrDir,smp2))
NetChopOut <- read.delim(file =paste0(smp,".NetChop.myout"),sep="\n",quote ="",stringsAsFactors=F,fill=T,skip = 17,header=F)
Netfasta <- read.table(file =paste0(smp,".for.NetChop.pep"),sep="\n",quote ="",stringsAsFactors=F,fill=T,header=F)
Npep <- nrow(Netfasta)/2
lenpep <- sapply(seq(Npep),function(x) nchar(Netfasta[2*x,1]))
lenpep_extend <- lenpep + 6
nstart <- c(0,cumsum(lenpep_extend[-length(lenpep_extend)]))
headerpep <- sapply(seq(Npep),function(x) Netfasta[2*x-1,1])
headerpep_s <- sapply(strsplit(headerpep,"_"),function(x) x[4])
headerpep_sn <- as.numeric(gsub("[A-Z]","",sub("p.","",headerpep_s,fixed=T)))
posM <- rep(14,length(headerpep_sn))
posM[headerpep_sn < 14] <- headerpep_sn[headerpep_sn < 14] -1 
NetChopOut_list <- lapply(seq(Npep),function(x){
	mt_ori <- NetChopOut[nstart[x]+4:c(lenpep[x]+3),1]
	mt <- t(sapply(strsplit(mt_ori,"  "),function(y) y[c(2:5)])) 
	mt[,2] <- sub(" ","",mt[,2])
	if((sum(mt[c(posM[x]+1):lenpep[x],3]=="S")>0))
	{
			upS <- c(0,which(mt[1:posM[x],3]=="S"))
			downS <- which(mt[c(posM[x]+1):lenpep[x],3]=="S")+posM[x]
			ud_combn <- matrix(as.numeric(unlist(strsplit(outer(upS,downS,paste)," "))),ncol=2,byrow=T)
			len <- ud_combn[,2]-ud_combn[,1]
			index <- which((len>=8)&(len<=14))
			if(length(index)>0){
				pep <- sapply(index,function(y) paste(mt[(ud_combn[y,1]+1):ud_combn[y,2],2],collapse=""))
				pep_df <- data.frame(header=Netfasta[2*(x-1)+1,1],peptide=pep,stringsAsFactors=F)
				return(pep_df)
			}
	}
})
NetChopOut_pep <- do.call(rbind,NetChopOut_list)
NetChopOut_pepU <- unique(unlist(NetChopOut_pep[,2]))

#write the result
setwd(paste0(chrDir,smp2))
write.table(NetChopOut_pep, file=paste0(smp,"NetChopOut_pep.txt"), quote=FALSE, row.names=F, col.names=T ) 
write.table(matrix(NetChopOut_pepU,ncol=1), file=paste0(smp,"NetChopcut.pep"), quote=FALSE, row.names=F, col.names=F ) 

#########################################################
#human leukocyte antigen (HLA) haplotyping by polysolver#
#########################################################
#because the reference files of polysover are under hg19,
#therefore matched normal cell is alse aligned to hg19
#${hg19RefFile} is the reference file for hg19
#${fastqFile1} is the fastq file for read 1
#${fastqFile2} is the fastq file for read 2
#${smp} is the sampleID
#${outputDir} is the output directory 
#${DIR} is the working directory for ploysolver
#${NAME} is the ID for matched normal sample
#${BAM} is the bam file of chromosome 6

#align by BWA 
bwa mem -o smp.sam -t 32 -M ${hg19RefFile}  ${fastqFile1} ${fastqFile2}; 
samtools view -@ 16 -bS ${smp}.sam > ${smp}.bam; 
samtools sort -@ 16 -o ${smp}.sort.bam ${smp}.bam; 
samtools index ${smp}.sort.bam; 
rm ${smp}.bam ${smp}.sam

#markduplicates by picards
picard.jar MarkDuplicates I=${smp}.sort.bam O=${smp}.sort.dedup.bam \
CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=. \
METRICS_FILE=${smp}.sort.metrics ASSUME_SORTED=true MAX_RECORDS_IN_RAM=8000000

#BaseRecalibrator by GATK
gatk BaseRecalibrator -I ${smp}.sort.dedup.bam -R ${hg19RefFile} \
--known-sites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
--known-sites dbsnp_138.hg19.vcf \
--known-sites 1000G_phase1.snps.high_confidence.hg19.sites.vcf \
--output ${smp}_recalibration_report.grp
gatk ApplyBQSR -R ${hg19RefFile} -I ${smp}.sort.dedup.bam \
--bqsr-recal-file ${smp}_recalibration_report.grp \
--output ${smp}.sort.dedup.baserecal.bam

#extract the result of chromosome 6
samtools view -h ${smp}.sort.dedup.baserecal.bam chr6 | sed 's/chr//g' > ${smp}.sort.dedup.baserecal.dechr.chr6.sam
samtools view -@ 16 -bS ${smp}.sort.dedup.baserecal.dechr.chr6.sam > ${smp}.sort.dedup.baserecal.dechr.chr6.bam
samtools sort -@ 16 -o ${smp}.sort.dedup.baserecal.dechr.chr6.sort.bam ${smp}.sort.dedup.baserecal.dechr.chr6.bam
samtools index ${smp}.sort.dedup.baserecal.dechr.chr6.sort.bam
rm ${smp}.sort.dedup.baserecal.dechr.chr6.bam ${smp}.sort.dedup.baserecal.dechr.chr6.sam

#run polysolver
docker run -d -P --name $NAME -v $DIR:/home/docker sachet/polysolver:v4 bash /home/polysolver/scripts/shell_call_hla_type /home/docker/$BAM Unknown 1 hg19 STDFQ 0 /home/docker 

################################################
#predict peptide-MHC binding event by netMHCpan#
################################################
#${smp} is sample ID
#${chrDir} is the directory of output
#${smp2} is sample ID for matched normal sample

#generate the script for netMHCpan
setwd(paste0(chrDir,smp2))
mhc_file <- read.table(file ="winners.hla.txt",sep="\t",quote ="",stringsAsFactors=F,fill=T,header=F)
mhc_list <- lapply(1:3,function(x) {
	a1 <- paste0(mhc_file[x,1],strsplit(mhc_file[x,2],"_")[[1]][3],":",strsplit(mhc_file[x,2],"_")[[1]][4])
	a2 <- paste0(mhc_file[x,1],strsplit(mhc_file[x,3],"_")[[1]][3],":",strsplit(mhc_file[x,3],"_")[[1]][4])
	return(c(a1,a2))
})
mhc_v <- paste(unique(unlist(mhc_list)),collapse=",")
code <- paste0("/NAS/dyl/software/netMHCpan/netMHCpan-4.1/netMHCpan -p ",
smp,"NetChopcut.pep -BA -xls -a ",mhc_v," -xlsfile ",smp,"_NetChopcut2netMHCpan.xls")
write.table(code, file=paste0(smp,".NetChopcut2netMHCpan.sh"), quote=FALSE, row.names=F, col.names=F )

#run netMHCpan
cd ${chrDir}${smp2}; 
sh ${chrDir}${smp2}/${smp}.NetChopcut2netMHCpan.sh > ${smp}.NetChopcut2netMHCpan.log 

#############
#Run Pyclone#
#############
#${DIR} is the working directory
#${smp} is the sample ID

conda activate pyclone
cd ${DIR}/${smp}
mkdir output
pathI=${DIR}/${smp}
PyClone setup_analysis --in_files ${pathI}/HKJ.tsv \
--working_dir ${DIR}/${smp}/output \
--config_extras_file ${DIR}/${smp}/config.txt

cd ${DIR}/${smp}/output
PyClone run_analysis --config_file config.yaml

PyClone build_table --config_file config.yaml \
--out_file ${DIR}/${smp}/output/out.tsv \
--table_type loci
conda deactivate

###############################
#Run clonevol(R version 4.0.3)#
###############################
#${smp} is the sample ID 
#${tumor} is the tumor ID of certain sample
#${DIR} is the working directory
#${inputDir} is the directory of input file, also the result directory of pyclone
#${annovarDir} is the directory of ANNOVAR results from WES
#${cancerGeneConsus} is the records from cancer gene consus

#Set working directory
setwd(DIR)
dir.create(smp)
setwd(paste0(DIR,smp))

#Load required packages
library(clonevol)

#Load input files
input <- read.table(file =paste0(inputDir,smp,"/output/out.tsv"),sep="\t",stringsAsFactors=F, header=T)
count.input <- read.table(file =paste0(inputDir,
smp,"/",smp,".tsv"),sep="\t",stringsAsFactors=F, header=T)
sample.ref.count <- count.input[match(input[,1],count.input[,1]),"ref_counts"]
sample.var.count <- count.input[match(input[,1],count.input[,1]),"var_counts"]
sample.depth <- sample.ref.count+sample.var.count
annova <- read.table(file =paste0(annovarDirtumor,".hg38_multianno.txt"),sep="\t",stringsAsFactors=F, header=T)
load(cancerGeneConsus)

#Format of input file
gene_label <- paste(annova[,1],annova[,2],sep=":")
gene_cgc <- annova[(annova[,"ExonicFunc.refGene"]=="nonsynonymous SNV")&(annova[,"Gene.refGene"]%in%cgc),"Gene.refGene"] 
names(gene_cgc) <- gene_label[(annova[,"ExonicFunc.refGene"]=="nonsynonymous SNV")&(annova[,"Gene.refGene"]%in%cgc)]
gene <- rep("-",nrow(input))
names(gene) <- input[,1]
gene[names(gene_cgc)] <- gene_cgc
names(gene) <- NULL
is.driver <- rep(F,nrow(input))
is.driver[gene!="-"] <- T
cluster <- input[,"cluster_id"]
cluster1 <- table(input[,"cluster_id"])
cluster_sel <- names(cluster1)[cluster1>=2]
res_df <- data.frame(cluster=cluster,sample.vaf=input[,"variant_allele_frequency"]*100,
sample.ccf=input[,"cellular_prevalence"]*100,gene=gene,is.driver=is.driver,
sample.ref.count=sample.ref.count,sample.var.count=sample.var.count,
sample.depth=sample.depth,sample=input[,"variant_allele_frequency"]*100,stringsAsFactors=F)
res_df <- res_df[res_df[,"cluster"]%in%cluster_sel,]
cluster_df <- data.frame(old=cluster_sel,new=seq(length(cluster_sel)),stringsAsFactors=F)
res_df[,"cluster"] <- cluster_df[match(res_df[,"cluster"],cluster_df[,"old"]),"new"]

#Run infer.clonal.models
if(length(cluster_sel)==1)
{
 y <- res_df
 save(y,file=paste0(smp,".clonal.models.RData"))
}else{
	vaf.col.names <- "sample"
	sample.names <- "sample"
	sample.groups <- "sample"
	names(sample.groups) <- vaf.col.names

	res_df <- res_df[order(res_df$cluster),]

	clone.colors <- NULL

	vaf_list <- split(res_df[,"sample.vaf"],f=factor(as.character(res_df[,"cluster"])))
	vaf_mean <- sapply(vaf_list,mean)
	founding.cluster <- as.numeric(names(vaf_mean)[which(vaf_mean==max(vaf_mean))[1]])

	y = infer.clonal.models(variants = res_df,
	cluster.col.name = 'cluster',
	vaf.col.names = vaf.col.names,
	sample.groups = sample.groups,
	cancer.initiation.model='monoclonal',
	subclonal.test = 'bootstrap',
	subclonal.test.model = 'non-parametric',
	num.boots = 1000,
	founding.cluster = founding.cluster,
	cluster.center = 'mean',
	ignore.clusters = NULL,
	clone.colors = clone.colors,
	min.cluster.vaf = 0.001,
	# min probability that CCF(clone) is non-negative
	sum.p = 0.05,
	# alpha level in confidence interval estimate for CCF(clone)
	alpha = 0.05)

	if(length(unique(res_df[,"gene"]))>1)
	{
	#Mapping driver events onto the trees
	y <- transfer.events.to.consensus.trees(y,
	res_df[res_df$is.driver,],
	cluster.col.name = 'cluster',
	event.col.name = 'gene')
	}

	#Plotting multiple plots and trees together
	plot.clonal.models(y,
	# box plot parameters
	box.plot = FALSE,
	fancy.boxplot = FALSE,
	# bell plot parameters
	clone.shape = 'bell',
	bell.event = TRUE,
	bell.event.label.color = 'blue',
	bell.event.label.angle = 60,
	clone.time.step.scale = 1,
	bell.curve.step = 2,
	cell.plot = FALSE,
	# branch-based consensus tree parameters
	merged.tree.clone.as.branch = FALSE,
	merged.tree.plot = FALSE,
	# output figure parameters
	out.dir = 'output',
	out.format = 'pdf',
	overwrite.output = TRUE,
	width = 7,
	height = 7)

	save(y,file=paste0(smp,".clonal.models.RData"))
}